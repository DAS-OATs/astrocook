"""The session log (``HistoryLogV2``): what is saved, how it is read back,
how it follows undo and redo, and how it replays.

Characterization written before phase 6-ter, PR 5a (the reasons of the
recipes and the notes of the agents in the log): the format of the log in
the ``.acs2`` archives and the relation between the states of a
``SessionHistory`` and the entries of its log must not change for the
logs that have neither.
"""
import json

import astropy.units as au
import numpy as np

from astrocook.core.session import SessionV2, load_session_from_file
from astrocook.core.session_manager import SessionHistory
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, HistoryLogV2, SpectrumDataV2

FLUX = au.erg / au.cm**2 / au.s / au.nm


def _session(name="main"):
    """A noisy power law, z_em = 2.5."""
    rng = np.random.default_rng(0)
    x = np.arange(400.0, 440.0, 0.01)
    cont = (x / 420.0) ** -1.5
    y = cont + 0.02 * cont * rng.standard_normal(len(x))
    spec = SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - 0.005, au.nm),
        xmax=DataColumnV2(x + 0.005, au.nm), y=DataColumnV2(y, FLUX),
        dy=DataColumnV2(0.02 * cont, FLUX), aux_cols={}, z_em=2.5))
    return SessionV2(name=name, gui=None, spec=spec)


def _log(*names):
    log = HistoryLogV2()
    for n in names:
        log.add_entry(n, {"p": n})
    return log


def _recipes(log):
    return [e.recipe_name for e in log.entries]


# --- The format ----------------------------------------------------------

def test_to_dict_of_a_log():
    log = HistoryLogV2()
    log.add_entry("smooth", {"sigma_kms": 50.0})
    log.add_entry("rebin", {"dx": 20.0, "xunit": "km/s"})
    assert log.to_dict() == {
        "entries": [{"recipe_name": "smooth", "params": {"sigma_kms": 50.0}},
                    {"recipe_name": "rebin", "params": {"dx": 20.0, "xunit": "km/s"}}],
        "current_index": 1}


def test_to_dict_of_an_empty_log():
    assert HistoryLogV2().to_dict() == {"entries": [], "current_index": -1}


def test_from_dict_reads_name_and_params_only():
    # What the reader of 2.5.0 does with keys it does not know: ignores them
    # (an archive written by a later version opens)
    log = HistoryLogV2.from_dict({
        "entries": [{"recipe_name": "smooth", "params": {"sigma_kms": 50.0},
                     "something_new": "x"}],
        "current_index": 0, "another_key": [1, 2]})
    assert _recipes(log) == ["smooth"] and log.current_index == 0
    assert log.entries[0].params == {"sigma_kms": 50.0}


def test_from_dict_of_nothing():
    for data in (None, {}, {"entries": None}):
        log = HistoryLogV2.from_dict(data)
        assert log.entries == [] and log.current_index == -1


# --- Undo and redo --------------------------------------------------------

def test_undo_redo_move_the_index():
    log = _log("a", "b", "c")
    assert log.undo() and log.undo()
    assert log.current_index == 0
    assert log.redo() and log.current_index == 1
    assert log.to_dict()["entries"] == [{"recipe_name": "a", "params": {"p": "a"}},
                                        {"recipe_name": "b", "params": {"p": "b"}}]


def test_undo_to_the_start_and_beyond():
    log = _log("a")
    assert log.undo() and log.current_index == -1
    assert not log.undo()
    assert log.to_dict() == {"entries": [], "current_index": -1}


def test_redo_at_the_end():
    log = _log("a", "b")
    assert not log.redo() and log.current_index == 1


def test_a_new_entry_after_undo_drops_the_redo():
    log = _log("a", "b", "c")
    log.undo()
    log.undo()
    log.add_entry("d", {})
    assert _recipes(log) == ["a", "d"] and log.current_index == 1
    assert not log.redo()


def test_the_history_and_its_log_move_together():
    session = _session()
    history = SessionHistory(session, HistoryLogV2())
    for name in ("a", "b"):
        history.log_manager.add_entry(name, {})
        history.add_state(history.current_state.with_new_spectrum(
            history.current_state.spec))
    assert history.current_index == 2 and history.log_manager.current_index == 1
    history.undo()
    assert history.current_index == 1 and history.log_manager.current_index == 0
    history.undo()
    assert history.current_index == 0 and history.log_manager.current_index == -1
    assert history.undo() is None and history.log_manager.current_index == -1
    history.redo()
    assert history.current_index == 1 and history.log_manager.current_index == 0
    # A new state after an undo: the log is truncated by add_entry, the
    # states by add_state
    history.log_manager.add_entry("c", {})
    history.add_state(history.current_state.with_new_spectrum(history.current_state.spec))
    assert _recipes(history.log_manager) == ["a", "c"]
    assert history.current_index == 2 and history.log_manager.current_index == 1
    assert all(s.log_manager is history.log_manager for s in history.states)


# --- Saved and opened ------------------------------------------------------

def _reload(session, tmp_path, name="log"):
    path = tmp_path / f"{name}.acs2"
    session.save(str(path))
    return load_session_from_file(str(path), name, None, "auto")


def _saved_log(path):
    import tarfile
    with tarfile.open(path) as tar:
        member = next(m for m in tar.getmembers() if m.name.endswith(".json"))
        return json.load(tar.extractfile(member))["log_history_json"]


def test_the_archive_holds_to_dict(tmp_path):
    session = _session()
    log = _log("smooth", "rebin")
    session.log_manager = log
    path = tmp_path / "s.acs2"
    session.save(str(path))
    assert _saved_log(path) == log.to_dict()


def test_a_log_reopened_is_the_same(tmp_path):
    session = _session()
    session.log_manager = _log("smooth", "rebin")
    reloaded = _reload(session, tmp_path)
    assert reloaded.log_manager.to_dict() == session.log_manager.to_dict()
    assert reloaded.log_manager.current_index == 1


def test_the_saved_entries_replay(tmp_path):
    session = _session()
    log = HistoryLogV2()
    log.add_entry("smooth", {"sigma_kms": 50.0})
    log.add_entry("rebin", {"dx": 20.0, "xunit": "km/s"})
    session.log_manager = log
    reloaded = _reload(session, tmp_path)
    replayed = _session()
    for e in reloaded.log_manager.to_dict()["entries"]:
        replayed = getattr(replayed.flux, e["recipe_name"])(**e["params"])
    expected = _session().flux.smooth(sigma_kms=50.0).flux.rebin(dx=20.0, xunit="km/s")
    assert np.array_equal(replayed.spec.y.value, expected.spec.y.value, equal_nan=True)


# --- Shown in the GUI ------------------------------------------------------

def test_the_log_scripter_shows_the_entries(qtbot):
    from astrocook.gui.log_scripter_dialog import LogScripterDialog
    log = HistoryLogV2()
    log.add_entry("set_properties", {"z_em": 2.5, "resol": None})
    log.add_entry("smooth", {"sigma_kms": 50.0})
    log.add_entry("rebin", {"xunit": "km/s"})
    log.undo()
    dialog = LogScripterDialog(log, "main", None, {})
    qtbot.addWidget(dialog)
    # The entries undone are commented out (Run All skips them)
    assert dialog.scripter_edit.toPlainText().splitlines() == [
        "set_properties(z_em=2.5, resol=None)",
        "smooth(sigma_kms=50.0)",
        "# rebin(xunit='km/s')"]
