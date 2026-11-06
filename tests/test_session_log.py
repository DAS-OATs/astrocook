"""The session log (``HistoryLogV2``): what is saved, how it is read back,
how it follows undo and redo, and how it replays.

Characterization written before phase 6-ter, PR 5a (the reasons of the
recipes and the notes of the agents in the log): the format of the log in
the ``.acs2`` archives and the relation between the states of a
``SessionHistory`` and the entries of its log must not change for the
logs that have neither.
"""
import json

import pytest

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


# --- Reasons and notes (phase 6-ter, PR 5a) ---------------------------------

def test_a_reason_is_saved_with_its_entry():
    log = HistoryLogV2()
    log.add_entry("smooth", {"sigma_kms": 50.0}, reason="  Noisy: smooth first. ")
    log.add_entry("rebin", {"dx": 20.0}, reason="   ")  # no reason
    assert [e.reason for e in log.entries] == ["Noisy: smooth first.", None]
    assert log.to_dict()["entries"] == [
        {"recipe_name": "smooth", "params": {"sigma_kms": 50.0},
         "reason": "Noisy: smooth first."},
        {"recipe_name": "rebin", "params": {"dx": 20.0}}]
    assert "notes" not in log.to_dict()


def test_notes_are_not_entries():
    log = _log("a")
    note = log.add_note("z_em from the Ly-a peak, computed outside Astrocook.")
    assert note.after == 1 and note.kind == "note"
    assert _recipes(log) == ["a"] and log.current_index == 0
    assert log.to_dict()["notes"] == [
        {"text": "z_em from the Ly-a peak, computed outside Astrocook.",
         "kind": "note", "after": 1}]


def test_a_note_needs_a_text_and_a_known_kind():
    log = HistoryLogV2()
    for text in ("", "  ", None):
        with pytest.raises(ValueError, match="needs a text"):
            log.add_note(text)
    with pytest.raises(ValueError, match="Unknown kind"):
        log.add_note("x", kind="gap")
    assert log.notes == []


def test_the_timeline_puts_notes_where_they_were_written():
    log = HistoryLogV2()
    log.add_note("before anything")
    log.add_entry("a", {})
    log.add_note("after a")
    log.add_note("after a, again")
    log.add_entry("b", {})
    items = [(getattr(i, "recipe_name", None) or i.text) for i, _ in log.timeline()]
    assert items == ["before anything", "a", "after a", "after a, again", "b"]
    assert all(active for _, active in log.timeline())


def test_notes_are_undone_and_redone_with_their_entry():
    log = _log("a", "b")
    log.add_note("on b")
    log.undo()
    assert "notes" not in log.to_dict()
    assert [(getattr(i, "recipe_name", None) or i.text, active)
            for i, active in log.timeline(include_undone=True)] == [
        ("a", True), ("b", False), ("on b", False)]
    log.redo()
    assert [n["text"] for n in log.to_dict()["notes"]] == ["on b"]


def test_a_note_after_undo_keeps_the_redo():
    log = _log("a", "b")
    log.undo()
    log.add_note("on a")
    assert log.redo() and _recipes(log) == ["a", "b"]
    assert [(getattr(i, "recipe_name", None) or i.text) for i, _ in log.timeline()] == [
        "a", "on a", "b"]


def test_a_new_entry_after_undo_drops_the_notes_undone():
    log = _log("a", "b")
    log.add_note("on b")
    log.undo()
    log.add_note("on a")
    log.add_entry("c", {})
    assert [n.text for n in log.notes] == ["on a"]
    assert [(getattr(i, "recipe_name", None) or i.text) for i, _ in log.timeline()] == [
        "a", "on a", "c"]


def test_notes_do_not_move_the_states_of_a_history():
    session = _session()
    history = SessionHistory(session, HistoryLogV2())
    history.log_manager.add_note("as opened")
    history.log_manager.add_entry("a", {})
    history.add_state(session.with_new_spectrum(session.spec))
    history.log_manager.add_note("on a")
    assert history.current_index == 1 and history.log_manager.current_index == 0
    history.undo()
    assert history.log_manager.current_index == -1
    assert [n.text for n, _ in history.log_manager.timeline()] == ["as opened"]
    history.redo()
    assert history.log_manager.current_index == 0
    assert len(history.log_manager.timeline()) == 3


def test_reasons_and_notes_are_saved_and_reopened(tmp_path):
    session = _session()
    log = HistoryLogV2()
    log.add_note("Opened as given.")
    log.add_entry("smooth", {"sigma_kms": 50.0}, reason="Too noisy to read.")
    log.add_note("Looked at the smoothed flux:\nno line above 3 sigma.")
    log.add_entry("rebin", {"dx": 20.0, "xunit": "km/s"})
    session.log_manager = log
    reloaded = _reload(session, tmp_path).log_manager
    assert reloaded.to_dict() == log.to_dict()
    assert [(getattr(i, "recipe_name", None) or i.text) for i, _ in reloaded.timeline()] == [
        "Opened as given.", "smooth",
        "Looked at the smoothed flux:\nno line above 3 sigma.", "rebin"]
    assert reloaded.entries[0].reason == "Too noisy to read."
    # The entries replay as before: the notes are not among them
    replayed = _session()
    for e in reloaded.to_dict()["entries"]:
        replayed = getattr(replayed.flux, e["recipe_name"])(**e["params"])
    expected = _session().flux.smooth(sigma_kms=50.0).flux.rebin(dx=20.0, xunit="km/s")
    assert np.array_equal(replayed.spec.y.value, expected.spec.y.value, equal_nan=True)


def test_notes_of_a_later_version_are_read():
    # A kind not known here, and an 'after' out of range, are kept as read
    # and brought within the log
    log = HistoryLogV2.from_dict({
        "entries": [{"recipe_name": "a", "params": {}}], "current_index": 0,
        "notes": [{"text": "late", "kind": "gap", "after": 7},
                  {"text": "early", "after": 0}]})
    assert [(n.text, n.kind, n.after) for n in log.notes] == [
        ("early", "note", 0), ("late", "gap", 1)]


def test_a_malformed_note_gives_an_empty_log(tmp_path, caplog):
    # As a malformed entry: the session opens, with an empty log and a warning
    session = _session()
    session.log_manager = _log("smooth")
    path = tmp_path / "bad.acs2"
    session.save(str(path))
    import tarfile
    with tarfile.open(path) as tar:
        members = [(m, tar.extractfile(m).read() if m.isfile() else None)
                   for m in tar.getmembers()]
    bad = tmp_path / "bad2.acs2"
    with tarfile.open(bad, "w:gz") as tar:
        for m, data in members:
            if m.name.endswith(".json"):
                meta = json.loads(data)
                meta["log_history_json"]["notes"] = [{"after": "x", "text": "t"}]
                data = json.dumps(meta).encode()
                m.size = len(data)
            import io
            tar.addfile(m, io.BytesIO(data) if data is not None else None)
    with caplog.at_level("WARNING"):
        reloaded = load_session_from_file(str(bad), "bad", None, "auto")
    assert reloaded and reloaded.log_manager is None
    assert "Could not restore the session log" in caplog.text


def test_the_log_scripter_shows_reasons_and_notes(qtbot):
    from astrocook.gui.log_scripter_dialog import LogScripterDialog
    log = HistoryLogV2()
    log.add_note("Opened as given.")
    log.add_entry("smooth", {"sigma_kms": 50.0}, reason="Too noisy.\nSee the plot.")
    log.add_entry("rebin", {"xunit": "km/s"}, reason="Coarser.")
    log.add_note("Looked at it.")
    log.undo()
    dialog = LogScripterDialog(log, "main", None, {})
    qtbot.addWidget(dialog)
    assert dialog.scripter_edit.toPlainText().splitlines() == [
        "# note: Opened as given.",
        "# reason: Too noisy.",
        "#         See the plot.",
        "smooth(sigma_kms=50.0)",
        "# # reason: Coarser.",
        "# rebin(xunit='km/s')",
        "# # note: Looked at it."]


def test_the_comments_of_a_script():
    from astrocook.gui.log_script import ScriptComments, log_to_script
    log = HistoryLogV2()
    log.add_note("First.\n  Indented line.")
    log.add_note("Second.")
    log.add_entry("smooth", {"sigma_kms": 50.0}, reason="Why.\nBecause.")
    log.add_entry("rebin", {})
    log.add_note("Last.")
    lines = log_to_script(log) + ["# a comment of the user", "", "#   not a note"]
    comments, read = ScriptComments(), []
    for line in lines:
        if line.startswith("#") or not line:
            comments.read(line)
        else:
            read.append((comments.notes(), comments.reason(), line))
    read.append((comments.notes(), comments.reason(), None))
    assert read == [
        ([("note", "First.\nIndented line."), ("note", "Second.")], "Why.\nBecause.",
         "smooth(sigma_kms=50.0)"),
        ([], None, "rebin()"),
        ([("note", "Last.")], None, None)]
