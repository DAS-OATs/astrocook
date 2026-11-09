"""Each saved archive is a node of the analysis (F109, the rule of the
sessions): it holds A (the original data, ``*_orig``), B (the current state)
and the log from A to B. When an archive is opened, B becomes the new A and
the log starts empty: the log of the archive is shown read only while the
session is open, never replayed nor saved again.

Characterization tests were written first (the previous commit): the log
of an opened archive was continued, so that an archive opened and saved again
had B as its original data and the old log with it, and a replay ran the old
recipes a second time. They are now the tests of the rule.
"""
import asyncio
import json
import tarfile

import numpy as np
import pytest
from PySide6.QtWidgets import QMessageBox

from astrocook.core.session import load_session_from_file
from astrocook.core.session_manager import SessionHistory
from astrocook.core.structures import HistoryLogV2
from astrocook.gui.log_script import log_to_script
from astrocook.gui.main_window import MainWindowV2
from astrocook.methods_text import methods_text

from test_gui_smoke import SESSION_FILE, _run_recipe


@pytest.fixture
def node(tmp_path):
    """An archive saved after one recipe that is not idempotent: A is the
    spectrum of the tests, B has its flux doubled, the log has the recipe."""
    a = load_session_from_file(str(SESSION_FILE), "node", None, "auto")
    history = SessionHistory(a, HistoryLogV2())
    params = {"target_col": "y", "expression": "y * 2"}
    history.log_manager.add_entry("apply_expression", params, reason="To double it.")
    b = a.edit.apply_expression(**params)
    history.add_state(b)
    path = tmp_path / "node.acs2"
    b.save(str(path), initial_session=a)
    return path, a.spec.y.value.copy()


def _original_flux(path, tmp_path):
    """The flux of the original data (``*_spec_orig.fits``) of an archive."""
    from astropy.table import Table
    with tarfile.open(path) as tar:
        name = next(n for n in tar.getnames() if n.endswith("_spec_orig.fits"))
        tar.extract(name, tmp_path / "orig", filter="data")
    return np.asarray(Table.read(tmp_path / "orig" / name)["y"])


def _saved_log(path):
    with tarfile.open(path) as tar:
        name = next(n for n in tar.getnames() if n.endswith("_meta.json"))
        meta = json.load(tar.extractfile(name))
    return [e["recipe_name"] for e in meta["log_history_json"].get("entries", [])]


@pytest.fixture
def window(qtbot, monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    for kind in ("warning", "critical", "information", "question"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(
            lambda *a, **k: QMessageBox.StandardButton.Ok))
    win = MainWindowV2(None, None)
    qtbot.addWidget(win)
    monkeypatch.setattr(win, "_show_custom_message",
                        lambda *a, **k: QMessageBox.StandardButton.Yes)
    win.show()
    return win


# --- The reader of an archive: faithful to the file (unchanged) -------------------------------

def test_the_reader_gives_the_log_of_the_archive(node):
    path, _ = node
    session = load_session_from_file(str(path), "node", None, "auto")
    assert [e.recipe_name for e in session.log_manager.entries] == ["apply_expression"]


# --- An archive opened: B is the new A, and the log starts empty --------------------------------
# (characterization tests before the fix: the log was continued, the archive
# saved again had B as its original data and the old log with it, and Run All
# doubled the flux twice)

def test_the_gui_starts_a_new_log(qtbot, window, node, tmp_path):
    path, a_flux = node
    window.open_session_from_path(str(path))
    history = window.active_history
    log = history.log_manager
    assert log.entries == []
    assert [e.recipe_name for e in log.previous.entries] == ["apply_expression"]
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072"})
    saved = tmp_path / "again.acs2"
    history.current_state.save(str(saved), initial_session=history.states[0])
    # A node: the original data are B (the flux doubled), and the log goes
    # from them to the state: replayed on them, it gives the state
    assert np.allclose(_original_flux(saved, tmp_path), 2 * a_flux, equal_nan=True)
    assert _saved_log(saved) == ["set_properties"]


def test_the_log_scripter_shows_the_previous_log(window, node):
    path, _ = node
    window.open_session_from_path(str(path))
    lines = log_to_script(window.active_history.log_manager)
    assert lines[0].startswith("## The log of the archive this session was opened from")
    assert "## # reason: To double it." in lines
    assert any(line.startswith("## apply_expression(") for line in lines)
    assert not any(line.startswith("apply_expression(") for line in lines)


def test_run_all_does_not_replay_the_archive(qtbot, window, node):
    path, a_flux = node
    window.open_session_from_path(str(path))
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072"})
    before = window.active_history
    window.run_script("\n".join(log_to_script(before.log_manager)))
    qtbot.waitUntil(lambda: window.active_history is not before, timeout=20000)
    history = window.active_history
    assert len(history.states) == 2
    flux = history.current_state.spec.y.value
    assert np.allclose(flux, 2 * a_flux, equal_nan=True)  # doubled once
    assert history.current_state.spec._data.z_em == 4.072
    # The log of the archive stays, read only
    assert [e.recipe_name for e in history.log_manager.entries] == ["set_properties"]
    assert [e.recipe_name for e in history.log_manager.previous.entries] == [
        "apply_expression"]


def test_the_server_starts_a_new_log(node, tmp_path):
    pytest.importorskip("mcp")
    from mcp import Client
    from astrocook.mcp_server.server import build_server
    path, a_flux = node
    saved = tmp_path / "server.acs2"

    async def go():
        async with Client(build_server()) as client:
            await client.call_tool("session_open", {"path": str(path)})
            await client.call_tool("edit_set_properties",
                                   {"session_id": "node", "z_em": 4.072})
            described = await client.call_tool("session_describe", {"session_id": "node"})
            history = await client.call_tool("session_history", {"session_id": "node"})
            await client.call_tool("session_save", {"session_id": "node",
                                                    "path": str(saved)})
            return described, history
    described, history = asyncio.run(go())
    before = described.structured_content["log_before"]
    assert [e["recipe"] for e in before] == ["apply_expression"]
    assert before[0]["reason"] == "To double it."
    states = history.structured_content["states"]
    assert [s["recipe"] for s in states] == [None, "set_properties"]
    assert np.allclose(_original_flux(saved, tmp_path), 2 * a_flux, equal_nan=True)
    assert _saved_log(saved) == ["set_properties"]


def test_methods_text_of_an_opened_archive(window, node):
    path, _ = node
    window.open_session_from_path(str(path))
    text = methods_text(window.active_history).text
    assert "To double it." in text
    assert "For the first step, read from the log of the archive" in text


def test_the_original_data_have_their_metadata(node):
    path, _ = node
    with tarfile.open(path) as tar:
        name = next(n for n in tar.getnames() if n.endswith("_meta_orig.json"))
        meta = json.load(tar.extractfile(name))
    assert meta["log_history_json"] == {} and meta["log_type"] == "unknown"
    assert "constraints_by_uuid" in meta and "component_uuids" in meta
    # The reader takes _meta.json, not the metadata of the original data
    session = load_session_from_file(str(path), "node", None, "auto")
    assert [e.recipe_name for e in session.log_manager.entries] == ["apply_expression"]


def test_the_previous_log_follows_copies_and_branches(qtbot, window, node, monkeypatch):
    path, _ = node
    window.open_session_from_path(str(path))
    source = window.active_history
    # Duplicate: a new session, its own empty log, no previous
    window._on_duplicate_session(source)
    assert window.active_history.log_manager.previous is None
    # A branch keeps the previous log (a copy of the log of its parent)
    window.active_history = source
    _run_recipe(qtbot, window, "edit", "split", {"expression": "(x > 500) & (x < 520)"})
    branch = window.active_history
    assert branch is not source
    assert [e.recipe_name for e in branch.log_manager.entries] == ["split"]
    assert [e.recipe_name for e in branch.log_manager.previous.entries] == [
        "apply_expression"]


def test_opened_from_keeps_only_the_entries_done():
    log = HistoryLogV2()
    log.add_entry("a", {})
    log.add_note("after a")
    log.add_entry("b", {})
    log.add_note("after b")
    log.undo()
    log.previous = HistoryLogV2.opened_from(HistoryLogV2())  # not carried over
    new = HistoryLogV2.opened_from(log)
    assert new.entries == [] and new.current_index == -1 and new.notes == []
    assert [e.recipe_name for e in new.previous.entries] == ["a"]
    assert [n.text for n in new.previous.notes] == ["after a"]
    assert new.previous.previous is None
    assert HistoryLogV2.opened_from(HistoryLogV2()).previous is None
    assert HistoryLogV2.opened_from(None).previous is None
    # Never saved
    assert "previous" not in new.to_dict()
