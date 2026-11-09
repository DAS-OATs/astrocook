"""Each saved archive is a node of the analysis (F109, the rule of the
sessions): it holds A (the original data, ``*_orig``), B (the current state)
and the log from A to B. When an archive is opened, B becomes the new A and
the log starts empty: the log of the archive is shown read only while the
session is open, never replayed nor saved again.

The tests marked "characterization" record the behaviour before the fix:
the log of an opened archive was continued, so that an archive opened and
saved again had B as its original data and the old log with it, and a replay
ran the old recipes a second time.
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


# --- Characterization: the log of an opened archive is continued ------------------------------

def test_characterization_the_gui_continues_the_log(qtbot, window, node, tmp_path):
    path, a_flux = node
    window.open_session_from_path(str(path))
    history = window.active_history
    assert [e.recipe_name for e in history.log_manager.entries] == ["apply_expression"]
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072"})
    saved = tmp_path / "again.acs2"
    history.current_state.save(str(saved), initial_session=history.states[0])
    # The original data of the new archive are B (flux doubled), and its log
    # still has the recipe that made B: a replay doubles the flux again
    assert np.allclose(_original_flux(saved, tmp_path), 2 * a_flux, equal_nan=True)
    assert _saved_log(saved) == ["apply_expression", "set_properties"]


def test_characterization_run_all_replays_the_archive(qtbot, window, node):
    path, a_flux = node
    window.open_session_from_path(str(path))
    script = "\n".join(log_to_script(window.active_history.log_manager))
    assert "apply_expression(" in script
    window.run_script(script)
    qtbot.waitUntil(lambda: len(window.active_history.states) == 2, timeout=20000)
    flux = window.active_history.current_state.spec.y.value
    assert np.allclose(flux, 4 * a_flux, equal_nan=True)  # doubled twice


def test_characterization_the_server_saves_the_archive_log(node, tmp_path):
    pytest.importorskip("mcp")
    from mcp import Client
    from astrocook.mcp_server.server import build_server
    path, a_flux = node
    saved = tmp_path / "server.acs2"

    async def go():
        async with Client(build_server()) as client:
            opened = await client.call_tool("session_open", {"path": str(path)})
            described = await client.call_tool("session_describe", {"session_id": "node"})
            await client.call_tool("edit_set_properties",
                                   {"session_id": "node", "z_em": 4.072})
            await client.call_tool("session_save", {"session_id": "node",
                                                    "path": str(saved)})
            return opened, described
    _, described = asyncio.run(go())
    before = described.structured_content["log_before"]
    assert [e["recipe"] for e in before] == ["apply_expression"]
    assert before[0]["reason"] == "To double it."
    assert np.allclose(_original_flux(saved, tmp_path), 2 * a_flux, equal_nan=True)
    assert _saved_log(saved) == ["apply_expression", "set_properties"]


def test_characterization_methods_text_of_an_opened_archive(window, node):
    path, _ = node
    window.open_session_from_path(str(path))
    text = methods_text(window.active_history).text
    assert "To double it." in text
    assert "For the first step, read from the log of the archive" in text


def test_characterization_the_original_data_have_no_metadata(node):
    path, _ = node
    with tarfile.open(path) as tar:
        names = tar.getnames()
    assert any(n.endswith("_spec_orig.fits") for n in names)
    assert not any(n.endswith("_meta_orig.json") for n in names)
