"""A session followed by the GUI (phase 6-ter, PR 5-bis; decision 7 of the
phase, option b): the MCP server saves the session to an archive after every
change and the GUI reloads it, keeping the view; the followed session is read
only in the GUI, with a banner and an editable copy.

The first tests are characterization tests of the behaviour of an archive
opened in the GUI, written before the follow was added: an archive without a
marker must keep opening as before.
"""
import shutil

import pytest
from PySide6.QtWidgets import QMessageBox

from astrocook.gui.main_window import MainWindowV2

from test_gui_smoke import SESSION_FILE, _run_recipe


def _new_window(qtbot, monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    for kind in ("warning", "critical", "information", "question"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(
            lambda *a, **k: QMessageBox.StandardButton.Ok))
    win = MainWindowV2(None, None)
    qtbot.addWidget(win)
    monkeypatch.setattr(win, "_show_custom_message",
                        lambda *a, **k: QMessageBox.StandardButton.Ok)
    win.show()
    return win


@pytest.fixture
def archive(tmp_path):
    """A copy of the J0103-1305 archive of the tests, in its own folder."""
    path = tmp_path / "data" / "J0103.acs2"
    path.parent.mkdir()
    shutil.copy(SESSION_FILE, path)
    return path


@pytest.fixture
def window(qtbot, monkeypatch, tmp_path):
    return _new_window(qtbot, monkeypatch, tmp_path)


# --- Characterization: an archive without a marker -------------------------------------------

def test_an_archive_opens_editable(qtbot, window, archive):
    window.open_session_from_path(str(archive))
    history = window.active_history
    assert len(history.states) == 1
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072"})
    assert len(history.states) == 2
    assert window.undo_action.isEnabled()
    window._undo_last_action()
    assert history.current_index == 0


def test_a_duplicate_starts_with_a_new_log(qtbot, window, archive):
    window.open_session_from_path(str(archive))
    history = window.active_history
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072"})
    assert len(history.log_manager.entries) == 1
    window._on_duplicate_session(history)
    copy = window.active_history
    assert copy is not history and copy.log_manager.entries == []
    assert copy.current_state.name == f"{history.current_state.name}_copy1"
