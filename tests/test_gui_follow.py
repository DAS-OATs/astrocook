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
from PySide6.QtWidgets import QFileDialog, QMessageBox

from astrocook import follow
from astrocook.session.session import load_session_from_file
from astrocook.gui.main_window import MainWindow

from test_gui_smoke import SESSION_FILE, _run_recipe


def _new_window(qtbot, monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    for kind in ("warning", "critical", "information", "question"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(
            lambda *a, **k: QMessageBox.StandardButton.Ok))
    win = MainWindow(None, None)
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


# --- A followed archive ------------------------------------------------------------------------

def _rewrite(archive, z_em):
    """Save the archive again, as the server does after a recipe: the session
    with z_em changed and the recipe in its log."""
    from astrocook.session.history import SessionHistory
    from astrocook.core.structures import HistoryLog
    session = load_session_from_file(str(SESSION_FILE), "J0103", None, "auto")
    history = SessionHistory(session, HistoryLog())
    history.log_manager.add_entry("set_properties", {"z_em": z_em})
    history.add_state(session.with_new_spectrum(session.spec.with_properties(z_em=z_em)))
    follow.write_archive(history, archive)


@pytest.fixture
def followed(qtbot, window, archive):
    """The window with the archive opened while a writer (this process)
    follows it."""
    follow.write_marker(archive, "J0103")
    window.open_session_from_path(str(archive))
    qtbot.wait(50)
    return window


def _z_em(window):
    return window.active_history.current_state.spec._data.z_em


def _asked(window, monkeypatch, button="Cancel"):
    asked = []

    def show(*args, **kwargs):
        asked.append(kwargs.get("header"))
        return getattr(QMessageBox.StandardButton, button)
    monkeypatch.setattr(window, "_show_custom_message", show)
    return asked


def test_a_marked_archive_is_followed(followed):
    assert followed.is_followed()
    assert not followed.follow_banner.isHidden()
    assert "Followed session." in followed.follow_banner_label.text()
    assert not followed.undo_action.isEnabled() and not followed.redo_action.isEnabled()
    # The cards float below the banner
    followed._reposition_floating_widgets()
    assert followed.left_sidebar_widget.y() >= followed._banner_room() > 0


def test_the_banner_names_the_writer(followed):
    """The pid and the version of Astrocook of the marker: which process
    writes the archive (a server of an older version left running)."""
    import os
    from astrocook import __version__
    text = followed.follow_banner_label.text()
    assert f"pid {os.getpid()}, Astrocook {__version__}" in text
    assert f"pid {os.getpid()}" in followed.follow_banner_label.toolTip()


def test_the_banner_of_a_marker_without_the_new_fields_is_as_before(
        qtbot, window, archive):
    follow.marker_path(archive).write_text('{"session_id": "J0103"}')
    window.open_session_from_path(str(archive))
    qtbot.wait(50)
    assert window.is_followed()
    assert window.follow_banner_label.text() == (
        "<b>Followed session.</b> Updated by an AI agent after every change. "
        "Read only.")
    assert "program that writes" not in window.follow_banner_label.toolTip()


def test_the_banner_follows_the_session_shown(qtbot, followed):
    followed.open_session_from_path(str(SESSION_FILE))  # not followed
    assert followed.follow_banner.isHidden() and not followed.is_followed()
    followed.active_history = followed.session_histories[0]
    followed._update_view_for_session(followed.active_history.current_state)
    assert not followed.follow_banner.isHidden()


def test_the_archive_is_reloaded_keeping_the_view(qtbot, followed, archive):
    history = followed.active_history
    followed._on_view_system_inspector()
    inspector = followed.system_inspector
    inspector.table_view.selectRow(5)
    selected = inspector._get_selected_components()[0].uuid
    followed.xmin_input.setText("500"); followed.xmax_input.setText("520")
    followed._on_set_custom_limits()
    xlim = followed.plot_viewer.canvas.axes.get_xlim()

    _rewrite(archive, 4.072)
    qtbot.waitUntil(lambda: _z_em(followed) == 4.072, timeout=10000)
    assert followed.active_history is history and len(history.states) == 1
    # The log of the agent, read only: the state read is the start (F109)
    assert history.log_manager.entries == []
    assert [e.recipe_name for e in history.log_manager.previous.entries] == ["set_properties"]
    assert history.current_state.log_manager is history.log_manager
    assert followed.plot_viewer.canvas.axes.get_xlim() == pytest.approx(xlim)
    assert inspector.isVisible()
    assert inspector._get_selected_components()[0].uuid == selected
    assert "updated by the agent" in followed.status_recipe.text()
    assert followed.is_followed()


def test_only_the_archive_reloads_and_once_per_burst(qtbot, followed, archive, monkeypatch):
    reloads = []
    reload = followed._reload_followed
    monkeypatch.setattr(followed, "_reload_followed",
                        lambda f, stamp: (reloads.append(stamp), reload(f, stamp)))
    # A temporary file of the writer, in the same folder
    (archive.parent / ".J0103.acs2.tmp").write_bytes(b"half")
    qtbot.wait(800)
    assert reloads == []
    for z in (4.0, 4.05, 4.072):
        _rewrite(archive, z)
    qtbot.waitUntil(lambda: _z_em(followed) == 4.072, timeout=10000)
    qtbot.wait(500)
    assert len(reloads) == 1


def test_a_failed_load_keeps_the_state(qtbot, followed, archive, caplog):
    import logging
    before = followed.active_history.current_state
    with caplog.at_level(logging.WARNING):
        archive.write_bytes(b"not an archive")
        qtbot.waitUntil(lambda: any("could not be reloaded" in r.getMessage()
                                    for r in caplog.records), timeout=10000)
    assert followed.active_history.current_state is before
    # The next archive is read
    _rewrite(archive, 4.072)
    qtbot.waitUntil(lambda: _z_em(followed) == 4.072, timeout=10000)


def test_a_followed_session_is_read_only(qtbot, followed, monkeypatch):
    history = followed.active_history
    asked = _asked(followed, monkeypatch)
    followed._on_recipe_requested("edit", "set_properties", {"z_em": "4.0"}, {})
    followed._launch_recipe_dialog("absorbers", "refit_all")
    assert followed.active_recipe_dialog is None
    followed._toggle_continuum_editor()
    assert not getattr(followed.plot_viewer, "_edit_mode_active", False)
    followed._undo_last_action()
    followed._redo_last_action()
    followed.run_script("pass")
    assert len(asked) == 6 and all(a == "'J0103' is read only." for a in asked)
    qtbot.wait(200)
    assert len(history.states) == 1 and len(followed.session_histories) == 1
    # Renaming it in the list: the name goes back
    index = followed.session_model.index(0)
    followed.session_model.setData(index, "renamed")
    assert followed.session_model.data(index) == "J0103"
    assert history.current_state.name == "J0103"


def test_saving_exporting_and_branching_are_allowed(qtbot, followed, monkeypatch, tmp_path):
    asked = _asked(followed, monkeypatch)
    assert not followed._refuse_recipe_if_followed("edit", "split")
    assert not followed._refuse_recipe_if_followed("file", "export_ascii")
    saved = tmp_path / "saved.acs2"
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(saved), "")))
    assert followed._save_session(followed.active_history) and saved.is_file()
    _run_recipe(qtbot, followed, "edit", "split", {"expression": "(x > 500) & (x < 520)"})
    assert asked == []
    branch = followed.active_history
    assert branch is followed.session_histories[1] and not followed.is_followed(branch)
    assert followed.follow_banner.isHidden()


def test_the_editable_copy(qtbot, followed, archive, monkeypatch):
    _rewrite(archive, 4.072)
    qtbot.waitUntil(lambda: _z_em(followed) == 4.072, timeout=10000)
    source = followed.active_history
    asked = _asked(followed, monkeypatch, "Yes")
    followed._on_recipe_requested("edit", "set_properties", {"z_em": "4.0"}, {})
    assert len(asked) == 1
    copy = followed.active_history
    assert copy is not source and not followed.is_followed(copy)
    assert copy.current_state.name == "J0103_copy1"
    # With the log of the agent, read only and not shared; its own starts
    # empty (F109: a new node, from the state of the agent)
    assert copy.log_manager.entries == []
    assert [e.recipe_name for e in copy.log_manager.previous.entries] == ["set_properties"]
    assert copy.log_manager is not source.log_manager
    assert followed.follow_banner.isHidden()
    _run_recipe(qtbot, followed, "edit", "set_properties", {"z_em": "4.0"})
    assert len(copy.states) == 2 and len(source.states) == 1
    assert [e.recipe_name for e in copy.log_manager.entries] == ["set_properties"]
    # The button of the banner makes another one
    followed.active_history = source
    followed.follow_copy_button.click()
    assert followed.active_history.current_state.name == "J0103_copy2"


def test_the_follow_ends_with_the_marker(qtbot, followed, archive):
    follow.remove(archive)
    qtbot.waitUntil(lambda: "No longer followed." in followed.follow_banner_label.text(),
                    timeout=10000)
    assert followed.is_followed()  # still read only
    assert followed._follow_watcher.files() == []
    assert "no longer followed" in followed.status_recipe.text()


def test_opened_again_it_is_shown(qtbot, followed, archive):
    followed.open_session_from_path(str(SESSION_FILE))
    followed.open_session_from_path(str(archive))
    assert len(followed.session_histories) == 2
    assert followed.active_history is followed.session_histories[0]
    assert not followed.follow_banner.isHidden()


def test_a_marker_of_a_dead_writer(qtbot, window, archive, monkeypatch):
    monkeypatch.setattr(follow.FollowMarker, "writer_alive", lambda self: False)
    follow.write_marker(archive, "J0103")
    window.open_session_from_path(str(archive))
    assert window.is_followed()
    assert "No longer followed." in window.follow_banner_label.text()
    assert window._follow_watcher.files() == []


def test_closing_a_followed_session(qtbot, followed):
    followed._on_close_session_requested()
    assert followed._followed == []
    assert followed._follow_watcher.files() == [] and followed._follow_watcher.directories() == []
    assert followed.follow_banner.isHidden()


def test_followed_from_the_command_line(qtbot, monkeypatch, tmp_path, archive):
    """The sessions given to the window at the start (the launcher sets the
    path they were read from)."""
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    follow.write_marker(archive, "J0103")
    session = load_session_from_file(str(archive), "J0103", None, "auto")
    session._source_path = str(archive)
    win = MainWindow([session], None)
    qtbot.addWidget(win)
    assert win.is_followed() and not win.follow_banner.isHidden()
    assert "_source_path" not in session.__dict__
