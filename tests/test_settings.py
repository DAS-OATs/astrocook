"""The *Settings...* dialog and what each setting does (phase 7-quater, U2).

The defaults keep the behaviour of the previous version; a setting is about
the interface only (no number, no script, no MCP tool, no log depends on it).
"""
import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QEvent
from PySide6.QtGui import QAction
from PySide6.QtWidgets import QApplication, QDialogButtonBox, QFileDialog, QMessageBox

from astrocook.gui import pyside_plot, software_update as su, theme
from astrocook.gui import settings as S
from astrocook.gui.data_inspector import format_value
from astrocook.gui.main_window import MainWindow
from astrocook.gui.settings import app_settings
from astrocook.gui.settings_dialog import SettingsDialog

from test_software_update import RELEASE_JSON, fake_urlopen

SESSION_FILE = Path(__file__).parent / "data" / "J0103-1305_metallines_lya_new.acs2"


@pytest.fixture(autouse=True)
def _no_stale_windows(qtbot):
    """Windows that other test files leave behind (an inspector whose main
    window is gone) would receive the font and palette change that Apply
    broadcasts, and fail in their ``resizeEvent`` inside one of these tests."""
    for widget in QApplication.topLevelWidgets():
        widget.close()
        widget.deleteLater()
    QApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


@pytest.fixture
def window(qtbot, monkeypatch):
    for kind in ("warning", "critical", "information", "question"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(
            lambda *a, **k: QMessageBox.StandardButton.Ok))
    win = MainWindow(None, None)
    qtbot.addWidget(win)
    win.show()
    return win


@pytest.fixture
def loaded(qtbot, window):
    window.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(100)
    return window


# --- Defaults and readers (characterization: the values of before) ---

def test_defaults_are_the_values_of_before():
    assert S.max_recent_files() == 10
    assert S.confirm_close_unsaved() is True
    assert S.initial_x_unit() == "nm"
    assert S.table_decimals() == 5
    assert S.line_width_factor() == 1.0
    assert S.get_bool(S.CHECK_UPDATES_KEY) is False  # no network unless asked
    assert S.downsample_threshold(pyside_plot.DECIMATION_THRESHOLD) == 1000000
    assert S.DEFAULTS[S.DOWNSAMPLE_KEY] == pyside_plot.DECIMATION_THRESHOLD


def test_the_old_keys_are_unchanged():
    assert (S.APPEARANCE_KEY, S.PLOT_COLOURS_KEY, S.RECENT_FILES_KEY, S.LAST_DIR_KEY) == (
        "appearance", "plot_colours", "recent_files", "last_working_dir")
    assert (S.ORGANIZATION, S.APPLICATION) == ("Astrocook", "AstrocookV2")


@pytest.mark.parametrize("stored, expected", [
    ("true", True), ("false", False), (True, True), (False, False), ("maybe", True), (3, True)])
def test_a_bool_comes_back_from_the_ini_file(stored, expected):
    app_settings().setValue(S.CONFIRM_CLOSE_KEY, stored)
    assert S.confirm_close_unsaved() is expected


def test_wrong_values_fall_back_to_the_defaults_or_the_limits():
    s = app_settings()
    s.setValue(S.MAX_RECENT_KEY, "many")
    s.setValue(S.DECIMALS_KEY, 99)
    s.setValue(S.INITIAL_X_UNIT_KEY, "parsec")
    s.setValue(S.LINE_WIDTH_KEY, "0.1")
    assert S.max_recent_files() == 10
    assert S.table_decimals() == S.DECIMALS_RANGE[1]
    assert S.initial_x_unit() == "nm"
    assert S.line_width_factor() == S.LINE_WIDTH_RANGE[0]


def test_starting_folder(tmp_path):
    s = app_settings()
    assert S.starting_folder()  # the working folder, as before
    s.setValue(S.LAST_DIR_KEY, str(tmp_path))
    assert S.starting_folder() == str(tmp_path)
    fixed = tmp_path / "fixed"
    fixed.mkdir()
    s.setValue(S.DIR_MODE_KEY, "fixed")
    assert S.starting_folder() == str(tmp_path)  # no folder chosen yet
    s.setValue(S.FIXED_DIR_KEY, str(fixed))
    assert S.starting_folder() == str(fixed)
    s.setValue(S.FIXED_DIR_KEY, str(tmp_path / "gone"))  # no longer there
    assert S.starting_folder() == str(tmp_path)


@pytest.mark.parametrize("frequency, hours, due", [
    ("start", 0, True), ("daily", 23, False), ("daily", 25, True),
    ("weekly", 24 * 6, False), ("weekly", 24 * 8, True)])
def test_update_check_due(frequency, hours, due):
    s = app_settings()
    now = datetime(2026, 11, 25, 12, 0)
    assert not S.update_check_due(now)  # off by default
    s.setValue(S.CHECK_UPDATES_KEY, True)
    s.setValue(S.CHECK_FREQUENCY_KEY, frequency)
    assert S.update_check_due(now)  # never checked
    S.record_update_check(now - timedelta(hours=hours))
    assert S.update_check_due(now) is due


# --- The dialog ---

def test_the_menu_entry(window):
    assert window.settings_action.menuRole() == QAction.MenuRole.PreferencesRole
    assert window.settings_action.text() == "&Settings..."


def test_the_dialog_shows_the_defaults(window):
    dlg = SettingsDialog(window)
    assert [dlg.tabs.tabText(i) for i in range(3)] == ["General", "Appearance", "Performance"]
    assert not dlg.check_updates.isChecked()
    assert dlg.max_recent.value() == 10 and dlg.confirm_close.isChecked()
    assert dlg.dir_mode.currentData() == "last" and not dlg.fixed_row.isEnabled()
    assert dlg.theme.currentData() == "system" and dlg.plot_colours.currentData() == "standard"
    assert dlg.font_scale.currentData() == "normal" and dlg.line_width.value() == 1.0
    assert dlg.x_unit.currentData() == "nm" and dlg.decimals.value() == 5
    assert dlg.downsample.value() == pyside_plot.DECIMATION_THRESHOLD


def test_changes_are_saved_and_the_next_dialog_reads(window, tmp_path):
    dlg = SettingsDialog(window)
    dlg.check_updates.setChecked(True)
    dlg.check_frequency.setCurrentIndex(dlg.check_frequency.findData("weekly"))
    dlg.max_recent.setValue(3)
    dlg.confirm_close.setChecked(False)
    dlg.dir_mode.setCurrentIndex(dlg.dir_mode.findData("fixed"))
    dlg.fixed_dir.setText(str(tmp_path))
    dlg.theme.setCurrentIndex(dlg.theme.findData("dark"))
    dlg.plot_colours.setCurrentIndex(dlg.plot_colours.findData("colorblind"))
    dlg.font_scale.setCurrentIndex(dlg.font_scale.findData("large"))
    dlg.line_width.setValue(2.0)
    dlg.x_unit.setCurrentIndex(dlg.x_unit.findData("Angstrom"))
    dlg.decimals.setValue(3)
    dlg.downsample.setValue(500_000)
    s = app_settings()
    assert S.get_bool(S.CHECK_UPDATES_KEY) and S.get_choice(S.CHECK_FREQUENCY_KEY, S.CHECK_FREQUENCIES) == "weekly"
    assert (S.max_recent_files(), S.confirm_close_unsaved()) == (3, False)
    assert S.starting_folder() == str(tmp_path)
    assert s.value("appearance") == "dark" and s.value("plot_colours") == "colorblind"
    assert (S.line_width_factor(), S.initial_x_unit(), S.table_decimals()) == (2.0, "Angstrom", 3)
    assert S.downsample_threshold(1) == 500_000
    again = SettingsDialog(window)
    assert again.max_recent.value() == 3 and again.theme.currentData() == "dark"
    assert again.fixed_row.isEnabled() and again.check_frequency.isEnabled()
    assert again.font_scale.currentData() == "large" and again.downsample.value() == 500_000


def test_the_buttons(window):
    dlg = SettingsDialog(window)
    SB = QDialogButtonBox.StandardButton
    assert dlg.buttons.standardButtons() == SB.Close | SB.RestoreDefaults


def test_a_change_applies_without_ok(window):
    dlg = SettingsDialog(window)
    assert app_settings().value("max_recent_files") is None
    dlg.max_recent.setValue(2)
    assert S.max_recent_files() == 2
    dlg.confirm_close.setChecked(False)
    assert not S.confirm_close_unsaved()
    dlg.x_unit.setCurrentIndex(dlg.x_unit.findData("micron"))
    assert S.initial_x_unit() == "micron"


def test_loading_the_dialog_saves_nothing(window):
    SettingsDialog(window)
    assert app_settings().value("max_recent_files") is None


def test_a_typed_number_applies_when_finished(window):
    dlg = SettingsDialog(window)
    assert not dlg.max_recent.keyboardTracking()
    dlg.fixed_dir.setText("/somewhere")  # typing: not yet
    assert app_settings().value(S.FIXED_DIR_KEY) is None
    dlg.fixed_dir.editingFinished.emit()
    assert app_settings().value(S.FIXED_DIR_KEY) == "/somewhere"


def test_restore_defaults_applies(window):
    app_settings().setValue(S.MAX_RECENT_KEY, 4)
    dlg = SettingsDialog(window)
    assert dlg.max_recent.value() == 4
    dlg.restore_defaults()
    assert dlg.max_recent.value() == 10 and dlg.downsample.value() == 1000000
    assert S.max_recent_files() == 10


def test_the_theme_stays_in_sync_with_the_menu(window):
    dlg = SettingsDialog(window)
    dlg.theme.setCurrentIndex(dlg.theme.findData("light"))
    dlg.plot_colours.setCurrentIndex(dlg.plot_colours.findData("colorblind"))
    dlg.apply()
    assert window.appearance_actions["light"].isChecked()
    assert window.plot_colours_actions["colorblind"].isChecked()
    window._set_appearance("system")  # from the View menu
    assert SettingsDialog(window).theme.currentData() == "system"


@pytest.fixture
def fonts(monkeypatch):
    """Fresh state of the font scaling, and the fonts back afterwards."""
    app = QApplication.instance()
    saved = {None: app.font(), **{n: app.font(n) for n in theme._FONT_CLASSES}}
    monkeypatch.setattr(theme, "_original_fonts", None)
    yield app
    for name, font in saved.items():
        app.setFont(font) if name is None else app.setFont(font, name)


def test_normal_never_calls_setfont(window, fonts, monkeypatch):
    calls = []
    monkeypatch.setattr(QApplication, "setFont", lambda *a, **k: calls.append(a))
    theme.apply_font_scale(1.0)
    window._apply_settings()
    window._restore_view_settings()
    assert calls == []


def test_the_font_scale(window, fonts):
    sizes = {n: fonts.font(n).pointSizeF() for n in ("QLabel", "QToolButton", "QHeaderView")}
    base = fonts.font().pointSizeF()
    dlg = SettingsDialog(window)
    dlg.font_scale.setCurrentIndex(dlg.font_scale.findData("large"))
    assert app_settings().value(S.FONT_SCALE_KEY) == "large"
    assert fonts.font().pointSizeF() == pytest.approx(base * 1.2)
    for name, size in sizes.items():  # the relations between classes stay
        assert fonts.font(name).pointSizeF() == pytest.approx(size * 1.2)
    dlg.font_scale.setCurrentIndex(dlg.font_scale.findData("small"))
    assert fonts.font("QLabel").pointSizeF() == pytest.approx(sizes["QLabel"] * 0.85)
    dlg.font_scale.setCurrentIndex(dlg.font_scale.findData("normal"))
    assert fonts.font().pointSizeF() == base
    for name, size in sizes.items():
        assert fonts.font(name).pointSizeF() == size


def test_the_headers_and_the_labels_scale_together(window, fonts):
    from PySide6.QtWidgets import QLabel
    from astrocook.gui.floating_cards import CollapsibleSection
    theme.apply_font_scale(1.2)
    section = CollapsibleSection("Group")
    section.show()
    label = QLabel("x")
    assert section.header.font().pointSizeF() == pytest.approx(label.font().pointSizeF())


def test_the_old_font_size_key_is_ignored(window, fonts):
    app_settings().setValue("font_size", 14)
    assert S.font_scale_factor() == 1.0
    window._apply_settings()
    assert SettingsDialog(window).font_scale.currentData() == "normal"
    app_settings().setValue(S.FONT_SCALE_KEY, "huge")  # not a choice either
    assert S.font_scale_factor() == 1.0


# --- The effect of each setting ---

def test_recent_files_limit(window, tmp_path):
    files = []
    for i in range(5):
        f = tmp_path / f"s{i}.fits"
        f.write_text("x")
        files.append(str(f))
    for f in files:
        window._add_recent_file(f)
    assert len(window._recent_files()) == 5
    app_settings().setValue(S.MAX_RECENT_KEY, 2)
    window._add_recent_file(files[0])
    assert window._recent_files() == [files[0], files[4]]


def _dirty(window, monkeypatch):
    history = window.active_history
    history.add_state(history.current_state)
    asked = []
    monkeypatch.setattr(window, "_show_custom_message",
                        lambda **k: asked.append(k["header"]) or QMessageBox.StandardButton.Cancel)
    return history, asked


def test_the_question_before_closing_a_changed_session(loaded, monkeypatch):
    history, asked = _dirty(loaded, monkeypatch)
    loaded._on_close_session_requested(history)
    assert len(asked) == 1 and history in loaded.session_histories  # Cancel


def test_closing_without_the_question(loaded, monkeypatch):
    history, asked = _dirty(loaded, monkeypatch)
    app_settings().setValue(S.CONFIRM_CLOSE_KEY, False)
    loaded._on_close_session_requested(history)
    assert asked == [] and history not in loaded.session_histories


def test_the_file_dialogs_start_in_the_fixed_folder(window, monkeypatch, tmp_path):
    seen = []
    monkeypatch.setattr(QFileDialog, "getOpenFileNames",
                        staticmethod(lambda parent, title, folder, *a: seen.append(folder) or ([], "")))
    app_settings().setValue(S.LAST_DIR_KEY, str(tmp_path / "last"))
    app_settings().setValue(S.DIR_MODE_KEY, "fixed")
    app_settings().setValue(S.FIXED_DIR_KEY, str(tmp_path))
    window._on_open_spectrum()
    assert seen == [str(tmp_path)]


def test_the_initial_x_unit(qtbot, monkeypatch):
    app_settings().setValue(S.INITIAL_X_UNIT_KEY, "micron")
    win = MainWindow(None, None)
    qtbot.addWidget(win)
    assert win.x_unit_combo.currentText() == "micron"


def test_decimals_of_the_data_inspector():
    assert format_value("x", 600.123456789) == "600.12346"
    app_settings().setValue(S.DECIMALS_KEY, 2)
    assert format_value("x", 600.123456789) == "600.12"
    assert format_value("y", 1.23456789) == "1.2346"  # the other columns do not follow


def _plotted_lines(window):
    window.plot_viewer.plot_spectrum(window.active_history.current_state)
    return {line.get_label(): line for line in window.plot_viewer.canvas.axes.get_lines()}


def test_line_width(loaded):
    base = _plotted_lines(loaded)["Spectrum"].get_linewidth()
    app_settings().setValue(S.LINE_WIDTH_KEY, 2.0)
    assert _plotted_lines(loaded)["Spectrum"].get_linewidth() == pytest.approx(2 * base)


def test_the_downsample_threshold(loaded):
    assert "Spectrum" in _plotted_lines(loaded)
    app_settings().setValue(S.DOWNSAMPLE_KEY, 10_000)
    assert "Spectrum (decimated)" in _plotted_lines(loaded)


# --- The check for updates at the start ---

def test_no_check_when_it_is_off(window):
    assert su.start_update_check(window, fake_urlopen({})) is None


@pytest.fixture
def pip_install(monkeypatch):
    """The tests run from a git checkout, which the check at start skips."""
    monkeypatch.setattr(su, "install_kind", lambda *a, **k: "pip")


def test_the_check_at_start_skips_development_installs(window, monkeypatch):
    app_settings().setValue(S.CHECK_UPDATES_KEY, True)
    monkeypatch.setattr(su, "install_kind", lambda *a, **k: "editable")
    uo = fake_urlopen({})  # any call would raise KeyError
    assert su.start_update_check(window, uo) is None
    assert uo.calls == [] and getattr(window, "_update_dialog", None) is None


def test_the_check_at_start_opens_the_dialog_for_a_newer_release(qtbot, window, pip_install):
    app_settings().setValue(S.CHECK_UPDATES_KEY, True)
    uo = fake_urlopen({su.RELEASES_API: json.dumps(RELEASE_JSON).encode()})
    worker = su.start_update_check(window, uo)
    assert worker is not None
    qtbot.waitUntil(lambda: getattr(window, "_update_dialog", None) is not None, timeout=5000)
    assert "v2.9.0" in window._update_dialog.status.text()
    assert app_settings().value(S.LAST_CHECK_KEY)  # recorded: a daily check is not due again
    assert su.start_update_check(window, uo) is None


def test_the_check_at_start_is_silent_without_network(qtbot, window, pip_install):
    app_settings().setValue(S.CHECK_UPDATES_KEY, True)
    worker = su.start_update_check(window, fake_urlopen({su.RELEASES_API: OSError("offline")}))
    qtbot.waitUntil(worker.isFinished, timeout=5000)
    qtbot.wait(50)
    assert getattr(window, "_update_dialog", None) is None
    assert not app_settings().value(S.LAST_CHECK_KEY)  # tried again at the next start


# --- Scripts, MCP and the log do not read the settings ---

def test_the_package_outside_the_gui_does_not_import_the_settings():
    root = Path(S.__file__).parents[1]
    for sub in ("core", "recipes", "fitting", "io", "session", "mcp_server"):
        for path in (root / sub).rglob("*.py"):
            assert "gui.settings" not in path.read_text(), path


def test_the_plot_arrays_do_not_depend_on_the_settings(loaded):
    spec = loaded.active_history.current_state.spec
    before = np.array(spec.y.value)
    app_settings().setValue(S.DOWNSAMPLE_KEY, 10_000)
    app_settings().setValue(S.LINE_WIDTH_KEY, 3.0)
    loaded.plot_viewer.plot_spectrum(None)
    np.testing.assert_array_equal(loaded.active_history.current_state.spec.y.value, before)
