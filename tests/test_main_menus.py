"""The menus of the main window: their entries and shortcuts, the actions
enabled with a session, the recipe dialog each entry opens, undo and
redo, the windows of the View menu, export and import, the documentation.

Characterization tests written before the menus were moved to a mixin
(phase 5, split of the GUI).
"""
import sys

import pytest
from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import QMessageBox

from astrocook.gui.main_window import MainWindowV2

from test_gui_smoke import SESSION_FILE, _run_recipe

SEP = "---"

MENUS = {
    "File": ["Open Session...", "Export to ASCII...", "Import from ASCII...", "Save Session...",
             SEP, "Close Session"],
    "View": ["Toggle Session Panel", "Toggle Plot Controls", "Show Toolbar", "Appearance", SEP,
             "View System Inspector", SEP,
             "View Session Log", "View Script Console", "View Identifications"],
    "Edit": ["Undo", "Redo", SEP, "Set Properties...", SEP, "Apply Expression...",
             "Mask by Expression...", "Smooth Column...", "Import Systems...",
             "Import Systems from ASCII...", "Delete Elements...", SEP, "Split Out Region...",
             "Extract Preset Region...", SEP, "Co-add Spectra..."],
    "Flux": ["Calculate Running StdDev...", SEP, "Smooth Spectrum...", "Rebin Spectrum...",
             "Equalize to Reference...", "Resample on Grid...", SEP, "Flux Calibrate..."],
    "Continuum": ["Auto-estimate Continuum...", SEP, "Find Absorbed Regions...",
                  "Fit Continuum to Mask...", SEP, "Fit Power-Law..."],
    "Absorbers": ["Auto-detect Doublets...", "Auto-fit Ly-alpha Forest...", SEP,
                  "Identify Absorption Lines...", "Refit All Systems...",
                  "Clean Negligible Components..."],
    "Help": ["Online Documentation", SEP, "About Astrocook..."],
}

#: The recipe dialog opened by each entry (category, recipe)
RECIPE_ENTRIES = {
    "Set Properties...": ("edit", "set_properties"),
    "Apply Expression...": ("edit", "apply_expression"),
    "Mask by Expression...": ("edit", "mask_expression"),
    "Smooth Column...": ("edit", "smooth_column"),
    "Import Systems...": ("edit", "import_systems"),
    "Delete Elements...": ("edit", "delete"),
    "Split Out Region...": ("edit", "split"),
    "Extract Preset Region...": ("edit", "extract_preset"),
    "Co-add Spectra...": ("edit", "coadd"),
    "Calculate Running StdDev...": ("flux", "calculate_running_std"),
    "Smooth Spectrum...": ("flux", "smooth"),
    "Rebin Spectrum...": ("flux", "rebin"),
    "Equalize to Reference...": ("flux", "equalize"),
    "Resample on Grid...": ("flux", "resample"),
    "Flux Calibrate...": ("flux", "calibrate_from_magnitudes"),
    "Auto-estimate Continuum...": ("continuum", "estimate_auto"),
    "Find Absorbed Regions...": ("continuum", "find_absorbed"),
    "Fit Continuum to Mask...": ("continuum", "fit_continuum"),
    "Fit Power-Law...": ("continuum", "fit_powerlaw"),
    "Auto-detect Doublets...": ("absorbers", "doublets_auto"),
    "Auto-fit Ly-alpha Forest...": ("absorbers", "lya_auto"),
    "Identify Absorption Lines...": ("absorbers", "identify_lines"),
    "Refit All Systems...": ("absorbers", "refit_all"),
    "Clean Negligible Components...": ("absorbers", "clean_negligible"),
    "Export to ASCII...": ("file", "export_ascii"),
}


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
def empty(qtbot, monkeypatch, tmp_path):
    return _new_window(qtbot, monkeypatch, tmp_path)


@pytest.fixture
def window(qtbot, monkeypatch, tmp_path):
    win = _new_window(qtbot, monkeypatch, tmp_path)
    win.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(50)
    return win


def _menus(window):
    """{menu title: [entry text, ...]} without the accelerators."""
    out = {}
    for menu_action in window.menuBar().actions():
        menu = menu_action.menu()
        out[menu.title().replace("&", "")] = [
            SEP if a.isSeparator() else a.text().replace("&", "") for a in menu.actions()]
    return out


def _actions(window):
    return {a.text().replace("&", ""): a
            for m in window.menuBar().actions() for a in m.menu().actions()
            if not a.isSeparator()}


def _module(window, method):
    return sys.modules[getattr(type(window), method).__module__]


# --- Entries -------------------------------------------------------------------------------

def test_menus_and_their_entries(empty):
    assert _menus(empty) == MENUS


def test_shortcuts(empty):
    actions = _actions(empty)
    assert actions["Open Session..."].shortcut().toString() == "Ctrl+O"
    assert actions["Save Session..."].shortcut().toString() == "Ctrl+S"
    assert actions["Close Session"].shortcut().toString() == "Ctrl+W"
    assert actions["Undo"].shortcut() == QKeySequence(QKeySequence.StandardKey.Undo)
    assert actions["Redo"].shortcut() == QKeySequence(QKeySequence.StandardKey.Redo)
    assert (actions["Online Documentation"].shortcut()
            == QKeySequence(QKeySequence.StandardKey.HelpContents))


def test_actions_enabled_with_a_session(empty, qtbot):
    # Also with no session: the toolbar and the appearance (phase 5-bis)
    always = {"Open Session...", "Online Documentation", "About Astrocook...",
              "Show Toolbar", "Appearance"}
    for text, action in _actions(empty).items():
        assert action.isEnabled() == (text in always), text
    empty.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(50)
    # Undo/redo follow the history, View Identifications needs identifications
    off = {"Undo", "Redo", "View Identifications"}
    for text, action in _actions(empty).items():
        assert action.isEnabled() == (text not in off), text
    # The two pipelines kept out of the menu are enabled too
    assert empty.forest_auto_action.isEnabled() and empty.metals_auto_action.isEnabled()


def test_recipe_entries_open_their_dialog(window, monkeypatch):
    opened = []
    monkeypatch.setattr(window, "_launch_recipe_dialog",
                        lambda cat, name, initial_params=None: opened.append(
                            (cat, name, initial_params)))
    actions = _actions(window)
    for text, (cat, name) in RECIPE_ENTRIES.items():
        actions[text].trigger()
        assert opened[-1][:2] == (cat, name), text
    actions["Co-add Spectra..."].trigger()
    assert opened[-1][2] == {"session_names": window.active_history.display_name}
    for action, (cat, name) in ((window.forest_auto_action, ("absorbers", "forest_auto")),
                                (window.metals_auto_action, ("absorbers", "metals_auto"))):
        action.trigger()
        assert opened[-1][:2] == (cat, name)


# --- Undo and redo ---------------------------------------------------------------------------

def test_undo_and_redo(qtbot, window):
    history = window.active_history
    actions = _actions(window)
    assert not actions["Undo"].isEnabled() and not actions["Redo"].isEnabled()
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072"})
    new_state = history.current_state
    assert actions["Undo"].isEnabled() and not actions["Redo"].isEnabled()
    actions["Undo"].trigger()
    assert history.current_state is not new_state and history.current_state.spec._data.z_em != 4.072
    assert not actions["Undo"].isEnabled() and actions["Redo"].isEnabled()
    actions["Undo"].trigger()  # nothing older: no change
    actions["Redo"].trigger()
    assert history.current_state is new_state
    assert actions["Undo"].isEnabled() and not actions["Redo"].isEnabled()


def test_undo_without_a_session(empty):
    empty._undo_last_action()
    empty._redo_last_action()
    assert not empty.undo_action.isEnabled() and not empty.redo_action.isEnabled()


# --- View --------------------------------------------------------------------------------------

def test_view_windows(qtbot, window):
    actions = _actions(window)
    actions["View System Inspector"].trigger()
    assert window.system_inspector.isVisible()
    actions["View Script Console"].trigger()
    assert window.script_dialog.isVisible()
    console = window.script_dialog
    actions["View Script Console"].trigger()
    assert window.script_dialog is console  # the same window
    actions["View Session Log"].trigger()
    assert window.log_scripter_dialog.isVisible()


def test_identifications(qtbot, window, monkeypatch):
    launched = []
    monkeypatch.setattr(window, "_launch_identification_viewer", launched.append)
    window._on_view_identifications()
    assert launched == [window.active_history]


def test_view_without_a_session(empty, caplog):
    empty._on_view_log()
    empty._on_view_identifications()
    assert "View Log called with no active history." in caplog.text
    assert "View Identifications called with no active history." in caplog.text


# --- File and Help -----------------------------------------------------------------------------

def test_import_ascii(window, monkeypatch, tmp_path):
    requests = []
    monkeypatch.setattr(window, "_on_recipe_requested",
                        lambda cat, name, params, extra: requests.append((cat, name, params)))
    dialog = _module(window, "_on_import_requested").QFileDialog
    monkeypatch.setattr(dialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (str(tmp_path / "cols.csv"), "")))
    _actions(window)["Import from ASCII..."].trigger()
    assert requests == [("file", "import_ascii", {"file_path": str(tmp_path / "cols.csv")})]
    # Cancelled at the warning: nothing
    monkeypatch.setattr(window, "_show_custom_message",
                        lambda *a, **k: QMessageBox.StandardButton.Cancel)
    _actions(window)["Import from ASCII..."].trigger()
    assert len(requests) == 1
    # No file chosen: nothing
    monkeypatch.setattr(window, "_show_custom_message",
                        lambda *a, **k: QMessageBox.StandardButton.Ok)
    monkeypatch.setattr(dialog, "getOpenFileName", staticmethod(lambda *a, **k: ("", "")))
    _actions(window)["Import from ASCII..."].trigger()
    assert len(requests) == 1


def test_import_systems_ascii(window, monkeypatch, tmp_path):
    requests = []
    monkeypatch.setattr(window, "_on_recipe_requested",
                        lambda cat, name, params, extra: requests.append((cat, name, params)))
    dialog = _module(window, "_import_systems_ascii").QFileDialog
    monkeypatch.setattr(dialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (str(tmp_path / "systs.csv"), "")))
    _actions(window)["Import Systems from ASCII..."].trigger()
    assert requests == [("edit", "import_systems_ascii",
                         {"file_path": str(tmp_path / "systs.csv"), "append": "True",
                          "refit": "False"})]


def test_documentation(empty, monkeypatch):
    opened = []
    services = _module(empty, "_open_docs").QDesktopServices
    monkeypatch.setattr(services, "openUrl", staticmethod(lambda url: opened.append(url.toString())))
    _actions(empty)["Online Documentation"].trigger()
    assert opened == ["https://das-oats.github.io/astrocook/"]


def test_placeholder_dialogs_are_gone(empty):
    # Placeholders not reachable from the menus, removed (findings register F77)
    for name in ("_launch_x_convert_dialog", "_launch_y_convert_dialog",
                 "_launch_rebin_dialog"):
        assert not hasattr(empty, name)
