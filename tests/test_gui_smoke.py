"""Smoke tests of the GUI (offscreen, pytest-qt).

They check that the main paths of the interface run without errors, not
what is drawn: the main window starts, a session is opened from a file,
every recipe dialog of the menus can be built, recipes run through the GUI
(dialog parameters as strings, worker thread, result processing, history),
errors are reported to the user without closing the window, undo and redo
move through the history, the tool windows open and a session can be
closed.

Message boxes are replaced by stubs that record what would have been shown
(modal dialogs would block an offscreen test).
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest
from PySide6.QtWidgets import QMessageBox

from astrocook.recipes.absorbers import ABSORBERS_RECIPES_SCHEMAS
from astrocook.recipes.continuum import CONTINUUM_RECIPES_SCHEMAS
from astrocook.recipes.edit import EDIT_RECIPES_SCHEMAS
from astrocook.recipes.file import FILE_RECIPES_SCHEMAS
from astrocook.recipes.flux import FLUX_RECIPES_SCHEMAS

# Small real archive (J0103-1305, 2156 pixels, 764 components), quick to load
SESSION_FILE = Path(__file__).parent / "data" / "J0103-1305_metallines_lya_new_3.acs2"

SCHEMAS = {"absorbers": ABSORBERS_RECIPES_SCHEMAS,
           "continuum": CONTINUUM_RECIPES_SCHEMAS,
           "edit": EDIT_RECIPES_SCHEMAS, "file": FILE_RECIPES_SCHEMAS,
           "flux": FLUX_RECIPES_SCHEMAS}
# Recipes with a dialog in the menus (the hidden ones are called by the
# system inspector or by other recipes with their own parameters)
DIALOG_RECIPES = [(cat, name) for cat, schemas in SCHEMAS.items()
                  for name, schema in schemas.items()
                  if not schema.get("gui_hidden")]


@pytest.fixture
def messages(monkeypatch):
    """Record the message boxes instead of showing them."""
    shown = []

    def record(kind):
        def stub(*args, **kwargs):
            shown.append((kind, args[1] if len(args) > 1 else kwargs.get("title")))
            return QMessageBox.StandardButton.Ok
        return staticmethod(stub)

    for kind in ("warning", "critical", "information", "question"):
        monkeypatch.setattr(QMessageBox, kind, record(kind))
    return shown


@pytest.fixture
def window(qtbot, monkeypatch, tmp_path, messages):
    # Keep the settings of the test run away from the user's ones
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    # Import only once qtbot has created the QApplication: the GUI selects the
    # QtAgg matplotlib backend at import time, which fails without one (F7).
    from astrocook.gui.main_window import MainWindowV2
    win = MainWindowV2(None, None)
    qtbot.addWidget(win)

    def custom_message(*args, **kwargs):
        messages.append(("custom", kwargs.get("header") or kwargs.get("title")))
        return QMessageBox.StandardButton.Ok

    monkeypatch.setattr(win, "_show_custom_message", custom_message)
    win.show()
    return win


@pytest.fixture
def loaded(window):
    window.open_session_from_path(str(SESSION_FILE))
    assert window.active_history is not None
    return window


def _run_recipe(qtbot, win, category, name, params):
    """Run a recipe as the dialog does and wait for the worker."""
    win._on_recipe_requested(category, name, params, {})
    qtbot.waitUntil(lambda: win.thread_pool.activeThreadCount() == 0, timeout=60000)
    qtbot.wait(50)  # let the result signals be processed


# --- Start and open -------------------------------------------------------

def test_main_window_starts_empty(window):
    assert window.windowTitle() == "Astrocook"
    assert window.active_history is None
    assert not window.doublets_auto_action.isEnabled()
    assert not window.undo_action.isEnabled()


def test_debug_plots_of_the_core_reach_the_gui(window, monkeypatch):
    # The core requests debug plots through a hook, not by importing the GUI
    # (F6); the window sends them to the main thread through its signal
    import matplotlib.pyplot as plt
    from astrocook.core import debug_hooks
    from astrocook.gui.debug_utils import GLOBAL_PLOTTER
    received = []
    GLOBAL_PLOTTER.plot_requested.connect(received.append)
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: None)
    plot = {"title": "CIV z=3.0000 (R2=0.990)", "v_compare": [0.0, 1.0],
            "Y_data": [0.1, 0.2], "Y_model": [0.1, 0.2]}
    try:
        debug_hooks.request_debug_plot(plot)
    finally:
        GLOBAL_PLOTTER.plot_requested.disconnect(received.append)
        debug_hooks.set_debug_plot_handler(None)
        plt.close("all")
    assert received == [plot]


def test_main_window_has_icon(window):
    # QIcon of a missing file is silently null (F34)
    assert not window.windowIcon().isNull()


def test_open_session(loaded, messages):
    history = loaded.active_history
    assert history.display_name == SESSION_FILE.stem
    assert SESSION_FILE.stem in loaded.windowTitle()
    assert len(loaded.session_histories) == 1
    assert len(history.current_state.systs.components) == 764
    assert loaded.doublets_auto_action.isEnabled()
    assert messages == []


def test_open_missing_or_invalid_file(window, messages, tmp_path):
    window.open_session_from_path(str(tmp_path / "missing.acs2"))
    invalid = tmp_path / "invalid.acs2"
    invalid.write_text("not an archive")
    window.open_session_from_path(str(invalid))
    assert window.session_histories == []
    assert [kind for kind, _ in messages] in (["warning"], ["critical"])


# --- Recipe dialogs -------------------------------------------------------

@pytest.mark.parametrize("category, name", DIALOG_RECIPES,
                         ids=[f"{c}.{n}" for c, n in DIALOG_RECIPES])
def test_recipe_dialog_builds(loaded, category, name):
    from astrocook.gui.recipe_dialog import RecipeDialog
    dialog = RecipeDialog(category, name, loaded.active_history.current_state, loaded)
    dialog.show()
    visible = [p["name"] for p in SCHEMAS[category][name]["params"]
               if not p.get("gui_hidden")]
    assert set(visible) <= set(dialog.input_widgets)
    dialog.close()


def test_launch_recipe_dialog_from_menu(loaded):
    loaded._launch_recipe_dialog("flux", "smooth")
    assert loaded.active_recipe_dialog is not None
    assert loaded.active_recipe_dialog.isVisible()
    loaded.active_recipe_dialog.close()


# --- Recipes through the GUI ----------------------------------------------

def test_set_properties(qtbot, loaded, messages):
    history = loaded.active_history
    _run_recipe(qtbot, loaded, "edit", "set_properties",
                {"z_em": "4.072", "resol": "45000"})
    assert history.current_index == 1
    assert history.current_state.spec._data.z_em == 4.072
    assert loaded.undo_action.isEnabled()
    assert messages == []


def test_recipe_requiring_z_em_asks_for_it(qtbot, loaded, messages):
    # The archive has no z_em: the GUI asks to set it instead of running
    history = loaded.active_history
    _run_recipe(qtbot, loaded, "continuum", "estimate_auto", {})
    assert history.current_index == 0
    assert messages == [("custom", "Emission Redshift (z_em) is required.")]


@pytest.fixture
def ready(qtbot, loaded):
    """Loaded session with z_em and R set (state 1 of the history)."""
    _run_recipe(qtbot, loaded, "edit", "set_properties",
                {"z_em": "4.072", "resol": "45000"})
    assert loaded.active_history.current_index == 1
    return loaded


@pytest.mark.parametrize("category, name, params", [
    ("continuum", "estimate_auto", {}),
    ("flux", "rebin", {"dx": "5.0", "xunit": "km/s"}),
    ("flux", "smooth", {"sigma_kms": "50.0"}),
])
def test_recipe_adds_a_state(qtbot, ready, messages, category, name, params):
    _run_recipe(qtbot, ready, category, name, params)
    assert ready.active_history.current_index == 2
    assert messages == []


def test_clean_negligible_after_continuum(qtbot, ready, messages):
    # Without a continuum the components cannot be deleted (F25)
    history = ready.active_history
    n_before = len(history.current_state.systs.components)
    _run_recipe(qtbot, ready, "continuum", "estimate_auto", {})
    _run_recipe(qtbot, ready, "absorbers", "clean_negligible", {})
    assert history.current_index == 3
    assert len(history.current_state.systs.components) < n_before
    assert messages == []


def test_recipe_returning_a_new_session(qtbot, ready, messages):
    # extract_preset opens the extracted region as a new session
    loaded = ready
    _run_recipe(qtbot, loaded, "edit", "extract_preset", {"region": "red_side"})
    assert len(loaded.session_histories) == 2
    x = loaded.session_histories[1].current_state.spec.x.value
    assert x.min() > 121.567 * (1 + 4.072) - 0.1
    assert messages == []


@pytest.mark.parametrize("name, params, title", [
    ("smooth", {"sigma_kms": "abc"}, "Recipe Aborted"),  # recipe returns 0
    ("rebin", {"dx": "abc"}, "Recipe Aborted"),   # returns 0 since the F23 fix
])
def test_recipe_errors_are_reported(qtbot, loaded, messages, name, params, title):
    history = loaded.active_history
    _run_recipe(qtbot, loaded, "flux", name, params)
    assert history.current_index == 0
    assert ("custom", title) in messages
    assert loaded.isVisible()


def test_undo_redo(qtbot, loaded):
    history = loaded.active_history
    _run_recipe(qtbot, loaded, "flux", "smooth", {"sigma_kms": "50.0"})
    smoothed = history.current_state
    loaded._undo_last_action()
    assert history.current_index == 0
    assert loaded.redo_action.isEnabled()
    loaded._redo_last_action()
    assert history.current_index == 1
    assert history.current_state is smoothed


# --- Tool windows, close --------------------------------------------------

def test_tool_windows_open(loaded):
    loaded._on_view_system_inspector()
    assert loaded.system_inspector.isVisible()
    loaded._on_view_log()
    assert loaded.log_scripter_dialog.isVisible()
    loaded._on_view_script_console()
    assert loaded.script_dialog.isVisible()


@pytest.mark.parametrize("installed", [True, False])
def test_inspector_bayesian_entries(loaded, monkeypatch, installed):
    # emcee and corner are optional (extra "bayes"): without them the menu
    # entries are disabled and their tooltip names the extra
    from PySide6.QtWidgets import QMenu
    from astrocook.fitting import bayesian_fitter
    from astrocook.gui import system_inspector
    if not installed:
        monkeypatch.setattr(bayesian_fitter, "emcee", None)
        monkeypatch.setattr(system_inspector, "_corner_installed", lambda: False)
    menus = []

    class RecordedMenu(QMenu):
        def exec(self, *args):  # record the menu instead of showing it
            menus.append(self)

    monkeypatch.setattr(system_inspector, "QMenu", RecordedMenu)

    loaded._on_view_system_inspector()
    inspector = loaded.system_inspector
    pos = inspector.table_view.visualRect(inspector.proxy_model.index(0, 0)).center()
    inspector._on_context_menu(pos)
    actions = {a.text(): a for a in menus[0].actions()}
    mcmc = next(a for t, a in actions.items() if t.startswith("Refit Bayesian (MCMC)"))
    corner = actions["Show Corner Plot"]
    assert mcmc.isEnabled() == installed
    assert not corner.isEnabled()  # no Bayesian results in this session
    for action in (mcmc, corner):
        assert ('astrocook[bayes]' in action.toolTip()) == (not installed)


def test_gui_starts_without_corner():
    # corner is optional (extra "bayes"): the GUI must import without it.
    # In a subprocess, where importing corner fails.
    code = ("import sys; sys.modules['corner'] = None\n"
            "from PySide6.QtWidgets import QApplication\n"
            "app = QApplication([])\n"
            "import astrocook.gui.main_window\n")
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=120,
                            env={**os.environ, "QT_QPA_PLATFORM": "offscreen"})
    assert result.returncode == 0, result.stderr[-2000:]


def test_close_session(loaded):
    loaded._close_session(loaded.active_history)
    assert loaded.session_histories == []
    assert loaded.active_history is None


# --- Session logs ---------------------------------------------------------
# Characterization for phase 3 (removal of the V1 GUILog): the log of every
# history is a HistoryLogV2, shared by all its states; a branch gets a copy.

def _recipes(log):
    return [e.recipe_name for e in log.entries]


def test_history_log_is_v2(ready):
    from astrocook.core.structures import HistoryLogV2
    history = ready.active_history
    assert isinstance(history.log_manager, HistoryLogV2)
    assert _recipes(history.log_manager) == ["set_properties"]
    assert all(s.log_manager is history.log_manager for s in history.states)


def test_branch_copies_the_log(qtbot, ready):
    from astrocook.core.structures import HistoryLogV2
    original = ready.active_history
    _run_recipe(qtbot, ready, "edit", "extract_preset", {"region": "red_side"})
    branch = ready.session_histories[1]
    assert isinstance(branch.log_manager, HistoryLogV2)
    assert branch.log_manager is not original.log_manager
    assert _recipes(branch.log_manager) == ["set_properties", "extract_preset"]
    assert all(s.log_manager is branch.log_manager for s in branch.states)
    assert len(original.states) == 2  # the parent session does not change


def test_branch_does_not_change_the_parent_log(qtbot, ready):
    # The parent session does not change, so neither does its log (F41)
    original = ready.active_history
    _run_recipe(qtbot, ready, "edit", "extract_preset", {"region": "red_side"})
    assert _recipes(original.log_manager) == ["set_properties"]
    assert original.log_manager.current_index == original.current_index - 1


def test_branch_after_undo_keeps_the_parent_redo(qtbot, ready):
    # Logging the branch on the parent also truncated its redo entries (F41)
    original = ready.active_history
    _run_recipe(qtbot, ready, "flux", "smooth", {"sigma_kms": "50.0"})
    ready._undo_last_action()
    _run_recipe(qtbot, ready, "edit", "extract_preset", {"region": "red_side"})
    assert _recipes(original.log_manager) == ["set_properties", "smooth"]
    assert original.can_redo() and original.log_manager.current_index == 0
    branch = ready.session_histories[1]
    assert _recipes(branch.log_manager) == ["set_properties", "extract_preset"]
    assert len(branch.states) == 3


def test_each_branch_logs_the_recipe(qtbot, ready):
    original = ready.active_history
    _run_recipe(qtbot, ready, "edit", "extract_preset",
                {"region": "lya_forest,red_side"})
    assert len(ready.session_histories) == 3
    assert _recipes(original.log_manager) == ["set_properties"]
    for branch in ready.session_histories[1:]:
        assert _recipes(branch.log_manager) == ["set_properties", "extract_preset"]
        assert branch.log_manager.current_index == branch.current_index - 1


def test_script_run_builds_a_v2_history(ready):
    from astrocook.core.structures import HistoryLogV2
    from astrocook.gui.qt_workers import ScriptWorker
    script = ("set_properties(z_em='4.072', resol='45000')\n"
              "# a comment\n"
              "smooth(sigma_kms='50.0')\n")
    worker = ScriptWorker(script, ready.active_history.states[0], ready)
    results, errors = [], []
    worker.signals.finished.connect(results.append)
    worker.signals.error.connect(errors.append)
    worker.run()  # synchronously, in this thread
    assert errors == []
    (history,) = results
    assert isinstance(history.log_manager, HistoryLogV2)
    assert _recipes(history.log_manager) == ["set_properties", "smooth"]
    assert len(history.states) == 3


def test_open_v1_archive(window, messages):
    # V1 .acs archives are read by io/v1_archive.py (phase 3)
    path = Path(__file__).parent / "data" / "v1_metals_constraints.acs"
    window.open_session_from_path(str(path))
    history = window.active_history
    assert history is not None and history.display_name == path.stem
    assert len(history.current_state.systs.components) == 10
    # The V1 links are shown by the System Inspector (F43)
    window._on_view_system_inspector()
    assert window.system_inspector.isVisible()
    assert messages == []
