"""The recipes and scripts of the main window: the requirements (z_em,
resolution) and the recipe resumed after Set Properties, the recipe
dialog, the checks before a recipe runs (export, resolution overwrite,
co-add with mixed resolutions), the progress dialog and the worker, the
result (log, warnings, cancellation), errors, Stop, the scripts and the
log scripter.

Characterization tests written before the recipes and scripts were moved
to mixins (phase 5, split of the GUI). The workers are not started: the
recipes run in ``test_gui_smoke.py``.
"""
import logging
import sys

import numpy as np
import pytest
from PySide6.QtWidgets import QMessageBox, QProgressDialog

from astrocook.core.structures import HistoryLogV2
from astrocook.gui.main_window import MainWindowV2
from astrocook.gui.recipe_progress import RecipeProgressDialog

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
def window(qtbot, monkeypatch, tmp_path):
    win = _new_window(qtbot, monkeypatch, tmp_path)
    win.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(50)
    return win


@pytest.fixture
def ready(qtbot, window):
    """The archive with z_em and the resolution."""
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072", "resol": "45000"})
    return window


@pytest.fixture
def workers(window, monkeypatch):
    """The workers started (not run)."""
    started = []
    monkeypatch.setattr(window.thread_pool, "start", started.append)
    return started


def _module(window, method):
    return sys.modules[getattr(type(window), method).__module__]


def _messages(window, monkeypatch, answer="Ok", checked=False):
    """The messages shown (title, header, text, kwargs), answered."""
    shown = []

    def show(*args, **kwargs):
        shown.append(kwargs)
        button = getattr(QMessageBox.StandardButton, answer)
        return (button, checked) if "checkbox_text" in kwargs else button
    monkeypatch.setattr(window, "_show_custom_message", show)
    return shown


def _dialogs(window, monkeypatch):
    opened = []
    monkeypatch.setattr(window, "_launch_recipe_dialog",
                        lambda cat, name, initial_params=None: opened.append(
                            (cat, name, initial_params)))
    return opened


# --- Requirements -------------------------------------------------------------------------

@pytest.mark.parametrize("answer", ["Yes", "Cancel"])
def test_recipe_needs_z_em(window, monkeypatch, answer):
    shown = _messages(window, monkeypatch, answer)
    opened = _dialogs(window, monkeypatch)
    assert not window._check_and_handle_requirements("continuum", "estimate_auto", {"a": 1},
                                                     mode="direct")
    assert shown[-1]["header"] == "Emission Redshift (z_em) is required."
    if answer == "Yes":
        assert opened == [("edit", "set_properties", None)]
        assert window._pending_recipe == {
            "category": "continuum", "name": "estimate_auto", "params": {"a": 1},
            "mode": "direct", "after": "set_properties"}
    else:
        assert opened == [] and window._pending_recipe is None


def test_recipe_needs_the_resolution(qtbot, window, monkeypatch):
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072"})
    shown = _messages(window, monkeypatch, "Cancel")
    assert not window._check_and_handle_requirements("absorbers", "lya_auto")
    assert shown[-1]["header"] == "Instrument Resolution is missing."
    assert window._check_and_handle_requirements("continuum", "estimate_auto")  # z_em only
    assert window._check_and_handle_requirements("flux", "smooth")  # nothing


def test_recipe_resumed_after_set_properties(ready, monkeypatch, qtbot):
    window = ready
    opened = _dialogs(window, monkeypatch)
    window._set_pending_action("continuum", "estimate_auto", {"a": 1}, "dialog")
    state = window.active_history.current_state
    window._process_recipe_result((state.with_new_spectrum(state.spec), "set_properties",
                                   {"z_em": "4.072"}))
    qtbot.waitUntil(lambda: bool(opened))
    assert opened == [("continuum", "estimate_auto", {"a": 1})]
    assert window._pending_recipe is None
    requested = []
    monkeypatch.setattr(window, "_on_recipe_requested",
                        lambda cat, name, params, extra: requested.append((cat, name, params)))
    window._set_pending_action("continuum", "estimate_auto", {"b": 2}, "direct")
    window._process_recipe_result((state.with_new_spectrum(state.spec), "set_properties",
                                   {"z_em": "4.072"}))
    qtbot.waitUntil(lambda: bool(requested))
    assert requested == [("continuum", "estimate_auto", {"b": 2})]


# --- The continuum (phase 6-bis, PR 5b) ---------------------------------------------------

def _requested(window, monkeypatch):
    requested = []
    monkeypatch.setattr(window, "_on_recipe_requested",
                        lambda cat, name, params, extra: requested.append((cat, name, params)))
    return requested


def test_recipe_needs_the_continuum_estimate(ready, monkeypatch):
    # The archive has z_em and the resolution, and no continuum
    window = ready
    assert window.active_history.current_state.spec.norm is None
    shown = _messages(window, monkeypatch, "Yes")
    requested = _requested(window, monkeypatch)
    assert not window._check_and_handle_requirements("absorbers", "refit_all", {"a": 1},
                                                     mode="direct")
    assert shown[-1]["header"] == "The spectrum has no continuum."
    assert shown[-1]["text"].startswith("Refit all systems needs a continuum.")
    assert shown[-1]["button_texts"][QMessageBox.StandardButton.Open] == \
        "Open Auto-estimate Continuum…"
    assert requested == [("continuum", "estimate_auto", {})]
    assert window._pending_recipe == {"category": "absorbers", "name": "refit_all",
                                      "params": {"a": 1}, "mode": "direct",
                                      "after": "estimate_auto"}


def test_recipe_needs_the_continuum_open_or_cancel(ready, monkeypatch):
    window = ready
    _messages(window, monkeypatch, "Open")
    opened = _dialogs(window, monkeypatch)
    assert not window._check_and_handle_requirements("absorbers", "refit_all")
    assert opened == [("continuum", "estimate_auto", None)]
    assert window._pending_recipe["after"] == "estimate_auto"
    window._pending_recipe = None
    _messages(window, monkeypatch, "Cancel")
    requested = _requested(window, monkeypatch)
    assert not window._check_and_handle_requirements("absorbers", "refit_all")
    assert requested == [] and window._pending_recipe is None
    # The other gaps come first, in the order z_em, resolution, continuum
    assert window._check_and_handle_requirements("flux", "smooth")


def test_recipe_resumed_after_the_continuum(ready, monkeypatch, qtbot):
    window = ready
    requested = _requested(window, monkeypatch)
    window._set_pending_action("absorbers", "refit_all", {"b": 2}, "direct",
                               after="estimate_auto")
    state = window.active_history.current_state
    # Another recipe does not resume it
    window._process_recipe_result((state.with_new_spectrum(state.spec), "smooth", {}))
    qtbot.wait(20)
    assert requested == [] and window._pending_recipe is not None
    window._process_recipe_result((state.with_new_spectrum(state.spec), "estimate_auto", {}))
    qtbot.waitUntil(lambda: bool(requested))
    assert requested == [("absorbers", "refit_all", {"b": 2})]
    assert window._pending_recipe is None


def test_a_failed_estimate_drops_the_pending_recipe(ready):
    window = ready
    window._set_pending_action("absorbers", "refit_all", {}, "direct", after="estimate_auto")
    window._on_recipe_error(("Recipe Failed", "no unabsorbed pixels", None))
    assert window._pending_recipe is None


def test_cancelled_estimate_drops_the_pending_recipe(ready):
    window = ready
    window._set_pending_action("absorbers", "refit_all", {}, "dialog", after="estimate_auto")
    window._launch_recipe_dialog("continuum", "estimate_auto")
    window.active_recipe_dialog.reject()
    assert window._pending_recipe is None


# --- The recipe dialog ----------------------------------------------------------------------

def test_recipe_dialog(ready, monkeypatch):
    window = ready
    viewer = window.plot_viewer
    window._launch_recipe_dialog("flux", "smooth")
    dialog = window.active_recipe_dialog
    assert dialog.recipe_name == "smooth" and dialog.isVisible()
    # No tool is on while the dialog is open (before, pan(False) and
    # zoom(False) toggled them and left zoom on: findings register F81)
    assert viewer.toolbar.mode == ""
    dialog.reject()
    # A tool that is on is switched off
    viewer.toolbar.pan()
    assert viewer.toolbar.mode == "pan/zoom"
    window._launch_recipe_dialog("flux", "smooth")
    dialog = window.active_recipe_dialog
    assert viewer.toolbar.mode == ""
    # A second dialog is not opened while one is
    window._launch_recipe_dialog("flux", "rebin")
    assert window.active_recipe_dialog is dialog
    dialog.reject()
    assert window.active_recipe_dialog is None


def test_the_dialogs_keep_a_reasonable_width(ready):
    # The descriptions are sentences since phase 6-bis (PR 6a): their labels
    # wrap, and the dialogs did grow wider than the window without it
    from astrocook.recipes.registry import all_recipes
    window = ready
    for spec in all_recipes():
        if spec.gui_hidden:
            continue
        window._launch_recipe_dialog(spec.category, spec.name)
        dialog = window.active_recipe_dialog
        if dialog is None:
            continue
        assert dialog.sizeHint().width() < 900, spec.name
        dialog.reject()


def test_set_properties_dialog_shows_the_current_values(ready):
    window = ready
    window._launch_recipe_dialog("edit", "set_properties", {"object": "QSO"})
    widgets = window.active_recipe_dialog.input_widgets
    assert widgets["name"].text() == window.active_history.display_name
    assert widgets["z_em"].text() == "4.07200"
    assert widgets["resol"].text() == "45000"
    assert widgets["object"].text() == "QSO"  # given: it wins
    window.active_recipe_dialog.reject()


def test_cancelled_set_properties_drops_the_pending_recipe(window, monkeypatch):
    window._set_pending_action("continuum", "estimate_auto", {}, "dialog")
    window._hidden_stitch_sessions = "b, c"
    window._launch_recipe_dialog("edit", "set_properties")
    window.active_recipe_dialog.reject()
    assert window._pending_recipe is None
    assert not hasattr(window, "_hidden_stitch_sessions")


def test_recipe_dialog_without_a_session(qtbot, monkeypatch, tmp_path):
    empty = _new_window(qtbot, monkeypatch, tmp_path)
    warned = []
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda *a: warned.append(a[1])))
    empty._launch_recipe_dialog("flux", "smooth")
    assert warned == ["No Session"] and empty.active_recipe_dialog is None


def test_split_from_the_view(window, monkeypatch):
    opened = _dialogs(window, monkeypatch)
    ax = window.plot_viewer.canvas.axes
    ax.set_xlim(500.0, 520.5)
    window.launch_split_from_current_view()
    assert opened == [("edit", "split", {"expression": "(λ > 500.0000) & (λ < 520.5000)"})]


# --- Before a recipe runs --------------------------------------------------------------------

def test_short_and_long_recipes(ready, workers, monkeypatch):
    window = ready
    # The progress dialogs, not the prerequisites (the archive has no
    # continuum: the GUI would ask for it first, PR 5b)
    monkeypatch.setattr(window, "_check_and_handle_requirements", lambda *a, **k: True)
    window._on_recipe_requested("flux", "smooth", {"sigma_kms": "10"}, {})
    assert isinstance(window.progress_dialog, QProgressDialog)
    assert window.qt_log_handler is None
    worker = workers[-1]
    assert (worker.recipe_name, worker.params) == ("smooth", {"sigma_kms": "10"})
    assert window._last_attempted_recipe == {"category": "flux", "recipe_name": "smooth"}
    for name, params, indeterminate in (("lya_auto", {}, False),
                                        ("add_component", {}, True),
                                        ("refit_all", {}, False),
                                        ("refit_all", {"selected_uuids": "x"}, True)):
        window._on_recipe_requested("absorbers", name, params, {})
        dialog = window.progress_dialog
        assert isinstance(dialog, RecipeProgressDialog), name
        assert (dialog.progress_bar.maximum() == 0) == indeterminate, (name, params)
        assert window.qt_log_handler in logging.getLogger().handlers
        window._cleanup_progress_ui()
    window._on_recipe_requested("edit", "import_systems", {"refit": "True"}, {})
    assert isinstance(window.progress_dialog, RecipeProgressDialog)
    window._cleanup_progress_ui()


def _save_as(window, monkeypatch, tmp_path, chosen=None):
    """The file dialog answered with the suggested name (or ``chosen``)."""
    dialog = _module(window, "_on_recipe_requested").QFileDialog
    suggested = []

    def save(parent, title, name, filters):
        suggested.append(name)
        return str(tmp_path / (chosen or name)), ""
    monkeypatch.setattr(dialog, "getSaveFileName", staticmethod(save))
    return suggested


def test_export_ascii(qtbot, ready, monkeypatch, tmp_path):
    window = ready
    _run_recipe(qtbot, window, "continuum", "estimate_auto", {})
    shown = _messages(window, monkeypatch, "Yes")
    suggested = _save_as(window, monkeypatch, tmp_path)
    name = window.active_history.current_state.name
    # The default targets: the spectrum (cont) and the systems
    window._on_recipe_requested("file", "export_ascii", {"targets": "F, cont, systems"}, {})
    assert suggested == [name]
    assert (tmp_path / f"{name}_spec.csv").exists() and (tmp_path / f"{name}_systems.csv").exists()
    assert shown[-1]["header"] == "Data exported successfully"
    assert shown[-1]["text"] == (f"Saved in {tmp_path}:\n\n• {name}_spec.csv\n"
                                 f"• {name}_systems.csv")
    # The suggested name (exported in other folders: no conflicts)
    for targets, expected, folder in (("y, dy", f"{name}_spec.csv", "a"),
                                      ("systems", f"{name}_systems.csv", "b")):
        (tmp_path / folder).mkdir()
        suggested = _save_as(window, monkeypatch, tmp_path / folder)
        window._on_recipe_requested("file", "export_ascii", {"targets": targets}, {})
        assert suggested == [expected]
        assert (tmp_path / folder / expected).exists()


@pytest.mark.parametrize("answer", ["Yes", "No"])
def test_export_over_existing_files(ready, monkeypatch, tmp_path, answer):
    # Asked before overwriting (F83: the question raised UnboundLocalError,
    # QMessageBox was imported again further down in the same method)
    window = ready
    name = window.active_history.current_state.name
    (tmp_path / f"{name}_spec.csv").write_text("old")
    shown = _messages(window, monkeypatch, answer)
    _save_as(window, monkeypatch, tmp_path, f"{name}_spec.csv")
    window._on_recipe_requested("file", "export_ascii", {"targets": "y, dy"}, {})
    assert shown[0]["header"] == "Files already exist"
    assert shown[0]["text"].startswith(f"The following files will be overwritten:\n\n• {name}_spec.csv")
    assert ((tmp_path / f"{name}_spec.csv").read_text() == "old") == (answer == "No")


@pytest.mark.parametrize("targets", ["F, dF", "λ, F", "F, cont, systems", "spec"])
def test_export_with_the_names_of_the_gui(ready, monkeypatch, tmp_path, targets):
    # Targets named as in the GUI (λ, F, dF), or the whole spectrum: the
    # spectrum file is expected, and an existing one is not overwritten
    # without asking (F82)
    window = ready
    name = window.active_history.current_state.name
    (tmp_path / f"{name}_spec.csv").write_text("old")
    shown = _messages(window, monkeypatch, "No")
    suggested = _save_as(window, monkeypatch, tmp_path)
    window._on_recipe_requested("file", "export_ascii", {"targets": targets}, {})
    expected = name if "systems" in targets else f"{name}_spec.csv"
    assert suggested == [expected]
    assert [m["header"] for m in shown] == ["Files already exist"]
    assert (tmp_path / f"{name}_spec.csv").read_text() == "old"


def test_variable_resolution_overwrite_asked(ready, monkeypatch, workers):
    window = ready
    history = window.active_history
    spec = history.current_state.spec
    n = len(spec.x.value)
    spec = spec.update_column("resol", np.linspace(40000, 50000, n))
    history.add_state(history.current_state.with_new_spectrum(spec))
    shown = _messages(window, monkeypatch, "Yes")
    params = {"resol": "30000"}
    window._on_recipe_requested("edit", "set_properties", params, {})
    assert shown[0]["header"] == "The resolution varies with wavelength."
    assert params["overwrite_resol_col"] is True
    window._cleanup_progress_ui()
    params = {"resol": "_current_"}
    window._on_recipe_requested("edit", "set_properties", params, {})
    assert "overwrite_resol_col" not in params


@pytest.mark.parametrize("answer", ["Yes", "No"])
def test_coadd_with_mixed_resolutions(qtbot, ready, monkeypatch, workers, answer):
    window = ready
    window.open_session_from_path(str(SESSION_FILE))  # no resolution
    qtbot.wait(50)
    shown = _messages(window, monkeypatch, answer)
    names = ", ".join(h.display_name for h in window.session_histories)
    window._on_recipe_requested("edit", "coadd", {"session_names": names}, {})
    assert shown[0]["header"] == "Missing Resolution Data"
    assert "Only 1 out of 2" in shown[0]["text"]
    assert len(workers) == (1 if answer == "Yes" else 0)


def test_equalize_and_stitch_gets_the_hidden_sessions(ready, workers):
    window = ready
    window._hidden_stitch_sessions = "b, c"
    params = {}
    window._on_recipe_requested("edit", "equalize_and_stitch", params, {})
    assert params["other_sessions"] == "b, c" and not hasattr(window, "_hidden_stitch_sessions")


# --- Results --------------------------------------------------------------------------------

def test_result_logged_once(ready):
    window = ready
    history = window.active_history
    log = history.log_manager
    n = len(log.entries)
    state = history.current_state
    window._process_recipe_result((state.with_new_spectrum(state.spec), "smooth",
                                   {"sigma_kms": "10", "alias_map": {}}))
    assert len(log.entries) == n + 1
    entry = log.entries[-1]
    assert entry.recipe_name == "smooth" and entry.params == {"sigma_kms": 10.0}  # typed (6-bis)


def test_branching_result_logged_in_the_branch(ready):
    window = ready
    parent = window.active_history
    n = len(parent.log_manager.entries)
    state = parent.current_state
    window._process_recipe_result(([state.with_new_spectrum(state.spec),
                                    state.with_new_spectrum(state.spec)], "split",
                                   {"expression": "x > 0"}))
    # Two branches, both of the parent; logged in each branch, not in the
    # parent (F41)
    assert len(window.session_histories) == 3
    assert len(parent.log_manager.entries) == n
    base = parent.states[0].name
    assert [h.display_name for h in window.session_histories[1:]] == [f"{base}-1", f"{base}-2"]
    for branch in window.session_histories[1:]:
        assert len(branch.log_manager.entries) == n + 1
        assert branch.log_manager.entries[-1].recipe_name == "split"


def test_stopped_recipe(ready):
    # The worker says it was stopped (it was a _stop_flag written in the
    # result, phase 6-ter PR 7-quater): the history does not change
    window = ready
    history = window.active_history
    n = len(history.states)
    window._pending_recipe = ("absorbers", "lya_auto", {}, {})
    window._on_recipe_stopped("lya_auto")
    assert len(history.states) == n and window._pending_recipe is None


def test_refit_warning(ready, monkeypatch):
    window = ready
    shown = _messages(window, monkeypatch, "Ok", checked=True)
    assert window._suppress_refit_warning  # off by default
    window._suppress_refit_warning = False
    state = window.active_history.current_state
    window._process_recipe_result((state.with_new_spectrum(state.spec), "update_component", {}))
    assert shown[-1]["header"] == "Component parameters updated."
    assert window._suppress_refit_warning  # "Don't show this again"
    window._process_recipe_result((state.with_new_spectrum(state.spec), "update_component", {}))
    assert len(shown) == 1


def test_captured_warnings_shown(ready, monkeypatch):
    window = ready
    shown = _messages(window, monkeypatch)
    window._captured_warnings = ["b at the bound", "b at the bound", "Recipe Aborted"]
    state = window.active_history.current_state
    window._process_recipe_result((state.with_new_spectrum(state.spec), "smooth", {}))
    assert shown[-1]["header"] == "Process completed with warnings"
    assert shown[-1]["text"] == ("The recipe finished, but generated 1 warning(s):"
                                 "<ul><li>b at the bound</li></ul>")


# --- Errors and Stop ------------------------------------------------------------------------

def test_internal_error(ready, monkeypatch):
    window = ready
    shown = _messages(window, monkeypatch)
    window._captured_warnings = []
    window._on_recipe_error(("Recipe Error: smooth", "boom", "Traceback ..."))
    assert shown[-1]["title"] == "Internal Error"
    assert shown[-1]["detailed_text"].startswith("Traceback:\n")


@pytest.mark.parametrize("answer", ["Retry", "Cancel"])
def test_user_error(ready, monkeypatch, answer):
    window = ready
    shown = _messages(window, monkeypatch, answer)
    opened = _dialogs(window, monkeypatch)
    window._last_attempted_recipe = {"category": "flux", "recipe_name": "smooth"}
    window._on_recipe_error(("Recipe Failed: smooth", "bad sigma", None))
    assert shown[-1]["title"] == "Recipe Failed" and shown[-1]["detailed_text"] is None
    assert opened == ([("flux", "smooth", None)] if answer == "Retry" else [])


def test_stop(ready, workers, monkeypatch):
    window = ready
    monkeypatch.setattr(window, "_check_and_handle_requirements", lambda *a, **k: True)
    window._on_recipe_requested("absorbers", "lya_auto", {}, {})
    worker = workers[-1]
    stopped = []
    worker.stop = lambda: stopped.append(True)
    window._on_stop_recipe()
    assert stopped == [True]
    assert window.progress_dialog.stop_btn.text() == "Stopping..."
    assert window.progress_dialog.status_label.text() == "Aborting process..."
    window._current_worker = None
    window._on_stop_recipe()  # no worker: the dialog is closed
    assert window.progress_dialog is None


# --- Scripts and the log scripter -------------------------------------------------------------

@pytest.mark.parametrize("answer", ["Yes", "Cancel"])
def test_run_script(ready, monkeypatch, workers, answer):
    window = ready
    shown = _messages(window, monkeypatch, answer)
    window.run_script("smooth(sigma_kms=10)")
    assert shown[-1]["header"] == "Run script and replace history?"
    if answer == "Yes":
        assert len(workers) == 1 and type(workers[-1]).__name__ == "ScriptWorker"
        assert isinstance(window.progress_dialog, QProgressDialog)
        window._on_script_progress("line 1")
        assert window.progress_dialog.labelText() == "line 1"
    else:
        assert workers == []


def test_script_result_replaces_the_history(qtbot, ready):
    window = ready
    from astrocook.core.session_manager import SessionHistory
    old = window.active_history
    new = SessionHistory(old.states[0], HistoryLogV2())
    window._process_script_result(new)
    assert window.session_histories == [new] and window.active_history is new


def test_script_error(ready, monkeypatch):
    window = ready
    errors = []
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a: errors.append(a[2])))
    window._on_script_error((3, "smooth(x)", "bad"))
    assert errors == ["Failed on line 3:\nsmooth(x)\n\nError: bad"]


def test_log_scripter(ready):
    window = ready
    history = window.active_history
    window._launch_log_scripter(history)
    dialog = window.log_scripter_dialog
    assert dialog.isVisible() and dialog.log_object is history.log_manager
    window._launch_log_scripter(history)
    assert window.log_scripter_dialog is dialog
    assert dialog.windowTitle() == f"Session Log: {history.display_name}"
    # The log is code: a fixed-width font (visual review of phase 5-bis)
    from PySide6.QtGui import QFontInfo
    assert QFontInfo(dialog.scripter_edit.font()).fixedPitch()
    dialog.close()
    dialog.done(0)
    assert window.log_scripter_dialog is None


def test_dead_code_is_gone(window):
    # Called by nothing, removed (findings register F80)
    assert not hasattr(window, "_show_refit_warning_with_checkbox")
    assert not hasattr(window, "_SCRIPT_LINE_REGEX")
