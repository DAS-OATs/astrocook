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
        assert window._pending_recipe_on_properties_set == {
            "category": "continuum", "name": "estimate_auto", "params": {"a": 1},
            "mode": "direct"}
    else:
        assert opened == [] and window._pending_recipe_on_properties_set is None


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
    assert window._pending_recipe_on_properties_set is None
    requested = []
    monkeypatch.setattr(window, "_on_recipe_requested",
                        lambda cat, name, params, extra: requested.append((cat, name, params)))
    window._set_pending_action("continuum", "estimate_auto", {"b": 2}, "direct")
    window._process_recipe_result((state.with_new_spectrum(state.spec), "set_properties",
                                   {"z_em": "4.072"}))
    qtbot.waitUntil(lambda: bool(requested))
    assert requested == [("continuum", "estimate_auto", {"b": 2})]


# --- The recipe dialog ----------------------------------------------------------------------

def test_recipe_dialog(ready, monkeypatch):
    window = ready
    viewer = window.plot_viewer
    window._launch_recipe_dialog("flux", "smooth")
    dialog = window.active_recipe_dialog
    assert dialog.recipe_name == "smooth" and dialog.isVisible()
    # pan(False) and zoom(False) were meant to switch the tools off, but
    # they toggle them: zoom ends up on (findings register F81)
    assert viewer.toolbar.mode == "zoom rect"
    # A second dialog is not opened while one is
    window._launch_recipe_dialog("flux", "rebin")
    assert window.active_recipe_dialog is dialog
    dialog.reject()
    assert window.active_recipe_dialog is None


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
    assert window._pending_recipe_on_properties_set is None
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

def test_short_and_long_recipes(ready, workers):
    window = ready
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


def test_export_over_existing_files(ready, monkeypatch, tmp_path):
    # The question before overwriting fails: QMessageBox is imported again
    # further down in the same method, so it is a local name not yet bound
    # (findings register F83). Nothing is written
    window = ready
    name = window.active_history.current_state.name
    (tmp_path / f"{name}_spec.csv").write_text("old")
    _messages(window, monkeypatch, "Yes")
    _save_as(window, monkeypatch, tmp_path, f"{name}_spec.csv")
    with pytest.raises(UnboundLocalError):
        window._on_recipe_requested("file", "export_ascii", {"targets": "y, dy"}, {})
    assert (tmp_path / f"{name}_spec.csv").read_text() == "old"


def test_export_with_the_names_of_the_gui(ready, monkeypatch, tmp_path):
    # Targets named only as in the GUI (λ, F, dF) are not recognized as
    # columns of the spectrum: no _spec suggested, and an existing file is
    # overwritten without asking (findings register F82)
    window = ready
    name = window.active_history.current_state.name
    (tmp_path / f"{name}_spec.csv").write_text("old")
    shown = _messages(window, monkeypatch, "No")
    suggested = _save_as(window, monkeypatch, tmp_path)
    window._on_recipe_requested("file", "export_ascii", {"targets": "F, dF"}, {})
    assert suggested == [name]
    assert "Files already exist" not in [m["header"] for m in shown]
    assert (tmp_path / f"{name}_spec.csv").read_text() != "old"


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
    assert entry.recipe_name == "smooth" and entry.params == {"sigma_kms": "10"}


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


def test_cancelled_result(ready):
    window = ready
    history = window.active_history
    state = history.current_state
    cancelled = state.with_new_spectrum(state.spec)
    cancelled._stop_flag = True
    n = len(history.states)
    window._process_recipe_result((cancelled, "lya_auto", {}))
    assert len(history.states) == n


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


def test_stop(ready, workers):
    window = ready
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
    assert dialog.windowTitle() == f"Log Scripter: {history.display_name}"
    dialog.close()
    dialog.done(0)
    assert window.log_scripter_dialog is None


def test_dead_code(window):
    # Not called by anything (findings register F80)
    assert callable(window._show_refit_warning_with_checkbox)
    assert window._SCRIPT_LINE_REGEX.match("smooth(sigma=1)").groups() == ("smooth", "sigma=1")
