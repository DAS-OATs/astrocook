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
from PySide6.QtWidgets import QApplication, QMessageBox

from astrocook.gui.main_window import MainWindowV2
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
    """Run a recipe as the dialog does and wait until the window has
    processed its result (or its error).

    The window processes the result with a QTimer after the worker ends: a
    fixed wait after the worker was not always enough on a loaded machine
    (tests in parallel, coverage), so the test waits for the processing
    itself."""
    done = []
    for attr in ("_process_recipe_result", "_on_recipe_error"):
        method = getattr(win, attr)

        def wrapped(*args, _method=method, **kwargs):
            try:
                return _method(*args, **kwargs)
            finally:
                done.append(True)

        setattr(win, attr, wrapped)
    started = []
    start = win.thread_pool.start
    win.thread_pool.start = lambda *args: (started.append(True), start(*args))
    try:
        win._on_recipe_requested(category, name, params, {})
        if started:  # not when the GUI asks for something first (e.g. z_em)
            qtbot.waitUntil(lambda: bool(done), timeout=60000)
    finally:
        for attr in ("_process_recipe_result", "_on_recipe_error"):
            delattr(win, attr)
        del win.thread_pool.start
    qtbot.wait(10)  # let the signals queued by the processing be delivered


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


def _columns_fit_their_contents(inspector):
    view = inspector.table_view
    header = view.horizontalHeader()
    return all(
        header.sectionSize(col) == max(view.sizeHintForColumn(col),
                                       header.sectionSizeHint(col))
        for col in range(inspector.proxy_model.columnCount())
        if not view.isColumnHidden(col))  # btur, hidden when not in use


def test_inspector_columns_fit_their_contents(qtbot, loaded):
    """The columns of the inspector table fit the rows shown: all of them,
    then only those of the group (characterization, before the table was
    made faster in phase 5)."""
    loaded._on_view_system_inspector()
    inspector = loaded.system_inspector
    header = inspector.table_view.horizontalHeader()
    qtbot.waitUntil(lambda: _columns_fit_their_contents(inspector))
    wide = header.sectionSize(1)
    inspector.group_cb.setChecked(True)
    # The group has shorter series than the longest of the session
    qtbot.waitUntil(lambda: _columns_fit_their_contents(inspector)
                    and header.sectionSize(1) < wide)
    inspector.group_cb.setChecked(False)
    qtbot.waitUntil(lambda: _columns_fit_their_contents(inspector)
                    and header.sectionSize(1) == wide)


def test_inspector_boxes_above_the_table(loaded):
    """The boxes above the table do not touch, and the group box says what
    it shows (findings register F85: no spacing in their row, and the
    label "Group View (Selected Component Only)")."""
    loaded._on_view_system_inspector()
    inspector = loaded.system_inspector
    inspector.resize(1000, 700)
    QApplication.processEvents()
    group, flagged, button = (inspector.group_cb, inspector.flagged_cb,
                              inspector.select_flagged_btn)
    assert group.text() == "Only group" and "fitted together" in group.toolTip()
    assert group.geometry().right() + 10 < flagged.geometry().left()
    assert flagged.geometry().right() + 10 < button.geometry().left()


def test_inspector_fields_in_groups(loaded):
    # One row above the plot, as the boxes above the table, in three groups
    # (the view, a new component, the residuals) with tooltips (visual
    # review of phase 5-bis)
    from PySide6.QtWidgets import QFrame
    loaded._on_view_system_inspector()
    inspector = loaded.system_inspector
    inspector.resize(1200, 700)
    QApplication.processEvents()
    fields = (inspector.trans_in, inspector.z_in, inspector.vmin_in, inspector.vmax_in,
              inspector.logn_in, inspector.b_in, inspector.fit_window_in, inspector.resid_cb)
    assert all(f.toolTip() for f in fields)
    assert "right click" in inspector.logn_in.toolTip()
    tops = {f.mapTo(inspector, f.rect().topLeft()).y() // 10 for f in fields}
    assert len(tops) == 1  # one row
    lines = [w for w in inspector.findChildren(QFrame)
             if w.frameShape() == QFrame.Shape.VLine and w.isVisible()]
    xs = sorted(w.mapTo(inspector, w.rect().topLeft()).x() for w in lines)
    x_of = lambda w: w.mapTo(inspector, w.rect().topLeft()).x()  # noqa: E731
    assert len(xs) == 2
    assert x_of(inspector.vmax_in) < xs[0] < x_of(inspector.logn_in)
    assert x_of(inspector.fit_window_in) < xs[1] < x_of(inspector.resid_cb)


def test_inspector_btur_column(loaded, monkeypatch):
    # btur, almost always a frozen 0, hidden unless a component uses it;
    # the header menu shows or hides a column, remembered (visual review)
    import dataclasses
    from PySide6.QtCore import QSettings
    from astrocook.gui import system_inspector
    from astrocook.gui.inspector_table import COLUMNS
    loaded._on_view_system_inspector()
    inspector = loaded.system_inspector
    view, col = inspector.table_view, COLUMNS.index("btur")
    assert view.isColumnHidden(col)
    # A component with btur: shown
    model = inspector.table_model
    comps = list(model._components)
    comps[0] = dataclasses.replace(comps[0], btur=5.0)
    model.update_data(comps, model._constraints_map)
    inspector._update_column_visibility()
    assert not view.isColumnHidden(col)
    # The user hides it: hidden, and remembered
    inspector._set_column_shown("btur", False)
    assert view.isColumnHidden(col)
    assert "btur" in QSettings("Astrocook", "AstrocookV2").value(
        system_inspector.HIDDEN_COLUMNS_KEY)
    inspector._set_column_shown("logN", False)
    assert view.isColumnHidden(COLUMNS.index("logN"))
    # Default Columns: btur shown again (in use), logN back
    inspector._reset_columns()
    assert not view.isColumnHidden(col) and not view.isColumnHidden(COLUMNS.index("logN"))
    # The menu of the header: one entry per column but ID
    from PySide6.QtWidgets import QMenu
    shown = []

    class RecordedMenu(QMenu):
        def exec(self, *args):
            shown.append([x.text() for x in self.actions() if not x.isSeparator()])

    monkeypatch.setattr(system_inspector, "QMenu", RecordedMenu)
    inspector._on_header_menu(view.horizontalHeader().rect().center())
    assert shown[0][:len(COLUMNS) - 1] == COLUMNS[1:]
    assert shown[0][-1] == "Default Columns"


def test_inspector_fit_window(loaded, monkeypatch):
    # #518 (F4): the velocity window of the manual fits is visible and
    # editable; empty keeps the default of each recipe
    loaded._on_view_system_inspector()
    inspector = loaded.system_inspector
    calls = []
    monkeypatch.setattr(loaded, "_on_recipe_requested",
                        lambda menu, name, params, extra: calls.append((name, params)))
    assert inspector.fit_window_in.text() == ""
    assert inspector.fit_window_in.placeholderText() == "def."
    assert inspector.fit_window_params() == {}
    comps = list(loaded.active_history.current_state.systs.components)[:2]
    inspector._refit_list(comps)
    inspector._refit_all()
    assert "z_window_kms" not in calls[0][1] and calls[1][1] == {}
    inspector.fit_window_in.setText("40")
    inspector._refit_list(comps)
    inspector._refit_all()
    assert calls[2][1]["z_window_kms"] == "40.0"
    assert calls[2][1]["group_depth"] == "2"
    assert calls[3] == ("refit_all", {"z_window_kms": "40.0"})


def test_inspector_unfreeze_btur_with_its_column_hidden(loaded, monkeypatch):
    # With the btur column hidden, the menu of the table frees btur of the
    # selected components and shows the column (Guido, phase 5-bis)
    from PySide6.QtWidgets import QMenu
    from astrocook.gui import system_inspector
    from astrocook.gui.inspector_table import COLUMNS
    menus, calls = [], []

    class RecordedMenu(QMenu):
        def exec(self, *args):
            menus.append(self)

    monkeypatch.setattr(system_inspector, "QMenu", RecordedMenu)
    monkeypatch.setattr(loaded, "_on_recipe_requested",
                        lambda menu, name, params, extra: calls.append((name, params)))
    loaded._on_view_system_inspector()
    inspector = loaded.system_inspector
    view, col = inspector.table_view, COLUMNS.index("btur")
    assert view.isColumnHidden(col)
    view.selectRow(0)
    pos = view.visualRect(inspector.proxy_model.index(0, 0)).center()
    inspector._on_context_menu(pos)
    action = next(a for a in menus[-1].actions() if a.text() == "Unfreeze 'btur'")
    action.trigger()
    assert not view.isColumnHidden(col)
    uuid = inspector._get_selected_components()[0].uuid
    assert calls and calls[-1][0] == "update_constraint"
    assert calls[-1][1]["uuid"] == uuid and calls[-1][1]["param"] == "btur"
    # The column shown: the entry is not needed (the cells of btur have it)
    inspector._on_context_menu(pos)
    assert "Unfreeze 'btur'" not in [a.text() for a in menus[-1].actions()]


def test_status_bar_outcome_of_a_recipe(qtbot, loaded):
    # A recipe that succeeded said nothing: one line in the status bar,
    # also in the System Inspector (phase 5-bis)
    loaded._on_view_system_inspector()
    _run_recipe(qtbot, loaded, "edit", "set_properties", {"z_em": "4.072"})
    text = loaded.status_recipe.text()
    assert text.startswith("Set ") and ": done in " in text and text.endswith(" s")
    assert "components" not in text  # the number did not change
    assert loaded.system_inspector.status_recipe.text() == text


def test_status_bar_components_and_warnings(loaded):
    import time
    state = loaded.active_history.current_state
    n = len(state.systs.components)
    loaded._recipe_started = (time.perf_counter(), "absorbers", "clean_negligible", n + 3)
    loaded._captured_warnings = ["one", "two", "one"]
    loaded._report_recipe_outcome("clean_negligible", state)
    text = loaded.status_recipe.text()
    assert f"{n + 3} \u2192 {n} components" in text and text.endswith("with 2 warnings")
    # Only the recipe that was started
    loaded._recipe_started = (time.perf_counter(), "edit", "set_properties", n)
    loaded._report_recipe_outcome("smooth", state)
    assert loaded.status_recipe.text() == text


def test_coordinates_in_the_status_bars(loaded):
    # The coordinates under the mouse, from the toolbars to the status bars:
    # they stay when the toolbar is hidden (phase 5-bis)
    toolbar = loaded.plot_viewer.toolbar
    assert not toolbar.locLabel.isVisible()
    toolbar.set_message("λ = 500.0000 nm")
    assert loaded.status_coords.text() == "λ = 500.0000 nm"
    loaded.show_toolbar_action.trigger()  # hidden
    toolbar.set_message("λ = 501.0000 nm")
    assert loaded.status_coords.text() == "λ = 501.0000 nm"
    loaded.show_toolbar_action.trigger()
    loaded._on_view_system_inspector()
    inspector = loaded.system_inspector
    # As tall as the status bar of the window (it took the height left by
    # the splitter, reported by Guido)
    inspector.resize(1200, 900)
    QApplication.processEvents()
    assert inspector.status_bar.height() <= loaded.statusBar().height() + 4
    inspector.table_view.selectRow(0)
    QApplication.processEvents()
    vel = inspector.vel_plot
    text = vel.axes[0].format_coord(-12.0, 0.5)  # as in the main window
    assert text.startswith("Δv = -12.0 km/s, ƒ = 0.500, z = ")
    vel.toolbar.set_message(text)
    assert inspector.status_coords.text() == text


def test_inspector_keys(loaded, monkeypatch):
    # Return refits and Delete deletes the selected components, in the
    # table; Undo and Redo in the inspector too (phase 5-bis)
    from PySide6.QtGui import QKeySequence, QShortcut
    loaded._on_view_system_inspector()
    inspector = loaded.system_inspector
    calls = []
    monkeypatch.setattr(inspector, "_refit_list", lambda comps: calls.append(("refit", len(comps))))
    monkeypatch.setattr(inspector, "_delete_list", lambda comps: calls.append(("delete", len(comps))))
    monkeypatch.setattr(loaded, "_undo_last_action", lambda: calls.append(("undo",)))
    monkeypatch.setattr(loaded, "_redo_last_action", lambda: calls.append(("redo",)))
    inspector.table_view.selectRow(0)
    shortcuts = {s.key().toString(): s for s in inspector.findChildren(QShortcut)}

    def press(key):
        shortcuts[QKeySequence(key).toString()].activated.emit()

    press(QKeySequence.StandardKey.Undo)
    press(QKeySequence.StandardKey.Redo)
    for name in ("Return", "Enter"):
        press(name)
    assert calls[:2] == [("undo",), ("redo",)]
    assert calls[2:] == []  # no refit without a continuum (the archive has none)
    monkeypatch.setattr(inspector, "has_continuum", lambda: True)
    press("Return")
    press("Del")
    press("Backspace")
    assert calls[2:] == [("refit", 1), ("delete", 1), ("delete", 1)]
    # In the table only: not while typing in the boxes above the plot
    from PySide6.QtCore import Qt
    for name in ("Return", "Del"):
        assert shortcuts[name].context() == Qt.ShortcutContext.WidgetWithChildrenShortcut
        assert shortcuts[name].parent() is inspector.table_view


def _inspector_state(window, continuum=True, systems=True):
    """The session of the window with or without a continuum and systems."""
    import numpy as np
    from astrocook.core.structures import SystemListDataV2
    from astrocook.core.system_list import SystemListV2
    history = window.active_history
    state = history.current_state
    spec = state.spec
    if continuum:
        spec = spec.update_column("cont", np.full(len(spec.x), np.nanmedian(spec.y.value)))
    state = state.with_new_spectrum(spec)
    if not systems:
        state = state.with_new_system_list(SystemListV2(SystemListDataV2()))
    history.add_state(state)
    window._on_view_system_inspector()
    inspector = window.system_inspector
    inspector.set_session(state)
    QApplication.processEvents()
    return inspector


def test_inspector_without_a_continuum(loaded, monkeypatch):
    # A message instead of empty panels, with a button to estimate the
    # continuum; the fits disabled, they failed with an error of the fitter
    # (visual review of phase 5-bis)
    from PySide6.QtWidgets import QMenu
    from astrocook.gui import system_inspector
    inspector = _inspector_state(loaded, continuum=False)
    assert loaded.active_history.current_state.spec.norm is None
    vel = inspector.vel_plot
    inspector.table_view.selectRow(0)
    QApplication.processEvents()
    assert vel.plot_stack.currentWidget() is vel.message
    assert vel.message_label.text() == vel.NO_CONTINUUM
    assert not vel.continuum_button.isHidden()
    opened = []
    monkeypatch.setattr(loaded, "_launch_recipe_dialog",
                        lambda cat, name, initial_params=None: opened.append((cat, name)))
    vel.continuum_button.click()
    assert opened == [("continuum", "estimate_auto")]
    menus = []

    class RecordedMenu(QMenu):
        def exec(self, *args):
            menus.append(self)

    monkeypatch.setattr(system_inspector, "QMenu", RecordedMenu)
    pos = inspector.table_view.visualRect(inspector.proxy_model.index(0, 0)).center()
    inspector._on_context_menu(pos)
    actions = {a.text().split(" (")[0]: a for a in menus[-1].actions()}
    for name in ("Refit Selected", "Refit All", "Refit Bayesian"):
        assert not actions[name].isEnabled(), name
        assert actions[name].toolTip() == system_inspector.NO_CONTINUUM_TIP
    assert actions["Delete"].isEnabled()


def test_inspector_without_systems(loaded):
    # With a continuum and no systems: a message says how to start; the
    # transitions and the redshift typed above show the spectrum, also
    # without a component (it showed nothing)
    inspector = _inspector_state(loaded, continuum=True, systems=False)
    vel = inspector.vel_plot
    vel.plot_spectrum()
    QApplication.processEvents()
    assert vel.plot_stack.currentWidget() is vel.message
    assert vel.message_label.text() == vel.NO_SYSTEMS
    assert vel.continuum_button.isHidden()
    assert not inspector.empty_label.isHidden()
    inspector.trans_in.setText("CIV")
    inspector._apply()
    inspector.z_in.setText("2.0")
    inspector._apply_z()
    QApplication.processEvents()
    assert vel.plot_stack.currentWidget() is vel.canvas
    # The panels of a page, as many as fit the height of the window
    assert vel._all_transitions == ["CIV_1548", "CIV_1550"] and vel.axes
    assert {v[0] for v in vel._panel_map.values()} <= set(vel._all_transitions)


def test_inspector_with_a_continuum_and_systems(loaded):
    inspector = _inspector_state(loaded, continuum=True)
    inspector.table_view.selectRow(0)
    QApplication.processEvents()
    assert inspector.vel_plot.plot_stack.currentWidget() is inspector.vel_plot.canvas
    assert inspector.empty_label.isHidden()
    assert inspector.has_continuum()


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
    # The archive has no continuum, without which the fits are disabled
    monkeypatch.setattr(inspector, "has_continuum", lambda: True)
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


def test_gui_imports_after_pyplot():
    # F7: pyplot imported first, without a QApplication (a script, a test)
    # fixed a headless backend, and the GUI then failed at import, asking
    # matplotlib for QtAgg. The GUI no longer picks the backend of pyplot
    # (the launcher does); its canvases are QtAgg anyway
    code = ("import matplotlib.pyplot\n"
            "import astrocook.gui.main_window as mw\n"
            "from PySide6.QtWidgets import QApplication\n"
            "app = QApplication([])\n"
            "window = mw.MainWindowV2(None, None)\n")
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=120,
                            env={**os.environ, "QT_QPA_PLATFORM": "offscreen",
                                 "MPLBACKEND": "Agg"})
    assert result.returncode == 0, result.stderr[-2000:]


def test_plot_colours():
    # The colours of the plots (characterization, before pyplot was
    # dropped from pyside_plot.py for F7)
    import matplotlib
    from astrocook.gui.pyside_plot import get_color_cycle
    tab10 = matplotlib.colormaps["tab10"]
    assert get_color_cycle(4, cmap="tab10") == [tab10(i / 4) for i in range(4)]
    cycle = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
    assert get_color_cycle(12) == [cycle[i % len(cycle)] for i in range(12)]
    assert get_color_cycle(3, cmap="no_such_map") == [cycle[i] for i in range(3)]


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


def test_inspector_velocity_plot_is_drawn_up_to_date(qtbot, loaded):
    """Once the events are processed, the velocity plot shows its final
    state: the same pixels as a synchronous draw (characterization, before
    its draws were deferred in phase 5)."""
    import numpy as np
    loaded._on_view_system_inspector()
    inspector = loaded.system_inspector
    qtbot.wait(100)
    canvas = inspector.vel_plot.canvas
    assert inspector.vel_plot.axes  # the panels of the selected component
    shown = np.array(canvas.buffer_rgba()).copy()
    canvas.draw()
    assert np.array_equal(shown, np.array(canvas.buffer_rgba()))


def test_main_plot_highlights_follow_the_selection(qtbot, loaded):
    """After each selection of a sequence of rows of the inspector, the main
    plot shows the same pixels as a full rebuild of the plot for that
    selection (characterization, before the highlights were drawn without
    a full redraw in phase 5)."""
    import numpy as np
    loaded._on_view_system_inspector()
    inspector = loaded.system_inspector
    viewer = loaded.plot_viewer
    canvas = viewer.canvas
    qtbot.wait(100)
    shown = {}
    for row in (3, 40, 41, 0, 2):  # one selection after the other, as a user does
        inspector.table_view.selectRow(row)
        qtbot.wait(100)
        uuid = inspector._get_selected_components()[0].uuid
        shown[uuid] = np.array(canvas.buffer_rgba()).copy()
    for uuid, pixels in shown.items():
        # The same component again (by uuid: the sorting may move the rows)
        inspector.focus_on_component(uuid)
        qtbot.wait(100)
        assert inspector._get_selected_components()[0].uuid == uuid
        viewer.plot_spectrum(loaded.active_history.current_state)
        canvas.draw()
        qtbot.wait(100)
        assert np.array_equal(pixels, np.array(canvas.buffer_rgba())), uuid
