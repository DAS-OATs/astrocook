"""The sessions of the main window: opening (dialog, path, drag and drop),
the list (switching, renaming, the context menu), the view of a session
(the auxiliary columns offered, the empty view), a recipe result added to
the history or as a new branch, duplicating, stitching and co-adding,
the session inspector, closing (with or without saving).

Characterization tests written before the sessions were moved to a mixin
(phase 5, split of the GUI). Saving is also in ``test_gui_save_session.py``.
"""
import logging
import sys

import numpy as np
import pytest
from PySide6.QtCore import QMimeData, QPoint, QUrl
from PySide6.QtWidgets import QDialog, QLabel, QMenu, QMessageBox

from astrocook.core.structures import HistoryLogV2
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
def empty(qtbot, monkeypatch, tmp_path):
    return _new_window(qtbot, monkeypatch, tmp_path)


@pytest.fixture
def window(qtbot, monkeypatch, tmp_path):
    win = _new_window(qtbot, monkeypatch, tmp_path)
    win.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(50)
    return win


@pytest.fixture
def two(qtbot, window):
    window.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(50)
    return window


def _module(window, method):
    return sys.modules[getattr(type(window), method).__module__]


def _names(window):
    return [h.display_name for h in window.session_histories]


def _requests(window, monkeypatch):
    requests = []
    monkeypatch.setattr(window, "_on_recipe_requested",
                        lambda cat, name, params, extra: requests.append((cat, name, params)))
    return requests


def _answer(window, monkeypatch, button):
    asked = []

    def show(*args, **kwargs):
        asked.append(kwargs.get("header"))
        return getattr(QMessageBox.StandardButton, button)
    monkeypatch.setattr(window, "_show_custom_message", show)
    return asked


# --- Opening --------------------------------------------------------------------------------

def test_open_from_a_path(empty):
    empty.open_session_from_path(str(SESSION_FILE))
    assert _names(empty) == [SESSION_FILE.stem]
    assert empty.active_history is empty.session_histories[0]
    assert empty.central_stack.currentIndex() == 0
    assert empty.session_list_view.selectionModel().currentIndex().row() == 0
    assert SESSION_FILE.stem in empty.windowTitle()


def test_open_a_missing_or_bad_file(empty, monkeypatch, tmp_path, caplog):
    empty.open_session_from_path(str(tmp_path / "missing.acs2"))
    assert "does not exist" in caplog.text and empty.session_histories == []
    bad = tmp_path / "bad.fits"
    bad.write_text("not a spectrum")
    shown = []
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a: shown.append(a[1])))
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda *a: shown.append(a[1])))
    empty.open_session_from_path(str(bad))
    assert empty.session_histories == [] and len(shown) == 1


def test_open_with_the_dialog(empty, monkeypatch, tmp_path):
    from PySide6.QtCore import QSettings
    dialog = _module(empty, "_on_open_spectrum").QFileDialog
    bad = tmp_path / "bad.fits"
    bad.write_text("not a spectrum")
    monkeypatch.setattr(dialog, "getOpenFileNames", staticmethod(
        lambda *a, **k: ([str(SESSION_FILE), str(bad), str(SESSION_FILE)], "")))
    empty._on_open_spectrum()
    assert _names(empty) == [SESSION_FILE.stem] * 2  # the bad file skipped
    assert QSettings("Astrocook", "AstrocookV2").value("last_working_dir") == str(SESSION_FILE.parent)


def test_open_while_editing_the_continuum(qtbot, window, monkeypatch):
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072"})
    _run_recipe(qtbot, window, "continuum", "estimate_auto", {})
    window.edit_cont_button.click()
    _answer(window, monkeypatch, "Cancel")
    dialog = _module(window, "_on_open_spectrum").QFileDialog
    asked = []
    monkeypatch.setattr(dialog, "getOpenFileNames",
                        staticmethod(lambda *a, **k: asked.append(True) or ([], "")))
    window._on_open_spectrum()
    assert asked == []  # the editor is kept open, no file dialog


def test_drag_and_drop(empty, monkeypatch):
    opened = []
    monkeypatch.setattr(empty, "open_session_from_path", opened.append)

    class Event:
        def __init__(self, urls):
            self.data = QMimeData()
            if urls:
                self.data.setUrls([QUrl.fromLocalFile(u) for u in urls])
            self.accepted = self.ignored = False

        def mimeData(self):
            return self.data

        def acceptProposedAction(self):
            self.accepted = True

        def ignore(self):
            self.ignored = True

    event = Event(["/a/one.fits", "/a/two.acs2"])
    empty.dragEnterEvent(event)
    assert event.accepted
    empty.dropEvent(event)
    assert opened == ["/a/one.fits", "/a/two.acs2"]
    event = Event([])
    empty.dragEnterEvent(event)
    assert event.ignored


# --- Logs of new sessions ---------------------------------------------------------------------

def test_log_of_a_new_session(window):
    history = window.active_history
    log = history.log_manager
    assert isinstance(log, HistoryLogV2)
    # A duplicate does not share the log of its source
    window._on_duplicate_session(history)
    assert window.active_history.log_manager is not log


# --- The list -----------------------------------------------------------------------------------

def test_switch_session(two):
    first, second = two.session_histories
    assert two.active_history is second
    two._on_session_switched(two.session_model.index(0))
    assert two.active_history is first
    two._on_session_switched(two.session_model.index(5))  # not a session
    assert two.active_history is first


def test_switch_while_editing_the_continuum(qtbot, two, monkeypatch):
    _run_recipe(qtbot, two, "edit", "set_properties", {"z_em": "4.072"})
    _run_recipe(qtbot, two, "continuum", "estimate_auto", {})
    two.edit_cont_button.click()
    _answer(two, monkeypatch, "Cancel")
    second = two.active_history
    two._on_session_switched(two.session_model.index(0))
    assert two.active_history is second and two.plot_viewer._edit_mode_active
    assert two.session_list_view.selectionModel().currentIndex().row() == 1


def test_rename_in_the_list(window, monkeypatch):
    requests = _requests(window, monkeypatch)
    window.session_model.setData(window.session_model.index(0), "renamed")
    assert window.active_history.display_name == "renamed"
    assert requests == [("edit", "set_properties", {"name": "renamed", "object": "_current_",
                                                     "z_em": "_current_",
                                                     "resol": "_current_"})]


@pytest.fixture
def menus(window, monkeypatch):
    shown = []

    class RecordedMenu(QMenu):
        def exec(self, *args):
            shown.append({a.text(): a for a in self.actions() if a.text()})

    monkeypatch.setattr(_module(window, "_on_session_list_context_menu"), "QMenu", RecordedMenu)
    return shown


def test_context_menu_of_a_session(window, menus, monkeypatch):
    history = window.active_history
    called = []
    for name in ("_on_session_info", "_launch_log_scripter", "_launch_identification_viewer",
                 "_on_duplicate_session", "_on_save_session_context",
                 "_on_close_session_requested"):
        monkeypatch.setattr(window, name, lambda item, n=name: called.append((n, item)))
    window._on_session_list_context_menu(QPoint(5, 5))
    actions = menus[-1]
    assert list(actions) == ["View Info", "View Log", "View Identifications", "Duplicate",
                             "Save as...", "Close"]
    assert not actions["View Identifications"].isEnabled()  # no identifications
    for text in ("View Info", "View Log", "Duplicate", "Save as...", "Close"):
        actions[text].trigger()
    assert [n for n, _ in called] == ["_on_session_info", "_launch_log_scripter",
                                      "_on_duplicate_session", "_on_save_session_context",
                                      "_on_close_session_requested"]
    assert all(item is history for _, item in called)


def test_context_menu_of_several_sessions(two, monkeypatch):
    shown = []

    class RecordedMenu(QMenu):
        def exec(self, *args):
            shown.append({a.text(): a for a in self.actions() if a.text()})

    monkeypatch.setattr(_module(two, "_on_session_list_context_menu"), "QMenu", RecordedMenu)
    from PySide6.QtCore import QItemSelectionModel
    selection = two.session_list_view.selectionModel()
    selection.select(two.session_model.index(0), QItemSelectionModel.SelectionFlag.Select)
    called = []
    for name in ("_on_stitch_requested", "_on_equalized_stitch_dialog_requested",
                 "_on_coadd_requested"):
        monkeypatch.setattr(two, name, lambda primary, others, n=name: called.append(
            (n, primary, others)))
    monkeypatch.setattr(two, "_on_close_multiple_sessions", lambda rows: called.append(rows))
    two._on_session_list_context_menu(QPoint(5, 5))
    assert list(shown[-1]) == ["Stitch 2 sessions...", "Equalize and Stitch 2 Arms...",
                               "Co-add 2 sessions...", "Close 2 sessions"]
    for action in shown[-1].values():
        action.trigger()
    first = two.session_histories[0]
    assert called[:3] == [(n, first, [SESSION_FILE.stem]) for n in (
        "_on_stitch_requested", "_on_equalized_stitch_dialog_requested", "_on_coadd_requested")]
    assert called[3] == [0, 1]


# --- The view of a session ------------------------------------------------------------------

def test_auxiliary_columns_offered(window):
    history = window.active_history
    spec = history.current_state.spec
    n = len(spec.x.value)
    for name in ("zeta", "alpha", "running_std"):
        spec = spec.update_column(name, np.ones(n))
    history.add_state(history.current_state.with_new_spectrum(spec))
    combo = window.aux_col_combo
    window._update_view_for_session(history.current_state)
    items = [combo.itemText(i) for i in range(combo.count())]
    assert items[0] == "None" and items[1:] == sorted(items[1:])
    assert {"alpha", "running_std", "zeta"} <= set(items) and "x" not in items
    assert combo.currentText() == "None"
    window._update_view_for_session(history.current_state, auto_show_aux="zeta")
    assert combo.currentText() == "zeta"
    window._update_view_for_session(history.current_state)  # the choice is kept
    assert combo.currentText() == "zeta"


def test_empty_view(window):
    window._update_view_for_session(None)
    assert window.central_stack.currentIndex() == 1
    assert window.aux_col_combo.count() == 1


# --- Results of recipes ---------------------------------------------------------------------

def test_linear_result(window):
    history = window.active_history
    state = history.current_state
    new = state.with_new_spectrum(state.spec)
    new.name = "renamed"
    n = len(history.states)
    window.update_gui_session_state(new, 0, is_branching=False)
    assert len(history.states) == n + 1 and history.current_state is new
    assert history.display_name == "renamed" and _names(window) == ["renamed"]
    assert new._gui is window and window.undo_action.isEnabled()


def test_branching_result(window):
    history = window.active_history
    state = history.current_state
    window.update_gui_session_state(state.with_new_spectrum(state.spec), 0, is_branching=False)
    new = state.with_new_spectrum(state.spec)
    window.update_gui_session_state(new, 0, is_branching=True)
    base = history.states[0].name
    assert _names(window) == [history.display_name, f"{base}-1"]
    branch = window.active_history
    assert branch is window.session_histories[1] and branch.current_state is new
    # The states up to the current one are copied, with the name of the branch
    assert len(branch.states) == 3 and all(s.name == f"{base}-1" for s in branch.states)
    assert branch.log_manager is not history.log_manager
    # A branch of the branch is named after the branch
    window.update_gui_session_state(new.with_new_spectrum(new.spec), 0, is_branching=True)
    assert _names(window)[-1] == f"{base}-1-1"
    # Another branch of the first session takes the next number
    window.active_history = history
    window.update_gui_session_state(state.with_new_spectrum(state.spec), 0, is_branching=True)
    assert _names(window)[-1] == f"{base}-2"


# --- Duplicate, stitch, co-add ------------------------------------------------------------

def test_duplicate(window):
    history = window.active_history
    base = history.current_state.name
    window._on_duplicate_session(history)
    window._on_duplicate_session(history)
    window._on_duplicate_session(window.session_histories[1])  # a copy of a copy
    assert _names(window) == [base, f"{base}_copy1", f"{base}_copy2", f"{base}_copy3"]
    assert window.active_history is window.session_histories[-1]


def test_stitch_and_coadd(two, monkeypatch):
    requests = _requests(two, monkeypatch)
    opened = []
    monkeypatch.setattr(two, "_launch_recipe_dialog",
                        lambda cat, name, initial_params=None: opened.append(
                            (cat, name, initial_params)))
    first = two.session_histories[0]
    two._on_stitch_requested(first, ["other"])
    assert two.active_history is first
    assert requests == [("edit", "stitch", {"other_sessions": ["other"]})]
    two._on_coadd_requested(two.session_histories[1], ["a", "b"])
    assert two.active_history is two.session_histories[1]
    assert opened[-1] == ("edit", "coadd",
                          {"session_names": f"{two.session_histories[1].display_name}, a, b"})


def test_equalize_and_stitch(two, monkeypatch):
    from astrocook.recipes.edit import EDIT_RECIPES_SCHEMAS
    opened = []
    monkeypatch.setattr(two, "_launch_recipe_dialog",
                        lambda cat, name, initial_params=None: opened.append((cat, name)))
    params = EDIT_RECIPES_SCHEMAS["equalize_and_stitch"]["params"]
    docs = {p["name"]: p["doc"] for p in params}
    try:
        two._on_equalized_stitch_dialog_requested(two.session_histories[0], ["b", "c"])
        assert opened == [("edit", "equalize_and_stitch")]
        assert two._hidden_stitch_sessions == "b, c"
        # The documentation of the parameters is rewritten in the shared
        # schema (findings register F79)
        now = {p["name"]: p["doc"] for p in params}
        assert now["stitch_wavelengths"] == "Comma-separated cut-offs (requires 2 values)."
        two._on_equalized_stitch_dialog_requested(two.session_histories[0], ["b"])
        now = {p["name"]: p["doc"] for p in params}
        assert now["manual_factors"] == "Fixed multipliers (requires 1 value) or 'auto'."
    finally:
        for p in params:
            p["doc"] = docs[p["name"]]


# --- Session inspector -------------------------------------------------------------------------

@pytest.fixture
def inspector(window, monkeypatch):
    """The session inspector dialog, recorded instead of shown."""
    shown = []

    class RecordedDialog(QDialog):
        def exec(self):
            shown.append(self)

    monkeypatch.setattr(_module(window, "_on_session_info"), "QDialog", RecordedDialog)
    return shown


def test_session_inspector(window, inspector):
    history = window.active_history
    spec = history.current_state.spec
    window._on_session_info(history)
    dialog = inspector[-1]
    assert dialog.windowTitle() == f"Session Info: {history.display_name}"
    assert window.info_name_edit.text() == history.display_name
    assert window.info_object_edit.text() == str(spec.meta.get("OBJECT") or "undefined")
    assert window.info_zem_edit.text() == f"{spec._data.z_em:.5f}"
    stats = next(label.text() for label in dialog.findChildren(QLabel) if "Statistics" in label.text())
    x = spec.x.value
    assert f"{np.min(x):.2f} - {np.max(x):.2f} nm" in stats
    assert "<td>Reduced χ²:</td><td>N/A (no model)</td>" in stats
    assert f"<td>Pixels:</td><td>{len(x)}</td>" in stats
    assert f"<td>Components:</td><td>{len(history.current_state.systs.components)}</td>" in stats
    # The step of a logarithmic grid in velocity first (visual review)
    dv = 299792.458 * np.log(10) * np.median(np.diff(np.log10(x)))
    assert f"<td>Median step:</td><td>{dv:.2f} km/s (logarithmic grid, " in stats


def test_session_inspector_chi2(window, inspector):
    history = window.active_history
    state = history.current_state
    spec = state.spec.with_properties(z_em=4.072)
    model = np.convolve(spec.y.value, np.ones(5) / 5, mode="same")
    spec = spec.update_column("model", model)
    history.add_state(state.with_new_spectrum(spec))
    window._on_session_info(history)
    stats = next(label.text() for label in inspector[-1].findChildren(QLabel)
                 if "Statistics" in label.text())
    x, y, dy = spec.x.value, spec.y.value, spec.dy.value
    valid = (dy > 0) & np.isfinite(y) & np.isfinite(model) & np.isfinite(dy)
    chi2 = ((y - model) / dy)[valid] ** 2
    red = x[valid] >= 121.567 * (1 + 4.072)
    assert (f"<td>{chi2[~red].mean():.2f} / {chi2[red].mean():.2f}</td>") in stats


def test_session_inspector_save(window, inspector, monkeypatch):
    requests = _requests(window, monkeypatch)
    window._on_session_info(window.active_history)
    window.info_name_edit.setText("new name")
    window.info_object_edit.setText("undefined")
    window.info_zem_edit.setText("3,5")
    window._apply_session_info_changes(inspector[-1], window.active_history)
    # The resolution not edited is left as it is
    assert requests == [("edit", "set_properties", {"name": "new name", "object": "",
                                                     "z_em": "3.5", "resol": "_current_"})]
    window._on_session_info(window.active_history)
    assert window.info_resol_edit.text() == "undefined"  # no resolution
    window.info_resol_edit.setText("")
    window._apply_session_info_changes(inspector[-1], window.active_history)
    assert requests[-1][2]["resol"] == "0.0"
    window.info_resol_edit.setText("40000")
    window._apply_session_info_changes(inspector[-1], window.active_history)
    assert requests[-1][2]["resol"] == "40000"


# --- Saving and closing ---------------------------------------------------------------------

def test_save_actions(window, monkeypatch):
    saved = []
    monkeypatch.setattr(window, "_save_session", saved.append)
    window._on_save_session()
    other = window.active_history
    window._on_save_session_context(other)
    assert saved == [window.active_history, other]


def test_close_a_clean_session(two):
    first, second = two.session_histories
    two._on_close_session_requested(False)  # from the menu: the active one
    assert two.session_histories == [first] and two.active_history is first
    two._on_close_session_requested(first)
    assert two.session_histories == [] and two.active_history is None
    assert two.central_stack.currentIndex() == 1


@pytest.mark.parametrize("answer, saved_ok, closed", [
    ("Save", True, True), ("Save", False, False), ("Discard", None, True),
    ("Cancel", None, False)])
def test_close_a_changed_session(qtbot, window, monkeypatch, answer, saved_ok, closed):
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072"})
    history = window.active_history
    asked = _answer(window, monkeypatch, answer)
    saves = []
    monkeypatch.setattr(window, "_save_session", lambda h: saves.append(h) or saved_ok)
    window._on_close_session_requested(history)
    assert asked == [f"Do you want to save changes to '{history.display_name}'?"]
    assert saves == ([history] if answer == "Save" else [])
    assert (history not in window.session_histories) == closed


def test_close_a_background_session(two):
    first, second = two.session_histories
    two._close_session(first)
    assert two.session_histories == [second] and two.active_history is second


def test_close_several_sessions(qtbot, two, monkeypatch):
    two.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(50)
    closed = []
    monkeypatch.setattr(two, "_on_close_session_requested", closed.append)
    histories = list(two.session_histories)
    two._on_close_multiple_sessions([0, 2])
    assert closed == [histories[2], histories[0]]  # from the last


def test_close_without_a_session(empty, caplog):
    with caplog.at_level(logging.WARNING):
        empty._on_close_session_requested()
    assert "Close request with no session specified or active." in caplog.text


def test_close_the_active_session_in_the_middle(qtbot, two):
    two.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(50)
    first, middle, last = two.session_histories
    two._on_session_switched(two.session_model.index(1))
    two._close_session(middle)
    # The session now at the same row becomes active
    assert two.session_histories == [first, last] and two.active_history is last
    assert two.session_list_view.selectionModel().currentIndex().row() == 1
