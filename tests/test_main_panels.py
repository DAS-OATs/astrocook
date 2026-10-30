"""The panels of the main window: the plot controls (each control redraws
the plot), the redshift cursor, the axis limits typed in the boxes and
shown from the plot, normalized and SNR views, the error columns of the
SNR view, the sidebars opened and closed, the state of the window with
and without a session, and the welcome view.

Characterization tests written before the panels were moved to a mixin
(phase 5, split of the GUI). The cards (their height, fading, sections
and swatches) are in ``test_gui_cards.py``.
"""
import numpy as np
import pytest
from PySide6.QtWidgets import QLabel, QMessageBox

from astrocook.gui.main_panels import LEFT_SIDEBAR_WIDTH, RIGHT_SIDEBAR_WIDTH
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
    win.resize(1400, 900)
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
def replots(window, monkeypatch):
    """The plots asked to the viewer (the session state of each)."""
    calls = []
    monkeypatch.setattr(window.plot_viewer, "plot_spectrum",
                        lambda session_state=None, **k: calls.append(session_state))
    return calls


# --- Plot controls -----------------------------------------------------------------------

def test_each_control_redraws_the_plot(window, replots):
    state = window.active_history.current_state
    for box in ("error_checkbox", "error_line_checkbox", "continuum_checkbox", "model_checkbox",
                "systems_checkbox", "strong_lines_checkbox", "log_x_checkbox", "log_y_checkbox"):
        before = len(replots)
        getattr(window, box).toggle()
        assert len(replots) == before + 1 and replots[-1] is state, box
    before = len(replots)
    window.x_unit_combo.setCurrentIndex(1)
    window.aux_col_combo.addItem("xyz")
    window.aux_col_combo.setCurrentIndex(window.aux_col_combo.findText("xyz"))
    assert len(replots) == before + 2


def test_isomag_grid(window, monkeypatch):
    toggled = []
    monkeypatch.setattr(window.plot_viewer, "toggle_isomag_grid", toggled.append)
    window.isomag_checkbox.setChecked(True)
    window.isomag_checkbox.setChecked(False)
    assert toggled == [True, False]


def test_error_style_redraws_the_inspector(qtbot, window, monkeypatch):
    window._on_view_system_inspector()
    updates = []
    monkeypatch.setattr(window.system_inspector.vel_plot, "_update_plot",
                        lambda: updates.append(True))
    window.error_line_checkbox.toggle()
    assert updates == [True]
    window.system_inspector.hide()
    window.error_line_checkbox.toggle()
    assert updates == [True]  # not when hidden


def test_redraw_while_editing_the_continuum(window, monkeypatch, replots):
    refreshed = []
    viewer = window.plot_viewer
    monkeypatch.setattr(viewer, "refresh_edit_layer", lambda: refreshed.append(True))
    viewer._edit_mode_active = True
    window._trigger_replot()
    viewer._edit_mode_active = False
    window._trigger_replot()
    assert refreshed == [True] and len(replots) == 1


def test_log_f_is_in_the_panel(window):
    # Below "Logarithmic λ" (before, it was not in the panel: findings
    # register F78)
    boxes = window.right_sidebar_widget.findChildren(type(window.log_y_checkbox))
    assert window.log_y_checkbox in boxes
    layout = window.log_x_checkbox.parentWidget().layout()
    assert (layout.indexOf(window.log_y_checkbox)
            == layout.indexOf(window.log_x_checkbox) + 1)


def test_custom_limits_with_log_f(window, caplog):
    import logging
    ax = window.plot_viewer.canvas.axes
    ax.set_yscale("log")
    ax.set_ylim(1.0, 100.0)
    window.log_y_checkbox.blockSignals(True)
    window.log_y_checkbox.setChecked(True)
    window.log_y_checkbox.blockSignals(False)
    window.ymin_input.setText("-5")
    window.ymax_input.setText("50")
    with caplog.at_level(logging.WARNING):
        window._on_set_custom_limits()
    # A lower limit that is not positive is not applied
    assert ax.get_ylim() == pytest.approx((1.0, 50.0))
    assert "must be positive" in caplog.text


# --- Normalized and SNR ------------------------------------------------------------------

def test_normalized_and_snr_exclude_each_other(window, replots):
    window.norm_y_checkbox.setChecked(True)
    assert window.lbl_flux.text() == "ƒ:"
    window.snr_checkbox.setChecked(True)
    assert not window.norm_y_checkbox.isChecked() and window.lbl_flux.text() == "F:"
    window.norm_y_checkbox.setChecked(True)
    assert not window.snr_checkbox.isChecked()
    window.norm_y_checkbox.setChecked(False)
    assert window.lbl_flux.text() == "F:"
    assert len(replots) == 4


def test_error_columns_of_the_snr_view(window):
    history = window.active_history
    spec = history.current_state.spec
    n = len(spec.x.value)
    for name in ("running_std", "y_err", "ERRsky", "cont"):
        spec = spec.update_column(name, np.ones(n))
    history.add_state(history.current_state.with_new_spectrum(spec))
    window._populate_snr_combo()
    combo = window.snr_col_combo
    assert [combo.itemText(i) for i in range(combo.count())] == ["dF", "running_std", "y_err",
                                                                 "ERRsky"]
    combo.setCurrentIndex(2)
    window._populate_snr_combo()
    assert combo.currentText() == "y_err"  # kept
    history.undo()
    window._populate_snr_combo()
    assert combo.count() == 1 and combo.currentText() == "dF"


def test_flux_normalized_needs_a_continuum(qtbot, window):
    assert not window.norm_y_checkbox.isEnabled()
    assert window.norm_y_checkbox.toolTip() == "Requires a continuum to normalize F."
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072"})
    _run_recipe(qtbot, window, "continuum", "estimate_auto", {})
    assert window.norm_y_checkbox.isEnabled()
    assert window.norm_y_checkbox.toolTip() == "Plot F / Continuum and set F limits."


# --- Redshift cursor -----------------------------------------------------------------------

def test_cursor_inputs(window, replots):
    window.cursor_z_input.setText("3.0")
    window.cursor_z_input.editingFinished.emit()
    assert replots == []  # cursor hidden
    window.cursor_series_input.returnPressed.emit()
    assert window.cursor_show_checkbox.isChecked() and len(replots) == 1
    window.cursor_z_input.editingFinished.emit()
    assert len(replots) == 2


@pytest.mark.parametrize("answer", ["Save", "Discard", "Cancel"])
def test_cursor_while_editing_the_continuum(window, monkeypatch, replots, answer):
    viewer = window.plot_viewer
    viewer.start_continuum_edit(initial_stride=500)
    window.edit_cont_button.setText("Save")
    saved = []
    monkeypatch.setattr(window, "_toggle_continuum_editor", lambda: saved.append(True))
    monkeypatch.setattr(window, "_show_custom_message",
                        lambda *a, **k: getattr(QMessageBox.StandardButton, answer))
    window.cursor_show_checkbox.setChecked(True)
    if answer == "Save":
        assert saved == [True] and window.cursor_show_checkbox.isChecked()
    elif answer == "Discard":
        assert not viewer._edit_mode_active and window.edit_cont_button.text() == "Start"
        assert window.cursor_show_checkbox.isChecked()
    else:
        assert viewer._edit_mode_active and not window.cursor_show_checkbox.isChecked()
        assert replots == []


# --- Axis limits ---------------------------------------------------------------------------

@pytest.mark.parametrize("unit, factor", [("nm", 1.0), ("Angstrom", 10.0), ("micron", 1e-3)])
def test_limits_typed_in_the_boxes(window, unit, factor):
    window.x_unit_combo.blockSignals(True)
    window.x_unit_combo.setCurrentIndex(window.x_unit_combo.findText(unit))
    window.x_unit_combo.blockSignals(False)
    ax = window.plot_viewer.canvas.axes
    for box, value in (("xmin_input", 500.0 * factor), ("xmax_input", 520.0 * factor),
                       ("ymin_input", 0.0), ("ymax_input", 2e-17)):
        getattr(window, box).setText(str(value))
    window._on_set_custom_limits()
    assert ax.get_xlim() == pytest.approx((500.0, 520.0)) and ax.get_ylim() == (0.0, 2e-17)
    # The boxes show the limits of the plot, in the unit shown
    window.xmin_input.setText("0")
    window._update_limit_boxes_from_plot()
    assert window.xmin_input.text() == f"{500.0 * factor:.4g}"
    assert window.xmax_input.text() == f"{520.0 * factor:.4g}"
    assert window.ymax_input.text() == f"{2e-17:.4g}"


def test_invalid_limits_are_ignored(window, caplog):
    ax = window.plot_viewer.canvas.axes
    before = ax.get_xlim(), ax.get_ylim()
    window.xmin_input.setText("abc")
    window._on_set_custom_limits()
    window.xmin_input.setText("600")
    window.xmax_input.setText("500")
    window._on_set_custom_limits()
    assert (ax.get_xlim(), ax.get_ylim()) == before
    assert "Invalid custom limit values entered." in caplog.text
    assert "Invalid custom limit range (max <= min)." in caplog.text


# --- Sidebars ------------------------------------------------------------------------------

def _settle(qtbot, window):
    for group in (window.left_animation_group, window.right_animation_group):
        if group is not None:
            qtbot.waitUntil(lambda g=group: g.state() != g.State.Running, timeout=2000)
    qtbot.wait(20)


def test_sidebars_open_and_close(qtbot, window, monkeypatch):
    _settle(qtbot, window)
    left, right = window.left_sidebar_widget, window.right_sidebar_widget
    open_left, open_right = left.width() > 0, right.width() > 0
    replots = []
    monkeypatch.setattr(window, "_trigger_replot", lambda *a: replots.append(True))
    window._toggle_left_sidebar()
    window._toggle_right_sidebar()
    _settle(qtbot, window)
    assert (left.width() > 0) == (not open_left) and (right.width() > 0) == (not open_right)
    assert left.isVisible() == (not open_left) and right.isVisible() == (not open_right)
    window._toggle_left_sidebar()
    window._toggle_right_sidebar()
    _settle(qtbot, window)
    assert left.width() == (LEFT_SIDEBAR_WIDTH if open_left else 0)
    assert right.width() == (RIGHT_SIDEBAR_WIDTH if open_right else 0)
    # The tab of each card is next to its edge; the plot is not redrawn
    # (it was, to move the legend out of the card: the plot has no legend
    # now, F86)
    assert window.session_collapse_button.x() == left.width()
    assert window.plot_controls_collapse_button.x() == window.width() - right.width() - 28
    assert window.session_collapse_button.toolTip() == (
        "Collapse Session Panel" if open_left else "Expand Session Panel")
    assert replots == []


def test_restack_keeps_the_size(qtbot, window):
    size = window.size()
    window._force_restack_floating_widgets()
    assert window.size() == size


# --- Empty and loaded window ----------------------------------------------------------------

def test_welcome_view(empty):
    assert empty.central_stack.currentIndex() == 1
    texts = [label.text() for label in empty.central_stack.currentWidget().findChildren(QLabel)]
    assert "Welcome to" in texts and "Use 'File > Open Session...'" in texts
    assert not empty.session_collapse_button.isVisible()
    assert not empty.plot_controls_collapse_button.isVisible()
    assert not empty.left_sidebar_widget.isVisible() and not empty.right_sidebar_widget.isVisible()


def test_state_with_and_without_a_session(window):
    assert window.central_stack.currentIndex() == 0
    assert window.session_collapse_button.isVisible()
    assert window.plot_controls_collapse_button.isVisible()
    window._update_ui_state(False)
    assert window.central_stack.currentIndex() == 1
    assert not window.session_collapse_button.isVisible()
    assert not window.left_sidebar_widget.isVisible() and window.left_sidebar_widget.width() == 0
    assert not window.save_action.isEnabled() and not window.norm_y_checkbox.isEnabled()
    window._update_ui_state(True)
    assert window.central_stack.currentIndex() == 0 and window.save_action.isEnabled()


def test_the_drop_down_menus_line_up(qtbot, window):
    # Aux. Column, λ Unit and Error Col.: the same left edge and width (F88:
    # three form layouts, each with the width of its own label)
    if window.right_sidebar_widget.width() == 0:
        window._toggle_right_sidebar()
    for section in window.right_sidebar_widget.findChildren(type(window.plot_elements_section)):
        section.setExpanded(True)
    qtbot.wait(400)
    combos = [window.aux_col_combo, window.x_unit_combo, window.snr_col_combo]
    lefts = {c.mapTo(window, c.rect().topLeft()).x() for c in combos}
    widths = {c.width() for c in combos}
    assert len(lefts) == 1 and len(widths) == 1, (lefts, widths)
