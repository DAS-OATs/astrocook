"""The mouse on the main plot: the tooltip after a pause (component, region
of an identification, data near the mouse), the menu of the right click,
and the knots of the continuum editor moved, added and removed with the
mouse.

Characterization tests written before the mouse handlers and the tooltip
were moved to a mixin (phase 5, split of the GUI). The menus and the
tooltips are recorded instead of shown; the names they use (QMenu,
QToolTip, QApplication) are replaced in the module that defines the
handlers, wherever it is.
"""
import inspect
import json

import numpy as np
import pytest
from matplotlib.backend_bases import MouseEvent
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMenu, QMessageBox

from astrocook.gui.main_window import MainWindowV2

from test_gui_smoke import SESSION_FILE, _run_recipe


@pytest.fixture
def window(qtbot, monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    for kind in ("warning", "critical", "information", "question"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(
            lambda *a, **k: QMessageBox.StandardButton.Ok))
    win = MainWindowV2(None, None)
    qtbot.addWidget(win)
    monkeypatch.setattr(win, "_show_custom_message", lambda *a, **k: None)
    win.show()
    win.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(100)
    return win


def _module(viewer, method):
    return inspect.getmodule(getattr(type(viewer), method))


@pytest.fixture
def modifiers(window, monkeypatch):
    """The keyboard modifiers seen by the handlers (none by default)."""
    held = [Qt.KeyboardModifier.NoModifier]
    viewer = window.plot_viewer

    class Application:
        keyboardModifiers = staticmethod(lambda: held[0])

    for method in ("on_press", "show_data_tooltip"):
        monkeypatch.setattr(_module(viewer, method), "QApplication", Application)
    return held


@pytest.fixture
def tips(window, monkeypatch):
    """The text of each tooltip shown (None when hidden)."""
    shown = []

    class Tip:
        showText = staticmethod(lambda pos, text, widget=None: shown.append(text))
        hideText = staticmethod(lambda: shown.append(None))

    monkeypatch.setattr(_module(window.plot_viewer, "_update_tooltip_display"), "QToolTip", Tip)
    return shown


@pytest.fixture
def menus(window, monkeypatch):
    """The actions of each menu that would be shown, by text (and the mode
    of the toolbar while it is shown, in ``tool_modes``)."""
    shown = []
    toolbar = window.plot_viewer.toolbar

    class RecordedMenu(QMenu):
        def exec(self, *args):
            shown.append({a.text(): a for a in self.actions() if a.text()})
            tool_modes.append(toolbar.mode)

    tool_modes = []
    monkeypatch.setattr(_module(window.plot_viewer, "on_press"), "QMenu", RecordedMenu)
    RecordedMenu.tool_modes = tool_modes
    return shown


def _event(canvas, name, x, y, button=None):
    x_disp, y_disp = canvas.axes.transData.transform((x, y))
    return MouseEvent(name, canvas, x_disp, y_disp, button=button)


def _spec(window):
    return window.active_history.current_state.spec


def _zoom(window, qtbot, lo, hi):
    """Zoom on [lo, hi] nm (the plot is rebuilt after the change of limits)."""
    window.plot_viewer.canvas.axes.set_xlim(lo, hi)
    qtbot.wait(100)


def _free_x(viewer, fraction=0.5):
    """A wavelength of the view away from every tick."""
    ax = viewer.canvas.axes
    lo, hi = ax.get_xlim()
    ticks = np.array([x for a, x, _ in viewer._tick_visuals if a is ax])
    grid = np.linspace(lo, hi, 2000)[100:-100]
    grid = np.r_[grid[int(len(grid) * fraction):], grid[:int(len(grid) * fraction)]]
    for x in grid:
        # The menu looks for a tick within 1.5% of the view
        if ticks.size == 0 or np.min(np.abs(ticks - x)) > 0.02 * (hi - lo):
            return float(x)
    raise AssertionError("no free wavelength")


# --- Tooltip ---------------------------------------------------------------------------

def test_tooltip_of_a_component(window, tips, modifiers):
    viewer = window.plot_viewer
    ax, x, comp = viewer._tick_visuals[0]
    _, y_hi = ax.get_ylim()
    viewer.show_data_tooltip(_event(viewer.canvas, "motion_notify_event", x, 0.9 * y_hi))
    assert tips[-1] == (
        "<div style='font-size: 13px'>"
        f"<b>ID:</b> {comp.id} | <b>{comp.series}</b><br>"
        f"<b>z:</b> {comp.z:.5f}<br>"
        "<b>χ²:</b> N/A | <b>R:</b> spectrum"
        "</div>")


def test_tooltip_of_the_flux_when_zoomed(window, qtbot, tips, modifiers):
    viewer = window.plot_viewer
    spec = _spec(window)
    x, y = spec.x.value, spec.y.value
    i = 1000
    _zoom(window, qtbot, x[i - 20], x[i + 20])
    viewer.show_data_tooltip(_event(viewer.canvas, "motion_notify_event", x[i], y[i]))
    unit = str(spec.y.unit)
    rows = (f"<tr><td style='padding-right:8px'><b>λ:</b></td><td>{x[i]:.4f} nm</td></tr>"
            f"<tr><td style='padding-right:8px'><b>dλ:</b></td>"
            f"<td>{spec.xmax.value[i] - spec.xmin.value[i]:.4f} nm</td></tr>"
            f"<tr><td style='padding-right:8px'><b>F:</b></td><td>{y[i]:.4e} {unit}</td></tr>"
            f"<tr><td style='padding-right:8px'><b>dF:</b></td>"
            f"<td>{spec.dy.value[i]:.4e} {unit}</td></tr>")
    assert tips[-1].startswith(f"<table style='font-size: 13px'>{rows}</table>")


def test_tooltip_far_from_the_data(window, qtbot, tips, modifiers):
    viewer = window.plot_viewer
    x = _spec(window).x.value
    _zoom(window, qtbot, x[980], x[1020])
    lo, hi = viewer.canvas.axes.get_ylim()
    xf = _free_x(viewer)
    viewer.show_data_tooltip(_event(viewer.canvas, "motion_notify_event", xf, hi - 0.01 * (hi - lo)))
    assert viewer._tip_data_rows == [] and viewer._tip_component_html is None
    assert tips[-1] is None  # nothing to show: hidden


def test_tooltip_of_the_whole_spectrum(window, tips, modifiers):
    # The whole spectrum is visible (2156 pixels, below the limit of 20000):
    # the data rows are given
    viewer = window.plot_viewer
    spec = _spec(window)
    i = 1000
    viewer.show_data_tooltip(_event(viewer.canvas, "motion_notify_event",
                                    spec.x.value[i], spec.y.value[i]))
    assert [r[0] for r in viewer._tip_data_rows][:3] == ["λ", "dλ", "F"]


def test_tooltip_of_the_continuum_normalized(qtbot, window, tips, modifiers):
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072", "resol": "45000"})
    _run_recipe(qtbot, window, "continuum", "estimate_auto", {})
    viewer = window.plot_viewer
    window.norm_y_checkbox.setChecked(True)
    qtbot.wait(100)
    spec = _spec(window)
    i = 1500
    _zoom(window, qtbot, spec.x.value[i - 20], spec.x.value[i + 20])
    viewer.show_data_tooltip(_event(viewer.canvas, "motion_notify_event", spec.x.value[i], 1.0))
    labels = [r[0] for r in viewer._tip_data_rows]
    assert "Cont" in labels
    assert viewer._tip_data_rows[labels.index("Cont")] == ("Cont", f"{1.0:.4e}", "")
    # The flux, normalized, near the mouse
    y = spec.y.value[i] / spec.norm.value[i]
    viewer.show_data_tooltip(_event(viewer.canvas, "motion_notify_event", spec.x.value[i], y))
    assert ("ƒ", f"{y:.4e}", "") in viewer._tip_data_rows
    assert ("dƒ", f"{spec.dy.value[i] / spec.norm.value[i]:.4e}", "") in viewer._tip_data_rows


def test_tooltip_of_a_region(window, qtbot, tips, modifiers):
    viewer = window.plot_viewer
    spec = _spec(window)
    i = 1000
    ids = np.zeros(len(spec.x.value))
    ids[i - 5:i + 5] = 7
    spec._data.aux_cols['abs_ids'] = ids
    spec._data.meta['region_identifications'] = json.dumps(
        json.dumps({"7": [["CIV", 0.83], ["Ly_a", 0.0]]}))  # encoded twice
    _zoom(window, qtbot, spec.x.value[i - 20], spec.x.value[i + 20])
    viewer.show_data_tooltip(_event(viewer.canvas, "motion_notify_event",
                                    spec.x.value[i], spec.y.value[i]))
    region = ("<table style='font-size: 13px'><tr><td style='padding-right:8px'>"
              "<b>Likely:</b></td><td>CIV (0.83), Ly_a</td></tr></table>")
    assert f"</table><hr style='margin: 4px 0px;'>{region}" in tips[-1]
    # Outside the region: no line
    viewer.show_data_tooltip(_event(viewer.canvas, "motion_notify_event",
                                    spec.x.value[i + 10], spec.y.value[i + 10]))
    assert "Likely" not in (tips[-1] or "")


def test_no_tooltip_with_ctrl_or_the_cursor(window, tips, modifiers):
    viewer = window.plot_viewer
    spec = _spec(window)
    event = _event(viewer.canvas, "motion_notify_event", spec.x.value[1000], spec.y.value[1000])
    modifiers[0] = Qt.KeyboardModifier.ControlModifier
    viewer.show_data_tooltip(event)
    modifiers[0] = Qt.KeyboardModifier.NoModifier
    window.cursor_show_checkbox.setChecked(True)
    viewer.show_data_tooltip(event)
    assert tips == []


def test_tooltip_order(window, tips):
    viewer = window.plot_viewer
    viewer._tip_data_rows = [("F", "1.0", "u")]
    viewer._tip_region_html = "<tr><td>R</td></tr>"
    viewer._tip_component_html = "<div>C</div>"
    viewer._update_tooltip_display()
    sep = "<hr style='margin: 4px 0px;'>"
    assert tips[-1] == ("<table style='font-size: 13px'><tr><td style='padding-right:8px'>"
                        "<b>F:</b></td><td>1.0 u</td></tr></table>" + sep
                        + "<table style='font-size: 13px'><tr><td>R</td></tr></table>" + sep
                        + "<div>C</div>")
    viewer._tip_data_rows = []
    viewer._tip_region_html = None
    viewer._update_tooltip_display()
    assert tips[-1] == "<div>C</div>"  # no separator
    viewer._tip_component_html = None
    viewer._update_tooltip_display()
    assert tips[-1] is None


# --- Right click -------------------------------------------------------------------------

def test_right_click_on_a_component(window, menus, modifiers, monkeypatch):
    viewer = window.plot_viewer
    ax, x, comp = viewer._tick_visuals[0]
    opened = []
    monkeypatch.setattr(window, "open_inspector_on_component", opened.append)
    _, y_hi = ax.get_ylim()
    assert viewer.on_press(_event(viewer.canvas, "button_press_event", x, 0.9 * y_hi, button=3))
    texts = list(menus[-1])
    spec = _spec(window)
    idx = int(np.clip(np.searchsorted(spec.x.value, x), 0, len(spec.x.value) - 1))
    assert texts == [f"Inspect {comp.series} (its group)", "Split from Current Zoom",
                     f"Inspect Data (near {spec.x.value[idx]:.4f} nm)"]
    menus[-1][texts[0]].trigger()
    assert opened == [comp.uuid]


def test_right_click_with_the_cursor(window, menus, modifiers, monkeypatch):
    viewer = window.plot_viewer
    requests, splits, inspected = [], [], []
    monkeypatch.setattr(window, "_on_recipe_requested",
                        lambda cat, name, params, extra: requests.append((cat, name, params)))
    monkeypatch.setattr(window, "launch_split_from_current_view", lambda: splits.append(True))
    monkeypatch.setattr(viewer, "_open_data_inspector", lambda s, i: inspected.append(i))
    window.cursor_series_input.setText("CIV")
    window.cursor_show_checkbox.setChecked(True)
    x = _free_x(viewer)
    lo, hi = viewer.canvas.axes.get_ylim()
    assert viewer.on_press(_event(viewer.canvas, "button_press_event", x, 0.5 * (lo + hi), button=3))
    texts = list(menus[-1])
    z = viewer.calculate_z_from_x(x)
    assert texts[:3] == [f"Add CIV at z={z:.5f}", f"Set Emission Redshift to z={z:.5f}",
                         "Split from Current Zoom"]
    assert texts[3].startswith("Inspect Data (near ")
    for text in texts:
        menus[-1][text].trigger()
    assert requests[0][:2] == ("absorbers", "add_component")
    assert requests[0][2]["series"] == "CIV" and requests[0][2]["z"] == pytest.approx(z)
    assert requests[1][:2] == ("edit", "set_properties")
    assert float(requests[1][2]["z_em"]) == pytest.approx(z, abs=1e-12)
    assert splits == [True]
    spec = _spec(window)
    assert inspected == [int(np.clip(np.searchsorted(spec.x.value, x), 0, len(spec.x.value) - 1))]


def test_right_click_with_a_tool(window, menus, modifiers):
    viewer = window.plot_viewer
    viewer.toolbar.zoom()
    x = _free_x(viewer)
    lo, hi = viewer.canvas.axes.get_ylim()
    event = _event(viewer.canvas, "button_press_event", x, 0.5 * (lo + hi), button=3)
    assert viewer.on_press(event) is False  # left to the tool
    assert menus == []
    # With Ctrl: the menu, the tool paused while it is shown and restored
    modifiers[0] = Qt.KeyboardModifier.ControlModifier
    assert viewer.on_press(event)
    assert len(menus) == 1
    assert _module(viewer, "on_press").QMenu.tool_modes == [""]
    assert viewer.toolbar.mode == "zoom rect"
    viewer.toolbar.zoom()


def test_left_click_does_nothing(window, menus, modifiers):
    viewer = window.plot_viewer
    lo, hi = viewer.canvas.axes.get_ylim()
    assert viewer.on_press(_event(viewer.canvas, "button_press_event", _free_x(viewer),
                                  0.5 * (lo + hi), button=1)) is False
    assert menus == []


# --- Knots of the continuum editor -----------------------------------------------------

@pytest.fixture
def editing(window, modifiers):
    viewer = window.plot_viewer
    viewer.start_continuum_edit(initial_stride=500)
    return viewer


def test_drag_a_knot(editing, modifiers):
    viewer = editing
    canvas = viewer.canvas
    i = len(viewer._knots_x) // 2
    x0, y0 = viewer._knots_x[i], viewer._knots_y[i]
    assert viewer.on_press(_event(canvas, "button_press_event", x0, y0, button=1))
    assert viewer._active_knot_idx == i
    x1, y1 = x0 + 0.05, 1.3 * y0
    viewer.on_mouse_move(_event(canvas, "motion_notify_event", x1, y1))
    assert viewer._knots_x[i] == pytest.approx(x1) and viewer._knots_y[i] == pytest.approx(y1)
    assert viewer._is_anchor[i]
    np.testing.assert_array_equal(viewer.knot_artist.get_xdata(), viewer._knots_x)
    viewer.on_release(_event(canvas, "button_release_event", x1, y1, button=1))
    assert viewer._active_knot_idx is None
    # Moving without a knot grabbed changes nothing
    before = viewer._knots_y.copy()
    viewer.on_mouse_move(_event(canvas, "motion_notify_event", x1 + 1, y1))
    np.testing.assert_array_equal(viewer._knots_y, before)


def test_drag_pauses_the_tool(editing):
    viewer = editing
    viewer.toolbar.pan()
    i = len(viewer._knots_x) // 2
    event = _event(viewer.canvas, "button_press_event",
                   viewer._knots_x[i], viewer._knots_y[i], button=1)
    assert viewer.on_press(event)
    assert viewer.toolbar.mode == "" and viewer._paused_tool_mode == "pan/zoom"
    viewer.on_release(event)
    assert viewer.toolbar.mode == "pan/zoom" and viewer._paused_tool_mode is None
    viewer.toolbar.pan()


def test_add_a_knot_with_shift(editing, modifiers):
    viewer = editing
    n = len(viewer._knots_x)
    x_full = viewer.main_window.active_history.current_state.spec.x.value
    x_new = 0.5 * (viewer._knots_x[3] + viewer._knots_x[4])
    expected_y = np.interp(x_new, x_full, viewer._frozen_spline_y)
    modifiers[0] = Qt.KeyboardModifier.ShiftModifier
    assert viewer.on_press(_event(viewer.canvas, "button_press_event", x_new, 0.0, button=1))
    assert len(viewer._knots_x) == n + 1
    assert viewer._knots_x[4] == pytest.approx(x_new) and viewer._is_anchor[4]
    assert viewer._knots_y[4] == pytest.approx(expected_y)  # on the curve
    assert np.all(np.diff(viewer._knots_x) > 0)


def test_remove_a_knot(editing, menus):
    viewer = editing
    n = len(viewer._knots_x)
    x5 = viewer._knots_x[5]
    assert viewer.on_press(_event(viewer.canvas, "button_press_event",
                                  x5, viewer._knots_y[5], button=3))
    assert len(viewer._knots_x) == n - 1 and x5 not in viewer._knots_x
    assert menus == []  # no menu over a knot
    # Two knots are always kept
    viewer._knots_x, viewer._knots_y = viewer._knots_x[:2], viewer._knots_y[:2]
    viewer._is_anchor = viewer._is_anchor[:2]
    assert viewer.on_press(_event(viewer.canvas, "button_press_event",
                                  viewer._knots_x[0], viewer._knots_y[0], button=3))
    assert len(viewer._knots_x) == 2
