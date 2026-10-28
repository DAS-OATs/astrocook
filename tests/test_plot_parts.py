"""The parts of the main plot around SpectrumPlotWidget: the style helpers
and the decimation, the canvas (cursor, selection of a region, tooltip
after a pause, blitting), the toolbar and the data inspector.

Characterization tests written before these parts were moved out of
``pyside_plot.py`` (phase 5, split of the GUI).
"""
import numpy as np
import pytest
from matplotlib.backend_bases import MouseEvent
from PySide6.QtWidgets import QMessageBox

from astrocook.gui.main_window import MainWindowV2
from astrocook.gui import plot_style

from test_gui_smoke import SESSION_FILE


# --- Style and decimation -------------------------------------------------------------

def test_style_colours():
    cycle = ["c0", "c1", "c2"]
    for key, style in plot_style.PLOT_STYLE.items():
        colour = plot_style.get_style_color(key, cycle)
        if "color" in style:
            assert colour == style["color"]
        else:
            assert colour == cycle[style.get("color_idx", 0) % 3]


def test_z_from_x():
    assert plot_style.z_convert_inverse(486.268, 121.567) == pytest.approx(3.0)
    assert plot_style.z_convert_inverse("bad", 121.567) is None


def test_min_max_decimation():
    y = np.array([3.0, 1.0, 2.0, 5.0, 4.0, 0.0, 9.0])
    x = np.arange(7.0)
    np.testing.assert_array_equal(plot_style.decimate_y_min_max(y, 3),
                                  [1.0, 3.0, 0.0, 5.0])  # the last pixel left out
    np.testing.assert_array_equal(plot_style.decimate_x_for_min_max(x, 3),
                                  [0.0, 2.0, 3.0, 5.0])
    # Too short, or nothing: unchanged
    assert plot_style.decimate_y_min_max(y[:2], 3) is not None
    np.testing.assert_array_equal(plot_style.decimate_y_min_max(y[:2], 3), y[:2])
    np.testing.assert_array_equal(plot_style.decimate_x_for_min_max(x[:2], 3), x[:2])
    assert plot_style.decimate_y_min_max(None, 3) is None


# --- The main window ---------------------------------------------------------------------

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


def _event(canvas, name, x, y=1.0, button=None):
    x_disp, y_disp = canvas.axes.transData.transform((x, y))
    return MouseEvent(name, canvas, x_disp, y_disp, button=button)


def _x_inside(canvas, fraction):
    lo, hi = canvas.axes.get_xlim()
    return lo + fraction * (hi - lo)


def _y_inside(canvas):
    lo, hi = canvas.axes.get_ylim()
    return 0.5 * (lo + hi)


# --- Canvas ------------------------------------------------------------------------------

def test_cursor_follows_the_mouse(qtbot, window):
    canvas = window.plot_viewer.canvas
    window.cursor_z_input.setText("3.0")
    window.cursor_show_checkbox.setChecked(True)
    qtbot.wait(100)
    assert canvas.cursor_artists  # the lines of the cursor
    x = _x_inside(canvas, 0.5)
    canvas.on_motion(_event(canvas, "motion_notify_event", x, _y_inside(canvas)))
    z = float(window.cursor_z_input.text())
    assert z == pytest.approx(window.plot_viewer.calculate_z_from_x(x), abs=1e-7)
    positions = window.plot_viewer.get_cursor_line_positions_at_z(
        window.active_history.current_state, z)
    assert [a.get_xdata()[0] for a in canvas.cursor_artists] == pytest.approx(positions)


def test_region_selection_is_not_reachable(window):
    # The selection of a region with the mouse (toggle_select_mode ->
    # toggle_region_selector, then on_press/on_motion/on_release) is dead
    # code: the toolbar has no 'Select' action, nothing else calls it, and
    # its last step calls main_window.launch_split_from_region, which does
    # not exist (findings register F73). Splitting by the current view is in
    # the context menu (launch_split_from_current_view)
    viewer = window.plot_viewer
    assert "Select" not in [a.text() for a in viewer.toolbar.actions()]
    assert not hasattr(window, "launch_split_from_region")
    assert hasattr(window, "launch_split_from_current_view")
    assert not viewer._is_selecting_region


def test_selection_span_is_drawn_and_removed(qtbot, window):
    # The canvas part of the (unreachable) selection: the span follows the
    # mouse and is removed by cleanup_selection
    viewer = window.plot_viewer
    canvas = viewer.canvas
    viewer.toggle_region_selector()
    y = _y_inside(canvas)
    canvas.selection_start_x = _x_inside(canvas, 0.6)
    canvas.on_motion(_event(canvas, "motion_notify_event", _x_inside(canvas, 0.3), y))
    assert canvas.selection_artist is not None
    canvas.cleanup_selection()
    assert canvas.selection_artist is None and canvas.selection_start_x is None
    viewer.toggle_region_selector()
    assert not viewer._is_selecting_region


def test_tooltip_after_a_pause(window, monkeypatch):
    viewer = window.plot_viewer
    canvas = viewer.canvas
    shown = []
    monkeypatch.setattr(viewer, "show_data_tooltip", shown.append)
    event = _event(canvas, "motion_notify_event", _x_inside(canvas, 0.5), _y_inside(canvas))
    canvas._last_mouse_event = event
    canvas._on_hover_timeout()
    assert shown == [event]


def test_blitting_needs_a_background(qtbot, window, monkeypatch):
    canvas = window.plot_viewer.canvas
    canvas.background = None
    idle = []
    monkeypatch.setattr(canvas, "draw_idle", lambda: idle.append(True))
    canvas.draw_animated()
    assert idle == [True]


def test_changing_the_limits_replots(qtbot, window, monkeypatch):
    viewer = window.plot_viewer
    canvas = viewer.canvas
    qtbot.wait(100)
    assert canvas.background is not None
    calls = []
    monkeypatch.setattr(viewer, "plot_spectrum", lambda **k: calls.append(k))
    lo, hi = canvas.axes.get_xlim()
    canvas.axes.set_xlim(lo, 0.5 * (lo + hi))
    assert canvas.background is None
    qtbot.waitUntil(lambda: bool(calls))


# --- Toolbar ------------------------------------------------------------------------------

def test_toolbar(window, monkeypatch):
    viewer = window.plot_viewer
    texts = [a.text() for a in viewer.toolbar.actions() if a.text()]
    for kept in ("Home", "Back", "Forward", "Pan", "Zoom", "Save"):
        assert kept in texts
    assert "Customize" not in texts and "Subplots" not in texts
    calls = []
    monkeypatch.setattr(viewer, "plot_spectrum", lambda **k: calls.append(k))
    viewer.toolbar.home()
    assert calls == [{"session_state": window.active_history.current_state,
                      "force_autoscale": True}]


# --- Data inspector -----------------------------------------------------------------------

def test_data_inspector(qtbot, window):
    viewer = window.plot_viewer
    session = window.active_history.current_state
    viewer._open_data_inspector(session, 1000)
    dialog = viewer._data_inspector_dialog
    table = dialog.table
    assert dialog.windowTitle().startswith("Data Inspector (x = ")
    headers = [table.horizontalHeaderItem(c).text() for c in range(table.columnCount())]
    assert headers[:8] == ["Index", "λ", "λmin", "λmax", "F", "dF", "cont", "model"]
    assert table.rowCount() == 31  # window of 15 on each side
    center = 1000 - dialog.current_start
    assert table.item(center, 0).text() == "1000"
    assert table.item(center, 0).font().bold()  # the clicked pixel highlighted
    # Scrolling to the bottom loads the next chunk
    vbar = table.verticalScrollBar()
    vbar.setValue(vbar.maximum())
    assert dialog.current_end == 1016 + dialog.chunk_size
    assert table.rowCount() == 31 + dialog.chunk_size
    # ... and to the top the previous one
    vbar.setValue(vbar.minimum())
    assert dialog.current_start == 985 - dialog.chunk_size
    # A new state keeps the window and refreshes the columns
    dialog.refresh(session)
    assert table.rowCount() == dialog.current_end - dialog.current_start
    dialog.close()
    qtbot.waitUntil(lambda: viewer._data_inspector_dialog is None)


def test_data_inspector_keeps_a_window_of_rows(qtbot, window):
    viewer = window.plot_viewer
    viewer._open_data_inspector(window.active_history.current_state, 1000)
    dialog = viewer._data_inspector_dialog
    vbar = dialog.table.verticalScrollBar()
    for _ in range(12):
        vbar.setValue(vbar.maximum())
    assert dialog.table.rowCount() <= dialog.max_rows
    assert dialog.table.rowCount() == dialog.current_end - dialog.current_start
    dialog.close()
