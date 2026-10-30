"""The parts of the main plot around SpectrumPlotWidget: the style helpers
and the decimation, the canvas (cursor, selection of a region, tooltip
after a pause, blitting), the toolbar and the data inspector.

Characterization tests written before these parts were moved out of
``pyside_plot.py`` (phase 5, split of the GUI).
"""
import numpy as np
import pytest
from matplotlib.backend_bases import MouseEvent
from PySide6.QtWidgets import QMenu, QMessageBox

from astrocook.gui.main_window import MainWindowV2
from astrocook.gui import plot_mouse, plot_style

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


def test_cursor_label(qtbot, window):
    # The redshift of the cursor at the top of its first line in view: the
    # box of the card fades while the mouse is on the plot (F92)
    canvas = window.plot_viewer.canvas
    window.cursor_z_input.setText("3.0")
    window.cursor_show_checkbox.setChecked(True)
    qtbot.wait(100)
    label = canvas.cursor_label
    assert label is not None and label.get_animated()
    assert label.get_text().strip() == "z = 3.00000"
    x = _x_inside(canvas, 0.5)
    canvas.on_motion(_event(canvas, "motion_notify_event", x, _y_inside(canvas)))
    z = float(window.cursor_z_input.text())
    assert label.get_text().strip() == f"z = {z:.5f}"
    lo, hi = canvas.axes.get_xlim()
    in_view = [a.get_xdata()[0] for a in canvas.cursor_artists
               if lo <= a.get_xdata()[0] <= hi]
    assert label.get_visible() and label.get_position()[0] == pytest.approx(in_view[0])
    # On the side of the line away from the card of the plot controls
    canvas.update_cursor_label(z, [lo + 0.2 * (hi - lo)])
    assert label.get_horizontalalignment() == "left"
    canvas.update_cursor_label(z, [lo + 0.8 * (hi - lo)])
    assert label.get_horizontalalignment() == "right"
    # No line in view: no label
    canvas.update_cursor_label(z, [hi + 1.0])
    assert not label.get_visible()
    window.cursor_show_checkbox.setChecked(False)
    qtbot.wait(100)
    assert canvas.cursor_label is None


def test_cursor_after_a_resize(qtbot, window):
    # The pixels kept for blitting the cursor were those before the resize:
    # moving the cursor after View > Show Toolbar drew the old plot over the
    # new one (reported by Guido in phase 5-bis)
    canvas = window.plot_viewer.canvas
    window.cursor_z_input.setText("3.0")
    window.cursor_show_checkbox.setChecked(True)
    qtbot.wait(100)
    before = canvas.background.get_extents()
    window.show_toolbar_action.trigger()
    qtbot.wait(100)
    canvas.on_motion(_event(canvas, "motion_notify_event", _x_inside(canvas, 0.5),
                            _y_inside(canvas)))
    bbox = canvas.axes.bbox
    extents = canvas.background.get_extents()
    assert extents != before
    assert extents[3] - extents[1] == pytest.approx(bbox.height, abs=2)
    window.show_toolbar_action.trigger()


def test_select_tool(window):
    # The Select tool (F73: the selection of a region was dead code, with no
    # button): in the main toolbar after Zoom, exclusive with Pan and Zoom
    viewer = window.plot_viewer
    toolbar = viewer.toolbar
    texts = [a.text() for a in toolbar.actions()]
    assert texts[texts.index("Zoom") + 1] == "Select"
    toolbar.pan()
    toolbar.select_action.trigger()
    assert viewer._is_selecting_region and toolbar.select_action.isChecked()
    assert toolbar.mode == ""
    toolbar.zoom()
    assert not viewer._is_selecting_region and not toolbar.select_action.isChecked()
    assert toolbar.mode == "zoom rect"
    toolbar.select_action.trigger()
    assert viewer._is_selecting_region and toolbar.mode == ""
    toolbar.select_action.trigger()
    assert not viewer._is_selecting_region
    # Not in the other plots
    from astrocook.gui.plot_toolbar import AstrocookToolbar
    other = AstrocookToolbar(viewer.canvas, viewer)
    assert other.select_action is None and "Select" not in [a.text() for a in other.actions()]


def _drag(viewer, x0, x1):
    canvas = viewer.canvas
    y = _y_inside(canvas)
    viewer.on_press(_event(canvas, "button_press_event", x0, y, button=1))
    canvas.on_motion(_event(canvas, "motion_notify_event", 0.5 * (x0 + x1), y))
    viewer.on_release(_event(canvas, "button_release_event", x1, y, button=1))


@pytest.mark.parametrize("choice, recipe, params", [
    ("Split region", ("edit", "split"), "expression"),
    ("Mask region...", ("edit", "mask_expression"), "expression"),
    ("Rebin region...", ("flux", "rebin"), "range"),
    ("Zoom to region", None, None)])
def test_region_menu(window, monkeypatch, choice, recipe, params):
    viewer = window.plot_viewer
    canvas = viewer.canvas
    launched, menus = [], []
    monkeypatch.setattr(window, "_launch_recipe_dialog",
                        lambda cat, name, initial_params=None: launched.append(
                            ((cat, name), initial_params)))

    class RecordedMenu(QMenu):
        def exec(self, *args):
            # The span is drawn while the menu is open
            assert canvas.selection_artist is not None
            actions = {a.text(): a for a in self.actions() if not a.isSeparator()}
            menus.append(list(actions))
            actions[choice].trigger()
    monkeypatch.setattr(plot_mouse, "QMenu", RecordedMenu)

    viewer.set_region_selection(True)
    x0, x1 = _x_inside(canvas, 0.6), _x_inside(canvas, 0.3)  # dragged leftwards
    _drag(viewer, x0, x1)
    assert menus == [["Split region", "Mask region...", "Rebin region...", "Zoom to region"]]
    lo, hi = x1, x0
    if recipe is None:
        assert canvas.axes.get_xlim() == pytest.approx((lo, hi))
        assert launched == []
    elif params == "expression":
        assert launched == [(recipe, {"expression": f"(λ > {lo:.4f}) & (λ < {hi:.4f})"})]
    else:
        assert launched == [(recipe, {"xstart": f"{lo:.4f}", "xend": f"{hi:.4f}"})]
    # The span is gone after the menu; the tool stays on
    assert canvas.selection_artist is None and viewer._is_selecting_region


def test_a_click_is_not_a_region(window, monkeypatch):
    menus = []

    class RecordedMenu(QMenu):
        def exec(self, *args):
            menus.append(self)
    monkeypatch.setattr(plot_mouse, "QMenu", RecordedMenu)
    viewer = window.plot_viewer
    viewer.set_region_selection(True)
    x = _x_inside(viewer.canvas, 0.5)
    _drag(viewer, x, x)
    assert menus == [] and viewer.canvas.selection_artist is None
    # Not in selection mode, a drag selects nothing
    viewer.set_region_selection(False)
    _drag(viewer, _x_inside(viewer.canvas, 0.2), _x_inside(viewer.canvas, 0.4))
    assert menus == []


def test_select_needs_a_wavelength_axis(window):
    import astropy.units as au
    viewer = window.plot_viewer
    viewer.set_region_selection(True)
    state = window.active_history.current_state
    in_velocity = state.with_new_spectrum(state.spec.x_convert(z_ref=3.0, xunit=au.km / au.s))
    viewer.plot_spectrum(in_velocity)
    assert not viewer.toolbar.select_action.isEnabled() and not viewer._is_selecting_region
    viewer.plot_spectrum(state)
    assert viewer.toolbar.select_action.isEnabled()


def test_selection_span_is_drawn_and_removed(qtbot, window):
    # The span follows the mouse and is removed by cleanup_selection
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


# --- The toolbar (redesigned in phase 5-bis) --------------------------------------

def test_toolbar_icons_and_tips(window):
    from astrocook.gui.plot_toolbar import TOOL_ICONS
    toolbar = window.plot_viewer.toolbar
    assert set(toolbar._tool_actions) == set(TOOL_ICONS)
    for key, action in toolbar._tool_actions.items():
        assert not action.icon().isNull(), key
    tips = {k: a.toolTip() for k, a in toolbar._tool_actions.items()}
    assert tips["pan"].endswith("(P)") and tips["zoom"].endswith("(Z)")
    assert tips["select"].endswith("(X)") and tips["home"].endswith("(H)")
    assert toolbar.iconSize().width() == 24
    # The coordinates under the mouse are still written by the toolbar
    assert toolbar.locLabel is not None
    # The other plots get the same look, without Select
    from astrocook.gui.plot_toolbar import AstrocookToolbar
    other = AstrocookToolbar(window.plot_viewer.canvas, None)
    assert "select" not in other._tool_actions and not other._tool_actions["pan"].icon().isNull()


def test_toolbar_shortcuts(qtbot, window):
    from PySide6.QtCore import Qt
    viewer = window.plot_viewer
    toolbar = viewer.toolbar
    window.activateWindow()
    viewer.canvas.setFocus()
    qtbot.waitUntil(lambda: viewer.canvas.hasFocus())
    qtbot.keyClick(viewer.canvas, Qt.Key.Key_P)
    assert toolbar.mode == "pan/zoom"
    qtbot.keyClick(viewer.canvas, Qt.Key.Key_Z)
    assert toolbar.mode == "zoom rect"
    qtbot.keyClick(viewer.canvas, Qt.Key.Key_X)
    assert viewer._is_selecting_region and toolbar.mode == ""
    qtbot.keyClick(viewer.canvas, Qt.Key.Key_Escape)
    assert not viewer._is_selecting_region and toolbar.mode == ""
    # Not while typing in a field of the window
    window.cursor_z_input.setFocus()
    qtbot.waitUntil(lambda: window.cursor_z_input.hasFocus())
    qtbot.keyClick(window.cursor_z_input, Qt.Key.Key_P)
    assert toolbar.mode == ""


def test_a_click_gives_the_plot_the_focus(qtbot, window):
    canvas = window.plot_viewer.canvas
    window.activateWindow()
    window.cursor_z_input.setFocus()
    canvas.on_press(_event(canvas, "button_press_event", _x_inside(canvas, 0.5),
                           _y_inside(canvas), button=1))
    qtbot.waitUntil(lambda: canvas.hasFocus())
