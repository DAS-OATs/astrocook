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


@pytest.mark.parametrize("key, label", [
    ("Ly_a", "Lyα"), ("Ly_b", "Lyβ"), ("Ly_6", "Ly6"), ("Ly_lim", "Ly limit"),
    ("H_a", "Hα"), ("D_b", "D Lyβ"), ("CIV_1548", "C IV"), ("CIII_1908", "C III]"),
    ("OIII_5008", "[O III]"), ("CII_1334", "C II 1334"), ("CII*_1335", "C II* 1335"),
    ("NI_1134a", "N I 1134a"), ("SiIV_1402", "Si IV 1402"), ("unknown", "unknown")])
def test_emission_labels(key, label):
    # The names of the emission lines in spectroscopic notation (visual
    # review of phase 5-bis: the plot showed the keys, e.g. NV_1238)
    assert plot_style.emission_label(key) == label


def test_every_line_has_a_label():
    from astrocook.core.atomic_data import xem_d
    assert all(plot_style.emission_label(k) != k for k in xem_d)


def test_flux_unit_labels():
    import astropy.units as au
    flux = au.erg / (au.cm**2 * au.s * au.Angstrom)
    assert plot_style.flux_unit_label(flux) == "erg s⁻¹ cm⁻² Å⁻¹"
    assert plot_style.flux_unit_label(au.erg / (au.cm**2 * au.s * au.Hz)) == "erg s⁻¹ cm⁻² Hz⁻¹"
    assert plot_style.flux_unit_label(au.Unit("1e-17 erg / (s cm2 Angstrom)")) == "1×10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹"
    assert plot_style.flux_unit_label(au.Jy) == "Jy"


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


def test_floating_toolbar(qtbot, window):
    # The toolbar floats over the plot as the side cards (phase 5-bis,
    # Guido's request): flush with the bottom left corner of the plot area,
    # fading, with a tab at its right edge that shows or hides it, as View >
    # Show Toolbar does
    from PySide6.QtCore import QSettings
    viewer = window.plot_viewer
    toolbar, tab = viewer.toolbar, viewer.toolbar.tab
    assert viewer.main_layout.indexOf(toolbar) == -1
    # As wide as the side cards, and all its tools fit (Guido)
    from astrocook.gui.floating_cards import CARD_WIDTH
    assert toolbar.width() == CARD_WIDTH
    assert toolbar.sizeHint().width() <= CARD_WIDTH

    def in_the_corner():
        g, t = toolbar.geometry(), tab.geometry()
        bottom = viewer.rect().bottom()
        return (g.left() == 0 and g.bottom() == bottom
                and t.left() == g.right() + 1 and t.bottom() == bottom)
    assert in_the_corner()
    window.resize(window.width() - 100, window.height() - 80)
    qtbot.wait(100)
    assert in_the_corner()
    # It fades as the cards, with its tab
    toolbar._fader.set_faded(True, animate=False)
    assert toolbar._fader.effect.opacity() < 0.5
    toolbar._fader.set_faded(False, animate=False)
    # The tab hides it (and stays in the corner), as View > Show Toolbar
    tab.click()
    assert not toolbar.isVisible() and tab.isVisible()
    assert tab.geometry().left() == 0
    assert not window.show_toolbar_action.isChecked()
    assert QSettings("Astrocook", "AstrocookV2").value("show_toolbar") in (False, "false")
    window.show_toolbar_action.trigger()
    assert toolbar.isVisible() and tab.geometry().left() == toolbar.width()
    tab.click()
    tab.click()
    assert toolbar.isVisible() and window.show_toolbar_action.isChecked()


def test_session_card_stops_above_the_toolbar(qtbot, window, monkeypatch):
    # With many sessions the card of the sessions reached the floating
    # toolbar in the bottom left corner and covered it (Guido): it stops
    # above the toolbar, or above its tab when the toolbar is hidden
    toolbar = window.plot_viewer.toolbar
    monkeypatch.setattr(window, "_session_list_height", lambda: 5000)
    window._reposition_floating_widgets()

    def top_of(widget):
        return widget.mapTo(window, widget.rect().topLeft()).y()
    card = window.left_sidebar_widget
    assert card.geometry().bottom() < top_of(toolbar)
    window.show_toolbar_action.trigger()  # hidden: only the tab
    assert card.geometry().bottom() < top_of(toolbar.tab)
    assert card.height() > 0
    window.show_toolbar_action.trigger()
    assert card.geometry().bottom() < top_of(toolbar)


def test_inspector_toolbar_remembered(qtbot, window):
    # The toolbar of the velocity plot, shown or hidden by its tab, as the
    # user left it (Guido)
    from PySide6.QtCore import QSettings
    from astrocook.gui.velocity_plot import INSPECTOR_TOOLBAR_KEY
    window._on_view_system_inspector()
    toolbar = window.system_inspector.vel_plot.toolbar
    assert toolbar.isVisible()
    toolbar.tab.click()
    assert not toolbar.isVisible()
    assert QSettings("Astrocook", "AstrocookV2").value(INSPECTOR_TOOLBAR_KEY) in (False, "false")
    from astrocook.gui.system_inspector import SystemInspector
    other = SystemInspector(window)
    qtbot.addWidget(other)
    assert other.vel_plot.toolbar.isHidden()
    toolbar.tab.click()
    assert QSettings("Astrocook", "AstrocookV2").value(INSPECTOR_TOOLBAR_KEY) in (True, "true")


def test_status_bars_look_the_same(qtbot, window):
    # The status line of the inspector looked different from the status bar
    # of the window (Guido): one style, and the full width of the window
    from astrocook.gui.theme import status_bar_style
    window._on_view_system_inspector()
    inspector = window.system_inspector
    inspector.resize(1200, 700)
    qtbot.wait(50)
    assert status_bar_style() in window.styleSheet()
    assert inspector.status_bar.styleSheet() == status_bar_style()
    bar = inspector.status_bar
    assert bar.mapTo(inspector, bar.rect().topLeft()).x() == 0
    assert bar.width() == inspector.width()


def test_axes_on_whole_pixels(qtbot, window):
    # The constrained layout moved the axes by ~1e-13 px between draws,
    # which changed the simplified path of the spectrum: the blitted
    # highlights did not match a full draw at some window heights
    canvas = window.plot_viewer.canvas
    for dh in (0, 7, 13):
        window.resize(window.width(), window.height() - dh)
        qtbot.wait(50)
        canvas.draw()
        bounds = canvas.axes.bbox.bounds
        assert all(abs(v - round(v)) < 1e-9 for v in bounds), bounds


def test_cursor_after_a_resize(qtbot, window):
    # The pixels kept for blitting the cursor were those before the resize:
    # moving the cursor after View > Show Toolbar drew the old plot over the
    # new one (reported by Guido in phase 5-bis). The toolbar now floats
    # over the plot: a resize of the window instead
    canvas = window.plot_viewer.canvas
    window.cursor_z_input.setText("3.0")
    window.cursor_show_checkbox.setChecked(True)
    qtbot.wait(100)
    before = canvas.background.get_extents()
    window.resize(window.width(), window.height() - 60)
    qtbot.wait(100)
    canvas.on_motion(_event(canvas, "motion_notify_event", _x_inside(canvas, 0.5),
                            _y_inside(canvas)))
    box = canvas.blit_box()
    extents = canvas.background.get_extents()
    assert extents != before
    assert extents[3] - extents[1] == pytest.approx(box.height, abs=2)


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
    # 22 px floating over the plot, to fit the width of the cards (phase 5-bis)
    assert toolbar.iconSize().width() == 22
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
