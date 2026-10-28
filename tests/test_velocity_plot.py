"""The velocity plot of the System Inspector: the tooltip of a tick, the
profile under the cursor and the menu of the right click (focus on a
component, add a component, a subset of a multiplet or a linked system).

Characterization tests written before the velocity plot was moved to its
own module (phase 5, split of the GUI). The menus are recorded instead of
shown (``QMenu.exec`` would block an offscreen test).
"""
import sys

import numpy as np
import pytest
from matplotlib.backend_bases import MouseEvent
from PySide6.QtWidgets import QMenu, QMessageBox

from astrocook.gui.main_window import MainWindowV2

from test_gui_smoke import SESSION_FILE


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
    return win


@pytest.fixture
def civ(qtbot, window):
    """The inspector on the first CIV component of the archive."""
    window._on_view_system_inspector()
    inspector = window.system_inspector
    comp = next(c for c in window.active_history.current_state.systs.components
                if c.series == "CIV")
    inspector.focus_on_component(comp.uuid)
    qtbot.wait(100)
    return inspector, comp


def _plot_module(inspector):
    """The module of the velocity plot, whose names (QMenu, QToolTip) the
    tests replace: found from the class, wherever it is defined."""
    return sys.modules[type(inspector.vel_plot).__module__]


@pytest.fixture
def menus(civ, monkeypatch):
    """The texts of the actions of each menu that would be shown, and the
    actions themselves, to trigger them."""
    shown = []

    class RecordedMenu(QMenu):
        def exec(self, *args):  # record the menu instead of showing it
            shown.append({a.text(): a for a in self.actions() if a.text()})

    monkeypatch.setattr(_plot_module(civ[0]), "QMenu", RecordedMenu)
    return shown


def _event(plot, name, ax, v, y=0.5, button=None):
    x_disp, y_disp = ax.transData.transform((v, y))
    return MouseEvent(name, plot.canvas, x_disp, y_disp, button=button)


def _free_v(plot, ax):
    """A velocity of the panel away from every tick (a click on the
    background, not on a component)."""
    lo, hi = plot._xlim
    tol = (hi - lo) * 0.015
    ticks = [v for a, v, _ in plot._tick_map_visuals if a is ax]
    for v in np.linspace(lo, hi, 200)[10:-10]:
        if all(abs(v - t) > 2 * tol for t in ticks):
            return float(v)
    raise AssertionError("no free velocity")


def _with_transitions(qtbot, inspector, text):
    inspector.trans_in.setText(text)
    inspector._apply()
    qtbot.wait(100)


# --- Hover ---------------------------------------------------------------------------------

def test_tooltip_of_a_tick(civ, monkeypatch):
    inspector, comp = civ
    plot = inspector.vel_plot
    tips = []

    class RecordedTip:
        showText = staticmethod(lambda *a, **k: tips.append(a[1]))
        hideText = staticmethod(lambda: None)

    monkeypatch.setattr(_plot_module(inspector), "QToolTip", RecordedTip)
    ax, v, tick_comp = plot._tick_map_visuals[0]
    plot.on_mouse_move(_event(plot, "motion_notify_event", ax, v))
    assert tips and f"<b>ID:</b> {tick_comp.id} | <b>{tick_comp.series}</b>" in tips[-1]
    assert f"<b>z:</b> {tick_comp.z:.5f}" in tips[-1]


def test_profile_under_the_cursor(civ):
    inspector, comp = civ
    plot = inspector.vel_plot
    ax = next(iter(plot._panel_map))
    plot.on_mouse_move(_event(plot, "motion_notify_event", ax, 37.0))
    for trans, line_main, _ in plot._panel_map.values():
        v, flux = line_main.get_data()
        assert len(v) == 300 and np.all(np.isfinite(flux)), trans
        assert np.min(flux) < 1.0  # the line of the cursor (logN 13.5, b 10)


# --- Right click ----------------------------------------------------------------------

def test_right_click_on_a_tick_and_on_the_background(civ, menus, monkeypatch):
    inspector, comp = civ
    plot = inspector.vel_plot
    window = inspector.main_window
    requests = []
    monkeypatch.setattr(window, "_on_recipe_requested",
                        lambda cat, name, params, extra: requests.append((cat, name, params)))
    ax, v_tick, tick_comp = plot._tick_map_visuals[0]
    plot.on_press(_event(plot, "button_press_event", ax, v_tick, button=3))
    assert f"Focus on Component {tick_comp.id}" in menus[-1]

    trans = plot._panel_map[ax][0]
    v_free = _free_v(plot, ax)
    plot.on_press(_event(plot, "button_press_event", ax, v_free, button=3))
    texts = list(menus[-1])
    z_new = (1 + plot._plot_center_z) * (1 + v_free / 299792.458) - 1
    assert texts == [f"Add {trans} (Single Line) at z={z_new:.5f}",
                     f"Add CIV (Standard) at z={z_new:.5f}"]
    menus[-1][texts[1]].trigger()
    cat, name, params = requests[-1]
    assert (cat, name) == ("absorbers", "add_component")
    assert params["series"] == "CIV" and params["z"] == pytest.approx(z_new)
    assert (params["logN"], params["b"]) == (13.5, 10.0)


def test_right_click_with_a_subset_of_a_large_multiplet(qtbot, civ, menus):
    # A FeII component of the archive: three lines of a large multiplet
    inspector, _ = civ
    feii = next(c for c in inspector.main_window.active_history.current_state.systs.components
                if c.series == "FeII_2344,FeII_2374,FeII_2382")
    inspector.focus_on_component(feii.uuid)
    qtbot.wait(100)
    plot = inspector.vel_plot
    assert plot._all_transitions == ["FeII_2344", "FeII_2374", "FeII_2382"]
    ax = next(iter(plot._panel_map))
    plot.on_press(_event(plot, "button_press_event", ax, _free_v(plot, ax), button=3))
    texts = list(menus[-1])
    assert len(texts) == 2 and texts[0].startswith("Add FeII_2344 (Single Line)")
    assert texts[1].startswith("Add FeII_2344,FeII_2374,FeII_2382 (Subset)")


def test_right_click_with_mixed_species(qtbot, civ, menus, monkeypatch):
    inspector, comp = civ
    window = inspector.main_window
    requests = []
    monkeypatch.setattr(window, "_on_recipe_requested",
                        lambda cat, name, params, extra: requests.append((cat, name, params)))
    _with_transitions(qtbot, inspector, "CIV,SiIV")
    plot = inspector.vel_plot
    ax = next(iter(plot._panel_map))
    plot.on_press(_event(plot, "button_press_event", ax, _free_v(plot, ax), button=3))
    linked = [t for t in menus[-1] if "(Linked)" in t]
    assert len(linked) == 1
    # The label is cut at 32 characters
    assert linked[0].startswith("Add CIV_1548,CIV_1550,SiIV_1393,SiIV... (Linked)")
    menus[-1][linked[0]].trigger()
    cat, name, params = requests[-1]
    assert name == "add_linked_system"
    assert params["series_list"] == "CIV_1548,CIV_1550,SiIV_1393,SiIV_1402"


def test_left_click_shows_no_menu(civ, menus):
    inspector, comp = civ
    plot = inspector.vel_plot
    ax = next(iter(plot._panel_map))
    plot.on_press(_event(plot, "button_press_event", ax, 20.0, button=1))
    assert menus == []


# --- Pixels ----------------------------------------------------------------------------------

def _pixels(qtbot, plot):
    qtbot.wait(100)
    plot.canvas.draw()
    return np.array(plot.canvas.buffer_rgba()).copy()


def test_the_plot_is_drawn_the_same(qtbot, civ):
    """The pixels of the velocity plot for a sequence of selections, with
    and without the residuals, compared with a full rebuild with the
    constrained layout of matplotlib (characterization, before the
    rebuild and the layout were made faster in phase 5)."""
    inspector, comp = civ
    plot = inspector.vel_plot
    for row in (0, 5, 6, 40, 5):
        for resid in (False, True):
            inspector.resid_cb.setChecked(resid)
            inspector.table_view.selectRow(row)
            shown = _pixels(qtbot, plot)
            # The reference: a new figure laid out from scratch
            getattr(plot, "_layout_cache", {}).clear()
            plot.fig.clear()
            plot.fig.set_layout_engine(layout='constrained')
            plot._update_plot()
            assert np.array_equal(shown, _pixels(qtbot, plot)), (row, resid)
