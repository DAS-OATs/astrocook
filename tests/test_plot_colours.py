"""View > Plot Colours (phase 5-bis): the standard palette of the plots and
the colour-blind safe one, in the main plot, in the velocity plot and in
the cues of the flags of the System Inspector, saved in the settings."""
import matplotlib.colors as mcolors
import pytest
from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import QMessageBox

from astrocook.gui import plot_style
from astrocook.gui.inspector_table import flag_cue
from astrocook.gui.main_window import MainWindow

from test_gui_smoke import SESSION_FILE
from test_inspector_flags import _cell, _inspector


def _hex(colour):
    return mcolors.to_hex(colour)


@pytest.fixture
def window(qtbot, monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    for kind in ("warning", "critical", "information", "question"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(
            lambda *a, **k: QMessageBox.StandardButton.Ok))
    win = MainWindow(None, None)
    qtbot.addWidget(win)
    monkeypatch.setattr(win, "_show_custom_message", lambda *a, **k: None)
    win.show()
    win.open_session_from_path(str(SESSION_FILE))
    yield win
    # The palette is global and saved: back to the standard one for the
    # other tests of the worker
    win._set_plot_colours("standard")


def _main_lines(window):
    return {line.get_label(): line for line in window.plot_viewer.canvas.axes.get_lines()}


def test_palettes_have_the_same_keys():
    keys = [set(p) for p in plot_style.PLOT_PALETTES.values()]
    assert all(k == keys[0] for k in keys)
    assert set(plot_style.PLOT_PALETTE_NAMES) == set(plot_style.PLOT_PALETTES)


def test_standard_palette_is_the_old_one_with_a_thinner_model():
    # The colours of 2.2.0: the first and second of the tab20 cycle, the
    # cursor its fourth, the HI ticks and the flags as before
    standard = plot_style.PLOT_PALETTES["standard"]
    cycle = [_hex(c) for c in plot_style.get_color_cycle(5, cmap="tab20")]
    assert standard["flux"] == cycle[0] and standard["model"] == cycle[1]
    assert standard["cursor"] == cycle[3]
    assert standard["hi_ticks"] == plot_style.HI_TICK_COLOR
    from astrocook.gui import theme
    assert standard["flag_bound"] == theme.FLAG_BOUND_COLOR
    assert standard["flag_nan"] == theme.FLAG_NAN_COLOR
    assert standard["model_lw"] < 1.6


def test_unknown_palette_is_standard():
    try:
        assert plot_style.set_plot_palette("no such palette") == "standard"
        assert plot_style.plot_colour("model") == plot_style.PLOT_PALETTES["standard"]["model"]
    finally:
        plot_style.set_plot_palette("standard")


def test_menu_switches_the_main_plot_and_is_saved(window):
    actions = window.plot_colours_actions
    assert list(actions) == ["standard", "colorblind"]
    assert actions["standard"].isChecked()
    lines = _main_lines(window)
    assert _hex(lines["Spectrum"].get_color()) == "#1f77b4"
    actions["colorblind"].trigger()
    assert actions["colorblind"].isChecked()
    assert QSettings("Astrocook", "AstrocookV2").value("plot_colours") == "colorblind"
    lines = _main_lines(window)
    palette = plot_style.PLOT_PALETTES["colorblind"]
    assert _hex(lines["Spectrum"].get_color()) == palette["flux"]
    model = lines.get("Absorption model")
    if model is not None:
        assert _hex(model.get_color()) == palette["model"]
        assert model.get_linewidth() == palette["model_lw"]
    # A new window starts with the palette the user chose
    other = MainWindow(None, None)
    try:
        assert other.plot_colours_actions["colorblind"].isChecked()
    finally:
        other.close()
    actions["standard"].trigger()
    assert _hex(_main_lines(window)["Spectrum"].get_color()) == "#1f77b4"


def test_flag_cues_follow_the_palette(qtbot):
    inspector = _inspector(with_spectrum=True)
    qtbot.addWidget(inspector)
    try:
        plot_style.set_plot_palette("colorblind")
        palette = plot_style.PLOT_PALETTES["colorblind"]
        inspector.refresh_plot_colours()
        assert _cell(inspector, 0, "b", Qt.BackgroundRole).name() == palette["flag_bound"]
        assert _cell(inspector, 1, "ID", Qt.BackgroundRole).name() == palette["flag_nan"]
        comps = inspector.table_model._components
        assert flag_cue(comps[0], "b")[0] == palette["flag_bound"]
    finally:
        plot_style.set_plot_palette("standard")


def test_velocity_plot_follows_the_palette(qtbot, window, monkeypatch):
    window._on_view_system_inspector()
    inspector = window.system_inspector
    # The archive of the tests has no continuum, without which the plot is
    # replaced by a message: as with one
    monkeypatch.setattr(inspector, "has_continuum", lambda: True)
    comp = next(c for c in window.active_history.current_state.systs.components
                if c.series == "CIV")
    inspector.resid_cb.setChecked(True)
    inspector.focus_on_component(comp.uuid)
    qtbot.wait(100)
    plot = inspector.vel_plot
    for name in ("colorblind", "standard"):
        window.plot_colours_actions[name].trigger()
        palette = plot_style.PLOT_PALETTES[name]
        main, resid = plot.axes[0], plot.axes[1]
        colours = {_hex(line.get_color()) for line in main.get_lines()}
        assert palette["flux"] in colours
        # The archive has no model: the residual axis shows its zero line
        # and the line of the cursor profile, in the colour of the model
        resid_colours = {_hex(line.get_color()) for line in resid.get_lines()}
        assert resid_colours == {palette["model"]}
