"""The ticks of the absorption systems on the main plot.

Characterization test written before the positions of the ticks were
computed in one conversion instead of one per transition (phase 5): the
positions (in nm and in km/s), the components they belong to and the
three groups drawn (HI, metals, selected) are compared with a snapshot.
"""
import numpy as np
import pytest
from PySide6.QtWidgets import QMessageBox

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


def _ticks(viewer, state, selected):
    viewer._highlighted_uuids = set(selected)
    viewer.plot_spectrum(state, initial_draw=True, force_autoscale=True)
    lines = {line.get_label(): np.asarray(line.get_xdata(), dtype=float)
             for line in viewer.canvas.axes.get_lines()
             if line.get_label() in ("HI Systems", "Metal Systems", "Selected")}
    return {
        "x": np.array([x for _, x, _ in viewer._tick_visuals]),
        "id": np.array([comp.id for _, _, comp in viewer._tick_visuals]),
        "hi": lines.get("HI Systems", np.array([])),
        "metal": lines.get("Metal Systems", np.array([])),
        "selected": lines.get("Selected", np.array([])),
    }


def _selection(state):
    """A CIV doublet and a Ly-alpha line."""
    comps = state.systs.components
    return [next(c for c in comps if c.series == "CIV").uuid,
            next(c for c in comps if c.series == "Ly_a").uuid]


def test_system_ticks_in_nm(window, snapshot):
    state = window.active_history.current_state
    ticks = _ticks(window.plot_viewer, state, _selection(state))
    assert len(ticks["x"]) > 700 and len(ticks["selected"]) == 3
    snapshot.check("main_plot_system_ticks_nm", ticks)


def test_system_ticks_in_velocity(window, snapshot):
    # The archive has no z_em: the reference of the ticks is 0, as the
    # z_ref of the conversion
    state = window.active_history.current_state.edit.x_convert(
        z_ref="0", xunit="km/s")
    ticks = _ticks(window.plot_viewer, state, _selection(state))
    assert len(ticks["x"]) > 700 and len(ticks["selected"]) == 3
    snapshot.check("main_plot_system_ticks_kms", ticks)


def test_system_ticks_in_velocity_with_z_em(window):
    """With z_em, the axis converted at z_ref = z_em (the default): the
    ticks are referred to z_em too, all on the plot (F66: they were
    referred to 0 and fell outside it)."""
    from astrocook.core.spectrum_operations import convert_axis_velocity
    import astropy.units as au
    state = window.active_history.current_state
    plain = _ticks(window.plot_viewer, state.edit.x_convert(z_ref="0", xunit="km/s"),
                   _selection(state))
    z_em = 4.072
    with_z_em = state.edit.set_properties(z_em=str(z_em))
    converted = with_z_em.edit.x_convert(xunit="km/s")  # z_ref = z_em
    ticks = _ticks(window.plot_viewer, converted, _selection(with_z_em))
    np.testing.assert_array_equal(ticks["id"], plain["id"])
    shift = convert_axis_velocity(np.array([1.0]) * au.nm, z_em, au.km / au.s).value[0] - \
        convert_axis_velocity(np.array([1.0]) * au.nm, 0.0, au.km / au.s).value[0]
    np.testing.assert_allclose(ticks["x"], plain["x"] + shift, rtol=0, atol=1e-6)
    assert len(ticks["selected"]) == 3
    lo, hi = window.plot_viewer.canvas.axes.get_xlim()
    assert lo <= ticks["x"].min() and ticks["x"].max() <= hi
