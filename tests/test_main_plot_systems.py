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
    # The ticks are converted with z_ref = 0 whatever the z_ref of the
    # conversion (F66): only z_ref = 0 puts them on the plot today
    state = window.active_history.current_state.edit.x_convert(
        z_ref="0", xunit="km/s")
    ticks = _ticks(window.plot_viewer, state, _selection(state))
    assert len(ticks["x"]) > 700 and len(ticks["selected"]) == 3
    snapshot.check("main_plot_system_ticks_kms", ticks)
