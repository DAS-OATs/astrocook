"""The continuum editor of the main plot: the knots spaced in velocity, the
draft curve frozen on the continuum, a knot moved and the curve patched
around it, the stride of the knots, reset, save and cancel.

Characterization tests written before the editor was moved to a mixin
(phase 5, split of the GUI): no test reached it.
"""
import numpy as np
import pytest
from PySide6.QtWidgets import QMessageBox
from scipy.interpolate import PchipInterpolator

from astrocook.gui.main_window import MainWindowV2

from test_gui_smoke import SESSION_FILE, _run_recipe

C_KMS = 299792.458


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


@pytest.fixture
def editing(window):
    """The editor started on the archive (no continuum: the knots on the
    flux), with knots every 500 km/s."""
    viewer = window.plot_viewer
    viewer.start_continuum_edit(initial_stride=500)
    return viewer


def _spec(viewer):
    return viewer.main_window.active_history.current_state.spec


# --- Knots -----------------------------------------------------------------------------

def test_knots_spaced_in_velocity(window):
    viewer = window.plot_viewer
    x = 500.0 * np.exp(np.arange(1000) * 2.0 / C_KMS)  # 2 km/s pixels
    idx = viewer._generate_velocity_knots(x, 100.0)
    v = C_KMS * np.log(x / x[0])
    targets = np.arange(0, v[-1], 100.0)
    # Each knot is the first pixel at or beyond its velocity (k * 100 km/s)
    assert len(idx) == len(targets) + 1
    assert np.all(v[idx[:-1]] >= targets) and np.all(v[idx[1:-1] - 1] < targets[1:])
    assert idx[-1] == 999  # the last pixel always included
    # A stride below the pixel: every pixel once
    np.testing.assert_array_equal(viewer._generate_velocity_knots(x, 0.5), np.arange(1000))


def test_start(editing):
    viewer = editing
    spec = _spec(viewer)
    x, y = spec.x.value, spec.y.value
    idx = viewer._generate_velocity_knots(x, 500.0)
    np.testing.assert_array_equal(viewer._knots_x, x[idx])
    np.testing.assert_array_equal(viewer._knots_y, y[idx])  # no continuum: the flux
    assert not viewer._is_anchor.any()
    np.testing.assert_array_equal(viewer._frozen_spline_y, y)
    assert viewer._edit_mode_active
    assert viewer.knot_artist.get_animated() and viewer.spline_artist.get_animated()
    assert viewer.background is not None  # captured at the draw


def test_moved_knot_patches_its_neighbourhood(editing):
    viewer = editing
    x_full = _spec(viewer).x.value
    before = viewer._frozen_spline_y.copy()
    i = len(viewer._knots_x) // 2
    viewer._knots_y[i] *= 1.5
    viewer._apply_local_patch(viewer._knots_x[i])
    xk, yk = viewer._get_sanitized_knots()
    local = PchipInterpolator(xk[i - 1:i + 2], yk[i - 1:i + 2])
    inside = (x_full >= xk[i - 1]) & (x_full <= xk[i + 1])
    np.testing.assert_allclose(viewer._frozen_spline_y[inside], local(x_full[inside]))
    np.testing.assert_array_equal(viewer._frozen_spline_y[~inside], before[~inside])
    np.testing.assert_array_equal(viewer.spline_artist.get_ydata(), viewer._frozen_spline_y)


def test_stride_keeps_the_anchors_and_the_curve(editing):
    viewer = editing
    x_full = _spec(viewer).x.value
    viewer._is_anchor[3] = True
    anchor = viewer._knots_x[3]
    frozen = viewer._frozen_spline_y.copy()
    viewer.begin_stride_interaction()
    assert viewer._stride_backup_x is not None
    viewer.update_continuum_stride(2000)
    viewer.end_stride_interaction()
    assert viewer._stride_backup_x is None
    floating = x_full[viewer._generate_velocity_knots(x_full, 2000.0)]
    np.testing.assert_array_equal(viewer._knots_x, np.unique(np.r_[anchor, floating]))
    np.testing.assert_array_equal(viewer._knots_y,
                                  np.interp(viewer._knots_x, x_full, frozen))
    assert viewer._is_anchor.sum() == 1 and viewer._knots_x[viewer._is_anchor][0] == anchor
    np.testing.assert_array_equal(viewer._frozen_spline_y, frozen)  # unchanged


def test_reset(editing):
    viewer = editing
    y = _spec(viewer).y.value
    viewer._knots_y[2] *= 2.0
    viewer._frozen_spline_y[:] = 0.0
    viewer.reset_continuum_to_original(1000)
    x = _spec(viewer).x.value
    np.testing.assert_array_equal(viewer._knots_x, x[viewer._generate_velocity_knots(x, 1000.0)])
    np.testing.assert_array_equal(viewer._frozen_spline_y, y)


def test_refresh_edit_layer(editing):
    viewer = editing
    viewer._knots_y = viewer._knots_y * 1.1
    viewer.refresh_edit_layer()
    np.testing.assert_array_equal(viewer.knot_artist.get_ydata(), viewer._knots_y)


def test_save(editing):
    viewer = editing
    viewer._frozen_spline_y = viewer._frozen_spline_y * 1.01
    expected = viewer._frozen_spline_y.tolist()
    x, y = viewer.stop_continuum_edit(save=True)
    assert x == _spec(viewer).x.value.tolist() and y == expected  # every pixel
    assert not viewer._edit_mode_active
    assert viewer.knot_artist is None and viewer.spline_artist is None


def test_stop_without_saving(editing):
    assert editing.stop_continuum_edit(save=False) == (None, None)


def test_cancel(editing, monkeypatch):
    viewer = editing
    replots = []
    monkeypatch.setattr(viewer, "plot_spectrum", lambda **k: replots.append(k))
    viewer.cancel_continuum_edit()
    assert not viewer._edit_mode_active and viewer.knot_artist is None
    assert viewer.canvas.background is None
    assert replots == [{"session_state": viewer.main_window.active_history.current_state}]
    viewer.cancel_continuum_edit()  # nothing to cancel
    assert len(replots) == 1


# --- From the main window ------------------------------------------------------------

def test_edit_and_save_from_the_window(qtbot, window):
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072", "resol": "45000"})
    _run_recipe(qtbot, window, "continuum", "estimate_auto", {})
    viewer = window.plot_viewer
    history = window.active_history
    cont = _spec(viewer).cont.value.copy()
    window._toggle_continuum_editor()  # start
    assert viewer._edit_mode_active and window.edit_cont_button.text() == "Save"
    np.testing.assert_array_equal(viewer._frozen_spline_y, cont)  # on the continuum
    i = len(viewer._knots_x) // 2
    viewer._knots_y[i] *= 1.2
    viewer._apply_local_patch(viewer._knots_x[i])
    edited = viewer._frozen_spline_y.copy()
    n_states = len(history.states)
    window._toggle_continuum_editor()  # save: a new state, without a worker
    assert not viewer._edit_mode_active
    assert len(history.states) == n_states + 1
    # The continuum through every pixel of the draft curve
    np.testing.assert_allclose(history.current_state.spec.cont.value, edited, rtol=1e-12)
