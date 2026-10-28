"""The controls of the continuum editor in the plot controls: Start/Save,
Cancel and Reset, the spacing of the knots (slider on a log scale, and
typed), the cursor switched off while editing, the continuum missing, the
model renormalized or not when saving.

Characterization tests written before the controls were moved to a mixin
(phase 5, split of the GUI). The editor on the plot is in
``test_continuum_editor.py``.
"""
import logging

import numpy as np
import pytest
from PySide6.QtWidgets import QMessageBox

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
    monkeypatch.setattr(win, "_show_custom_message",
                        lambda *a, **k: QMessageBox.StandardButton.Ok)
    win.show()
    win.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(50)
    return win


@pytest.fixture
def with_continuum(qtbot, window):
    _run_recipe(qtbot, window, "edit", "set_properties", {"z_em": "4.072"})
    _run_recipe(qtbot, window, "continuum", "estimate_auto", {})
    return window


@pytest.fixture
def editing(with_continuum):
    with_continuum.edit_cont_button.click()
    assert with_continuum.plot_viewer._edit_mode_active
    return with_continuum


def _answer(window, monkeypatch, button):
    asked = []

    def show(*args, **kwargs):
        asked.append(kwargs.get("header"))
        return getattr(QMessageBox.StandardButton, button)
    monkeypatch.setattr(window, "_show_custom_message", show)
    return asked


# --- Spacing of the knots ------------------------------------------------------------------

def test_spacing_on_a_log_scale(window):
    assert window._get_log_stride(0) == pytest.approx(100.0)
    assert window._get_log_stride(100) == pytest.approx(20000.0)
    assert window._get_log_stride(50) == pytest.approx(100.0 * 200 ** 0.5)
    for stride, position in ((100.0, 0), (20000.0, 100), (1500.0, 51), (50.0, 0), (1e6, 100)):
        window._set_slider_from_stride(stride)
        assert window.stride_slider.value() == position, stride
    # The default: 1500 km/s, as a position of the slider
    assert window._get_log_stride(51) == pytest.approx(1491.16, abs=0.01)


# --- Start -------------------------------------------------------------------------------------

@pytest.mark.parametrize("answer", ["Yes", "Cancel"])
def test_start_without_a_continuum(window, monkeypatch, answer):
    asked = _answer(window, monkeypatch, answer)
    launched = []
    monkeypatch.setattr(window, "_launch_recipe_dialog",
                        lambda cat, name, initial_params=None: launched.append((cat, name)))
    window.edit_cont_button.click()
    assert asked == ["No continuum found to edit."]
    assert launched == ([("continuum", "estimate_auto")] if answer == "Yes" else [])
    assert not window.plot_viewer._edit_mode_active
    assert window.edit_cont_button.text() == "Start"
    # The hints are shown before the check, and stay shown
    assert not window.edit_hints_label.isHidden()


def test_start(with_continuum):
    window = with_continuum
    window.cursor_series_input.setText("CIV")
    window.cursor_z_input.setText("3.0")
    window.cursor_show_checkbox.setChecked(True)
    assert window.plot_viewer.canvas.cursor_artists
    window.edit_cont_button.click()
    viewer = window.plot_viewer
    assert viewer._edit_mode_active
    # The cursor is switched off and its lines removed (not in the background)
    assert not window.cursor_show_checkbox.isChecked() and viewer.canvas.cursor_artists == []
    assert window.edit_cont_button.text() == "Save"
    for widget in (window.cancel_cont_button, window.stride_slider, window.stride_entry,
                   window.reset_cont_button):
        assert widget.isEnabled()
    stride = window._get_log_stride(window.stride_slider.value())
    assert window.stride_entry.text() == f"{stride:.0f}"
    assert window._default_continuum_stride == stride
    x = window.active_history.current_state.spec.x.value
    np.testing.assert_array_equal(viewer._knots_x, x[viewer._generate_velocity_knots(x, stride)])


# --- Spacing while editing -----------------------------------------------------------------

def test_slider_changes_the_spacing(editing, monkeypatch):
    window = editing
    strides = []
    monkeypatch.setattr(window.plot_viewer, "update_continuum_stride", strides.append)
    window.stride_slider.setValue(80)
    assert strides == [pytest.approx(window._get_log_stride(80))]
    assert window.stride_entry.text() == f"{window._get_log_stride(80):.0f}"


def test_typed_spacing(editing, monkeypatch):
    window = editing
    strides = []
    monkeypatch.setattr(window.plot_viewer, "update_continuum_stride", strides.append)
    for typed, used in (("1250", 1250.0), ("50", 100.0), ("99999", 20000.0)):
        window.stride_entry.setText(typed)
        window._on_stride_text_changed()
        assert strides[-1] == used and window.stride_entry.text() == f"{used:.0f}"
    assert len(strides) == 3  # the slider moved without changing the spacing again
    window.stride_entry.setText("1250")
    window._on_stride_text_changed()
    assert window.stride_slider.value() == 47
    window.stride_entry.setText("abc")
    window._on_stride_text_changed()
    assert len(strides) == 4


def test_slider_pressed_and_released(editing):
    window = editing
    viewer = window.plot_viewer
    window._on_stride_slider_pressed()
    assert viewer._stride_backup_x is not None
    window._on_stride_slider_released()
    assert viewer._stride_backup_x is None


def test_reset(editing, monkeypatch):
    window = editing
    resets = []
    monkeypatch.setattr(window.plot_viewer, "reset_continuum_to_original", resets.append)
    strides = []
    monkeypatch.setattr(window.plot_viewer, "update_continuum_stride", strides.append)
    default = window._default_continuum_stride
    window.stride_slider.setValue(90)
    window.reset_cont_button.click()
    assert resets == [default] and len(strides) == 1  # no update from the reset itself
    assert window.stride_entry.text() == f"{default:.0f}"
    assert window.stride_slider.value() == 51


# --- Save, cancel ----------------------------------------------------------------------------

@pytest.mark.parametrize("renorm", [True, False])
def test_save(editing, monkeypatch, renorm):
    window = editing
    history = window.active_history
    viewer = window.plot_viewer
    i = len(viewer._knots_x) // 2
    viewer._knots_y[i] *= 1.2
    viewer._apply_local_patch(viewer._knots_x[i])
    edited = viewer._frozen_spline_y.copy()
    asked = []

    def ask(name, params):
        asked.append((name, sorted(params)))
        return {**params, "renorm_model": str(renorm)}
    monkeypatch.setattr(window, "_ask_to_renormalize_model", ask)
    n = len(history.states)
    window.edit_cont_button.click()
    assert asked == [("update_from_knots", ["knots_x", "knots_y", "renorm_model"])]
    assert len(history.states) == n + 1
    np.testing.assert_allclose(history.current_state.spec.cont.value, edited, rtol=1e-12)
    assert window.edit_cont_button.text() == "Start" and not window.cancel_cont_button.isEnabled()
    assert window.edit_hints_label.isHidden() and window.stride_entry.text() == ""
    assert not hasattr(window, "_default_continuum_stride")


def test_failed_save_is_logged(editing, monkeypatch, caplog):
    window = editing
    recipes = type(window.active_history.current_state.continuum)

    def fail(self, **params):
        raise ValueError("bad knots")
    monkeypatch.setattr(recipes, "update_from_knots", fail)
    with caplog.at_level(logging.ERROR):
        window.edit_cont_button.click()
    assert "Failed to save continuum: bad knots" in caplog.text


@pytest.mark.parametrize("answer", ["Save", "Discard", "Cancel"])
def test_cancel_button(editing, monkeypatch, answer):
    window = editing
    asked = _answer(window, monkeypatch, answer)
    saved = []
    monkeypatch.setattr(window, "_toggle_continuum_editor", lambda: saved.append(True))
    window.cancel_cont_button.click()
    assert asked == ["Unsaved Continuum Edits"]
    if answer == "Save":
        assert saved == [True]
    elif answer == "Discard":
        assert not window.plot_viewer._edit_mode_active
        assert window.edit_cont_button.text() == "Start"
    else:
        assert window.plot_viewer._edit_mode_active


def test_confirm_without_editing(window, monkeypatch):
    asked = _answer(window, monkeypatch, "Cancel")
    assert window._confirm_discard_continuum_edits() is True and asked == []
