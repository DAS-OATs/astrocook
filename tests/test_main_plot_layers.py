"""The layers of the main plot (``SpectrumPlotWidget.plot_spectrum``): what
each view setting draws, the arrays normalized or divided by the error, the
zoom kept from one plot to the next, the units of the x axis, the
rest-frame axis and the emission lines, the caches of the arrays and the
decimation of long spectra.

Characterization tests written before ``plot_spectrum`` was split into its
parts (phase 5, split of the GUI); the split was also checked pixel by
pixel on a sequence of view settings (``ROADMAP.md``, phase 5).
"""
import logging

import numpy as np
import pytest
from PySide6.QtWidgets import QMessageBox

from astrocook.gui import pyside_plot
from astrocook.gui.main_window import MainWindowV2

from test_gui_smoke import SESSION_FILE
from test_recipes_edit_flux import _spec

Z_EM = 4.072


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
def columns(window):
    """The archive with z_em, a continuum, a model, a running error, the
    ids of detected regions and a boolean mask."""
    history = window.active_history
    state = history.current_state
    spec = state.spec.with_properties(z_em=Z_EM, resol=45000)
    y = spec.y.value
    n = len(y)
    cont = np.convolve(y, np.ones(51) / 51, mode="same")
    spec = spec.update_column("cont", cont)
    spec = spec.update_column("model", cont * (0.9 + 0.1 * np.cos(np.arange(n) / 7.0)))
    spec = spec.update_column("running_std", np.abs(spec.dy.value) * 1.3)
    spec = spec.update_column("abs_ids", ((np.arange(n) // 40) % 3).astype(float))
    spec = spec.update_column("my_mask", (np.arange(n) % 17) < 5)
    history.add_state(state.with_new_spectrum(spec))
    _replot(window, force_autoscale=True)
    return spec


def _replot(window, **kwargs):
    viewer = window.plot_viewer
    viewer.plot_spectrum(window.active_history.current_state, **kwargs)
    return viewer.canvas.axes


def _set(window, **boxes):
    """Set check boxes without their signals (one plot afterwards)."""
    for name, value in boxes.items():
        box = getattr(window, name)
        box.blockSignals(True)
        box.setChecked(value)
        box.blockSignals(False)


def _lines(ax):
    return {line.get_label(): line for line in ax.get_lines()}


def _labels(ax):
    return [a.get_label() for a in ax.get_children()
            if isinstance(a.get_label(), str) and not a.get_label().startswith("_")]


# --- Layers ------------------------------------------------------------------------------

def test_layers_follow_the_check_boxes(window, columns):
    _set(window, error_checkbox=True, error_line_checkbox=False, continuum_checkbox=True,
         model_checkbox=True, systems_checkbox=True)
    ax = _replot(window)
    lines = _lines(ax)
    spec = columns
    np.testing.assert_array_equal(lines["Spectrum"].get_xdata(), spec.x.value)
    assert lines["Spectrum"].get_drawstyle() == "steps-mid"
    # The error as a band around the flux
    np.testing.assert_allclose(lines["1-sigma error"].get_ydata(), spec.y.value - spec.dy.value)
    np.testing.assert_array_equal(lines["Normalization"].get_ydata(), spec.cont.value)
    np.testing.assert_array_equal(lines["Absorption model"].get_ydata(), spec.model.value)
    assert {"HI Systems", "Metal Systems"} <= set(lines)

    _set(window, error_line_checkbox=True, continuum_checkbox=False, model_checkbox=False,
         systems_checkbox=False)
    lines = _lines(_replot(window))
    np.testing.assert_array_equal(lines["1-sigma error"].get_ydata(), spec.dy.value)
    assert not {"Normalization", "Absorption model", "HI Systems", "Metal Systems"} & set(lines)
    _set(window, error_checkbox=False)
    assert "1-sigma error" not in _lines(_replot(window))


def test_auxiliary_columns(window, columns):
    for col, label in (("running_std", "running_std"), ("abs_ids", "Detected Regions"),
                       ("my_mask", "Mask: my_mask")):
        window.aux_col_combo.blockSignals(True)
        if window.aux_col_combo.findText(col) < 0:
            window.aux_col_combo.addItem(col)
        window.aux_col_combo.setCurrentIndex(window.aux_col_combo.findText(col))
        window.aux_col_combo.blockSignals(False)
        ax = _replot(window)
        assert label in _labels(ax), col
        if col == "running_std":  # a line
            np.testing.assert_array_equal(_lines(ax)[col].get_ydata(),
                                          columns.get_column(col).value)


def test_cursor_lines_are_animated(window, columns):
    window.cursor_series_input.setText("CIV")
    window.cursor_z_input.setText("3.0")
    _set(window, cursor_show_checkbox=True)
    _replot(window)
    viewer = window.plot_viewer
    lines = viewer.canvas.cursor_artists
    assert [line.get_xdata()[0] for line in lines] == pytest.approx(
        viewer.get_cursor_line_positions_at_z(window.active_history.current_state, 3.0))
    assert all(line.get_animated() for line in lines)
    assert lines[0].get_label() == "Cursor (z=3.0000)" and lines[1].get_label().startswith("_")


def test_scatter_style(window, columns):
    columns._data.meta["plot_style"] = "scatter"
    lines = _lines(_replot(window, force_autoscale=True))
    assert lines["Spectrum"].get_linestyle() == "None" and lines["Spectrum"].get_marker() == "."
    assert lines["Spectrum"].get_markersize() == 6.0  # 2156 pixels
    columns._data.meta["plot_style"] = "step"


def test_other_selected_sessions(qtbot, window, columns):
    from PySide6.QtCore import QItemSelectionModel
    window.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(100)
    view = window.session_list_view
    for row in range(view.model().rowCount()):
        view.selectionModel().select(view.model().index(row, 0),
                                     QItemSelectionModel.SelectionFlag.Select)
    lines = _lines(_replot(window))
    ghost = [label for label in lines if label.startswith("Ghost: ")]
    assert len(ghost) == 1
    other = window.session_histories[0].current_state.spec
    np.testing.assert_array_equal(lines[ghost[0]].get_ydata(), other.y.value)
    # In the SNR view they are divided by their error too (before, the
    # default "dF" was looked up as a column and they stayed in flux, F75)
    _set(window, snr_checkbox=True)
    lines = _lines(_replot(window))
    np.testing.assert_allclose(lines[ghost[0]].get_ydata(), other.y.value / other.dy.value)


# --- Normalized and SNR ----------------------------------------------------------------

def test_normalized(window, columns):
    spec = columns
    _set(window, norm_y_checkbox=True, continuum_checkbox=True)
    ax = _replot(window)
    lines = _lines(ax)
    np.testing.assert_allclose(lines["Spectrum"].get_ydata(), spec.y.value / spec.cont.value)
    np.testing.assert_array_equal(lines["Normalization"].get_ydata(), 1.0)
    assert ax.get_ylim() == (-0.3, 1.3)  # switched on: the default range
    assert ax.get_ylabel() == "Normalized Flux (ƒ)"


def test_signal_to_noise(window, columns, caplog):
    spec = columns
    _set(window, snr_checkbox=True)
    with caplog.at_level(logging.ERROR):
        ax = _replot(window)
    np.testing.assert_allclose(_lines(ax)["Spectrum"].get_ydata(), spec.y.value / spec.dy.value)
    assert ax.get_ylabel() == "Signal-to-Noise Ratio (F / dy)"
    # The default "dF" of the box is dy: no error in the log (before, it was
    # looked up as a column, findings register F75), and the arrays are
    # cached for the next plot
    assert "Could not get SNR error column" not in caplog.text
    assert window.plot_viewer._cached_snr_error_col == "dy"
    assert "snr" in window.plot_viewer._plot_cache
    window.snr_col_combo.blockSignals(True)
    if window.snr_col_combo.findText("running_std") < 0:
        window.snr_col_combo.addItem("running_std")
    window.snr_col_combo.setCurrentIndex(window.snr_col_combo.findText("running_std"))
    window.snr_col_combo.blockSignals(False)
    ax = _replot(window)
    np.testing.assert_allclose(_lines(ax)["Spectrum"].get_ydata(),
                               spec.y.value / spec.get_column("running_std").value)
    assert ax.get_ylabel() == "Signal-to-Noise Ratio (F / running_std)"


# --- Zoom ---------------------------------------------------------------------------------

def test_zoom_is_kept(window, columns):
    spec = columns
    x = spec.x.value
    ax = window.plot_viewer.canvas.axes
    ax.set_xlim(x[900], x[1100])
    ax.set_ylim(0.0, 1e-16)
    ax = _replot(window)
    assert ax.get_xlim() == (x[900], x[1100]) and ax.get_ylim() == (0.0, 1e-16)
    # The arrays cut to the view, with 50 pixels on each side
    np.testing.assert_array_equal(_lines(ax)["Spectrum"].get_xdata(), x[850:1151])
    # Home: the whole spectrum
    ax = _replot(window, force_autoscale=True)
    assert len(_lines(ax)["Spectrum"].get_xdata()) == len(x)


# --- Axes ----------------------------------------------------------------------------------

def test_units_of_the_x_axis(window, columns):
    for unit, symbol, factor in (("nm", "nm", 1.0), ("Angstrom", "Å", 10.0),
                                 ("micron", "µm", 1e-3)):
        window.x_unit_combo.blockSignals(True)
        window.x_unit_combo.setCurrentIndex(window.x_unit_combo.findText(unit))
        window.x_unit_combo.blockSignals(False)
        ax = _replot(window)
        assert ax.get_xlabel() == f"Wavelength (λ, {symbol})"
        text = ax.format_coord(500.0, 2e-17)
        rest = 500.0 / (1 + Z_EM)
        assert text == (f"λ = {500.0 * factor:.4f} {symbol}, F = {2e-17:.3e} "
                        f"{columns.y.unit:s}".replace("Angstrom", "Å").replace("**", "^")
                        + f", λᵣ = {rest * factor:.4f} {symbol}")
        formatter = ax.xaxis.get_major_formatter()
        if unit != "nm":
            assert formatter(500.0, 0) == f"{500.0 * factor:g}"


def test_rest_frame_axis_and_emission_lines(window, columns):
    _set(window, strong_lines_checkbox=True)
    ax = _replot(window, force_autoscale=True)
    assert len(ax.child_axes) == 1
    assert ax.child_axes[0].get_xlabel() == "Rest-frame (λᵣ, nm)"
    strong = {t.get_text() for t in ax.texts}
    assert strong and strong <= set(pyside_plot.STRONG_EMISSION_LINES)
    _set(window, strong_lines_checkbox=False)
    ax = _replot(window, force_autoscale=True)
    assert {t.get_text() for t in ax.texts} > strong


def test_no_rest_frame_without_z_em(window):
    ax = _replot(window, force_autoscale=True)  # the archive has no z_em
    assert len(ax.child_axes) == 0 and len(ax.texts) == 0


def test_no_session(window):
    viewer = window.plot_viewer
    viewer.plot_spectrum(None)
    assert viewer.canvas.axes.get_title() == "No Spectrum Data Loaded"
    assert not viewer.canvas.axes.has_data()


def test_legend_without_the_swatches(window, columns):
    # With the Plot Elements section open, the entries shown there as
    # swatches are left out of the legend
    window.plot_elements_section.setExpanded(True)
    legend = _replot(window).get_legend()
    texts = {t.get_text() for t in legend.get_texts()} if legend else set()
    panel_open = window.right_sidebar_widget.width() > 1
    assert not (panel_open and texts & set(pyside_plot.PANEL_LEGEND_LABELS))


# --- Decimation -----------------------------------------------------------------------------

def test_decimation_of_long_spectra(window, monkeypatch):
    """A spectrum longer than the threshold is decimated (min-max) in the
    full view, only at the first plot after the axes were cleared: then the
    limits are kept (was_zoomed) and the decimated arrays stay in the
    cache, drawn as steps (findings register F76)."""
    monkeypatch.setattr(pyside_plot, "DECIMATION_THRESHOLD", 1000)
    long = _spec(x0=600.0, x1=620.0, step=0.01)  # 2000 pixels
    history = window.active_history
    history.add_state(history.current_state.with_new_spectrum(long))
    viewer = window.plot_viewer
    viewer.plot_spectrum(None)
    lines = _lines(_replot(window))
    decimated = lines["Spectrum (decimated)"]
    assert len(decimated.get_xdata()) == 2 * (2000 // 10)
    np.testing.assert_array_equal(decimated.get_ydata(),
                                  pyside_plot.decimate_y_min_max(long.y.value, 10))
    lines = _lines(_replot(window))
    assert "Spectrum (decimated)" not in lines
    assert len(lines["Spectrum"].get_xdata()) == 400  # still the decimated arrays
    # Home: the cache is emptied, and the full view is not decimated
    assert len(_lines(_replot(window, force_autoscale=True))["Spectrum"].get_xdata()) == 2000
