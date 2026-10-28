"""The flux methods of SpectrumV2: scaling, running standard deviation,
smoothing, calibration.

Characterization tests written before the methods were moved to their own
module (phase 5, split of SpectrumV2): the paths that the recipe tests
(``test_recipes_edit_flux.py``) and the snapshots do not reach.
"""
import logging

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.photometry import generate_calibration_curve
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.spectrum_operations import running_std, smooth_spectrum
from astrocook.core.structures import DataColumnV2, SpectrumDataV2

FLUX = au.erg / au.cm**2 / au.s / au.nm
FLAM = au.erg / (au.cm**2 * au.s * au.Angstrom)
SIGMA_KMS = 50.0


@pytest.fixture
def spec():
    rng = np.random.default_rng(1)
    x = np.linspace(400.0, 700.0, 3001)
    y = 1e-16 * (2.0 + 0.001 * (x - 400.0)) + 1e-18 * rng.normal(size=x.size)
    aux = {
        "cont": DataColumnV2(1e-16 * (2.0 + 0.001 * (x - 400.0)), FLUX, "Continuum"),
        "model": DataColumnV2(1e-16 * np.full(x.size, 1.9), FLUX, "Model"),
        "cont_pl": DataColumnV2(1e-16 * np.full(x.size, 2.1), FLUX, "Power law"),
        "telluric_model": DataColumnV2(np.linspace(0.9, 1.0, x.size),
                                       au.dimensionless_unscaled, "Transmission"),
        "flag": DataColumnV2(np.zeros(x.size, dtype=bool), au.dimensionless_unscaled),
    }
    data = SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - 0.05, au.nm),
        xmax=DataColumnV2(x + 0.05, au.nm), y=DataColumnV2(y, FLUX),
        dy=DataColumnV2(np.full(x.size, 1e-18), FLUX), aux_cols=aux, z_em=3.0)
    return SpectrumV2(data, history=["start"])


def _col(spec, name):
    return getattr(spec, name) if name in ("y", "dy") else spec.get_column(name)


# --- scale_flux ---------------------------------------------------------------------

def test_scale_flux_by_a_number(spec):
    new = spec.scale_flux(2.0)
    for name in ("y", "dy"):
        np.testing.assert_array_equal(getattr(new, name).value, 2.0 * getattr(spec, name).value)
    for name in ("cont", "model"):
        np.testing.assert_array_equal(_col(new, name).value, 2.0 * _col(spec, name).value)
        assert new._data.aux_cols[name].description == spec._data.aux_cols[name].description
    # Not scaled: the power law and the transmission
    for name in ("cont_pl", "telluric_model"):
        np.testing.assert_array_equal(_col(new, name).value, _col(spec, name).value)
    assert new.history[-1] == "Scaled flux by 2.000"


def test_scale_flux_by_an_array(spec):
    factor = np.linspace(1.0, 3.0, spec.x.size)
    new = spec.scale_flux(factor)
    np.testing.assert_array_equal(new.y.value, factor * spec.y.value)
    assert new.history[-1] == "Scaled flux by array(mean=2.000)"


# --- calculate_running_std -----------------------------------------------------------

def test_running_std(spec):
    new = spec.calculate_running_std("y", "rms", 11)
    col = _col(new, "rms")
    np.testing.assert_array_equal(col.value, running_std(y=spec.y.value, h=5))
    assert col.unit == FLUX
    assert new._data.aux_cols["rms"].description == "Running StdDev (window=11pix) of y"
    assert new.history[-1] == "Calculated running RMS (window=11pix) on y"


def test_running_std_errors(spec):
    with pytest.raises(ValueError, match="not found"):
        spec.calculate_running_std("missing", "rms", 11)
    with pytest.raises(TypeError, match="non-numerical"):
        spec.calculate_running_std("flag", "rms", 11)


# --- smooth, smooth_column ----------------------------------------------------------

def _smoothed(spec, values):
    return smooth_spectrum(spec._data.x.quantity, values, SIGMA_KMS, spec._data.z_em)


def test_smooth(spec):
    new = spec.smooth(SIGMA_KMS)
    np.testing.assert_array_equal(new.y.value, _smoothed(spec, spec.y.value))
    np.testing.assert_array_equal(new.dy.value, _smoothed(spec, spec.dy.value))
    assert new._data.y.description == "Smoothed Flux"
    for name in ("cont", "model", "telluric_model"):
        np.testing.assert_array_equal(_col(new, name).value,
                                      _smoothed(spec, _col(spec, name).value))
    # Not smoothed: the power law, and the columns that are not numbers
    np.testing.assert_array_equal(_col(new, "cont_pl").value, _col(spec, "cont_pl").value)
    assert new.history[-1] == f"Smoothed spectrum (sigma={SIGMA_KMS} km/s)"


def test_smooth_leaves_a_column_that_is_not_numbers(spec):
    labels = np.array(["x"] * spec.x.size, dtype=object)
    new = spec.update_column("model", labels).smooth(SIGMA_KMS)
    np.testing.assert_array_equal(new._data.aux_cols["model"].values, labels)
    np.testing.assert_array_equal(_col(new, "cont").value,
                                  _smoothed(spec, _col(spec, "cont").value))


@pytest.mark.parametrize("name", ["y", "cont"])
def test_smooth_column(spec, name):
    new = spec.smooth_column(name, SIGMA_KMS)
    np.testing.assert_array_equal(_col(new, name).value,
                                  _smoothed(spec, _col(spec, name).value))
    # The model is kept without renorm_model
    np.testing.assert_array_equal(_col(new, "model").value, _col(spec, "model").value)
    assert new.history[-1] == f"Smoothed column: {name} (sigma={SIGMA_KMS} km/s)"
    assert new.meta == spec.meta and new._data.z_em == spec._data.z_em


def test_smooth_column_renormalizes_the_model(spec):
    new = spec.smooth_column("cont", SIGMA_KMS, renorm_model=True)
    expected = (_col(new, "cont").value * _col(spec, "model").value
                / _col(spec, "cont").value)
    np.testing.assert_allclose(_col(new, "model").value, expected, rtol=1e-14)


def test_smooth_column_errors(spec):
    with pytest.raises(ValueError, match="not found"):
        spec.smooth_column("missing", SIGMA_KMS)
    with pytest.raises(TypeError, match="non-numerical"):
        spec.smooth_column("flag", SIGMA_KMS)


# --- calibrate -----------------------------------------------------------------------

MAGS = "SDSS_g=17.5, SDSS_r=17.0"


def test_calibrate(spec):
    curve = generate_calibration_curve(spec.x, spec.y, MAGS)
    new = spec.calibrate(MAGS)
    for name in ("y", "dy"):
        col = getattr(new, name)
        np.testing.assert_array_equal(col.value, getattr(spec, name).value * curve)
        assert col.unit == FLAM
    for name in ("cont", "model", "cont_pl"):
        np.testing.assert_array_equal(_col(new, name).value, _col(spec, name).value * curve)
        assert _col(new, name).unit == FLAM
    assert new.history[-1] == f"Flux calibrated with: {MAGS}"


def test_calibrate_also_scales_the_transmission(spec):
    # The telluric model is a dimensionless transmission (see
    # SpectrumV2.norm), yet it is multiplied by the calibration curve and
    # given a flux unit (findings register F69)
    curve = generate_calibration_curve(spec.x, spec.y, MAGS)
    new = spec.calibrate(MAGS)
    np.testing.assert_array_equal(_col(new, "telluric_model").value,
                                  _col(spec, "telluric_model").value * curve)
    assert _col(new, "telluric_model").unit == FLAM


def test_calibrate_with_invalid_magnitudes(spec, caplog):
    with caplog.at_level(logging.ERROR), pytest.raises(ValueError, match="Invalid format"):
        spec.calibrate("SDSS_g 17.5")
    assert "Calibration failed" in caplog.text
