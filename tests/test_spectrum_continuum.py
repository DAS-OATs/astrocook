"""The continuum methods of Spectrum on a small synthetic spectrum.

Characterization tests written before the methods were moved to their own
module (phase 5, split of Spectrum): the paths that the snapshots of the
continuum (``test_snapshots_continuum.py``) do not reach, i.e. the
continuum from manual knots, the errors and the renormalization of the
model.
"""
import logging

import astropy.units as au
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, SpectrumData

FLUX = au.erg / au.cm**2 / au.s / au.nm


def _spec(model=False, cont=False):
    x = np.linspace(500.0, 510.0, 201)
    y = 2.0 + 0.01 * (x - 500.0)
    y[90:110] *= 0.4  # an absorption
    aux = {"abs_mask": DataColumn(
        (np.arange(x.size) >= 90) & (np.arange(x.size) < 110),
        au.dimensionless_unscaled)}
    if cont:
        aux["cont"] = DataColumn(np.full(x.size, 2.0), FLUX)
    if model:
        aux["model"] = DataColumn(np.where(aux["abs_mask"].values, 1.0, 2.0), FLUX)
    data = SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x - 0.025, au.nm),
        xmax=DataColumn(x + 0.025, au.nm), y=DataColumn(y, FLUX),
        dy=DataColumn(np.full(x.size, 0.05), FLUX), aux_cols=aux, z_em=3.0)
    return Spectrum(data, history=["start"])


# --- update_continuum_from_knots -------------------------------------------

KNOTS_X = [510.0, 500.0, 505.0]  # not sorted: sorted by the method
KNOTS_Y = [2.2, 2.0, 2.05]


def test_knots_give_a_pchip_continuum():
    spec = _spec()
    new = spec.update_continuum_from_knots(KNOTS_X, KNOTS_Y)
    expected = PchipInterpolator([500.0, 505.0, 510.0], [2.0, 2.05, 2.2])(spec.x.value)
    np.testing.assert_array_equal(new.get_column("cont").value, expected)
    assert new.get_column("cont").unit == FLUX
    assert new.history == ["start", "Updated column 'cont'", "Updated continuum from knots"]
    assert spec.history == ["start"]  # the original is not changed


def test_knots_renormalize_the_model():
    spec = _spec(model=True, cont=True)
    new = spec.update_continuum_from_knots(KNOTS_X, KNOTS_Y)
    cont = new.get_column("cont").value
    # The old normalized model (model / old cont) times the new continuum
    expected = cont * spec.get_column("model").value / 2.0
    np.testing.assert_allclose(new.get_column("model").value, expected, rtol=1e-15)
    kept = spec.update_continuum_from_knots(KNOTS_X, KNOTS_Y, renorm_model=False)
    np.testing.assert_array_equal(kept.get_column("model").value,
                                  spec.get_column("model").value)


def test_knots_with_a_model_and_no_continuum_keep_the_model():
    # Without an old continuum the model cannot be renormalized: it is kept
    # as it is, as in fit_continuum (it raised AttributeError, findings
    # register F67)
    spec = _spec(model=True)
    for renorm in (True, False):
        new = spec.update_continuum_from_knots(KNOTS_X, KNOTS_Y,
                                               renorm_model=renorm)
        assert new.has_aux_column("cont")
        np.testing.assert_array_equal(new.get_column("model").value,
                                      spec.get_column("model").value)


# Invalid knots raise a ValueError with a message (phase 6, F95): they
# returned the spectrum unchanged, and the edit of the continuum was lost
@pytest.mark.parametrize("knots_x, knots_y, message", [
    ([], [], "no knots"), (None, None, "no knots"),
    ([500.0, 510.0], [2.0], "2 wavelengths and 1 fluxes"),
    ([500.0], [2.0], "At least two knots"),
    ([500.0, 500.0], [1.0, 2.0], "Two knots have the same wavelength"),
])
def test_invalid_knots_raise(knots_x, knots_y, message):
    with pytest.raises(ValueError, match=message):
        _spec().update_continuum_from_knots(knots_x, knots_y)


def test_unsorted_knots_build_the_same_continuum():
    spec = _spec()
    a = spec.update_continuum_from_knots([500.0, 505.0, 510.0], [2.0, 2.1, 2.2])
    b = spec.update_continuum_from_knots([510.0, 500.0, 505.0], [2.2, 2.0, 2.1])
    np.testing.assert_array_equal(a.get_column("cont").value, b.get_column("cont").value)


# --- fit_continuum -------------------------------------------------------------

def test_fit_continuum_without_the_mask():
    with pytest.raises(ValueError, match="Run 'find_absorbed' first"):
        _spec().fit_continuum(1.0, 100.0, mask_col="missing")


def test_fit_continuum_without_smoothing_interpolates():
    spec = _spec()
    new = spec.fit_continuum(1.0, 0.0)
    cont = new.get_column("cont").value
    # Through the unabsorbed pixels, interpolated across the absorption
    unabs = ~spec.get_column("abs_mask").value.astype(bool)
    np.testing.assert_allclose(cont[unabs], spec.y.value[unabs], rtol=1e-12)
    np.testing.assert_allclose(cont[100], 2.0 + 0.01 * (spec.x.value[100] - 500.0),
                               rtol=1e-12)
    assert new.history[-1] == ("Fitted continuum (fudge=1.000, smooth=0.0 km/s, "
                               "template=False)")


def test_fit_continuum_with_everything_absorbed(caplog):
    spec = _spec()
    spec = spec.update_column("abs_mask", np.ones(spec.x.size, dtype=bool),
                              unit=au.dimensionless_unscaled)
    # It reaches the recipe (phase 6, F95): the spectrum was returned
    # unchanged, with an ERROR in the log only
    with pytest.raises(ValueError, match="fewer than 2 unabsorbed pixels"):
        spec.fit_continuum(1.0, 0.0)


def test_fit_continuum_with_an_invalid_fudge(caplog):
    spec = _spec()
    with caplog.at_level(logging.WARNING):
        new = spec.fit_continuum("bad", 0.0)
    assert "The fudge 'bad' is not a number: 1.0 is used." in caplog.text
    np.testing.assert_array_equal(new.get_column("cont").value,
                                  spec.fit_continuum(1.0, 0.0).get_column("cont").value)


def test_fit_continuum_scales_with_the_fudge():
    spec = _spec()
    base = spec.fit_continuum(1.0, 0.0).get_column("cont").value
    np.testing.assert_allclose(spec.fit_continuum(1.1, 0.0).get_column("cont").value,
                               1.1 * base, rtol=1e-15)


def test_fit_continuum_renormalizes_the_model():
    spec = _spec(model=True, cont=True)
    new = spec.fit_continuum(1.0, 0.0, renorm_model=True)
    cont = new.get_column("cont").value
    expected = cont * spec.get_column("model").value / 2.0
    np.testing.assert_allclose(new.get_column("model").value, expected, rtol=1e-15)
    kept = spec.fit_continuum(1.0, 0.0)
    np.testing.assert_array_equal(kept.get_column("model").value,
                                  spec.get_column("model").value)


def test_fit_continuum_renormalizes_the_model_with_the_tellurics():
    spec = _spec(model=True, cont=True)
    telluric = np.linspace(0.9, 1.0, spec.x.size)
    spec = spec.update_column("telluric_model", telluric, unit=au.dimensionless_unscaled)
    new = spec.fit_continuum(1.0, 0.0, renorm_model=True)
    cont = new.get_column("cont").value
    # The normalized model is taken with the tellurics, old and new
    expected = (cont * telluric) * spec.get_column("model").value / (2.0 * telluric)
    np.testing.assert_allclose(new.get_column("model").value, expected, rtol=1e-14)


# --- find_absorbed, fit_powerlaw -------------------------------------------------

def test_find_absorbed_adds_the_mask():
    spec = _spec()
    new = spec.find_absorbed(5000.0 * au.km / au.s, 400.0 * au.km / au.s, 2.0, False)
    mask = new.get_column("abs_mask")
    assert mask.unit == au.dimensionless_unscaled and mask.value.dtype == bool
    assert new.history[-1] == "Found absorbed (kappa=2.0)"


def test_fit_powerlaw_adds_its_column():
    spec = _spec()
    new = spec.fit_powerlaw("125.0-126.0")
    assert new.has_aux_column("cont_pl")
    assert new.get_column("cont_pl").unit == FLUX
    assert new.history[-1] == "Fitted power-law (regions=125.0-126.0)"
