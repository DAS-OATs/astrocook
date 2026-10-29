"""The methods of SpectrumV2 that combine or regrid spectra: stitch, rebin,
coadd, resampling.

Characterization tests written before the methods were moved to their own
module (phase 5, split of SpectrumV2): the paths that the recipe tests and
the snapshots do not reach.
"""
import logging

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, SpectrumDataV2

FLUX = au.erg / au.cm**2 / au.s / au.nm


def _spec(lo, hi, n=201, aux=None, meta=None):
    x = np.linspace(lo, hi, n)
    half = 0.5 * (x[1] - x[0])
    data = SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - half, au.nm),
        xmax=DataColumnV2(x + half, au.nm),
        y=DataColumnV2(1.0 + 0.01 * (x - lo), FLUX),
        dy=DataColumnV2(np.full(n, 0.05), FLUX),
        aux_cols=aux(x) if aux else {}, meta=meta or {}, z_em=3.0)
    return SpectrumV2(data, history=["start"])


# --- stitch -----------------------------------------------------------------------------

def test_stitch_pads_missing_columns():
    a = _spec(500.0, 510.0, aux=lambda x: {"cont": DataColumnV2(np.ones(x.size), FLUX)})
    b = _spec(490.0, 499.0)
    new = a.stitch([b])
    assert np.all(np.diff(new.x.value) > 0)  # sorted by default
    cont = new.get_column("cont").value
    assert np.isnan(cont[new.x.value < 500.0]).all()
    assert (cont[new.x.value >= 500.0] == 1.0).all()


def test_stitch_without_sorting():
    a, b = _spec(500.0, 510.0), _spec(490.0, 499.0)
    new = a.stitch([b], sort=False)
    np.testing.assert_array_equal(new.x.value, np.concatenate([a.x.value, b.x.value]))


# --- rebin --------------------------------------------------------------------------------

def test_rebin_in_velocity():
    spec = _spec(500.0, 510.0)
    new = spec.rebin(None, None, 10.0 * au.km / au.s, None, np.nan)
    assert new.x.unit == au.nm
    dv = np.diff(np.log(new.x.value)) * 299792.458
    np.testing.assert_allclose(dv, 10.0, rtol=1e-6)


def test_rebin_drops_the_masks(caplog):
    # A boolean mask (e.g. abs_mask) is not interpolated: dropped, with a
    # warning; the numerical columns are kept
    spec = _spec(500.0, 510.0, aux=lambda x: {
        "cont": DataColumnV2(np.ones(x.size), FLUX),
        "abs_mask": DataColumnV2(x > 505.0, au.dimensionless_unscaled),
    })
    with caplog.at_level(logging.WARNING):
        new = spec.rebin(None, None, 0.1 * au.nm, None, np.nan)
    cont = new.get_column("cont").value
    assert np.isfinite(cont).sum() > 0.9 * cont.size  # NaN only beyond the edges
    np.testing.assert_allclose(cont[np.isfinite(cont)], 1.0)
    assert not new.has_aux_column("abs_mask")
    assert "Failed to interpolate aux col 'abs_mask'" in caplog.text


# --- coadd --------------------------------------------------------------------------------

def test_coadd_without_resolution_everywhere(caplog):
    a = _spec(500.0, 510.0, aux=lambda x: {
        "resol": DataColumnV2(np.full(x.size, 45000.0), au.dimensionless_unscaled)},
        meta={"plot_style": "step"})
    b = _spec(505.0, 515.0)
    with caplog.at_level(logging.WARNING):
        new = a.coadd([b], None, None, 0.05 * au.nm, equalize_order=-1)
    assert "Resolution mismatch" in caplog.text
    assert "plot_style" not in new.meta and new.meta["ORIGIN"] == "coadd"
    assert new.x.value.min() < 500.5 and new.x.value.max() > 514.5


# --- resample_on_grid, get_resampled_column ---------------------------------------------

def test_resample_drops_the_masks(caplog):
    spec = _spec(500.0, 510.0, aux=lambda x: {
        "cont": DataColumnV2(np.ones(x.size), FLUX),
        "abs_mask": DataColumnV2(x > 505.0, au.dimensionless_unscaled),
    })
    grid = _spec(502.0, 508.0, n=61)
    with caplog.at_level(logging.WARNING):
        new = spec.resample_on_grid(grid)
    np.testing.assert_array_equal(new.x.value, grid.x.value)
    np.testing.assert_allclose(new.y.value, 1.0 + 0.01 * (grid.x.value - 500.0), rtol=1e-12)
    assert new.has_aux_column("cont") and not new.has_aux_column("abs_mask")
    assert "Cannot resample non-numerical column 'abs_mask'" in caplog.text


def test_resampled_column():
    spec = _spec(500.0, 510.0)
    grid = np.array([499.0, 502.5, 505.0])
    out = spec.get_resampled_column("y", grid)
    assert np.isnan(out[0])  # outside: the fill value
    np.testing.assert_allclose(out[1:], 1.0 + 0.01 * (grid[1:] - 500.0), rtol=1e-12)


def test_resampled_column_errors():
    spec = _spec(500.0, 510.0, aux=lambda x: {
        "abs_mask": DataColumnV2(x > 505.0, au.dimensionless_unscaled)})
    # The message once referred to self.name, which SpectrumV2 does not
    # have: AttributeError instead of ValueError (findings register F71)
    with pytest.raises(ValueError, match="'missing' not found"):
        spec.get_resampled_column("missing", np.array([505.0]))
    with pytest.raises(TypeError, match="non-numerical"):
        spec.get_resampled_column("abs_mask", np.array([505.0]))
