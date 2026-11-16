"""The methods of Spectrum that combine or regrid spectra: stitch, rebin,
coadd, resampling.

Characterization tests written before the methods were moved to their own
module (phase 5, split of Spectrum): the paths that the recipe tests and
the snapshots do not reach.
"""
import logging

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, SpectrumData

FLUX = au.erg / au.cm**2 / au.s / au.nm


def _spec(lo, hi, n=201, aux=None, meta=None):
    x = np.linspace(lo, hi, n)
    half = 0.5 * (x[1] - x[0])
    data = SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x - half, au.nm),
        xmax=DataColumn(x + half, au.nm),
        y=DataColumn(1.0 + 0.01 * (x - lo), FLUX),
        dy=DataColumn(np.full(n, 0.05), FLUX),
        aux_cols=aux(x) if aux else {}, meta=meta or {}, z_em=3.0)
    return Spectrum(data, history=["start"])


# --- stitch -----------------------------------------------------------------------------

def test_stitch_pads_missing_columns():
    a = _spec(500.0, 510.0, aux=lambda x: {"cont": DataColumn(np.ones(x.size), FLUX)})
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


def test_rebin_masks_follow_the_new_pixels(caplog):
    # F72: a boolean mask (e.g. abs_mask) follows the new grid, a new pixel
    # being masked if any old pixel in it was; the numerical columns are
    # still interpolated
    spec = _spec(500.0, 510.0, aux=lambda x: {
        "cont": DataColumn(np.ones(x.size), FLUX),
        "abs_mask": DataColumn(x > 505.0, au.dimensionless_unscaled),
    })
    with caplog.at_level(logging.WARNING):
        new = spec.rebin(None, None, 0.1 * au.nm, None, np.nan)
    assert caplog.text == ""
    cont = new.get_column("cont").value
    assert np.isfinite(cont).sum() > 0.9 * cont.size  # NaN only beyond the edges
    np.testing.assert_allclose(cont[np.isfinite(cont)], 1.0)
    mask = new.get_column("abs_mask").value
    assert mask.dtype == bool
    # The old pixels > 505 start at 505.025: the new pixel [505.0, 505.1] has
    # one in it, the one before ([504.9, 505.0]) has none
    assert not mask[new.xmax.value <= 505.0 + 1e-9].any()
    assert mask[new.xmin.value >= 505.0 - 1e-9].all()


def test_rebin_masks_of_a_coarse_grid():
    # Few masked pixels, a coarse grid: any new pixel that holds one is masked
    spec = _spec(500.0, 510.0, aux=lambda x: {
        "abs_mask": DataColumn(np.abs(x - 505.0) < 0.03, au.dimensionless_unscaled)})
    new = spec.rebin(500.0 * au.nm, 510.0 * au.nm, 1.0 * au.nm, None, np.nan)
    mask = new.get_column("abs_mask").value
    # The masked pixels straddle the edge at 505.0: both new pixels are masked
    np.testing.assert_allclose(new.x.value[mask], [504.5, 505.5])


def test_rebin_masks_in_velocity():
    spec = _spec(500.0, 510.0, aux=lambda x: {
        "abs_mask": DataColumn(x > 505.0, au.dimensionless_unscaled)})
    new = spec.rebin(None, None, 10.0 * au.km / au.s, None, np.nan)
    mask = new.get_column("abs_mask").value
    assert not mask[new.xmax.value < 505.0].any()
    assert mask[new.xmin.value > 505.0].all()


def test_follow_mask_edges_and_gaps():
    from astrocook.core.spectrum_operations import follow_mask
    lo = np.array([0.0, 1.0, 2.0, 5.0]); hi = lo + 1.0
    mask = np.array([False, True, False, True])
    new_lo = np.array([-2.0, 0.5, 1.5, 3.0, 6.0, 6.0])
    new_hi = new_lo + 1.0
    # beyond the data (twice); touching only at an edge ([0.5, 1.5] overlaps
    # [1, 2]; [1.5, 2.5] overlaps it too); in a gap; on the masked last one
    np.testing.assert_array_equal(
        follow_mask(lo, hi, mask, new_lo, new_hi),
        [False, True, True, False, False, False])
    # Touching at an edge only is not an overlap
    np.testing.assert_array_equal(
        follow_mask(lo, hi, mask, np.array([2.0, 0.0]), np.array([3.0, 1.0])),
        [False, False])
    # A mask with nothing masked, and an empty one
    assert not follow_mask(lo, hi, np.zeros(4, bool), new_lo, new_hi).any()


def test_rebin_interpolates_the_integer_columns():
    # F72, characterization: an integer column (e.g. the ids of the
    # absorbers) is interpolated as if it were a number, so the result is
    # float, with values between the ids (not an id any more)
    spec = _spec(500.0, 510.0, aux=lambda x: {
        "ids": DataColumn(np.where(x > 505.0, 2, 1), au.dimensionless_unscaled)})
    new = spec.rebin(None, None, 0.07 * au.nm, None, np.nan)
    ids = new.get_column("ids").value
    assert ids.dtype.kind == "f"
    inner = ids[np.isfinite(ids)]
    assert np.any((inner > 1.0) & (inner < 2.0))


def test_coadd_warns_that_it_drops_all_columns_but_the_resolution(caplog):
    # F72: no auxiliary column of the inputs (a mask included) is carried to
    # the co-added spectrum (the masks are made again on the result): said
    mask = lambda x: {"abs_mask": DataColumn(x > 505.0, au.dimensionless_unscaled),
                      "cont": DataColumn(np.ones(x.size), FLUX)}
    a, b = _spec(500.0, 510.0, aux=mask), _spec(505.0, 515.0, aux=mask)
    with caplog.at_level(logging.WARNING):
        new = a.coadd([b], None, None, 0.05 * au.nm, equalize_order=-1)
    assert not new.has_aux_column("abs_mask") and not new.has_aux_column("cont")
    assert "['abs_mask', 'cont'] of the spectra are not carried" in caplog.text


# --- resample_on_grid, get_resampled_column ---------------------------------------------

def test_resample_masks_follow_the_target_pixels(caplog):
    spec = _spec(500.0, 510.0, aux=lambda x: {
        "cont": DataColumn(np.ones(x.size), FLUX),
        "abs_mask": DataColumn(x > 505.0, au.dimensionless_unscaled),
    })
    grid = _spec(502.0, 508.0, n=61)
    with caplog.at_level(logging.WARNING):
        new = spec.resample_on_grid(grid)
    assert caplog.text == ""
    np.testing.assert_array_equal(new.x.value, grid.x.value)
    np.testing.assert_allclose(new.y.value, 1.0 + 0.01 * (grid.x.value - 500.0), rtol=1e-12)
    assert new.has_aux_column("cont")
    mask = new.get_column("abs_mask").value
    assert mask.dtype == bool
    assert not mask[grid.xmax.value < 505.0].any()
    assert mask[grid.xmin.value > 505.0].all()
    # Outside the data nothing is masked
    wide = spec.resample_on_grid(_spec(495.0, 515.0, n=201))
    wmask = wide.get_column("abs_mask").value
    assert not wmask[wide.x.value < 500.0].any() and wmask[wide.x.value > 506.0].any()


def test_resampled_column():
    spec = _spec(500.0, 510.0)
    grid = np.array([499.0, 502.5, 505.0])
    out = spec.get_resampled_column("y", grid)
    assert np.isnan(out[0])  # outside: the fill value
    np.testing.assert_allclose(out[1:], 1.0 + 0.01 * (grid[1:] - 500.0), rtol=1e-12)


def test_resampled_column_errors():
    spec = _spec(500.0, 510.0, aux=lambda x: {
        "abs_mask": DataColumn(x > 505.0, au.dimensionless_unscaled)})
    # The message once referred to self.name, which Spectrum does not
    # have: AttributeError instead of ValueError (findings register F71)
    with pytest.raises(ValueError, match="'missing' not found"):
        spec.get_resampled_column("missing", np.array([505.0]))
    with pytest.raises(TypeError, match="non-numerical"):
        spec.get_resampled_column("abs_mask", np.array([505.0]))
