"""Expressions, masks and regions of SpectrumV2.

Characterization tests written before the methods were moved to their own
module (phase 5, split of SpectrumV2): the paths that the recipe tests
(``test_recipes_edit_flux.py`` and others) do not reach.
"""
import logging

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, SpectrumDataV2

FLUX = au.erg / au.cm**2 / au.s / au.nm


def _spec(lo=500.0, hi=510.0, n=101, mask=False):
    x = np.linspace(lo, hi, n)
    aux = {"cont": DataColumnV2(np.full(n, 2.0), FLUX)}
    if mask:
        aux["abs_mask"] = DataColumnV2(np.ones(n, dtype=bool), au.dimensionless_unscaled)
    data = SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - 0.05, au.nm),
        xmax=DataColumnV2(x + 0.05, au.nm), y=DataColumnV2(1.0 + 0.1 * (x - lo), FLUX),
        dy=DataColumnV2(np.full(n, 0.05), FLUX), aux_cols=aux, z_em=3.0)
    return SpectrumV2(data, history=["start"])


# --- prepare_expression_context --------------------------------------------------------

def test_context_resamples_the_other_spectra():
    spec, other = _spec(), _spec(495.0, 515.0, n=201)
    expr, extra = spec.prepare_expression_context("y - s1.y + s1.y + s2.y", {"s1": other})
    # A known alias becomes a variable (once); an unknown one is left as it is
    assert expr == "y - s1_y + s1_y + s2.y"
    np.testing.assert_allclose(extra["s1_y"], 1.0 + 0.1 * (spec.x.value - 495.0), rtol=1e-12)
    assert list(extra) == ["s1_y"]


def test_context_with_a_missing_column():
    # get_resampled_column fails on its own message (F71): the error says so
    with pytest.raises(ValueError, match="Failed to resample s1.nope"):
        _spec().prepare_expression_context("s1.nope", {"s1": _spec()})


# --- apply_expression, mask_expression -------------------------------------------------

def test_expression_with_a_scalar_result():
    spec = _spec()
    new = spec.apply_expression("level", "max(y)")
    np.testing.assert_array_equal(new.get_column("level").value,
                                  np.full(spec.x.size, spec.y.value.max()))


@pytest.mark.parametrize("target", ["y", "cont"])
def test_mask_expression(target):
    spec = _spec()
    new = spec.mask_expression(target, "x > 505")
    values = getattr(new, target).value if target == "y" else new.get_column(target).value
    above = spec.x.value > 505.0
    assert np.isnan(values[above]).all() and np.isfinite(values[~above]).all()
    assert new.history[-1] == f"Mask: {target} where x > 505"


def test_mask_expression_with_extra_variables():
    spec = _spec()
    new = spec.mask_expression("y", "lim > 0", extra_vars={"lim": spec.x.value - 505.0})
    assert np.isnan(new.y.value[spec.x.value > 505.0]).all()


def test_mask_expression_must_be_boolean(caplog):
    with caplog.at_level(logging.ERROR), pytest.raises(TypeError, match="boolean"):
        _spec().mask_expression("y", "x * 2")


# --- split ---------------------------------------------------------------------------------

def test_split_with_extra_variables():
    spec = _spec()
    new = spec.split("keep", extra_vars={"keep": spec.x.value < 503.0})
    assert new.x.value.max() < 503.0 and new.x.size == np.sum(spec.x.value < 503.0)


def test_split_errors(caplog):
    spec = _spec()
    with caplog.at_level(logging.ERROR), pytest.raises(TypeError, match="boolean"):
        spec.split("x * 2")
    assert "Numexpr evaluation failed for split" in caplog.text
    with pytest.raises(ValueError, match="empty spectrum"):
        spec.split("x > 1000")


# --- mask_region, get_region_bounds -----------------------------------------------------

def test_mask_region():
    spec = _spec(mask=True)
    new = spec.mask_region("abs_mask", 502.0 * au.nm, 5080.0 * au.Angstrom)
    mask = new.get_column("abs_mask").value.astype(bool)
    x = spec.x.value
    np.testing.assert_array_equal(mask, (x >= 502.0) & (x <= 508.0))


def test_mask_region_without_the_mask(caplog):
    spec = _spec()
    with caplog.at_level(logging.WARNING):
        assert spec.mask_region("abs_mask", 502.0 * au.nm, 508.0 * au.nm) is spec
    assert "Mask column 'abs_mask' not found" in caplog.text


def test_region_bounds():
    lo, hi = _spec().get_region_bounds("lya_forest")
    assert lo.to_value(au.nm) == pytest.approx(102.57 * 4.0)
    assert hi.to_value(au.nm) == pytest.approx(121.567 * 4.0)
    with pytest.raises(ValueError, match="Unknown region"):
        _spec().get_region_bounds("nowhere")
