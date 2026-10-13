"""V1/V2 parity tests (ROADMAP.md, phase 1, item 11).

The V1 code in ``astrocook/legacy/`` is the validated reference behind the
published papers. Each test runs V1 and V2 on the same input: where the
algorithm is the same the outputs must agree; known differences are strict
xfails pointing to the findings register (ROADMAP.md, section 3.2) or
characterizations of a documented difference.
"""
from types import SimpleNamespace

import numpy as np
import astropy.units as au
import pytest

from astrocook.core.spectrum import SpectrumV2
from astrocook.core.spectrum_operations import rebin_spectrum, smooth_spectrum
from astrocook.core.structures import DataColumnV2, SpectrumDataV2
from astrocook.legacy.spectrum import Spectrum as SpectrumV1

FLUX_UNIT = au.erg / (au.cm**2 * au.s * au.nm)


def _arrays(start, stop, step, scale=1.0, noise=0.02, seed=0, spikes=()):
    """Sloped continuum with absorption lines, heteroscedastic noise."""
    rng = np.random.default_rng(seed)
    x = np.arange(start, stop, step)
    flux = 1.0 + 0.03 * (x - 500.0)
    for centre, depth, width in [(501.0, 0.7, 0.02), (503.2, 0.4, 0.05),
                                 (506.5, 0.9, 0.015)]:
        flux = flux - depth * np.exp(-0.5 * ((x - centre) / width) ** 2)
    dy = noise * (0.5 + flux)
    y = scale * (flux + dy * rng.standard_normal(len(x)))
    for i in spikes:
        y[i] += 1.0
    return x, step, y, scale * dy


def _v1(x, step, y, dy):
    return SpectrumV1(x * au.nm, (x - step / 2) * au.nm, (x + step / 2) * au.nm,
                      y * FLUX_UNIT, dy * FLUX_UNIT, xunit=au.nm, yunit=FLUX_UNIT)


def _v2(x, step, y, dy):
    return SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - step / 2, au.nm),
        xmax=DataColumnV2(x + step / 2, au.nm),
        y=DataColumnV2(y, FLUX_UNIT), dy=DataColumnV2(dy, FLUX_UNIT)))


def _rebin_both(arrays, xstart, xend, dx, kappa):
    """Rebin with V1 and V2 on the same bins.

    V1 takes the centre of the first bin as start, V2 its left edge: V1 is
    given ``xstart + dx/2`` so that the bins coincide. The first and last
    V1 bins are half-width (edges from midpoints) and are dropped.
    """
    s1 = _v1(*arrays)
    v1 = s1._rebin(xstart + dx / 2, xend, dx, au.nm, s1.y, s1.dy, kappa=kappa)
    s2 = _v2(*arrays)._data
    x2, _, _, y2, dy2, _ = rebin_spectrum(
        s2.x.quantity, s2.xmin.quantity, s2.xmax.quantity, s2,
        xstart=xstart * au.nm, xend=xend * au.nm, dx=dx * au.nm,
        calc_unit=au.nm, kappa=kappa, filling=np.nan)
    n = min(len(v1.x), len(x2))
    inner = slice(1, n - 1)
    np.testing.assert_allclose(v1.x.value[inner], x2.value[inner], rtol=0, atol=1e-9)
    return (v1.y.value[inner], v1.dy.value[inner]), (y2.value[inner], dy2.value[inner])


@pytest.mark.parametrize("dx", [0.01, 0.02])
def test_rebin_flux(dx):
    (y1, _), (y2, _) = _rebin_both(_arrays(500.0, 508.0, 0.004), 501.0, 507.0, dx, None)
    np.testing.assert_allclose(y2, y1, rtol=1e-10)


@pytest.mark.xfail(strict=True, reason=(
    "V2 rebinned errors do not decrease with the number of pixels per bin; "
    "V1 gives ~sigma/sqrt(N). Findings register F1, issue #490."))
def test_rebin_errors():
    (_, dy1), (_, dy2) = _rebin_both(_arrays(500.0, 508.0, 0.004), 501.0, 507.0, 0.02, None)
    np.testing.assert_allclose(dy2, dy1, rtol=1e-6)


@pytest.mark.xfail(strict=True, reason=(
    "With 3 pixels per bin, V2 clips against a first-pass mean that includes "
    "the outlier and can empty the bin (NaN); V1 keeps the bin. Note that V1 "
    "does not remove the spike either (sigma_clip on 3-5 values cannot). "
    "Findings register F2."))
def test_rebin_clipping_keeps_bins_with_a_spike():
    arrays = _arrays(500.0, 508.0, 0.004, spikes=(400, 1200))
    (y1, _), (y2, _) = _rebin_both(arrays, 501.0, 507.0, 0.012, 3.0)
    (y1_raw, _), _ = _rebin_both(arrays, 501.0, 507.0, 0.012, None)
    # V1: the bins are kept, but the spikes are not removed
    assert np.all(np.isfinite(y1))
    np.testing.assert_allclose(y1, y1_raw)
    # V2 must at least keep the bins
    assert np.all(np.isfinite(y2))


def test_coadd_flux():
    # V1 coadd: combine the exposures (append the bins) and rebin.
    # V2 coadd without equalization and clipping, on the same bins.
    a = _arrays(499.0, 508.0, 0.004, noise=0.02, seed=1)
    b = _arrays(499.5, 508.5, 0.005, noise=0.03, seed=2)
    dx, xstart, xend = 0.01, 500.0, 507.0

    s1 = _v1(*a)
    s1._append(_v1(*b))
    v1 = s1._rebin(xstart + dx / 2, xend, dx, au.nm, s1.y, s1.dy, kappa=None)

    v2 = _v2(*a).coadd([_v2(*b)], xstart=xstart * au.nm, xend=xend * au.nm,
                       dx=dx * au.nm, kappa=None, equalize_order=-1)
    n = min(len(v1.x), len(v2.x))
    inner = slice(1, n - 1)
    np.testing.assert_allclose(v1.x.value[inner], v2.x.value[inner], atol=1e-9)
    np.testing.assert_allclose(v2.y.value[inner], v1.y.value[inner], rtol=1e-10)


def test_smoothing_width_differs_by_sqrt2():
    # Characterization (findings register F16): for the same "std" parameter,
    # the V1 kernel is exp(-(dv/std)^2), i.e. sigma = std/sqrt(2); V2 uses a
    # Gaussian with sigma = std.
    x = np.linspace(500.0, 510.0, 4001)
    step = x[1] - x[0]
    y = np.zeros_like(x)
    y[2000] = 1.0
    s1 = _v1(x, step, y, np.full_like(x, 0.1))
    s1._gauss_convolve(std=20 * au.km / au.s, verb=False)
    v1 = np.array(s1._t["conv"])
    v2 = smooth_spectrum(x * au.nm, y, 20.0)

    def sigma_kms(p):
        v = (x - 505.0) / 505.0 * 299792.458
        p = p / p.sum()
        mean = (p * v).sum()
        return np.sqrt((p * (v - mean) ** 2).sum())

    assert sigma_kms(v1) == pytest.approx(20.0 / np.sqrt(2), rel=1e-2)
    assert sigma_kms(v2) == pytest.approx(20.0, rel=1e-2)


def test_equalization_factor():
    # V1: ratio of the medians of y in a window. V2: equalize_to_reference
    # (order 0). On a featureless overlap both must recover the scale.
    from astrocook.legacy.cookbook_general import CookbookGeneral

    def flat(start, stop, step, scale, seed):
        rng = np.random.default_rng(seed)
        x = np.arange(start, stop, step)
        y = scale * (1 + 0.02 * rng.standard_normal(len(x)))
        return x, step, y, np.full_like(x, 0.02 * scale)

    a, b = flat(500.0, 506.0, 0.004, 1.0, 1), flat(503.0, 509.0, 0.005, 0.8, 2)
    a1, b1 = _v1(*a), _v1(*b)
    cookbook = CookbookGeneral()
    cookbook.sess = SimpleNamespace(_gui=SimpleNamespace(_sess_list=[
        SimpleNamespace(spec=a1), SimpleNamespace(spec=b1)]))
    cookbook.equalize(503.0, 506.0, _sel=[0, 1])
    factor_v1 = np.median(b1.y.value) / np.median(b[2])

    b2 = _v2(*b).equalize_to_reference(_v2(*a), order=0)
    factor_v2 = np.median(b2.y.value) / np.median(b[2])

    assert factor_v1 == pytest.approx(1.25, rel=5e-3)
    assert factor_v2 == pytest.approx(factor_v1, rel=5e-3)
