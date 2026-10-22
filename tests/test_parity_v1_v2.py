"""V1/V2 parity tests (ROADMAP.md, phase 1, item 11).

The V1 code is the validated reference behind the published papers; the
routines used here are verbatim copies of ``astrocook/legacy/`` in
``tests/reference/v1/`` (phase 3), so that they outlive the removal of the
V1 code from the package. Each test runs V1 and V2 on the same input: where the
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
from reference.v1.spectrum import Spectrum as SpectrumV1

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


def test_rebin_errors():
    (_, dy1), (_, dy2) = _rebin_both(_arrays(500.0, 508.0, 0.004), 501.0, 507.0, 0.02, None)
    np.testing.assert_allclose(dy2, dy1, rtol=1e-6)


def test_rebin_clipping_keeps_bins_with_a_spike():
    # F2 (findings register): V2 used to clip against a first-pass mean that
    # included the outlier and could empty the bin (NaN). V1 keeps the bin
    # but cannot remove a spike with 3 pixels per bin; V2 now clips against
    # a running median (upward only) and removes it.
    spikes = (400, 1200)
    arrays = _arrays(500.0, 508.0, 0.004, spikes=spikes)
    clean = _arrays(500.0, 508.0, 0.004)
    (y1, _), (y2, dy2) = _rebin_both(arrays, 501.0, 507.0, 0.012, 3.0)
    (y1_raw, _), (y2_clean, _) = _rebin_both(clean, 501.0, 507.0, 0.012, None)
    (y1_spiked, _), _ = _rebin_both(arrays, 501.0, 507.0, 0.012, None)
    # V1: the bins are kept, but the spikes are not removed
    assert np.all(np.isfinite(y1))
    np.testing.assert_allclose(y1, y1_spiked)
    # V2: the bins are kept and the spikes are gone
    assert np.all(np.isfinite(y2))
    assert np.all(np.abs(y2 - y2_clean) < 5 * dy2)


def test_rebin_clipping_keeps_narrow_absorption():
    # F2: a two-sided clipping against a running median removed the cores of
    # narrow lines; the clipping is upward only, so absorption is untouched
    x = np.arange(500.0, 508.0, 0.004)
    sigma = 2.5 * 0.004 / 2.3548  # FWHM 2.5 pixels
    flux = 1.0 - 0.8 * np.exp(-0.5 * ((x - 504.0) / sigma) ** 2)
    dy = np.full_like(x, 0.02)
    y = flux + dy * np.random.default_rng(1).standard_normal(len(x))
    arrays = (x, 0.004, y, dy)
    (_, _), (y_clip, _) = _rebin_both(arrays, 501.0, 507.0, 0.016, 3.0)
    (_, _), (y_raw, _) = _rebin_both(arrays, 501.0, 507.0, 0.016, None)
    # Bins within 5 sigma of the line (elsewhere a few upward noise
    # fluctuations beyond 3 sigma are clipped, as expected)
    x_out = 501.0 + 0.016 * (np.arange(len(y_raw)) + 1.5)
    line = np.abs(x_out - 504.0) < 5 * sigma
    assert line.sum() >= 2
    np.testing.assert_array_equal(y_clip[line], y_raw[line])


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
    from reference.v1.cookbook_general import CookbookGeneral

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


# --- Voigt model ---------------------------------------------------------

def _voigt_v1(x, step, series, z, logN, b, resol):
    """V1 model as in SystModel: lines_voigt convolved with spec.psf_gauss."""
    from reference.v1.functions import convolve_simple, lines_voigt
    one = np.ones_like(x)
    spec = _v1(x, step, one, 0.02 * one)
    model = lines_voigt(x, z, logN, b, 0.0, series)
    return model if resol is None else convolve_simple(model, spec.psf_gauss(x, resol))


def _voigt_v2(x, step, series, z, logN, b, resol=None, resol_col=None):
    from astrocook.core.structures import ComponentDataV2, SystemListDataV2
    from astrocook.core.system_list import SystemListV2
    from astrocook.fitting.voigt_fitter import VoigtFitterV2
    one = np.ones_like(x)
    unit = au.dimensionless_unscaled
    aux = {"cont": DataColumnV2(one, unit)}
    if resol_col is not None:
        aux["resol"] = DataColumnV2(resol_col, unit)
    spec = SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - step / 2, au.nm),
        xmax=DataColumnV2(x + step / 2, au.nm), y=DataColumnV2(one, unit),
        dy=DataColumnV2(0.02 * one, unit), aux_cols=aux, z_em=4.0,
        resol=resol or 0.0))
    comp = ComponentDataV2(id=1, z=z, dz=None, logN=logN, dlogN=None, b=b,
                           db=None, series=series)
    _, model = VoigtFitterV2(
        spec, SystemListV2(SystemListDataV2(components=[comp]))).compute_model_flux()
    return model


VOIGT_X = np.arange(612.0, 620.0, 0.002)  # CIV at z = 2.97, Ly-alpha at z = 4.06


@pytest.fixture
def v1_civ_wavelengths(monkeypatch):
    """The CIV wavelengths of V1 (Morton 2003) in the V2 table, to compare
    the physics: V2 takes those of Petitjean & Aracil (2004) since #568."""
    from astrocook.core.atomic_data import ATOM_DATA
    for name, wave in (("CIV_1548", 1548.204), ("CIV_1550", 1550.781)):
        monkeypatch.setitem(ATOM_DATA[name], "wave", wave)


@pytest.mark.parametrize("resol", [None, 45000.0, 20000.0])
@pytest.mark.parametrize("series, z, logN, b, ew_rtol", [
    ("CIV", 2.97, 13.6, 12.0, 1e-3),
    # V1 uses f = 0.416 for Ly-alpha, V2 0.4164: 0.1% more absorption in V2
    ("Ly_a", 4.059, 13.8, 25.0, 2e-3),
])
def test_voigt_model(series, z, logN, b, ew_rtol, resol, v1_civ_wavelengths):
    # Same physics: the models agree within 5e-4 in normalised flux. A
    # constant 0.05% difference in equivalent width remains for CIV, with
    # identical atomic data; not traced further (2e-4 dex in logN).
    v1 = _voigt_v1(VOIGT_X, 0.002, series, z, logN, b, resol)
    v2 = _voigt_v2(VOIGT_X, 0.002, series, z, logN, b, resol)
    assert np.max(np.abs(v2 - v1)) < 5e-4 + 0.5 * ew_rtol
    assert np.sum(1 - v2) == pytest.approx(np.sum(1 - v1), rel=ew_rtol)


def test_voigt_model_with_resolution_column(v1_civ_wavelengths):
    # Regression test for F10 (V2 left the pixels with a non-round
    # resolution at 0); V1 takes the resolution at the system centre
    resol = 45050.0
    v1 = _voigt_v1(VOIGT_X, 0.002, "CIV", 2.97, 13.6, 12.0, resol)
    v2 = _voigt_v2(VOIGT_X, 0.002, "CIV", 2.97, 13.6, 12.0,
                   resol_col=np.full(len(VOIGT_X), resol))
    assert np.max(np.abs(v2 - v1)) < 1e-3


# --- Continuum (J0103-1305 window) -------------------------------------------

@pytest.fixture(scope="module")
def continua(j0103_window):
    """V1 clip_flux and V2 find_absorbed + fit_continuum on the same window.

    The V2 smoothing is set to the effective V1 one (std 400 km/s in V1 is
    sigma = 283 km/s, see F16), to compare the algorithms rather than the
    smoothing convention.
    """
    from reference.v1.cookbook_continuum import CookbookContinuum
    from conftest import J0103_Z_EM

    w = j0103_window

    def v1_continuum(fudge):
        spec = SpectrumV1(w.x, w.xmin, w.xmax, w.y, w.dy,
                          xunit=au.nm, yunit=w.y.unit)
        cookbook = CookbookContinuum()
        cookbook.sess = SimpleNamespace(spec=spec, lines=None, nodes=None)
        cookbook.nodes_extract = lambda *a, **k: 0  # knots for the GUI editor
        cookbook.clip_flux(zem=J0103_Z_EM, fudge=fudge)
        return (np.array(spec._t["cont"]),
                np.array(spec._t["mask_unabs"]).astype(bool))

    cont1, unabs1 = v1_continuum("auto")
    cont1_nofudge, _ = v1_continuum("1.0")
    absorbed = w.find_absorbed(5000 * au.km / au.s, 400 * au.km / au.s, 2.0, False)
    unabs2 = ~absorbed.get_column("abs_mask").value.astype(bool)
    cont2 = absorbed.fit_continuum("auto", 400.0 / np.sqrt(2)).get_column("cont").value
    forest = w.x.value < 121.567 * (1 + J0103_Z_EM)
    return SimpleNamespace(cont1=cont1, cont1_nofudge=cont1_nofudge,
                           unabs1=unabs1, cont2=cont2, unabs2=unabs2,
                           forest=forest, red=~forest)


def test_continuum_red_side(continua):
    c = continua
    # Masks agree except a few pixels at the red edge of the window
    assert np.mean(c.unabs2[c.red] == c.unabs1[c.red]) > 0.999
    ratio = c.cont2[c.red] / c.cont1[c.red]
    assert np.nanmedian(ratio) == pytest.approx(1.0, abs=0.005)


def test_v1_auto_fudge_is_small(continua):
    # The V1 automatic fudge is effective but small on this spectrum: it
    # corrects the running mean V1 builds the continuum from (5th-95th
    # percentile 0.998-1.002 in the forest, 0.999-1.000 on the red side)
    c = continua
    ratio = c.cont1 / c.cont1_nofudge
    assert np.nanmedian(ratio[c.forest]) == pytest.approx(1.0, abs=1e-3)
    assert np.nanmedian(ratio[c.red]) == pytest.approx(1.0, abs=1e-3)


def test_v2_auto_fudge_is_one(j0103_window):
    # Characterization (findings register F3, accepted): V2 interpolates the
    # unabsorbed pixels, so the residuals the automatic fudge is computed
    # from are zero and the fudge is exactly 1. A V1-like local fudge was
    # tried and does not reproduce V1 (the V2 continuum is not a running
    # mean): see the F3 entry of the register.
    absorbed = j0103_window.find_absorbed(5000 * au.km / au.s, 400 * au.km / au.s,
                                          2.0, False)
    auto = absorbed.fit_continuum("auto", 400.0).get_column("cont").value
    one = absorbed.fit_continuum("1.0", 400.0).get_column("cont").value
    np.testing.assert_array_equal(auto, one)


@pytest.mark.xfail(strict=True, reason=(
    "In the Lyman-alpha forest V2 marks more pixels as unabsorbed than V1 "
    "(57.6% vs 54.8% of the window), and its continuum is ~1.2% lower "
    "(median; up to 5%). Findings register F17."))
def test_continuum_forest(continua):
    c = continua
    assert np.mean(c.unabs2[c.forest] == c.unabs1[c.forest]) > 0.99
    ratio = c.cont2[c.forest] / c.cont1[c.forest]
    assert np.nanmedian(ratio) == pytest.approx(1.0, abs=0.005)
