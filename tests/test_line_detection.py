"""The detection of the lines by matched filters (phase 7, PR 7c.4b; F57)."""
import numpy as np
import pytest

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.line_detection import (
    C_KMS, DEFAULT_B_KMS, LINE_SIGMA, MULTIPLET_SIGMA, detect_lines,
    filter_sigma_kms, filter_sums,
    joint_significance, matched_filter_detection, multiplet_lines,
    multiplet_seeds,
)
from astrocook.mocks import MockComponent, MockTruth, generate

R = 45000.0
DV = 2.5


def _grid(n, dv=DV, lam0=400.0):
    v = np.arange(n) * dv
    return v, lam0 * np.exp(v / C_KMS)


def _gauss(v, centre, sigma):
    return np.exp(-0.5 * ((v - centre) / sigma) ** 2)


def test_on_pure_noise_s_is_a_standard_normal_for_each_width():
    rng = np.random.default_rng(1)
    _, x = _grid(80000)
    err = np.full(len(x), 0.05)
    det = matched_filter_detection(x, 1 + rng.normal(0, 0.05, len(x)), err, R,
                                   multiplets=())
    for s in det.significance:
        assert abs(np.mean(s)) < 0.1
        assert abs(np.std(s) - 1) < 0.1
    # The one-sided tail beyond 3 sigma, 0.00135, on the narrowest filter
    # (the broader ones have too few independent values)
    assert 0.0008 < np.mean(det.significance[0] > 3) < 0.002


def test_on_a_grid_linear_in_wavelength_too():
    rng = np.random.default_rng(2)
    x = np.linspace(400.0, 600.0, 60000)
    det = matched_filter_detection(x, 1 + rng.normal(0, 0.05, len(x)),
                                   np.full(len(x), 0.05), 8900.0, multiplets=())
    for s in det.significance:
        assert abs(np.std(s) - 1) < 0.1


def test_the_significance_of_a_matched_line_is_its_amplitude_over_its_error():
    v, x = _grid(4000)
    b = 15.0
    sigma = filter_sigma_kms(np.array([R]), b)[0]
    err = np.full(len(x), 0.02)
    amp = 0.05
    flux = 1 - amp * _gauss(v, v[2000], sigma)
    det = matched_filter_detection(x, flux, err, R, b_kms=(b,), multiplets=())
    # S = A sqrt(sum w^2) / err for the filter of the line, at its centre
    expected = amp * np.sqrt(np.sum(_gauss(v, v[2000], sigma) ** 2)) / 0.02
    assert det.significance[0, 2000] == pytest.approx(expected, rel=1e-3)
    assert np.argmax(det.significance[0]) == 2000


def test_the_filter_combines_the_instrumental_profile_and_b():
    sigma = filter_sigma_kms(np.array([C_KMS / 10.0]), 2.0)
    assert sigma[0] == pytest.approx(np.sqrt((10 / 2.3548200450309493) ** 2 + 2.0))


def test_a_line_and_its_region_by_the_two_thresholds():
    v, x = _grid(4000)
    sigma = filter_sigma_kms(np.array([R]), 15.0)[0]
    err = np.full(len(x), 0.02)
    flux = 1 - 0.04 * _gauss(v, v[2000], sigma)
    det = matched_filter_detection(x, flux, err, R, multiplets=())
    best = det.max_significance
    assert det.n_regions == 1
    assert det.lines[2000] and det.mask[2000]
    # The line where S > 3, its region where S > 1, around it
    np.testing.assert_array_equal(det.lines, best > 3)
    np.testing.assert_array_equal(det.mask, best > 1)
    assert det.mask.sum() > det.lines.sum()


def test_a_weak_line_below_the_threshold_gives_no_region():
    v, x = _grid(4000)
    sigma = filter_sigma_kms(np.array([R]), 15.0)[0]
    flux = 1 - 0.01 * _gauss(v, v[2000], sigma)
    det = matched_filter_detection(x, flux, np.full(len(x), 0.02), R, multiplets=())
    assert 1 < np.nanmax(det.max_significance) < 3
    assert det.n_regions == 0 and not det.mask.any()


def _doublet(s_strong, n=6000, b=15.0, err=0.02, multiplet="CIV"):
    """A noiseless thin doublet at the shape of the filter of b, its
    stronger line at S = s_strong, the weaker scaled by f lambda."""
    v, x = _grid(n, lam0=450.0)
    lags, ratios = multiplet_lines(multiplet)
    sigma = filter_sigma_kms(np.array([R]), b)[0]
    norm = np.sqrt(np.sum(_gauss(v, v[n // 4], sigma) ** 2)) / err
    i0 = n // 4
    depth = np.zeros(n)
    for lag, r in zip(lags, ratios, strict=True):
        depth += s_strong / norm * r * _gauss(v, v[i0] + lag, sigma)
    return v, x, 1 - depth, np.full(n, err), i0, lags, ratios


def test_the_joint_significance_of_a_doublet():
    v, x, flux, err, i0, lags, ratios = _doublet(2.8)
    num, den = filter_sums(v, 1 - flux, 1 / err ** 2,
                           filter_sigma_kms(np.full(len(v), R), 15.0))
    j, s_lines = joint_significance(v, num, den, lags, ratios)
    # J = S_1 sqrt(1 + r^2) for lines at the shape of the filter
    r = min(ratios)
    assert j[i0] == pytest.approx(2.8 * np.sqrt(1 + r ** 2), rel=2e-3)
    assert s_lines[:, i0] == pytest.approx([2.8, 2.8 * r], rel=2e-3)
    # Consistent with the ratio: chi2 = sum S^2 - J^2 = 0
    assert multiplet_seeds(j, s_lines, 3.0, 1.0)[i0]
    # The example of section 5: lines at 3 and 1.5 sigma give about 3.35
    assert (3 + 0.5 * 1.5) / np.sqrt(1 + 0.25) == pytest.approx(3.354, abs=1e-3)


def test_a_doublet_below_the_threshold_in_each_line_is_found_by_its_joint_significance():
    v, x, flux, err, i0, lags, ratios = _doublet(2.8)
    alone = matched_filter_detection(x, flux, err, R, multiplets=())
    assert alone.n_regions == 0
    # Found jointly at the threshold of the design (3 sigma), not at the
    # default (4 sigma: J = 3.13)
    det = matched_filter_detection(x, flux, err, R, multiplets=("CIV",),
                                   multiplet_sigma=LINE_SIGMA)
    assert det.n_regions == 2
    assert matched_filter_detection(x, flux, err, R, multiplets=("CIV",)).n_regions == 0
    for lag in lags:
        i = int(np.argmin(np.abs(v - (v[i0] + lag))))
        assert det.mask[i] and det.lines[i]
    assert np.nanargmax(det.joint["CIV"][1]) == i0


def test_at_the_default_threshold_the_doublet_gives_its_region_to_the_weaker_line():
    # The stronger line at 5 sigma, the weaker at 2.5: J = 5.6 > 4
    assert MULTIPLET_SIGMA == 4.0
    v, x, flux, err, i0, lags, ratios = _doublet(5.0)
    alone = matched_filter_detection(x, flux, err, R, multiplets=())
    assert alone.n_regions == 1
    det = matched_filter_detection(x, flux, err, R, multiplets=("CIV",))
    assert det.n_regions == 2
    i_weak = int(np.argmin(np.abs(v - (v[i0] + lags[1]))))
    assert det.mask[i_weak] and not alone.mask[i_weak]


def test_z_max_bounds_the_multiplets():
    v, x, flux, err, i0, _, _ = _doublet(2.8)
    z = x[i0] / (ATOM_DATA["CIV_1548"]["wave"] / 10.0) - 1
    det = matched_filter_detection(x, flux, err, R, multiplets=("CIV",), z_max=z - 0.001)
    assert det.n_regions == 0
    assert np.all(np.isnan(det.joint["CIV"][:, i0:]))


def test_the_lines_of_a_multiplet_outside_the_spectrum():
    v, x, flux, err, i0, _, _ = _doublet(5.0)
    det = matched_filter_detection(x, flux, err, R, multiplets=("CIV",))
    lag = multiplet_lines("CIV")[0][1]
    assert np.all(np.isnan(det.joint["CIV"][:, v + lag > v[-1]]))


def test_invalid_pixels_are_left_out():
    v, x = _grid(4000)
    sigma = filter_sigma_kms(np.array([R]), 15.0)[0]
    flux = 1 - 0.1 * _gauss(v, v[2000], sigma)
    err = np.full(len(x), 0.02)
    ref = matched_filter_detection(x, flux, err, R, multiplets=())
    flux_bad, err_bad = flux.copy(), err.copy()
    flux_bad[100:200] = np.nan
    err_bad[300:400] = 0.0
    flux_bad[3000:3100] = -50.0
    err_bad[3000:3100] = np.inf
    det = matched_filter_detection(x, flux_bad, err_bad, R, multiplets=())
    np.testing.assert_allclose(det.significance[:, 1500:2500], ref.significance[:, 1500:2500])
    assert np.all(np.isnan(det.significance[:, 150]))
    empty = matched_filter_detection(x, np.full(len(x), np.nan), err, R)
    assert np.all(np.isnan(empty.significance)) and empty.n_regions == 0


def test_the_resolution_is_required():
    _, x = _grid(100)
    with pytest.raises(ValueError, match="resolution"):
        matched_filter_detection(x, np.ones(100), np.full(100, 0.1), np.nan)


def test_detect_lines_on_a_spectrum():
    truth = MockTruth(components=[MockComponent("CIV", z=2.0, logN=13.0, b=10.0)],
                      lam_min_nm=455.0, lam_max_nm=470.0, dv_kms=DV, resol=R,
                      snr=50.0, z_em=2.5, seed=3)
    spec = generate(truth)
    det = detect_lines(spec)
    assert det.b_kms == DEFAULT_B_KMS
    x = spec.x.to_value("nm")
    cont = spec.cont.value
    ref = matched_filter_detection(x, spec.y.value / cont, spec.dy.value / cont,
                                   spec.resolution())
    np.testing.assert_array_equal(det.significance, ref.significance)
    for t in ("CIV_1548", "CIV_1550"):
        i = int(np.argmin(np.abs(x - ATOM_DATA[t]["wave"] / 10 * 3.0)))
        assert det.mask[i]


def test_a_strong_line_makes_no_multiplet_with_the_noise_at_its_lags():
    rng = np.random.default_rng(4)
    v, x = _grid(20000, lam0=450.0)
    sigma = filter_sigma_kms(np.array([R]), 15.0)[0]
    err = np.full(len(x), 0.02)
    flux = 1 - 0.5 * _gauss(v, v[10000], sigma) + rng.normal(0, 0.02, len(x))
    alone = matched_filter_detection(x, flux, err, R, multiplets=())
    det = matched_filter_detection(x, flux, err, R, multiplets=("CIV", "SiIV", "MgII"))
    # The joint significance is high at every lag of the strong line ...
    for name in ("CIV", "SiIV", "MgII"):
        for lag in multiplet_lines(name)[0]:
            i = int(np.argmin(np.abs(v - (v[10000] - lag))))
            assert np.nanmax(det.joint[name][:, i]) > 10
    # ... but its lines are not consistent with the ratio, and in its wings
    # J is not at a peak: no region at the lags (within 60 km/s)
    for name in ("CIV", "SiIV", "MgII"):
        lag = multiplet_lines(name)[0][1]
        for sign in (-1, 1):
            near = np.abs(v - v[10000] - sign * lag) < 60
            np.testing.assert_array_equal(det.mask[near], alone.mask[near])


def test_a_doublet_inconsistent_with_its_ratio_is_not_found_jointly():
    v, x, flux, err, i0, lags, ratios = _doublet(2.8)
    # The weaker line removed: J = 2.8 / sqrt(1 + r^2) < 3 anyway; the
    # stronger one at 2.9 and the weaker at 2.9 too (r = 0.5 expects 1.45):
    # J = 3.89, chi2 = 2 2.9^2 - 3.89^2 = 1.7, consistent within 3 sigma
    depth = 1 - flux
    depth_weak = np.where(v > v[i0] + lags[1] / 2, depth / ratios[1], depth)
    flux2 = 1 - depth_weak
    det = matched_filter_detection(x, flux2, err, R, multiplets=("CIV",),
                                   multiplet_sigma=LINE_SIGMA)
    assert det.n_regions == 2
    # A weaker line much stronger than its ratio (a blend, a saturated
    # doublet) is not a joint detection, but alone above 3 sigma it is found
    depth_blend = np.where(v > v[i0] + lags[1] / 2, depth / ratios[1] * 5, depth)
    det = matched_filter_detection(x, 1 - depth_blend, err, R, multiplets=("CIV",),
                                   multiplet_sigma=LINE_SIGMA)
    i_strong = i0
    assert not det.mask[i_strong]
    assert det.n_regions == 1


def test_on_pure_noise_the_doublets_add_regions_by_their_own_threshold():
    rng = np.random.default_rng(5)
    _, x = _grid(100000)
    err = np.full(len(x), 0.05)
    flux = 1 + rng.normal(0, 0.05, len(x))
    alone = matched_filter_detection(x, flux, err, R, multiplets=())
    det = matched_filter_detection(x, flux, err, R, multiplet_sigma=LINE_SIGMA)
    assert alone.n_regions < det.n_regions
    default = matched_filter_detection(x, flux, err, R)
    assert alone.n_regions <= default.n_regions < det.n_regions
    high = matched_filter_detection(x, flux, err, R, multiplet_sigma=np.inf)
    np.testing.assert_array_equal(high.mask, alone.mask)
