"""Edge cases of the pure functions of ``core/spectrum_operations.py``.

Written after mutation testing (``mutmut``, phase 1, item 6 of
``ROADMAP.md``): the snapshots exercise these functions on regular input,
but changing the handling of bad pixels, thresholds or fallbacks did not
make any test fail. These tests pin the current behaviour on degenerate
input: invalid pixels, windows as long as the array, exact thresholds,
no overlap, no smoothing. One of them records a V1/V2 difference found on
the way (F30, findings register): the running mean and standard deviation
pad the edges with zeros.
"""
import numpy as np
import astropy.units as au
import pytest

from astrocook.core.spectrum_operations import (
    compute_auto_fudge, compute_flux_scaling, interpolate_continuum_mask,
    rebin_spectrum, running_mean, running_std, smooth_spectrum,
)
from astrocook.core.structures import DataColumnV2, SpectrumDataV2


def _noisy(n=200, seed=0):
    rng = np.random.default_rng(seed)
    return 1.0 + 0.1 * rng.standard_normal(n)


# --- smooth_spectrum ------------------------------------------------------

def test_smooth_without_sigma_returns_the_input():
    x = np.linspace(600.0, 601.0, 200) * au.nm
    y = _noisy()
    assert smooth_spectrum(x, y, 0.0) is y
    assert smooth_spectrum(x, y, -5.0) is y


def test_smooth_small_sigma_still_smooths():
    # sigma below 1 km/s (pixels of ~2.5 km/s) is applied, not skipped
    x = np.linspace(600.0, 601.0, 200) * au.nm
    y = _noisy()
    out = smooth_spectrum(x, y, 0.5)
    assert not np.array_equal(out, y)
    assert np.std(out) < np.std(y)


def test_smooth_default_reference_redshift_is_zero():
    x = np.linspace(600.0, 601.0, 200) * au.nm
    y = _noisy()
    np.testing.assert_array_equal(smooth_spectrum(x, y, 20.0),
                                  smooth_spectrum(x, y, 20.0, z_ref=0.0))


# --- running_mean, running_std --------------------------------------------

@pytest.mark.parametrize("func, numpy_func", [(running_mean, np.mean),
                                              (running_std, np.nanstd)])
def test_running_window_longer_than_array(func, numpy_func):
    y = _noisy(11)
    np.testing.assert_allclose(func(y, 6), numpy_func(y))  # window 13 > 11


@pytest.mark.parametrize("func", [running_mean, running_std])
def test_running_window_as_long_as_array(func):
    # Window of 11 on 11 values: computed with the running kernel, not as a
    # single value for the whole array
    y = _noisy(11)
    out = func(y, 5)
    assert not np.allclose(out, out[0])


def test_running_std_of_a_constant_is_zero():
    # <y^2> - <y>^2 can be slightly negative: clipped to 0, not NaN
    y = np.full(50, 0.1)
    out = running_std(y, 3)
    assert np.all(np.isfinite(out))
    np.testing.assert_array_less(out[3:-3], 1e-6)


def test_running_mean_and_std_pad_the_edges_with_zeros():
    # Characterization of F30: np.convolve(mode='same') pads with zeros, so
    # the first and last h values are biased towards 0 (V1 repeats the first
    # and last full-window value instead)
    y = np.ones(50)
    h = 3
    mean = running_mean(y, h)
    np.testing.assert_allclose(mean[:h + 1], [4 / 7, 5 / 7, 6 / 7, 1.0])
    np.testing.assert_allclose(mean[h:-h], 1.0)
    assert running_std(y, h)[0] > 0.4


# --- compute_auto_fudge, interpolate_continuum_mask -----------------------

def test_auto_fudge_value():
    y = np.full(20, 1.1)
    model = np.full(20, 1.0)
    mask_abs = np.zeros(20, dtype=bool)
    assert compute_auto_fudge(y, mask_abs, model) == pytest.approx(1.1)


@pytest.mark.parametrize("n_unabs, expected", [(9, 1.0), (10, 1.1)])
def test_auto_fudge_needs_ten_unabsorbed_pixels(n_unabs, expected):
    y = np.full(20, 1.1)
    model = np.full(20, 1.0)
    mask_abs = np.arange(20) >= n_unabs
    assert compute_auto_fudge(y, mask_abs, model) == pytest.approx(expected)


@pytest.mark.parametrize("model_value", [0.0, np.nan])
def test_auto_fudge_without_a_model_level(model_value):
    y = np.full(20, 1.1)
    model = np.full(20, model_value)
    assert compute_auto_fudge(y, np.zeros(20, dtype=bool), model) == 1.0


def test_interpolate_continuum_with_two_points():
    x = np.arange(5.0)
    y = np.array([1.0, 0.0, 0.0, 0.0, 3.0])
    mask_abs = np.array([False, True, True, True, False])
    np.testing.assert_allclose(interpolate_continuum_mask(x, y, mask_abs),
                               [1.0, 1.5, 2.0, 2.5, 3.0])


def test_interpolate_continuum_with_one_point_raises():
    mask_abs = np.array([False, True, True])
    with pytest.raises(ValueError, match="^Cannot fit continuum: < 2 unabsorbed points found.$"):
        interpolate_continuum_mask(np.arange(3.0), np.ones(3), mask_abs)


# --- rebin_spectrum -------------------------------------------------------

def _spec_data(x, y, dy, half=0.0025):
    unit = au.dimensionless_unscaled
    return SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - half, au.nm),
        xmax=DataColumnV2(x + half, au.nm), y=DataColumnV2(y, unit),
        dy=DataColumnV2(dy, unit))


def _rebin(data, xstart=None, xend=None):
    out = rebin_spectrum(data.x.quantity, data.xmin.quantity, data.xmax.quantity,
                         data, xstart, xend, 0.02 * au.nm, au.nm, None, np.nan)
    return out[3].value, out[4].value


@pytest.mark.parametrize("bad", ["nan_y", "zero_dy", "negative_dy", "empty_bin"])
def test_rebin_ignores_bad_pixels(bad):
    # A bad pixel must not change the result: same as without it
    x = np.arange(600.0, 601.0, 0.005)
    y, dy = _noisy(len(x)), np.full(len(x), 0.1)
    data = _spec_data(x, y.copy(), dy.copy())
    i = 100
    if bad == "nan_y":
        data.y.values[i] = np.nan
    elif bad == "zero_dy":
        data.dy.values[i] = 0.0
    elif bad == "negative_dy":
        data.dy.values[i] = -0.1
    else:
        data.xmax.values[i] = data.xmin.values[i]
    keep = np.arange(len(x)) != i
    clean = SpectrumDataV2(
        x=DataColumnV2(x[keep], au.nm),
        xmin=DataColumnV2(data.xmin.values[keep], au.nm),
        xmax=DataColumnV2(data.xmax.values[keep], au.nm),
        y=DataColumnV2(y[keep], au.dimensionless_unscaled),
        dy=DataColumnV2(dy[keep], au.dimensionless_unscaled))
    start, end = 600.0 * au.nm, 600.995 * au.nm
    for got, expected in zip(_rebin(data, start, end), _rebin(clean, start, end)):
        np.testing.assert_allclose(got, expected, rtol=1e-12)


# --- compute_flux_scaling -------------------------------------------------

def _scaling_inputs(scale=2.0, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    x = np.linspace(600.0, 640.0, n)
    y_ref = scale * (x / 620.0) ** -1.5 * (1 + 0.01 * rng.standard_normal(n))
    y_in = (x / 620.0) ** -1.5 * (1 + 0.01 * rng.standard_normal(n))
    return x, y_ref, y_in


def test_flux_scaling_default_is_a_scalar():
    x, y_ref, y_in = _scaling_inputs()
    out = compute_flux_scaling(x, y_ref, x, y_in)
    assert np.all(out == out[0])
    assert out[0] == pytest.approx(2.0, rel=0.01)


@pytest.mark.parametrize("shift", [40.0, 50.0])  # touching, disjoint
def test_flux_scaling_without_overlap(shift):
    x, y_ref, y_in = _scaling_inputs()
    out = compute_flux_scaling(x, y_ref, x + shift, y_in)
    np.testing.assert_array_equal(out, np.ones_like(y_in))


def test_flux_scaling_ignores_bad_pixels():
    # Negative and NaN fluxes in the input are left out of the comparison
    x, y_ref, y_in = _scaling_inputs()
    spoiled = y_in.copy()
    spoiled[::7] = -1.0
    spoiled[3::7] = np.nan
    good = compute_flux_scaling(x, y_ref, x, y_in)[0]
    assert compute_flux_scaling(x, y_ref, x, spoiled)[0] == pytest.approx(good, rel=1e-3)


def test_flux_scaling_with_too_few_valid_pixels():
    x, y_ref, y_in = _scaling_inputs()
    y_in = np.full_like(y_in, -1.0)
    y_in[:9] = 1.0
    np.testing.assert_array_equal(compute_flux_scaling(x, y_ref, x, y_in),
                                  np.ones_like(y_in))
