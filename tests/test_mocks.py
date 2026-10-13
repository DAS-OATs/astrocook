"""Validation of the mock spectrum generator (astrocook.mocks).

Analytic checks where possible: equivalent width of an optically thin line,
conservation of the equivalent width by the convolution, width of the
instrumental profile, noise statistics, convergence with oversampling.
"""
import numpy as np
import pytest
from astropy import constants as const

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.mocks import (
    MockComponent, MockTruth, convolve_velocity, generate, make_grid,
    normalized_flux,
)
from astrocook.mocks.physics import C_KMS, FWHM_TO_SIGMA

Z = 2.97
CIV = ATOM_DATA["CIV_1548"]


def _truth(components, **kwargs):
    base = dict(lam_min_nm=612.0, lam_max_nm=618.0, dv_kms=2.5,
                resol=45000.0, snr=None, z_em=3.1, seed=1)
    base.update(kwargs)
    return MockTruth(components=components, **base)


def _ew_ang(truth):
    x, flux = normalized_flux(truth)
    _, xmin, xmax = make_grid(truth.lam_min_nm, truth.lam_max_nm, truth.dv_kms)
    return np.sum((1.0 - flux) * (xmax - xmin)) * 10.0


def test_grid_is_uniform_in_velocity():
    x, xmin, xmax = make_grid(600.0, 610.0, 2.5)
    np.testing.assert_allclose(np.diff(np.log(x)), np.log1p(2.5 / C_KMS), rtol=1e-8)
    np.testing.assert_allclose(xmin[1:], xmax[:-1], rtol=0, atol=0)
    assert np.all((x > xmin) & (x < xmax))


@pytest.mark.parametrize("resol", [None, 45000.0])
def test_thin_line_equivalent_width(resol):
    # Optically thin limit: W_rest = (pi e^2 / m_e c^2) f lambda_0^2 N
    logN = 11.0
    truth = _truth([MockComponent("CIV_1548", Z, logN, 15.0)], resol=resol)
    r_e = (const.e.esu**2 / (const.m_e * const.c**2)).cgs.value  # cm
    w_rest_cm = np.pi * r_e * CIV["f"] * (CIV["wave"] * 1e-8) ** 2 * 10**logN
    expected = w_rest_cm * 1e8 * (1 + Z)
    assert _ew_ang(truth) == pytest.approx(expected, rel=5e-3)


@pytest.mark.parametrize("logN", [12.5, 14.0, 15.0])
def test_convolution_conserves_equivalent_width(logN):
    comps = [MockComponent("CIV", Z, logN, 12.0)]
    unconvolved = _ew_ang(_truth(comps, resol=None))
    assert _ew_ang(_truth(comps, resol=20000.0)) == pytest.approx(unconvolved, rel=1e-3)


@pytest.mark.parametrize("resol", [45000.0, np.linspace(40000.0, 50000.0, 4001)])
def test_instrumental_profile_width(resol):
    # A delta function is turned into a Gaussian with FWHM = c / R
    dv = 0.5
    delta = np.zeros(4001)
    delta[2000] = 1.0
    out = convolve_velocity(1.0 - delta, dv, resol)
    profile = 1.0 - out
    v = (np.arange(4001) - 2000) * dv
    sigma = np.sqrt(np.sum(profile * v**2) / np.sum(profile))
    r_centre = resol if np.isscalar(resol) else resol[2000]
    assert sigma == pytest.approx(C_KMS / r_centre * FWHM_TO_SIGMA, rel=1e-2)


def test_equivalent_width_grows_with_column_density():
    ews = [_ew_ang(_truth([MockComponent("CIV_1548", Z, logN, 10.0)]))
           for logN in (12.0, 13.0, 14.0, 15.0, 16.0)]
    assert np.all(np.diff(ews) > 0)


def test_oversampling_converges_for_narrow_lines():
    # b = 1 km/s with 2.5 km/s pixels: the pixel average matters
    comps = [MockComponent("CIV_1548", Z, 13.0, 1.0)]
    ew = {k: _ew_ang(_truth(comps, resol=None, oversample=k)) for k in (1, 20, 40)}
    assert abs(ew[1] / ew[40] - 1) > 0.01            # centre sampling is biased
    assert ew[20] == pytest.approx(ew[40], rel=1e-3)  # converged


def test_noise_statistics_and_reproducibility():
    truth = _truth([], lam_min_nm=500.0, lam_max_nm=520.0, snr=20.0, seed=7,
                   cont_amplitude=3.0, cont_slope=-1.5)
    spec = generate(truth)
    pull = (spec.y.value - spec.get_column("cont").value) / spec.dy.value
    assert np.mean(pull) == pytest.approx(0.0, abs=0.03)
    assert np.std(pull) == pytest.approx(1.0, abs=0.03)
    np.testing.assert_array_equal(generate(truth).y.value, spec.y.value)
    other = MockTruth.from_dict({**truth.to_dict(), "seed": 8})
    assert not np.array_equal(generate(other).y.value, spec.y.value)


def test_spectrum_properties():
    spec = generate(_truth([MockComponent("CIV", Z, 13.5, 10.0)], snr=50.0))
    assert spec._data.z_em == 3.1
    assert spec._data.resol == 45000.0
    assert spec.meta["resol"] == 45000.0
    np.testing.assert_allclose(spec.get_column("mock_flux").value,
                               spec.get_column("cont").value * normalized_flux(
                                   _truth([MockComponent("CIV", Z, 13.5, 10.0)]))[1])


def test_truth_json_round_trip(tmp_path):
    truth = _truth([MockComponent("CIV", Z, 13.5, 10.0, btur=2.0),
                    MockComponent("Ly_a", 4.05, 14.0, 25.0)], snr=30.0)
    path = tmp_path / "truth.json"
    truth.to_json(str(path))
    assert MockTruth.from_json(str(path)) == truth


def test_agrees_with_v2_model_for_resolved_lines():
    # Cross-check with the V2 fitting code for well-sampled lines (b well
    # above the pixel size); the two implementations are independent
    from astrocook.core.structures import ComponentDataV2, SystemListDataV2
    from astrocook.core.system_list import SystemListV2
    from astrocook.fitting.voigt_fitter import VoigtFitterV2

    truth = _truth([MockComponent("CIV", Z, 13.6, 12.0)], oversample=10)
    spec = generate(truth)
    comp = ComponentDataV2(id=1, z=Z, dz=None, logN=13.6, dlogN=None, b=12.0,
                           db=None, series="CIV")
    _, v2 = VoigtFitterV2(
        spec, SystemListV2(SystemListDataV2(components=[comp]))).compute_model_flux()
    np.testing.assert_allclose(v2, normalized_flux(truth)[1], atol=3e-3)
