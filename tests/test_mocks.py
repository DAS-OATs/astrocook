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


# --- Populations, emission lines, spikes ------------------------------------

from astrocook.mocks import (  # noqa: E402
    QSO_EMISSION_LINES, EmissionLine, continuum, expected_forest_lines,
    lya_forest, mean_transmission_fg08, metal_systems, spike_indices,
)


def test_forest_number_of_lines():
    rng = np.random.default_rng(3)
    comps = lya_forest(rng, 2.5, 4.0)
    expected = expected_forest_lines(2.5, 4.0)
    assert abs(len(comps) - expected) < 4 * np.sqrt(expected)
    assert all(c.series == "Ly_a" and 2.5 <= c.z <= 4.0 for c in comps)


def test_forest_column_density_distribution():
    # Fraction above logN = 13 for f(N) ~ N^-1.5 between 10^12 and 10^17
    comps = lya_forest(np.random.default_rng(4), 2.0, 5.0)
    logN = np.array([c.logN for c in comps])
    e = 1.0 - 1.5
    p = (1e13**e - 1e17**e) / (1e12**e - 1e17**e)
    n = len(logN)
    assert abs(np.mean(logN > 13.0) - p) < 4 * np.sqrt(p * (1 - p) / n)
    assert logN.min() >= 12.0 and logN.max() <= 17.0


def test_forest_doppler_parameters():
    comps = lya_forest(np.random.default_rng(5), 2.0, 5.0)
    b = np.array([c.b for c in comps])
    assert b.min() >= 10.0
    assert np.median(b) == pytest.approx(25.0, rel=0.03)


def test_forest_redshift_distribution():
    # p(z) ~ (1+z)^2.5: fraction of lines in the upper half of the range
    comps = lya_forest(np.random.default_rng(6), 2.0, 4.0)
    z = np.array([c.z for c in comps])
    p = (5.0**3.5 - 4.0**3.5) / (5.0**3.5 - 3.0**3.5)
    assert abs(np.mean(z > 3.0) - p) < 4 * np.sqrt(p * (1 - p) / len(z))


@pytest.mark.parametrize("z", [3.0, 3.8])
def test_forest_mean_transmission(z):
    # Calibration of the default line density against the observed mean
    # transmission (Faucher-Giguere et al. 2008), averaged over 3 seeds
    lam = 121.567 * (1 + np.array([z - 0.1, z + 0.1]))
    fluxes = []
    for seed in range(3):
        comps = lya_forest(np.random.default_rng(seed), z - 0.12, z + 0.12)
        truth = MockTruth(components=comps, lam_min_nm=lam[0], lam_max_nm=lam[1],
                          dv_kms=5.0, resol=None, snr=None, oversample=2)
        fluxes.append(np.mean(normalized_flux(truth)[1]))
    assert np.mean(fluxes) == pytest.approx(mean_transmission_fg08(z), abs=0.06)


def test_metal_systems():
    comps = metal_systems(np.random.default_rng(7), 20, 2.0, 3.0,
                          species=("CIV", "MgII"), dv_kms=50.0)
    assert 20 <= len(comps) <= 80
    assert {c.series for c in comps} <= {"CIV", "MgII"}
    assert all(1.99 < c.z < 3.01 and 4.0 <= c.b <= 15.0 for c in comps)


def test_emission_lines_in_the_continuum():
    truth = _truth([], lam_min_nm=580.0, lam_max_nm=640.0, z_em=4.0,
                   emission_lines=[EmissionLine(121.567, 2.0, 4000.0)])
    x, _, _ = make_grid(580.0, 640.0, 2.5)
    cont = continuum(truth, x)
    peak = np.argmin(np.abs(x - 121.567 * 5.0))
    assert cont[peak] == pytest.approx(3.0, rel=1e-3)
    assert cont[0] == pytest.approx(1.0, abs=1e-3)


def test_spikes():
    truth = _truth([], snr=20.0, n_spikes=5, spike_amplitude=10.0, seed=2)
    spec = generate(truth)
    idx = spike_indices(truth, len(spec.x.value))
    assert len(idx) == 5
    pull = (spec.y.value - spec.get_column("mock_flux").value) / spec.dy.value
    assert np.all(pull[idx] > 6.0)
    others = np.delete(pull, idx)
    assert np.max(np.abs(others)) < 6.0
    np.testing.assert_array_equal(spike_indices(truth, len(spec.x.value)), idx)


def test_full_mock_spectrum(tmp_path):
    # Forest, metals and emission lines together: the forest absorbs only
    # blueward of the Lyman-alpha emission
    rng = np.random.default_rng(11)
    z_em = 4.0
    lya_em = 121.567 * (1 + z_em)
    comps = (lya_forest(rng, 3.6, z_em)
             + metal_systems(rng, 5, 2.8, 3.3, species=("CIV",)))
    truth = MockTruth(components=comps, lam_min_nm=560.0, lam_max_nm=640.0,
                      dv_kms=5.0, resol=40000.0, snr=30.0, z_em=z_em, seed=11,
                      oversample=3, emission_lines=QSO_EMISSION_LINES)
    spec = generate(truth)
    x = spec.x.value
    norm = spec.get_column("mock_flux").value / spec.get_column("cont").value
    assert np.mean(norm[x < lya_em - 2]) < np.mean(norm[x > lya_em + 2]) - 0.1
    path = tmp_path / "truth.json"
    truth.to_json(str(path))
    assert MockTruth.from_json(str(path)) == truth


# --- Strong H I absorbers and the Lyman limit (phase 4, validation) --------

from scipy.integrate import quad  # noqa: E402

from astrocook.mocks import (  # noqa: E402
    DLA_LOGN, LLS_LOGN, LYMAN_SERIES, strong_hi_systems,
)
from astrocook.mocks.physics import (  # noqa: E402
    LYMAN_LIMIT_ANG, SIGMA_LYMAN_LIMIT, is_hydrogen,
)

Z_LLS = 3.0


def _lls_truth(logN, lyman_limit):
    # 880-1000 A rest frame, noiseless, no instrumental profile
    return MockTruth(
        components=[MockComponent(LYMAN_SERIES, Z_LLS, logN, 20.0)],
        lam_min_nm=88.0 * (1 + Z_LLS), lam_max_nm=100.0 * (1 + Z_LLS),
        dv_kms=5.0, resol=None, snr=None, lyman_limit=lyman_limit)


def _rest_ang(x_nm):
    return x_nm * 10.0 / (1 + Z_LLS)


def test_lyman_series_covers_the_table():
    names = LYMAN_SERIES.split(",")
    assert names[0] == "Ly_a" and names[-1] == "Ly_15" and "Ly_lim" not in names
    assert is_hydrogen(LYMAN_SERIES) and is_hydrogen("Ly_a")
    assert not is_hydrogen("CIV")


def test_lyman_limit_transmission():
    # Below the edge (no lines there) the flux is exp(-N sigma (l/l_LL)^3)
    logN = 17.5  # tau = 2 at the edge
    x, flux = normalized_flux(_lls_truth(logN, lyman_limit=True))
    lam = _rest_ang(x)
    below = (lam > 885.0) & (lam < 908.0)
    expected = np.exp(-10**logN * SIGMA_LYMAN_LIMIT * (lam / LYMAN_LIMIT_ANG) ** 3)
    assert flux[below] == pytest.approx(expected[below], rel=1e-3)
    assert expected[below].min() < 0.2


def test_lyman_limit_is_off_by_default():
    _, with_edge = normalized_flux(_lls_truth(17.5, lyman_limit=True))
    x, without = normalized_flux(_lls_truth(17.5, lyman_limit=False))
    lam = _rest_ang(x)
    assert MockTruth(components=[], lam_min_nm=1, lam_max_nm=2,
                     dv_kms=1).lyman_limit is False
    assert without[(lam > 885.0) & (lam < 908.0)] == pytest.approx(1.0, abs=1e-3)
    redwards = lam > LYMAN_LIMIT_ANG + 1.0
    assert with_edge[redwards] == pytest.approx(without[redwards])


def test_dla_damping_wings():
    # Equivalent width of a DLA in Ly-a against the Lorentzian wings,
    # tau = C / dlambda^2 with C = N (pi e^2 / m_e c) f gamma lambda0^4 /
    # (4 pi^2 c^2), integrated over the same window
    logN, half = 20.5, 60.0  # half-width of the window, Angstrom (rest)
    lya = ATOM_DATA["Ly_a"]
    lam0, z = lya["wave"], 3.0
    truth = MockTruth(components=[MockComponent("Ly_a", z, logN, 20.0)],
                      lam_min_nm=(lam0 - half) / 10 * (1 + z),
                      lam_max_nm=(lam0 + half) / 10 * (1 + z),
                      dv_kms=5.0, resol=None, snr=None)
    w_mock = _ew_ang(truth) / (1 + z)
    pi_e2_mec = (np.pi * const.e.esu**2 / (const.m_e * const.c)).cgs.value
    c_cms = const.c.cgs.value
    lam0_cm = lam0 * 1e-8
    c_wing = (10**logN * pi_e2_mec * lya["f"] * lya["gamma"] * lam0_cm**4
              / (4 * np.pi**2 * c_cms**2)) * 1e16  # cm^2 -> Angstrom^2
    w_wing = 2 * quad(lambda d: 1 - np.exp(-c_wing / d**2), 0, half)[0]
    assert w_mock == pytest.approx(w_wing, rel=0.02)
    assert w_wing > 10.0  # damped regime


def test_strong_hi_systems():
    rng = np.random.default_rng(3)
    comps = strong_hi_systems(rng, 200, 2.0, 3.0, logN_range=DLA_LOGN)
    assert len(comps) == 200
    assert all(c.series == LYMAN_SERIES for c in comps)
    logN = np.array([c.logN for c in comps])
    assert DLA_LOGN[0] <= logN.min() and logN.max() <= DLA_LOGN[1]
    assert np.mean(logN < 21.0) > 0.5  # f(N) falls with N
    assert [c.z for c in comps] == sorted(c.z for c in comps)
    again = strong_hi_systems(np.random.default_rng(3), 200, 2.0, 3.0,
                              logN_range=DLA_LOGN)
    assert again == comps
    lls = strong_hi_systems(rng, 50, 2.0, 3.0, logN_range=LLS_LOGN)
    assert all(LLS_LOGN[0] <= c.logN <= LLS_LOGN[1] for c in lls)


def test_metal_systems_with_hi():
    def draw(**kwargs):
        return metal_systems(np.random.default_rng(7), 20, 2.0, 3.0,
                             species=("CIV", "MgII"), dv_kms=50.0, **kwargs)
    plain, with_hi = draw(), draw(hi_logN_range=(14.0, 16.0))
    hi = [c for c in with_hi if c.series == LYMAN_SERIES]
    metals = [c for c in with_hi if c.series != LYMAN_SERIES]
    assert len(hi) == 20
    assert all(14.0 <= c.logN <= 16.0 and 15.0 <= c.b <= 35.0 for c in hi)
    for c in hi:  # each H I sits at its system, within the metal spread
        assert min(abs(m.z - c.z) / (1 + c.z) * C_KMS for m in metals) <= 50.0
    assert not any(c.series == LYMAN_SERIES for c in plain)


def test_forest_series():
    rng = np.random.default_rng(1)
    comps = lya_forest(rng, 2.9, 3.0, series=LYMAN_SERIES)
    assert comps and all(c.series == LYMAN_SERIES for c in comps)


def test_generate_without_truth_columns():
    truth = _truth([MockComponent("CIV", Z, 13.5, 10.0)], snr=30.0)
    observed = generate(truth, truth_columns=False)
    full = generate(truth)
    assert not observed.has_aux_column("cont")
    assert not observed.has_aux_column("mock_flux")
    assert observed.y.value == pytest.approx(full.y.value)
    assert observed.dy.value == pytest.approx(full.dy.value)
    assert observed._data.z_em == full._data.z_em


def test_truth_round_trip_with_and_without_lyman_limit():
    truth = _lls_truth(17.5, lyman_limit=True)
    assert MockTruth.from_dict(truth.to_dict()) == truth
    old = truth.to_dict()
    del old["lyman_limit"]  # a truth saved before the option existed
    assert MockTruth.from_dict(old).lyman_limit is False


def test_noise_follows_the_flux():
    # sigma = (cont / snr) sqrt(a + (1 - a) F / cont): the continuum value
    # where unabsorbed, sqrt(a) of it in a black core
    snr, floor = 20.0, 0.1
    comps = [MockComponent("CIV", Z, 15.5, 10.0)]
    flat = generate(_truth(comps, snr=snr))
    truth = _truth(comps, snr=snr, noise_floor=floor)
    spec = generate(truth)
    cont = spec.get_column("cont").value
    true_flux = spec.get_column("mock_flux").value
    dy = spec.dy.value
    assert flat.dy.value == pytest.approx(cont / snr)  # default: constant
    clean = true_flux / cont > 0.999
    black = true_flux / cont < 1e-3
    assert clean.any() and black.any()
    assert dy[clean] == pytest.approx(cont[clean] / snr, rel=1e-3)
    assert dy[black] == pytest.approx(cont[black] / snr * np.sqrt(floor), rel=1e-2)
    # Same random draws, scaled: the residuals in units of the error match
    pulls = (spec.y.value - true_flux) / dy
    assert pulls == pytest.approx((flat.y.value - true_flux) / flat.dy.value)


@pytest.mark.parametrize("floor", [-0.1, 1.5])
def test_noise_floor_out_of_range(floor):
    with pytest.raises(ValueError, match="noise_floor"):
        generate(_truth([], snr=20.0, noise_floor=floor))


def test_forest_strong_lines_get_their_series():
    plain = lya_forest(np.random.default_rng(5), 2.8, 3.0)
    mixed = lya_forest(np.random.default_rng(5), 2.8, 3.0,
                       strong_series=LYMAN_SERIES, strong_logN=15.0)
    assert [(c.z, c.logN, c.b) for c in mixed] == [(c.z, c.logN, c.b) for c in plain]
    strong = [c for c in mixed if c.logN >= 15.0]
    assert strong and all(c.series == LYMAN_SERIES for c in strong)
    assert all(c.series == "Ly_a" for c in mixed if c.logN < 15.0)
    assert len(strong) < 0.1 * len(mixed)
