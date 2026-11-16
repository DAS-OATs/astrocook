"""The cases of the mocks added for the robustness grid (phase 7, PR 7a.1):
BAL troughs, emission lines with shifts and asymmetries, tellurics, the
proximity zone and the associated absorbers, the truth map."""
import numpy as np
import pytest

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.mocks import (
    CATEGORIES, LYMAN_SERIES, QSO_EMISSION_LINES_FULL, BALTrough,
    EmissionLine, MockComponent, MockTruth, TelluricModel,
    associated_absorbers, bal_optical_depth, component_category, continuum,
    generate, lya_forest, make_grid, normalized_flux, proximity_omega,
    telluric_optical_depth, truth_map,
)
from astrocook.mocks import telluric as tel
from astrocook.mocks.physics import C_KMS, outflow_velocity, trough_profile
from astrocook.mocks.truthmap import BITS

Z_EM = 3.5


def _truth(components=(), **kwargs):
    base = dict(lam_min_nm=600.0, lam_max_nm=720.0, dv_kms=5.0, resol=30000.0,
                snr=None, z_em=Z_EM, seed=1, oversample=2)
    base.update(kwargs)
    return MockTruth(components=list(components), **base)


# --- Defaults: the mocks made before are unchanged ----------------------

def test_new_fields_are_off_by_default():
    truth = _truth([MockComponent("CIV", 3.3, 13.5, 10.0)])
    assert truth.bal_troughs == [] and truth.telluric is None
    assert truth.proximity_kms is None
    line = EmissionLine(121.567, 2.0, 4000.0)
    assert line.shift_kms == 0.0 and line.asymmetry == 0.0


def test_old_truth_dict_still_loads():
    truth = _truth([MockComponent("CIV", 3.3, 13.5, 10.0)],
                   emission_lines=[EmissionLine(154.9, 0.8, 5000.0)])
    old = truth.to_dict()
    for key in ("bal_troughs", "telluric", "proximity_kms"):
        del old[key]
    for line in old["emission_lines"]:
        del line["shift_kms"], line["asymmetry"]
    assert MockTruth.from_dict(old) == truth


def test_round_trip_with_the_new_fields(tmp_path):
    truth = _truth(bal_troughs=[BALTrough(3000.0, 9000.0, 2.0, column_ratios=(1, 0.3, 1, 2))],
                   telluric=TelluricModel(1.5, 4.0, 1.0), proximity_kms=2000.0,
                   emission_lines=QSO_EMISSION_LINES_FULL)
    path = tmp_path / "truth.json"
    truth.to_json(str(path))
    assert MockTruth.from_json(str(path)) == truth


# --- BAL troughs ---------------------------------------------------------

def test_outflow_velocity():
    lam0 = ATOM_DATA["CIV_1548"]["wave"]
    beta = 0.03
    lam = lam0 * (1 + Z_EM) * np.sqrt((1 - beta) / (1 + beta))
    assert outflow_velocity(np.array([lam]), lam0, Z_EM)[0] == pytest.approx(beta * C_KMS)
    assert outflow_velocity(np.array([lam0 * (1 + Z_EM)]), lam0, Z_EM)[0] == pytest.approx(0)


def test_trough_profile():
    v = np.array([0.0, 2000.0, 2500.0, 3000.0, 5000.0, 7000.0, 8000.0, 9000.0])
    p = trough_profile(v, 2000.0, 8000.0, 1000.0)
    assert p[[0, 1, 6, 7]] == pytest.approx(0.0)
    assert p[2] == pytest.approx(0.5)
    assert p[[3, 4, 5]] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        trough_profile(v, 3000.0, 2000.0, 100.0)


def _tau_at(line, v_kms, trough):
    lam0 = ATOM_DATA[line]["wave"]
    beta = v_kms / C_KMS
    lam = lam0 * (1 + Z_EM) * np.sqrt((1 - beta) / (1 + beta))
    return lam, bal_optical_depth(np.array([lam]), Z_EM, trough)[0]


def test_bal_optical_depth_in_the_ratio_of_f_lambda():
    # A trough far from the others, so that the lines do not overlap at
    # the velocity probed: only CIV 1548 and SiIV 1394
    trough = BALTrough(5000.0, 6000.0, 2.0, edge_kms=100.0, ions=("CIV_1548", "SiIV_1393"))
    _, tau_civ = _tau_at("CIV_1548", 5500.0, trough)
    _, tau_si = _tau_at("SiIV_1393", 5500.0, trough)
    assert tau_civ == pytest.approx(2.0)
    f_lam = {k: ATOM_DATA[k]["f"] * ATOM_DATA[k]["wave"] for k in ("CIV_1548", "SiIV_1393")}
    assert tau_si / tau_civ == pytest.approx(f_lam["SiIV_1393"] / f_lam["CIV_1548"])
    halved = BALTrough(5000.0, 6000.0, 2.0, edge_kms=100.0,
                       ions=("CIV_1548", "SiIV_1393"), column_ratios=(1.0, 0.5))
    assert _tau_at("SiIV_1393", 5500.0, halved)[1] == pytest.approx(tau_si / 2)
    with pytest.raises(ValueError):
        bal_optical_depth(np.array([600.0]), Z_EM, BALTrough(1.0, 2.0, 1.0, column_ratios=(1.0,)))


def test_bal_trough_in_the_mock():
    trough = BALTrough(3000.0, 8000.0, 3.0)
    x, flux = normalized_flux(_truth(bal_troughs=[trough], lam_min_nm=680.0,
                                     lam_max_nm=700.0))
    v = outflow_velocity(x * 10.0, ATOM_DATA["CIV_1548"]["wave"], Z_EM)
    # Saturated on the flat part, untouched outside the troughs of CIV
    # (the SiIV one is far to the blue)
    assert np.all(flux[(v > 4500) & (v < 6000)] < 0.1)
    assert np.allclose(flux[(v < 1500) & (v > -2000)], 1.0)


# --- Emission lines ------------------------------------------------------

def test_emission_line_shift_and_asymmetry():
    line = EmissionLine(154.906, 1.0, 5000.0, shift_kms=-1000.0, asymmetry=0.3)
    truth = _truth(emission_lines=[line])
    x, _, _ = make_grid(truth.lam_min_nm, truth.lam_max_nm, 1.0)
    cont = continuum(truth, x)
    peak = x[np.argmax(cont)]
    centre = 154.906 * (1 + Z_EM) * (1 - 1000.0 / C_KMS)
    assert (peak - centre) / centre * C_KMS == pytest.approx(0.0, abs=1.0)
    # Half maximum on each side: the blue half is 1.3, the red 0.7 times
    # the half FWHM
    half = x[cont >= 1.5]
    blue = (centre - half[0]) / centre * C_KMS
    red = (half[-1] - centre) / centre * C_KMS
    assert blue == pytest.approx(2500.0 * 1.3, rel=0.01)
    assert red == pytest.approx(2500.0 * 0.7, rel=0.01)
    with pytest.raises(ValueError):
        continuum(_truth(emission_lines=[EmissionLine(154.9, 1.0, 5000.0, asymmetry=1.0)]), x)


def test_full_set_of_emission_lines():
    rest = [line.rest_nm for line in QSO_EMISSION_LINES_FULL]
    assert rest == sorted(rest) and len(rest) == 7
    civ = [line for line in QSO_EMISSION_LINES_FULL if abs(line.rest_nm - 154.9) < 0.1]
    assert civ[0].shift_kms < 0 and civ[0].asymmetry > 0


# --- Tellurics -----------------------------------------------------------

def test_telluric_template():
    lam, tau, meta = tel.template()
    assert lam[0] == pytest.approx(550.0) and lam[-1] == pytest.approx(1000.0)
    assert np.all(tau >= 0)
    assert meta["airmass"] == tel.TEMPLATE_AIRMASS and meta["pwv_mm"] == tel.TEMPLATE_PWV_MM
    # The A band is saturated at its core, the window around 600 nm clear
    a_band = (lam > 759) & (lam < 772)
    assert tau[a_band].max() > 10
    assert np.exp(-tau[(lam > 600) & (lam < 610)]).min() > 0.9


def test_telluric_scaling():
    lam = np.array([690.0, 694.0, 760.5, 763.0, 720.0, 723.0, 930.0, 940.0])
    ref = telluric_optical_depth(lam, TelluricModel())
    o2 = tel.is_o2(lam)
    assert o2.any() and (~o2).any()
    double_airmass = telluric_optical_depth(lam, TelluricModel(airmass=2.2))
    assert double_airmass == pytest.approx(2 * ref)
    double_pwv = telluric_optical_depth(lam, TelluricModel(pwv_mm=5.0))
    assert double_pwv[o2] == pytest.approx(ref[o2])
    assert double_pwv[~o2] == pytest.approx(2 * ref[~o2])
    # Outside the template: no absorption
    assert telluric_optical_depth(np.array([500.0, 1100.0]), TelluricModel()) == pytest.approx(0)


def test_telluric_shift():
    lam_t, tau_t, _ = tel.template()
    i = np.argmax(tau_t * ((lam_t > 760) & (lam_t < 770)))
    shift = 2.0
    lam = lam_t[i] * (1 + shift / C_KMS)
    assert telluric_optical_depth(np.array([lam]), TelluricModel(shift_kms=shift))[0] == \
        pytest.approx(tau_t[i], rel=1e-6)


def test_tellurics_before_the_convolution():
    kwargs = dict(lam_min_nm=758.0, lam_max_nm=772.0, z_em=4.5,
                  telluric=TelluricModel())
    _, sharp = normalized_flux(_truth(resol=None, **kwargs))
    _, smooth = normalized_flux(_truth(resol=8000.0, **kwargs))
    assert sharp.min() < 1e-3 < smooth.min()
    # The convolution conserves the absorbed area (away from the edges)
    assert np.sum(1 - sharp) == pytest.approx(np.sum(1 - smooth), rel=0.02)


# --- Proximity zone and associated absorbers -----------------------------

def test_proximity_omega():
    z = Z_EM - np.array([1000.0, 2000.0, 4000.0]) / C_KMS * (1 + Z_EM)
    assert proximity_omega(z, Z_EM, 2000.0) == pytest.approx([4.0, 1.0, 0.25])
    assert proximity_omega(np.array([Z_EM]), Z_EM, 2000.0)[0] == np.inf


def test_proximity_zone_of_the_forest():
    def forest(**kw):
        return lya_forest(np.random.default_rng(3), 3.3, Z_EM, **kw)
    plain = forest()
    near = forest(z_em=Z_EM, proximity_kms=2000.0)
    assert len(plain) == len(near)
    assert [c.z for c in plain] == [c.z for c in near]
    assert [c.b for c in plain] == [c.b for c in near]
    dlogN = np.array([p.logN - n.logN for p, n in zip(plain, near, strict=True)])
    omega = proximity_omega(np.array([c.z for c in plain]), Z_EM, 2000.0)
    assert dlogN == pytest.approx(np.log10(1 + omega))
    with pytest.raises(ValueError):
        forest(proximity_kms=2000.0)


def test_proximity_zone_lowers_the_density_of_lines():
    rng_args = dict(z_em=Z_EM, proximity_kms=3000.0)
    z_lo = Z_EM - 0.02
    plain = lya_forest(np.random.default_rng(4), 3.0, Z_EM)
    near = lya_forest(np.random.default_rng(4), 3.0, Z_EM, **rng_args)
    def strong(comps, lo, hi):
        return sum(1 for c in comps if lo <= c.z < hi and c.logN > 13)
    assert strong(near, z_lo, Z_EM) < 0.6 * strong(plain, z_lo, Z_EM)
    assert strong(near, 3.0, 3.1) >= 0.95 * strong(plain, 3.0, 3.1)


def test_associated_absorbers():
    rng = np.random.default_rng(5)
    comps = associated_absorbers(rng, 20, Z_EM, v_range_kms=(3500.0, 6000.0))
    metals = [c for c in comps if c.series != LYMAN_SERIES]
    assert {c.series for c in metals} <= {"CIV", "NV", "SiIV"}
    v = [(c.z - Z_EM) / (1 + Z_EM) * C_KMS for c in metals]
    assert min(v) > 3500.0 - 60 and max(v) < 6000.0 + 60
    assert sum(c.series == LYMAN_SERIES for c in comps) == 20


# --- Truth map -----------------------------------------------------------

def test_component_category():
    assert component_category(MockComponent("Ly_ab", 3.0, 14.0, 25.0)) == "forest"
    assert component_category(MockComponent(LYMAN_SERIES, 3.0, 20.5, 25.0)) == "strong_hi"
    assert component_category(MockComponent("CIV", 3.0, 13.0, 10.0)) == "metal"


def test_truth_map_of_one_category_is_the_flux():
    truth = _truth([MockComponent("CIV", 3.3, 13.5, 10.0)], lam_min_nm=660.0,
                   lam_max_nm=670.0)
    tmap = truth_map(truth)
    _, flux = normalized_flux(truth)
    assert tmap.depth["metal"] == pytest.approx(1 - flux)
    for name in ("forest", "strong_hi", "bal", "telluric"):
        assert not tmap.depth[name].any()
    assert set(np.unique(tmap.mask())) == {0, BITS["metal"]}


def test_truth_map_of_every_category():
    comps = [MockComponent("Ly_a", 3.0, 14.0, 25.0),          # 486.7 nm
             MockComponent("Ly_a", 3.02, 20.5, 25.0),         # strong
             MockComponent("CIV", 2.2, 14.0, 10.0)]           # 495.5 nm
    truth = _truth(comps, lam_min_nm=480.0, lam_max_nm=700.0, lyman_limit=True,
                   bal_troughs=[BALTrough(30000.0, 35000.0, 2.0, ions=("CIV",))],
                   telluric=TelluricModel())
    tmap = truth_map(truth)
    assert set(tmap.depth) == set(CATEGORIES)
    for name in CATEGORIES:
        assert tmap.absorbed(name).any(), name
    # Where one category absorbs alone, the flux is its transmission
    _, flux = normalized_flux(truth)
    mask = tmap.mask(1e-6)
    for name in CATEGORIES:
        alone = mask == BITS[name]
        if alone.any():
            assert flux[alone] == pytest.approx(1 - tmap.depth[name][alone], abs=1e-5)


def test_generate_with_every_case():
    truth = _truth([MockComponent("CIV", 3.3, 13.5, 10.0)], snr=30.0,
                   bal_troughs=[BALTrough(3000.0, 8000.0, 1.0)],
                   telluric=TelluricModel(1.3, 3.0, 0.5), proximity_kms=2000.0,
                   emission_lines=QSO_EMISSION_LINES_FULL)
    spec = generate(truth)
    assert np.all(np.isfinite(spec.y.value))
    assert spec.get_column("mock_flux").value.min() < 0.5 * spec.get_column("cont").value.max()
