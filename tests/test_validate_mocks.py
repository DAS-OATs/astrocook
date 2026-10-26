"""The validation benchmark on mocks (benchmarks/validate_mocks.py, phase 4):
its mocks, its detection limit and its metrics, on cases with a known
answer, and a run on the smallest grid."""
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy import constants as const

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "benchmarks"))

import validate_mocks as vm  # noqa: E402

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS  # noqa: E402
from astrocook.core.structures import ComponentDataV2  # noqa: E402
from astrocook.mocks import LYMAN_SERIES, MockComponent, generate  # noqa: E402
from astrocook.mocks.physics import C_KMS  # noqa: E402


def _z_offset(z, dv_kms):
    return z + dv_kms / C_KMS * (1 + z)


# --- Mocks ------------------------------------------------------------------

@pytest.mark.parametrize("resol", [8900.0, 45000.0, 140000.0])
def test_mock_sampling(resol):
    truth, _, _ = vm.make_truth(30.0, resol, 3.5, 2)
    fwhm = C_KMS / resol
    assert truth.dv_kms == pytest.approx(fwhm / vm.PIXELS_PER_FWHM)
    assert truth.dv_kms / truth.oversample <= vm.MAX_SUBPIXEL_KMS
    assert truth.lam_min_nm == pytest.approx(104.0 * 4.5)
    assert truth.lyman_limit and truth.noise_floor == vm.NOISE_FLOOR


def test_metal_systems_on_the_red_side():
    z_em = 3.5
    truth, systems, _ = vm.make_truth(30.0, 45000.0, z_em, 2)
    assert {s["species"]: 0 for s in systems}.keys() == set(vm.METALS)
    assert len(systems) == sum(vm.METAL_SYSTEMS.values())
    lya_obs = vm.LYA_NM * (1 + z_em)
    for s in systems:
        for c in s["components"]:
            for t in STANDARD_MULTIPLETS[s["species"]]:
                lam = ATOM_DATA[t]["wave"] / 10 * (1 + c.z)
                assert lya_obs < lam < truth.lam_max_nm, (s["species"], c.z)
    # Their H I is in the truth too
    hi = [c for c in truth.components if c.series == LYMAN_SERIES and c.logN < 17]
    assert len(hi) >= len(systems)


def test_strong_absorbers_in_seeds_0_and_1():
    for seed, kind, (lo, hi) in [(0, "DLA", (20.3, 22.0)), (1, "subDLA", (19.0, 20.3))]:
        truth, _, strong = vm.make_truth(30.0, 45000.0, 3.0, seed)
        assert strong == kind
        assert any(c.series == LYMAN_SERIES and lo <= c.logN <= hi
                   for c in truth.components)
    assert vm.make_truth(30.0, 45000.0, 3.0, 2)[2] is None


def test_mocks_are_reproducible():
    a = vm.make_truth(10.0, 8900.0, 2.5, 3)
    b = vm.make_truth(10.0, 8900.0, 2.5, 3)
    assert a[0] == b[0]
    assert vm.make_truth(10.0, 8900.0, 2.5, 4)[0] != a[0]


# --- Detection limit -----------------------------------------------------

def test_weak_line():
    assert vm.weak_line("CIV") == "CIV_1550"
    assert vm.weak_line("MgII") == "MgII_2803"


def test_equivalent_width_of_a_thin_line():
    # W = (pi e^2 / m_e c^2) N f lambda^2 for an optically thin line
    line, z, logN = "CIV_1550", 2.0, 11.0
    atom = ATOM_DATA[line]
    r_e = (const.e.esu**2 / (const.m_e * const.c**2)).cgs.value
    lam_cm = atom["wave"] * 1e-8 * (1 + z)
    expected = np.pi * r_e * 10**logN * atom["f"] * lam_cm**2 / (1 + z) * 1e8
    assert vm.equivalent_width_ang(line, z, logN, 10.0) == pytest.approx(expected, rel=0.01)


def test_detection_limit():
    truth, _, _ = vm.make_truth(30.0, 45000.0, 3.5, 2)
    def system(logN):
        return dict(species="CIV", z=3.0,
                    components=[MockComponent("CIV", 3.0, logN, 10.0)])
    strong, ew_s, sigma = vm.is_detectable(system(13.5), truth)
    weak, ew_w, _ = vm.is_detectable(system(11.0), truth)
    assert strong and not weak
    assert ew_s > ew_w > 0 and sigma > 0


# --- Metrics -------------------------------------------------------------

class _Spec:
    def __init__(self, candidates):
        self.meta = {"series_map_json": {str(i): s for i, (s, _) in enumerate(candidates)},
                     "z_map_json": {str(i): z for i, (_, z) in enumerate(candidates)}}


def _systems():
    return [dict(species="CIV", z=2.0, components=[MockComponent("CIV", 2.0, 14.0, 10.0)]),
            dict(species="CIV", z=2.3, components=[MockComponent("CIV", 2.3, 14.0, 10.0)]),
            dict(species="SiIV", z=2.1, components=[MockComponent("SiIV", 2.1, 13.5, 8.0)])]


def test_identification_metrics():
    truth, _, _ = vm.make_truth(100.0, 45000.0, 2.5, 2)
    cands = [("CIV", _z_offset(2.0, 5.0)), ("CIV", _z_offset(2.0, 6.0)),  # duplicate
             ("CIV", _z_offset(2.15, 0.0)),  # false
             ("CIV", 2.6),  # false, beyond z_em
             ("SiIV", _z_offset(2.1, -15.0)), ("Ly_a", 2.4)]
    m = vm.identification_metrics(truth, _systems(), _Spec(cands))
    assert m["n_lya"] == 1
    civ = m["CIV"]
    assert civ["n_true"] == 2 and civ["n_candidates"] == 4
    assert civ["n_found"] == 1 and civ["n_duplicates"] == 1
    assert civ["n_false"] == 2 and civ["n_beyond_z_em"] == 1
    assert m["SiIV"]["n_found"] == 1 and m["SiIV"]["n_false"] == 0
    assert m["MgII"]["n_candidates"] == 0


def _comp(series, z, logN, b, dz=1e-5, dlogN=0.05, db=0.5, chi2=1.0):
    return ComponentDataV2(id=1, z=z, dz=dz, logN=logN, dlogN=dlogN, b=b,
                           db=db, series=series, chi2=chi2)


class _Systs:
    def __init__(self, comps):
        self.components = comps


def test_fit_metrics():
    fitted = [_comp("CIV", _z_offset(2.0, 3.0), 14.1, 11.0),
              _comp("CIV", _z_offset(2.0, 1.0), 12.0, 1.0, dlogN=2.0),
              _comp("SiIV", _z_offset(2.1, 50.0), 13.5, 8.0),  # too far
              _comp("Ly_a", 2.4, 14.0, 20.0)]  # not a metal: ignored
    m = vm.fit_metrics(None, _systems(), _Systs(fitted))
    assert m["n_true"] == 3 and m["n_fitted"] == 3
    assert m["n_matched"] == 1 and m["n_spurious"] == 2 and m["n_missed"] == 2
    # The closest pair wins: the true CIV goes with the weak component 1 km/s
    # away, and the good one 3 km/s away counts as spurious (a limit of the
    # metric when a fit splits a line)
    assert m["diffs"]["logN"] == pytest.approx([-2.0])
    assert m["flags"]["absurd_errors"] == 1 and m["flags"]["b_at_bound"] == 1
    assert m["n_flagged"] == 1


def test_fit_metrics_pulls():
    true = [dict(species="CIV", z=2.0, components=[MockComponent("CIV", 2.0, 14.0, 10.0)])]
    dz = 2e-5
    fitted = [_comp("CIV", 2.0 + 2 * dz, 14.1, 11.0, dz=dz, dlogN=0.05, db=0.5)]
    m = vm.fit_metrics(None, true, _Systs(fitted))
    assert m["pulls"]["z"] == pytest.approx([2.0])
    assert m["pulls"]["logN"] == pytest.approx([2.0])
    assert m["pulls"]["b"] == pytest.approx([2.0])
    assert m["diffs"]["z"] == pytest.approx([2 * dz / 3.0 * C_KMS])


def test_continuum_metrics_of_the_true_continuum():
    truth, _, _ = vm.make_truth(30.0, 8900.0, 2.5, 2)
    spec = generate(truth)  # with the true continuum in "cont"
    m = vm.continuum_metrics(truth, spec)
    for region in ("forest", "red"):
        assert m[region]["median"] == pytest.approx(1.0)
        assert m[region]["p05"] == pytest.approx(1.0)


def test_schema_params_are_the_gui_defaults():
    params = vm.schema_params("absorbers", "doublets_auto")
    assert params == {"multiplets": "CIV,SiIV,MgII", "score_threshold": "0.5",
                      "z_window_kms": "20.0", "max_b": "50.0"}
    assert "debug_rating" not in vm.schema_params("absorbers", "identify_lines")
    assert vm.schema_params("continuum", "estimate_auto")["fudge"] == "auto"


# --- A run --------------------------------------------------------------

def test_tiny_grid():
    result = vm.run("tiny")
    (record,) = result["records"]
    assert record["n_pixels"] > 1000
    assert set(record["times"]) == {"mock", "continuum", "identify", "doublets"}
    s = result["summary"]["all"]
    assert s["n_mocks"] == 1
    assert 0.5 < s["continuum_red"]["median_of_medians"] < 1.5
    for species in vm.METALS:
        assert 0 <= s["identification"][species]["n_found"] <= vm.METAL_SYSTEMS[species]
    assert s["fit"]["n_true"] > 0
    assert set(result["summary"]["by_snr"]) == {"30.0"}


def test_robust_stats():
    values = list(np.random.default_rng(0).standard_normal(10000)) + [1e9]
    st = vm._stats(values)
    assert st["sigma"] == pytest.approx(1.0, abs=0.03)
    assert st["median"] == pytest.approx(0.0, abs=0.03)
    assert st["frac_gt3"] == pytest.approx(0.0027, abs=0.002)
    assert vm._stats([np.nan]) is None


def test_pull_categories():
    true = MockComponent("CIV", 2.0, 14.0, 10.0)
    alone = [_comp("CIV", _z_offset(2.0, 2.0), 14.0, 10.0),
             _comp("CIV", _z_offset(2.0, 80.0), 13.0, 10.0),  # beyond CROWD_KMS
             _comp("SiIV", 2.0, 13.0, 10.0)]  # another series
    crowded = alone + [_comp("CIV", _z_offset(2.0, -10.0), 12.5, 5.0)]
    assert vm.pull_category(true, 1, alone) == "single_clean"
    assert vm.pull_category(true, 1, crowded) == "single_crowded"
    assert vm.pull_category(true, 3, alone) == "multi_clean"
    assert vm.pull_category(true, 3, crowded) == "multi_crowded"


def test_fit_metrics_pulls_by_category():
    systems = [dict(species="CIV", z=2.0, components=[MockComponent("CIV", 2.0, 14.0, 10.0)])]
    fitted = [_comp("CIV", 2.0, 14.1, 10.0, dlogN=0.1)]
    m = vm.fit_metrics(None, systems, _Systs(fitted))
    assert m["pulls_by_category"]["single_clean"]["logN"] == pytest.approx([1.0])
    assert m["pulls_by_category"]["multi_crowded"]["logN"] == []


def test_completeness_by_significance():
    records = [{"identification": {"systems": [
        dict(species="CIV", significance=0.5, found=False),
        dict(species="CIV", significance=2.5, found=False),
        dict(species="CIV", significance=2.7, found=True),
        dict(species="SiIV", significance=20.0, found=True)]}}]
    c = vm.completeness_by_significance(records)
    assert c["0-1"] == dict(n=1, found=0.0)
    assert c["2-3"] == dict(n=2, found=0.5)
    assert c["10-inf"] == dict(n=1, found=1.0)
    assert c["3-5"] == dict(n=0, found=None)


def test_system_column():
    system = dict(species="CIV", z=2.0, components=[
        MockComponent("CIV", 2.0, 13.0, 10.0),
        MockComponent("CIV", _z_offset(2.0, 30.0), 13.0, 10.0)])
    # One fitted component carrying the whole column, within the range
    fitted = [_comp("CIV", _z_offset(2.0, 15.0), np.log10(2e13), 15.0, dlogN=0.05),
              _comp("CIV", _z_offset(2.0, 200.0), 12.0, 10.0),  # outside
              _comp("SiIV", 2.0, 13.0, 10.0)]  # another species
    c = vm.system_column(system, fitted)
    assert c["n_true"] == 2 and c["n_fitted"] == 1
    assert c["fit_logN"] == pytest.approx(c["true_logN"])
    assert c["dlogN"] == pytest.approx(0.05)
    assert vm.system_column(system, fitted[1:])["fit_logN"] is None


# --- Regression guard on the reduced grid (slow) ----------------------------

@pytest.mark.slow
def test_reduced_grid_does_not_regress():
    """The chain on the reduced grid (16 mocks, ~2 min with 4 workers)
    against loose thresholds set from the baseline of phase 4 (version
    2.1.0; ROADMAP.md, "Baseline"): it catches regressions, not the
    improvements the phase aims at. Values at the baseline in comments,
    and after an intended change of the phase where it moved them."""
    import os
    result = vm.run("reduced", workers=min(4, os.cpu_count() or 1))
    s = result["summary"]["all"]
    assert 0.99 < s["continuum_red"]["median_of_medians"] < 1.01      # 1.003
    assert s["continuum_forest"]["median_of_medians"] > 0.85           # 0.929
    for species in vm.METALS:
        t = s["identification"][species]
        assert t["purity"] > 0.9, species                              # 0.98-1.00
        assert t["completeness"] > 0.6, species                        # 0.72-0.82
        assert t["n_beyond_z_em"] == 0, species
    f = s["fit"]
    # One candidate per doublet (F49): far fewer spurious components, but
    # the systems of several close components are fitted with fewer (F58),
    # so fewer true components are matched and the pulls are biased
    assert f["n_matched"] / f["n_true"] > 0.25                         # 0.46; 0.34 (F49)
    assert f["n_spurious"] / f["n_fitted"] < 0.1                       # 0.15; 0.02 (F49)
    # b: 1.26 after F61. The damping wings, (1 + z) too strong, broadened
    # the profiles and offset the positive bias of the systems fitted with
    # fewer components (F49, F58); single components are less biased now
    # (full grid: -0.25 -> -0.13)
    for par in ("z", "logN", "b"):
        assert abs(f[f"pull_{par}"]["median"]) < 1.5, par              # |0.14|; b 0.84 (F49); b 1.26 (F61)
    column = f["system_column"]
    assert abs(column["diff"]["median"]) < 0.05                        # 0.000; -0.022 (F49)
    assert column["diff"]["sigma"] < 0.4                               # 0.236; 0.175 (F49)
    assert column["pull"]["sigma"] < 3.0                               # 1.21; 2.00 (F49)
