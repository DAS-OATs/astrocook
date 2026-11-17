"""The line methods of Spectrum: detection of the absorbed regions,
identification, velocity bounds, redshift of a doublet across two spectra.

Characterization tests written before the methods were moved to their own
module (phase 5, split of Spectrum): the paths that the identification
tests (``test_identify_lines.py``, ``test_snapshots_identification.py``)
do not reach. ``detect_absorptions``, ``get_velocity_bounds`` and
``detect_doublet_z`` had no test.
"""
import json

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.spectrum_operations import detect_regions
from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.mocks.physics import C_KMS

Z_EM = 3.0


def _mock(components, lam_min_nm=490.0, lam_max_nm=700.0, seed=1):
    return generate(MockTruth(components=components, lam_min_nm=lam_min_nm,
                              lam_max_nm=lam_max_nm, dv_kms=2.5, resol=45000.0,
                              snr=100.0, z_em=Z_EM, seed=seed))


@pytest.fixture(scope="module")
def civ():
    """A CIV doublet on the red side, with its mask of absorbed pixels."""
    spec = _mock([MockComponent("CIV", 2.5, 14.0, 10.0)])
    return spec.find_absorbed(5000.0 * au.km / au.s, 400.0 * au.km / au.s, 2.0, False)


def _systems(spec):
    series = json.loads(spec.meta["series_map_json"])
    z = json.loads(spec.meta["z_map_json"])
    return sorted((series[k], z[k]) for k in series)


# --- detect_absorptions ---------------------------------------------------------

def test_detect_absorptions_labels_the_regions(civ):
    new = civ.detect_absorptions(min_pix=3)
    mask = civ.get_column("abs_mask").value.astype(bool)
    expected, n = detect_regions(mask, min_pix=3)
    ids = new.get_column("region_id")
    np.testing.assert_array_equal(ids.value, expected)
    assert ids.unit == au.dimensionless_unscaled and n > 0
    assert new.history[-1] == f"Detected {n} absorbers from abs_mask"


def test_detect_absorptions_casts_an_integer_mask(civ):
    mask = civ.get_column("abs_mask").value.astype(bool)
    as_int = civ.update_column("abs_mask", mask.astype(int), unit=au.dimensionless_unscaled)
    np.testing.assert_array_equal(as_int.detect_absorptions().get_column("region_id").value,
                                  civ.detect_absorptions().get_column("region_id").value)


def test_detect_absorptions_refuses_a_float_mask(civ):
    as_float = civ.update_column("abs_mask", civ.get_column("abs_mask").value * 1.0,
                                 unit=au.dimensionless_unscaled)
    with pytest.raises(TypeError, match="expected boolean or integer mask"):
        as_float.detect_absorptions()


def test_detect_absorptions_without_the_mask(civ):
    with pytest.raises(ValueError, match="Run 'find_absorbed' first"):
        civ.detect_absorptions(mask_col="missing")


# --- identify_lines -----------------------------------------------------------------

def test_identify_lines_needs_z_em(civ):
    with pytest.raises(ValueError, match="z_em must be set"):
        civ.with_properties(z_em=0.0).identify_lines(["CIV"])


def test_identify_lines_makes_the_mask_if_missing(civ):
    raw = _mock([MockComponent("CIV", 2.5, 14.0, 10.0)])
    assert not raw.has_aux_column("abs_mask")
    found = raw.identify_lines(["CIV"])
    expected = civ.identify_lines(["CIV"])  # the mask of find_absorbed's defaults
    assert _systems(found) == _systems(expected)
    np.testing.assert_array_equal(found.get_column("abs_ids").value,
                                  expected.get_column("abs_ids").value)
    assert not found.has_aux_column("abs_mask")  # made for the search only


def test_identify_lines_writes_its_results(civ):
    new = civ.identify_lines(["CIV"])
    systems = _systems(new)
    assert [s for s, _ in systems] == ["CIV"]
    assert abs(systems[0][1] - 2.5) / 3.5 * C_KMS < 5.0
    for key in ("series_map_json", "z_map_json", "logN_map_json", "b_map_json",
                "region_identifications"):
        json.loads(new.meta[key])
    assert new.history[-1] == "Identified 1 systems (de-blended)"
    assert civ.meta.get("series_map_json") is None  # the original is not changed


@pytest.mark.parametrize("bypass_scoring", [False, True])
def test_no_metal_doublet_in_the_forest(bypass_scoring):
    # A CIV doublet in the Lyman-alpha forest (Ly-a at 486 nm): the metal
    # doublets are searched on the red side only, so its regions are taken
    # as Ly-a, even without the scoring
    spec = _mock([MockComponent("CIV", 2.0, 14.0, 10.0)], lam_min_nm=440.0,
                 lam_max_nm=700.0)
    systems = _systems(spec.identify_lines(["CIV"], bypass_scoring=bypass_scoring))
    assert systems and all(s == "Ly_a" for s, _ in systems)


def _candidates(spec):
    meta = {k: json.loads(spec.meta[f"{k}_map_json"]) for k in ("series", "z", "score")}
    return sorted((meta["series"][k], meta["z"][k], meta["score"][k]) for k in meta["series"])


def test_ly_ab_pairs_are_rated():
    # A Ly-a/Ly-b pair in the forest. Until phase 7, PR 7c.1 the score was
    # always 0 and the Ly-b line became a spurious Ly-a at z = 2.29
    # (findings register F68): rate_doublet_candidate measured velocities
    # in a linear frame (c dlambda / lambda) and the lags of the lines in a
    # logarithmic one, ~4100 km/s apart for Ly-a/Ly-b. Now the pair is
    # recognised as with bypass_scoring, and gives a Ly-a system
    spec = _mock([MockComponent("Ly_ab", 2.9, 14.0, 25.0)], lam_min_nm=390.0,
                 lam_max_nm=520.0)
    z_lyb_as_lya = 3.9 * ATOM_DATA["Ly_b"]["wave"] / ATOM_DATA["Ly_a"]["wave"] - 1.0
    scored = spec.identify_lines(["Ly_ab"])
    bypassed = spec.identify_lines(["Ly_ab"], bypass_scoring=True)
    assert _systems(scored) == _systems(bypassed)
    pair = [c for c in _candidates(scored) if abs(c[1] - 2.9) < 1e-4]
    assert len(pair) == 1 and pair[0][0] == "Ly_a" and pair[0][2] > 0.9
    assert not any(abs(z - z_lyb_as_lya) < 1e-3 for _, z in _systems(scored))
    assert not any(s == "Ly_b" for s, _ in _systems(scored))


@pytest.mark.parametrize(("series", "z", "logN", "b", "lam_min_nm"), [
    ("CIV", 2.5, 14.0, 10.0, 490.0),   # 0.996 before PR 7c.1
    ("MgII", 1.2, 13.0, 5.0, 490.0),   # 0.960
    ("SiIV", 2.6, 13.0, 5.0, 490.0),   # 0.205: dropped
    ("SiIV", 2.6, 13.0, 3.0, 490.0),   # 0: dropped
    ("SiIV", 2.6, 12.9, 9.0, 490.0),   # 0.559 (the system of F119)
])
def test_clean_doublets_score_high(series, z, logN, b, lam_min_nm):
    # The lines of a clean doublet aligned in the frame of their lags (F68):
    # the linear frame put them 0.4 (CIV), 1 (MgII) and 6 km/s (SiIV)
    # apart, enough to drop a narrow SiIV (the scores before in the
    # comments, S/N 100, R 45000)
    from astrocook.core import spectrum_lines

    spec = _mock([MockComponent(series, z, logN, b)], lam_min_nm=lam_min_nm)
    spec = spec.find_absorbed(5000.0 * au.km / au.s, 400.0 * au.km / au.s, 2.0, False)
    scores = []
    original = spectrum_lines.rate_doublet_candidate

    def record(*args, **kwargs):
        scores.append(original(*args, **kwargs))
        return scores[-1]

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(spectrum_lines, "rate_doublet_candidate", record)
        systems = _systems(spec.identify_lines([series]))
    assert len(scores) == 1 and scores[0] > 0.95
    assert len(systems) == 1 and systems[0][0] == series
    assert abs(systems[0][1] - z) / (1 + z) * C_KMS < 5.0


def test_doublet_score_frame_is_logarithmic():
    # An exact doublet, F_1 = F_2 ** gamma with the lines at c ln(l_2/l_1)
    # of each other, built on the pixels without the mocks: R^2 = 1 to the
    # linear
    # interpolation. In the linear frame of the scoring before PR 7c.1 the
    # Ly-a/Ly-b pair scored 0 (F68)
    from astrocook.core.atomic_data import STANDARD_MULTIPLETS
    from astrocook.core.spectrum import Spectrum
    from astrocook.core.spectrum_operations import rate_doublet_candidate
    from astrocook.core.structures import DataColumn, SpectrumData

    z = 2.9
    line_1, line_2 = STANDARD_MULTIPLETS["Ly_ab"]
    lam = {t: ATOM_DATA[t]["wave"] / 10.0 * (1 + z) for t in (line_1, line_2)}
    x = np.exp(np.arange(np.log(390.0), np.log(480.0), 2.5 / C_KMS))
    gamma = (ATOM_DATA[line_1]["f"] * ATOM_DATA[line_1]["wave"]
             / (ATOM_DATA[line_2]["f"] * ATOM_DATA[line_2]["wave"]))
    tau_2 = 0.8 * np.exp(-0.5 * (C_KMS * np.log(x / lam[line_2]) / 20.0) ** 2)
    tau_1 = gamma * 0.8 * np.exp(-0.5 * (C_KMS * np.log(x / lam[line_1]) / 20.0) ** 2)
    flux = np.exp(-tau_1 - tau_2)
    one = np.ones_like(x)
    spec = Spectrum(SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x, au.nm), xmax=DataColumn(x, au.nm),
        y=DataColumn(flux, au.dimensionless_unscaled),
        dy=DataColumn(0.01 * one, au.dimensionless_unscaled),
        aux_cols={"cont": DataColumn(one, au.dimensionless_unscaled)}, z_em=3.0))
    # The regions as the kinematic search gives them: the first is that of
    # the second line of the multiplet (Ly-b), the second that of the first
    near = lambda t: np.abs(C_KMS * np.log(x / lam[t])) < 60.0  # noqa: E731
    score = rate_doublet_candidate(spec, near(line_2), near(line_1), "Ly_ab", z)
    assert score == pytest.approx(1.0, abs=1e-4)  # the linear interpolation


def test_identify_lines_casts_an_integer_mask(civ):
    mask = civ.get_column("abs_mask").value.astype(bool)
    as_int = civ.update_column("abs_mask", mask.astype(int), unit=au.dimensionless_unscaled)
    assert _systems(as_int.identify_lines(["CIV"])) == _systems(civ.identify_lines(["CIV"]))


# --- get_velocity_bounds, detect_doublet_z -----------------------------------------

def test_velocity_bounds(civ):
    x = civ.x.value
    lam = ATOM_DATA["CIV_1548"]["wave"] / 10.0 * 3.5
    v_min, v_max = civ.get_velocity_bounds(2.5, "CIV_1548")
    assert v_min == pytest.approx(C_KMS * (x.min() - lam) / lam, rel=1e-12)
    assert v_max == pytest.approx(C_KMS * (x.max() - lam) / lam, rel=1e-12)
    assert civ.get_velocity_bounds(2.5, "not_a_line") == (-np.inf, np.inf)


def test_doublet_redshift_across_two_spectra():
    z = 2.5
    # The two lines of a CIV doublet in two spectra that do not overlap
    blue = _mock([MockComponent("CIV", z, 14.0, 10.0)], 540.0, 542.5, seed=2)
    red = _mock([MockComponent("CIV", z, 14.0, 10.0)], 542.5, 545.0, seed=3)
    blue = blue.update_column("y", blue.y.value / blue.get_column("cont").value)
    red = red.update_column("y", red.y.value / red.get_column("cont").value)
    found = blue.detect_doublet_z(red, "CIV_1548", "CIV_1550")
    assert abs(found - z) / (1 + z) * C_KMS < 5.0


def test_doublet_redshift_without_overlap():
    blue = _mock([], 540.0, 541.0)
    red = _mock([], 600.0, 601.0)
    # A ValueError since phase 6-quater, PR 10 (F145: 0.0, as a redshift)
    with pytest.raises(ValueError, match="cover no redshift in common"):
        blue.detect_doublet_z(red, "CIV_1548", "CIV_1550")
