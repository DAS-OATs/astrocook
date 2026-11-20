"""The identification as a comparison between hypotheses (phase 7, PR 7c.4c):
``identify_lines(method="hypotheses")``, :mod:`astrocook.fitting.hypotheses`."""
import json

import numpy as np
import pytest

from astrocook.core import identification as ident
from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.spectrum_operations import rate_doublet_candidate
from astrocook.fitting import hypotheses as H
from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.mocks.physics import C_KMS
from astrocook.recipes.errors import InvalidParameter
from astrocook.session import Session

Z_EM = 3.0

#: Metal systems on the red side (Ly-a at 486 nm), one of two components
COMPONENTS = [
    MockComponent("CIV", 2.50, 14.0, 10.0),
    MockComponent("CIV", 2.70, 13.8, 12.0),
    MockComponent("CIV", 2.70 + 40.0 / C_KMS * 3.70, 13.6, 10.0),
    MockComponent("SiIV", 2.65, 13.5, 10.0),
    MockComponent("MgII", 1.20, 13.5, 8.0),
]


def _mock(components, lam=(490.0, 700.0), snr=100.0, resol=45000.0, seed=1):
    return generate(MockTruth(components=components, lam_min_nm=lam[0], lam_max_nm=lam[1],
                              dv_kms=2.5, resol=resol, snr=snr, z_em=Z_EM, seed=seed))


@pytest.fixture(scope="module")
def result():
    return H.identify_hypotheses(_mock(COMPONENTS), ["CIV", "SiIV", "MgII"], screen=None)


def _systems(spec):
    series = json.loads(spec.meta["series_map_json"])
    z = json.loads(spec.meta["z_map_json"])
    return [(series[k], z[k]) for k in series]


def _dv(z1, z2):
    return abs(z1 - z2) / (1 + z2) * C_KMS


@pytest.mark.parametrize("comp", COMPONENTS, ids=lambda c: f"{c.series}-{c.z:.4f}")
def test_every_true_component_has_a_candidate(result, comp):
    found = [z for s, z in _systems(result.spectrum) if s == comp.series]
    assert min(_dv(z, comp.z) for z in found) < 5.0


def test_no_false_candidates(result):
    for series, z in _systems(result.spectrum):
        assert any(c.series == series and _dv(z, c.z) < 20.0 for c in COMPONENTS)


def test_the_structure(result):
    st = result.spectrum.meta[ident.META_KEY]
    assert st.method == "hypotheses"
    assert st.params["criterion"] == "spic"
    compared = [g for g in st.groups if len(g.hypotheses) > 1]
    assert len(compared) == 4  # CIV twice, SiIV, MgII
    for g in compared:
        weights = [h.weight for h in g.hypotheses if h.status != ident.Status.DROPPED]
        assert sum(weights) == pytest.approx(1.0)
        accepted = [h for h in g.hypotheses if h.status == ident.Status.ACCEPTED]
        assert len(accepted) == 1
        assert accepted[0].weight == max(weights)
        assert accepted[0].stage == ident.Stage.FIT
        assert accepted[0].anchors  # a pair, not the default labels
        assert {"spic", "bic", "aic", "aicc"} <= set(accepted[0].criteria)
        for h in g.hypotheses:
            assert all(0 <= a < len(g.edges) for a in h.anchors)
        # The hypothesis of the default labels (unidentified on the red side)
        assert any(not h.anchors for h in g.hypotheses)
    # Every region is in abs_ids, the unidentified ones too
    ids = set(np.unique(result.spectrum.get_column("abs_ids").value).astype(int)) - {0}
    assert ids == {r.id for r in st.regions}
    assert any(r.kind == ident.RegionKind.UNIDENTIFIED for r in st.regions)
    metal = {t for r in st.regions if r.kind == ident.RegionKind.METAL for t in r.transitions}
    assert {"CIV_1548", "CIV_1550", "SiIV_1393", "SiIV_1402", "MgII_2796",
            "MgII_2803"} <= metal


def test_the_candidates_carry_the_weight_and_the_fitted_parameters(result):
    meta = result.spectrum.meta
    scores = json.loads(meta["score_map_json"])
    logn = json.loads(meta["logN_map_json"])
    assert all(0.5 < s <= 1.0 for s in scores.values())
    series = json.loads(meta["series_map_json"])
    for k, s in series.items():
        true = [c for c in COMPONENTS if c.series == s]
        assert any(abs(logn[k] - c.logN) < 0.3 for c in true)


def test_json_round_trip(result):
    st = result.spectrum.meta[ident.META_KEY]
    back = ident.Identification.from_dict(json.loads(json.dumps(st.to_dict())))
    assert back == st


def test_the_ly_a_only_when_asked():
    # F51: a line of the forest is Ly-a only when Ly-a is asked (assumed:
    # no Ly-b links it); otherwise it is unidentified
    comps = [MockComponent("Ly_a", 2.90, 13.5, 25.0), MockComponent("CIV", 2.50, 14.0, 10.0)]
    spec = _mock(comps, lam=(450.0, 560.0))
    without = H.identify_hypotheses(spec, ["CIV"]).spectrum
    assert {s for s, _ in _systems(without)} == {"CIV"}
    with_lya = H.identify_hypotheses(spec, ["CIV", "Ly_a"]).spectrum
    lya = [z for s, z in _systems(with_lya) if s == "Ly_a"]
    assert min(_dv(z, 2.90) for z in lya) < 5.0
    st = with_lya.meta[ident.META_KEY]
    kinds = {r.kind for r in st.regions if r.transitions == ("Ly_a",)}
    assert kinds == {ident.RegionKind.FOREST_ASSUMED}
    scores = json.loads(with_lya.meta["score_map_json"])
    series = json.loads(with_lya.meta["series_map_json"])
    assert all(scores[k] == 0.0 for k, s in series.items() if s == "Ly_a")


def test_a_ly_a_with_its_ly_b_is_a_ly_ab_system():
    # Both the Ly-a (474 nm) and the Ly-b (400 nm) in the forest of z_em 3
    comps = [MockComponent("Ly_ab", 2.90, 14.5, 25.0)]
    spec = _mock(comps, lam=(380.0, 480.0), snr=50.0)
    out = H.identify_hypotheses(spec, ["Ly_ab"]).spectrum
    lyab = [z for s, z in _systems(out) if s == "Ly_ab"]
    assert len(lyab) >= 1 and min(_dv(z, 2.90) for z in lyab) < 5.0
    st = out.meta[ident.META_KEY]
    labelled = {t for r in st.regions if r.kind == ident.RegionKind.FOREST for t in r.transitions}
    assert labelled == {"Ly_a", "Ly_b"}


def test_two_unrelated_lines_at_the_lag_are_not_a_doublet():
    # The "1550" twice as strong as the "1548": no doublet ratio
    z = 2.50
    lam2 = ATOM_DATA["CIV_1550"]["wave"] * (1 + z)
    comps = [MockComponent("Ly_a", ATOM_DATA["CIV_1548"]["wave"] * (1 + z)
                           / ATOM_DATA["Ly_a"]["wave"] - 1, 13.0, 20.0),
             MockComponent("Ly_a", lam2 / ATOM_DATA["Ly_a"]["wave"] - 1, 13.8, 20.0)]
    # (prefilter 0: the R^2 of the pair is 0, the pre-filter would drop it)
    out = H.identify_hypotheses(_mock(comps), ["CIV"], prefilter=0.0)
    st = out.spectrum.meta[ident.META_KEY]
    group = [g for g in st.groups if any(_dv(e.z, z) < 20 for e in g.edges)]
    assert len(group) == 1  # the pair is compared...
    accepted = [h for h in group[0].hypotheses if h.status == ident.Status.ACCEPTED]
    assert not accepted[0].anchors  # ...and left unidentified
    assert not any(_dv(zc, z) < 20 for s, zc in _systems(out.spectrum) if s == "CIV")


def test_the_pair_score_is_that_of_the_current_method():
    spec = _mock([MockComponent("CIV", 2.50, 13.8, 12.0)], lam=(530.0, 560.0))
    x = spec.x.to_value("nm")
    flux = spec.y.value / spec.cont.value
    v = C_KMS * np.log(x)
    z = 2.50
    m1 = np.abs(C_KMS * np.log(x / (ATOM_DATA["CIV_1548"]["wave"] / 10 * (1 + z)))) < 40
    m2 = np.abs(C_KMS * np.log(x / (ATOM_DATA["CIV_1550"]["wave"] / 10 * (1 + z)))) < 40
    i1, i2 = np.flatnonzero(m1), np.flatnonzero(m2)
    ours = H.pair_score(v, flux, (i1[0], i1[-1]), (i2[0], i2[-1]),
                        ("CIV_1548", "CIV_1550"), z)
    # rate_doublet_candidate: the region of the redder line first
    theirs = rate_doublet_candidate(spec, m2, m1, "CIV", z)
    assert ours == pytest.approx(theirs, abs=1e-3)
    assert ours > 0.9


def _region(rid):
    return H._Region(id=rid, lo=10 * rid, hi=10 * rid + 5, zone="red", significance=5.0,
                     peaks=(H._Peak(10 * rid + 2, 20.0),))


def test_a_group_too_large_is_split_at_its_weakest_edges():
    regions = [_region(i) for i in range(1, 9)]
    scores = [0.9, 0.8, 0.3, 0.7, 0.6, 0.5, 0.4]
    edges = [H._Candidate(i, "CIV", ("CIV_1548", "CIV_1550"), (i + 1, i + 2), 2.5, s)
             for i, s in enumerate(scores)]
    groups, cut = H._groups(regions, edges, max_group=4)
    # 8 regions: the weakest edge (0.3) cut, then the weakest of the group
    # of 5 left (0.4)
    assert [e.score for e in cut] == [0.3, 0.4]
    assert [g[0] for g in groups] == [[1, 2, 3], [4, 5, 6, 7], [8]]
    assert [len(g[1]) for g in groups] == [2, 3, 0]


def test_a_group_too_large_for_the_fit_keeps_its_default_labels(monkeypatch):
    # MAX_FIT_K: the time of the short fit; a group with a hypothesis above
    # it is not compared, its default labels accepted with the AOD
    monkeypatch.setattr(H, "MAX_FIT_K", 2)
    res = H.identify_hypotheses(_mock(COMPONENTS[:1]), ["CIV"], screen=None)
    st = res.spectrum.meta[ident.META_KEY]
    groups = [g for g in st.groups if len(g.hypotheses) > 1]
    assert groups
    for g in groups:
        accepted = [h for h in g.hypotheses if h.status == ident.Status.ACCEPTED]
        assert len(accepted) == 1 and not accepted[0].anchors
        assert "not compared" in accepted[0].reason
        assert all(h.status == ident.Status.DROPPED for h in g.hypotheses if h is not accepted[0])
    assert not [s for s, _ in _systems(res.spectrum) if s == "CIV"]


def test_the_species():
    assert H.parse_species(["all"])[0] == list(H.ALL_METALS)
    assert H.parse_species(["CIV", "Ly_a"]) == (["CIV"], True)
    assert H.parse_species(["CIV", "Ly_ab"]) == (["CIV", "Ly_ab"], True)
    assert H.parse_species(["CIV"]) == (["CIV"], False)
    with pytest.raises(ValueError, match="No anchor"):
        H.parse_species(["OI"])
    assert ("FeII", ("FeII_2382", "FeII_2600")) in H.anchor_pairs(["FeII"])


# --- The recipes -----------------------------------------------------------

@pytest.fixture(scope="module")
def session():
    return Session("mock", None, spec=_mock(COMPONENTS))


def test_the_recipe(session):
    out = session.absorbers.identify_lines(multiplets="CIV,SiIV,MgII", method="hypotheses")
    assert out.spec.meta[ident.META_KEY].method == "hypotheses"
    assert {s for s, _ in _systems(out.spec)} == {"CIV", "SiIV", "MgII"}


def test_the_recipe_refuses_an_unknown_species(session):
    with pytest.raises(InvalidParameter, match="No anchor"):
        session.absorbers.identify_lines(multiplets="CIV,XYZ", method="hypotheses")


def test_the_current_method_is_the_default(session):
    out = session.absorbers.identify_lines(multiplets="CIV,SiIV,MgII")
    assert out.spec.meta[ident.META_KEY].method == "doublets"


def test_the_pipeline_with_the_hypotheses(session):
    out = session.absorbers.doublets_auto(method="hypotheses")
    fitted = [(c.series, c.z, c.logN) for c in out.systs.components]
    assert len(fitted) == len(COMPONENTS)
    for comp in COMPONENTS:
        assert any(s == comp.series and _dv(z, comp.z) < 5.0 and abs(n - comp.logN) < 0.2
                   for s, z, n in fitted)


def test_populate_starts_from_the_short_fit(session):
    ident_session = session.absorbers.identify_lines(multiplets="CIV,SiIV,MgII",
                                                     method="hypotheses")
    meta = ident_session.spec.meta
    out = ident_session.absorbers.populate_from_identification(
        series_map_json=meta["series_map_json"], z_map_json=meta["z_map_json"])
    logn = json.loads(meta["logN_map_json"])
    assert sorted(round(c.logN, 6) for c in out.systs.components) == sorted(
        round(v, 6) for v in logn.values())


# --- PR 7c.4c-bis ------------------------------------------------------------

def _arrays(spec):
    flux = spec.y.value / spec.cont.value
    err = spec.dy.value / spec.cont.value
    return C_KMS * np.log(spec.x.to_value("nm")), flux, err


def _at(v, line, z, dv=0.0):
    return C_KMS * np.log(ATOM_DATA[line]["wave"] / 10 * (1 + z)) + dv


def test_the_score_of_a_pair_is_on_its_peaks():
    # F169: a CIV doublet whose partner region also holds an unrelated line
    # scores low on the regions, high on the peaks it matches
    z = 2.50
    lam = ATOM_DATA["CIV_1550"]["wave"] * (1 + z) * (1 + 90.0 / C_KMS)
    comps = [MockComponent("CIV", z, 13.8, 12.0),
             MockComponent("Ly_a", lam / ATOM_DATA["Ly_a"]["wave"] - 1, 14.0, 20.0)]
    spec = _mock(comps, lam=(530.0, 560.0))
    v, flux, err = _arrays(spec)
    fwhm = np.full(len(v), C_KMS / 45000.0)

    def region(rid, lo_kms, hi_kms, line):
        lo, hi = np.searchsorted(v, [_at(v, line, z, lo_kms), _at(v, line, z, hi_kms)])
        return H._Region(rid, int(lo), int(hi), "red", 10.0, H._peaks(v, flux, err, lo, hi))

    a, b = region(1, -40, 40, "CIV_1548"), region(2, -40, 160, "CIV_1550")
    assert len(b.peaks) == 2
    c = H._Candidate(0, "CIV", ("CIV_1548", "CIV_1550"), (1, 2), z, 0.0)
    by_id = {1: a, 2: b}
    assert H.anchor_peaks(c, by_id, v, fwhm) == [(0, 0)]
    on_regions = H.pair_score(v, flux, (a.lo, a.hi), (b.lo, b.hi), c.lines, z)
    on_peaks = H.peak_score(c, by_id, v, flux, fwhm)
    assert on_peaks > 0.9
    assert on_regions < on_peaks - 0.2


def test_a_region_is_split_at_its_gaps():
    # Two lines 300 km/s apart, the flux back to the continuum between them:
    # two regions of the identification, one each
    comps = [MockComponent("CIV", 2.50, 14.0, 10.0),
             MockComponent("CIV", 2.50 + 300.0 / C_KMS * 3.5, 14.0, 10.0)]
    spec = _mock(comps, lam=(530.0, 560.0))
    v, flux, err = _arrays(spec)
    lo, hi = np.searchsorted(v, [_at(v, "CIV_1548", 2.50, -60),
                                 _at(v, "CIV_1548", 2.50, 360)])
    labels = np.zeros(len(v), dtype=int)
    labels[lo:hi + 1] = 1
    sig = np.full(len(v), 10.0)
    regions, new = H._regions(spec, labels, sig, v, flux, err, Z_EM, None)
    assert len(regions) == 2
    assert regions[0].hi + 1 == regions[1].lo
    assert set(np.unique(new[lo:hi + 1])) == {1, 2}
    assert all(len(r.peaks) == 1 for r in regions)


def test_a_region_is_not_split_within_a_blend(monkeypatch):
    comps = [MockComponent("CIV", 2.50, 14.0, 10.0),
             MockComponent("CIV", 2.50 + 30.0 / C_KMS * 3.5, 14.0, 10.0)]
    spec = _mock(comps, lam=(530.0, 560.0))
    v, flux, err = _arrays(spec)
    lo, hi = np.searchsorted(v, [_at(v, "CIV_1548", 2.50, -60),
                                 _at(v, "CIV_1548", 2.50, 90)])
    labels = np.zeros(len(v), dtype=int)
    labels[lo:hi + 1] = 1
    regions, _ = H._regions(spec, labels, np.full(len(v), 10.0), v, flux, err, Z_EM, None)
    assert len(regions) == 1 and len(regions[0].peaks) == 2
    monkeypatch.setattr(H, "SPLIT_GAPS", False)
    spec2 = _mock(comps[:1] + [MockComponent("CIV", 2.50 + 300.0 / C_KMS * 3.5, 14.0, 10.0)],
                  lam=(530.0, 560.0))
    v, flux, err = _arrays(spec2)
    hi = np.searchsorted(v, _at(v, "CIV_1548", 2.50, 360))
    labels[lo:hi + 1] = 1
    regions, _ = H._regions(spec2, labels, np.full(len(v), 10.0), v, flux, err, Z_EM, None)
    assert len(regions) == 1


def test_the_window_of_a_group_leaves_out_the_other_regions():
    spec = _mock(COMPONENTS[:2])
    v, flux, err = _arrays(spec)
    systs = H._empty().add_component("CIV", 2.50, logN=14.0, b=10.0)
    lo, hi = np.searchsorted(v, [_at(v, "CIV_1548", 2.50, -40), _at(v, "CIV_1548", 2.50, 40)])
    labels = np.zeros(len(v), dtype=int)
    labels[lo:hi + 1] = 1
    # A region of another group inside the window of the first
    o_lo, o_hi = np.searchsorted(v, [_at(v, "CIV_1548", 2.50, 45), _at(v, "CIV_1548", 2.50, 52)])
    labels[o_lo:o_hi + 1] = 2
    ctx = H._Context(spec, None, "spic", None, False, False, labels, (0.0, 0.0))
    from astrocook.fitting.evaluation_window import EvaluationWindow
    full = EvaluationWindow.around(spec, [systs])
    assert full.mask[o_lo:o_hi + 1].all()
    win = H._window(ctx, [systs], [1])
    assert not win.mask[o_lo:o_hi + 1].any()
    assert win.mask[lo:hi + 1].all()
    assert (win.mask | (labels == 2)).sum() == full.mask.sum()


def _added(res):
    st = res.spectrum.meta[ident.META_KEY]
    return st, [len(h.components) for g in st.groups for h in g.hypotheses]


def test_the_refinement_is_an_option(monkeypatch):
    # F170: on by default; off, every hypothesis keeps the components of its
    # peaks; MAX_ADDED caps it
    spec = _mock(COMPONENTS[1:3])
    on = H.identify_hypotheses(spec, ["CIV"])
    off = H.identify_hypotheses(spec, ["CIV"], refine=False)
    assert on.spectrum.meta[ident.META_KEY].params["refine"] is True
    assert off.spectrum.meta[ident.META_KEY].params["refine"] is False
    n_on, n_off = _added(on)[1], _added(off)[1]
    assert sum(n_on) > sum(n_off)
    monkeypatch.setattr(H, "MAX_ADDED", 1)
    capped = H.identify_hypotheses(spec, ["CIV"])
    st = capped.spectrum.meta[ident.META_KEY]
    assert st.params["max_added"] == 1
    for g_off, g_cap in zip(off.spectrum.meta[ident.META_KEY].groups, st.groups, strict=True):
        for h_off, h_cap in zip(g_off.hypotheses, g_cap.hypotheses, strict=False):
            assert len(h_cap.components) <= len(h_off.components) + 1


def test_the_recipe_takes_refine(session):
    out = session.absorbers.identify_lines(multiplets="CIV", method="hypotheses", refine=False)
    assert out.spec.meta[ident.META_KEY].params["refine"] is False
