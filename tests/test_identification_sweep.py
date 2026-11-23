"""The sweep of the anchored systems (phase 7, PR 7c.4d):
:mod:`astrocook.fitting.hypotheses_sweep`, the second pass of
``identify_lines(method="hypotheses")``."""
import json

import pytest

from astrocook.core import identification as ident
from astrocook.core.atomic_data import series_lines
from astrocook.fitting import hypotheses as H
from astrocook.fitting import hypotheses_sweep as S
from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.mocks.physics import C_KMS
from astrocook.session import Session

Z_EM = 3.0
Z_SYS = 2.6
#: The forest from 430 nm (Ly-a emission at 486 nm), without its lines: the
#: H I and the low ions of the system are the only lines there
LAM = (430.0, 700.0)


def _at(dv):
    return Z_SYS + dv / C_KMS * (1 + Z_SYS)


def _rich(with_ions=True):
    """A CIV system of two components with its H I and, with ``with_ions``,
    SiII (forest and red side), SiIII (forest) and AlII (red side); a second
    CIV system alone."""
    comps = [MockComponent("Ly_a,Ly_b", Z_SYS, 15.5, 25.0),
             MockComponent("CIV", 2.8, 13.8, 10.0)]
    for dv, logn in ((0.0, 14.0), (25.0, 13.6)):
        comps.append(MockComponent("CIV", _at(dv), logn, 8.0))
        if with_ions:
            comps.append(MockComponent(S.SWEEP_SPECIES["SiII"], _at(dv), logn - 0.6, 6.0))
            comps.append(MockComponent("SiIII_1206", _at(dv), logn - 0.8, 7.0))
            comps.append(MockComponent("AlII_1670", _at(dv), logn - 1.0, 6.0))
    return generate(MockTruth(components=sorted(comps, key=lambda c: c.z), lam_min_nm=LAM[0],
                              lam_max_nm=LAM[1], dv_kms=2.5, resol=45000.0, snr=50.0,
                              z_em=Z_EM, seed=1))


@pytest.fixture(scope="module")
def rich():
    return H.identify_hypotheses(_rich(), ["CIV", "SiIV", "MgII", "Ly_a"]).spectrum


def _systems(spec):
    series = json.loads(spec.meta["series_map_json"])
    z = json.loads(spec.meta["z_map_json"])
    score = json.loads(spec.meta["score_map_json"])
    return [(series[k], z[k], score[k]) for k in series]


def _dv(z1, z2):
    return abs(z1 - z2) / (1 + z2) * C_KMS


def _sweep_groups(spec):
    return [g for g in spec.meta[ident.META_KEY].groups if g.anchored_on is not None]


@pytest.mark.parametrize("species", ["HI", "SiII", "SiIII", "AlII"])
def test_the_species_of_the_system_are_found(rich, species):
    series = S.SWEEP_SPECIES[species]
    found = [(z, w) for s, z, w in _systems(rich) if s == series]
    assert found and min(_dv(z, Z_SYS) for z, _ in found) < 10.0
    # Their score is the weight of the test against the labels
    assert all(0 < w <= 1 for _, w in found)


def test_no_species_at_the_lone_system(rich):
    # The CIV at z = 2.8 has no other ion: nothing is added to it
    added = [(s, z) for s, z, _ in _systems(rich) if s != "CIV" and _dv(z, 2.8) < 200]
    assert not added


def test_the_assumed_ly_a_becomes_the_h_i_of_the_system(rich):
    lya = [z for s, z, _ in _systems(rich) if s == "Ly_a" and _dv(z, Z_SYS) < 50]
    assert not lya
    st = rich.meta[ident.META_KEY]
    hi = next(g for g in _sweep_groups(rich) if g.hypotheses[0].label.startswith("HI"))
    assert hi.hypotheses[0].criteria["replaced_lya"] >= 1
    # The region of its Ly-a is labelled as identified forest, with Ly_a
    rid = hi.edges[0].regions[1]
    region = st.region(rid)
    assert region.kind == ident.RegionKind.FOREST and "Ly_a" in region.transitions


def test_the_structure_of_the_sweep(rich):
    st = rich.meta[ident.META_KEY]
    first = {g.id: g for g in st.groups if g.anchored_on is None}
    groups = _sweep_groups(rich)
    assert groups
    for g in groups:
        # Each test a group of one hypothesis, anchored on a group of the
        # first pass whose accepted hypothesis has a pair
        assert len(g.hypotheses) == 1
        anchor = first[g.anchored_on]
        assert any(h.anchors for h in anchor.hypotheses if h.status == ident.Status.ACCEPTED)
        h = g.hypotheses[0]
        assert "(sweep of CIV)" in h.label
        assert 0 <= h.weight <= 1
        if h.status == ident.Status.ACCEPTED:
            # An accepted species is an edge from the anchor's strong line
            assert h.anchors == (0,) and len(g.edges) == 1
            e = g.edges[0]
            assert e.transitions[0] == "CIV_1548"
            assert e.transitions[1] in series_lines(e.series)
            assert set(e.regions) <= set(g.regions)
            assert h.criteria["spic"] <= h.criteria["labels_spic"] + S.SWEEP_TIE
        else:
            assert not g.edges
    accepted = {g.hypotheses[0].label.split()[0] for g in groups
                if g.hypotheses[0].status == ident.Status.ACCEPTED}
    assert {"HI", "SiII", "SiIII", "AlII"} <= accepted
    # The regions an accepted metal explains are labelled with its lines
    sii = next(g for g in groups if g.hypotheses[0].label.startswith("SiII"))
    for rid in sii.regions[1:]:
        region = st.region(rid)
        assert region.kind == ident.RegionKind.METAL
        assert any(t.startswith("SiII") for t in region.transitions)


def test_without_the_ions_nothing_is_added():
    spec = H.identify_hypotheses(_rich(with_ions=False), ["CIV", "SiIV", "MgII", "Ly_a"]).spectrum
    added = {s for s, z, _ in _systems(spec) if not s.startswith("Ly_") and s != "CIV"}
    assert not added


def test_the_sweep_is_an_option():
    spec = _rich()
    off = H.identify_hypotheses(spec, ["CIV", "SiIV", "MgII", "Ly_a"], sweep=False).spectrum
    st = off.meta[ident.META_KEY]
    assert st.params["sweep"] is False
    assert not _sweep_groups(off)
    assert {s for s, _, _ in _systems(off)} == {"CIV", "Ly_a"}


def test_the_recipe_takes_sweep():
    session = Session("mock", None, spec=_rich())
    out = session.absorbers.identify_lines(multiplets="CIV,Ly_a", method="hypotheses",
                                           sweep=False)
    assert out.spec.meta[ident.META_KEY].params["sweep"] is False


def test_the_species_of_the_sweep_are_known():
    # Every line of a species of the sweep is in the table, and the
    # excited fine-structure lines are left out (F173)
    from astrocook.core.atomic_data import ATOM_DATA
    for series in S.SWEEP_SPECIES.values():
        lines = series_lines(series)
        assert lines and all(t in ATOM_DATA for t in lines)
        assert not any("*" in t for t in lines)
