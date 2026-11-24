"""The structure of the identification and the structured data of the
archive (phase 7, PR 7c.3; ``design/identification.md``, section 9; F28,
F93 item 1)."""
import dataclasses
import json
import logging
import tarfile

import numpy as np
import pytest

from astrocook.core import identification as ident
from astrocook.core.identification import (
    Component, Edge, Group, Hypothesis, Identification, Region, RegionKind, Stage, Status,
)
from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.mocks.physics import C_KMS
from astrocook.session.session import Session

from full_chain import run_chain

Z_EM = 3.0
COMPONENTS = [
    MockComponent("CIV", 2.50, 14.0, 10.0),
    MockComponent("CIV", 2.70, 13.8, 12.0),
    MockComponent("CIV", 2.70 + 40.0 / C_KMS * 3.70, 13.6, 10.0),
    MockComponent("SiIV", 2.65, 13.5, 10.0),
    MockComponent("MgII", 1.20, 13.5, 8.0),
    MockComponent("Ly_a", 2.85, 14.0, 25.0),
]


@pytest.fixture(scope="module")
def session():
    truth = MockTruth(components=COMPONENTS, lam_min_nm=440.0, lam_max_nm=700.0,
                      dv_kms=2.5, resol=45000.0, snr=100.0, z_em=Z_EM, seed=1)
    return Session("mock", None, spec=generate(truth))


@pytest.fixture(scope="module")
def identified(session):
    return session.absorbers.identify_lines(multiplets="CIV,SiIV,MgII,Ly_ab",
                                            score_threshold=0.5)


def _structure(s) -> Identification:
    return s.spec._data.meta[ident.META_KEY]


# --- The structure of the current method -----------------------------------

def test_the_current_method_fills_the_structure(identified):
    st = _structure(identified)
    assert st.method == "doublets"
    assert st.params["score_threshold"] == 0.5
    assert st.params["multiplets"] == ["CIV", "SiIV", "MgII", "Ly_ab"]
    # Every region of abs_ids is in the structure, with the same ID
    ids = identified.spec.get_column("abs_ids").value
    in_column = set(np.unique(ids[ids > 0]).astype(int).tolist())
    assert in_column <= {r.id for r in st.regions}
    assert len({r.id for r in st.regions}) == len(st.regions)
    # F93 item 1: the pairs are stored, as edges, with their score
    pairs = {e.series for g in st.groups for e in g.edges}
    assert {"CIV", "SiIV", "MgII"} <= pairs
    for g in st.groups:
        for e in g.edges:
            assert set(e.regions) <= set(g.regions)
            assert e.prefilter is not None and e.prefilter > 0.5
            assert e.lag_kms is not None and e.lag_kms > 0


def test_the_kinds_of_the_regions(identified):
    st = _structure(identified)
    for r in st.regions:
        if r.kind == RegionKind.METAL:
            assert any(t.split("_")[0] in ("CIV", "SiIV", "MgII") for t in r.transitions)
        if r.kind == RegionKind.FOREST_ASSUMED:
            assert r.transitions == ("Ly_a",)
        if r.kind == RegionKind.UNIDENTIFIED:
            assert r.transitions == ()
    kinds = {r.kind for r in st.regions}
    assert RegionKind.METAL in kinds and RegionKind.FOREST_ASSUMED in kinds
    # The labels of region_identifications and of the structure agree
    labels = json.loads(identified.spec.meta["region_identifications"])
    for rid, cands in labels.items():
        assert st.region(int(rid)).transitions == tuple(dict.fromkeys(c[0] for c in cands))


def test_every_candidate_system_is_in_one_accepted_hypothesis(identified):
    st = _structure(identified)
    series = {int(k): v for k, v in json.loads(identified.spec.meta["series_map_json"]).items()}
    z = {int(k): v for k, v in json.loads(identified.spec.meta["z_map_json"]).items()}
    seen = [c for h in st.hypotheses(Status.ACCEPTED) for c in h.components]
    assert sorted(c.candidate for c in seen) == sorted(series)
    for c in seen:
        assert (c.series, c.z) == (series[c.candidate], z[c.candidate])
    for h in st.hypotheses():
        assert h.stage == Stage.SCORE and h.status == Status.ACCEPTED and h.reason
    # The two CIV components 40 km/s apart: one pair, two components
    civ = [h for h in st.hypotheses() if h.label.startswith("CIV at z = 2.70")]
    assert len(civ) == 1 and len(civ[0].components) == 2


def test_the_dropped_candidates_are_in_the_structure(session):
    out = session.absorbers.identify_lines(multiplets="CIV,SiIV,MgII", score_threshold=0.999)
    st = _structure(out)
    assert st.dropped
    from astrocook.recipes.absorbers import _identification_extra
    extra = _identification_extra(out.spec)["dropped_candidates"]
    assert [(d.series, d.z, d.score) for d in st.dropped] == \
           [(d["series"], d["z"], d["score"]) for d in extra]


def test_the_numbers_of_the_identification_are_unchanged(session, identified):
    # The maps of the header are left as they were (PR 7c.3, decision c):
    # the structure is added next to them
    meta = identified.spec.meta
    for key in ("series_map_json", "z_map_json", "logN_map_json", "b_map_json",
                "region_map_json", "score_map_json", "rejected_candidates_json",
                "region_identifications"):
        assert isinstance(meta[key], str)


def test_populate_gives_the_uuids_of_the_components(identified):
    populated = identified.absorbers.populate_from_identification(
        series_map_json=identified.spec.meta["series_map_json"],
        z_map_json=identified.spec.meta["z_map_json"])
    st = _structure(populated)
    uuids = {c.uuid for c in populated.systs.components}
    comps = [c for h in st.hypotheses(Status.ACCEPTED) for c in h.components]
    assert comps and all(c.uuid in uuids for c in comps)
    # One component per candidate system
    assert len({c.uuid for c in comps}) == len({c.candidate for c in comps})


# --- JSON ---------------------------------------------------------------

def test_json_round_trip(identified):
    st = _structure(identified)
    d = json.loads(json.dumps(st.to_dict()))
    assert d["schema"] == ident.SCHEMA_NAME and d["schema_version"] == ident.SCHEMA_VERSION
    assert Identification.from_dict(d) == st


def _full():
    comp = Component("CIV", 2.5, 13.5, 10.0, 1e-5, 0.05, 1.0, candidate=1, uuid="u1")
    hyps = (Hypothesis("CIV at z = 2.5", (1, 2), (comp,), Stage.FIT, Status.ACCEPTED,
                       "SpIC_H lower by 12", 0.97, 101.5, 3, 120,
                       {"aic": 107.5, "bic": 115.9, "spic_h": 110.2}),
            Hypothesis("two Ly_a", (1, 2), (), Stage.AOD, Status.ALTERNATIVE,
                       "", 0.03))
    group = Group(1, (1, 2), (Edge((1, 2), ("CIV_1548", "CIV_1550"), "CIV", 2.5,
                                   497.9, 0.9),), hyps)
    # A test of the sweep of the anchored systems (PR 7c.4d)
    swept = Group(2, (1, 3), (Edge((1, 3), ("CIV_1548", "AlII_1670"), "AlII_1670", 2.5),),
                  (Hypothesis("AlII at z = 2.50000 (sweep of CIV)", (1, 3), (), Stage.FIT,
                              Status.ACCEPTED, "spic -3.0", 0.82, anchors=(0,),
                              criteria={"spic": 50.0, "replaced_unidentified": 1.0}),),
                  anchored_on=1)
    return Identification("hypotheses", {"criterion": "spic"},
                          (Region(1, 541.0, 541.2, RegionKind.METAL, ("CIV_1548",), 5.2),
                           Region(2, 542.0, 542.2, RegionKind.METAL, ("CIV_1550",), 3.1),
                           Region(3, 400.0, 400.1, RegionKind.TELLURIC)),
                          (group, swept),
                          (ident.DroppedCandidate("SiIV", 2.9, 0.3, 0.5, 530.0, 530.1),))


def test_json_round_trip_of_every_field():
    st = _full()
    assert Identification.from_dict(json.loads(json.dumps(st.to_dict()))) == st


def test_a_field_unknown_to_this_version_is_ignored():
    d = _full().to_dict()
    d["new_field"] = 1
    d["regions"][0]["new_field"] = 2
    assert Identification.from_dict(d) == _full()


@pytest.mark.parametrize("change, match", [
    ({"schema_version": ident.SCHEMA_VERSION + 1}, "is not known"),
    ({"schema": "something else"}, "Not an identification"),
])
def test_a_later_or_other_schema_is_refused(change, match):
    with pytest.raises(ValueError, match=match):
        Identification.from_dict({**_full().to_dict(), **change})


def test_the_vocabulary_maps_to_the_truth_map():
    from astrocook.mocks.truthmap import CATEGORIES
    assert set(ident.TRUTH_CATEGORY) == set(RegionKind)
    assert {c for c in ident.TRUTH_CATEGORY.values() if c} == set(CATEGORIES)


# --- Merge of the steps of a pipeline ---------------------------------------

def test_merged_follows_merge_region_ids():
    from astrocook.recipes.absorbers import merge_region_ids
    old_ids = np.array([0, 1, 1, 0, 2, 2, 0, 0, 0, 0])
    new_ids = np.array([0, 0, 0, 0, 1, 1, 0, 2, 2, 0])
    old = Identification("doublets", {}, (
        Region(1, 1, 2, RegionKind.METAL, ("CIV_1548",)),
        Region(2, 4, 5, RegionKind.METAL, ("CIV_1550",)),
        Region(3, 8, 9, RegionKind.UNIDENTIFIED)),
        (Group(1, (1, 2), (Edge((1, 2), ("CIV_1548", "CIV_1550"), "CIV", 2.0),), ()),))
    new = Identification("doublets", {"run": 2}, (
        Region(1, 4, 5, RegionKind.FOREST_ASSUMED, ("Ly_a",)),
        Region(2, 7, 8, RegionKind.FOREST_ASSUMED, ("Ly_a",)),
        Region(3, 9.5, 9.6, RegionKind.UNIDENTIFIED)),
        (Group(1, (2,), (), (Hypothesis("Ly_a (assumed)", (2,), (), Stage.SCORE,
                                        Status.ACCEPTED),)),
         Group(2, (2,), (), (), anchored_on=1)))
    final_ids, _ = merge_region_ids(old_ids, {}, new_ids, {})
    merged = old.merged(new, int(old_ids.max()), np.unique(old_ids), np.unique(new_ids))
    ids = [r.id for r in merged.regions]
    assert len(set(ids)) == len(ids)
    # The regions of the column keep the IDs of the merged column
    by_id = {r.id: r for r in merged.regions}
    for rid in np.unique(final_ids[final_ids > 0]):
        assert int(rid) in by_id
    assert by_id[4].x_min == 7  # region 2 of the new run, shifted by 2
    assert merged.groups[1].regions == (4,) and merged.groups[1].id == 2
    # A group of the sweep follows the group it was anchored on (PR 7c.4d)
    assert merged.groups[2].id == 3 and merged.groups[2].anchored_on == 2
    assert merged.groups[0].anchored_on is None
    assert merged.groups[0].edges[0].regions == (1, 2)
    # The unlabelled regions are numbered past the others
    assert sorted(r.id for r in merged.regions if r.kind == RegionKind.UNIDENTIFIED) == [5, 6]
    assert merged.params == {"run": 2}


@pytest.fixture(scope="module")
def chain(session):
    return run_chain(session)


def test_the_structure_after_the_pipelines(chain):
    for state in chain[1:]:
        st = _structure(state)
        ids = state.spec.get_column("abs_ids").value
        in_column = set(np.unique(ids[ids > 0]).astype(int).tolist())
        regions = {r.id: r for r in st.regions}
        assert len(regions) == len(st.regions)
        assert in_column <= set(regions)
        for rid in in_column:
            assert regions[rid].kind != RegionKind.UNIDENTIFIED
        for g in st.groups:
            assert set(g.regions) <= set(regions)
            for h in g.hypotheses:
                assert set(h.regions) <= set(g.regions)
    # The doublets of the first step are still there after the forest
    assert any(e.series == "CIV" for g in _structure(chain[2]).groups for e in g.edges)


# --- The archive -----------------------------------------------------------

def test_the_identification_in_the_archive(identified, tmp_path):
    path = tmp_path / "s.acs2"
    identified.save(str(path), initial_session=identified)
    with tarfile.open(path) as tar:
        names = tar.getnames()
    assert {"s_identification.json", "s_identification_orig.json"} <= set(names)
    reloaded = Session.open_new(str(path), "s", None, "auto")
    assert _structure(reloaded) == _structure(identified)
    # Not in the FITS header
    from astropy.io import fits
    with tarfile.open(path) as tar:
        tar.extract("s_spec.fits", tmp_path, filter="data")
    assert ident.META_KEY not in fits.getheader(tmp_path / "s_spec.fits", 1)


def test_an_unreadable_file_is_left_out_with_a_warning(identified, tmp_path, caplog):
    path = tmp_path / "s.acs2"
    identified.save(str(path))
    work = tmp_path / "work"
    with tarfile.open(path) as tar:
        tar.extractall(work, filter="data")
        names = tar.getnames()
    d = json.loads((work / "s_identification.json").read_text())
    d["schema_version"] = ident.SCHEMA_VERSION + 1
    (work / "s_identification.json").write_text(json.dumps(d))
    with tarfile.open(path, "w:gz") as tar:
        for n in names:
            tar.add(work / n, arcname=n)
    with caplog.at_level(logging.WARNING):
        reloaded = Session.open_new(str(path), "s", None, "auto")
    assert reloaded != 0
    assert ident.META_KEY not in reloaded.spec._data.meta
    assert "the identification of the regions, but it cannot be read" in caplog.text
    np.testing.assert_array_equal(reloaded.spec.get_column("abs_ids").value,
                                  identified.spec.get_column("abs_ids").value)


def test_an_archive_without_the_files_opens_as_before(identified, tmp_path):
    session = identified.with_new_spectrum(identified.spec.remove_meta_keys([ident.META_KEY]))
    path = tmp_path / "old.acs2"
    session.save(str(path))
    with tarfile.open(path) as tar:
        assert not any("identification" in n for n in tar.getnames())
    reloaded = Session.open_new(str(path), "old", None, "auto")
    assert ident.META_KEY not in reloaded.spec._data.meta
    assert reloaded.spec.meta["series_map_json"] == session.spec.meta["series_map_json"]


def test_the_structure_follows_the_undo(identified):
    from astrocook.session.history import SessionHistory
    from astrocook.core.structures import HistoryLog
    history = SessionHistory(identified, HistoryLog())
    later = identified.with_new_spectrum(identified.spec.remove_meta_keys([ident.META_KEY]))
    history.add_state(later)
    assert ident.META_KEY not in history.current_state.spec._data.meta
    history.undo()
    assert history.current_state.spec._data.meta[ident.META_KEY] == _structure(identified)


def test_the_structure_is_frozen(identified):
    st = _structure(identified)
    with pytest.raises(dataclasses.FrozenInstanceError):
        st.method = "other"  # type: ignore[misc]


def test_a_second_save_has_no_spurious_warning(identified, tmp_path, caplog):
    # F165: after a round trip region_identifications is a dict, and the
    # save warned that it could not be saved (it is, in _meta.json)
    identified.save(str(tmp_path / "a.acs2"))
    reloaded = Session.open_new(str(tmp_path / "a.acs2"), "a", None, "auto")
    with caplog.at_level(logging.WARNING):
        reloaded.save(str(tmp_path / "b.acs2"))
    assert "cannot be saved" not in caplog.text
    again = Session.open_new(str(tmp_path / "b.acs2"), "b", None, "auto")
    assert again.spec.meta["region_identifications"] == \
           reloaded.spec.meta["region_identifications"]
    assert _structure(again) == _structure(identified)
