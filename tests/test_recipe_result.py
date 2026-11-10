"""The structured result of a recipe (phase 6-bis, PR 4a).

``run`` runs a recipe of the registry and returns a RecipeResult: the
session, what changed (by a generic comparison of the session before and
after, Guido's choice), the warnings, the recipes run as steps, the values
of the parameters used.
"""
import json
import logging
import warnings

import numpy as np
import pytest

from astrocook.core.session import SessionV2
from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.recipes.errors import InvalidParameter, RecipeFailed
from astrocook.recipes.result import RecipeResult, compare_sessions, run

TRUTH = MockTruth(components=[MockComponent("CIV", 1.93, 13.8, 10.0)],
                  lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5, resol=45000.0,
                  snr=50.0, z_em=2.7, seed=1)


@pytest.fixture(scope="module")
def sessions():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw = SessionV2("mock", None, spec=generate(TRUTH, truth_columns=False))
        cont = raw.continuum.estimate_auto()
        fitted = cont.absorbers.doublets_auto()
    return raw, cont, fitted


# --- The comparison ----------------------------------------------------------

def test_nothing_changed(sessions):
    raw, _, _ = sessions
    changes = compare_sessions(raw, raw)
    assert changes.empty and changes.describe() == []


def test_columns_and_components(sessions):
    raw, cont, fitted = sessions
    changes = compare_sessions(raw, cont)
    assert changes.columns_added == ("abs_mask", "cont")
    assert not changes.columns_changed and not changes.components_added
    changes = compare_sessions(cont, fitted)
    (added,) = changes.components_added
    assert added.series == "CIV" and added.z == pytest.approx(1.93, abs=1e-4)
    assert type(added.z) is float  # plain values, for JSON
    back = compare_sessions(fitted, cont)
    assert [c.uuid for c in back.components_removed] == [added.uuid]


def test_changed_components_and_properties(sessions):
    _, _, fitted = sessions
    uuid = fitted.systs.components[0].uuid
    edited = fitted.absorbers.update_component(uuid, b=20.0)
    edited = edited.edit.set_properties(z_em=2.8, name="other")
    changes = compare_sessions(fitted, edited)
    (change,) = changes.components_changed
    assert change.component.uuid == uuid
    assert change.changed["b"][1] == 20.0
    assert changes.properties_changed["z_em"] == (2.7, 2.8)
    assert changes.properties_changed["name"] == ("mock", "other")
    assert "1 component changed" in changes.describe()


def test_a_refit_that_changed_nothing_changes_nothing(sessions):
    # Trial 3 of PR 7: a refit returns the values of the start to ~1e-10,
    # and every component of the group was listed as changed
    _, _, fitted = sessions
    uuid = fitted.systs.components[0].uuid
    c = fitted.systs.components[0]
    nudged = fitted.absorbers.update_component(
        uuid, z=c.z * (1 + 1e-10), logN=c.logN * (1 + 1e-10), b=c.b * (1 + 1e-10))
    assert not compare_sessions(fitted, nudged).components_changed
    # a real change, however small, is listed
    moved = fitted.absorbers.update_component(uuid, b=c.b * (1 + 1e-6))
    (change,) = compare_sessions(fitted, moved).components_changed
    assert set(change.changed) == {"b"}


def test_another_grid_is_said_by_the_pixels(sessions):
    raw, _, _ = sessions
    part = raw.edit.split("(λ > 445) & (λ < 450)")
    changes = compare_sessions(raw, part)
    assert not changes.columns_changed  # every column differs: not listed
    old, new = changes.properties_changed["pixels"]
    assert new < old


# --- run -------------------------------------------------------------------

def test_run_gives_the_session_the_changes_and_the_values_used(sessions):
    raw, _, _ = sessions
    result = run(raw, "continuum.estimate_auto", kappa="3")
    assert isinstance(result, RecipeResult)
    assert result.recipe == "continuum.estimate_auto" and not result.branching
    assert result.session.spec.has_aux_column("cont")
    assert result.params["kappa"] == 3.0 and result.params["smooth_fwhm"] == 800.0
    assert result.summary == ("Auto-estimate continuum (single-click): columns "
                              "'abs_mask', 'cont' added.")
    assert run(raw, "estimate_auto").recipe == "continuum.estimate_auto"  # by name


def test_the_steps_of_a_pipeline(sessions):
    _, cont, _ = sessions
    result = run(cont, "absorbers.doublets_auto")
    assert result.steps == ("absorbers.components_auto", "continuum.find_absorbed",
                            "absorbers.identify_lines",
                            "absorbers.populate_from_identification",
                            "absorbers.refit_all")
    assert len(result.changes.components_added) == 1
    # Outside run, the steps are not collected
    assert cont.absorbers.doublets_auto().systs is not None


def test_warnings_are_in_the_result(sessions):
    raw, _, _ = sessions
    result = run(raw, "absorbers.refit_all")
    assert result.session is raw
    assert result.warnings == ("There are no components to fit.",)
    assert result.summary == "Refit all systems: nothing changed, with 1 warning."


def test_several_sessions_and_a_value(sessions, tmp_path):
    raw, _, fitted = sessions
    result = run(raw, "edit.extract_preset", region="lya_forest, red_side")
    assert result.branching and len(result.sessions) == 2
    assert len(result.changes) == 2
    assert result.summary.startswith("Extract region by preset: 2 sessions")
    written = run(fitted, "file.export_ascii", targets="systems", path=str(tmp_path),
                  prefix="p")
    assert written.session is None and written.value is True and written.changes is None


def test_a_failure_is_the_exception_of_the_recipe(sessions):
    raw, _, _ = sessions
    with pytest.raises(InvalidParameter, match="kappa"):
        run(raw, "continuum.estimate_auto", kappa="x")
    with pytest.raises(RecipeFailed):
        run(raw, "flux.resample", target_session="missing")
    with pytest.raises(KeyError):
        run(raw, "no_such_recipe")


def test_the_fields_of_phase_7_are_there_and_empty(sessions):
    raw, _, _ = sessions
    result = run(raw, "flux.smooth", sigma_kms=10.0)
    assert (result.continuum_uncertainty, result.credibility, result.alternatives,
            dict(result.extra)) == (None, None, (), {})


def test_the_changes_can_be_written_as_json(sessions):
    import dataclasses
    _, cont, _ = sessions
    result = run(cont, "absorbers.doublets_auto")
    json.dumps(dataclasses.asdict(result.changes))


# --- The GUI --------------------------------------------------------------

def test_the_gui_worker_keeps_the_result(qtbot, sessions, caplog):
    from astrocook.gui.qt_workers import RecipeWorker
    raw, _, _ = sessions
    worker = RecipeWorker(raw, "flux", "smooth", {"sigma_kms": "10"}, {})
    done = []
    worker.signals.finished.connect(done.append)
    with caplog.at_level(logging.ERROR):
        worker.run()
    (state, name, _), = done
    assert worker.result.recipe == "flux.smooth" and worker.result.session is state
    # dy is constant on the mock, and stays the same once smoothed
    assert worker.result.changes.columns_changed == ("y",)
    assert np.all(np.isfinite(state.spec.y.value))


# --- What a recipe knows (PR 4b) ----------------------------------------------

def test_the_identification_keeps_region_and_score(sessions):
    # F49: the candidate systems kept neither their region nor their score
    _, cont, _ = sessions
    result = run(cont, "absorbers.identify_lines")
    (civ,) = [c for c in result.extra["candidates"] if c["series"] == "CIV"]
    assert civ["z"] == pytest.approx(1.93, abs=1e-4)
    assert civ["score"] > 0.9 and civ["region"] > 0
    regions = {r["region"]: r for r in result.extra["regions"]}
    region = regions[civ["region"]]
    assert region["x_unit"] == "nm"
    assert region["x_min"] < 154.8204 * (1 + civ["z"]) < region["x_max"]  # CIV 1548
    assert region["identifications"][0]["transition"].startswith("CIV")
    # Kept in the session too, with the other maps of the identification
    meta = result.session.spec.meta
    assert "region_map_json" in meta and "score_map_json" in meta


def test_the_reasons_of_what_is_removed(sessions):
    _, _, fitted = sessions
    uuid = fitted.systs.components[0].uuid
    deleted = run(fitted, "absorbers.delete_component", uuid=uuid)
    assert deleted.extra["removed"] == {uuid: "deleted"}
    cleaned = run(fitted, "absorbers.clean_negligible", min_logN=14.0)
    assert cleaned.extra["removed"] == {uuid: "logN < min_logN = 14.0 (too weak)"}
    assert [c.uuid for c in cleaned.changes.components_removed] == [uuid]
    nothing = run(fitted, "absorbers.clean_negligible")
    assert "removed" not in nothing.extra


def test_what_a_pipeline_makes_and_drops_is_discarded():
    # The narrow lines that lya_auto finds and removes were not in the
    # session: they are 'discarded', not 'removed'
    truth = MockTruth(components=[MockComponent("Ly_a", 2.65, 13.5, 25.0),
                                  MockComponent("Ly_a", 2.64, 14.0, 30.0)],
                      lam_min_nm=440.0, lam_max_nm=450.0, dv_kms=2.5, resol=45000.0,
                      snr=50.0, z_em=2.7, seed=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cont = SessionV2("m", None, spec=generate(truth, truth_columns=False)
                         ).continuum.estimate_auto()
        result = run(cont, "absorbers.lya_auto", min_b=40.0)
    assert "removed" not in result.extra
    assert set(result.extra["discarded"].values()) == {"narrower than min_b = 40.0 km/s"}


def test_the_groups_that_could_not_be_fitted_are_flagged(sessions, monkeypatch):
    # Agreed with Guido in phase 6: the flag not_fitted, set by refit_all on
    # the components of a group whose fit failed (they keep their values)
    from astrocook.fitting.voigt_fitter import VoigtFitterV2
    from astrocook.gui.inspector_table import flag_cue
    _, _, fitted = sessions

    def fit(self, *args, **kwargs):
        raise ValueError("x0 is infeasible")
    monkeypatch.setattr(VoigtFitterV2, "fit", fit)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = run(fitted, "absorbers.refit_all")
    (comp,) = result.session.systs.components
    assert "not_fitted" in comp.flags
    assert result.extra["not_fitted"] == [comp.uuid]
    assert comp.z == fitted.systs.components[0].z  # the values of before
    assert any("could not be fitted" in w for w in result.warnings)
    _, messages = flag_cue(comp)
    assert any("The fit of its group failed" in m for m in messages)


def test_report_outside_run_does_nothing():
    from astrocook.recipes.result import report, report_removed
    report("x", 1)
    report_removed(["u"], "why")  # no error, nothing kept


# --- What the user is told (PR 4c) ----------------------------------------------

def test_the_notes_of_a_result(sessions):
    _, cont, fitted = sessions
    uuid = fitted.systs.components[0].uuid
    cleaned = run(fitted, "absorbers.clean_negligible", min_logN=14.0)
    assert cleaned.notes == ("1 component removed: logN < min_logN = 14.0 (too weak)",)
    assert cleaned.counts == ("1 removed",)
    result = RecipeResult("absorbers.refit_all", "Refit all systems", None, False, None, {},
                          extra={"not_fitted": [uuid, "u2"],
                                 "discarded": {"a": "why", "b": "why", "c": "other"}})
    assert result.notes == ("2 components made and discarded by the recipe: why",
                            "1 component made and discarded by the recipe: other",
                            "2 components not fitted: the fit of their group failed, "
                            "the values are those from before")
    assert result.counts == ("3 discarded", "2 not fitted")
    identified = run(cont, "absorbers.identify_lines")
    n = len(identified.extra["candidates"])
    assert f"{n} candidates" in identified.counts
    assert run(cont, "flux.smooth", sigma_kms=10.0).notes == ()


def test_lya_auto_does_not_refit_the_red_side():
    # F102 (phase 6-bis, PR 4c, Guido's decision): the refit after the
    # removal of the narrow lines is limited to the forest, as the pipeline;
    # it refitted the CIV of the red side too, which changed slightly
    truth = MockTruth(components=[MockComponent("Ly_a", 2.65, 13.5, 25.0),
                                  MockComponent("Ly_a", 2.64, 14.0, 30.0),
                                  MockComponent("CIV", 1.93, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5, resol=45000.0,
                      snr=50.0, z_em=2.7, seed=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cont = SessionV2("m", None, spec=generate(truth, truth_columns=False)
                         ).continuum.estimate_auto()
        civ = cont.absorbers.doublets_auto()
        result = run(civ, "absorbers.lya_auto", min_b=40.0)
    assert result.extra["discarded"]  # the cleanup and its refit did run
    before = [(c.uuid, c.z, c.logN, c.b) for c in civ.systs.components]
    after = [(c.uuid, c.z, c.logN, c.b) for c in result.session.systs.components
             if c.series == "CIV"]
    assert after == before


def test_the_candidates_dropped_are_in_the_result():
    # F119 (trial 3 of PR 7): a doublet scored below the threshold was
    # dropped without a word; its score and threshold are in the result now
    truth = MockTruth(components=[MockComponent("CIV", 1.93, 13.8, 10.0),
                                  MockComponent("SiIV", 1.93, 12.4, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5, resol=45000.0,
                      snr=30.0, z_em=2.7, seed=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw = SessionV2("mock", None, spec=generate(truth, truth_columns=False))
        cont = raw.continuum.estimate_auto()
        result = run(cont, "identify_lines", score_threshold=0.999)
    dropped = result.extra["dropped_candidates"]
    assert dropped
    assert all(d["score"] <= d["threshold"] == 0.999 for d in dropped)
    assert all(set(d) == {"series", "z", "score", "threshold", "x_nm"} for d in dropped)


def test_optimize_system_says_what_it_tried_and_why_it_stopped(sessions):
    # F118 (trial 3 of PR 7): "1 component changed" and nothing else
    _, _, fitted = sessions
    uuid = fitted.systs.components[0].uuid
    result = run(fitted, "optimize_system", uuid=uuid)
    opt = result.extra["optimization"]
    assert isinstance(opt["start_aic"], float) and isinstance(opt["tried"], list)
    assert "no residual deeper than threshold_sigma = 2.5" in opt["stopped"]
    # a threshold low enough makes it try, and each trial says how it ended
    result = run(fitted, "optimize_system", uuid=uuid, threshold_sigma=0.1,
                 max_components=2)
    tried = result.extra["optimization"]["tried"]
    assert tried and all({"z", "residual_sigma", "outcome"} <= set(t) for t in tried)
    assert all(t["outcome"] == "added" or t["outcome"].startswith("rejected")
               or t["outcome"].endswith("not fitted") for t in tried)
