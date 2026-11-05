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
