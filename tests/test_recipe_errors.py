"""The exceptions of the recipes (phase 6-bis, PR 3).

A recipe that cannot do what was asked raises a RecipeError, which names
the recipe (and the steps of a pipeline) and the reason; it logs nothing,
and the error is logged once where it is caught: the worker of the GUI, the
worker of the scripts, the MCP server (Guido's decisions at the start of
PR 3). The errors of the core that mean "this data cannot be used"
(ValueError, KeyError, OSError) become RecipeFailed at the boundary of the
recipe; the others (bugs) go on, with a note naming the recipe.
"""
import logging
from types import SimpleNamespace

import pytest

from astrocook.session.session import Session
from astrocook.recipes.errors import (InvalidParameter, MissingPrerequisite, RecipeError,
                                      RecipeFailed, reason_of)
from astrocook.recipes.params import typed_recipe
from test_recipes_edit_flux import _session, _spec


# --- The message ---------------------------------------------------------

def test_message_without_and_with_the_recipes():
    error = RecipeFailed("the region 'nope' is unknown")
    assert str(error) == "the region 'nope' is unknown."
    assert error.recipe is None
    error.enter("absorbers", "refit_all", "Refit all systems")
    error.enter("absorbers", "doublets_auto", "Auto-fit metal doublets (Red Side)")
    assert str(error) == ("Auto-fit metal doublets (Red Side) → Refit all systems: "
                          "the region 'nope' is unknown.")
    assert error.recipe == "absorbers.refit_all"  # where it happened
    assert [c[1] for c in error.chain] == ["doublets_auto", "refit_all"]
    assert str(RecipeFailed("Done!")) == "Done!"  # no second full stop


def test_the_kinds():
    assert issubclass(InvalidParameter, RecipeError) and issubclass(InvalidParameter, ValueError)
    assert issubclass(RecipeFailed, RecipeError) and not issubclass(RecipeFailed, ValueError)
    assert issubclass(MissingPrerequisite, RecipeError)
    assert reason_of(KeyError("abs_mask")) == "'abs_mask' not found"
    assert reason_of(ValueError("bad")) == "bad"
    assert reason_of(ValueError()) == "ValueError"


# --- The boundary of the recipes -----------------------------------------

def test_errors_of_the_core_become_recipe_failed():
    # An error of the core (here, a column that the spectrum does not have)
    session = _session()
    with pytest.raises(RecipeFailed) as info:
        session.flux.calculate_running_std(input_col="nope")
    assert str(info.value).startswith("Calculate running standard deviation: ")
    assert info.value.recipe == "flux.calculate_running_std"
    assert info.value.__cause__ is not None  # the original error, for the traceback


def test_a_bug_goes_on_with_a_note(monkeypatch):
    from astrocook.core.spectrum import Spectrum

    def broken(self, sigma_kms):
        raise AttributeError("a bug")
    monkeypatch.setattr(Spectrum, "smooth", broken)
    with pytest.raises(AttributeError) as info:
        _session().flux.smooth(sigma_kms=10.0)
    assert info.value.__notes__ == ["In the recipe flux.smooth (Smooth spectrum)."]


def test_a_pipeline_names_its_step(monkeypatch):
    # The error of a recipe called by another names both, outermost first
    from astrocook.recipes import continuum

    @typed_recipe
    def estimate_auto(self):
        raise RecipeFailed("no unabsorbed pixels")
    monkeypatch.setattr(continuum.RecipeContinuum, "estimate_auto", estimate_auto)
    spec = _spec(scale=100.0).remove_columns(["cont", "model"])
    session = Session(name="s", gui=None, spec=spec)
    with pytest.raises(RecipeFailed) as info:
        session.absorbers.populate_from_identification(
            series_map_json='{"1": "CIV"}', z_map_json='{"1": 2.97}')
    assert str(info.value) == ("Populate components from identification → "
                               "Auto-estimate continuum (single-click): no unabsorbed pixels.")
    assert info.value.recipe == "continuum.estimate_auto"


def test_a_missing_prerequisite():
    spec = _spec(z_em=0.0)
    with pytest.raises(MissingPrerequisite, match="^Fit power-law continuum: The emission "
                                                   "redshift"):
        Session(name="s", gui=None, spec=spec).continuum.fit_powerlaw()


def test_the_recipe_logs_nothing(caplog):
    with caplog.at_level(logging.DEBUG, logger="astrocook.recipes"):
        with pytest.raises(RecipeFailed):
            _session(gui=SimpleNamespace(session_histories=[])).flux.resample("missing")
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


# --- Where it is caught --------------------------------------------------

def test_the_gui_worker_shows_and_logs_it_once(qtbot, caplog):
    from astrocook.gui.qt_workers import RecipeWorker
    worker = RecipeWorker(_session(), "flux", "smooth", {"sigma_kms": "fifty"}, {})
    errors, done = [], []
    worker.signals.error.connect(errors.append)
    worker.signals.finished.connect(done.append)
    with caplog.at_level(logging.ERROR):
        worker.run()  # synchronously, in this thread
    assert done == []
    assert errors == [("Recipe Aborted",
                       "Smooth spectrum: sigma_kms: 'fifty' is not a number.", "")]
    logged = [r.getMessage() for r in caplog.records if r.levelno >= logging.ERROR]
    assert logged == ["Smooth spectrum: sigma_kms: 'fifty' is not a number."]


def test_the_mcp_server_passes_it_to_the_agent():
    pytest.importorskip("mcp")
    from mcp.server.mcpserver.exceptions import ToolError
    from astrocook.mcp_server.errors import run_recipe
    session = _session()
    with pytest.raises(ToolError) as info:
        run_recipe("flux_smooth", session.flux.smooth, sigma_kms="fifty")
    assert str(info.value) == "Smooth spectrum: sigma_kms: 'fifty' is not a number."


# --- Absorbers (PR 3c): nothing to do, a pipeline whose fit fails --------

def test_nothing_to_do_is_not_a_failure(caplog):
    # Guido's decision: the session back, with a WARNING (it was a failure, 0)
    session = _session()
    with caplog.at_level(logging.WARNING):
        assert session.absorbers.refit_all() is session
        assert session.absorbers.metals_auto() is session
        assert session.absorbers.clean_negligible() is session
    warned = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert "There are no components to fit." in warned
    assert "There are no systems to search metals at: nothing is added." in warned
    assert "There are no components to clean." in warned


@pytest.mark.parametrize("call, message", [
    (lambda s: s.absorbers.fit_component("x"), "Fit component: There are no components to fit."),
    (lambda s: s.absorbers.optimize_system("x"), "There are no components to optimize."),
    (lambda s: s.absorbers.update_constraints_column("b", False), "There are no components."),
    # "No components are selected" before phase 6-quater, PR 10-bis: the
    # session is checked first, as by the other fits
    (lambda s: s.absorbers.fit_bayesian([]), "Fit Bayesian: There are no components to fit."),
])
def test_without_components_the_reason_is_given(call, message):
    # It returned 0, logging nothing (Guido, phase 6-bis, PR 3)
    with pytest.raises(RecipeFailed) as info:
        call(_session())
    assert message in str(info.value)


def test_a_pipeline_whose_fit_fails_fails(monkeypatch):
    # The 0 of a failed refit was taken for "no components to fit", and the
    # pipeline gave back the components unfitted, as a result (Guido, PR 3c)
    from astrocook.mocks import MockComponent, MockTruth, generate
    from astrocook.recipes import absorbers

    @typed_recipe
    def refit_all(self, **kwargs):
        raise RecipeFailed("the fit did not converge")
    monkeypatch.setattr(absorbers.RecipeAbsorbers, "refit_all", refit_all)
    truth = MockTruth(components=[MockComponent("CIV", 1.93, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5, resol=45000.0,
                      snr=50.0, z_em=2.7, seed=1)
    session = Session("mock", None, spec=generate(truth, truth_columns=False))
    with pytest.raises(RecipeFailed) as info:
        session.continuum.estimate_auto().absorbers.doublets_auto()
    assert str(info.value) == ("Auto-fit metal doublets (Red Side) → Auto-detect and fit "
                               "components → Refit all systems: the fit did not converge.")
