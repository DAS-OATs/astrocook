"""The prerequisites of the recipes (phase 6-bis, PR 5).

The registry declares what each recipe requires (z_em, the resolution, the
continuum) and what it provides; ``recipes/prerequisites.py`` says what is
missing and fills it at the boundary of the recipe, as the policy of the
script says: ``provide`` (the default: the continuum estimated first, with
a WARNING; z_em and the resolution left to the recipe, as before) or
``raise``. The step run to provide is reported in the result, and is a
state and a log entry of its own in the GUI.
"""
import warnings

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.session import SessionV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, SpectrumDataV2
from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.recipes import prerequisites
from astrocook.recipes.errors import MissingPrerequisite
from astrocook.recipes.prerequisites import missing, policy, providers_of
from astrocook.recipes.registry import REQUIREMENTS, all_recipes, find_recipe
from astrocook.recipes.result import run

FLUX = au.erg / au.cm**2 / au.s / au.nm
TRUTH = MockTruth(components=[MockComponent("CIV", 1.93, 13.8, 10.0)],
                  lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5, resol=45000.0,
                  snr=50.0, z_em=2.7, seed=1)


def _physical(mock):
    """The mock in physical units, on a power law: a spectrum without a
    continuum (the mock itself is normalized, so its continuum is 1)."""
    d = mock._data
    k = 1e-16 * (mock.x.value / 450.0) ** -1.5
    return SpectrumV2(SpectrumDataV2(
        x=d.x, xmin=d.xmin, xmax=d.xmax, y=DataColumnV2(mock.y.value * k, FLUX),
        dy=DataColumnV2(mock.dy.value * k, FLUX), aux_cols=dict(d.aux_cols), z_em=d.z_em))


@pytest.fixture(scope="module")
def raw():
    return SessionV2("raw", None, spec=_physical(generate(TRUTH, truth_columns=False)))


def _components(session):
    return [(c.series, c.z, c.logN, c.b) for c in session.systs.components]


# --- The registry ---------------------------------------------------------------

def test_every_requirement_has_a_provider():
    for name in REQUIREMENTS:
        assert providers_of(name), name
        automatic = prerequisites.PREREQUISITES[name].automatic
        if automatic is not None:
            assert providers_of(name)[0] == automatic
            category, recipe = automatic.split('.')
            assert name in find_recipe(recipe).provides
            # It runs with its defaults: no parameter without one
            assert all(p.default not in ("", None) or p.optional
                       for p in find_recipe(recipe).params)
    assert providers_of("z_em") == ("edit.set_properties",)


def test_what_needs_the_continuum():
    # The recipes that estimated it themselves (populate_from_identification,
    # fit_component, optimize_system, refit_all), those that fit, and the
    # identification with the pipelines (F103, Guido)
    needs = {s.name for s in all_recipes() if "continuum" in s.requires}
    assert needs == {"identify_lines", "populate_from_identification", "components_auto",
                     "doublets_auto", "lya_auto", "forest_auto", "metals_auto",
                     "add_component", "add_linked_system", "fit_component",
                     "optimize_system", "fit_bayesian", "refit_all"}
    # Every recipe that fits needs the resolution too (F99)
    for name in needs - {"identify_lines", "populate_from_identification"}:
        assert "resolution" in find_recipe(name).requires, name


# --- What is missing ----------------------------------------------------------------

def test_missing_says_what_and_who(raw):
    (gap,) = missing(raw, find_recipe("refit_all"))
    assert gap.name == "continuum" and gap.automatic == "continuum.estimate_auto"
    assert gap.message == ("The spectrum has no continuum: Auto-estimate continuum "
                           "(single-click) (continuum.estimate_auto) provides it with its "
                           "defaults")
    assert "continuum.fit_continuum" in gap.providers
    no_z = raw.edit.set_properties(z_em=0.0)
    gaps = missing(no_z, find_recipe("doublets_auto"))
    assert [g.name for g in gaps] == ["z_em", "continuum"]
    assert gaps[0].message == ("The emission redshift (z_em) is not set: set it with "
                               "Set session properties (edit.set_properties)")
    assert missing(raw, find_recipe("smooth")) == ()


def test_a_normalized_spectrum_has_its_continuum():
    # A spectrum with a flux of ~1 has a continuum of 1 (SpectrumV2.norm)
    session = SessionV2("norm", None, spec=generate(TRUTH, truth_columns=False))
    assert missing(session, find_recipe("refit_all")) == ()


# --- The policy "provide" -----------------------------------------------------

def test_provided_first_as_the_replay_will(raw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = run(raw, "absorbers.doublets_auto")
        replay = raw.continuum.estimate_auto().absorbers.doublets_auto()
    (step,) = result.provided
    assert step.recipe == "continuum.estimate_auto"
    assert step.reason == "the spectrum has no continuum"
    assert step.params["kappa"] == 2.0  # its defaults, for the log
    assert step.session.spec.norm is not None
    assert result.warnings[0] == ("The spectrum has no continuum: it is estimated first "
                                  "with Auto-estimate continuum (single-click), with its "
                                  "defaults.")
    assert result.notes[0] == ("Auto-estimate continuum (single-click) run first: the "
                               "spectrum has no continuum")
    # The log has the two steps: their replay gives the same result
    assert _components(result.session) == _components(replay)
    assert len(result.session.systs.components) == 1  # the CIV, found (F103)


def test_only_the_outermost_recipe_provides(raw):
    # The steps of a pipeline find the continuum already there
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = run(raw, "absorbers.lya_auto")
    assert len(result.provided) == 1
    assert result.steps.count("continuum.estimate_auto") == 1


def test_without_z_em_the_recipe_goes_on_as_before(raw):
    # No recipe can set z_em by itself: with "provide" the recipe goes on,
    # with its own WARNING (Guido: the numbers do not change)
    no_z = raw.edit.set_properties(z_em=0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = run(no_z, "continuum.estimate_auto")
    assert result.provided == ()
    assert any("z_em" in w for w in result.warnings)


# --- The policy "raise" --------------------------------------------------------

def test_raise_names_what_is_missing(raw):
    with pytest.raises(MissingPrerequisite) as info:
        run(raw, "absorbers.refit_all", on_missing="raise")
    assert str(info.value).startswith("Refit all systems: The spectrum has no continuum: "
                                      "Auto-estimate continuum")
    with policy("raise"), pytest.raises(MissingPrerequisite, match="z_em"):
        raw.edit.set_properties(z_em=0.0).continuum.estimate_auto()
    # Outside the block, "provide" again
    assert prerequisites.current_policy() == "provide"


def test_raise_lets_a_complete_session_run(raw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cont = raw.continuum.estimate_auto()
        with policy("raise"):
            fitted = cont.absorbers.doublets_auto()
    assert len(fitted.systs.components) == 1


def test_an_unknown_policy():
    with pytest.raises(ValueError, match="Unknown policy 'ask'"):
        with policy("ask"):
            pass


def test_the_region_of_a_refit_needs_z_em(raw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fitted = raw.continuum.estimate_auto().absorbers.doublets_auto()
    with pytest.raises(MissingPrerequisite, match="The region 'red_side' needs the emission"):
        fitted.edit.set_properties(z_em=0.0).absorbers.refit_all(region_limit="red_side")


# --- The script console -----------------------------------------------------------

def test_the_script_worker_logs_the_step_provided(raw, qtbot):
    from astrocook.gui.qt_workers import ScriptWorker

    class Window:
        recipe_category_map = {"refit_all": "absorbers"}
        active_history = None

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fitted = raw.continuum.estimate_auto().absorbers.doublets_auto()
    start = raw.with_new_system_list(fitted.systs)
    worker = ScriptWorker("refit_all()", start, Window())
    done, errors = [], []
    worker.signals.finished.connect(done.append)
    worker.signals.error.connect(errors.append)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        worker.run()
    assert errors == []
    (history,) = done
    assert [e.recipe_name for e in history.log_manager.entries] == ["estimate_auto",
                                                                    "refit_all"]
    assert len(history.states) == 3
    assert history.states[1].spec.norm is not None
    assert np.array_equal(history.states[1].spec.cont.value,
                          history.current_state.spec.cont.value)
