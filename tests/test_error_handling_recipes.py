"""Errors of the recipes reported once, with their context (phase 6, PR 5).

PR 5a: the recipes catch only the expected exceptions inside, their entry
points log an unexpected error once, with the traceback, and their
messages follow the rule of phase 6 (``ROADMAP.md``).
"""
import logging

import pytest
from astrocook.recipes.errors import InvalidParameter, RecipeFailed  # phase 6-bis, PR 3
from astrocook.recipes.params import typed_recipe

from astrocook.core.session import SessionV2

from test_recipes_absorbers import _systs
from test_recipes_edit_flux import _spec


def _two():
    session = SessionV2("s", None, spec=_spec(),
                        systs=_systs(("CIV", 2.0, 13.0, 10.0), ("SiIV", 2.0001, 12.5, 8.0)))
    return session, [c.uuid for c in session.systs.components]


def test_linked_value_uses_the_functions_of_the_fit(caplog):
    # F98: the update of the linked value offered 'sqrt' and not 'ln', the
    # fit the opposite; both use LINK_FUNCTIONS now
    session, (a, b) = _two()
    with caplog.at_level(logging.WARNING):
        new = session.absorbers.update_constraint(
            b, "logN", is_free=False, expression=f"ln(p['{a}'].logN)", target_uuid=a)
    assert new != 0
    assert new.systs.get_component_by_uuid(b).logN == pytest.approx(2.564949357)
    assert caplog.text == ""


def test_link_with_sqrt_is_refused(caplog):
    session, (a, b) = _two()
    with pytest.raises(RecipeFailed, match="Update constraint: The link of b") as info:
        session.absorbers.update_constraint(
            b, "b", is_free=False, expression=f"sqrt(p['{a}'].b)", target_uuid=a)
    assert "NameError" in str(info.value)


# --- PR 5b: the rows D of the recipes -------------------------------------------------
# As decided with Guido (ROADMAP_HISTORY.md, phase 6, "Rows D of the recipes"); each
# test reproduces the old fallback, which went on with a result that
# looked valid.

from types import SimpleNamespace  # noqa: E402

import astrocook.recipes.absorbers as absorbers  # noqa: E402
from astrocook.core.system_list import SystemListV2  # noqa: E402
from astrocook.fitting.voigt_fitter import VoigtFitterV2  # noqa: E402

from test_recipes_edit_flux import _gui, _session  # noqa: E402


def _errors(caplog):
    return [r.getMessage() for r in caplog.records if r.levelno >= logging.ERROR]


def _warnings(caplog):
    return [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]


@pytest.mark.parametrize("recipe", ["components_auto", "refit_all"])
def test_unknown_region_fails(recipe, caplog):
    # It was skipped with a WARNING: the whole spectrum was used; a
    # RecipeFailed in phase 6, refused as a parameter since phase 6-bis
    # (PR 6c: a list of choices)
    from astrocook.recipes.errors import InvalidParameter
    session, _ = _two()
    with pytest.raises(InvalidParameter, match="region_limit: 'nowhere' is not one of"):
        getattr(session.absorbers, recipe)(region_limit="nowhere")


def test_failed_cleaning_keeps_the_result_of_lya_auto(monkeypatch, caplog):
    # Option (b), Guido's choice: the result is kept, the WARNING says what
    # it still contains
    session = _session(systs=_systs(("Ly_a", 4.0, 13.0, 5.0), ("Ly_a", 4.001, 13.5, 30.0)))
    monkeypatch.setattr(absorbers.RecipeAbsorbersV2, "components_auto",
                        lambda self, **kwargs: session)

    def broken(*args, **kwargs):
        raise RuntimeError("cleaning")
    monkeypatch.setattr(absorbers, "_remove_narrow", broken)
    with caplog.at_level(logging.WARNING):
        out = session.absorbers.lya_auto(min_b="10.0")
    assert out != 0 and len(out.systs.components) == 2
    assert ("The lines narrower than 10.0 km/s could not be removed "
            "(RuntimeError: cleaning): 1 of them are in the result.") in _warnings(caplog)


def test_unreadable_min_b_of_lya_auto_fails(monkeypatch, caplog):
    session = _session(systs=_systs(("Ly_a", 4.0, 13.0, 5.0)))
    monkeypatch.setattr(absorbers.RecipeAbsorbersV2, "components_auto",
                        lambda self, **kwargs: session)
    with pytest.raises(InvalidParameter, match="min_b: 'ten' is not a number"):
        session.absorbers.lya_auto(min_b="ten")


def _failing_fit(monkeypatch, error):
    def fit(self, *args, **kwargs):
        raise error
    monkeypatch.setattr(VoigtFitterV2, "fit", fit)


def test_failed_groups_of_refit_all_are_counted(monkeypatch, caplog):
    session, _ = _two()
    _failing_fit(monkeypatch, ValueError("x0 is infeasible"))
    with caplog.at_level(logging.WARNING):
        out = session.absorbers.refit_all()
    assert out != 0
    [warning] = [m for m in _warnings(caplog) if "could not be fitted" in m]
    assert warning == ("2 of 2 groups could not be fitted (z = 2.00000, 2.00010): "
                       "their components keep their previous values.")
    assert [c.z for c in out.systs.components] == [2.0, 2.0001]


def test_unexpected_error_in_refit_all_fails_the_recipe(monkeypatch):
    session, _ = _two()
    # A KeyError is turned into RecipeFailed at the boundary (phase 6-bis, PR 3)
    _failing_fit(monkeypatch, KeyError("bug"))
    with pytest.raises(RecipeFailed, match="Refit all systems: 'bug' not found"):
        session.absorbers.refit_all()


def test_failed_trial_fits_of_metals_auto_are_counted(monkeypatch, caplog):
    session = _session(systs=_systs(("CIV", 2.97, 13.6, 12.0)))
    _failing_fit(monkeypatch, ValueError("x0 is infeasible"))
    with caplog.at_level(logging.INFO):
        out = session.absorbers.metals_auto(ions="SiIV,NV")
    assert out != 0 and len(out.systs.components) == 1
    assert "2 trial fits failed: those candidates were not added." in _warnings(caplog)


def test_unexpected_error_in_metals_auto_reaches_the_caller(monkeypatch):
    # metals_auto has no boundary of its own: the error reaches the caller
    # (in the GUI, the worker, which shows it)
    # A KeyError is one of the errors of the core that the boundary of the
    # recipes turns into RecipeFailed (Guido, phase 6-bis, PR 3), with the
    # original as its cause; another error (a bug) goes on as it is, with a
    # note naming the recipe
    session = _session(systs=_systs(("CIV", 2.97, 13.6, 12.0)))
    _failing_fit(monkeypatch, KeyError("bug"))
    with pytest.raises(RecipeFailed, match="Check for associated metals: 'bug' not found") as info:
        session.absorbers.metals_auto(ions="SiIV")
    assert isinstance(info.value.__cause__, KeyError)
    _failing_fit(monkeypatch, AttributeError("bug"))
    with pytest.raises(AttributeError) as info:
        session.absorbers.metals_auto(ions="SiIV")
    assert "In the recipe absorbers.metals_auto" in info.value.__notes__[-1]


def test_failed_delete_component_fails(monkeypatch):
    # It returned the session unchanged, which the GUI took as a success;
    # since phase 6-bis (PR 3) the error goes on, naming the recipe (an
    # error that is not one of the core's is a bug, shown with its traceback)
    session, (a, _) = _two()

    def broken(self, **kwargs):
        raise RuntimeError("delete")
    monkeypatch.setattr(SystemListV2, "delete_component", broken)
    with pytest.raises(RuntimeError) as info:
        session.absorbers.delete_component(uuid=a)
    assert "In the recipe absorbers.delete_component" in info.value.__notes__[-1]


def test_failed_refit_of_an_import_keeps_the_import(monkeypatch, caplog):
    # Option (b), Guido's choice
    source = SimpleNamespace(systs=_systs(("CIV", 2.97, 13.6, 12.0)))
    session = _session(gui=_gui(src=source))
    def failing(self, **kwargs):
        raise RecipeFailed("the fit did not converge")
    monkeypatch.setattr(absorbers.RecipeAbsorbersV2, "refit_all", failing)
    with caplog.at_level(logging.WARNING):
        out = session.edit.import_systems("src", refit=True)
    assert len(out.systs.components) == 1
    assert ("The imported systems could not be refitted (the fit did not converge): "
            "they keep their imported values.") in _warnings(caplog)


@pytest.mark.parametrize("kwargs, message", [
    (dict(manual_factors="one"), "manual_factors: give numbers separated by commas"),
    (dict(manual_factors="1.0, 2.0"), "manual_factors: give 1, one per other session"),
    (dict(equalize_ranges="600"), "equalize_ranges: write them 'min-max, min-max'"),
    (dict(equalize_ranges="600-610, 620-630"), "equalize_ranges: give 1, one per other"),
])
def test_equalize_and_stitch_refuses_parameters_it_cannot_read(kwargs, message):
    # They were replaced by 'auto', with a WARNING; since phase 6-bis (PR 3)
    # an InvalidParameter that names them
    session = _session(gui=_gui(red=SimpleNamespace(spec=_spec(620.0, 628.0))))
    with pytest.raises(InvalidParameter) as info:
        session.edit.equalize_and_stitch("red", "620.0", **kwargs)
    assert message in str(info.value)


def test_trim_common_without_a_common_range_fails(caplog):
    # A WARNING and the session unchanged
    far = SimpleNamespace(spec=_spec(300.0, 301.0, 0.005))
    session = _session(gui=_gui(far=far))
    with pytest.raises(RecipeFailed, match="The sessions have no range in common: "
                                           "nothing is trimmed."):
        session.edit.trim_common("far", "2.97", "CIV_1548", "Ly_a")


# --- PR 5c: the prerequisites of the fits ---------------------------------------------
# The continuum estimated first when missing (and said so), its failure
# handled; a fit without a resolution said at WARNING (it was at DEBUG).

import numpy as np  # noqa: E402
import astropy.units as au  # noqa: E402

import astrocook.recipes.continuum as continuum  # noqa: E402
from astrocook.recipes.messages import msg_no_resolution  # noqa: E402

# Said when the continuum is estimated first (phase 6-bis, PR 5: by the
# boundary of the recipe, at WARNING; it was msg_continuum_estimated, INFO)
CONTINUUM_ESTIMATED = ("The spectrum has no continuum: it is estimated first with "
                       "Auto-estimate continuum (single-click), with its defaults.")

FIT_RECIPES = {
    "refit_all": lambda s, uuid: s.absorbers.refit_all(),
    "fit_component": lambda s, uuid: s.absorbers.fit_component(uuid),
    "optimize_system": lambda s, uuid: s.absorbers.optimize_system(uuid, max_components="1"),
}


def _without_continuum():
    # A physical flux (median far from 1): without 'cont' there is no norm
    spec = _spec(scale=100.0).remove_columns(["cont", "model"])
    assert spec.norm is None
    session = _session(spec, systs=_systs(("CIV", 2.97, 13.6, 12.0)))
    return session, session.systs.components[0].uuid


@typed_recipe
def estimate_auto(self):
    raise RecipeFailed("no unabsorbed pixels")


_failing_estimate = estimate_auto


@pytest.mark.parametrize("name", FIT_RECIPES)
def test_failed_estimate_of_the_continuum_stops_the_fit(name, monkeypatch, caplog):
    # The 0 of the failed estimate went on as the session: "'int' object has
    # no attribute 'systs'"
    session, uuid = _without_continuum()
    monkeypatch.setattr(continuum.RecipeContinuumV2, "estimate_auto", _failing_estimate)
    with caplog.at_level(logging.INFO), pytest.raises(RecipeFailed) as info:
        FIT_RECIPES[name](session, uuid)
    assert CONTINUUM_ESTIMATED in [r.getMessage() for r in caplog.records]
    # The recipe, the step and the reason (phase 6-bis, PR 3)
    assert str(info.value).endswith(" → Auto-estimate continuum (single-click): "
                                    "no unabsorbed pixels.")
    assert info.value.chain[0][1] == name


def test_failed_estimate_of_the_continuum_stops_populate(monkeypatch):
    session, _ = _without_continuum()
    monkeypatch.setattr(continuum.RecipeContinuumV2, "estimate_auto", _failing_estimate)
    with pytest.raises(RecipeFailed, match="Populate components from identification → "
                                           "Auto-estimate continuum"):
        session.absorbers.populate_from_identification(
            series_map_json='{"1": "CIV"}', z_map_json='{"1": 2.97}')


def test_refit_all_says_it_estimates_the_continuum(caplog):
    # It estimated it without saying so
    session, _ = _without_continuum()
    with caplog.at_level(logging.INFO):
        assert session.absorbers.refit_all() != 0
    assert CONTINUUM_ESTIMATED in [r.getMessage() for r in caplog.records]


@pytest.mark.parametrize("name", FIT_RECIPES)
def test_fit_without_a_resolution_is_a_warning(name, caplog):
    session, (uuid, _) = _two()
    with caplog.at_level(logging.WARNING):
        if name == "fit_component":
            # The components of _two are outside the spectrum: refused
            # since F111 (phase 6-ter, PR 7-quater), after the warning
            with pytest.raises(RecipeFailed, match="No valid pixels around the components"):
                FIT_RECIPES[name](session, uuid)
        else:
            FIT_RECIPES[name](session, uuid)
    assert msg_no_resolution in _warnings(caplog)


def test_fit_with_a_resolution_has_no_such_warning(caplog):
    session, (uuid, _) = _two()
    spec = session.spec.update_column("resol", np.full(session.spec.x.size, 45000.0),
                                      unit=au.dimensionless_unscaled)
    session = session.with_new_spectrum(spec)
    with caplog.at_level(logging.WARNING):
        session.absorbers.refit_all()
    assert msg_no_resolution not in _warnings(caplog)


# --- Trial 3 of the MCP server (phase 6-ter, PR 7; F117) ---------------------------
# The example of the description, ``p["z"]``, cannot be evaluated, and with
# ``is_free`` left out the link was stored unchecked (the check ran only with
# ``is_free=False``): the fit then kept the parameter frozen, with a line of
# the log that nobody read

def test_link_without_is_free_is_checked(caplog):
    session, (a, b) = _two()
    with pytest.raises(RecipeFailed, match="Update constraint: The link of z") as info:
        session.absorbers.update_constraint(
            b, "z", expression='p["z"]', target_uuid=a)
    assert "KeyError" in str(info.value)
    # nothing is stored: the session is the one it was
    assert (b, "z") not in session.systs.parsed_constraints
    assert caplog.text == ""


def test_link_without_is_free_that_works_is_kept():
    # unchanged: a valid link is a link, free or not given
    session, (a, b) = _two()
    new = session.absorbers.update_constraint(
        b, "z", expression=f"p['{a}'].z", target_uuid=a)
    assert new.systs.get_component_by_uuid(b).z == pytest.approx(2.0)
    assert new.systs.parsed_constraints[(b, "z")].expression == f"p['{a}'].z"


def test_unlinking_does_not_check_the_old_expression():
    # is_free=True wipes the link: no check, whatever was passed
    session, (a, b) = _two()
    linked = session.absorbers.update_constraint(
        b, "z", is_free=False, expression=f"p['{a}'].z", target_uuid=a)
    freed = linked.absorbers.update_constraint(b, "z", is_free=True)
    assert freed.systs.parsed_constraints[(b, "z")].expression is None
