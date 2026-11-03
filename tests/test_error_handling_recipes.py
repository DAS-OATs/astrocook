"""Errors of the recipes reported once, with their context (phase 6, PR 5).

PR 5a: the recipes catch only the expected exceptions inside, their entry
points log an unexpected error once, with the traceback, and their
messages follow the rule of phase 6 (``ROADMAP.md``).
"""
import logging

import pytest

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
    with caplog.at_level(logging.ERROR):
        assert session.absorbers.update_constraint(
            b, "b", is_free=False, expression=f"sqrt(p['{a}'].b)", target_uuid=a) == 0
    assert "Could not set the constraint" in caplog.text
    assert "NameError" in caplog.text


# --- PR 5b: the rows D of the recipes -------------------------------------------------
# As decided with Guido (ROADMAP.md, phase 6, "Rows D of the recipes"); each
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
    # It was skipped with a WARNING: the whole spectrum was used
    session, _ = _two()
    with caplog.at_level(logging.ERROR):
        assert getattr(session.absorbers, recipe)(region_limit="nowhere") == 0
    assert any("The region 'nowhere' cannot be used" in m for m in _errors(caplog))


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
    assert session.absorbers.lya_auto(min_b="ten") == 0


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
    _failing_fit(monkeypatch, KeyError("bug"))
    assert session.absorbers.refit_all() == 0


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
    session = _session(systs=_systs(("CIV", 2.97, 13.6, 12.0)))
    _failing_fit(monkeypatch, KeyError("bug"))
    with pytest.raises(KeyError):
        session.absorbers.metals_auto(ions="SiIV")


def test_failed_delete_component_returns_0(monkeypatch):
    # It returned the session unchanged, which the GUI took as a success
    session, (a, _) = _two()

    def broken(self, **kwargs):
        raise RuntimeError("delete")
    monkeypatch.setattr(SystemListV2, "delete_component", broken)
    assert session.absorbers.delete_component(uuid=a) == 0


def test_failed_refit_of_an_import_keeps_the_import(monkeypatch, caplog):
    # Option (b), Guido's choice
    source = SimpleNamespace(systs=_systs(("CIV", 2.97, 13.6, 12.0)))
    session = _session(gui=_gui(src=source))
    monkeypatch.setattr(absorbers.RecipeAbsorbersV2, "refit_all", lambda self, **kw: 0)
    with caplog.at_level(logging.WARNING):
        out = session.edit.import_systems("src", refit=True)
    assert out != 0 and len(out.systs.components) == 1
    assert ("The imported systems could not be refitted: they keep their imported "
            "values.") in _warnings(caplog)


@pytest.mark.parametrize("kwargs, message", [
    (dict(manual_factors="one"), "The manual factors must be numbers"),
    (dict(manual_factors="1.0, 2.0"), "Give 1 manual factors"),
    (dict(equalize_ranges="600"), "The equalize ranges must be written"),
    (dict(equalize_ranges="600-610, 620-630"), "Give 1 equalize ranges"),
])
def test_equalize_and_stitch_refuses_parameters_it_cannot_read(kwargs, message, caplog):
    # They were replaced by 'auto', with a WARNING
    session = _session(gui=_gui(red=SimpleNamespace(spec=_spec(620.0, 628.0))))
    with caplog.at_level(logging.ERROR):
        assert session.edit.equalize_and_stitch("red", "620.0", **kwargs) == 0
    assert any(message in m for m in _errors(caplog))


def test_trim_common_without_a_common_range_fails(caplog):
    # A WARNING and the session unchanged
    far = SimpleNamespace(spec=_spec(300.0, 301.0, 0.005))
    session = _session(gui=_gui(far=far))
    with caplog.at_level(logging.ERROR):
        assert session.edit.trim_common("far", "2.97", "CIV_1548", "Ly_a") == 0
    assert "The sessions have no range in common: nothing is trimmed." in _errors(caplog)


# --- PR 5c: the prerequisites of the fits ---------------------------------------------
# The continuum estimated first when missing (and said so), its failure
# handled; a fit without a resolution said at WARNING (it was at DEBUG).

import numpy as np  # noqa: E402
import astropy.units as au  # noqa: E402

import astrocook.recipes.continuum as continuum  # noqa: E402
from astrocook.recipes.messages import (msg_continuum_estimated,  # noqa: E402
                                        msg_no_resolution)

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


@pytest.mark.parametrize("name", FIT_RECIPES)
def test_failed_estimate_of_the_continuum_stops_the_fit(name, monkeypatch, caplog):
    # The 0 of the failed estimate went on as the session: "'int' object has
    # no attribute 'systs'"
    session, uuid = _without_continuum()
    monkeypatch.setattr(continuum.RecipeContinuumV2, "estimate_auto", lambda self: 0)
    with caplog.at_level(logging.INFO):
        assert FIT_RECIPES[name](session, uuid) == 0
    messages = [r.getMessage() for r in caplog.records]
    assert msg_continuum_estimated in messages
    assert "The continuum could not be estimated: the recipe is not run." in messages
    assert not any("'int' object" in m for m in messages)


def test_failed_estimate_of_the_continuum_stops_populate(monkeypatch):
    session, _ = _without_continuum()
    monkeypatch.setattr(continuum.RecipeContinuumV2, "estimate_auto", lambda self: 0)
    assert session.absorbers.populate_from_identification(
        series_map_json='{"1": "CIV"}', z_map_json='{"1": 2.97}') == 0


def test_refit_all_says_it_estimates_the_continuum(caplog):
    # It estimated it without saying so
    session, _ = _without_continuum()
    with caplog.at_level(logging.INFO):
        assert session.absorbers.refit_all() != 0
    assert msg_continuum_estimated in [r.getMessage() for r in caplog.records]


@pytest.mark.parametrize("name", FIT_RECIPES)
def test_fit_without_a_resolution_is_a_warning(name, caplog):
    session, (uuid, _) = _two()
    with caplog.at_level(logging.WARNING):
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
