"""Characterization tests of the recipe layer of ``recipes/absorbers.py``
and of the system list operations behind it.

The automatic pipelines (identification, populate, ``refit_all``,
``doublets_auto``, ``lya_auto``) are covered by the identification and
full-chain snapshots; here the target is the single-component recipes of
the system inspector and of the Absorbers menu: parameter parsing (string
inputs as from the GUI dialogs), what they do to the system list and to
the model column, and what is returned on bad input.

Input: the synthetic normalised spectrum of ``test_snapshots_voigt.py``
(a blended CIV doublet with two components and a Ly-alpha line, R = 45000,
noise 0.02). Fitted values are compared with the known input within loose
tolerances, not stored as snapshots: they depend on the machine (see
``test_snapshots_full_chain.py``).

These tests record the current behaviour. F24 (``fit_bayesian`` with at
most 200 MCMC steps overwrote the component parameters with NaN) and F25
(components could not be deleted from a spectrum without a continuum) are
fixed and covered by regression tests (findings register, ``ROADMAP.md``,
section 3.2).
"""
import numpy as np
import pytest

from astrocook.core.atomic_data import STANDARD_MULTIPLETS
from astrocook.core.session import SessionV2
from astrocook.core.structures import ComponentDataV2, SystemListDataV2
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2

from test_snapshots_voigt import (
    R_SYNTH, TRUE_COMPONENTS, _model, _synthetic_spectrum, _system_list,
)

C_KMS = 299792.458
CIV_Z = 2.9700


@pytest.fixture(scope="module")
def data():
    """Noisy synthetic spectrum with the components of TRUE_COMPONENTS."""
    clean = _synthetic_spectrum(resol=R_SYNTH)
    true_model = _model(clean, _system_list(TRUE_COMPONENTS))
    rng = np.random.default_rng(42)
    return _synthetic_spectrum(
        y=true_model + rng.normal(0.0, 0.02, len(true_model)), resol=R_SYNTH)


def _session(spec, systs=None):
    return SessionV2(name="synth", gui=None, spec=spec, systs=systs)


def _systs(*params):
    """System list from (series, z, logN, b) tuples."""
    comps = [ComponentDataV2(id=i + 1, z=z, dz=None, logN=logN, dlogN=None,
                             b=b, db=None, series=series)
             for i, (series, z, logN, b) in enumerate(params)]
    return SystemListV2(SystemListDataV2(components=comps))


@pytest.fixture(scope="module")
def one_civ(data):
    """Session with one CIV component added (and fitted) by the recipe."""
    return _session(data).absorbers.add_component("CIV", str(CIV_Z + 1e-4),
                                                  "13.5", "10.0")


def _dv(z1, z2):
    return abs(z1 - z2) / (1 + z2) * C_KMS


@pytest.fixture
def restore_multiplets():
    saved = dict(STANDARD_MULTIPLETS)
    yield
    STANDARD_MULTIPLETS.clear()
    STANDARD_MULTIPLETS.update(saved)


# --- Add, update, delete --------------------------------------------------

def test_add_component_fits_it(one_civ):
    (comp,) = one_civ.systs.components
    assert comp.series == "CIV"
    # One component for the two blended ones (logN 13.6 + 13.1, b 12 and 6)
    assert _dv(comp.z, CIV_Z) < 3.0
    assert comp.logN == pytest.approx(13.72, abs=0.1)
    assert comp.dz is not None and np.isfinite(comp.dz)
    # The model column is updated
    model = one_civ.spec.get_column("model").value
    assert model.min() < 0.7


@pytest.mark.parametrize("kwargs", [dict(z="abc"), dict(logN=""),
                                    dict(z_window_kms="x")])
def test_add_component_bad_input_returns_0(data, kwargs):
    assert _session(data).absorbers.add_component("CIV", **{"z": "2.97", **kwargs}) == 0


def test_add_component_custom_multiplet(data, restore_multiplets):
    series = "CIV_1548,SiIV_1393"
    out = _session(data).absorbers.add_component(series, "2.97")
    assert STANDARD_MULTIPLETS[series] == ["CIV_1548", "SiIV_1393"]
    assert out.systs.components[0].series == series


def test_update_component(one_civ):
    uuid = one_civ.systs.components[0].uuid
    out = one_civ.absorbers.update_component(uuid, logN="14.0", b="20")
    comp = out.systs.components[0]
    assert (comp.logN, comp.b) == (14.0, 20.0)
    assert comp.z == one_civ.systs.components[0].z  # 'None' leaves it as is
    assert one_civ.absorbers.update_component(uuid, z="abc") == 0
    # Unknown UUID: logged, list unchanged
    same = one_civ.absorbers.update_component("unknown", logN="14.0")
    assert same.systs.components == one_civ.systs.components


def test_update_components(data):
    session = _session(data, _systs(("CIV", 2.97, 13.0, 10.0),
                                    ("Ly_a", 4.059, 13.0, 20.0)))
    uuids = [c.uuid for c in session.systs.components]
    out = session.absorbers.update_components(uuids, b="15.0")
    assert [c.b for c in out.systs.components] == [15.0, 15.0]
    assert session.absorbers.update_components(uuids) is session


def test_delete_component_recomputes_model(data):
    session = _session(data, _systs(("CIV", 2.97, 13.6, 12.0),
                                    ("Ly_a", 4.059, 13.8, 25.0)))
    civ, lya = session.systs.components
    out = session.absorbers.delete_component(uuid=civ.uuid)
    assert [c.uuid for c in out.systs.components] == [lya.uuid]
    _, expected = VoigtFitterV2(data, out.systs).compute_model_flux()
    np.testing.assert_allclose(out.spec.get_column("model").value, expected)
    both = session.absorbers.delete_component(uuids=[civ.uuid, lya.uuid])
    assert both.systs.components == []


def test_delete_component_without_continuum(data):
    # Regression test for F25: without a continuum the model could not be
    # recomputed, the error was logged and the component was not deleted
    import astropy.units as au
    from test_recipes_edit_flux import _spec
    # Physical flux (the fitter takes a median flux between 0.1 and 10 as
    # already normalised)
    session = _session(_spec(scale=1e-16), _systs(("CIV", 2.97, 13.6, 12.0),
                                                   ("CIV", 2.98, 11.0, 12.0)))
    session = session.edit.delete("cont")
    assert session.spec.y.unit != au.dimensionless_unscaled
    assert session.spec.norm is None
    session = session.with_new_spectrum(
        session.spec.update_model(np.ones_like(session.spec.y.value)))
    out = session.absorbers.delete_component(uuid=session.systs.components[0].uuid)
    assert [c.z for c in out.systs.components] == [2.98]
    # The model, which still included the deleted component, is dropped
    assert not out.spec.has_aux_column("model")
    # "Clean negligible" goes through delete_component
    cleaned = session.absorbers.clean_negligible(min_logN="11.5")
    assert [c.z for c in cleaned.systs.components] == [2.97]


def test_delete_component_resets_links_to_it(data):
    session = _session(data, _systs(("CIV", 2.97, 13.6, 12.0),
                                    ("SiIV", 2.97, 13.0, 12.0)))
    civ, siiv = session.systs.components
    linked = session.absorbers.update_constraint(
        siiv.uuid, "z", is_free=False, expression=f"p['{civ.uuid}'].z",
        target_uuid=civ.uuid)
    out = linked.absorbers.delete_component(uuid=civ.uuid)
    constraint = out.systs._data.v2_constraints_map[siiv.uuid]["z"]
    assert constraint.is_free and constraint.expression is None


# --- Fit, constraints -----------------------------------------------------

@pytest.mark.parametrize("kwargs", [dict(max_nfev="x"), dict(z_window_kms="x"),
                                    dict(group_depth="x")])
def test_fit_component_bad_input_returns_0(one_civ, kwargs):
    uuid = one_civ.systs.components[0].uuid
    assert one_civ.absorbers.fit_component(uuid, **kwargs) == 0


def test_frozen_parameter_is_not_fitted(one_civ):
    uuid = one_civ.systs.components[0].uuid
    moved = one_civ.absorbers.update_component(uuid, b="20.0")
    frozen = moved.absorbers.update_constraint(uuid, "b", is_free=False)
    out = frozen.absorbers.fit_component(uuid)
    comp = out.systs.components[0]
    assert comp.b == 20.0
    assert comp.logN != one_civ.systs.components[0].logN


def test_link_updates_the_value(data):
    session = _session(data, _systs(("CIV", 2.97, 13.6, 12.0),
                                    ("SiIV", 2.98, 13.0, 12.0)))
    civ, siiv = session.systs.components
    out = session.absorbers.update_constraint(
        siiv.uuid, "z", is_free=False, expression=f"p['{civ.uuid}'].z",
        target_uuid=civ.uuid)
    assert out.systs.components[1].z == 2.97
    # Unlinking clears the expression
    free = out.absorbers.update_constraint(siiv.uuid, "z", is_free=True)
    constraint = free.systs._data.v2_constraints_map[siiv.uuid]["z"]
    assert constraint.is_free and constraint.expression is None


def test_update_constraints_column(data):
    session = _session(data, _systs(("CIV", 2.97, 13.6, 12.0),
                                    ("Ly_a", 4.059, 13.8, 25.0)))
    out = session.absorbers.update_constraints_column("b", is_free=False)
    for comp in out.systs.components:
        assert not out.systs._data.v2_constraints_map[comp.uuid]["b"].is_free


def test_add_linked_system(data):
    out = _session(data).absorbers.add_linked_system("CIV, SiIV", CIV_Z)
    civ, siiv = out.systs.components
    assert (civ.series, siiv.series) == ("CIV", "SiIV")
    assert siiv.z == civ.z  # linked to the primary, also after the fit
    assert _dv(civ.z, CIV_Z) < 3.0
    constraint = out.systs._data.v2_constraints_map[siiv.uuid]["z"]
    assert constraint.target_uuid == civ.uuid and not constraint.is_free
    assert _session(data).absorbers.add_linked_system("", CIV_Z) == 0


# --- Optimize, metals, clean ----------------------------------------------

def test_optimize_system_adds_the_blended_component(one_civ):
    uuid = one_civ.systs.components[0].uuid
    out = one_civ.absorbers.optimize_system(uuid, threshold_sigma="1.5")
    civ = sorted((c for c in out.systs.components if c.series == "CIV"),
                 key=lambda c: c.z)
    assert len(civ) == 2
    for comp, true in zip(civ, TRUE_COMPONENTS[:2]):
        assert _dv(comp.z, true["z"]) < 3.0
        assert comp.logN == pytest.approx(true["logN"], abs=0.2)
    # With the default threshold (2.5 sigma) no component is added
    default = one_civ.absorbers.optimize_system(uuid)
    assert len(default.systs.components) == 1
    assert one_civ.absorbers.optimize_system(uuid, max_components="x") == 0


def test_metals_auto(one_civ, data):
    # No other ion is kept at the redshift of the CIV system
    out = one_civ.absorbers.metals_auto()
    assert [c.series for c in out.systs.components] == ["CIV"]
    assert _session(data).absorbers.metals_auto() == 0
    assert one_civ.absorbers.metals_auto(ions="CIV", aic_penalty="x") == 0


@pytest.mark.parametrize("combined, removed", [
    ("True", ["narrow", "weak", "broad_weak"]),
    ("False", ["narrow", "weak", "broad_weak", "broad_strong"]),
])
def test_clean_negligible(data, combined, removed):
    params = {"narrow": ("CIV", 2.97, 13.5, 2.0),
              "weak": ("CIV", 2.971, 11.0, 10.0),
              "broad_weak": ("Ly_a", 4.05, 12.5, 100.0),
              "broad_strong": ("Ly_a", 4.06, 14.0, 100.0),
              "good": ("Ly_a", 4.059, 13.8, 25.0)}
    session = _session(data, _systs(*params.values()))
    out = session.absorbers.clean_negligible(combined_check=combined)
    kept = {c.z for c in out.systs.components}
    assert kept == {p[1] for name, p in params.items() if name not in removed}


def test_clean_negligible_nothing_to_do(data):
    session = _session(data, _systs(("Ly_a", 4.059, 13.8, 25.0)))
    assert session.absorbers.clean_negligible() is session
    assert session.absorbers.clean_negligible(min_logN="x") == 0


# --- Reproducibility of refit_all (F19) ----------------------------------

def test_refit_all_does_not_depend_on_the_uuids(data):
    # Regression test for F19: with group_depth >= 1 the seeds of the groups
    # were taken from a set of random UUIDs, so the depth-limited groups (and
    # the fitted values) changed with the UUIDs and from run to run.
    import dataclasses
    # A chain of blended Ly-alpha lines around the one in the spectrum:
    # neighbours overlap, so a depth-1 group depends on its seed
    z0, dz = 4.0590, 25.0 / C_KMS * (1 + 4.0590)
    params = [("Ly_a", z0 + k * dz, 13.0, 20.0) for k in range(-4, 5)]
    results = []
    for labels in range(6):
        rng = np.random.default_rng(labels)
        comps = [dataclasses.replace(c, uuid=f"{rng.integers(1 << 62):x}")
                 for c in _systs(*params).components]
        systs = SystemListV2(SystemListDataV2(components=comps))
        out = _session(data, systs).absorbers.refit_all(
            max_nfev="30", group_depth="1", max_group_size="3")
        results.append(sorted((c.z, c.logN, c.b) for c in out.systs.components))
    assert all(r == results[0] for r in results[1:])


# --- group_depth (F11) ----------------------------------------------------

def _lya_chain(n=7):
    # 60 km/s apart, b = 20 km/s: only adjacent lines overlap (threshold
    # max(4 b, 30 km/s) = 80 km/s), so each depth step adds two neighbours
    z0, dz = 4.0590, 60.0 / C_KMS * (1 + 4.0590)
    return _systs(*[("Ly_a", z0 + k * dz, 13.0, 20.0) for k in range(-(n // 2), n // 2 + 1)])


def test_set_active_components_honours_group_depth():
    # Regression test for F11: the depth was ignored and 2 was always used
    systs = _lya_chain()
    centre = systs.components[3].uuid
    model = systs.constraint_model
    sizes = {}
    for depth in (0, 1, 2):
        model.set_active_components([centre], group_depth=depth)
        sizes[depth] = len(model._active_uuids)
    model.set_active_components(None)
    assert sizes[0] == 1
    assert sizes[0] < sizes[1] < sizes[2]


def test_fit_component_group_depth_0_freezes_neighbours(data):
    session = _session(data, _lya_chain())
    centre = session.systs.components[3]
    out = session.absorbers.fit_component(centre.uuid, group_depth="0")
    changed = [c.uuid for c, old in zip(out.systs.components, session.systs.components)
               if (c.z, c.logN, c.b) != (old.z, old.logN, old.b)]
    assert changed == [centre.uuid]
    # The default (2) fits the neighbours too
    out2 = session.absorbers.fit_component(centre.uuid)
    changed2 = [c for c, old in zip(out2.systs.components, session.systs.components)
                if (c.z, c.logN, c.b) != (old.z, old.logN, old.b)]
    assert len(changed2) > 1


def test_refit_all_errors_come_from_the_final_model(data):
    # Issue #534: refit_all fits each unit while its neighbours still have
    # their initial values; the errors and chi2 are recomputed on the final
    # model, so they match an evaluation at the final parameters
    import dataclasses
    from astrocook.fitting.voigt_fitter import VoigtFitterV2
    from astrocook.recipes.absorbers import FIT_EXPANSION_DEPTH
    out = _session(data, _lya_chain()).absorbers.refit_all(
        max_nfev="50", group_depth="1", max_group_size="3")
    systs = out.systs
    target = systs.components[3]
    local = SystemListV2(dataclasses.replace(systs._data, components=list(systs.components)))
    with local.fitting_context([target.uuid], group_depth=FIT_EXPANSION_DEPTH):
        evaluated = VoigtFitterV2(out.spec, local).evaluate_uncertainties(100.0)
    again = next(c for c in evaluated.components if c.uuid == target.uuid)
    assert again.chi2 == pytest.approx(target.chi2, rel=0.2)
    assert again.dlogN == pytest.approx(target.dlogN, rel=0.2)


# --- Upper bound on b (F32) ---------------------------------------------

def test_refit_all_max_b(one_civ):
    out = one_civ.absorbers.refit_all(max_b="8.0")
    assert all(c.b <= 8.0 + 1e-9 for c in out.systs.components)
    assert one_civ.absorbers.refit_all(max_b="x") == 0


def test_pipelines_bound_b_and_manual_fits_do_not():
    import inspect
    from astrocook.recipes.absorbers import RecipeAbsorbersV2, ABSORBERS_RECIPES_SCHEMAS
    for name in ("components_auto", "doublets_auto", "lya_auto", "forest_auto",
                 "metals_auto"):
        sig = inspect.signature(getattr(RecipeAbsorbersV2, name))
        assert sig.parameters["max_b"].default == "50.0", name
        schema = {p["name"]: p for p in ABSORBERS_RECIPES_SCHEMAS[name]["params"]}
        assert schema["max_b"]["default"] == 50.0, name
    sig = inspect.signature(RecipeAbsorbersV2.refit_all)
    assert sig.parameters["max_b"].default == "300.0"
    assert "max_b" not in inspect.signature(RecipeAbsorbersV2.fit_component).parameters


# --- Bayesian, anchor -----------------------------------------------------

def test_fit_bayesian_mcmc(one_civ):
    uuid = one_civ.systs.components[0].uuid
    out = one_civ.absorbers.fit_bayesian([uuid], nsteps="300")
    comp = out.systs.components[0]
    ref = one_civ.systs.components[0]
    assert _dv(comp.z, ref.z) < 1.0
    assert comp.logN == pytest.approx(ref.logN, abs=0.05)
    assert np.isfinite(comp.dlogN) and comp.dlogN > 0
    stored = out.systs._data.meta["bayesian_results"][uuid]
    assert stored["engine"] == "mcmc" and len(stored["samples"]) > 0


@pytest.mark.parametrize("nsteps", ["200", "10"])
def test_fit_bayesian_few_steps_keeps_parameters(one_civ, nsteps):
    # Regression test for F24: the burn-in was a fixed 200 steps, so a chain
    # of at most 200 steps was empty and the NaN medians overwrote the
    # components and the model.
    uuid = one_civ.systs.components[0].uuid
    out = one_civ.absorbers.fit_bayesian([uuid], nsteps=nsteps)
    comp = out.systs.components[0]
    assert np.isfinite([comp.z, comp.logN, comp.b]).all()
    assert np.isfinite(out.spec.model.value).all()
    stored = out.systs._data.meta["bayesian_results"][uuid]
    assert len(stored["samples"]) > 0


def test_fit_bayesian_non_finite_result_returns_0(one_civ, monkeypatch):
    from astrocook.fitting.bayesian_fitter import BayesianVoigtFitter

    def nan_mcmc(self, **kwargs):
        return {"samples": np.empty((0, self.ndim)),
                "medians": np.full(self.ndim, np.nan),
                "std": np.full(self.ndim, np.nan)}

    monkeypatch.setattr(BayesianVoigtFitter, "run_mcmc", nan_mcmc)
    uuid = one_civ.systs.components[0].uuid
    assert one_civ.absorbers.fit_bayesian([uuid], nsteps="500") == 0


def test_run_mcmc_rejects_an_empty_chain(one_civ):
    from astrocook.fitting.bayesian_fitter import BayesianVoigtFitter
    uuid = one_civ.systs.components[0].uuid
    with one_civ.systs.fitting_context([uuid], group_depth=0):
        b_fitter = BayesianVoigtFitter(one_civ.spec, one_civ.systs)
        with pytest.raises(ValueError):
            b_fitter.run_mcmc(nsteps=0, progress=False)
        with pytest.raises(ValueError):
            b_fitter.run_mcmc(nsteps=100, burn_in=100, progress=False)


def test_fit_bayesian_bad_input_returns_0(one_civ):
    uuid = one_civ.systs.components[0].uuid
    assert one_civ.absorbers.fit_bayesian(None) == 0
    assert one_civ.absorbers.fit_bayesian([uuid], nsteps="x") == 0
    assert one_civ.absorbers.fit_bayesian([uuid], nsteps="0") == 0


def test_detect_anchor_without_sibling(one_civ):
    assert one_civ.absorbers.detect_anchor("missing", "CIV_1548", "CIV_1550") == 0.0
