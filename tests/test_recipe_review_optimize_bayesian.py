"""Review of ``optimize_system`` and ``fit_bayesian`` (phase 6-quater, PR 10-bis).

The behaviour kept and the behaviour changed, each test saying what it
was. The values of the valid results were taken on master before the
change (``tests/snapshots/optimize_bayesian_*.npz``), with the tolerances
of PR 10 (1e-4 relative for the values and the model, 1e-2 for the
errors): the numbers of ``optimize_system`` and of the MCMC do not move.
emcee draws from the global generator of NumPy, so the MCMC cases fix it
with ``np.random.seed`` (the generator a recipe without ``seed`` uses, as
on master). Nested sampling crashed on master (F139): its numbers are new,
and its posterior is the weighted one (F140).
"""
import warnings

import numpy as np
import pytest

from test_recipe_review_absorbers_manual import _fingerprint, _ids, _values


def _quiet(run):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run()


def _check(snapshot, case, session):
    values = _values(session)
    errors = {"errors": values.pop("errors")}
    results = session.systs._data.meta.get("bayesian_results", {})
    if results:
        # The same group for each of its components: one is enough
        first = results[next(iter(results))]
        values["samples_shape"] = np.array(first["samples"].shape, dtype=float)
        values["samples_mean"] = np.mean(first["samples"], axis=0)
        values["labels"] = np.array(first["labels"], dtype=str)
    snapshot.check(f"optimize_bayesian_{case}", values, rtol=1e-4, atol=1e-12)
    snapshot.check(f"optimize_bayesian_{case}_errors", errors, rtol=1e-2, atol=1e-12)


def _seeded(seed, run):
    np.random.seed(seed)
    return run()


#: Valid calls, on the sessions of test_recipe_inputs_absorbers (the mock
#: of a CIV doublet and two Ly-a, fitted by components_auto)
KEPT = {
    "optimize_system": lambda S, u, lya: S["fitted"].absorbers.optimize_system(u),
    "optimize_system_options": lambda S, u, lya: S["fitted"].absorbers.optimize_system(
        lya[0], max_components=2, threshold_sigma=1.5, aic_penalty=0.0,
        z_window_kms=60.0, min_dv=5.0, group_depth=0, patience=1),
    "optimize_system_strings": lambda S, u, lya: S["fitted"].absorbers.optimize_system(
        u, "3", "2.0", "1.0", "80.0"),
    "optimize_system_civ": lambda S, u, lya: S["fitted"].absorbers.optimize_system(
        u, threshold_sigma=1.0, aic_penalty=0.0, min_dv=3.0, patience=3),
    "mcmc": lambda S, u, lya: _seeded(0, lambda: S["fitted"].absorbers.fit_bayesian(
        [u], nsteps=100)),
    "mcmc_defaults": lambda S, u, lya: _seeded(1, lambda: S["fitted"].absorbers.fit_bayesian(
        [u])),
    "mcmc_group": lambda S, u, lya: _seeded(2, lambda: S["fitted"].absorbers.fit_bayesian(
        lya, nsteps=60, group_depth=0)),
}


@pytest.mark.parametrize("case", sorted(KEPT))
def test_valid_results_unchanged(case, snapshot):
    sessions, u, lya = _ids()
    before = _fingerprint(sessions["fitted"])
    result = _quiet(lambda: KEPT[case](sessions, u, lya))
    _check(snapshot, case, result)
    assert _fingerprint(sessions["fitted"]) == before


def test_max_nfev_default_is_the_old_fixed_value(snapshot):
    # max_nfev was fixed at 1000 inside optimize_hierarchy: the default
    # keeps it
    sessions, u, lya = _ids()
    result = _quiet(lambda: sessions["fitted"].absorbers.optimize_system(
        lya[0], max_components=2, threshold_sigma=1.5, aic_penalty=0.0,
        z_window_kms=60.0, min_dv=5.0, group_depth=0, patience=1, max_nfev=1000))
    _check(snapshot, "optimize_system_options", result)


def test_a_fit_that_reaches_max_nfev_adds_nothing():
    # F149 (registered, not changed): a fit that does not converge within
    # max_nfev counts as failed, so its candidate is rejected; here the
    # starting fit fails too, and nothing can be added
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    result = _quiet(lambda: fitted.absorbers.optimize_system(
        lya[0], max_components=2, threshold_sigma=1.5, aic_penalty=0.0,
        z_window_kms=60.0, min_dv=5.0, group_depth=0, patience=1, max_nfev=1))
    assert len(result.systs.components) == len(fitted.systs.components)


# --- Refused before the recipe runs ---------------------------------------

UNKNOWN = "00000000-0000-0000-0000-000000000000"

#: What each did on master in the comment
REFUSED = {
    # A session with the model recomputed and nothing else (F142)
    "optimize_unknown": (lambda S, u, lya: S["fitted"].absorbers.optimize_system(UNKNOWN),
                         "uuid: there is no component"),
    "optimize_empty": (lambda S, u, lya: S["fitted"].absorbers.optimize_system(""),
                       "uuid: there is no component"),
    # RecipeFailed: "cond is not defined on empty arrays" (F142)
    "bayes_unknown": (lambda S, u, lya: S["fitted"].absorbers.fit_bayesian([u, UNKNOWN]),
                      "uuids: there is no component"),
    # A RecipeFailed: "No components are selected"
    "bayes_empty": (lambda S, u, lya: S["fitted"].absorbers.fit_bayesian([]),
                    "uuids: give at least one"),
    # nlive 1: IndexError; 3: RecipeFailed; 4 and 7: no end (3 parameters)
    "nlive_too_few": (lambda S, u, lya: S["fitted"].absorbers.fit_bayesian(
        [u], engine="nested", nlive=29, group_depth=0),
        "nlive: at least 10 live points per free parameter are needed, 30 for the 3"),
    "nlive_below_the_minimum": (lambda S, u, lya: S["fitted"].absorbers.fit_bayesian(
        [u], engine="nested", nlive=4), "nlive"),
}


@pytest.mark.parametrize("case", sorted(REFUSED))
def test_refused_before_the_recipe_runs(case):
    from astrocook.recipes.errors import InvalidParameter
    sessions, u, lya = _ids()
    run, message = REFUSED[case]
    before = _fingerprint(sessions["fitted"])
    with pytest.raises(InvalidParameter, match=message):
        _quiet(lambda: run(sessions, u, lya))
    assert _fingerprint(sessions["fitted"]) == before


def test_uuids_is_required():
    sessions, u, lya = _ids()
    with pytest.raises(TypeError, match="missing a required argument"):
        sessions["fitted"].absorbers.fit_bayesian()


# --- The samplers ---------------------------------------------------------

def _many(sessions, n=6):
    """The fitted session with n Ly-a components blended near the first
    one: 3 n free parameters."""
    fitted = sessions["fitted"]
    z0 = [c.z for c in fitted.systs.components if c.series == "Ly_a"][0]
    systs = fitted.systs
    for k in range(1, n):
        systs = systs.add_component(series="Ly_a", z=z0 + k * 1e-4, logN=12.5, b=15.0)
    ly = [c.uuid for c in systs.components if c.series == "Ly_a"]
    return fitted.with_new_system_list(systs), ly


def test_mcmc_with_more_than_16_free_parameters():
    # F141: 32 walkers, and emcee needs twice the free parameters: a
    # RuntimeError (not converted by the recipes) with more than 16. The
    # walkers are now max(32, 2 * ndim): 36 for the 18 here
    sessions, u, lya = _ids()
    session, ly = _many(sessions)
    out = _quiet(lambda: session.absorbers.fit_bayesian(
        ly[:6], nsteps=4, group_depth=0, seed=1))
    samples = out.systs._data.meta["bayesian_results"][ly[0]]["samples"]
    assert samples.shape == (36 * 2, 18)


def test_too_few_walkers_are_refused():
    from astrocook.fitting.bayesian_fitter import BayesianVoigtFitter
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    with fitted.systs.fitting_context([u], group_depth=0, spec=fitted.spec):
        b_fitter = BayesianVoigtFitter(fitted.spec, fitted.systs)
        with pytest.raises(ValueError, match="twice the free parameters"):
            b_fitter.run_mcmc(nwalkers=5, nsteps=10, progress=False)


def test_a_sampler_error_is_a_failure_of_the_recipe(monkeypatch):
    from astrocook.fitting.bayesian_fitter import BayesianVoigtFitter
    from astrocook.recipes.errors import RecipeFailed

    def broken(self, **kwargs):
        raise RuntimeError("It is unadvisable to use a red-blue move")

    monkeypatch.setattr(BayesianVoigtFitter, "run_mcmc", broken)
    sessions, u, lya = _ids()
    with pytest.raises(RecipeFailed, match="The sampler failed: It is unadvisable"):
        sessions["fitted"].absorbers.fit_bayesian([u])


# The time of nested sampling varies with the seed, by a factor of 5 on
# the same component (7-45 s here, not with nlive): the seeds of the tests
# are fast ones
@pytest.mark.parametrize("engine, kwargs", [("mcmc", {"nsteps": 40}),
                                            ("nested", {"nlive": 30})])
def test_the_seed_repeats_the_sampling(engine, kwargs):
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    state = np.random.get_state()[1].copy()
    runs = [_quiet(lambda: fitted.absorbers.fit_bayesian(
        [u], engine=engine, group_depth=0, seed=1, **kwargs)) for _ in range(2)]
    assert _fingerprint(runs[0]) == _fingerprint(runs[1])
    # The global generator of NumPy is not touched
    assert np.array_equal(np.random.get_state()[1], state)


def test_nested_sampling_runs_and_logs_its_progress(monkeypatch, caplog):
    # F139: the line of the log on the progress formatted a dict (dynesty
    # 3), a TypeError after 0.5 s; here every iteration is "after 0.5 s"
    import itertools
    import logging
    import time
    clock = itertools.count()
    monkeypatch.setattr(time, "time", lambda: float(next(clock)))
    sessions, u, lya = _ids()
    with caplog.at_level(logging.INFO, logger="astrocook.fitting.bayesian_fitter"):
        out = sessions["fitted"].absorbers.fit_bayesian(
            [u], engine="nested", nlive=30, group_depth=0, seed=1)
    assert "dlogZ" in caplog.text
    comp = [c for c in out.systs.components if c.uuid == u][0]
    assert np.isfinite([comp.z, comp.logN, comp.b, comp.dz, comp.dlogN, comp.db]).all()


def test_nested_posterior_is_weighted():
    # F140: the medians and standard deviations were those of the dead
    # points without their weights: on this component the errors came out
    # hundreds of times those of the MCMC. Now the two samplers agree
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    nested = _quiet(lambda: fitted.absorbers.fit_bayesian(
        [u], engine="nested", nlive=100, group_depth=0, seed=2))
    mcmc = _quiet(lambda: fitted.absorbers.fit_bayesian(
        [u], nsteps=500, group_depth=0, seed=2))
    n = [c for c in nested.systs.components if c.uuid == u][0]
    m = [c for c in mcmc.systs.components if c.uuid == u][0]
    for p in ("z", "logN", "b"):
        err = getattr(m, "d" + p)
        assert abs(getattr(n, p) - getattr(m, p)) < 2 * err, p
        assert 0.5 < getattr(n, "d" + p) / err < 2, p
    # The samples stored for the corner plot have equal weights
    stored = nested.systs._data.meta["bayesian_results"][u]
    assert stored["engine"] == "nested"
    assert np.std(stored["samples"][:, 1]) == pytest.approx(n.dlogN, rel=0.3)


def test_the_gui_passes_only_the_components_and_the_engine():
    # It passed "500", "250", "2": strings that copied the defaults
    from types import SimpleNamespace
    from unittest.mock import Mock
    from astrocook.gui.system_inspector import SystemInspector
    window = Mock()
    fake = SimpleNamespace(main_window=window)
    SystemInspector._fit_bayesian(fake, [SimpleNamespace(uuid="a"), SimpleNamespace(uuid="b")],
                                  "mcmc")
    window._on_recipe_requested.assert_called_once_with(
        "absorbers", "fit_bayesian", {"uuids": ["a", "b"], "engine": "mcmc"}, {})
