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
