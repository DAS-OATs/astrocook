"""Check of the Voigt-fit uncertainties against the reduced chi-squared
(phase 1, item 9 of ``ROADMAP.md``).

Monte Carlo on mock spectra (``astrocook.mocks``, independent of the fitting
code): one CIV component (z = 2.97, logN = 13.5, b = 12 km/s), R = 45000,
2.5 km/s pixels, S/N = 20, fitted by ``VoigtFitterV2`` from a nearby guess
in many noise realizations. The "pull" (fitted - true) / error must have a
standard deviation of 1 if the errors are right.

The errors are sqrt(diag((J^T J)^-1)) with residuals weighted by 1/dy,
scaled by sqrt(reduced chi2) as in V1 (lmfit, ``scale_covar=True``):
findings register F29, decided with Guido. Unscaled, they were right only
when the error spectrum was right, and wrong by the same factor as the
error spectrum otherwise; scaled, they are right in all three cases.
"""
import dataclasses

import numpy as np
import pytest

from astrocook.core.structures import ComponentDataV2, DataColumnV2, SystemListDataV2
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.mocks import MockComponent, MockTruth, generate

TRUE = dict(z=2.97, logN=13.5, b=12.0)
N_REALIZATIONS = 100   # pull std known to ~7%
PULL_TOL = 0.2


def _fit(seed, dy_factor):
    """Fit one noise realization; the error spectrum is multiplied by
    ``dy_factor`` (1: correct; 2: overestimated; 0.5: underestimated)."""
    truth = MockTruth(components=[MockComponent("CIV", **TRUE)],
                      lam_min_nm=613.5, lam_max_nm=616.5, dv_kms=2.5,
                      resol=45000.0, snr=20.0, z_em=3.1, seed=seed)
    spec = generate(truth)
    data = spec._data
    spec = type(spec)(dataclasses.replace(
        data, dy=DataColumnV2(data.dy.values * dy_factor, data.dy.unit)))
    guess = ComponentDataV2(id=1, z=TRUE["z"] + 2e-5, dz=None, logN=13.3,
                            dlogN=None, b=10.0, db=None, series="CIV")
    systs = SystemListV2(SystemListDataV2(components=[guess]))
    fitted, _, result = VoigtFitterV2(spec, systs).fit(max_nfev=2000,
                                                       z_window_kms=20.0)
    comp = fitted.components[0]
    red_chi2 = np.sum(result.fun ** 2) / (len(result.fun) - 3)
    return [comp.z, comp.dz, comp.logN, comp.dlogN, comp.b, comp.db, red_chi2]


@pytest.fixture(scope="module", params=[1.0, 2.0, 0.5], ids=["dy_ok", "dy_x2", "dy_x0.5"])
def monte_carlo(request):
    rows = np.array([_fit(seed, request.param) for seed in range(N_REALIZATIONS)])
    return request.param, rows


def _pulls(rows):
    z, dz, logN, dlogN, b, db, _ = rows.T
    return {"z": (z - TRUE["z"]) / dz,
            "logN": (logN - TRUE["logN"]) / dlogN,
            "b": (b - TRUE["b"]) / db}


def test_reduced_chi2_follows_the_error_spectrum(monte_carlo):
    dy_factor, rows = monte_carlo
    assert np.median(rows[:, -1]) == pytest.approx(dy_factor ** -2, rel=0.15)


def test_errors_are_right_in_all_cases(monte_carlo):
    # Regression test for F29: with the errors scaled by sqrt(reduced chi2)
    # the pulls have std 1 also when the error spectrum is off by 2 or 1/2
    _, rows = monte_carlo
    for name, pull in _pulls(rows).items():
        assert np.std(pull) == pytest.approx(1.0, rel=PULL_TOL), name
        assert abs(np.mean(pull)) < 0.35, name
