"""What the typed recipes of Absorbers do with each kind of input (phase 6-bis, PR 2c).

As test_recipe_inputs.py and test_recipe_inputs_edit.py: for each
value given to a parameter (with the other parameters in base), on a
small mock (a CIV doublet on the red side, two Ly-a lines in the forest,
z_em = 2.7, R = 45000), the value the recipe uses (the recipe must give
the same result as with that typed value: the data, the mask, the
components with their flags, the constraints and the identification),
how it fails (REFUSED, FAILS: ``recipe_outcomes.py``; exceptions since
PR 3c) or ACCEPTED (fit_bayesian,
whose sampling is not seeded: the result is not compared). The cases whose
outcome changed from 2.4.0 say what 2.4.0 did (the battery was run on
master before the code changed). Guido's decisions on the bounds: the
upper bounds of the fits within those of the fitter (max_logN in
(10, 23], max_b in (1, 300]); logN, b and btur of a component within them
(10-23, 1-300, 0-100), its z above -1, its resol positive; score_threshold
in [0, 1]; 0 allowed where it means "none" (merge_dv, min_b, aic_penalty,
min_dv, max_components, patience); counts at least 1 (min_pix_region,
max_nfev, max_group_size, nsteps, nlive), group_depth at least 0 (0 is
the "iterative single" strategy of refit_all, used by metals_auto);
windows and thresholds in sigma positive.
"""
import functools
import logging
import warnings

import numpy as np
import pytest

from astrocook.core.session import SessionV2
from astrocook.mocks import MockComponent, MockTruth, generate

# The outcomes of a failure (phase 6-bis, PR 3): see recipe_outcomes.py
from recipe_outcomes import FAILS, FAILURES, RAISES, REFUSED, check_failure  # noqa: E402,F401
ACCEPTED = "accepted"
NAN = float("nan")
TRUTH = MockTruth(components=[MockComponent("CIV", 1.93, 13.8, 10.0),
                              MockComponent("Ly_a", 2.65, 13.5, 25.0),
                              MockComponent("Ly_a", 2.64, 14.0, 30.0)],
                  lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5, resol=45000.0,
                  snr=50.0, z_em=2.7, seed=1)


class _Token:
    """A value known only when the sessions are built (the UUID of the CIV
    component, the identification)."""
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return self.name


UUID, SERIES_JSON, Z_JSON, LIST_AS_TEXT = (_Token(n) for n in
                                           ("UUID", "SERIES_JSON", "Z_JSON", "LIST_AS_TEXT"))


@functools.lru_cache(maxsize=1)
def _sessions():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = SessionV2("mock", None, spec=generate(TRUTH, truth_columns=False)
                         ).continuum.estimate_auto()
        fitted = base.absorbers.components_auto()
        ident = base.absorbers.identify_lines()
    uuid = [c.uuid for c in fitted.systs.components if c.series == "CIV"][0]
    values = {"UUID": uuid, "[UUID]": [uuid], "LIST_AS_TEXT": repr([uuid]),
              "SERIES_JSON": ident.spec.meta.get("series_map_json"),
              "Z_JSON": ident.spec.meta.get("z_map_json")}
    return {"base": base, "fitted": fitted, "ident": ident}, values


def _resolve(value):
    """The value, with the tokens replaced (also inside a list)."""
    values = _sessions()[1]
    if isinstance(value, _Token):
        return values[value.name]
    if isinstance(value, list) and len(value) == 1 and value[0] is UUID:
        return [values["UUID"]]
    return value


def _run(name, session, param, value, base):
    kwargs = {k: _resolve(v) for k, v in base.items()}
    kwargs[param] = _resolve(value)
    return getattr(_sessions()[0][session].absorbers, name)(**kwargs)


def _fingerprint(result):
    """The data, the mask and identification, the components with their
    flags, the constraints of a result."""
    columns = [np.asarray(result.spec.y.value)]
    for c in ("cont", "model", "abs_mask", "abs_ids"):
        if result.spec.has_aux_column(c):
            columns.append(np.round(np.nan_to_num(np.asarray(result.spec.t[c], dtype=float),
                                                  nan=-9e9), 6))
    comps = [(c.series, round(c.z, 7), round(c.logN, 4), round(c.b, 3), round(c.btur, 3),
              repr(getattr(c, "resol", None)), repr(c.flags))
             for c in (result.systs.components if result.systs else [])]
    constraints = repr(getattr(result.systs, "constraints", None)) if result.systs else ""
    return ([c.tobytes() for c in columns], comps, constraints,
            result.spec.meta.get("series_map_json"), result.spec.meta.get("z_map_json"))


@functools.lru_cache(maxsize=None)
def _reference(name, session, param, expected, base):
    """The result with the typed value (cached: several inputs share it)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = _run(name, session, param, _unfreeze(expected),
                      {k: _unfreeze(v) for k, v in base})
    assert isinstance(result, SessionV2)
    return _fingerprint(result)


def _freeze(value):
    return tuple(value) if isinstance(value, list) else value


def _unfreeze(value):
    return list(value) if isinstance(value, tuple) else value


CASES = [
    ('components_auto', 'base', 'merge_dv', '2.0', {}, 2.0),
    ('components_auto', 'base', 'merge_dv', 2.0, {}, 2.0),
    ('components_auto', 'base', 'merge_dv', 2, {}, 2.0),
    ('components_auto', 'base', 'merge_dv', ' 2 ', {}, 2.0),
    ('components_auto', 'base', 'merge_dv', '1e1', {}, 10.0),
    ('components_auto', 'base', 'merge_dv', 'abc', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('components_auto', 'base', 'merge_dv', '', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('components_auto', 'base', 'merge_dv', 'None', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('components_auto', 'base', 'merge_dv', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('components_auto', 'base', 'merge_dv', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('components_auto', 'base', 'merge_dv', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('components_auto', 'base', 'merge_dv', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('components_auto', 'base', 'merge_dv', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'merge_dv', '0', {}, 0.0),
    ('components_auto', 'base', 'merge_dv', '30', {}, 30.0),
    ('components_auto', 'base', 'z_window_kms', '2.0', {}, 2.0),
    ('components_auto', 'base', 'z_window_kms', 2.0, {}, 2.0),
    ('components_auto', 'base', 'z_window_kms', 2, {}, 2.0),
    ('components_auto', 'base', 'z_window_kms', ' 2 ', {}, 2.0),
    ('components_auto', 'base', 'z_window_kms', '1e1', {}, 10.0),
    ('components_auto', 'base', 'z_window_kms', 'abc', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'z_window_kms', '', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'z_window_kms', 'None', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'z_window_kms', None, {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'z_window_kms', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('components_auto', 'base', 'z_window_kms', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('components_auto', 'base', 'z_window_kms', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('components_auto', 'base', 'z_window_kms', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'z_window_kms', '0', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'z_window_kms', '30', {}, 30.0),
    ('components_auto', 'base', 'max_b', '2.0', {}, 2.0),
    ('components_auto', 'base', 'max_b', 2.0, {}, 2.0),
    ('components_auto', 'base', 'max_b', 2, {}, 2.0),
    ('components_auto', 'base', 'max_b', ' 2 ', {}, 2.0),
    ('components_auto', 'base', 'max_b', '1e1', {}, 10.0),
    ('components_auto', 'base', 'max_b', 'abc', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_b', '', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_b', 'None', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_b', None, {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_b', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('components_auto', 'base', 'max_b', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('components_auto', 'base', 'max_b', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('components_auto', 'base', 'max_b', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_b', '0', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_b', '30', {}, 30.0),
    ('components_auto', 'base', 'max_logN', '2.0', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_logN', 2.0, {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_logN', 2, {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_logN', ' 2 ', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_logN', '1e1', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_logN', 'abc', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_logN', '', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_logN', 'None', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_logN', None, {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_logN', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('components_auto', 'base', 'max_logN', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('components_auto', 'base', 'max_logN', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('components_auto', 'base', 'max_logN', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_logN', '0', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'max_logN', '30', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'score_threshold', '0.3', {}, 0.3),
    ('components_auto', 'base', 'score_threshold', 0.3, {}, 0.3),
    ('components_auto', 'base', 'score_threshold', '0', {}, 0.0),
    ('components_auto', 'base', 'score_threshold', '1', {}, 1.0),
    ('components_auto', 'base', 'score_threshold', '1.5', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'score_threshold', '-0.1', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'score_threshold', 'abc', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('components_auto', 'base', 'score_threshold', '', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('components_auto', 'base', 'score_threshold', 'None', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('components_auto', 'base', 'score_threshold', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('components_auto', 'base', 'score_threshold', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('components_auto', 'base', 'score_threshold', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('components_auto', 'base', 'min_pix_region', '2', {}, 2),
    ('components_auto', 'base', 'min_pix_region', 2, {}, 2),
    ('components_auto', 'base', 'min_pix_region', '2.0', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('components_auto', 'base', 'min_pix_region', 2.0, {}, 2),
    ('components_auto', 'base', 'min_pix_region', 2.5, {}, REFUSED),  # 2.4.0: accepted, truncated
    ('components_auto', 'base', 'min_pix_region', '', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('components_auto', 'base', 'min_pix_region', 'abc', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('components_auto', 'base', 'min_pix_region', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('components_auto', 'base', 'min_pix_region', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('components_auto', 'base', 'min_pix_region', 'None', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('components_auto', 'base', 'min_pix_region', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'min_pix_region', '0', {}, REFUSED),  # 2.4.0: accepted
    ('components_auto', 'base', 'min_pix_region', '1', {}, 1),
    ('components_auto', 'base', 'min_pix_region', '5', {}, 5),
    ('components_auto', 'base', 'region_limit', 'None', {}, 'None'),
    ('components_auto', 'base', 'region_limit', 'red_side', {}, 'red_side'),
    ('components_auto', 'base', 'region_limit', 'lya_forest', {}, 'lya_forest'),
    # Refused since phase 6-bis (PR 6c): a list of choices (they failed in the recipe)
    ('components_auto', 'base', 'region_limit', 'nope', {}, REFUSED),
    ('components_auto', 'base', 'region_limit', '', {}, REFUSED),
    ('components_auto', 'base', 'region_limit', None, {}, REFUSED),
    ('components_auto', 'base', 'multiplets', 'CIV', {}, 'CIV'),
    ('components_auto', 'base', 'multiplets', 'CIV,Ly_ab', {}, 'CIV,Ly_ab'),
    ('components_auto', 'base', 'multiplets', 'nope', {}, 'nope'),
    ('components_auto', 'base', 'multiplets', '', {}, ''),
    ('components_auto', 'base', 'multiplets', None, {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('doublets_auto', 'base', 'z_window_kms', '2.0', {}, 2.0),
    ('doublets_auto', 'base', 'z_window_kms', 2.0, {}, 2.0),
    ('doublets_auto', 'base', 'z_window_kms', 2, {}, 2.0),
    ('doublets_auto', 'base', 'z_window_kms', ' 2 ', {}, 2.0),
    ('doublets_auto', 'base', 'z_window_kms', '1e1', {}, 10.0),
    ('doublets_auto', 'base', 'z_window_kms', 'abc', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'z_window_kms', '', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'z_window_kms', 'None', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'z_window_kms', None, {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'z_window_kms', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('doublets_auto', 'base', 'z_window_kms', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('doublets_auto', 'base', 'z_window_kms', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('doublets_auto', 'base', 'z_window_kms', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'z_window_kms', '0', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'z_window_kms', '30', {}, 30.0),
    ('doublets_auto', 'base', 'max_b', '2.0', {}, 2.0),
    ('doublets_auto', 'base', 'max_b', 2.0, {}, 2.0),
    ('doublets_auto', 'base', 'max_b', 2, {}, 2.0),
    ('doublets_auto', 'base', 'max_b', ' 2 ', {}, 2.0),
    ('doublets_auto', 'base', 'max_b', '1e1', {}, 10.0),
    ('doublets_auto', 'base', 'max_b', 'abc', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_b', '', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_b', 'None', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_b', None, {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_b', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('doublets_auto', 'base', 'max_b', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('doublets_auto', 'base', 'max_b', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('doublets_auto', 'base', 'max_b', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_b', '0', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_b', '30', {}, 30.0),
    ('doublets_auto', 'base', 'max_logN', '2.0', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_logN', 2.0, {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_logN', 2, {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_logN', ' 2 ', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_logN', '1e1', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_logN', 'abc', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_logN', '', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_logN', 'None', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_logN', None, {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_logN', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('doublets_auto', 'base', 'max_logN', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('doublets_auto', 'base', 'max_logN', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('doublets_auto', 'base', 'max_logN', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_logN', '0', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'max_logN', '30', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'score_threshold', '0.3', {}, 0.3),
    ('doublets_auto', 'base', 'score_threshold', 0.3, {}, 0.3),
    ('doublets_auto', 'base', 'score_threshold', '0', {}, 0.0),
    ('doublets_auto', 'base', 'score_threshold', '1', {}, 1.0),
    ('doublets_auto', 'base', 'score_threshold', '1.5', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'score_threshold', '-0.1', {}, REFUSED),  # 2.4.0: accepted
    ('doublets_auto', 'base', 'score_threshold', 'abc', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('doublets_auto', 'base', 'score_threshold', '', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('doublets_auto', 'base', 'score_threshold', 'None', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('doublets_auto', 'base', 'score_threshold', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('doublets_auto', 'base', 'score_threshold', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('doublets_auto', 'base', 'score_threshold', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('lya_auto', 'base', 'min_b', '2.0', {}, 2.0),
    ('lya_auto', 'base', 'min_b', 2.0, {}, 2.0),
    ('lya_auto', 'base', 'min_b', 2, {}, 2.0),
    ('lya_auto', 'base', 'min_b', ' 2 ', {}, 2.0),
    ('lya_auto', 'base', 'min_b', '1e1', {}, 10.0),
    ('lya_auto', 'base', 'min_b', 'abc', {}, REFUSED),
    ('lya_auto', 'base', 'min_b', '', {}, REFUSED),
    ('lya_auto', 'base', 'min_b', 'None', {}, REFUSED),
    ('lya_auto', 'base', 'min_b', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('lya_auto', 'base', 'min_b', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('lya_auto', 'base', 'min_b', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('lya_auto', 'base', 'min_b', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('lya_auto', 'base', 'min_b', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'min_b', '0', {}, 0.0),
    ('lya_auto', 'base', 'min_b', '30', {}, 30.0),
    ('lya_auto', 'base', 'z_window_kms', '2.0', {}, 2.0),
    ('lya_auto', 'base', 'z_window_kms', 2.0, {}, 2.0),
    ('lya_auto', 'base', 'z_window_kms', 2, {}, 2.0),
    ('lya_auto', 'base', 'z_window_kms', ' 2 ', {}, 2.0),
    ('lya_auto', 'base', 'z_window_kms', '1e1', {}, 10.0),
    ('lya_auto', 'base', 'z_window_kms', 'abc', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'z_window_kms', '', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'z_window_kms', 'None', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'z_window_kms', None, {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'z_window_kms', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('lya_auto', 'base', 'z_window_kms', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('lya_auto', 'base', 'z_window_kms', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('lya_auto', 'base', 'z_window_kms', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'z_window_kms', '0', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'z_window_kms', '30', {}, 30.0),
    ('lya_auto', 'base', 'max_b', '2.0', {}, 2.0),
    ('lya_auto', 'base', 'max_b', 2.0, {}, 2.0),
    ('lya_auto', 'base', 'max_b', 2, {}, 2.0),
    ('lya_auto', 'base', 'max_b', ' 2 ', {}, 2.0),
    ('lya_auto', 'base', 'max_b', '1e1', {}, 10.0),
    ('lya_auto', 'base', 'max_b', 'abc', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_b', '', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_b', 'None', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_b', None, {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_b', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('lya_auto', 'base', 'max_b', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('lya_auto', 'base', 'max_b', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('lya_auto', 'base', 'max_b', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_b', '0', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_b', '30', {}, 30.0),
    ('lya_auto', 'base', 'max_logN', '2.0', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_logN', 2.0, {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_logN', 2, {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_logN', ' 2 ', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_logN', '1e1', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_logN', 'abc', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_logN', '', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_logN', 'None', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_logN', None, {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_logN', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('lya_auto', 'base', 'max_logN', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('lya_auto', 'base', 'max_logN', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('lya_auto', 'base', 'max_logN', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_logN', '0', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'max_logN', '30', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'score_threshold', '0.3', {}, 0.3),
    ('lya_auto', 'base', 'score_threshold', 0.3, {}, 0.3),
    ('lya_auto', 'base', 'score_threshold', '0', {}, 0.0),
    ('lya_auto', 'base', 'score_threshold', '1', {}, 1.0),
    ('lya_auto', 'base', 'score_threshold', '1.5', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'score_threshold', '-0.1', {}, REFUSED),  # 2.4.0: accepted
    ('lya_auto', 'base', 'score_threshold', 'abc', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('lya_auto', 'base', 'score_threshold', '', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('lya_auto', 'base', 'score_threshold', 'None', {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('lya_auto', 'base', 'score_threshold', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('lya_auto', 'base', 'score_threshold', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('lya_auto', 'base', 'score_threshold', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('forest_auto', 'base', 'min_b', '2.0', {}, 2.0),
    ('forest_auto', 'base', 'min_b', 2.0, {}, 2.0),
    ('forest_auto', 'base', 'min_b', 2, {}, 2.0),
    ('forest_auto', 'base', 'min_b', ' 2 ', {}, 2.0),
    ('forest_auto', 'base', 'min_b', '1e1', {}, 10.0),
    ('forest_auto', 'base', 'min_b', 'abc', {}, REFUSED),
    ('forest_auto', 'base', 'min_b', '', {}, REFUSED),
    ('forest_auto', 'base', 'min_b', 'None', {}, REFUSED),
    ('forest_auto', 'base', 'min_b', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('forest_auto', 'base', 'min_b', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('forest_auto', 'base', 'min_b', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('forest_auto', 'base', 'min_b', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('forest_auto', 'base', 'min_b', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('forest_auto', 'base', 'min_b', '0', {}, 0.0),
    ('forest_auto', 'base', 'min_b', '30', {}, 30.0),
    ('forest_auto', 'base', 'z_window_kms', '2.0', {}, 2.0),
    ('forest_auto', 'base', 'z_window_kms', 2.0, {}, 2.0),
    ('forest_auto', 'base', 'z_window_kms', 2, {}, 2.0),
    ('forest_auto', 'base', 'z_window_kms', ' 2 ', {}, 2.0),
    ('forest_auto', 'base', 'z_window_kms', '1e1', {}, 10.0),
    ('forest_auto', 'base', 'z_window_kms', 'abc', {}, REFUSED),
    ('forest_auto', 'base', 'z_window_kms', '', {}, REFUSED),
    ('forest_auto', 'base', 'z_window_kms', 'None', {}, REFUSED),
    ('forest_auto', 'base', 'z_window_kms', None, {}, REFUSED),
    ('forest_auto', 'base', 'z_window_kms', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('forest_auto', 'base', 'z_window_kms', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('forest_auto', 'base', 'z_window_kms', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('forest_auto', 'base', 'z_window_kms', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('forest_auto', 'base', 'z_window_kms', '0', {}, REFUSED),  # 2.4.0: accepted
    ('forest_auto', 'base', 'z_window_kms', '30', {}, 30.0),
    ('forest_auto', 'base', 'max_b', '2.0', {}, 2.0),
    ('forest_auto', 'base', 'max_b', 2.0, {}, 2.0),
    ('forest_auto', 'base', 'max_b', 2, {}, 2.0),
    ('forest_auto', 'base', 'max_b', ' 2 ', {}, 2.0),
    ('forest_auto', 'base', 'max_b', '1e1', {}, 10.0),
    ('forest_auto', 'base', 'max_b', 'abc', {}, REFUSED),
    ('forest_auto', 'base', 'max_b', '', {}, REFUSED),
    ('forest_auto', 'base', 'max_b', 'None', {}, REFUSED),
    ('forest_auto', 'base', 'max_b', None, {}, REFUSED),
    ('forest_auto', 'base', 'max_b', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('forest_auto', 'base', 'max_b', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('forest_auto', 'base', 'max_b', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('forest_auto', 'base', 'max_b', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('forest_auto', 'base', 'max_b', '0', {}, REFUSED),  # 2.4.0: accepted
    ('forest_auto', 'base', 'max_b', '30', {}, 30.0),
    ('metals_auto', 'fitted', 'aic_penalty', '2.0', {}, 2.0),
    ('metals_auto', 'fitted', 'aic_penalty', 2.0, {}, 2.0),
    ('metals_auto', 'fitted', 'aic_penalty', 2, {}, 2.0),
    ('metals_auto', 'fitted', 'aic_penalty', ' 2 ', {}, 2.0),
    ('metals_auto', 'fitted', 'aic_penalty', '1e1', {}, 10.0),
    ('metals_auto', 'fitted', 'aic_penalty', 'abc', {}, REFUSED),
    ('metals_auto', 'fitted', 'aic_penalty', '', {}, REFUSED),
    ('metals_auto', 'fitted', 'aic_penalty', 'None', {}, REFUSED),
    ('metals_auto', 'fitted', 'aic_penalty', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('metals_auto', 'fitted', 'aic_penalty', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('metals_auto', 'fitted', 'aic_penalty', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('metals_auto', 'fitted', 'aic_penalty', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('metals_auto', 'fitted', 'aic_penalty', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('metals_auto', 'fitted', 'aic_penalty', '0', {}, 0.0),
    ('metals_auto', 'fitted', 'aic_penalty', '30', {}, 30.0),
    ('metals_auto', 'fitted', 'z_window_kms', '2.0', {}, 2.0),
    ('metals_auto', 'fitted', 'z_window_kms', 2.0, {}, 2.0),
    ('metals_auto', 'fitted', 'z_window_kms', 2, {}, 2.0),
    ('metals_auto', 'fitted', 'z_window_kms', ' 2 ', {}, 2.0),
    ('metals_auto', 'fitted', 'z_window_kms', '1e1', {}, 10.0),
    ('metals_auto', 'fitted', 'z_window_kms', 'abc', {}, REFUSED),
    ('metals_auto', 'fitted', 'z_window_kms', '', {}, REFUSED),
    ('metals_auto', 'fitted', 'z_window_kms', 'None', {}, REFUSED),
    ('metals_auto', 'fitted', 'z_window_kms', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('metals_auto', 'fitted', 'z_window_kms', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('metals_auto', 'fitted', 'z_window_kms', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('metals_auto', 'fitted', 'z_window_kms', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('metals_auto', 'fitted', 'z_window_kms', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('metals_auto', 'fitted', 'z_window_kms', '0', {}, REFUSED),  # 2.4.0: accepted
    ('metals_auto', 'fitted', 'z_window_kms', '30', {}, 30.0),
    ('metals_auto', 'fitted', 'max_b', '2.0', {}, 2.0),
    ('metals_auto', 'fitted', 'max_b', 2.0, {}, 2.0),
    ('metals_auto', 'fitted', 'max_b', 2, {}, 2.0),
    ('metals_auto', 'fitted', 'max_b', ' 2 ', {}, 2.0),
    ('metals_auto', 'fitted', 'max_b', '1e1', {}, 10.0),
    ('metals_auto', 'fitted', 'max_b', 'abc', {}, REFUSED),
    ('metals_auto', 'fitted', 'max_b', '', {}, REFUSED),
    ('metals_auto', 'fitted', 'max_b', 'None', {}, REFUSED),
    ('metals_auto', 'fitted', 'max_b', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('metals_auto', 'fitted', 'max_b', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('metals_auto', 'fitted', 'max_b', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('metals_auto', 'fitted', 'max_b', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('metals_auto', 'fitted', 'max_b', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('metals_auto', 'fitted', 'max_b', '0', {}, REFUSED),  # 2.4.0: accepted
    ('metals_auto', 'fitted', 'max_b', '30', {}, 30.0),
    ('metals_auto', 'fitted', 'ions', 'NV', {}, 'NV'),
    ('metals_auto', 'fitted', 'ions', 'SiIV,NV', {}, 'SiIV,NV'),
    ('metals_auto', 'fitted', 'ions', 'nope', {}, 'nope'),
    ('metals_auto', 'fitted', 'ions', '', {}, ''),
    ('metals_auto', 'fitted', 'ions', None, {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('identify_lines', 'base', 'min_pix_region', '2', {}, 2),
    ('identify_lines', 'base', 'min_pix_region', 2, {}, 2),
    ('identify_lines', 'base', 'min_pix_region', '2.0', {}, REFUSED),
    ('identify_lines', 'base', 'min_pix_region', 2.0, {}, 2),
    ('identify_lines', 'base', 'min_pix_region', 2.5, {}, REFUSED),  # 2.4.0: accepted, truncated
    ('identify_lines', 'base', 'min_pix_region', '', {}, REFUSED),
    ('identify_lines', 'base', 'min_pix_region', 'abc', {}, REFUSED),
    ('identify_lines', 'base', 'min_pix_region', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('identify_lines', 'base', 'min_pix_region', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('identify_lines', 'base', 'min_pix_region', 'None', {}, REFUSED),
    ('identify_lines', 'base', 'min_pix_region', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('identify_lines', 'base', 'min_pix_region', '0', {}, REFUSED),  # 2.4.0: accepted
    ('identify_lines', 'base', 'min_pix_region', '1', {}, 1),
    ('identify_lines', 'base', 'min_pix_region', '5', {}, 5),
    ('identify_lines', 'base', 'merge_dv', '2.0', {}, 2.0),
    ('identify_lines', 'base', 'merge_dv', 2.0, {}, 2.0),
    ('identify_lines', 'base', 'merge_dv', 2, {}, 2.0),
    ('identify_lines', 'base', 'merge_dv', ' 2 ', {}, 2.0),
    ('identify_lines', 'base', 'merge_dv', '1e1', {}, 10.0),
    ('identify_lines', 'base', 'merge_dv', 'abc', {}, REFUSED),
    ('identify_lines', 'base', 'merge_dv', '', {}, REFUSED),
    ('identify_lines', 'base', 'merge_dv', 'None', {}, REFUSED),
    ('identify_lines', 'base', 'merge_dv', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('identify_lines', 'base', 'merge_dv', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('identify_lines', 'base', 'merge_dv', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('identify_lines', 'base', 'merge_dv', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('identify_lines', 'base', 'merge_dv', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('identify_lines', 'base', 'merge_dv', '0', {}, 0.0),
    ('identify_lines', 'base', 'merge_dv', '30', {}, 30.0),
    ('identify_lines', 'base', 'score_threshold', '0.3', {}, 0.3),
    ('identify_lines', 'base', 'score_threshold', 0.3, {}, 0.3),
    ('identify_lines', 'base', 'score_threshold', '0', {}, 0.0),
    ('identify_lines', 'base', 'score_threshold', '1', {}, 1.0),
    ('identify_lines', 'base', 'score_threshold', '1.5', {}, REFUSED),  # 2.4.0: accepted
    ('identify_lines', 'base', 'score_threshold', '-0.1', {}, REFUSED),  # 2.4.0: accepted
    ('identify_lines', 'base', 'score_threshold', 'abc', {}, REFUSED),
    ('identify_lines', 'base', 'score_threshold', '', {}, REFUSED),
    ('identify_lines', 'base', 'score_threshold', 'None', {}, REFUSED),
    ('identify_lines', 'base', 'score_threshold', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('identify_lines', 'base', 'score_threshold', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('identify_lines', 'base', 'score_threshold', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('identify_lines', 'base', 'bypass_scoring', 'True', {}, True),
    ('identify_lines', 'base', 'bypass_scoring', 'False', {}, False),
    ('identify_lines', 'base', 'bypass_scoring', 'true', {}, True),
    ('identify_lines', 'base', 'bypass_scoring', 'false', {}, False),
    ('identify_lines', 'base', 'bypass_scoring', '1', {}, True),
    ('identify_lines', 'base', 'bypass_scoring', '0', {}, False),
    ('identify_lines', 'base', 'bypass_scoring', 'yes', {}, REFUSED),  # 2.4.0: accepted as False
    ('identify_lines', 'base', 'bypass_scoring', True, {}, True),  # 2.4.0: crashed (AttributeError)
    ('identify_lines', 'base', 'bypass_scoring', False, {}, False),  # 2.4.0: crashed (AttributeError)
    ('identify_lines', 'base', 'bypass_scoring', 1, {}, True),  # 2.4.0: crashed (AttributeError)
    ('identify_lines', 'base', 'bypass_scoring', 0, {}, False),  # 2.4.0: crashed (AttributeError)
    ('identify_lines', 'base', 'bypass_scoring', '', {}, REFUSED),  # 2.4.0: accepted as False
    ('identify_lines', 'base', 'bypass_scoring', None, {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('identify_lines', 'base', 'bypass_scoring', 'None', {}, REFUSED),  # 2.4.0: accepted as False
    ('identify_lines', 'base', 'debug_rating', 'False', {}, False),
    ('identify_lines', 'base', 'debug_rating', False, {}, False),  # 2.4.0: crashed (AttributeError)
    ('identify_lines', 'base', 'debug_rating', 'false', {}, False),
    ('identify_lines', 'base', 'debug_rating', '0', {}, False),
    ('identify_lines', 'base', 'debug_rating', 0, {}, False),  # 2.4.0: crashed (AttributeError)
    ('identify_lines', 'base', 'debug_rating', 'no', {}, REFUSED),  # 2.4.0: accepted as False
    ('identify_lines', 'base', 'debug_rating', '', {}, REFUSED),  # 2.4.0: accepted as False
    ('identify_lines', 'base', 'debug_rating', None, {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('identify_lines', 'base', 'debug_rating', 'None', {}, REFUSED),  # 2.4.0: accepted as False
    ('populate_from_identification', 'ident', 'series_map_json', SERIES_JSON, {'z_map_json': Z_JSON}, SERIES_JSON),
    ('populate_from_identification', 'ident', 'series_map_json', '', {'z_map_json': Z_JSON}, FAILS),
    ('populate_from_identification', 'ident', 'series_map_json', 'abc', {'z_map_json': Z_JSON}, FAILS),
    ('populate_from_identification', 'ident', 'series_map_json', None, {'z_map_json': Z_JSON}, REFUSED),
    ('populate_from_identification', 'ident', 'series_map_json', '{}', {'z_map_json': Z_JSON}, '{}'),
    ('add_component', 'fitted', 'z', '2.0', {'series': 'CIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_component', 'fitted', 'z', 2.0, {'series': 'CIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_component', 'fitted', 'z', 2, {'series': 'CIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_component', 'fitted', 'z', ' 2 ', {'series': 'CIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_component', 'fitted', 'z', '1e1', {'series': 'CIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_component', 'fitted', 'z', 'abc', {'series': 'CIV'}, REFUSED),
    ('add_component', 'fitted', 'z', '', {'series': 'CIV'}, REFUSED),
    ('add_component', 'fitted', 'z', 'None', {'series': 'CIV'}, REFUSED),
    ('add_component', 'fitted', 'z', None, {'series': 'CIV'}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('add_component', 'fitted', 'z', 'nan', {'series': 'CIV'}, REFUSED),
    ('add_component', 'fitted', 'z', 'inf', {'series': 'CIV'}, REFUSED),
    ('add_component', 'fitted', 'z', True, {'series': 'CIV'}, REFUSED),  # 2.4.0: accepted as 1
    ('add_component', 'fitted', 'z', '-1', {'series': 'CIV'}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'z', '0', {'series': 'CIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_component', 'fitted', 'z', '1.931', {'series': 'CIV'}, 1.931),
    ('add_component', 'fitted', 'logN', '2.0', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'logN', 2.0, {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'logN', 2, {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'logN', ' 2 ', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'logN', '1e1', {'series': 'CIV', 'z': 1.931}, 10.0),
    ('add_component', 'fitted', 'logN', 'abc', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'logN', '', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'logN', 'None', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'logN', None, {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('add_component', 'fitted', 'logN', 'nan', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'logN', 'inf', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted as an infinity
    ('add_component', 'fitted', 'logN', True, {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted as 1
    ('add_component', 'fitted', 'logN', '-1', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'logN', '0', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'logN', '30', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'b', '2.0', {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'b', 2.0, {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'b', 2, {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'b', ' 2 ', {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'b', '1e1', {'series': 'CIV', 'z': 1.931}, 10.0),
    ('add_component', 'fitted', 'b', 'abc', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'b', '', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'b', 'None', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'b', None, {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('add_component', 'fitted', 'b', 'nan', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'b', 'inf', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted as an infinity
    ('add_component', 'fitted', 'b', True, {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted as 1
    ('add_component', 'fitted', 'b', '-1', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'b', '0', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'b', '30', {'series': 'CIV', 'z': 1.931}, 30.0),
    ('add_component', 'fitted', 'btur', '2.0', {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'btur', 2.0, {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'btur', 2, {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'btur', ' 2 ', {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'btur', '1e1', {'series': 'CIV', 'z': 1.931}, 10.0),
    ('add_component', 'fitted', 'btur', 'abc', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'btur', '', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'btur', 'None', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'btur', None, {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('add_component', 'fitted', 'btur', 'nan', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted as NaN
    ('add_component', 'fitted', 'btur', 'inf', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted as an infinity
    ('add_component', 'fitted', 'btur', True, {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted as 1
    ('add_component', 'fitted', 'btur', '-1', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'btur', '0', {'series': 'CIV', 'z': 1.931}, 0.0),
    ('add_component', 'fitted', 'btur', '30', {'series': 'CIV', 'z': 1.931}, 30.0),
    ('add_component', 'fitted', 'z_window_kms', '2.0', {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'z_window_kms', 2.0, {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'z_window_kms', 2, {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'z_window_kms', ' 2 ', {'series': 'CIV', 'z': 1.931}, 2.0),
    ('add_component', 'fitted', 'z_window_kms', '1e1', {'series': 'CIV', 'z': 1.931}, 10.0),
    ('add_component', 'fitted', 'z_window_kms', 'abc', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'z_window_kms', '', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'z_window_kms', 'None', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'z_window_kms', None, {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('add_component', 'fitted', 'z_window_kms', 'nan', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'z_window_kms', 'inf', {'series': 'CIV', 'z': 1.931}, REFUSED),
    ('add_component', 'fitted', 'z_window_kms', True, {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted as 1
    ('add_component', 'fitted', 'z_window_kms', '-1', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'z_window_kms', '0', {'series': 'CIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_component', 'fitted', 'z_window_kms', '30', {'series': 'CIV', 'z': 1.931}, 30.0),
    ('add_component', 'fitted', 'series', 'CIV', {'z': 1.931}, 'CIV'),
    ('add_component', 'fitted', 'series', 'Ly_a', {'z': 1.931}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_component', 'fitted', 'series', 'nope', {'z': 1.931}, FAILS),  # F111: no transition, so no pixel; fitted on the whole spectrum until PR 7-quater
    ('add_component', 'fitted', 'series', '', {'z': 1.931}, FAILS),  # F111: no transition, so no pixel; fitted on the whole spectrum until PR 7-quater
    ('add_component', 'fitted', 'series', None, {'z': 1.931}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('update_component', 'fitted', 'z', '2.0', {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'z', 2.0, {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'z', 2, {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'z', ' 2 ', {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'z', '1e1', {'uuid': UUID}, 10.0),
    ('update_component', 'fitted', 'z', 'abc', {'uuid': UUID}, REFUSED),
    ('update_component', 'fitted', 'z', '', {'uuid': UUID}, None),  # 2.4.0: refused (0)
    ('update_component', 'fitted', 'z', 'None', {'uuid': UUID}, None),
    ('update_component', 'fitted', 'z', None, {'uuid': UUID}, None),  # 2.4.0: refused (0)
    ('update_component', 'fitted', 'z', 'nan', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as NaN
    ('update_component', 'fitted', 'z', 'inf', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as an infinity
    ('update_component', 'fitted', 'z', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('update_component', 'fitted', 'z', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'z', '0', {'uuid': UUID}, 0.0),
    ('update_component', 'fitted', 'z', '1.9301', {'uuid': UUID}, 1.9301),
    ('update_components', 'fitted', 'z', '2.0', {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'z', 2.0, {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'z', 2, {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'z', ' 2 ', {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'z', '1e1', {'uuids': [UUID]}, 10.0),
    ('update_components', 'fitted', 'z', 'abc', {'uuids': [UUID]}, REFUSED),
    ('update_components', 'fitted', 'z', '', {'uuids': [UUID]}, None),  # 2.4.0: refused (0)
    ('update_components', 'fitted', 'z', 'None', {'uuids': [UUID]}, None),
    ('update_components', 'fitted', 'z', None, {'uuids': [UUID]}, None),  # 2.4.0: refused (0)
    ('update_components', 'fitted', 'z', 'nan', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as NaN
    ('update_components', 'fitted', 'z', 'inf', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as an infinity
    ('update_components', 'fitted', 'z', True, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as 1
    ('update_components', 'fitted', 'z', '-1', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'z', '0', {'uuids': [UUID]}, 0.0),
    ('update_components', 'fitted', 'z', '1.9301', {'uuids': [UUID]}, 1.9301),
    ('update_component', 'fitted', 'logN', '2.0', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'logN', 2.0, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'logN', 2, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'logN', ' 2 ', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'logN', '1e1', {'uuid': UUID}, 10.0),
    ('update_component', 'fitted', 'logN', 'abc', {'uuid': UUID}, REFUSED),
    ('update_component', 'fitted', 'logN', '', {'uuid': UUID}, None),  # 2.4.0: refused (0)
    ('update_component', 'fitted', 'logN', 'None', {'uuid': UUID}, None),
    ('update_component', 'fitted', 'logN', None, {'uuid': UUID}, None),  # 2.4.0: refused (0)
    ('update_component', 'fitted', 'logN', 'nan', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as NaN
    ('update_component', 'fitted', 'logN', 'inf', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as an infinity
    ('update_component', 'fitted', 'logN', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('update_component', 'fitted', 'logN', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'logN', '0', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'logN', '14', {'uuid': UUID}, 14.0),
    ('update_components', 'fitted', 'logN', '2.0', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'logN', 2.0, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'logN', 2, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'logN', ' 2 ', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'logN', '1e1', {'uuids': [UUID]}, 10.0),
    ('update_components', 'fitted', 'logN', 'abc', {'uuids': [UUID]}, REFUSED),
    ('update_components', 'fitted', 'logN', '', {'uuids': [UUID]}, None),  # 2.4.0: refused (0)
    ('update_components', 'fitted', 'logN', 'None', {'uuids': [UUID]}, None),
    ('update_components', 'fitted', 'logN', None, {'uuids': [UUID]}, None),  # 2.4.0: refused (0)
    ('update_components', 'fitted', 'logN', 'nan', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as NaN
    ('update_components', 'fitted', 'logN', 'inf', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as an infinity
    ('update_components', 'fitted', 'logN', True, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as 1
    ('update_components', 'fitted', 'logN', '-1', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'logN', '0', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'logN', '14', {'uuids': [UUID]}, 14.0),
    ('update_component', 'fitted', 'b', '2.0', {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'b', 2.0, {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'b', 2, {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'b', ' 2 ', {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'b', '1e1', {'uuid': UUID}, 10.0),
    ('update_component', 'fitted', 'b', 'abc', {'uuid': UUID}, REFUSED),
    ('update_component', 'fitted', 'b', '', {'uuid': UUID}, None),  # 2.4.0: refused (0)
    ('update_component', 'fitted', 'b', 'None', {'uuid': UUID}, None),
    ('update_component', 'fitted', 'b', None, {'uuid': UUID}, None),  # 2.4.0: refused (0)
    ('update_component', 'fitted', 'b', 'nan', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as NaN
    ('update_component', 'fitted', 'b', 'inf', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as an infinity
    ('update_component', 'fitted', 'b', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('update_component', 'fitted', 'b', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'b', '0', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'b', '14', {'uuid': UUID}, 14.0),
    ('update_components', 'fitted', 'b', '2.0', {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'b', 2.0, {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'b', 2, {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'b', ' 2 ', {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'b', '1e1', {'uuids': [UUID]}, 10.0),
    ('update_components', 'fitted', 'b', 'abc', {'uuids': [UUID]}, REFUSED),
    ('update_components', 'fitted', 'b', '', {'uuids': [UUID]}, None),  # 2.4.0: refused (0)
    ('update_components', 'fitted', 'b', 'None', {'uuids': [UUID]}, None),
    ('update_components', 'fitted', 'b', None, {'uuids': [UUID]}, None),  # 2.4.0: refused (0)
    ('update_components', 'fitted', 'b', 'nan', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as NaN
    ('update_components', 'fitted', 'b', 'inf', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as an infinity
    ('update_components', 'fitted', 'b', True, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as 1
    ('update_components', 'fitted', 'b', '-1', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'b', '0', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'b', '14', {'uuids': [UUID]}, 14.0),
    ('update_component', 'fitted', 'btur', '2.0', {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'btur', 2.0, {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'btur', 2, {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'btur', ' 2 ', {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'btur', '1e1', {'uuid': UUID}, 10.0),
    ('update_component', 'fitted', 'btur', 'abc', {'uuid': UUID}, REFUSED),
    ('update_component', 'fitted', 'btur', '', {'uuid': UUID}, None),  # 2.4.0: refused (0)
    ('update_component', 'fitted', 'btur', 'None', {'uuid': UUID}, None),
    ('update_component', 'fitted', 'btur', None, {'uuid': UUID}, None),  # 2.4.0: refused (0)
    ('update_component', 'fitted', 'btur', 'nan', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as NaN
    ('update_component', 'fitted', 'btur', 'inf', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as an infinity
    ('update_component', 'fitted', 'btur', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('update_component', 'fitted', 'btur', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'btur', '0', {'uuid': UUID}, 0.0),
    ('update_component', 'fitted', 'btur', '14', {'uuid': UUID}, 14.0),
    ('update_components', 'fitted', 'btur', '2.0', {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'btur', 2.0, {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'btur', 2, {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'btur', ' 2 ', {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'btur', '1e1', {'uuids': [UUID]}, 10.0),
    ('update_components', 'fitted', 'btur', 'abc', {'uuids': [UUID]}, REFUSED),
    ('update_components', 'fitted', 'btur', '', {'uuids': [UUID]}, None),  # 2.4.0: refused (0)
    ('update_components', 'fitted', 'btur', 'None', {'uuids': [UUID]}, None),
    ('update_components', 'fitted', 'btur', None, {'uuids': [UUID]}, None),  # 2.4.0: refused (0)
    ('update_components', 'fitted', 'btur', 'nan', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as NaN
    ('update_components', 'fitted', 'btur', 'inf', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as an infinity
    ('update_components', 'fitted', 'btur', True, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as 1
    ('update_components', 'fitted', 'btur', '-1', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'btur', '0', {'uuids': [UUID]}, 0.0),
    ('update_components', 'fitted', 'btur', '14', {'uuids': [UUID]}, 14.0),
    ('update_component', 'fitted', 'resol', '2.0', {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'resol', 2.0, {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'resol', 2, {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'resol', ' 2 ', {'uuid': UUID}, 2.0),
    ('update_component', 'fitted', 'resol', '1e1', {'uuid': UUID}, 10.0),
    ('update_component', 'fitted', 'resol', 'abc', {'uuid': UUID}, REFUSED),
    ('update_component', 'fitted', 'resol', '', {'uuid': UUID}, None),  # 2.4.0: refused (0)
    ('update_component', 'fitted', 'resol', 'None', {'uuid': UUID}, None),
    ('update_component', 'fitted', 'resol', None, {'uuid': UUID}, None),  # 2.4.0: refused (0)
    ('update_component', 'fitted', 'resol', 'nan', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as NaN
    ('update_component', 'fitted', 'resol', 'inf', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as an infinity
    ('update_component', 'fitted', 'resol', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('update_component', 'fitted', 'resol', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'resol', '0', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('update_component', 'fitted', 'resol', '30000', {'uuid': UUID}, 30000.0),
    ('update_components', 'fitted', 'resol', '2.0', {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'resol', 2.0, {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'resol', 2, {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'resol', ' 2 ', {'uuids': [UUID]}, 2.0),
    ('update_components', 'fitted', 'resol', '1e1', {'uuids': [UUID]}, 10.0),
    ('update_components', 'fitted', 'resol', 'abc', {'uuids': [UUID]}, REFUSED),
    ('update_components', 'fitted', 'resol', '', {'uuids': [UUID]}, None),  # 2.4.0: refused (0)
    ('update_components', 'fitted', 'resol', 'None', {'uuids': [UUID]}, None),
    ('update_components', 'fitted', 'resol', None, {'uuids': [UUID]}, None),  # 2.4.0: refused (0)
    ('update_components', 'fitted', 'resol', 'nan', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as NaN
    ('update_components', 'fitted', 'resol', 'inf', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as an infinity
    ('update_components', 'fitted', 'resol', True, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as 1
    ('update_components', 'fitted', 'resol', '-1', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'resol', '0', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('update_components', 'fitted', 'resol', '30000', {'uuids': [UUID]}, 30000.0),
    ('update_component', 'fitted', 'series', 'CIV', {'uuid': UUID}, 'CIV'),
    ('update_component', 'fitted', 'series', 'SiIV', {'uuid': UUID}, 'SiIV'),
    ('update_component', 'fitted', 'series', 'nope', {'uuid': UUID}, 'nope'),
    ('update_component', 'fitted', 'series', '', {'uuid': UUID}, ''),
    ('update_component', 'fitted', 'series', None, {'uuid': UUID}, None),  # 2.4.0: accepted, the series renamed 'None'
    ('update_component', 'fitted', 'series', 'None', {'uuid': UUID}, None),
    ('update_component', 'fitted', 'uuid', UUID, {}, UUID),
    ('update_component', 'fitted', 'uuid', 'nope', {}, 'nope'),
    ('update_component', 'fitted', 'uuid', '', {}, ''),
    ('update_component', 'fitted', 'uuid', None, {}, REFUSED),  # 2.4.0: accepted
    ('delete_component', 'fitted', 'uuid', UUID, {}, UUID),
    ('delete_component', 'fitted', 'uuid', 'nope', {}, 'nope'),
    ('delete_component', 'fitted', 'uuid', '', {}, ''),
    ('delete_component', 'fitted', 'uuid', None, {}, None),
    ('delete_component', 'fitted', 'uuids', [UUID], {}, [UUID]),
    ('delete_component', 'fitted', 'uuids', LIST_AS_TEXT, {}, [UUID]),  # 2.4.0: accepted, the text not read as a list
    ('delete_component', 'fitted', 'uuids', [], {}, []),
    ('delete_component', 'fitted', 'uuids', None, {}, None),
    ('delete_component', 'fitted', 'uuids', 'abc', {}, REFUSED),  # 2.4.0: accepted
    ('fit_component', 'fitted', 'max_nfev', '2', {'uuid': UUID}, 2),
    ('fit_component', 'fitted', 'max_nfev', 2, {'uuid': UUID}, 2),
    ('fit_component', 'fitted', 'max_nfev', '2.0', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'max_nfev', 2.0, {'uuid': UUID}, 2),
    ('fit_component', 'fitted', 'max_nfev', 2.5, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted, truncated
    ('fit_component', 'fitted', 'max_nfev', '', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'max_nfev', 'abc', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'max_nfev', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('fit_component', 'fitted', 'max_nfev', None, {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'max_nfev', 'None', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'max_nfev', '-1', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'max_nfev', '0', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'max_nfev', '1', {'uuid': UUID}, 1),
    ('fit_component', 'fitted', 'max_nfev', '3', {'uuid': UUID}, 3),
    ('fit_component', 'fitted', 'group_depth', '2', {'uuid': UUID}, 2),
    ('fit_component', 'fitted', 'group_depth', 2, {'uuid': UUID}, 2),
    ('fit_component', 'fitted', 'group_depth', '2.0', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'group_depth', 2.0, {'uuid': UUID}, 2),
    ('fit_component', 'fitted', 'group_depth', 2.5, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted, truncated
    ('fit_component', 'fitted', 'group_depth', '', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'group_depth', 'abc', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'group_depth', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('fit_component', 'fitted', 'group_depth', None, {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'group_depth', 'None', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'group_depth', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('fit_component', 'fitted', 'group_depth', '0', {'uuid': UUID}, 0),
    ('fit_component', 'fitted', 'group_depth', '1', {'uuid': UUID}, 1),
    ('fit_component', 'fitted', 'group_depth', '3', {'uuid': UUID}, 3),
    ('fit_component', 'fitted', 'z_window_kms', '2.0', {'uuid': UUID}, 2.0),
    ('fit_component', 'fitted', 'z_window_kms', 2.0, {'uuid': UUID}, 2.0),
    ('fit_component', 'fitted', 'z_window_kms', 2, {'uuid': UUID}, 2.0),
    ('fit_component', 'fitted', 'z_window_kms', ' 2 ', {'uuid': UUID}, 2.0),
    ('fit_component', 'fitted', 'z_window_kms', '1e1', {'uuid': UUID}, 10.0),
    ('fit_component', 'fitted', 'z_window_kms', 'abc', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'z_window_kms', '', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'z_window_kms', 'None', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'z_window_kms', None, {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'z_window_kms', 'nan', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'z_window_kms', 'inf', {'uuid': UUID}, REFUSED),
    ('fit_component', 'fitted', 'z_window_kms', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('fit_component', 'fitted', 'z_window_kms', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('fit_component', 'fitted', 'z_window_kms', '0', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('fit_component', 'fitted', 'z_window_kms', '30', {'uuid': UUID}, 30.0),
    ('update_constraint', 'fitted', 'is_free', 'True', {'uuid': UUID, 'param': 'b'}, True),
    ('update_constraint', 'fitted', 'is_free', 'False', {'uuid': UUID, 'param': 'b'}, False),
    ('update_constraint', 'fitted', 'is_free', 'true', {'uuid': UUID, 'param': 'b'}, True),
    ('update_constraint', 'fitted', 'is_free', 'false', {'uuid': UUID, 'param': 'b'}, False),
    ('update_constraint', 'fitted', 'is_free', '1', {'uuid': UUID, 'param': 'b'}, True),
    ('update_constraint', 'fitted', 'is_free', '0', {'uuid': UUID, 'param': 'b'}, False),
    ('update_constraint', 'fitted', 'is_free', 'yes', {'uuid': UUID, 'param': 'b'}, REFUSED),  # 2.4.0: accepted as False
    ('update_constraint', 'fitted', 'is_free', True, {'uuid': UUID, 'param': 'b'}, True),
    ('update_constraint', 'fitted', 'is_free', False, {'uuid': UUID, 'param': 'b'}, False),
    ('update_constraint', 'fitted', 'is_free', 1, {'uuid': UUID, 'param': 'b'}, True),
    ('update_constraint', 'fitted', 'is_free', 0, {'uuid': UUID, 'param': 'b'}, False),
    ('update_constraint', 'fitted', 'is_free', '', {'uuid': UUID, 'param': 'b'}, None),
    ('update_constraint', 'fitted', 'is_free', None, {'uuid': UUID, 'param': 'b'}, None),
    ('update_constraint', 'fitted', 'is_free', 'None', {'uuid': UUID, 'param': 'b'}, None),
    ('update_constraint', 'fitted', 'param', 'b', {'uuid': UUID, 'is_free': False}, 'b'),
    ('update_constraint', 'fitted', 'param', 'z', {'uuid': UUID, 'is_free': False}, 'z'),
    ('update_constraint', 'fitted', 'param', 'logN', {'uuid': UUID, 'is_free': False}, 'logN'),
    # Refused since phase 6-bis (PR 6c): a constraint on no parameter was stored
    ('update_constraint', 'fitted', 'param', 'nope', {'uuid': UUID, 'is_free': False}, REFUSED),
    ('update_constraint', 'fitted', 'param', '', {'uuid': UUID, 'is_free': False}, REFUSED),
    ('update_constraint', 'fitted', 'param', None, {'uuid': UUID, 'is_free': False}, REFUSED),  # 2.4.0: accepted
    ('update_constraints_column', 'fitted', 'is_free', 'True', {'param': 'b'}, True),
    ('update_constraints_column', 'fitted', 'is_free', 'False', {'param': 'b'}, False),
    ('update_constraints_column', 'fitted', 'is_free', 'true', {'param': 'b'}, True),
    ('update_constraints_column', 'fitted', 'is_free', 'false', {'param': 'b'}, False),
    ('update_constraints_column', 'fitted', 'is_free', '1', {'param': 'b'}, True),
    ('update_constraints_column', 'fitted', 'is_free', '0', {'param': 'b'}, False),
    ('update_constraints_column', 'fitted', 'is_free', 'yes', {'param': 'b'}, REFUSED),  # 2.4.0: accepted as False
    ('update_constraints_column', 'fitted', 'is_free', True, {'param': 'b'}, True),
    ('update_constraints_column', 'fitted', 'is_free', False, {'param': 'b'}, False),
    ('update_constraints_column', 'fitted', 'is_free', 1, {'param': 'b'}, True),
    ('update_constraints_column', 'fitted', 'is_free', 0, {'param': 'b'}, False),
    ('update_constraints_column', 'fitted', 'is_free', '', {'param': 'b'}, REFUSED),  # 2.4.0: accepted as False
    ('update_constraints_column', 'fitted', 'is_free', None, {'param': 'b'}, REFUSED),  # 2.4.0: accepted as False
    ('update_constraints_column', 'fitted', 'is_free', 'None', {'param': 'b'}, REFUSED),  # 2.4.0: accepted as False
    ('update_constraints_column', 'fitted', 'param', 'b', {'is_free': False}, 'b'),
    ('update_constraints_column', 'fitted', 'param', 'btur', {'is_free': False}, 'btur'),
    # Refused since phase 6-bis (PR 6c): as update_constraint
    ('update_constraints_column', 'fitted', 'param', 'nope', {'is_free': False}, REFUSED),
    ('update_constraints_column', 'fitted', 'param', '', {'is_free': False}, REFUSED),
    ('update_constraints_column', 'fitted', 'param', None, {'is_free': False}, REFUSED),  # 2.4.0: accepted
    ('add_linked_system', 'fitted', 'z', '2.0', {'series_list': 'CIV,SiIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_linked_system', 'fitted', 'z', 2.0, {'series_list': 'CIV,SiIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_linked_system', 'fitted', 'z', 2, {'series_list': 'CIV,SiIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_linked_system', 'fitted', 'z', ' 2 ', {'series_list': 'CIV,SiIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_linked_system', 'fitted', 'z', '1e1', {'series_list': 'CIV,SiIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_linked_system', 'fitted', 'z', 'abc', {'series_list': 'CIV,SiIV'}, REFUSED),  # 2.4.0: crashed (ValueError)
    ('add_linked_system', 'fitted', 'z', '', {'series_list': 'CIV,SiIV'}, REFUSED),  # 2.4.0: crashed (ValueError)
    ('add_linked_system', 'fitted', 'z', 'None', {'series_list': 'CIV,SiIV'}, REFUSED),  # 2.4.0: crashed (ValueError)
    ('add_linked_system', 'fitted', 'z', None, {'series_list': 'CIV,SiIV'}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('add_linked_system', 'fitted', 'z', 'nan', {'series_list': 'CIV,SiIV'}, REFUSED),  # 2.4.0: crashed (ValueError)
    ('add_linked_system', 'fitted', 'z', 'inf', {'series_list': 'CIV,SiIV'}, REFUSED),  # 2.4.0: crashed (ValueError)
    ('add_linked_system', 'fitted', 'z', True, {'series_list': 'CIV,SiIV'}, REFUSED),  # 2.4.0: accepted as 1
    ('add_linked_system', 'fitted', 'z', '-1', {'series_list': 'CIV,SiIV'}, REFUSED),  # 2.4.0: crashed (ValueError)
    ('add_linked_system', 'fitted', 'z', '0', {'series_list': 'CIV,SiIV'}, FAILS),  # F111: no pixel (outside the spectrum), fitted on the whole spectrum until PR 7-quater
    ('add_linked_system', 'fitted', 'z', '1.931', {'series_list': 'CIV,SiIV'}, 1.931),  # 2.4.0: crashed (ValueError)
    ('add_linked_system', 'fitted', 'logN', '2.0', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'logN', 2.0, {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_linked_system', 'fitted', 'logN', 2, {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_linked_system', 'fitted', 'logN', ' 2 ', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'logN', '1e1', {'series_list': 'CIV,SiIV', 'z': 1.931}, 10.0),  # 2.4.0: refused (0)
    ('add_linked_system', 'fitted', 'logN', 'abc', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'logN', '', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'logN', 'None', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'logN', None, {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_linked_system', 'fitted', 'logN', 'nan', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'logN', 'inf', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'logN', True, {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted as 1
    ('add_linked_system', 'fitted', 'logN', '-1', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'logN', '0', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'logN', '30', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'b', '2.0', {'series_list': 'CIV,SiIV', 'z': 1.931}, 2.0),  # 2.4.0: refused (0)
    ('add_linked_system', 'fitted', 'b', 2.0, {'series_list': 'CIV,SiIV', 'z': 1.931}, 2.0),
    ('add_linked_system', 'fitted', 'b', 2, {'series_list': 'CIV,SiIV', 'z': 1.931}, 2.0),
    ('add_linked_system', 'fitted', 'b', ' 2 ', {'series_list': 'CIV,SiIV', 'z': 1.931}, 2.0),  # 2.4.0: refused (0)
    ('add_linked_system', 'fitted', 'b', '1e1', {'series_list': 'CIV,SiIV', 'z': 1.931}, 10.0),  # 2.4.0: refused (0)
    ('add_linked_system', 'fitted', 'b', 'abc', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'b', '', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'b', 'None', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'b', None, {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'b', 'nan', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'b', 'inf', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'b', True, {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted as 1
    ('add_linked_system', 'fitted', 'b', '-1', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'b', '0', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'b', '30', {'series_list': 'CIV,SiIV', 'z': 1.931}, 30.0),  # 2.4.0: refused (0)
    ('add_linked_system', 'fitted', 'z_window_kms', '2.0', {'series_list': 'CIV,SiIV', 'z': 1.931}, 2.0),
    ('add_linked_system', 'fitted', 'z_window_kms', 2.0, {'series_list': 'CIV,SiIV', 'z': 1.931}, 2.0),
    ('add_linked_system', 'fitted', 'z_window_kms', 2, {'series_list': 'CIV,SiIV', 'z': 1.931}, 2.0),
    ('add_linked_system', 'fitted', 'z_window_kms', ' 2 ', {'series_list': 'CIV,SiIV', 'z': 1.931}, 2.0),
    ('add_linked_system', 'fitted', 'z_window_kms', '1e1', {'series_list': 'CIV,SiIV', 'z': 1.931}, 10.0),
    ('add_linked_system', 'fitted', 'z_window_kms', 'abc', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'z_window_kms', '', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'z_window_kms', 'None', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'z_window_kms', None, {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'z_window_kms', 'nan', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'z_window_kms', 'inf', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'z_window_kms', True, {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'z_window_kms', '-1', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_linked_system', 'fitted', 'z_window_kms', '0', {'series_list': 'CIV,SiIV', 'z': 1.931}, REFUSED),  # 2.4.0: accepted
    ('add_linked_system', 'fitted', 'z_window_kms', '30', {'series_list': 'CIV,SiIV', 'z': 1.931}, 30.0),
    ('add_linked_system', 'fitted', 'series_list', 'CIV', {'z': 1.931}, 'CIV'),
    ('add_linked_system', 'fitted', 'series_list', 'CIV,SiIV', {'z': 1.931}, 'CIV,SiIV'),
    ('add_linked_system', 'fitted', 'series_list', 'nope', {'z': 1.931}, FAILS),  # F111: no transition, so no pixel; fitted on the whole spectrum until PR 7-quater
    ('add_linked_system', 'fitted', 'series_list', '', {'z': 1.931}, REFUSED),
    ('add_linked_system', 'fitted', 'series_list', None, {'z': 1.931}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('optimize_system', 'fitted', 'max_components', '2', {'uuid': UUID}, 2),
    ('optimize_system', 'fitted', 'max_components', 2, {'uuid': UUID}, 2),
    ('optimize_system', 'fitted', 'max_components', '2.0', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'max_components', 2.0, {'uuid': UUID}, 2),
    ('optimize_system', 'fitted', 'max_components', 2.5, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted, truncated
    ('optimize_system', 'fitted', 'max_components', '', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'max_components', 'abc', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'max_components', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('optimize_system', 'fitted', 'max_components', None, {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'max_components', 'None', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'max_components', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('optimize_system', 'fitted', 'max_components', '0', {'uuid': UUID}, 0),
    ('optimize_system', 'fitted', 'max_components', '1', {'uuid': UUID}, 1),
    ('optimize_system', 'fitted', 'max_components', '3', {'uuid': UUID}, 3),
    ('optimize_system', 'fitted', 'group_depth', '2', {'uuid': UUID}, 2),
    ('optimize_system', 'fitted', 'group_depth', 2, {'uuid': UUID}, 2),
    ('optimize_system', 'fitted', 'group_depth', '2.0', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'group_depth', 2.0, {'uuid': UUID}, 2),
    ('optimize_system', 'fitted', 'group_depth', 2.5, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted, truncated
    ('optimize_system', 'fitted', 'group_depth', '', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'group_depth', 'abc', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'group_depth', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('optimize_system', 'fitted', 'group_depth', None, {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'group_depth', 'None', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'group_depth', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('optimize_system', 'fitted', 'group_depth', '0', {'uuid': UUID}, 0),
    ('optimize_system', 'fitted', 'group_depth', '1', {'uuid': UUID}, 1),
    ('optimize_system', 'fitted', 'group_depth', '3', {'uuid': UUID}, 3),
    ('optimize_system', 'fitted', 'patience', '2', {'uuid': UUID}, 2),
    ('optimize_system', 'fitted', 'patience', 2, {'uuid': UUID}, 2),
    ('optimize_system', 'fitted', 'patience', '2.0', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'patience', 2.0, {'uuid': UUID}, 2),
    ('optimize_system', 'fitted', 'patience', 2.5, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted, truncated
    ('optimize_system', 'fitted', 'patience', '', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'patience', 'abc', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'patience', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('optimize_system', 'fitted', 'patience', None, {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'patience', 'None', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'patience', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('optimize_system', 'fitted', 'patience', '0', {'uuid': UUID}, 0),
    ('optimize_system', 'fitted', 'patience', '1', {'uuid': UUID}, 1),
    ('optimize_system', 'fitted', 'patience', '3', {'uuid': UUID}, 3),
    ('optimize_system', 'fitted', 'threshold_sigma', '2.0', {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'threshold_sigma', 2.0, {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'threshold_sigma', 2, {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'threshold_sigma', ' 2 ', {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'threshold_sigma', '1e1', {'uuid': UUID}, 10.0),
    ('optimize_system', 'fitted', 'threshold_sigma', 'abc', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'threshold_sigma', '', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'threshold_sigma', 'None', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'threshold_sigma', None, {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'threshold_sigma', 'nan', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as NaN
    ('optimize_system', 'fitted', 'threshold_sigma', 'inf', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as an infinity
    ('optimize_system', 'fitted', 'threshold_sigma', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('optimize_system', 'fitted', 'threshold_sigma', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('optimize_system', 'fitted', 'threshold_sigma', '0', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('optimize_system', 'fitted', 'threshold_sigma', '30', {'uuid': UUID}, 30.0),
    ('optimize_system', 'fitted', 'aic_penalty', '2.0', {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'aic_penalty', 2.0, {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'aic_penalty', 2, {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'aic_penalty', ' 2 ', {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'aic_penalty', '1e1', {'uuid': UUID}, 10.0),
    ('optimize_system', 'fitted', 'aic_penalty', 'abc', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'aic_penalty', '', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'aic_penalty', 'None', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'aic_penalty', None, {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'aic_penalty', 'nan', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as NaN
    ('optimize_system', 'fitted', 'aic_penalty', 'inf', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as an infinity
    ('optimize_system', 'fitted', 'aic_penalty', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('optimize_system', 'fitted', 'aic_penalty', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('optimize_system', 'fitted', 'aic_penalty', '0', {'uuid': UUID}, 0.0),
    ('optimize_system', 'fitted', 'aic_penalty', '30', {'uuid': UUID}, 30.0),
    ('optimize_system', 'fitted', 'z_window_kms', '2.0', {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'z_window_kms', 2.0, {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'z_window_kms', 2, {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'z_window_kms', ' 2 ', {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'z_window_kms', '1e1', {'uuid': UUID}, 10.0),
    ('optimize_system', 'fitted', 'z_window_kms', 'abc', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'z_window_kms', '', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'z_window_kms', 'None', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'z_window_kms', None, {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'z_window_kms', 'nan', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as NaN
    ('optimize_system', 'fitted', 'z_window_kms', 'inf', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'z_window_kms', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('optimize_system', 'fitted', 'z_window_kms', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('optimize_system', 'fitted', 'z_window_kms', '0', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('optimize_system', 'fitted', 'z_window_kms', '30', {'uuid': UUID}, 30.0),
    ('optimize_system', 'fitted', 'min_dv', '2.0', {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'min_dv', 2.0, {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'min_dv', 2, {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'min_dv', ' 2 ', {'uuid': UUID}, 2.0),
    ('optimize_system', 'fitted', 'min_dv', '1e1', {'uuid': UUID}, 10.0),
    ('optimize_system', 'fitted', 'min_dv', 'abc', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'min_dv', '', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'min_dv', 'None', {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'min_dv', None, {'uuid': UUID}, REFUSED),
    ('optimize_system', 'fitted', 'min_dv', 'nan', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as NaN
    ('optimize_system', 'fitted', 'min_dv', 'inf', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as an infinity
    ('optimize_system', 'fitted', 'min_dv', True, {'uuid': UUID}, REFUSED),  # 2.4.0: accepted as 1
    ('optimize_system', 'fitted', 'min_dv', '-1', {'uuid': UUID}, REFUSED),  # 2.4.0: accepted
    ('optimize_system', 'fitted', 'min_dv', '0', {'uuid': UUID}, 0.0),
    ('optimize_system', 'fitted', 'min_dv', '30', {'uuid': UUID}, 30.0),
    ('fit_bayesian', 'fitted', 'nsteps', '', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'nsteps', 'abc', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'nsteps', None, {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'nsteps', 'None', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'nsteps', '-1', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'nsteps', '0', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'nsteps', 2.5, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted, truncated
    ('fit_bayesian', 'fitted', 'nsteps', True, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as 1
    ('fit_bayesian', 'fitted', 'nlive', '', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'nlive', 'abc', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'nlive', None, {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'nlive', 'None', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'nlive', '-1', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('fit_bayesian', 'fitted', 'nlive', '0', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('fit_bayesian', 'fitted', 'nlive', 2.5, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted, truncated
    ('fit_bayesian', 'fitted', 'nlive', True, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as 1
    ('fit_bayesian', 'fitted', 'group_depth', '', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'group_depth', 'abc', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'group_depth', None, {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'group_depth', 'None', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'group_depth', '-1', {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted
    ('fit_bayesian', 'fitted', 'group_depth', '0', {'uuids': [UUID]}, ACCEPTED),
    ('fit_bayesian', 'fitted', 'group_depth', 2.5, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted, truncated
    ('fit_bayesian', 'fitted', 'group_depth', True, {'uuids': [UUID]}, REFUSED),  # 2.4.0: accepted as 1
    # Refused since phase 6-bis (PR 6c): an unknown engine ran MCMC
    ('fit_bayesian', 'fitted', 'engine', '', {'uuids': [UUID]}, REFUSED),
    ('fit_bayesian', 'fitted', 'engine', None, {'uuids': [UUID]}, REFUSED),
    ('refit_all', 'fitted', 'max_nfev', '2', {}, 2),
    ('refit_all', 'fitted', 'max_nfev', 2, {}, 2),
    ('refit_all', 'fitted', 'max_nfev', '2.0', {}, REFUSED),
    ('refit_all', 'fitted', 'max_nfev', 2.0, {}, 2),
    ('refit_all', 'fitted', 'max_nfev', 2.5, {}, REFUSED),  # 2.4.0: accepted, truncated
    ('refit_all', 'fitted', 'max_nfev', '', {}, REFUSED),
    ('refit_all', 'fitted', 'max_nfev', 'abc', {}, REFUSED),
    ('refit_all', 'fitted', 'max_nfev', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('refit_all', 'fitted', 'max_nfev', None, {}, REFUSED),
    ('refit_all', 'fitted', 'max_nfev', 'None', {}, REFUSED),
    ('refit_all', 'fitted', 'max_nfev', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_nfev', '0', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_nfev', '1', {}, 1),
    ('refit_all', 'fitted', 'max_nfev', '3', {}, 3),
    ('refit_all', 'fitted', 'group_depth', '2', {}, 2),
    ('refit_all', 'fitted', 'group_depth', 2, {}, 2),
    ('refit_all', 'fitted', 'group_depth', '2.0', {}, REFUSED),
    ('refit_all', 'fitted', 'group_depth', 2.0, {}, 2),
    ('refit_all', 'fitted', 'group_depth', 2.5, {}, REFUSED),  # 2.4.0: accepted, truncated
    ('refit_all', 'fitted', 'group_depth', '', {}, REFUSED),
    ('refit_all', 'fitted', 'group_depth', 'abc', {}, REFUSED),
    ('refit_all', 'fitted', 'group_depth', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('refit_all', 'fitted', 'group_depth', None, {}, REFUSED),
    ('refit_all', 'fitted', 'group_depth', 'None', {}, REFUSED),
    ('refit_all', 'fitted', 'group_depth', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'group_depth', '0', {}, 0),
    ('refit_all', 'fitted', 'group_depth', '1', {}, 1),
    ('refit_all', 'fitted', 'group_depth', '3', {}, 3),
    ('refit_all', 'fitted', 'max_group_size', '2', {}, 2),
    ('refit_all', 'fitted', 'max_group_size', 2, {}, 2),
    ('refit_all', 'fitted', 'max_group_size', '2.0', {}, REFUSED),
    ('refit_all', 'fitted', 'max_group_size', 2.0, {}, 2),
    ('refit_all', 'fitted', 'max_group_size', 2.5, {}, REFUSED),  # 2.4.0: accepted, truncated
    ('refit_all', 'fitted', 'max_group_size', '', {}, REFUSED),
    ('refit_all', 'fitted', 'max_group_size', 'abc', {}, REFUSED),
    ('refit_all', 'fitted', 'max_group_size', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('refit_all', 'fitted', 'max_group_size', None, {}, REFUSED),
    ('refit_all', 'fitted', 'max_group_size', 'None', {}, REFUSED),
    ('refit_all', 'fitted', 'max_group_size', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_group_size', '0', {}, REFUSED),
    ('refit_all', 'fitted', 'max_group_size', '1', {}, 1),
    ('refit_all', 'fitted', 'max_group_size', '3', {}, 3),
    ('refit_all', 'fitted', 'z_window_kms', '2.0', {}, 2.0),
    ('refit_all', 'fitted', 'z_window_kms', 2.0, {}, 2.0),
    ('refit_all', 'fitted', 'z_window_kms', 2, {}, 2.0),
    ('refit_all', 'fitted', 'z_window_kms', ' 2 ', {}, 2.0),
    ('refit_all', 'fitted', 'z_window_kms', '1e1', {}, 10.0),
    ('refit_all', 'fitted', 'z_window_kms', 'abc', {}, REFUSED),
    ('refit_all', 'fitted', 'z_window_kms', '', {}, REFUSED),
    ('refit_all', 'fitted', 'z_window_kms', 'None', {}, REFUSED),
    ('refit_all', 'fitted', 'z_window_kms', None, {}, REFUSED),
    ('refit_all', 'fitted', 'z_window_kms', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('refit_all', 'fitted', 'z_window_kms', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('refit_all', 'fitted', 'z_window_kms', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('refit_all', 'fitted', 'z_window_kms', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'z_window_kms', '0', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'z_window_kms', '30', {}, 30.0),
    ('refit_all', 'fitted', 'max_b', '2.0', {}, 2.0),
    ('refit_all', 'fitted', 'max_b', 2.0, {}, 2.0),
    ('refit_all', 'fitted', 'max_b', 2, {}, 2.0),
    ('refit_all', 'fitted', 'max_b', ' 2 ', {}, 2.0),
    ('refit_all', 'fitted', 'max_b', '1e1', {}, 10.0),
    ('refit_all', 'fitted', 'max_b', 'abc', {}, REFUSED),
    ('refit_all', 'fitted', 'max_b', '', {}, REFUSED),
    ('refit_all', 'fitted', 'max_b', 'None', {}, REFUSED),
    ('refit_all', 'fitted', 'max_b', None, {}, REFUSED),
    ('refit_all', 'fitted', 'max_b', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('refit_all', 'fitted', 'max_b', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('refit_all', 'fitted', 'max_b', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('refit_all', 'fitted', 'max_b', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_b', '0', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_b', '30', {}, 30.0),
    ('refit_all', 'fitted', 'max_logN', '2.0', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_logN', 2.0, {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_logN', 2, {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_logN', ' 2 ', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_logN', '1e1', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_logN', 'abc', {}, REFUSED),
    ('refit_all', 'fitted', 'max_logN', '', {}, REFUSED),
    ('refit_all', 'fitted', 'max_logN', 'None', {}, REFUSED),
    ('refit_all', 'fitted', 'max_logN', None, {}, REFUSED),
    ('refit_all', 'fitted', 'max_logN', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('refit_all', 'fitted', 'max_logN', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('refit_all', 'fitted', 'max_logN', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('refit_all', 'fitted', 'max_logN', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_logN', '0', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'max_logN', '30', {}, REFUSED),  # 2.4.0: accepted
    ('refit_all', 'fitted', 'region_limit', 'None', {}, 'None'),
    ('refit_all', 'fitted', 'region_limit', 'red_side', {}, 'red_side'),
    ('refit_all', 'fitted', 'region_limit', 'lya_forest', {}, 'lya_forest'),
    # Refused since phase 6-bis (PR 6c): a list of choices (they failed in the recipe)
    ('refit_all', 'fitted', 'region_limit', 'nope', {}, REFUSED),
    ('refit_all', 'fitted', 'region_limit', '', {}, REFUSED),
    ('refit_all', 'fitted', 'region_limit', None, {}, REFUSED),
    ('refit_all', 'fitted', 'selected_uuids', [UUID], {}, [UUID]),
    ('refit_all', 'fitted', 'selected_uuids', LIST_AS_TEXT, {}, [UUID]),  # 2.4.0: accepted, the text not read as a list
    ('refit_all', 'fitted', 'selected_uuids', [], {}, []),
    ('refit_all', 'fitted', 'selected_uuids', None, {}, None),
    ('refit_all', 'fitted', 'selected_uuids', 'abc', {}, REFUSED),  # 2.4.0: accepted
    ('clean_negligible', 'fitted', 'min_logN', '2.0', {}, 2.0),
    ('clean_negligible', 'fitted', 'min_logN', 2.0, {}, 2.0),
    ('clean_negligible', 'fitted', 'min_logN', 2, {}, 2.0),
    ('clean_negligible', 'fitted', 'min_logN', ' 2 ', {}, 2.0),
    ('clean_negligible', 'fitted', 'min_logN', '1e1', {}, 10.0),
    ('clean_negligible', 'fitted', 'min_logN', 'abc', {}, REFUSED),
    ('clean_negligible', 'fitted', 'min_logN', '', {}, REFUSED),
    ('clean_negligible', 'fitted', 'min_logN', 'None', {}, REFUSED),
    ('clean_negligible', 'fitted', 'min_logN', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('clean_negligible', 'fitted', 'min_logN', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('clean_negligible', 'fitted', 'min_logN', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('clean_negligible', 'fitted', 'min_logN', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('clean_negligible', 'fitted', 'min_logN', '-1', {}, -1.0),
    ('clean_negligible', 'fitted', 'min_logN', '0', {}, 0.0),
    ('clean_negligible', 'fitted', 'min_logN', '13', {}, 13.0),
    ('clean_negligible', 'fitted', 'max_b', '2.0', {}, 2.0),
    ('clean_negligible', 'fitted', 'max_b', 2.0, {}, 2.0),
    ('clean_negligible', 'fitted', 'max_b', 2, {}, 2.0),
    ('clean_negligible', 'fitted', 'max_b', ' 2 ', {}, 2.0),
    ('clean_negligible', 'fitted', 'max_b', '1e1', {}, 10.0),
    ('clean_negligible', 'fitted', 'max_b', 'abc', {}, REFUSED),
    ('clean_negligible', 'fitted', 'max_b', '', {}, REFUSED),
    ('clean_negligible', 'fitted', 'max_b', 'None', {}, REFUSED),
    ('clean_negligible', 'fitted', 'max_b', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('clean_negligible', 'fitted', 'max_b', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('clean_negligible', 'fitted', 'max_b', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('clean_negligible', 'fitted', 'max_b', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('clean_negligible', 'fitted', 'max_b', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('clean_negligible', 'fitted', 'max_b', '0', {}, REFUSED),  # 2.4.0: accepted
    ('clean_negligible', 'fitted', 'max_b', '13', {}, 13.0),
    ('clean_negligible', 'fitted', 'min_b', '2.0', {}, 2.0),
    ('clean_negligible', 'fitted', 'min_b', 2.0, {}, 2.0),
    ('clean_negligible', 'fitted', 'min_b', 2, {}, 2.0),
    ('clean_negligible', 'fitted', 'min_b', ' 2 ', {}, 2.0),
    ('clean_negligible', 'fitted', 'min_b', '1e1', {}, 10.0),
    ('clean_negligible', 'fitted', 'min_b', 'abc', {}, REFUSED),
    ('clean_negligible', 'fitted', 'min_b', '', {}, REFUSED),
    ('clean_negligible', 'fitted', 'min_b', 'None', {}, REFUSED),
    ('clean_negligible', 'fitted', 'min_b', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('clean_negligible', 'fitted', 'min_b', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('clean_negligible', 'fitted', 'min_b', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('clean_negligible', 'fitted', 'min_b', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('clean_negligible', 'fitted', 'min_b', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('clean_negligible', 'fitted', 'min_b', '0', {}, 0.0),
    ('clean_negligible', 'fitted', 'min_b', '13', {}, 13.0),
    ('clean_negligible', 'fitted', 'combined_check', 'True', {}, True),
    ('clean_negligible', 'fitted', 'combined_check', 'False', {}, False),
    ('clean_negligible', 'fitted', 'combined_check', 'true', {}, True),
    ('clean_negligible', 'fitted', 'combined_check', 'false', {}, False),
    ('clean_negligible', 'fitted', 'combined_check', '1', {}, True),
    ('clean_negligible', 'fitted', 'combined_check', '0', {}, False),
    ('clean_negligible', 'fitted', 'combined_check', 'yes', {}, REFUSED),  # 2.4.0: accepted as False
    ('clean_negligible', 'fitted', 'combined_check', True, {}, True),  # 2.4.0: crashed (AttributeError)
    ('clean_negligible', 'fitted', 'combined_check', False, {}, False),  # 2.4.0: crashed (AttributeError)
    ('clean_negligible', 'fitted', 'combined_check', 1, {}, True),  # 2.4.0: crashed (AttributeError)
    ('clean_negligible', 'fitted', 'combined_check', 0, {}, False),  # 2.4.0: crashed (AttributeError)
    ('clean_negligible', 'fitted', 'combined_check', '', {}, REFUSED),  # 2.4.0: accepted as False
    ('clean_negligible', 'fitted', 'combined_check', None, {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('clean_negligible', 'fitted', 'combined_check', 'None', {}, REFUSED),  # 2.4.0: accepted as False
]


@pytest.mark.parametrize("name, session, param, value, base, expected", CASES,
                         ids=[f"{n}.{p}={v!r}" for n, _, p, v, _, _ in CASES])
def test_input(caplog, name, session, param, value, base, expected):
    # The outcome is characterized here, not the warnings of the libraries
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _check(caplog, name, session, param, value, base, expected)


def _check(caplog, name, session, param, value, base, expected):
    if expected in FAILURES:
        check_failure(lambda: _run(name, session, param, value, base), name, param,
                      expected)
        return
    with caplog.at_level(logging.ERROR):
        result = _run(name, session, param, value, base)
    assert isinstance(result, SessionV2), caplog.text
    if expected == ACCEPTED:
        return
    reference = _reference(name, session, param, _freeze(expected),
                           tuple(sorted((k, _freeze(v)) for k, v in base.items())))
    assert _fingerprint(result) == reference
