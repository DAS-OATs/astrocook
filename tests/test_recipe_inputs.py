"""What the typed recipes do with each kind of input (phase 6-bis, PR 2a).

A battery of edge cases, the same as the comparison with 2.4.0 made for
the PR: for each value given to a parameter, the value the recipe uses
(the recipe must give the same result as with that typed value), REFUSED
(the conversion stops the recipe: 0 and an error that names the
parameter) or FAILS (the value is converted, then the recipe itself
fails, e.g. on a session that is not open). The
cases whose outcome changed from 2.4.0 say what 2.4.0 did; every other
case behaved the same. Decisions of Guido: 1 and 0 are booleans; an empty
field is None only for a parameter whose default is None (an empty kappa
of Rebin is an error, not "no clipping"); widths and thresholds positive
and finite, and no infinity as a value (PR 2b).

The recipes of Edit are in ``test_recipe_inputs_edit.py``; those of
Absorbers are added when they are typed (PR 2c), each with its comparison
with 2.4.0.
"""
import logging
import warnings
from types import SimpleNamespace

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.session import SessionV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, SpectrumDataV2

REFUSED = "refused"  # by the conversion: 0, an error that names the parameter
FAILS = "fails"  # converted, then the recipe fails (e.g. a session not open): 0
RAISES = "raises"  # converted, then an exception leaves the recipe (PR 3 turns them into errors)
NAN = float("nan")
FLUX = au.erg / au.cm**2 / au.s / au.nm
CATEGORY = {"estimate_auto": "continuum", "find_absorbed": "continuum",
            "fit_continuum": "continuum", "fit_powerlaw": "continuum",
            "update_from_knots": "continuum", "calculate_running_std": "flux",
            "smooth": "flux", "rebin": "flux", "equalize": "flux", "resample": "flux"}


def _spec(scale=1.0):
    rng = np.random.default_rng(0)
    x = np.arange(380.0, 460.0, 0.02)
    cont = scale * (x / 420.0) ** -1.5 * (1 + 0.1 * (scale - 1) * (x - 420.0) / 40.0)
    y = (cont * (1 - 0.5 * np.exp(-0.5 * ((x - 420.0) / 0.05) ** 2))
         + 0.02 * cont * rng.standard_normal(len(x)))
    return SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - 0.01, au.nm),
        xmax=DataColumnV2(x + 0.01, au.nm), y=DataColumnV2(y, FLUX),
        dy=DataColumnV2(0.02 * cont, FLUX), aux_cols={}, z_em=2.5))


@pytest.fixture(scope="module")
def sessions():
    """The session each recipe runs on (fit_continuum needs a mask)."""
    plain = SessionV2(name="main", gui=None, spec=_spec())
    other = SessionV2(name="other", gui=None, spec=_spec(2.0))  # another level and slope
    plain._gui = SimpleNamespace(session_histories=[
        SimpleNamespace(display_name="other", current_state=other)])
    masked = plain.continuum.find_absorbed()
    masked._gui = plain._gui
    return {"plain": plain, "masked": masked}


def _run(sessions, name, param, value):
    session = sessions["masked" if name == "fit_continuum" else "plain"]
    kwargs = {param: value}
    if name == "update_from_knots":
        kwargs.update(knots_x=[390.0, 420.0, 450.0], knots_y=[1.0, 1.0, 1.0])
    if name == "equalize" and param == "order":
        kwargs["reference_session"] = "other"
    if name == "rebin" and param != "dx":
        kwargs.setdefault("dx", 20.0)
    return getattr(getattr(session, CATEGORY[name]), name)(**kwargs)


def _data(result):
    columns = [np.asarray(result.spec.y.value)]
    for c in ("cont", "abs_mask", "running_std", "cont_pl", "model"):
        if result.spec.has_aux_column(c):
            columns.append(np.asarray(result.spec.t[c], dtype=float))
    return columns


CASES = [
    ('estimate_auto', 'smooth_len_lya', '2.0', 2.0),
    ('estimate_auto', 'smooth_len_lya', 2.0, 2.0),
    ('estimate_auto', 'smooth_len_lya', 2, 2.0),
    ('estimate_auto', 'smooth_len_lya', ' 2 ', 2.0),
    ('estimate_auto', 'smooth_len_lya', '1e1', 10.0),
    ('estimate_auto', 'smooth_len_lya', 'abc', REFUSED),
    ('estimate_auto', 'smooth_len_lya', '', REFUSED),
    ('estimate_auto', 'smooth_len_lya', 'None', REFUSED),
    ('estimate_auto', 'smooth_len_lya', None, REFUSED),  # 2.4.0: crashed (TypeError)
    ('estimate_auto', 'smooth_len_lya', 'nan', REFUSED),
    ('estimate_auto', 'smooth_len_lya', 'inf', REFUSED),
    ('estimate_auto', 'smooth_len_lya', '-inf', REFUSED),
    ('estimate_auto', 'smooth_len_lya', True, REFUSED),  # 2.4.0: accepted as 1
    ('estimate_auto', 'smooth_len_lya', '-1', REFUSED),  # 2.4.0: accepted
    ('estimate_auto', 'smooth_len_lya', '0', REFUSED),  # 2.4.0: accepted
    ('estimate_auto', 'smooth_len_lya', '3,5', REFUSED),
    ('estimate_auto', 'smooth_len_out', '2.0', 2.0),
    ('estimate_auto', 'smooth_len_out', 2.0, 2.0),
    ('estimate_auto', 'smooth_len_out', 2, 2.0),
    ('estimate_auto', 'smooth_len_out', ' 2 ', 2.0),
    ('estimate_auto', 'smooth_len_out', '1e1', 10.0),
    ('estimate_auto', 'smooth_len_out', 'abc', REFUSED),
    ('estimate_auto', 'smooth_len_out', '', REFUSED),
    ('estimate_auto', 'smooth_len_out', 'None', REFUSED),
    ('estimate_auto', 'smooth_len_out', None, REFUSED),  # 2.4.0: crashed (TypeError)
    ('estimate_auto', 'smooth_len_out', 'nan', REFUSED),
    ('estimate_auto', 'smooth_len_out', 'inf', REFUSED),
    ('estimate_auto', 'smooth_len_out', '-inf', REFUSED),
    ('estimate_auto', 'smooth_len_out', True, REFUSED),  # 2.4.0: accepted as 1
    ('estimate_auto', 'smooth_len_out', '-1', REFUSED),  # 2.4.0: accepted
    ('estimate_auto', 'smooth_len_out', '0', REFUSED),  # 2.4.0: accepted
    ('estimate_auto', 'smooth_len_out', '3,5', REFUSED),
    ('estimate_auto', 'kappa', '2.0', 2.0),
    ('estimate_auto', 'kappa', 2.0, 2.0),
    ('estimate_auto', 'kappa', 2, 2.0),
    ('estimate_auto', 'kappa', ' 2 ', 2.0),
    ('estimate_auto', 'kappa', '1e1', 10.0),
    ('estimate_auto', 'kappa', 'abc', REFUSED),
    ('estimate_auto', 'kappa', '', REFUSED),
    ('estimate_auto', 'kappa', 'None', REFUSED),
    ('estimate_auto', 'kappa', None, REFUSED),  # 2.4.0: crashed (TypeError)
    ('estimate_auto', 'kappa', 'nan', REFUSED),  # 2.4.0: accepted as NaN
    ('estimate_auto', 'kappa', 'inf', REFUSED),  # 2.4.0: accepted as an infinity
    ('estimate_auto', 'kappa', '-inf', REFUSED),
    ('estimate_auto', 'kappa', True, REFUSED),  # 2.4.0: accepted as 1
    ('estimate_auto', 'kappa', '-1', REFUSED),  # 2.4.0: accepted
    ('estimate_auto', 'kappa', '0', REFUSED),  # 2.4.0: accepted
    ('estimate_auto', 'kappa', '3,5', REFUSED),
    ('estimate_auto', 'smooth_fwhm', '2.0', 2.0),
    ('estimate_auto', 'smooth_fwhm', 2.0, 2.0),
    ('estimate_auto', 'smooth_fwhm', 2, 2.0),
    ('estimate_auto', 'smooth_fwhm', ' 2 ', 2.0),
    ('estimate_auto', 'smooth_fwhm', '1e1', 10.0),
    ('estimate_auto', 'smooth_fwhm', 'abc', REFUSED),
    ('estimate_auto', 'smooth_fwhm', '', REFUSED),
    ('estimate_auto', 'smooth_fwhm', 'None', REFUSED),
    ('estimate_auto', 'smooth_fwhm', None, REFUSED),  # 2.4.0: crashed (TypeError)
    ('estimate_auto', 'smooth_fwhm', 'nan', REFUSED),  # 2.4.0: accepted as NaN
    ('estimate_auto', 'smooth_fwhm', 'inf', REFUSED),
    ('estimate_auto', 'smooth_fwhm', '-inf', REFUSED),  # 2.4.0: accepted as an infinity
    ('estimate_auto', 'smooth_fwhm', True, REFUSED),  # 2.4.0: accepted as 1
    ('estimate_auto', 'smooth_fwhm', '-1', REFUSED),  # 2.4.0: accepted
    ('estimate_auto', 'smooth_fwhm', '0', REFUSED),  # 2.4.0: accepted
    ('estimate_auto', 'smooth_fwhm', '3,5', REFUSED),
    ('estimate_auto', 'template', 'True', True),
    ('estimate_auto', 'template', 'False', False),
    ('estimate_auto', 'template', 'true', True),  # 2.4.0: accepted as False
    ('estimate_auto', 'template', 'false', False),
    ('estimate_auto', 'template', 'TRUE', True),  # 2.4.0: accepted as False
    ('estimate_auto', 'template', '1', True),  # 2.4.0: accepted as False
    ('estimate_auto', 'template', '0', False),
    ('estimate_auto', 'template', 'yes', REFUSED),  # 2.4.0: accepted as False
    ('estimate_auto', 'template', True, True),
    ('estimate_auto', 'template', False, False),
    ('estimate_auto', 'template', 1, True),  # 2.4.0: accepted as False
    ('estimate_auto', 'template', 0, False),
    ('estimate_auto', 'template', '', REFUSED),  # 2.4.0: accepted as False
    ('estimate_auto', 'template', None, REFUSED),  # 2.4.0: accepted as False
    ('estimate_auto', 'template', 'None', REFUSED),  # 2.4.0: accepted as False
    ('estimate_auto', 'renorm_model', 'True', True),
    ('estimate_auto', 'renorm_model', 'False', False),
    ('estimate_auto', 'renorm_model', 'true', True),
    ('estimate_auto', 'renorm_model', 'false', False),
    ('estimate_auto', 'renorm_model', 'TRUE', True),
    ('estimate_auto', 'renorm_model', '1', True),
    ('estimate_auto', 'renorm_model', '0', False),
    ('estimate_auto', 'renorm_model', 'yes', REFUSED),  # 2.4.0: accepted as False
    ('estimate_auto', 'renorm_model', True, True),
    ('estimate_auto', 'renorm_model', False, False),
    ('estimate_auto', 'renorm_model', 1, True),
    ('estimate_auto', 'renorm_model', 0, False),
    ('estimate_auto', 'renorm_model', '', REFUSED),  # 2.4.0: accepted as False
    ('estimate_auto', 'renorm_model', None, REFUSED),  # 2.4.0: accepted as False
    ('estimate_auto', 'renorm_model', 'None', REFUSED),  # 2.4.0: accepted as False
    ('estimate_auto', 'fudge', 'auto', 'auto'),
    ('estimate_auto', 'fudge', 'AUTO', 'AUTO'),
    ('estimate_auto', 'fudge', '1.0', '1.0'),
    ('estimate_auto', 'fudge', 1.0, '1.0'),  # 2.4.0: crashed (AttributeError)
    ('estimate_auto', 'fudge', 1, '1'),  # 2.4.0: crashed (AttributeError)
    ('estimate_auto', 'fudge', 'abc', FAILS),
    ('estimate_auto', 'fudge', None, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('estimate_auto', 'fudge', '', FAILS),
    ('find_absorbed', 'smooth_len_lya', '2.0', 2.0),
    ('find_absorbed', 'smooth_len_lya', 2.0, 2.0),
    ('find_absorbed', 'smooth_len_lya', 2, 2.0),
    ('find_absorbed', 'smooth_len_lya', ' 2 ', 2.0),
    ('find_absorbed', 'smooth_len_lya', '1e1', 10.0),
    ('find_absorbed', 'smooth_len_lya', 'abc', REFUSED),
    ('find_absorbed', 'smooth_len_lya', '', REFUSED),
    ('find_absorbed', 'smooth_len_lya', 'None', REFUSED),
    ('find_absorbed', 'smooth_len_lya', None, REFUSED),  # 2.4.0: crashed (TypeError)
    ('find_absorbed', 'smooth_len_lya', 'nan', REFUSED),
    ('find_absorbed', 'smooth_len_lya', 'inf', REFUSED),
    ('find_absorbed', 'smooth_len_lya', '-inf', REFUSED),
    ('find_absorbed', 'smooth_len_lya', True, REFUSED),  # 2.4.0: accepted as 1
    ('find_absorbed', 'smooth_len_lya', '-1', REFUSED),  # 2.4.0: accepted
    ('find_absorbed', 'smooth_len_lya', '0', REFUSED),  # 2.4.0: accepted
    ('find_absorbed', 'smooth_len_lya', '3,5', REFUSED),
    ('find_absorbed', 'kappa', '2.0', 2.0),
    ('find_absorbed', 'kappa', 2.0, 2.0),
    ('find_absorbed', 'kappa', 2, 2.0),
    ('find_absorbed', 'kappa', ' 2 ', 2.0),
    ('find_absorbed', 'kappa', '1e1', 10.0),
    ('find_absorbed', 'kappa', 'abc', REFUSED),
    ('find_absorbed', 'kappa', '', REFUSED),
    ('find_absorbed', 'kappa', 'None', REFUSED),
    ('find_absorbed', 'kappa', None, REFUSED),  # 2.4.0: crashed (TypeError)
    ('find_absorbed', 'kappa', 'nan', REFUSED),  # 2.4.0: accepted as NaN
    ('find_absorbed', 'kappa', 'inf', REFUSED),  # 2.4.0: accepted as an infinity
    ('find_absorbed', 'kappa', '-inf', REFUSED),  # 2.4.0: accepted as an infinity
    ('find_absorbed', 'kappa', True, REFUSED),  # 2.4.0: accepted as 1
    ('find_absorbed', 'kappa', '-1', REFUSED),  # 2.4.0: accepted
    ('find_absorbed', 'kappa', '0', REFUSED),  # 2.4.0: accepted
    ('find_absorbed', 'kappa', '3,5', REFUSED),
    ('find_absorbed', 'template', 'True', True),
    ('find_absorbed', 'template', 'False', False),
    ('find_absorbed', 'template', 'true', True),
    ('find_absorbed', 'template', 'false', False),
    ('find_absorbed', 'template', 'TRUE', True),
    ('find_absorbed', 'template', '1', True),
    ('find_absorbed', 'template', '0', False),
    ('find_absorbed', 'template', 'yes', REFUSED),  # 2.4.0: accepted as False
    ('find_absorbed', 'template', True, True),
    ('find_absorbed', 'template', False, False),
    ('find_absorbed', 'template', 1, True),
    ('find_absorbed', 'template', 0, False),
    ('find_absorbed', 'template', '', REFUSED),  # 2.4.0: accepted as False
    ('find_absorbed', 'template', None, REFUSED),  # 2.4.0: accepted as False
    ('find_absorbed', 'template', 'None', REFUSED),  # 2.4.0: accepted as False
    ('fit_continuum', 'smooth_fwhm', '2.0', 2.0),
    ('fit_continuum', 'smooth_fwhm', 2.0, 2.0),
    ('fit_continuum', 'smooth_fwhm', 2, 2.0),
    ('fit_continuum', 'smooth_fwhm', ' 2 ', 2.0),
    ('fit_continuum', 'smooth_fwhm', '1e1', 10.0),
    ('fit_continuum', 'smooth_fwhm', 'abc', REFUSED),
    ('fit_continuum', 'smooth_fwhm', '', REFUSED),
    ('fit_continuum', 'smooth_fwhm', 'None', REFUSED),
    ('fit_continuum', 'smooth_fwhm', None, REFUSED),  # 2.4.0: crashed (TypeError)
    ('fit_continuum', 'smooth_fwhm', 'nan', REFUSED),  # 2.4.0: accepted as NaN
    ('fit_continuum', 'smooth_fwhm', 'inf', REFUSED),
    ('fit_continuum', 'smooth_fwhm', '-inf', REFUSED),  # 2.4.0: accepted as an infinity
    ('fit_continuum', 'smooth_fwhm', True, REFUSED),  # 2.4.0: accepted as 1
    ('fit_continuum', 'smooth_fwhm', '-1', REFUSED),  # 2.4.0: accepted
    ('fit_continuum', 'smooth_fwhm', '0', REFUSED),  # 2.4.0: accepted
    ('fit_continuum', 'smooth_fwhm', '3,5', REFUSED),
    ('fit_continuum', 'renorm_model', 'True', True),
    ('fit_continuum', 'renorm_model', 'False', False),
    ('fit_continuum', 'renorm_model', 'true', True),
    ('fit_continuum', 'renorm_model', 'false', False),
    ('fit_continuum', 'renorm_model', 'TRUE', True),
    ('fit_continuum', 'renorm_model', '1', True),
    ('fit_continuum', 'renorm_model', '0', False),
    ('fit_continuum', 'renorm_model', 'yes', REFUSED),  # 2.4.0: accepted as False
    ('fit_continuum', 'renorm_model', True, True),
    ('fit_continuum', 'renorm_model', False, False),
    ('fit_continuum', 'renorm_model', 1, True),
    ('fit_continuum', 'renorm_model', 0, False),
    ('fit_continuum', 'renorm_model', '', REFUSED),  # 2.4.0: accepted as False
    ('fit_continuum', 'renorm_model', None, REFUSED),  # 2.4.0: accepted as False
    ('fit_continuum', 'renorm_model', 'None', REFUSED),  # 2.4.0: accepted as False
    ('fit_continuum', 'fudge', 'auto', 'auto'),
    ('fit_continuum', 'fudge', 'AUTO', 'AUTO'),
    ('fit_continuum', 'fudge', '1.0', '1.0'),
    ('fit_continuum', 'fudge', 1.0, '1.0'),  # 2.4.0: crashed (AttributeError)
    ('fit_continuum', 'fudge', 1, '1'),  # 2.4.0: crashed (AttributeError)
    ('fit_continuum', 'fudge', 'abc', FAILS),
    ('fit_continuum', 'fudge', None, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('fit_continuum', 'fudge', '', FAILS),
    ('fit_powerlaw', 'kappa', '2.0', 2.0),
    ('fit_powerlaw', 'kappa', 2.0, 2.0),
    ('fit_powerlaw', 'kappa', 2, 2.0),
    ('fit_powerlaw', 'kappa', ' 2 ', 2.0),
    ('fit_powerlaw', 'kappa', '1e1', 10.0),
    ('fit_powerlaw', 'kappa', 'abc', REFUSED),
    ('fit_powerlaw', 'kappa', '', REFUSED),
    ('fit_powerlaw', 'kappa', 'None', REFUSED),
    ('fit_powerlaw', 'kappa', None, REFUSED),
    ('fit_powerlaw', 'kappa', 'nan', REFUSED),  # 2.4.0: accepted as NaN
    ('fit_powerlaw', 'kappa', 'inf', REFUSED),  # 2.4.0: accepted as an infinity
    ('fit_powerlaw', 'kappa', '-inf', REFUSED),
    ('fit_powerlaw', 'kappa', True, REFUSED),  # 2.4.0: accepted as 1
    ('fit_powerlaw', 'kappa', '-1', REFUSED),
    ('fit_powerlaw', 'kappa', '0', REFUSED),  # 2.4.0: accepted
    ('fit_powerlaw', 'kappa', '3,5', REFUSED),
    ('update_from_knots', 'renorm_model', 'True', True),
    ('update_from_knots', 'renorm_model', 'False', False),
    ('update_from_knots', 'renorm_model', 'true', True),
    ('update_from_knots', 'renorm_model', 'false', False),
    ('update_from_knots', 'renorm_model', 'TRUE', True),
    ('update_from_knots', 'renorm_model', '1', True),
    ('update_from_knots', 'renorm_model', '0', False),
    ('update_from_knots', 'renorm_model', 'yes', REFUSED),  # 2.4.0: accepted as False
    ('update_from_knots', 'renorm_model', True, True),
    ('update_from_knots', 'renorm_model', False, False),
    ('update_from_knots', 'renorm_model', 1, True),
    ('update_from_knots', 'renorm_model', 0, False),
    ('update_from_knots', 'renorm_model', '', REFUSED),  # 2.4.0: accepted as False
    ('update_from_knots', 'renorm_model', None, REFUSED),  # 2.4.0: accepted as False
    ('update_from_knots', 'renorm_model', 'None', REFUSED),  # 2.4.0: accepted as False
    ('calculate_running_std', 'window_pix', '21', 21),
    ('calculate_running_std', 'window_pix', 21, 21),
    ('calculate_running_std', 'window_pix', '21.0', REFUSED),
    ('calculate_running_std', 'window_pix', 21.0, 21),
    ('calculate_running_std', 'window_pix', 21.5, REFUSED),  # 2.4.0: accepted, truncated
    ('calculate_running_std', 'window_pix', '', REFUSED),
    ('calculate_running_std', 'window_pix', 'abc', REFUSED),
    ('calculate_running_std', 'window_pix', True, REFUSED),  # 2.4.0: accepted as 1
    ('calculate_running_std', 'window_pix', False, REFUSED),  # 2.4.0: accepted as 0
    ('calculate_running_std', 'window_pix', '-1', REFUSED),  # 2.4.0: accepted
    ('calculate_running_std', 'window_pix', '0', REFUSED),  # 2.4.0: accepted
    ('calculate_running_std', 'window_pix', '2', 2),
    ('calculate_running_std', 'window_pix', None, REFUSED),  # 2.4.0: crashed (TypeError)
    ('calculate_running_std', 'window_pix', 'None', REFUSED),
    ('smooth', 'sigma_kms', '2.0', 2.0),
    ('smooth', 'sigma_kms', 2.0, 2.0),
    ('smooth', 'sigma_kms', 2, 2.0),
    ('smooth', 'sigma_kms', ' 2 ', 2.0),
    ('smooth', 'sigma_kms', '1e1', 10.0),
    ('smooth', 'sigma_kms', 'abc', REFUSED),
    ('smooth', 'sigma_kms', '', REFUSED),
    ('smooth', 'sigma_kms', 'None', REFUSED),
    ('smooth', 'sigma_kms', None, REFUSED),  # 2.4.0: crashed (TypeError)
    ('smooth', 'sigma_kms', 'nan', REFUSED),  # 2.4.0: accepted as NaN
    ('smooth', 'sigma_kms', 'inf', REFUSED),
    ('smooth', 'sigma_kms', '-inf', REFUSED),  # 2.4.0: accepted as an infinity
    ('smooth', 'sigma_kms', True, REFUSED),  # 2.4.0: accepted as 1
    ('smooth', 'sigma_kms', '-1', REFUSED),  # 2.4.0: accepted
    ('smooth', 'sigma_kms', '0', REFUSED),  # 2.4.0: accepted
    ('smooth', 'sigma_kms', '3,5', REFUSED),
    ('rebin', 'xstart', '2.0', 2.0),
    ('rebin', 'xstart', 2.0, 2.0),
    ('rebin', 'xstart', 2, 2.0),
    ('rebin', 'xstart', ' 2 ', 2.0),
    ('rebin', 'xstart', '1e1', 10.0),
    ('rebin', 'xstart', 'abc', REFUSED),
    ('rebin', 'xstart', '', None),  # 2.4.0: refused (0)
    ('rebin', 'xstart', 'None', None),
    ('rebin', 'xstart', None, None),  # 2.4.0: refused (0)
    ('rebin', 'xstart', 'nan', REFUSED),  # 2.4.0: crashed (ValueError)
    ('rebin', 'xstart', 'inf', REFUSED),  # 2.4.0: crashed (ValueError)
    ('rebin', 'xstart', '-inf', REFUSED),  # 2.4.0: crashed (ValueError)
    ('rebin', 'xstart', True, REFUSED),  # 2.4.0: accepted as 1
    ('rebin', 'xstart', '-1', RAISES),
    ('rebin', 'xstart', '0', RAISES),
    ('rebin', 'xstart', '3,5', REFUSED),
    ('rebin', 'xstart', '400', 400.0),
    ('rebin', 'xstart', 400, 400.0),
    ('rebin', 'xend', '2.0', 2.0),
    ('rebin', 'xend', 2.0, 2.0),
    ('rebin', 'xend', 2, 2.0),
    ('rebin', 'xend', ' 2 ', 2.0),
    ('rebin', 'xend', '1e1', 10.0),
    ('rebin', 'xend', 'abc', REFUSED),
    ('rebin', 'xend', '', None),  # 2.4.0: refused (0)
    ('rebin', 'xend', 'None', None),
    ('rebin', 'xend', None, None),  # 2.4.0: refused (0)
    ('rebin', 'xend', 'nan', REFUSED),  # 2.4.0: crashed (ValueError)
    ('rebin', 'xend', 'inf', REFUSED),  # 2.4.0: crashed (ValueError)
    ('rebin', 'xend', '-inf', REFUSED),  # 2.4.0: crashed (ValueError)
    ('rebin', 'xend', True, REFUSED),  # 2.4.0: accepted as 1
    ('rebin', 'xend', '-1', RAISES),
    ('rebin', 'xend', '0', RAISES),
    ('rebin', 'xend', '3,5', REFUSED),
    ('rebin', 'xend', '440', 440.0),
    ('rebin', 'xend', 440, 440.0),
    ('rebin', 'dx', '20', 20.0),
    ('rebin', 'dx', 20, 20.0),
    ('rebin', 'dx', 20.0, 20.0),
    ('rebin', 'dx', '', REFUSED),
    ('rebin', 'dx', 'None', REFUSED),
    ('rebin', 'dx', None, REFUSED),
    ('rebin', 'dx', 'abc', REFUSED),
    ('rebin', 'dx', '0', REFUSED),  # 2.4.0: crashed (ValueError)
    ('rebin', 'dx', '-5', REFUSED),  # 2.4.0: accepted
    ('rebin', 'dx', True, REFUSED),  # 2.4.0: accepted as 1
    ('rebin', 'dx', 'nan', REFUSED),  # 2.4.0: crashed (ValueError)
    ('rebin', 'kappa', '2.0', 2.0),
    ('rebin', 'kappa', 2.0, 2.0),
    ('rebin', 'kappa', 2, 2.0),
    ('rebin', 'kappa', ' 2 ', 2.0),
    ('rebin', 'kappa', '1e1', 10.0),
    ('rebin', 'kappa', 'abc', REFUSED),
    ('rebin', 'kappa', '', REFUSED),
    ('rebin', 'kappa', 'None', None),
    ('rebin', 'kappa', None, None),  # 2.4.0: refused (0)
    ('rebin', 'kappa', 'nan', REFUSED),  # 2.4.0: accepted as NaN
    ('rebin', 'kappa', 'inf', REFUSED),  # 2.4.0: accepted as an infinity
    ('rebin', 'kappa', '-inf', REFUSED),  # 2.4.0: accepted as an infinity
    ('rebin', 'kappa', True, REFUSED),  # 2.4.0: accepted as 1
    ('rebin', 'kappa', '-1', REFUSED),  # 2.4.0: accepted
    ('rebin', 'kappa', '0', REFUSED),  # 2.4.0: accepted
    ('rebin', 'kappa', '3,5', REFUSED),
    ('rebin', 'filling', '2.0', 2.0),
    ('rebin', 'filling', 2.0, 2.0),
    ('rebin', 'filling', 2, 2.0),
    ('rebin', 'filling', ' 2 ', 2.0),
    ('rebin', 'filling', '1e1', 10.0),
    ('rebin', 'filling', 'abc', REFUSED),
    ('rebin', 'filling', '', REFUSED),
    ('rebin', 'filling', 'None', REFUSED),
    ('rebin', 'filling', None, REFUSED),
    ('rebin', 'filling', 'nan', NAN),
    ('rebin', 'filling', 'inf', REFUSED),  # 2.4.0: accepted as an infinity
    ('rebin', 'filling', '-inf', REFUSED),  # 2.4.0: accepted as an infinity
    ('rebin', 'filling', True, REFUSED),  # 2.4.0: accepted as 1
    ('rebin', 'filling', '-1', -1.0),
    ('rebin', 'filling', '0', 0.0),
    ('rebin', 'filling', '3,5', REFUSED),
    ('rebin', 'xunit', 'km/s', 'km/s'),
    ('rebin', 'xunit', 'nm', 'nm'),
    ('rebin', 'xunit', 'bogus', FAILS),
    ('rebin', 'xunit', '', RAISES),
    ('rebin', 'xunit', None, REFUSED),
    ('rebin', 'xunit', 5, RAISES),
    ('equalize', 'order', '21', 21),
    ('equalize', 'order', 21, 21),
    ('equalize', 'order', '21.0', REFUSED),
    ('equalize', 'order', 21.0, 21),
    ('equalize', 'order', 21.5, REFUSED),  # 2.4.0: accepted, truncated
    ('equalize', 'order', '', REFUSED),
    ('equalize', 'order', 'abc', REFUSED),
    ('equalize', 'order', True, REFUSED),  # 2.4.0: accepted as 1
    ('equalize', 'order', False, REFUSED),  # 2.4.0: accepted as 0
    ('equalize', 'order', '-1', -1),
    ('equalize', 'order', '0', 0),
    ('equalize', 'order', '2', 2),
    ('equalize', 'order', None, REFUSED),
    ('equalize', 'order', 'None', REFUSED),
    ('equalize', 'reference_session', 'other', 'other'),
    ('equalize', 'reference_session', 'None', FAILS),
    ('equalize', 'reference_session', '', FAILS),
    ('equalize', 'reference_session', None, REFUSED),
    ('equalize', 'reference_session', 'missing', FAILS),
    ('resample', 'target_session', 'other', 'other'),
    ('resample', 'target_session', 'None', FAILS),
    ('resample', 'target_session', '', FAILS),
    ('resample', 'target_session', None, REFUSED),
    ('resample', 'target_session', 'missing', FAILS),
]


@pytest.mark.parametrize("name, param, value, expected", CASES,
                         ids=[f"{n}.{p}={v!r}" for n, p, v, _ in CASES])
def test_input(sessions, caplog, name, param, value, expected):
    # The outcome is characterized here, not the warnings of the libraries
    # on absurd values (e.g. the log of a negative kappa)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _check(sessions, caplog, name, param, value, expected)


def _check(sessions, caplog, name, param, value, expected):
    if expected == RAISES:
        with pytest.raises(Exception):
            _run(sessions, name, param, value)
        return
    with caplog.at_level(logging.ERROR):
        result = _run(sessions, name, param, value)
    if expected in (REFUSED, FAILS):
        assert isinstance(result, int) and result == 0
        if expected == REFUSED:
            assert f"{param}: " in caplog.text  # the error names the parameter
        return
    assert isinstance(result, SessionV2), caplog.text
    with caplog.at_level(logging.ERROR):
        reference = _run(sessions, name, param, expected)
    assert isinstance(reference, SessionV2)
    for got, want in zip(_data(result), _data(reference), strict=True):
        assert np.array_equal(got, want, equal_nan=True)
