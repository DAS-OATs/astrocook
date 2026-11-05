"""What the typed recipes of Edit do with each kind of input (phase 6-bis, PR 2b).

As test_recipe_inputs.py for continuum, flux and file: for each value
given to a parameter (with the other parameters in base), the value
the recipe uses (the recipe must give the same result as with that typed
value: same data, name, object, z_em, units and systems), or how it
fails (REFUSED, FAILS, RETURNS_0, RAISES: ``recipe_outcomes.py``; PR 3b
turns the RETURNS_0 into exceptions). The cases whose outcome
changed from 2.4.0 say what 2.4.0 did (the battery was run on master
before the code changed). Guido's decisions on the bounds (F46): a
redshift above -1 and finite; a resolution R or N pixels positive and
finite (0, negative and NaN were "no change", '0px' an infinite R);
widths and thresholds positive and finite.
"""
import hashlib
import logging
import warnings
from types import SimpleNamespace

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.session import SessionV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import (ComponentDataV2, DataColumnV2, SpectrumDataV2,
                                       SystemListDataV2)
from astrocook.core.system_list import SystemListV2

# The outcomes of a failure (phase 6-bis, PR 3): see recipe_outcomes.py
from recipe_outcomes import FAILS, FAILURES, RAISES, REFUSED, RETURNS_0, check_failure  # noqa: E402,F401
NAN = float("nan")
CSV = "<the CSV file of the fixture>"
FLUX = au.erg / au.cm**2 / au.s / au.nm


def _spec(x0=380.0, x1=460.0, scale=1.0):
    rng = np.random.default_rng(0)
    x = np.arange(x0, x1, 0.02)
    cont = scale * (x / 420.0) ** -1.5
    y = (cont * (1 - 0.5 * np.exp(-0.5 * ((x - 420.0) / 0.05) ** 2))
         + 0.02 * cont * rng.standard_normal(len(x)))
    return SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - 0.01, au.nm),
        xmax=DataColumnV2(x + 0.01, au.nm), y=DataColumnV2(y, FLUX),
        dy=DataColumnV2(0.02 * cont, FLUX), aux_cols={"cont": DataColumnV2(cont, FLUX)},
        z_em=2.5, resol=0.0))


def _systs(*zs):
    return SystemListV2(SystemListDataV2(components=[
        ComponentDataV2(id=i + 1, z=z, dz=None, logN=13.5, dlogN=None, b=10.0, db=None,
                        series="CIV") for i, z in enumerate(zs)]))


def _session():
    """The main session (380-460 nm, one CIV system) and three others in the
    GUI: blue (360-400 nm, a system inside the main range), red (440-500
    nm, twice the flux), same (the main range, 1.5 times the flux)."""
    s = SessionV2(name="main", gui=None, spec=_spec(), systs=_systs(1.70))
    others = {"blue": SessionV2(name="blue", gui=None, spec=_spec(360.0, 400.0, 1.0),
                                systs=_systs(1.40, 1.71)),
              "red": SessionV2(name="red", gui=None, spec=_spec(440.0, 500.0, 2.0)),
              "same": SessionV2(name="same", gui=None, spec=_spec(scale=1.5))}
    s._gui = SimpleNamespace(session_histories=[
        SimpleNamespace(display_name=k, current_state=v) for k, v in others.items()])
    return s


@pytest.fixture(scope="module")
def csv_file(tmp_path_factory):
    path = tmp_path_factory.mktemp("systems") / "s.csv"
    path.write_text("series,z,logN,b\nCIV,1.71,13.2,12.0\nCIV,1.72,13.5,8.0\n")
    return str(path)


def _fingerprint(result):
    """The data, name, object, z_em, units and systems of a result (or a
    list of results)."""
    if isinstance(result, list):
        return [_fingerprint(r) for r in result]
    columns = [np.asarray(result.spec.x.value), np.asarray(result.spec.y.value)]
    for c in ("cont", "resol", "dy"):
        if c in result.spec.t.colnames:
            columns.append(np.asarray(result.spec.t[c], dtype=float))
    meta = [result.name, str(result.spec.meta.get("OBJECT")), repr(result.spec._data.z_em),
            str(result.spec.x.unit), str(result.spec.y.unit),
            repr([(c.series, round(c.z, 6)) for c in
                  (result.systs.components if result.systs else [])])]
    digest = hashlib.md5()
    for v in columns:
        digest.update(np.round(np.nan_to_num(v, nan=-9e9), 9).tobytes())
    digest.update("|".join(meta).encode())
    return digest.hexdigest()


def _run(name, param, value, base, csv_file):
    kwargs = {k: (csv_file if v == CSV else v) for k, v in base.items()}
    kwargs[param] = csv_file if value == CSV else value
    if name == "coadd" and param != "dx":
        kwargs.setdefault("dx", 20.0)
    return getattr(_session().edit, name)(**kwargs)


CASES = [
    ('set_properties', 'z_em', '2.0', {}, 2.0),
    ('set_properties', 'z_em', 2.0, {}, 2.0),
    ('set_properties', 'z_em', 2, {}, 2.0),
    ('set_properties', 'z_em', ' 2 ', {}, 2.0),
    ('set_properties', 'z_em', '1e1', {}, 10.0),
    ('set_properties', 'z_em', 'abc', {}, REFUSED),
    ('set_properties', 'z_em', '', {}, REFUSED),
    ('set_properties', 'z_em', 'None', {}, REFUSED),
    ('set_properties', 'z_em', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('set_properties', 'z_em', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('set_properties', 'z_em', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('set_properties', 'z_em', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('set_properties', 'z_em', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('set_properties', 'z_em', '0', {}, 0.0),
    ('set_properties', 'z_em', '_current_', {}, '_current_'),
    ('set_properties', 'z_em', '3.1', {}, 3.1),
    ('set_properties', 'resol', '2.0', {}, '2.0'),
    ('set_properties', 'resol', 2.0, {}, '2.0'),
    ('set_properties', 'resol', 2, {}, '2'),
    ('set_properties', 'resol', ' 2 ', {}, ' 2 '),
    ('set_properties', 'resol', '1e1', {}, '1e1'),
    ('set_properties', 'resol', 'abc', {}, RETURNS_0),
    ('set_properties', 'resol', '', {}, RETURNS_0),
    ('set_properties', 'resol', 'None', {}, RETURNS_0),
    ('set_properties', 'resol', None, {}, REFUSED),
    ('set_properties', 'resol', 'nan', {}, RETURNS_0),  # 2.4.0: accepted as no change
    ('set_properties', 'resol', 'inf', {}, RETURNS_0),  # 2.4.0: accepted as an infinity
    ('set_properties', 'resol', True, {}, RETURNS_0),
    ('set_properties', 'resol', '-1', {}, RETURNS_0),  # 2.4.0: accepted as no change
    ('set_properties', 'resol', '0', {}, RETURNS_0),  # 2.4.0: accepted as no change
    ('set_properties', 'resol', '_current_', {}, '_current_'),
    ('set_properties', 'resol', '45000', {}, '45000'),
    ('set_properties', 'resol', 45000, {}, '45000'),
    ('set_properties', 'resol', '3px', {}, '3px'),
    ('set_properties', 'resol', '3 px', {}, '3 px'),
    ('set_properties', 'resol', 'px', {}, RETURNS_0),
    ('set_properties', 'resol', '-3px', {}, RETURNS_0),
    ('set_properties', 'resol', '0px', {}, RETURNS_0),  # 2.4.0: accepted as an infinite R
    ('set_properties', 'resol', '45000.5', {}, '45000.5'),
    ('set_properties', 'name', 'new', {}, 'new'),
    ('set_properties', 'name', '', {}, ''),
    ('set_properties', 'name', '_current_', {}, '_current_'),
    ('set_properties', 'name', None, {}, REFUSED),  # 2.4.0: accepted
    ('set_properties', 'name', 5, {}, '5'),  # 2.4.0: crashed (TypeError)
    ('set_properties', 'object', 'Q', {}, 'Q'),
    ('set_properties', 'object', '', {}, ''),
    ('set_properties', 'object', '_current_', {}, '_current_'),
    ('set_properties', 'object', None, {}, REFUSED),  # 2.4.0: accepted
    ('set_properties', 'object', 5, {}, '5'),
    ('set_properties', 'overwrite_resol_col', 'True', {'resol': '45000'}, True),
    ('set_properties', 'overwrite_resol_col', 'False', {'resol': '45000'}, False),
    ('set_properties', 'overwrite_resol_col', 'true', {'resol': '45000'}, True),
    ('set_properties', 'overwrite_resol_col', 'false', {'resol': '45000'}, False),
    ('set_properties', 'overwrite_resol_col', '1', {'resol': '45000'}, True),
    ('set_properties', 'overwrite_resol_col', '0', {'resol': '45000'}, False),
    ('set_properties', 'overwrite_resol_col', 'yes', {'resol': '45000'}, REFUSED),  # 2.4.0: accepted as False
    ('set_properties', 'overwrite_resol_col', True, {'resol': '45000'}, True),
    ('set_properties', 'overwrite_resol_col', False, {'resol': '45000'}, False),
    ('set_properties', 'overwrite_resol_col', 1, {'resol': '45000'}, True),
    ('set_properties', 'overwrite_resol_col', 0, {'resol': '45000'}, False),
    ('set_properties', 'overwrite_resol_col', '', {'resol': '45000'}, REFUSED),  # 2.4.0: accepted as False
    ('set_properties', 'overwrite_resol_col', None, {'resol': '45000'}, REFUSED),  # 2.4.0: accepted as False
    ('set_properties', 'overwrite_resol_col', 'None', {'resol': '45000'}, REFUSED),  # 2.4.0: accepted as False
    ('x_convert', 'z_ref', '2.0', {}, 2.0),
    ('x_convert', 'z_ref', 2.0, {}, 2.0),
    ('x_convert', 'z_ref', 2, {}, 2.0),
    ('x_convert', 'z_ref', ' 2 ', {}, 2.0),
    ('x_convert', 'z_ref', '1e1', {}, 10.0),
    ('x_convert', 'z_ref', 'abc', {}, REFUSED),
    ('x_convert', 'z_ref', '', {}, REFUSED),
    ('x_convert', 'z_ref', 'None', {}, REFUSED),
    ('x_convert', 'z_ref', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('x_convert', 'z_ref', 'nan', {}, REFUSED),  # 2.4.0: accepted as NaN
    ('x_convert', 'z_ref', 'inf', {}, REFUSED),  # 2.4.0: accepted as an infinity
    ('x_convert', 'z_ref', True, {}, REFUSED),  # 2.4.0: accepted as 1
    ('x_convert', 'z_ref', '-1', {}, REFUSED),  # 2.4.0: accepted
    ('x_convert', 'z_ref', '0', {}, 0.0),
    ('x_convert', 'z_ref', '_current_', {}, '_current_'),
    ('x_convert', 'z_ref', '2.5', {}, 2.5),
    ('x_convert', 'xunit', 'km/s', {}, 'km/s'),
    ('x_convert', 'xunit', 'nm', {}, 'nm'),
    ('x_convert', 'xunit', 'Angstrom', {}, 'Angstrom'),
    ('x_convert', 'xunit', 'bogus', {}, RETURNS_0),
    ('x_convert', 'xunit', '', {}, FAILS),
    ('x_convert', 'xunit', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('x_convert', 'xunit', 5, {}, FAILS),
    ('y_convert', 'yunit', 'erg/(nm s cm^2)', {}, 'erg/(nm s cm^2)'),
    ('y_convert', 'yunit', 'erg/(Angstrom s cm^2)', {}, 'erg/(Angstrom s cm^2)'),
    ('y_convert', 'yunit', 'Jy', {}, FAILS),
    ('y_convert', 'yunit', 'bogus', {}, RETURNS_0),
    ('y_convert', 'yunit', '', {}, FAILS),
    ('y_convert', 'yunit', None, {}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('y_convert', 'yunit', 5, {}, FAILS),
    ('apply_expression', 'expression', 'F / cont', {'target_col': 'F'}, 'F / cont'),
    ('apply_expression', 'expression', 'F * 2.0', {'target_col': 'F'}, 'F * 2.0'),
    ('apply_expression', 'expression', '', {'target_col': 'F'}, RETURNS_0),
    ('apply_expression', 'expression', '   ', {'target_col': 'F'}, RETURNS_0),
    ('apply_expression', 'expression', 'F +', {'target_col': 'F'}, RETURNS_0),
    ('apply_expression', 'expression', None, {'target_col': 'F'}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('apply_expression', 'expression', 5, {'target_col': 'F'}, '5'),  # 2.4.0: crashed (AttributeError)
    ('apply_expression', 'target_col', 'F', {'expression': 'F * 2'}, 'F'),
    ('apply_expression', 'target_col', 'new', {'expression': 'F * 2'}, 'new'),
    ('apply_expression', 'target_col', '', {'expression': 'F * 2'}, ''),
    ('apply_expression', 'target_col', None, {'expression': 'F * 2'}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('apply_expression', 'target_col', 5, {'expression': 'F * 2'}, '5'),  # 2.4.0: crashed (AttributeError)
    ('mask_expression', 'expression', 'F < 0.5', {'target_col': 'F'}, 'F < 0.5'),
    ('mask_expression', 'expression', '(λ > 400) & (λ < 410)', {'target_col': 'F'}, '(λ > 400) & (λ < 410)'),
    ('mask_expression', 'expression', '', {'target_col': 'F'}, RETURNS_0),
    ('mask_expression', 'expression', 'F <', {'target_col': 'F'}, RETURNS_0),
    ('mask_expression', 'expression', None, {'target_col': 'F'}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('smooth_column', 'sigma_kms', '2.0', {'target_col': 'F'}, 2.0),
    ('smooth_column', 'sigma_kms', 2.0, {'target_col': 'F'}, 2.0),
    ('smooth_column', 'sigma_kms', 2, {'target_col': 'F'}, 2.0),
    ('smooth_column', 'sigma_kms', ' 2 ', {'target_col': 'F'}, 2.0),
    ('smooth_column', 'sigma_kms', '1e1', {'target_col': 'F'}, 10.0),
    ('smooth_column', 'sigma_kms', 'abc', {'target_col': 'F'}, REFUSED),
    ('smooth_column', 'sigma_kms', '', {'target_col': 'F'}, REFUSED),
    ('smooth_column', 'sigma_kms', 'None', {'target_col': 'F'}, REFUSED),
    ('smooth_column', 'sigma_kms', None, {'target_col': 'F'}, REFUSED),  # 2.4.0: crashed (TypeError)
    ('smooth_column', 'sigma_kms', 'nan', {'target_col': 'F'}, REFUSED),  # 2.4.0: accepted as NaN
    ('smooth_column', 'sigma_kms', 'inf', {'target_col': 'F'}, REFUSED),
    ('smooth_column', 'sigma_kms', True, {'target_col': 'F'}, REFUSED),  # 2.4.0: accepted as 1
    ('smooth_column', 'sigma_kms', '-1', {'target_col': 'F'}, REFUSED),  # 2.4.0: accepted
    ('smooth_column', 'sigma_kms', '0', {'target_col': 'F'}, REFUSED),  # 2.4.0: accepted
    ('smooth_column', 'renorm_model', 'True', {'target_col': 'cont'}, True),
    ('smooth_column', 'renorm_model', 'False', {'target_col': 'cont'}, False),
    ('smooth_column', 'renorm_model', 'true', {'target_col': 'cont'}, True),
    ('smooth_column', 'renorm_model', 'false', {'target_col': 'cont'}, False),
    ('smooth_column', 'renorm_model', '1', {'target_col': 'cont'}, True),
    ('smooth_column', 'renorm_model', '0', {'target_col': 'cont'}, False),
    ('smooth_column', 'renorm_model', 'yes', {'target_col': 'cont'}, REFUSED),  # 2.4.0: accepted as False
    ('smooth_column', 'renorm_model', True, {'target_col': 'cont'}, True),
    ('smooth_column', 'renorm_model', False, {'target_col': 'cont'}, False),
    ('smooth_column', 'renorm_model', 1, {'target_col': 'cont'}, True),
    ('smooth_column', 'renorm_model', 0, {'target_col': 'cont'}, False),
    ('smooth_column', 'renorm_model', '', {'target_col': 'cont'}, REFUSED),  # 2.4.0: accepted as False
    ('smooth_column', 'renorm_model', None, {'target_col': 'cont'}, REFUSED),  # 2.4.0: accepted as False
    ('smooth_column', 'renorm_model', 'None', {'target_col': 'cont'}, REFUSED),  # 2.4.0: accepted as False
    ('smooth_column', 'target_col', 'F', {}, 'F'),
    ('smooth_column', 'target_col', 'dF', {}, 'dF'),
    ('smooth_column', 'target_col', 'cont', {}, 'cont'),
    ('smooth_column', 'target_col', 'nope', {}, RETURNS_0),
    ('smooth_column', 'target_col', '', {}, RETURNS_0),
    ('smooth_column', 'target_col', None, {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('import_systems', 'append', 'True', {'source_session': 'blue'}, True),
    ('import_systems', 'append', 'False', {'source_session': 'blue'}, False),
    ('import_systems', 'append', 'true', {'source_session': 'blue'}, True),  # 2.4.0: accepted as False
    ('import_systems', 'append', 'false', {'source_session': 'blue'}, False),
    ('import_systems', 'append', '1', {'source_session': 'blue'}, True),  # 2.4.0: accepted as False
    ('import_systems', 'append', '0', {'source_session': 'blue'}, False),
    ('import_systems', 'append', 'yes', {'source_session': 'blue'}, REFUSED),  # 2.4.0: accepted as False
    ('import_systems', 'append', True, {'source_session': 'blue'}, True),
    ('import_systems', 'append', False, {'source_session': 'blue'}, False),
    ('import_systems', 'append', 1, {'source_session': 'blue'}, True),  # 2.4.0: accepted as False
    ('import_systems', 'append', 0, {'source_session': 'blue'}, False),
    ('import_systems', 'append', '', {'source_session': 'blue'}, REFUSED),  # 2.4.0: accepted as False
    ('import_systems', 'append', None, {'source_session': 'blue'}, REFUSED),  # 2.4.0: accepted as False
    ('import_systems', 'append', 'None', {'source_session': 'blue'}, REFUSED),  # 2.4.0: accepted as False
    ('import_systems', 'refit', 'False', {'source_session': 'blue'}, False),
    ('import_systems', 'refit', False, {'source_session': 'blue'}, False),
    ('import_systems', 'refit', 'false', {'source_session': 'blue'}, False),
    ('import_systems', 'refit', '0', {'source_session': 'blue'}, False),
    ('import_systems', 'refit', 0, {'source_session': 'blue'}, False),
    ('import_systems', 'refit', 'no', {'source_session': 'blue'}, REFUSED),  # 2.4.0: accepted as False
    ('import_systems', 'refit', '', {'source_session': 'blue'}, REFUSED),  # 2.4.0: accepted as False
    ('import_systems', 'refit', None, {'source_session': 'blue'}, REFUSED),  # 2.4.0: accepted as False
    ('import_systems', 'source_session', 'blue', {}, 'blue'),
    ('import_systems', 'source_session', 'red', {}, 'red'),
    ('import_systems', 'source_session', 'missing', {}, RETURNS_0),
    ('import_systems', 'source_session', 'None', {}, RETURNS_0),
    ('import_systems', 'source_session', '', {}, RETURNS_0),
    ('import_systems', 'source_session', None, {}, REFUSED),
    ('import_systems_ascii', 'append', 'True', {'file_path': CSV}, True),
    ('import_systems_ascii', 'append', 'False', {'file_path': CSV}, False),
    ('import_systems_ascii', 'append', 'true', {'file_path': CSV}, True),  # 2.4.0: accepted as False
    ('import_systems_ascii', 'append', 'false', {'file_path': CSV}, False),
    ('import_systems_ascii', 'append', '1', {'file_path': CSV}, True),  # 2.4.0: accepted as False
    ('import_systems_ascii', 'append', '0', {'file_path': CSV}, False),
    ('import_systems_ascii', 'append', 'yes', {'file_path': CSV}, REFUSED),  # 2.4.0: accepted as False
    ('import_systems_ascii', 'append', True, {'file_path': CSV}, True),
    ('import_systems_ascii', 'append', False, {'file_path': CSV}, False),
    ('import_systems_ascii', 'append', 1, {'file_path': CSV}, True),  # 2.4.0: accepted as False
    ('import_systems_ascii', 'append', 0, {'file_path': CSV}, False),
    ('import_systems_ascii', 'append', '', {'file_path': CSV}, REFUSED),  # 2.4.0: accepted as False
    ('import_systems_ascii', 'append', None, {'file_path': CSV}, REFUSED),  # 2.4.0: accepted as False
    ('import_systems_ascii', 'append', 'None', {'file_path': CSV}, REFUSED),  # 2.4.0: accepted as False
    ('import_systems_ascii', 'file_path', CSV, {}, CSV),
    ('import_systems_ascii', 'file_path', '/nonexistent.csv', {}, RETURNS_0),
    ('import_systems_ascii', 'file_path', '', {}, RETURNS_0),
    ('import_systems_ascii', 'file_path', None, {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('delete', 'targets', 'cont', {}, 'cont'),
    ('delete', 'targets', 'systems', {}, 'systems'),
    ('delete', 'targets', 'cont, systems', {}, 'cont, systems'),
    ('delete', 'targets', '', {}, ''),
    ('delete', 'targets', 'nope', {}, 'nope'),
    ('delete', 'targets', None, {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('split', 'expression', '(λ > 400) & (λ < 410)', {}, '(λ > 400) & (λ < 410)'),
    ('split', 'expression', 'λ > 1000', {}, FAILS),
    ('split', 'expression', '', {}, RETURNS_0),
    ('split', 'expression', None, {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('split', 'expression', 'λ >', {}, RAISES),
    ('extract_preset', 'region', 'lya_forest', {}, 'lya_forest'),
    ('extract_preset', 'region', 'red_side', {}, 'red_side'),
    ('extract_preset', 'region', 'lya_forest, red_side', {}, 'lya_forest, red_side'),
    ('extract_preset', 'region', 'nope', {}, RETURNS_0),
    ('extract_preset', 'region', '', {}, RETURNS_0),
    ('extract_preset', 'region', None, {}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('trim_common', 'z_target', '2.0', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, RETURNS_0),
    ('trim_common', 'z_target', 2.0, {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, RETURNS_0),
    ('trim_common', 'z_target', 2, {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, RETURNS_0),
    ('trim_common', 'z_target', ' 2 ', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, RETURNS_0),
    ('trim_common', 'z_target', '1e1', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, RETURNS_0),
    ('trim_common', 'z_target', 'abc', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'z_target', '', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'z_target', 'None', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'z_target', None, {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'z_target', 'nan', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'z_target', 'inf', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'z_target', True, {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'z_target', '-1', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'z_target', '0', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, RETURNS_0),
    ('trim_common', 'z_target', '_current_', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, RETURNS_0),
    ('trim_common', 'z_target', '2.45', {'others_names': 'same', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, 2.45),
    ('trim_common', 'window_kms', '2.0', {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, RETURNS_0),
    ('trim_common', 'window_kms', 2.0, {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, RETURNS_0),
    ('trim_common', 'window_kms', 2, {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, RETURNS_0),
    ('trim_common', 'window_kms', ' 2 ', {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, RETURNS_0),
    ('trim_common', 'window_kms', '1e1', {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, 10.0),
    ('trim_common', 'window_kms', 'abc', {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'window_kms', '', {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'window_kms', 'None', {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'window_kms', None, {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'window_kms', 'nan', {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),  # 2.4.0: accepted as NaN
    ('trim_common', 'window_kms', 'inf', {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),  # 2.4.0: accepted as an infinity
    ('trim_common', 'window_kms', True, {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'window_kms', '-1', {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('trim_common', 'window_kms', '0', {'others_names': 'same', 'z_target': '2.45', 'trans_self': 'Ly_a', 'trans_others': 'Ly_a'}, REFUSED),
    ('equalize_and_stitch', 'stitch_wavelengths', '398', {'other_sessions': 'blue'}, '398'),
    ('equalize_and_stitch', 'stitch_wavelengths', 398, {'other_sessions': 'blue'}, '398'),  # 2.4.0: crashed (AttributeError)
    ('equalize_and_stitch', 'stitch_wavelengths', 398.0, {'other_sessions': 'blue'}, '398.0'),  # 2.4.0: crashed (AttributeError)
    ('equalize_and_stitch', 'stitch_wavelengths', '', {'other_sessions': 'blue'}, RETURNS_0),
    ('equalize_and_stitch', 'stitch_wavelengths', 'abc', {'other_sessions': 'blue'}, RETURNS_0),
    ('equalize_and_stitch', 'stitch_wavelengths', '398, 450', {'other_sessions': 'blue'}, RETURNS_0),
    ('equalize_and_stitch', 'stitch_wavelengths', None, {'other_sessions': 'blue'}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('equalize_and_stitch', 'manual_factors', 'auto', {'other_sessions': 'blue', 'stitch_wavelengths': '398'}, 'auto'),
    ('equalize_and_stitch', 'manual_factors', 'AUTO', {'other_sessions': 'blue', 'stitch_wavelengths': '398'}, 'AUTO'),
    ('equalize_and_stitch', 'manual_factors', '1.5', {'other_sessions': 'blue', 'stitch_wavelengths': '398'}, '1.5'),
    ('equalize_and_stitch', 'manual_factors', 1.5, {'other_sessions': 'blue', 'stitch_wavelengths': '398'}, '1.5'),  # 2.4.0: crashed (AttributeError)
    ('equalize_and_stitch', 'manual_factors', 'abc', {'other_sessions': 'blue', 'stitch_wavelengths': '398'}, RETURNS_0),
    ('equalize_and_stitch', 'manual_factors', '', {'other_sessions': 'blue', 'stitch_wavelengths': '398'}, RETURNS_0),
    ('equalize_and_stitch', 'manual_factors', None, {'other_sessions': 'blue', 'stitch_wavelengths': '398'}, REFUSED),  # 2.4.0: crashed (AttributeError)
    ('coadd', 'xstart', '2.0', {'session_names': 'main2, same'}, 2.0),
    ('coadd', 'xstart', 2.0, {'session_names': 'main2, same'}, 2.0),
    ('coadd', 'xstart', 2, {'session_names': 'main2, same'}, 2.0),
    ('coadd', 'xstart', ' 2 ', {'session_names': 'main2, same'}, 2.0),
    ('coadd', 'xstart', '1e1', {'session_names': 'main2, same'}, 10.0),
    ('coadd', 'xstart', 'abc', {'session_names': 'main2, same'}, REFUSED),
    ('coadd', 'xstart', '', {'session_names': 'main2, same'}, None),  # 2.4.0: refused (0)
    ('coadd', 'xstart', 'None', {'session_names': 'main2, same'}, None),
    ('coadd', 'xstart', None, {'session_names': 'main2, same'}, None),  # 2.4.0: refused (0)
    ('coadd', 'xstart', 'nan', {'session_names': 'main2, same'}, REFUSED),
    ('coadd', 'xstart', 'inf', {'session_names': 'main2, same'}, REFUSED),
    ('coadd', 'xstart', True, {'session_names': 'main2, same'}, REFUSED),  # 2.4.0: accepted as 1
    ('coadd', 'xstart', '-1', {'session_names': 'main2, same'}, RETURNS_0),
    ('coadd', 'xstart', '0', {'session_names': 'main2, same'}, RETURNS_0),
    ('coadd', 'xstart', '390', {'session_names': 'main2, same'}, 390.0),
    ('coadd', 'xstart', 390, {'session_names': 'main2, same'}, 390.0),
    ('coadd', 'xend', '2.0', {'session_names': 'same'}, 2.0),
    ('coadd', 'xend', 2.0, {'session_names': 'same'}, 2.0),
    ('coadd', 'xend', 2, {'session_names': 'same'}, 2.0),
    ('coadd', 'xend', ' 2 ', {'session_names': 'same'}, 2.0),
    ('coadd', 'xend', '1e1', {'session_names': 'same'}, 10.0),
    ('coadd', 'xend', 'abc', {'session_names': 'same'}, REFUSED),
    ('coadd', 'xend', '', {'session_names': 'same'}, None),  # 2.4.0: refused (0)
    ('coadd', 'xend', 'None', {'session_names': 'same'}, None),
    ('coadd', 'xend', None, {'session_names': 'same'}, None),  # 2.4.0: refused (0)
    ('coadd', 'xend', 'nan', {'session_names': 'same'}, REFUSED),
    ('coadd', 'xend', 'inf', {'session_names': 'same'}, REFUSED),
    ('coadd', 'xend', True, {'session_names': 'same'}, REFUSED),  # 2.4.0: accepted as 1
    ('coadd', 'xend', '-1', {'session_names': 'same'}, RETURNS_0),
    ('coadd', 'xend', '0', {'session_names': 'same'}, RETURNS_0),
    ('coadd', 'xend', '450', {'session_names': 'same'}, 450.0),
    ('coadd', 'xend', 450, {'session_names': 'same'}, 450.0),
    ('coadd', 'dx', '20', {'session_names': 'same'}, 20.0),
    ('coadd', 'dx', 20, {'session_names': 'same'}, 20.0),
    ('coadd', 'dx', 20.0, {'session_names': 'same'}, 20.0),
    ('coadd', 'dx', '', {'session_names': 'same'}, REFUSED),
    ('coadd', 'dx', 'None', {'session_names': 'same'}, REFUSED),
    ('coadd', 'dx', None, {'session_names': 'same'}, REFUSED),
    ('coadd', 'dx', 'abc', {'session_names': 'same'}, REFUSED),
    ('coadd', 'dx', 'nan', {'session_names': 'same'}, REFUSED),
    ('coadd', 'dx', True, {'session_names': 'same'}, REFUSED),  # 2.4.0: accepted as 1
    ('coadd', 'kappa', '2.0', {'session_names': 'same'}, 2.0),
    ('coadd', 'kappa', 2.0, {'session_names': 'same'}, 2.0),
    ('coadd', 'kappa', 2, {'session_names': 'same'}, 2.0),
    ('coadd', 'kappa', ' 2 ', {'session_names': 'same'}, 2.0),
    ('coadd', 'kappa', '1e1', {'session_names': 'same'}, 10.0),
    ('coadd', 'kappa', 'abc', {'session_names': 'same'}, REFUSED),
    ('coadd', 'kappa', '', {'session_names': 'same'}, REFUSED),
    ('coadd', 'kappa', 'None', {'session_names': 'same'}, None),
    ('coadd', 'kappa', None, {'session_names': 'same'}, None),  # 2.4.0: refused (0)
    ('coadd', 'kappa', 'nan', {'session_names': 'same'}, REFUSED),  # 2.4.0: accepted as NaN
    ('coadd', 'kappa', 'inf', {'session_names': 'same'}, REFUSED),  # 2.4.0: accepted as an infinity
    ('coadd', 'kappa', True, {'session_names': 'same'}, REFUSED),  # 2.4.0: accepted as 1
    ('coadd', 'kappa', '-1', {'session_names': 'same'}, REFUSED),  # 2.4.0: accepted
    ('coadd', 'kappa', '0', {'session_names': 'same'}, REFUSED),  # 2.4.0: accepted
    ('coadd', 'equalize_order', '0', {'session_names': 'same, blue'}, 0),
    ('coadd', 'equalize_order', 0, {'session_names': 'same, blue'}, 0),
    ('coadd', 'equalize_order', '1', {'session_names': 'same, blue'}, 1),
    ('coadd', 'equalize_order', 1, {'session_names': 'same, blue'}, 1),
    ('coadd', 'equalize_order', '-1', {'session_names': 'same, blue'}, -1),
    ('coadd', 'equalize_order', -1, {'session_names': 'same, blue'}, -1),
    ('coadd', 'equalize_order', '0.0', {'session_names': 'same, blue'}, REFUSED),
    ('coadd', 'equalize_order', 0.0, {'session_names': 'same, blue'}, 0),
    ('coadd', 'equalize_order', 0.5, {'session_names': 'same, blue'}, REFUSED),  # 2.4.0: accepted, truncated
    ('coadd', 'equalize_order', '', {'session_names': 'same, blue'}, REFUSED),
    ('coadd', 'equalize_order', 'abc', {'session_names': 'same, blue'}, REFUSED),
    ('coadd', 'equalize_order', True, {'session_names': 'same, blue'}, REFUSED),  # 2.4.0: accepted as 1
    ('coadd', 'equalize_order', None, {'session_names': 'same, blue'}, REFUSED),
    ('coadd', 'xunit', 'km/s', {'session_names': 'same'}, 'km/s'),
    ('coadd', 'xunit', 'nm', {'session_names': 'same'}, 'nm'),
    ('coadd', 'xunit', 'bogus', {'session_names': 'same'}, RETURNS_0),
    ('coadd', 'xunit', '', {'session_names': 'same'}, RETURNS_0),
    ('coadd', 'xunit', None, {'session_names': 'same'}, REFUSED),
]


@pytest.mark.parametrize("name, param, value, base, expected", CASES,
                         ids=[f"{n}.{p}={v!r}" for n, p, v, _, _ in CASES])
def test_input(caplog, csv_file, name, param, value, base, expected):
    # The outcome is characterized here, not the warnings of the libraries
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _check(caplog, csv_file, name, param, value, base, expected)


def _check(caplog, csv_file, name, param, value, base, expected):
    if expected in FAILURES:
        check_failure(lambda: _run(name, param, value, base, csv_file), name, param,
                      expected)
        return
    with caplog.at_level(logging.ERROR):
        result = _run(name, param, value, base, csv_file)
    assert isinstance(result, (SessionV2, list)), caplog.text
    reference = _run(name, param, expected, base, csv_file)
    assert _fingerprint(result) == _fingerprint(reference)
