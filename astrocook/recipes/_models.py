"""The external models that ``import_model`` brings into a spectrum (phase
7-ter, PR T2).

A model is a column of the spectrum with a role (agreed with Guido, 20
November 2026): ``import_model`` reads it from a file or from another
session, brings it onto the grid of the spectrum
(:func:`astrocook.core.model_grid.model_on_grid`: a copy on the same grid,
otherwise an average over the pixels) and writes it in the column of its
role. What changes with the role is in a :class:`ModelRole` of
:data:`MODEL_ROLES`: the column written, the names by which the column of
the model is found in a source, the value beyond the model, and the check
(and conversion) of its values. T2 has the role ``telluric``, T3 ``cont``
(a continuum in the unit of the flux, NaN beyond the model; the ``model``
of the spectrum follows it if the user wants, as in ``fit_continuum``), T4
``model`` (an absorption model, in the unit of the flux or normalized and
multiplied by the normalization ``cont × telluric`` at the import; the
normalization beyond it).
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, Optional

import astropy.units as au
import numpy as np

from astrocook.core.column_names import internal_column
from astrocook.core.model_grid import model_on_grid
from astrocook.core.telluric import (TELLURIC_COLUMN, TELLURIC_NAMES,
                                     is_telluric_corrected)
from astrocook.core.units import known_unit
from astrocook.recipes.errors import InvalidParameter, RecipeFailed

if TYPE_CHECKING:
    from astrocook.core.spectrum import Spectrum

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ModelSource:
    """The column of a model as read from its source, on its own grid."""
    #: Centres and edges of its pixels, as quantities
    x: Any
    xmin: Any
    xmax: Any
    #: Its values, and their unit (None if the source says none)
    values: np.ndarray
    unit: Optional[au.UnitBase]
    #: The name of the column in the source, and the source (for the log)
    column: str
    where: str


@dataclass(frozen=True)
class ModelRole:
    """What ``import_model`` does with a model of one role."""
    #: The column of the spectrum written
    column: str
    #: The names of the column of the model in a source, in order of
    #: preference (matched ignoring the case)
    names: tuple
    #: The value beyond the model (and where it does not cover a pixel)
    fill: float
    #: Words for the messages ("telluric model")
    label: str
    #: The values to write, in their unit, from the source and the spectrum;
    #: InvalidParameter or RecipeFailed if they cannot be this model
    prepare: Callable[[ModelSource, 'Spectrum'], tuple]
    #: What the role says after the import, if anything (e.g. that the flux
    #: is flagged as corrected), at WARNING
    after: Optional[Callable[['Spectrum'], Optional[str]]] = None
    #: A check of the values on the grid of the spectrum (the values, the
    #: pixels covered, the source, the spectrum); InvalidParameter or
    #: RecipeFailed if they cannot be this model
    check: Optional[Callable[[np.ndarray, np.ndarray, ModelSource, 'Spectrum'], None]] = None
    #: Whether the ``model`` of the spectrum may follow the new column (a
    #: continuum: ``renorm_model``, as in ``fit_continuum``)
    renormalizes_model: bool = False
    #: Whether the model may be given normalized (``normalized=True``: 0-1,
    #: multiplied by the normalization of the spectrum at the import; T4)
    normalizable: bool = False
    #: Whether the normalization of the spectrum (``cont × telluric``) is
    #: the value beyond the model, instead of ``fill`` (an absorption
    #: model: no absorption there; NaN where the spectrum has none; T4)
    beyond_is_norm: bool = False


#: How far a transmission may go beyond 0-1 without a warning (T2: a
#: Molecfit model of the request reaches 1 plus a little, from the numerics
#: of the model, and the warning said "from 0.0237 to 1"). Confirmed by
#: Guido (20 November 2026)
TRANSMISSION_TOLERANCE = 1e-3
#: A continuum without a unit is taken in the unit of the flux; it is
#: refused when its median and that of the flux, over the pixels it
#: covers, differ by more than this factor (T3, Guido, 20 November 2026: a
#: normalized continuum, ~1, under a flux of 1e-16, is not a continuum of
#: that flux)
CONTINUUM_SCALE_FACTOR = 10.0


def _prepare_telluric(source: ModelSource, spec: 'Spectrum'):
    """A transmission: dimensionless, between 0 and 1."""
    unit = source.unit
    if unit is not None and unit != au.dimensionless_unscaled:
        raise RecipeFailed(f"The column {source.column} of {source.where} is in {unit}: "
                           "a telluric transmission is dimensionless.")
    finite = source.values[np.isfinite(source.values)]
    if len(finite) and (finite.min() < -TRANSMISSION_TOLERANCE
                        or finite.max() > 1.0 + TRANSMISSION_TOLERANCE):
        logger.warning(f"The column {source.column} of {source.where} goes from "
                       f"{finite.min():.6g} to {finite.max():.6g}: a transmission is "
                       "between 0 and 1. Is it the column of the model?")
    return source.values, au.dimensionless_unscaled


def _after_telluric(spec: 'Spectrum') -> Optional[str]:
    if is_telluric_corrected(spec.meta):
        return ("The flux is flagged as already corrected for the tellurics: the "
                "transmission is kept but does not enter the normalization. "
                "Set Properties with telluric_corrected='no' to use it.")
    return None


#: The names of the column of a continuum in a source
CONTINUUM_NAMES = ('cont', 'continuum')


def _prepare_in_flux(source: ModelSource, spec: 'Spectrum', what: str):
    """A column in the unit of the flux (a continuum, an absorption model):
    a unit of the source is converted to that of the flux (a flux density
    per wavelength or per frequency, at the wavelengths of the model) and
    refused if it cannot be (e.g. a physical flux and ``arb_flux``);
    without a unit the values are taken in the unit of the flux (checked on
    the grid, :func:`_check_cont`, :func:`_check_model`)."""
    y_unit = spec.y.unit
    unit = source.unit
    if unit is None or unit == au.dimensionless_unscaled:
        return source.values, y_unit
    try:
        values = (source.values * unit).to(
            y_unit, equivalencies=au.spectral_density(source.x)).value
    except au.UnitConversionError as e:
        raise RecipeFailed(f"The column {source.column} of {source.where} is in {unit}, "
                           f"the flux of the spectrum in {y_unit}: {what} is in the "
                           "unit of the flux.") from e
    if unit != y_unit:
        logger.info(f"The column {source.column} of {source.where} was converted from "
                    f"{unit} to {y_unit}.")
    return values, y_unit


def _prepare_cont(source: ModelSource, spec: 'Spectrum'):
    """A continuum: in the unit of the flux (:func:`_prepare_in_flux`)."""
    return _prepare_in_flux(source, spec, "a continuum")


def _scale_medians(values: np.ndarray, covered: np.ndarray, source: ModelSource,
                   spec: 'Spectrum') -> Optional[tuple]:
    """The medians of a model without a unit and of the flux over the
    pixels the model covers, when they differ by more than
    :data:`CONTINUUM_SCALE_FACTOR`; None otherwise, and when the source
    gives a unit (it was converted) or a median is not positive."""
    if source.unit is not None and source.unit != au.dimensionless_unscaled:
        return None
    y = np.asarray(spec.y.value, dtype=float)
    ok = covered & np.isfinite(values) & np.isfinite(y)
    if not np.any(ok):
        return None
    m_model, m_flux = float(np.median(values[ok])), float(np.median(y[ok]))
    if m_flux <= 0 or m_model <= 0:
        return None
    if max(m_model / m_flux, m_flux / m_model) > CONTINUUM_SCALE_FACTOR:
        return m_model, m_flux
    return None


def _check_cont(values: np.ndarray, covered: np.ndarray, source: ModelSource,
                spec: 'Spectrum') -> None:
    """A continuum without a unit must be on the scale of the flux: its
    median and that of the flux, over the pixels it covers, within
    :data:`CONTINUUM_SCALE_FACTOR` (a normalized continuum under a
    physical flux is refused). Not checked when the source gives a unit (it
    was converted) or the median of the flux is not positive."""
    medians = _scale_medians(values, covered, source, spec)
    if medians is not None:
        raise RecipeFailed(
            f"The column {source.column} of {source.where} has no unit and a median of "
            f"{medians[0]:.4g}, the flux of the spectrum a median of {medians[1]:.4g} "
            f"{spec.y.unit} over the same pixels: a continuum is in the unit of the "
            "flux (a normalized continuum is not one). Give the continuum in the unit "
            "of the flux, or with its unit in the file.")


def _after_cont(spec: 'Spectrum') -> Optional[str]:
    cont = spec.get_column('cont')
    if cont is None:
        return None
    missing = int(np.count_nonzero(~np.isfinite(np.asarray(cont.value, dtype=float))))
    if missing:
        return (f"{missing} of {len(cont)} pixels have no continuum (beyond the model): "
                "there the flux is not normalized.")
    return None


#: The names of the column of an absorption model in a source
MODEL_NAMES = ('model',)


def _check_model(values: np.ndarray, covered: np.ndarray, source: ModelSource,
                 spec: 'Spectrum') -> None:
    """An absorption model in the unit of the flux (``normalized=False``),
    without a unit: on the scale of the flux, as a continuum
    (:func:`_check_cont`); a normalized model is refused with a pointer to
    ``normalized=True``."""
    medians = _scale_medians(values, covered, source, spec)
    if medians is not None:
        raise RecipeFailed(
            f"The column {source.column} of {source.where} has no unit and a median of "
            f"{medians[0]:.4g}, the flux of the spectrum a median of {medians[1]:.4g} "
            f"{spec.y.unit} over the same pixels: a model is in the unit of the flux. "
            "If the model is normalized (0-1), import it with normalized=True; "
            "otherwise give it in the unit of the flux, or with its unit in the file.")


def _prepare_model(source: ModelSource, spec: 'Spectrum'):
    """An absorption model in the unit of the flux (:func:`_prepare_in_flux`)."""
    return _prepare_in_flux(source, spec, "a model")


def _prepare_normalized(source: ModelSource) -> tuple:
    """A normalized model (T4, ``normalized=True``): dimensionless, between
    0 and 1 (as a transmission: values beyond, by more than
    :data:`TRANSMISSION_TOLERANCE`, a WARNING)."""
    unit = source.unit
    if unit is not None and unit != au.dimensionless_unscaled:
        raise RecipeFailed(f"The column {source.column} of {source.where} is in {unit}: "
                           "a normalized model is dimensionless (import it with "
                           "normalized=False to convert it to the unit of the flux).")
    finite = source.values[np.isfinite(source.values)]
    if len(finite) and (finite.min() < -TRANSMISSION_TOLERANCE
                        or finite.max() > 1.0 + TRANSMISSION_TOLERANCE):
        logger.warning(f"The column {source.column} of {source.where} goes from "
                       f"{finite.min():.6g} to {finite.max():.6g}: a normalized model "
                       "is between 0 and 1. Is it normalized?")
    return source.values, au.dimensionless_unscaled


def _after_model(spec: 'Spectrum') -> Optional[str]:
    model = spec.get_column('model')
    if model is None:
        return None
    missing = int(np.count_nonzero(~np.isfinite(np.asarray(model.value, dtype=float))))
    if missing:
        return (f"{missing} of {len(model)} pixels have no model (beyond the model "
                "the spectrum has no normalization, cont × telluric).")
    return None


#: The roles of the models, by the value of the parameter ``role``
MODEL_ROLES = {
    'telluric': ModelRole(column=TELLURIC_COLUMN, names=TELLURIC_NAMES,
                          fill=1.0, label="telluric model",
                          prepare=_prepare_telluric, after=_after_telluric),
    'cont': ModelRole(column='cont', names=CONTINUUM_NAMES, fill=np.nan,
                      label="continuum", prepare=_prepare_cont, after=_after_cont,
                      check=_check_cont, renormalizes_model=True),
    # T4: a model in the unit of the flux is prepared and checked as a
    # continuum (converted from a unit; on the scale of the flux without
    # one); a normalized one by _prepare_normalized, without the check
    'model': ModelRole(column='model', names=MODEL_NAMES, fill=np.nan,
                       label="absorption model", prepare=_prepare_model,
                       after=_after_model, check=_check_model,
                       normalizable=True, beyond_is_norm=True),
}


# --- The sources ------------------------------------------------------------------

def _not_found(role: ModelRole, column: str, available, where: str) -> RecipeFailed:
    asked = (f"no column {column}" if column else
             f"no column of a {role.label} ({', '.join(role.names)})")
    return RecipeFailed(f"{where} has {asked}: its columns are "
                        f"{', '.join(available)}. Name the column of the model.")


def source_from_file(path: str, role: ModelRole, column: str) -> ModelSource:
    """The model in a file (a table read by name, :mod:`astrocook.io.tables`)."""
    from astrocook.io import tables
    where = os.path.basename(path)
    if not os.path.exists(path):
        raise RecipeFailed(f"{path} does not exist.")
    try:
        table = tables.read_table(path, wanted=(column,) if column else role.names)
        x, xmin, xmax = tables.axis_of(table)
    except (ValueError, OSError) as e:
        raise RecipeFailed(f"{where} cannot be read as a table: {e}") from e
    if column:
        name = tables.find_column(table, [column])
    else:
        axis = {tables.find_column(table, n) for n in
                (tables.AXIS_NAMES, tables.XMIN_NAMES, tables.XMAX_NAMES)}
        name = tables.find_column(table, role.names, exclude=[a for a in axis if a])
    if name is None:
        raise _not_found(role, column, table.colnames, where)
    try:
        values = np.asarray(table[name], dtype=float)
    except (TypeError, ValueError) as e:
        raise RecipeFailed(f"The column {name} of {where} is not numeric.") from e
    unit = known_unit(table[name].unit)
    return ModelSource(x, xmin, xmax, values, unit, name, where)


def source_from_session(other, role: ModelRole, column: str) -> ModelSource:
    """The model in another session: a column of its spectrum."""
    spec = other.spec
    where = f"the session {other.name}"
    if spec is None:
        raise RecipeFailed(f"{where} has no spectrum.")
    cols = spec.t.colnames
    if column:
        name = internal_column(column)
        if name not in cols:
            raise _not_found(role, column, cols, where)
    else:
        lower = {c.lower(): c for c in cols}
        found = next((lower[n.lower()] for n in role.names if n.lower() in lower), None)
        if found is None:
            raise _not_found(role, column, cols, where)
        name = found
    q = spec.get_column(name) if name not in ('x', 'xmin', 'xmax', 'y', 'dy') \
        else getattr(spec, name)
    unit = getattr(q, 'unit', None)
    return ModelSource(spec.x, spec.xmin, spec.xmax,
                       np.asarray(getattr(q, 'value', q), dtype=float),
                       None if unit is None else unit, name, where)


# --- The import ---------------------------------------------------------------------

def _range(lo, hi, unit) -> str:
    return f"{lo:.6g}-{hi:.6g} {unit}"


def _renormalized_model(old: 'Spectrum', new: 'Spectrum') -> Optional['Spectrum']:
    """``new`` with its ``model`` following the new continuum, as
    ``fit_continuum`` does with ``renorm_model`` (the model divided by the
    old continuum, times the new one; 1 where the old continuum is 0);
    where either continuum is not finite (beyond an imported continuum)
    the model is kept as it was. None if there is nothing to follow (no
    model, or no old continuum)."""
    old_cont, new_cont = old.get_column('cont'), new.get_column('cont')
    model = old.get_column('model')
    if old_cont is None or new_cont is None or model is None:
        return None
    try:
        old_vals = old_cont.to(new_cont.unit).value
    except au.UnitConversionError:
        old_vals = np.asarray(old_cont.value, dtype=float)
    new_vals = np.asarray(new_cont.value, dtype=float)
    m = np.asarray(model.value, dtype=float)
    with np.errstate(invalid='ignore'):
        profile = np.divide(m, old_vals, out=np.ones_like(m), where=old_vals != 0)
    rescaled = new_vals * profile
    follow = np.isfinite(rescaled)
    values = np.where(follow, rescaled, m)
    kept = len(m) - int(np.count_nonzero(follow))
    logger.info("The model was renormalized to the new continuum"
                + (f" (kept as it was on {kept} pixels without a continuum)." if kept
                   else "."))
    description = old._data.aux_cols['model'].description
    return new.update_column('model', values, unit=model.unit, description=description)


def _normalization(spec: 'Spectrum') -> Optional[np.ndarray]:
    """The normalization of the spectrum (``cont × telluric``, 1 for a flux
    already normalized without a continuum: :attr:`Spectrum.norm`) in the
    unit of the flux, or None if it has none."""
    norm = spec.norm
    if norm is None:
        return None
    return np.asarray(au.Quantity(norm).to_value(spec.y.unit), dtype=float)


def import_model(spec: 'Spectrum', source: ModelSource, role_name: str,
                 renorm_model: bool = True, normalized: bool = False) -> 'Spectrum':
    """The spectrum with the model in the column of its role; for a
    continuum, the ``model`` of the spectrum follows it if
    ``renorm_model``; for an absorption model, ``normalized`` says that it
    is 0-1, multiplied by the normalization of the spectrum."""
    role = MODEL_ROLES.get(role_name)
    if role is None:
        raise InvalidParameter(f"role: {role_name!r} is not one of "
                               f"{', '.join(repr(r) for r in MODEL_ROLES)}")
    if normalized and not role.normalizable:
        raise InvalidParameter(f"normalized: only a model of the role 'model' can be "
                               f"normalized, not one of the role {role_name!r}.")
    unit = spec.x.unit
    try:
        m_x, m_xmin, m_xmax = (q.to(unit, equivalencies=au.spectral()).value
                               for q in (source.x, source.xmin, source.xmax))
    except au.UnitConversionError as e:
        raise RecipeFailed(f"The axis of {source.where} is in {source.x.unit}, that of "
                           f"the spectrum in {unit}: they cannot be compared.") from e
    if len(m_x) != len(source.values):
        raise RecipeFailed(f"{source.where}: the axis and the column {source.column} "
                           "have different lengths.")
    norm = _normalization(spec) if role.beyond_is_norm or normalized else None
    if normalized:
        if norm is None:
            raise RecipeFailed(
                f"The {role.label} of {source.where} is normalized: it is multiplied by "
                "the normalization of the spectrum (cont × telluric), and the spectrum "
                "has no continuum. Fit or import a continuum first.")
        values, _ = _prepare_normalized(source)
        # Normalized, no absorption beyond the model: 1 (a pixel covered in
        # part mixed with it, as a transmission), then times the normalization
        fill = 1.0
    else:
        values, values_unit = role.prepare(source, spec)
        fill = role.fill
    xmin, xmax = spec.xmin.to(unit).value, spec.xmax.to(unit).value
    result = model_on_grid(m_x, m_xmin, m_xmax, values, spec.x.value, xmin, xmax, fill)
    n, n_cov = len(result.values), int(np.count_nonzero(result.covered))
    if n_cov == 0:
        raise RecipeFailed(
            f"The {role.label} of {source.where} ({_range(np.nanmin(m_xmin), np.nanmax(m_xmax), unit)}) "
            f"does not overlap the spectrum ({_range(np.nanmin(xmin), np.nanmax(xmax), unit)}).")
    if role.check is not None and not normalized:
        role.check(result.values, result.covered, source, spec)
    grid_values = result.values
    if normalized:
        assert norm is not None
        grid_values = grid_values * norm
        values_unit = spec.y.unit
    elif role.beyond_is_norm:
        beyond = np.full(n, np.nan) if norm is None else norm
        grid_values = np.where(result.covered, grid_values, beyond)
    how = ("copied (the model is on the grid of the spectrum)" if result.method == 'copy'
           else "averaged over the pixels of the spectrum")
    if role.beyond_is_norm:
        fill_text = "the normalization (cont × telluric)" if norm is not None else "NaN"
    else:
        fill_text = "NaN" if np.isnan(role.fill) else f"{role.fill:g}"
    beyond_text = "" if n_cov == n else f"; beyond, {fill_text}"
    scaled = (" It was normalized: multiplied by the normalization (cont × telluric)."
              if normalized else "")
    replaced = " It replaces the previous one." if role.column in spec.t.colnames else ""
    logger.info(f"The {role.label} of {source.where} (column {source.column}) was {how} "
                f"into the column {role.column}: it covers {n_cov} of {n} pixels"
                f"{beyond_text}.{scaled}{replaced}")
    new_spec = spec.update_column(role.column, grid_values, unit=values_unit,
                                  description=f"{role.label} from {source.where}")
    if role.renormalizes_model and renorm_model:
        renormalized = _renormalized_model(spec, new_spec)
        if renormalized is not None:
            new_spec = renormalized
    if role.after is not None:
        message = role.after(new_spec)
        if message:
            logger.warning(message)
    return new_spec
