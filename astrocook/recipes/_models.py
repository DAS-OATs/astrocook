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
(and conversion) of its values. T2 has the role ``telluric``; ``cont``
(T3) and ``model`` (T4) are entries of the table, with no change of the
interface of the recipe.
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


#: How far a transmission may go beyond 0-1 without a warning (T2: a
#: Molecfit model of the request reaches 1 plus a little, from the numerics
#: of the model, and the warning said "from 0.0237 to 1")
TRANSMISSION_TOLERANCE = 1e-3


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


#: The roles of the models, by the value of the parameter ``role``
MODEL_ROLES = {
    'telluric': ModelRole(column=TELLURIC_COLUMN, names=TELLURIC_NAMES,
                          fill=1.0, label="telluric model",
                          prepare=_prepare_telluric, after=_after_telluric),
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


def import_model(spec: 'Spectrum', source: ModelSource, role_name: str) -> 'Spectrum':
    """The spectrum with the model in the column of its role."""
    role = MODEL_ROLES.get(role_name)
    if role is None:
        raise InvalidParameter(f"role: {role_name!r} is not one of "
                               f"{', '.join(repr(r) for r in MODEL_ROLES)}")
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
    values, values_unit = role.prepare(source, spec)
    xmin, xmax = spec.xmin.to(unit).value, spec.xmax.to(unit).value
    result = model_on_grid(m_x, m_xmin, m_xmax, values, spec.x.value, xmin, xmax,
                           role.fill)
    n, n_cov = len(result.values), int(np.count_nonzero(result.covered))
    if n_cov == 0:
        raise RecipeFailed(
            f"The {role.label} of {source.where} ({_range(np.nanmin(m_xmin), np.nanmax(m_xmax), unit)}) "
            f"does not overlap the spectrum ({_range(np.nanmin(xmin), np.nanmax(xmax), unit)}).")
    how = ("copied (the model is on the grid of the spectrum)" if result.method == 'copy'
           else "averaged over the pixels of the spectrum")
    beyond = "" if n_cov == n else f"; beyond, {role.fill:g}"
    replaced = " It replaces the previous one." if role.column in spec.t.colnames else ""
    logger.info(f"The {role.label} of {source.where} (column {source.column}) was {how} "
                f"into the column {role.column}: it covers {n_cov} of {n} pixels{beyond}."
                f"{replaced}")
    new_spec = spec.update_column(role.column, result.values, unit=values_unit,
                                  description=f"{role.label} from {source.where}")
    if role.after is not None:
        message = role.after(new_spec)
        if message:
            logger.warning(message)
    return new_spec
