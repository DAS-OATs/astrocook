"""The units of the spectra that a recipe combines (phase 6-quater, PR 9).

The rule agreed with Guido (F132-F134, F138): spectra are combined, or one
is scaled to another, in the units of a reference.

- An axis is converted between wavelengths (e.g. Å to nm); an axis in
  velocity against one in wavelength is refused (``x_convert`` first).
- A flux is converted when its unit is of the same kind as that of the
  reference (e.g. per nm to per Å). A flux in arbitrary units (``arb_flux``)
  against a physical one is combined only when it is scaled to it, and then
  takes the unit of the reference: the scaling calibrates it. A normalized
  flux (no unit) is never combined with one that is not.
- Nothing is done when the units are the same, so those results do not
  move.
"""
from __future__ import annotations

import logging
from typing import List

import astropy.units as au

from astrocook.core.units import has_no_unit, is_arbitrary_flux
from astrocook.recipes.errors import RecipeFailed

logger = logging.getLogger(__name__)


def unit_words(unit) -> str:
    """A unit of the flux as the messages say it."""
    if is_arbitrary_flux(unit):
        return "arbitrary units"
    if has_no_unit(unit):
        return "no unit (normalized)"
    return str(unit)


def check_axes(x_unit, ref_unit, name: str, ref_name: str) -> None:
    """RecipeFailed if the axis of ``name`` cannot be compared with that of
    ``ref_name``: the same unit, or two wavelengths."""
    x_u, x_ref = au.Unit(x_unit), au.Unit(ref_unit)
    if x_u != x_ref and not (x_u.is_equivalent(au.nm) and x_ref.is_equivalent(au.nm)):
        raise RecipeFailed(f"The axis of {name} is in {x_u} and that of {ref_name} in "
                           f"{x_ref}: convert them to wavelengths with x_convert.")


def flux_relation(y_unit, ref_unit, name: str, ref_name: str, scaled: bool) -> str:
    """How the flux of ``name`` joins that of ``ref_name``: ``'same'`` (the
    same unit), ``'convert'`` (a unit of the same kind), ``'declare'`` (an
    arbitrary flux and a physical one, ``scaled``: the scaled flux takes the
    unit of the reference); RecipeFailed otherwise."""
    if au.Unit(y_unit) == au.Unit(ref_unit):
        return 'same'
    if au.Unit(y_unit).is_equivalent(ref_unit):
        return 'convert'
    if scaled and (is_arbitrary_flux(y_unit) or is_arbitrary_flux(ref_unit)) \
            and not (has_no_unit(y_unit) or has_no_unit(ref_unit)):
        logger.info(f"The flux of {name}, in {unit_words(y_unit)}, takes the unit of "
                    f"{ref_name} ({unit_words(ref_unit)}) once scaled to it.")
        return 'declare'
    how = ("" if scaled else "; scaled to it (Equalize, Equalize and stitch arms, or "
           "Co-add with equalize_order), an arbitrary flux takes the unit of the other")
    raise RecipeFailed(f"The flux of {name} is in {unit_words(y_unit)} and that of {ref_name} "
                       f"in {unit_words(ref_unit)}: they cannot be combined{how}. "
                       "Set_properties declares the unit of a flux (yunit).")


def in_units_of(spec, ref, name: str, ref_name: str, scaled: bool):
    """``spec`` with the axis and the flux in the units of ``ref`` (see the
    rule of the module)."""
    check_axes(spec.x.unit, ref.x.unit, name, ref_name)
    if au.Unit(spec.x.unit) != au.Unit(ref.x.unit):
        spec = spec.x_convert(spec._data.z_em or 0.0, au.Unit(ref.x.unit))
    how = flux_relation(spec.y.unit, ref.y.unit, name, ref_name, scaled)
    if how == 'convert':
        return spec.y_convert(ref.y.unit)
    if how == 'declare':
        return spec.with_flux_unit(ref.y.unit)
    return spec


def in_units_of_first(specs: list, names: List[str], scaled: bool) -> list:
    """The spectra, the first as it is, the others in its units."""
    return [specs[0]] + [in_units_of(spec, specs[0], name, names[0], scaled)
                         for spec, name in zip(specs[1:], names[1:], strict=True)]
