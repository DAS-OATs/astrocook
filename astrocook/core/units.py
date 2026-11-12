"""The units of the flux (phase 6-quater, PR 8-bis, agreed with Guido,
12 November 2026).

A loader gives the flux a physical unit only when it reads one (``BUNIT``,
``TUNIT`` or a known descriptor of the format); otherwise the flux is in
:data:`ARBITRARY_FLUX`, the unit ``arb_flux`` of Astrocook. It is a unit of
its own, not dimensionless: a dimensionless flux is a normalized one, and an
arbitrary flux cannot be converted to a physical unit (``y_convert`` refuses
it; ``set_properties`` declares the unit, ``calibrate_from_magnitudes``
calibrates the flux).

``arb_flux`` is enabled in astropy when this module is imported (by
``core/structures.py``), so ``astropy.units.Unit('arb_flux')`` gives it. It
is not a unit of the FITS standard: the archives write it in ``TUNITn`` by
hand (astropy would drop it with a warning) and read it back with
:func:`known_unit`; other FITS readers show the string.
"""
from typing import Optional

import astropy.units as au

#: The unit of a flux whose physical unit is not known
ARBITRARY_FLUX = au.def_unit(
    'arb_flux', doc="Astrocook: a flux in arbitrary units (not calibrated)")
au.add_enabled_units([ARBITRARY_FLUX])

#: f_lambda in the unit that the calibration gives
FLAM = au.erg / (au.cm**2 * au.s * au.Angstrom)

#: The auxiliary columns in the unit of the flux, by name (F69, F130): the
#: unit cannot tell them, since a flux without a unit and a quality flag are
#: both dimensionless
FLUX_AUX_COLUMNS = ('cont', 'model', 'cont_pl')


def is_arbitrary_flux(unit) -> bool:
    """Whether ``unit`` is the arbitrary flux of Astrocook."""
    return unit is not None and unit == ARBITRARY_FLUX


def has_no_unit(unit) -> bool:
    """Whether ``unit`` is missing or dimensionless."""
    return unit is None or au.Unit(unit) == au.dimensionless_unscaled


def known_unit(unit):
    """``unit``, with an unrecognized unit (as astropy reads ``arb_flux``
    from a FITS header) replaced by the enabled unit of that name, if any."""
    if isinstance(unit, au.UnrecognizedUnit):
        try:
            return au.Unit(unit.name)
        except ValueError:
            return unit
    return unit


def parse_flux_unit(text) -> Optional[au.UnitBase]:
    """The unit of the flux written in a header (``BUNIT``, ``TUNITn``), or
    ``None`` if there is none, it cannot be read, or it is dimensionless.

    Read as a FITS unit first, then as an astropy one (``arb_flux``,
    ``erg/(s cm2 Angstrom)``)."""
    if text is None or not str(text).strip():
        return None
    text = str(text).strip()
    for fmt in ('fits', None):
        try:
            unit = au.Unit(text, format=fmt) if fmt else au.Unit(text)
        except ValueError:
            continue
        if unit == au.dimensionless_unscaled:
            return None
        return unit
    return None


def is_flux_density_per_wavelength(unit) -> bool:
    """Whether ``unit`` is a flux density per unit wavelength (f_lambda)."""
    return unit is not None and au.Unit(unit).is_equivalent(FLAM)
