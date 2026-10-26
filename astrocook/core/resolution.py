"""The spectral resolution: one source, the ``resol`` column (#526).

The ``resol`` column of a spectrum holds the resolution of each pixel. Its
unit says what it is: dimensionless for the resolving power R, ``km/s``
for the FWHM in velocity. A component of a system list can override it
with its own value (always R).

Before #526 the resolution could also live in the ``SpectrumDataV2.resol``
attribute and in the meta keys ``resol`` (R), ``RESOL`` (R, V1 archives)
and ``resol_fwhm`` (km/s), and its unit was guessed from its value (below
500: FWHM in km/s). These legacy sources are migrated into the column when
a spectrum is built (:func:`migrate_resolution`), and no longer read
anywhere else; the guess of the unit survives only there, for data saved
before #526. Since #526 (PR 4) nothing writes them: the attribute is
cleared by the migration, the meta keys are left as they came (they may be
keywords of the original file).
"""
import dataclasses
import logging
from typing import Optional, Tuple

import astropy.units as au
import numpy as np

C_KMS = 299792.458

#: Legacy rule: a resolution below this value was taken as a FWHM in km/s,
#: above as R. Used only to migrate data saved before #526.
LEGACY_KMS_THRESHOLD = 500.0

#: A component resolution within this fraction of the resolution of the
#: spectrum at its position was stamped by the fitter before #526 (it wrote
#: the resolution it used into every fitted component), not set by the
#: user: it is cleared at load, so that the component follows the column.
STAMPED_RTOL = 0.01


def column_to_resolving_power(values: np.ndarray, unit) -> np.ndarray:
    """Resolving power R from the values of a ``resol`` column.

    Parameters
    ----------
    values : np.ndarray
        The values of the column.
    unit : astropy.units.Unit
        Their unit: dimensionless (R) or a velocity (FWHM).

    Returns
    -------
    np.ndarray
        R for each pixel, NaN where the value is missing or not positive.
    """
    values = np.asarray(values, dtype=float)
    unit = au.Unit(unit) if unit is not None else au.dimensionless_unscaled
    with np.errstate(invalid='ignore', divide='ignore'):
        if unit.physical_type == 'speed':
            fwhm_kms = (values * unit).to_value(au.km / au.s)
            r = C_KMS / fwhm_kms
        else:
            r = values * unit.to(au.dimensionless_unscaled)
    return np.where(np.isfinite(r) & (r > 0), r, np.nan)


def _legacy_to_r(value: float) -> float:
    """A legacy scalar resolution, R or FWHM in km/s by its value."""
    return C_KMS / value if value < LEGACY_KMS_THRESHOLD else value


def legacy_resolution(data) -> Tuple[Optional[float], str]:
    """Resolving power from the legacy sources of a spectrum, in the order
    the fitter read them before #526: attribute, meta ``resol``, meta
    ``RESOL``, meta ``resol_fwhm``. Returns ``(None, '')`` if none is set."""
    def positive(value):
        try:
            value = float(value)
        except (TypeError, ValueError):
            return None
        return value if np.isfinite(value) and value > 0 else None

    value = positive(getattr(data, 'resol', None))
    if value is not None:
        return _legacy_to_r(value), "attribute"
    meta = data.meta or {}
    for key in ('resol', 'RESOL'):
        value = positive(meta.get(key))
        if value is not None:
            return _legacy_to_r(value), f"meta '{key}'"
    value = positive(meta.get('resol_fwhm'))
    if value is not None:
        return C_KMS / value, "meta 'resol_fwhm'"
    return None, ''


def migrate_resolution(data):
    """Move the resolution of a ``SpectrumDataV2`` into its ``resol`` column.

    A column with at least one valid value is kept (a dimensionless column
    with median below ``LEGACY_KMS_THRESHOLD`` was a FWHM in km/s before
    #526, and gets that unit). Otherwise the column is created, flat, from
    the legacy sources (:func:`legacy_resolution`), if any. The attribute
    is cleared either way, so that it cannot come back later (e.g. after
    the column is deleted).

    Returns
    -------
    SpectrumDataV2
        The same object if nothing changes, a new one otherwise.
    """
    from astrocook.core.structures import DataColumnV2

    if getattr(data, 'resol', 0.0):
        r_legacy, source = legacy_resolution(data)
        data = dataclasses.replace(data, resol=0.0)
        if r_legacy is not None and not _valid_column(data):
            return _with_flat_column(data, r_legacy, source)
    col = data.aux_cols.get('resol')
    if col is not None:
        values = np.asarray(col.values, dtype=float)
        valid = np.isfinite(values) & (values > 0)
        if valid.any():
            unit = au.Unit(col.unit) if col.unit is not None else au.dimensionless_unscaled
            if (unit == au.dimensionless_unscaled
                    and np.median(values[valid]) < LEGACY_KMS_THRESHOLD):
                logging.info("Resolution column read as FWHM in km/s "
                             "(dimensionless, median below %g).",
                             LEGACY_KMS_THRESHOLD)
                aux = dict(data.aux_cols)
                aux['resol'] = DataColumnV2(values, au.km / au.s,
                                            description="Resolution (FWHM km/s)")
                return dataclasses.replace(data, aux_cols=aux)
            return data

    r, source = legacy_resolution(data)
    if r is None:
        return data
    return _with_flat_column(data, r, source)


def _valid_column(data) -> bool:
    col = data.aux_cols.get('resol')
    if col is None:
        return False
    values = np.asarray(col.values, dtype=float)
    return bool((np.isfinite(values) & (values > 0)).any())


def _with_flat_column(data, r: float, source: str):
    from astrocook.core.structures import DataColumnV2
    aux = dict(data.aux_cols)
    aux['resol'] = DataColumnV2(np.full(len(data.x.values), r),
                                au.dimensionless_unscaled,
                                description=f"Resolving power (from {source})")
    return dataclasses.replace(data, aux_cols=aux)


def resolution_at(spec, x_values: np.ndarray) -> np.ndarray:
    """Resolving power of a spectrum at the given x (same unit as its x),
    interpolated; NaN where the spectrum has none."""
    r = spec.resolution()
    x = spec.x.value
    ok = np.isfinite(r)
    if not ok.any():
        return np.full(np.shape(x_values), np.nan)
    order = np.argsort(x[ok])
    return np.interp(x_values, x[ok][order], r[ok][order],
                     left=np.nan, right=np.nan)


def migrate_component_resolutions(systs, spec):
    """Component resolutions of a system list saved before #526.

    Legacy values below ``LEGACY_KMS_THRESHOLD`` were FWHM in km/s and are
    converted to R. Values within ``STAMPED_RTOL`` of the resolution of the
    spectrum at the component (its first transition) were stamped by the
    fitter and are cleared, so that the component follows the column;
    other values are kept as deliberate overrides.

    Returns
    -------
    SystemListV2
        The same object if nothing changes, a new one otherwise.
    """
    from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS

    if systs is None or spec is None or not systs.components:
        return systs
    x_unit = spec.x.unit
    changed, n_cleared, n_converted = False, 0, 0
    components = []
    for comp in systs.components:
        r = comp.resol
        if r is None or not np.isfinite(r) or r <= 0:
            components.append(comp)
            continue
        new_r = r
        if r < LEGACY_KMS_THRESHOLD:
            new_r = C_KMS / r
            n_converted += 1
        series = comp.series
        first = (STANDARD_MULTIPLETS[series][0] if series in STANDARD_MULTIPLETS
                 else series.split(',')[0].strip())
        atom = ATOM_DATA.get(first)
        if atom is not None:
            try:
                x_obs = (atom['wave'] * (1 + comp.z) * au.Angstrom).to_value(x_unit)
            except au.UnitConversionError:
                x_obs = None
            if x_obs is not None:
                r_spec = resolution_at(spec, np.array([x_obs]))[0]
                if np.isfinite(r_spec) and abs(new_r - r_spec) <= STAMPED_RTOL * r_spec:
                    new_r = None
                    n_cleared += 1
        if new_r != r:
            comp = dataclasses.replace(comp, resol=new_r)
            changed = True
        components.append(comp)
    if not changed:
        return systs
    if n_cleared or n_converted:
        logging.info("Component resolutions: %d stamped by an earlier fit cleared "
                     "(they follow the spectrum now), %d converted from km/s to R.",
                     n_cleared, n_converted)
    return systs.__class__(dataclasses.replace(systs._data, components=components))


def resolution_summary(spec) -> Tuple[float, bool]:
    """Median resolving power of a spectrum (0.0 if unknown) and whether it
    varies with wavelength, for displays and checks."""
    r = spec.resolution()
    r = r[np.isfinite(r)]
    if not r.size:
        return 0.0, False
    return float(np.median(r)), bool(not np.allclose(r, r[0], rtol=1e-6))
