"""The telluric transmission of a spectrum (phase 7-ter, PR T1).

The transmission is the auxiliary column ``telluric`` (dimensionless, 1
where there is no telluric absorption): the name of
``design/identification.md`` (section 7) and of phase 7, 7b. The core reads
it in one place, :func:`telluric_transmission`: the normalization
(``Spectrum.norm``) is ``cont × telluric``, unless the flux is already
corrected for the tellurics, which the spectrum says with the flag
``TELLURIC_CORRECTED`` of its metadata (then the normalization is ``cont``;
the column stays, for the identification).

:func:`migrate_telluric` is applied when a spectrum is built (as the
resolution, ``core/resolution.py``), so that every source (archives, files,
scripts) gives the same column:

- ``telluric_model``, the name read by the core before 2.8.0, becomes
  ``telluric``;
- a boolean ``telluric`` is the mask of V1's ``telluric_mask`` (True in the
  bands), not a transmission: it becomes ``telluric_mask``;
- a ``telluric`` already there when ``telluric_model`` is renamed (written
  by the ASCII telluric loader and read by nothing before 2.8.0) is kept as
  ``telluric_unused``: the normalization stays the one of
  ``telluric_model``.
"""
from __future__ import annotations

import dataclasses
import logging
from typing import Mapping, Optional

import numpy as np

logger = logging.getLogger(__name__)

#: The column of the transmission
TELLURIC_COLUMN = 'telluric'
#: Its name before 2.8.0, renamed when a spectrum is built
LEGACY_TELLURIC_COLUMN = 'telluric_model'
#: The V1 mask of the telluric bands (a boolean ``telluric`` column of V1)
V1_MASK_COLUMN = 'telluric_mask'
#: A ``telluric`` column met with ``telluric_model``, kept under this name
UNUSED_TELLURIC_COLUMN = 'telluric_unused'
#: The flag of the metadata: the flux is already corrected for the tellurics
TELLURIC_CORRECTED = 'TELLURIC_CORRECTED'
#: The names of a column of a telluric transmission in the files, in order
#: of preference (matched ignoring the case): Astrocook, before 2.8.0, the
#: Molecfit outputs (``trans_molecfit`` of the ASCII files of the request of
#: 7-ter, ``mtrans`` of its FITS tables), PypeIt (``TELLURIC``); T2
TELLURIC_NAMES = (TELLURIC_COLUMN, LEGACY_TELLURIC_COLUMN, 'trans_molecfit', 'mtrans')

_TRUE = {'t', 'true', 'yes', 'y', '1'}


def is_telluric_corrected(meta: Mapping) -> bool:
    """Whether the metadata say that the flux is already corrected for the
    tellurics (the flag ``TELLURIC_CORRECTED``: a bool, or ``T``/``F`` and
    the like as a FITS header or a user may write it); False without it."""
    value = meta.get(TELLURIC_CORRECTED)
    if value is None:
        # A HIERARCH card may be read back in another case
        value = next((v for k, v in meta.items()
                      if str(k).upper() == TELLURIC_CORRECTED), None)
    if value is None:
        return False
    if isinstance(value, str):
        return value.strip().lower() in _TRUE
    return bool(value)


def telluric_transmission(data) -> Optional[np.ndarray]:
    """The transmission that enters the normalization of a spectrum (its
    :class:`~astrocook.core.structures.SpectrumData`): the values of the
    ``telluric`` column, or None without it or when the flux is already
    corrected (``TELLURIC_CORRECTED``)."""
    col = data.aux_cols.get(TELLURIC_COLUMN)
    if col is None or is_telluric_corrected(data.meta):
        return None
    return col.values


def _free_name(base: str, names) -> str:
    name, i = base, 2
    while name in names:
        name, i = f"{base}_{i}", i + 1
    return name


def migrate_telluric(data):
    """``data`` with the telluric columns under their names (module
    docstring); the same object when there is nothing to rename. The order
    of the columns is kept."""
    cols = data.aux_cols
    renames: dict[str, str] = {}
    current = cols.get(TELLURIC_COLUMN)
    if current is not None and np.asarray(current.values).dtype == bool:
        renames[TELLURIC_COLUMN] = _free_name(V1_MASK_COLUMN, cols)
        logger.info(f"The boolean column '{TELLURIC_COLUMN}' (the mask of the "
                    f"telluric bands of V1) is renamed '{renames[TELLURIC_COLUMN]}': "
                    "it is not a transmission.")
    if LEGACY_TELLURIC_COLUMN in cols:
        if current is not None and TELLURIC_COLUMN not in renames:
            renames[TELLURIC_COLUMN] = _free_name(UNUSED_TELLURIC_COLUMN, cols)
            logger.warning(f"The column '{TELLURIC_COLUMN}' is renamed "
                           f"'{renames[TELLURIC_COLUMN]}': the transmission is "
                           f"'{LEGACY_TELLURIC_COLUMN}', which becomes '{TELLURIC_COLUMN}'.")
        renames[LEGACY_TELLURIC_COLUMN] = TELLURIC_COLUMN
    if not renames:
        return data
    new_cols = {renames.get(name, name): col for name, col in cols.items()}
    return dataclasses.replace(data, aux_cols=new_cols)
