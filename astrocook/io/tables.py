"""One reader of the tables of the files, by the names of their columns
(phase 7-ter, PR T2).

:func:`read_table` reads a file as an astropy ``Table`` whose columns have
the names of the file and, where the file says it, a unit:

- ECSV (the CSV of ``export_ascii``), with its units and metadata;
- a text table: columns separated by spaces, tabs or commas; the names
  from a header line (``x xmin xmax y dy``), also commented
  (``# x xmin xmax trans_molecfit dy``, the last comment line before the
  data with as many names as columns), or ``col1``, ``col2``, ... without
  one;
- FITS: a binary table (the files of Astrocook, the Molecfit outputs with
  ``mtrans``); a table of one row whose cells are arrays is flattened (the
  ESO Phase 3 format); the extension holding one of the columns asked for
  is chosen, and if it has no axis, the axis of another extension of the
  same length (the ``TELLURIC_DATA`` of the X-Shooter products).

The columns are found by name with :func:`find_column` (exact name first,
then ignoring the case). The axis (:func:`axis_of`) is in the unit of the
file, or without one in nm for a text file (as the ASCII loaders) and by
the rule of the FITS loaders for a FITS file (Å above 3000, nm below).

Used by ``import_model`` and by the loader of the ASCII spectra with named
columns (FT2), which share the names below.
"""
from __future__ import annotations

import logging
import re
from typing import Iterable, Optional, Sequence

import astropy.units as au
import numpy as np
from astropy.io import ascii, fits
from astropy.table import Column, Table

from astrocook.core.model_grid import auto_edges
from astrocook.core.units import known_unit

logger = logging.getLogger(__name__)

#: The names of the axis, of its edges, of the flux and of its error, in
#: order of preference: those of Astrocook (data, GUI in ASCII) first
AXIS_NAMES = ('x', 'lam', 'lambda', 'wave', 'wavelength')
XMIN_NAMES = ('xmin', 'lammin', 'wave_min', 'wmin')
XMAX_NAMES = ('xmax', 'lammax', 'wave_max', 'wmax')
FLUX_NAMES = ('y', 'F', 'flux')
ERROR_NAMES = ('dy', 'dF', 'err', 'error', 'sigma')

#: Where :func:`read_table` read the table (``table.meta['format']``)
FORMAT_KEY = 'format'

_NUMBER = re.compile(r'^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$|^[+-]?(nan|inf)$',
                     re.IGNORECASE)


def find_column(table: Table, names: Iterable[str],
                exclude: Sequence[str] = ()) -> Optional[str]:
    """The name in ``table`` of the first of ``names`` it has: the exact
    name first, then ignoring the case; None if it has none. Columns in
    ``exclude`` are not taken."""
    names = list(names)
    cols = [c for c in table.colnames if c not in exclude]
    for name in names:
        if name in cols:
            return name
    lower = {c.lower(): c for c in reversed(cols)}
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    return None


# --- Text files -----------------------------------------------------------------

def _is_number(token: str) -> bool:
    return bool(_NUMBER.match(token.strip()))


def _split(line: str, delimiter: Optional[str]) -> list[str]:
    if delimiter is None:
        return line.split()
    return [t.strip() for t in line.split(delimiter)]


def _delimiter(line: str) -> Optional[str]:
    """The separator of a line of data: comma, tab, or spaces (None)."""
    if ',' in line:
        return ','
    if '\t' in line.strip():
        return '\t'
    return None


def _text_layout(path: str):
    """The layout of a text table: the separator, the names of the columns
    (None without a header) and the number of lines before the data."""
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    comments: list[str] = []
    for i, raw in enumerate(lines):
        line = raw.strip()
        if not line:
            continue
        if line.startswith('#'):
            comments.append(line.lstrip('#').strip())
            continue
        delimiter = _delimiter(line)
        tokens = _split(line, delimiter)
        if not all(_is_number(t) for t in tokens):
            # The header, as a line of names
            return delimiter, tokens, i + 1
        names = None
        for comment in reversed(comments):
            parts = _split(comment, delimiter)
            if delimiter is not None and len(parts) != len(tokens):
                parts = comment.split()
            if len(parts) == len(tokens) and not any(_is_number(p) for p in parts):
                names = parts
                break
        return delimiter, names, i
    raise ValueError(f"{path} has no data.")


def text_column_names(path: str) -> Optional[list[str]]:
    """The names of the columns of a text table (None without a header),
    for the detection of the format; ValueError if it is not one."""
    _, names, _ = _text_layout(path)
    return names


def is_text_table(path: str) -> bool:
    """Whether the file is a text table: its first line of data has at
    least two numbers."""
    try:
        with open(path, 'r', encoding='utf-8', errors='strict') as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith('#'):
                    continue
                tokens = _split(line, _delimiter(line))
                if all(_is_number(t) for t in tokens):
                    return len(tokens) >= 2
                # A header: the next line decides
                continue
    except (OSError, UnicodeDecodeError):
        return False
    return False


def _is_ecsv(path: str) -> bool:
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        return f.readline().startswith('# %ECSV')


def _read_text(path: str) -> Table:
    if _is_ecsv(path):
        table = Table.read(path, format='ascii.ecsv')
        table.meta[FORMAT_KEY] = 'ecsv'
        return table
    delimiter, names, data_start = _text_layout(path)
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        lines = [ln for ln in f.readlines()[data_start:]
                 if ln.strip() and not ln.lstrip().startswith('#')]
    kwargs = {'format': 'no_header', 'guess': False, 'fast_reader': False}
    if delimiter is not None:
        kwargs['delimiter'] = delimiter
    table = ascii.read(lines, **kwargs)
    if names is not None:
        if len(names) != len(table.colnames):
            raise ValueError(f"{path}: {len(names)} names in the header and "
                             f"{len(table.colnames)} columns.")
        table = Table([table[c] for c in table.colnames], names=names)
    table.meta[FORMAT_KEY] = 'ascii'
    return table


# --- FITS files -----------------------------------------------------------------

def _fits_unit(text):
    if text is None or not str(text).strip():
        return None
    try:
        return known_unit(au.Unit(str(text).strip(), format='fits',
                                  parse_strict='silent'))
    except ValueError:
        return None


def _hdu_table(hdu) -> Table:
    """The columns of a binary table, flattened if it is one row of arrays."""
    data = hdu.data
    table = Table()
    if data is None:
        return table
    one_row = len(data) == 1 and any(np.ndim(data[c.name][0]) >= 1 for c in hdu.columns)
    for col in hdu.columns:
        values = np.asarray(data[col.name])
        if one_row:
            values = values.reshape(-1)
        elif values.ndim > 1:
            continue  # a column of arrays per pixel: not a column of a model
        table[col.name] = Column(values, unit=_fits_unit(col.unit))
    return table


def _read_fits(path: str, wanted: Sequence[str]) -> Table:
    with fits.open(path) as hdul:
        tables = [(hdu.name, _hdu_table(hdu)) for hdu in hdul
                  if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU))]
        header = dict(hdul[0].header)
    tables = [(n, t) for n, t in tables if len(t.colnames)]
    if not tables:
        raise ValueError(f"{path} has no table.")
    chosen = next(((n, t) for n, t in tables if find_column(t, wanted)), None) \
        or next(((n, t) for n, t in tables if find_column(t, AXIS_NAMES)), None) \
        or tables[0]
    name, table = chosen
    if find_column(table, AXIS_NAMES) is None:
        # The axis of another extension of the same length
        donor = next((t for n, t in tables if n != name and len(t) == len(table)
                      and find_column(t, AXIS_NAMES)), None)
        if donor is not None:
            for names in (AXIS_NAMES, XMIN_NAMES, XMAX_NAMES):
                col = find_column(donor, names)
                if col is not None and col not in table.colnames:
                    table[col] = donor[col]
    table.meta[FORMAT_KEY] = 'fits'
    table.meta['EXTNAME'] = name
    for key in ('ORIGIN', 'INSTRUME', 'OBJECT'):
        if key in header:
            table.meta[key] = header[key]
    return table


def _is_fits(path: str) -> bool:
    with open(path, 'rb') as f:
        return f.read(6) == b'SIMPLE'


def read_table(path: str, wanted: Sequence[str] = ()) -> Table:
    """The table of a file (module docstring); ``wanted`` are names of
    columns that choose the extension of a FITS file. ValueError if the
    file is not a table."""
    if _is_fits(path):
        return _read_fits(path, wanted)
    return _read_text(path)


def axis_of(table: Table):
    """The axis of a table and the edges of its pixels, as quantities:
    ``x``, ``xmin``, ``xmax`` (the edges halfway between the centres if the
    table has none); the first column of a text table without a header;
    ValueError without an axis. The unit is that of the
    file; without one, nm for a text file and the rule of the FITS loaders
    for a FITS file (Å above 3000, nm below)."""
    name = find_column(table, AXIS_NAMES)
    if name is None and table.colnames and table.colnames[0] == 'col1':
        name = 'col1'  # a text table without a header: the first column
    if name is None:
        raise ValueError(f"No column of the axis ({', '.join(AXIS_NAMES)}) among "
                         f"{', '.join(table.colnames)}.")
    x = np.asarray(table[name], dtype=float)
    unit = known_unit(table[name].unit)
    if unit is None or unit == au.dimensionless_unscaled:
        if table.meta.get(FORMAT_KEY) == 'fits' and np.nanmedian(x) > 3000:
            unit = au.Angstrom
        else:
            unit = au.nm
    lo_name = find_column(table, XMIN_NAMES, exclude=(name,))
    hi_name = find_column(table, XMAX_NAMES, exclude=(name,))
    if lo_name is not None and hi_name is not None:
        def edge(col):
            col_unit = known_unit(table[col].unit)
            q = np.asarray(table[col], dtype=float) * (col_unit if col_unit not in
                                                        (None, au.dimensionless_unscaled)
                                                        else unit)
            return q.to(unit, equivalencies=au.spectral())
        xmin, xmax = edge(lo_name), edge(hi_name)
    else:
        lo, hi = auto_edges(x)
        xmin, xmax = lo * unit, hi * unit
    return x * unit, xmin, xmax
