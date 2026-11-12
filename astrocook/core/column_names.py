"""The names of the columns of a spectrum (phase 6-quater, PR 8, F120).

One set of names for the user and the agents: those of the GUI, ``λ``,
``λmin``, ``λmax``, ``F``, ``dF`` for the columns that every spectrum has,
and the other columns by their own name (``cont``, ``model``, ``resol``,
...). The data keep their own names, ``x``, ``xmin``, ``xmax``, ``y``,
``dy`` (the archives, the CSV files of ``export_ascii``, the attributes
``Spectrum.x``, ``.y``, ``.dy``): accepted wherever a column is named, as
silent and permanent aliases, as ``lambda``, ``lambdamin`` and
``lambdamax``, the ASCII spelling of ``λ`` (Guido, 12 November 2026). So
the logs, the saved expressions and the scripts written with either set
still work, and nothing is translated in them.

The names of the parameters that name an *axis* and not a column keep
``x`` and ``y`` (``xunit``, ``xstart``, ``dx``, ``x_min`` of the tools,
``x_convert``).
"""
import re
from typing import Dict, Iterable, List

#: The names of the columns of every spectrum: name shown -> name of the data
GUI_COLUMN_NAMES: Dict[str, str] = {
    'λ': 'x', 'λmin': 'xmin', 'λmax': 'xmax', 'F': 'y', 'dF': 'dy'}

#: The ASCII spelling of the names of the axis
ASCII_COLUMN_NAMES: Dict[str, str] = {
    'lambda': 'x', 'lambdamin': 'xmin', 'lambdamax': 'xmax'}

#: Every name accepted for a column of every spectrum, with its name in the data
_TO_INTERNAL: Dict[str, str] = {**GUI_COLUMN_NAMES, **ASCII_COLUMN_NAMES}

_TO_DISPLAY: Dict[str, str] = {v: k for k, v in GUI_COLUMN_NAMES.items()}

#: The columns of every spectrum, by the names of the data
CORE_COLUMNS = ('x', 'xmin', 'xmax', 'y', 'dy')

#: The columns of the axis, by the names of the data
AXIS_COLUMNS = ('x', 'xmin', 'xmax')

# The longest first (λmin before λ), whole words only (λmin is not λ + min,
# F_corr and s1_F are not F)
_PATTERN = re.compile(
    r'\b(' + '|'.join(re.escape(n) for n in sorted(_TO_INTERNAL, key=len, reverse=True))
    + r')\b')


def internal_column(name: str) -> str:
    """The name in the data of the column called ``name`` (any spelling;
    a name of another column as it is)."""
    return _TO_INTERNAL.get(name, name)


def display_column(name: str) -> str:
    """The name shown of the column called ``name`` in the data."""
    return _TO_DISPLAY.get(name, name)


def display_columns(names: Iterable[str]) -> List[str]:
    """:func:`display_column` of each name."""
    return [display_column(n) for n in names]


def internal_expression(text: str) -> str:
    """``text`` (an expression, or a comma-separated list of names) with the
    names of the columns of every spectrum written as in the data."""
    return _PATTERN.sub(lambda m: _TO_INTERNAL[m.group(1)], text)


def is_core_column(name: str) -> bool:
    """Whether ``name`` (any spelling) is a column of every spectrum."""
    return internal_column(name) in CORE_COLUMNS


def is_axis_column(name: str) -> bool:
    """Whether ``name`` (any spelling) is a column of the axis (λ, λmin, λmax)."""
    return internal_column(name) in AXIS_COLUMNS
