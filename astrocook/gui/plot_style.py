"""Style of the plots: the strong emission lines, the colours of the
plotted elements, and the helpers (colour cycle, x to redshift, min-max
decimation of long spectra).

Moved out of ``pyside_plot.py`` in phase 5 (split of the GUI). The plots
stay white also in dark mode (Guido's decision): their colours do not
follow the theme of the interface, except the selection, shared with it
(``theme.SELECTION``).
"""
import logging
import re

import matplotlib
import numpy as np

from .theme import SELECTION

#: The error of the flux (band, line or bars), in every plot
ERROR_COLOR = '#aaaaaa'

#: The ticks of the HI components on the main plot: a softer red than pure
#: red, drawn translucent
HI_TICK_COLOR = '#d62728'


# Keys must match exactly what is in xem_d (core/atomic_data.py)
STRONG_EMISSION_LINES = {
    'Ly_a', 'Ly_b', 'NV_1238', 'SiIV_1393', 'CIV_1548', 
    'CIII_1908', 'MgII_2796', 'H_b', 'OIII_5008', 'H_a'
}

#: The names of the strong emission lines on the plot: ions in spectroscopic
#: notation, the forbidden and semi-forbidden lines with their brackets
EMISSION_LABELS = {
    'Ly_a': 'Ly\u03b1', 'Ly_b': 'Ly\u03b2', 'NV_1238': 'N V', 'SiIV_1393': 'Si IV',
    'CIV_1548': 'C IV', 'CIII_1908': 'C III]', 'MgII_2796': 'Mg II',
    'H_b': 'H\u03b2', 'OIII_5008': '[O III]', 'H_a': 'H\u03b1',
}

_GREEK = {'a': '\u03b1', 'b': '\u03b2', 'g': '\u03b3', 'd': '\u03b4', 'e': '\u03b5'}

#: The ions of the other lines: element, ionization (Roman numerals),
#: excited levels (*), rest wavelength, and a letter for close lines
_ION_KEY = re.compile(r'^([A-Z][a-z]?)([IVX]+)(\**)_(\d+)([a-z]?)$')


def emission_label(key: str) -> str:
    """The name of a line of ``xem_d`` as written on the plot (e.g.
    'CIV_1548' -> 'C IV', 'Ly_b' -> 'Ly\u03b2', 'CII_1334' -> 'C II 1334');
    the key itself if it follows no known pattern."""
    if key in EMISSION_LABELS:
        return EMISSION_LABELS[key]
    match = _ION_KEY.match(key)
    if match:
        element, ion, stars, wave, letter = match.groups()
        return f"{element} {ion}{stars} {wave}{letter}"
    series, _, member = key.partition('_')
    if series in ('Ly', 'H', 'D') and member:
        name = 'Ly' if series in ('Ly', 'D') else 'H'
        prefix = 'D ' if series == 'D' else ''
        if member == 'lim':
            return f"{prefix}{name} limit"
        return f"{prefix}{name}{_GREEK.get(member, member)}"
    return key


def flux_unit_label(unit) -> str:
    """The unit of the flux as written on the axis: the usual order for
    f_lambda and f_nu (erg s\u207b\u00b9 cm\u207b\u00b2 \u00c5\u207b\u00b9), Unicode for the
    others. Astropy writes erg \u00c5\u207b\u00b9 s\u207b\u00b9 cm\u207b\u00b2, and in mathtext the \u00c5
    was drawn out of place (visual review of phase 5-bis)."""
    import astropy.units as au
    per_s_cm2 = '\u207b\u00b9 cm\u207b\u00b2'
    scale = getattr(unit, 'scale', 1.0)
    base = unit / scale if scale != 1.0 else unit
    prefix = f"{au.Unit(scale).to_string('unicode')} " if scale != 1.0 else ""
    for per, name in ((au.Angstrom, '\u00c5'), (au.nm, 'nm'), (au.Hz, 'Hz')):
        if base == au.erg / (au.s * au.cm**2 * per):
            return f"{prefix}erg s{per_s_cm2} {name}\u207b\u00b9"
    return unit.to_string('unicode')


# --- Centralized Plot Style Configuration ---
# The span of a region selected with the Select tool (F73): the colour of
# the selection elsewhere (the highlighted group of the System Inspector)
SELECTION_STYLE = {'color': SELECTION, 'alpha': 0.25, 'linewidth': 0}

PLOT_STYLE = {
    'flux': {
        'color_idx': 0,       # Index in color cycle (tab20)
        'lw': 0.5,
        'alpha': 1.0,
        'label': 'Flux'
    },
    'model': {
        'color_idx': 1,       # Index 2 is usually Orange/Red in tab20
        'ls': '-',
        'lw': 0.8,
        'alpha': 1.0,
        'label': 'Model'
    },
    'error': {
        'color': ERROR_COLOR,   # Grey shading
        'lw': 0.3,
        'alpha': 1.0
    },
    'continuum': {
        'color': 'black',
        'ls': '--',
        'lw': 0.8
    },
    'cursor': {
        'color_idx': 3,
        'ls': '--',
        'lw': 1.0,
        'alpha': 0.8
    },
    'zero_line': {
        'color': 'gray',
        'ls': '-',
        'lw': 0.5,
        'alpha': 1.0
    },
    'center_line': {
        'color': 'gray',
        'ls': '-',
        'lw': 0.5,
        'alpha': 1.0
    }
}

def get_style_color(key: str, cycle: list):
    """Helper to get color from style dict or cycle."""
    if 'color' in PLOT_STYLE[key]:
        return PLOT_STYLE[key]['color']
    else:
        idx = PLOT_STYLE[key].get('color_idx', 0)
        return cycle[idx % len(cycle)]

# The x_plot it receives is *always* in the data's native units (nm)
def z_convert_inverse(x_plot_nm, xem_nm): # <<< Removed zem_spec, x_unit
    """Converts plot X coordinate (in data units, nm) back to redshift."""
    try:
        z = (x_plot_nm / xem_nm) - 1.0 # Simple calculation
        return z
    except Exception as e:
        logging.debug(f"Inverse z conversion failed: {e}")
        return None

# --- Helper to get colors ---
def get_color_cycle(n=10, cmap=None, fallback_cmap='tab20'):
    # ... (function definition as before) ...
    if cmap is not None:
        try:
            cmap_obj = matplotlib.colormaps.get_cmap(cmap)
            return [cmap_obj(i / n) for i in range(n)]
        except Exception as e:
            logging.warning(f"Could not get specified colormap '{cmap}': {e}. Using fallback.")
    try:
        colors = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']

        return [colors[i % len(colors)] for i in range(n)]
    except Exception as e:
        logging.warning(f"Could not get color cycle: {e}. Using fallback.")
        cmap = matplotlib.colormaps.get_cmap(fallback_cmap)
        return [cmap(i / n) for i in range(n)]

# ----------------------------

def decimate_y_min_max(y: np.ndarray, factor: int) -> np.ndarray:
    """
    Performs a min-max decimation on a Y-axis array.
    Returns an array of [y_min1, y_max1, y_min2, y_max2, ...]
    """
    if y is None:
        return None
    try:
        # 1. Trim array to be divisible by factor
        size = (len(y) // factor) * factor
        if size == 0:
            return y # Not enough data to decimate, return original
            
        y_trim = y[:size]
        
        # 2. Reshape
        y_reshaped = y_trim.reshape(-1, factor)
        
        # 3. Get min/max for y
        y_mins = y_reshaped.min(axis=1)
        y_maxs = y_reshaped.max(axis=1)
        
        # 4. Interleave
        return np.stack((y_mins, y_maxs), axis=1).ravel()
    except Exception as e:
        logging.warning(f"y-axis decimation failed: {e}")
        return y # Return original on failure
        
def decimate_x_for_min_max(x: np.ndarray, factor: int) -> np.ndarray:
    """
    Decimates an X-axis to match a min-max-decimated Y-axis.
    Returns an array of [x_start1, x_end1, x_start2, x_end2, ...]
    """
    try:
        # 1. Trim array to be divisible by factor
        size = (len(x) // factor) * factor
        if size == 0:
            return x # Not enough data to decimate, return original

        x_trim = x[:size]
        
        # 2. Reshape
        x_reshaped = x_trim.reshape(-1, factor)
        
        # 3. Get start/end for x
        x_starts = x_reshaped[:, 0]  # Get the start x of each bin
        x_ends = x_reshaped[:, -1] # Get the end x of each bin
        
        # 4. Interleave
        return np.stack((x_starts, x_ends), axis=1).ravel()
    except Exception as e:
        logging.warning(f"x-axis decimation failed: {e}")
        return x # Return original on failure
