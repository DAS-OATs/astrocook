"""Style of the main plot: the legend entries shown as swatches, the strong
emission lines, the colours of the plotted elements, and the helpers
(colour cycle, x to redshift, min-max decimation of long spectra).

Moved out of ``pyside_plot.py`` in phase 5 (split of the GUI).
"""
import logging
import matplotlib
import numpy as np


# Legend entries shown as swatches in the plot controls of the main window:
# while those are open they are left out of the legend, which would
# otherwise show through them
PANEL_LEGEND_LABELS = {
    'Spectrum', 'Spectrum (decimated)', '1-sigma error', 'bin size',
    'Normalization', 'Absorption model', 'HI Systems', 'Metal Systems',
}

# Keys must match exactly what is in xem_d (core/atomic_data.py)
STRONG_EMISSION_LINES = {
    'Ly_a', 'Ly_b', 'NV_1238', 'SiIV_1393', 'CIV_1548', 
    'CIII_1908', 'MgII_2796', 'H_b', 'OIII_5008', 'H_a'
}

# --- Centralized Plot Style Configuration ---
# The span of a region selected with the Select tool (F73): the colour of
# the selection elsewhere (the highlighted group of the System Inspector)
SELECTION_STYLE = {'color': '#f5a100', 'alpha': 0.25, 'linewidth': 0}

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
        'color': '#aaaaaa',   # Grey shading
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
