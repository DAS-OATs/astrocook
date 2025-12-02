import logging
import numpy as np
from astropy.table import Table
import astropy.units as au
from typing import Callable, Dict, Optional

# Import your V2 structures
from astrocook.core.structures import SpectrumDataV2, DataColumnV2 # Adjust import path as needed

# --- The Registry ---
_LOADER_REGISTRY: Dict[str, Callable] = {}

def register_loader(format_name: str):
    """Decorator to register a new loader function."""
    def decorator(func):
        _LOADER_REGISTRY[format_name] = func
        return func
    return decorator

def get_loader(format_name: str) -> Optional[Callable]:
    return _LOADER_REGISTRY.get(format_name)

# --- Loader Implementations ---

@register_loader("simple_csv")
def load_simple_csv(file_path: str, **kwargs) -> SpectrumDataV2:
    """
    Example V2 Loader: Reads a 3-column CSV (wave, flux, err).
    """
    try:
        data = Table.read(file_path, format='ascii.csv')
        
        # 1. Extract raw numpy arrays
        x = data['col1'].data
        y = data['col2'].data
        dy = data['col3'].data
        
        # 2. Define Units (You can pass these in kwargs if dynamic)
        x_unit = au.nm
        y_unit = au.erg / (au.cm**2 * au.s * au.Angstrom)
        
        # 3. Create Data Columns
        x_col = DataColumnV2(x, x_unit, "Wavelength")
        y_col = DataColumnV2(y, y_unit, "Flux")
        dy_col = DataColumnV2(dy, y_unit, "Error")
        
        # 4. Handle limits (can be computed if missing)
        # Simple approximation: midpoints
        xmin = x - (np.diff(x, append=x[-1])/2)
        xmax = x + (np.diff(x, append=x[-1])/2)
        xmin_col = DataColumnV2(xmin, x_unit)
        xmax_col = DataColumnV2(xmax, x_unit)

        # 5. Metadata
        meta = {'ORIGIN': 'simple_csv_loader', 'FILE': file_path}
        
        # 6. Return the Data Core
        return SpectrumDataV2(
            x=x_col, xmin=xmin_col, xmax=xmax_col,
            y=y_col, dy=dy_col,
            aux_cols={},
            meta=meta
        )
        
    except Exception as e:
        logging.error(f"Failed to load CSV: {e}")
        raise

@register_loader("xshooter_fits")
def load_new_user_format(file_path: str, **kwargs) -> SpectrumDataV2:
    """
    Placeholder for your user's new format.
    """
    # ... logic to open FITS, read specific extensions, handle headers ...
    pass