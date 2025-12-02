import logging
import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.units as au
from typing import Callable, Dict, Optional, Tuple

# Import V2 structures from the core
from astrocook.core.structures import SpectrumDataV2, DataColumnV2

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

# --- Helpers ---

def _auto_limits(x_arr: np.ndarray, x_unit: au.Unit) -> Tuple[DataColumnV2, DataColumnV2]:
    """
    Calculates xmin/xmax (pixel edges) from centers, assuming contiguous pixels.
    """
    # Calculate difference between pixels (gradient)
    # Prepend the first diff to handle the start edge
    if len(x_arr) > 1:
        diff = np.diff(x_arr, prepend=x_arr[0] - (x_arr[1] - x_arr[0]))
    else:
        diff = np.array([1.0]) # Fallback for single-pixel artifacts

    # Assume center is midpoint
    half_width = diff / 2.0
    
    xmin = x_arr - half_width
    xmax = x_arr + half_width
    
    return DataColumnV2(xmin, x_unit), DataColumnV2(xmax, x_unit)

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
        
        # 4. Use helper for limits
        xmin_col, xmax_col = _auto_limits(x, x_unit)

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
def load_xshooter_fits(file_path: str, **kwargs) -> SpectrumDataV2:
    """
    Loader for ESO X-Shooter 1D spectra (Binary Table FITS).
    Expects columns: WAVE, FLUX, ERR, QUAL (optional).
    """
    try:
        with fits.open(file_path) as hdul:
            # X-Shooter data is usually in extension 1
            data = hdul[1].data
            header = hdul[0].header
            
            # 1. Extract Arrays (Handle standard column naming variations)
            cols = data.columns.names
            
            # Wavelength
            if 'WAVE' in cols: x = data['WAVE']
            elif 'LAMBDA' in cols: x = data['LAMBDA']
            else: raise ValueError("Could not find wavelength column (WAVE/LAMBDA)")
            
            # Flux
            if 'FLUX' in cols: y = data['FLUX']
            else: raise ValueError("Could not find flux column (FLUX)")
            
            # Error
            if 'ERR' in cols: dy = data['ERR']
            elif 'ERROR' in cols: dy = data['ERROR']
            else: dy = np.zeros_like(y)
            
            # 2. Units (Hardcoded for X-Shooter standard, or parsed from header)
            x_unit = au.nm 
            y_unit = au.erg / (au.cm**2 * au.s * au.Angstrom)
            
            # 3. Build Columns
            x_col = DataColumnV2(x, x_unit)
            y_col = DataColumnV2(y, y_unit)
            dy_col = DataColumnV2(dy, y_unit)
            xmin_col, xmax_col = _auto_limits(x, x_unit)
            
            # 4. Handle Quality Flag (Aux Column)
            aux_cols = {}
            if 'QUAL' in cols:
                # Convert Bad Pixel Map to internal mask style if needed
                qual = data['QUAL']
                aux_cols['quality'] = DataColumnV2(qual, au.dimensionless_unscaled, "Quality Flag")

            # 5. Metadata
            meta = dict(header)
            meta['ORIGIN'] = 'xshooter_fits_loader'
            
            return SpectrumDataV2(
                x=x_col, xmin=xmin_col, xmax=xmax_col,
                y=y_col, dy=dy_col,
                aux_cols=aux_cols,
                meta=meta
            )

    except Exception as e:
        logging.error(f"Failed to load X-Shooter FITS: {e}")
        raise