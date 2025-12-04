from astropy.io import fits
from astropy.table import Table
import astropy.units as au
import logging
import numpy as np
import re
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


@register_loader("ascii_resvel_header")
def load_ascii_with_resvel_header(file_path: str, **kwargs) -> SpectrumDataV2:
    """
    Robust loader for ASCII files with 'RESVEL' headers and inline source tags.
    Handles files where data rows might be split across lines or mixed with tags.
    """
    try:
        meta = {'ORIGIN': 'ascii_resvel_loader', 'FILE': file_path}
        resol_values = []
        all_values = []

        # Use utf-8-sig to handle potential Byte Order Marks (BOM)
        with open(file_path, 'r', encoding='utf-8-sig', errors='replace') as f:
            for line in f:
                line = line.strip()
                if not line: continue

                # 1. Extract RESVEL metadata (and skip line if it's purely metadata)
                if 'RESVEL' in line:
                    match = re.search(r'RESVEL\s*[:=]?\s*([0-9\.]+)', line)
                    if match:
                        resol_values.append(float(match.group(1)))
                    # We assume lines with RESVEL don't contain spectral data points
                    continue
                
                # 2. Remove tags which might appear mid-stream
                # Regex removes "" but leaves the numbers around it
                line_clean = re.sub(r'\\', '', line)

                # 3. Extract numbers
                # This handles cases where X is on one line and Y DY are on the next
                try:
                    # Split by whitespace and convert to float
                    # We iterate so we can catch non-numeric garbage without losing the whole line
                    for token in line_clean.split():
                        try:
                            all_values.append(float(token))
                        except ValueError:
                            pass # Ignore non-numeric tokens
                except Exception:
                    continue

        # 4. Reshape Data
        # We expect a multiple of 3 (Wave, Flux, Err)
        total_points = len(all_values)
        if total_points == 0:
            raise ValueError(f"No numeric data found in {file_path}")
            
        if total_points % 3 != 0:
            logging.warning(f"Data length ({total_points}) is not divisible by 3. Truncating excess.")
            valid_len = (total_points // 3) * 3
            all_values = all_values[:valid_len]

        data = np.array(all_values).reshape(-1, 3)

        x = data[:, 0]
        y = data[:, 1]
        dy = data[:, 2]

        # 5. Units (Assumption: Angstrom for X, dimensionless for Y)
        # Heuristic: If median > 500, it's almost certainly Angstroms (UV/Optical)
        # Lya is 1215 A (121.5 nm). 
        if np.median(x) > 500:
            logging.info("Detected Angstroms (median > 500). Converting to nm.")
            x = x / 10.0
            x_unit = au.nm
        else:
            logging.info("Assuming nm (median < 500).")
            x_unit = au.nm

        y_unit = au.dimensionless_unscaled
        
        # 6. Build Columns
        x_col = DataColumnV2(x, x_unit, "Wavelength")
        y_col = DataColumnV2(y, y_unit, "Flux")
        dy_col = DataColumnV2(dy, y_unit, "Error")
        xmin_col, xmax_col = _auto_limits(x, x_unit)

        # 7. Handle Resolution
        aux_cols = {}
        if resol_values:
            # Use the median if multiple RESVELs found (e.g. multiple sources stitched in file)
            final_resol = np.median(resol_values)
            c_kms = 299792.458
            meta['resol_fwhm'] = final_resol
            meta['resol'] = c_kms / final_resol
            
            resol_arr = np.full_like(x, final_resol)
            aux_cols['resol'] = DataColumnV2(resol_arr, au.km/au.s, "Resolution (FWHM km/s)")
            logging.info(f"Loaded {len(x)} points. Resolution FWHM: {final_resol} km/s (R~{meta['resol']:.0f})")

        return SpectrumDataV2(
            x=x_col, xmin=xmin_col, xmax=xmax_col,
            y=y_col, dy=dy_col,
            aux_cols=aux_cols,
            meta=meta
        )

    except Exception as e:
        logging.error(f"Failed to load ASCII Resvel file: {e}")
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