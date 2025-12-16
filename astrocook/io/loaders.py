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

def detect_file_format(file_path: str) -> str:
    path_lower = file_path.lower()
    if path_lower.endswith(('.acs', '.acs2')):
        return "archive"

    try:
        with fits.open(file_path, mode='readonly', memmap=True) as hdul:
            ext_names = [h.name.upper() for h in hdul]
            
            # 1. PypeIt spec1d (Multi-Extension Order-by-Order)
            # Look for HDU names like 'OBJ...-ORDER...'
            if any('ORDER' in en for en in ext_names):
                return "pypeit_spec1d_multi"

            # 2. PypeIt Coadds/Stacks
            if 'SPECTRUM' in ext_names:
                cols = hdul['SPECTRUM'].columns.names
                if any(c == 'wave_stack' for c in cols):
                    return "pypeit_stack"
                if 'ivar' in cols or 'sigma' in cols:
                    return "pypeit_coadd"

            # 3. ESPRESSO S2D
            if 'SCIDATA' in ext_names and \
               any(w in ext_names for w in ['WAVEDATA_VAC_BARY', 'WAVEDATA_AIR']):
                return "espresso_s2d_fits"

            # 4. X-Shooter Binary Table (Standard ESO)
            instrume = hdul[0].header.get('INSTRUME', '').upper()
            if ('XSHOOTER' in instrume) and 'BinTableHDU' in str(type(hdul[-1])):
                return "xshooter_fits"
            
    except (OSError, fits.VerifyError):
        pass

    # CSV/ASCII checks
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            head = [next(f) for _ in range(10)]
        head_str = "".join(head).upper()
        if "RESVEL" in head_str: return "ascii_resvel_header"
        if path_lower.endswith('.csv'): return "simple_csv"
    except Exception:
        pass

    return "generic_spectrum"

# --- Loader Implementations ---

@register_loader("pypeit_coadd")
def load_pypeit_coadd(file_path: str, **kwargs) -> SpectrumDataV2:
    """
    Loader for PypeIt 'Coadd' 1D files.
    Features: Sorts wavelength and removes negative/bad values.
    """
    try:
        data = Table.read(file_path, hdu='SPECTRUM')
        
        x = data['wave'].data
        y = data['flux'].data
        
        # Error handling: prefer sigma, fallback to IVAR
        if 'sigma' in data.columns:
            dy = data['sigma'].data
        elif 'ivar' in data.columns:
            ivar = data['ivar'].data
            with np.errstate(divide='ignore', invalid='ignore'):
                dy = np.sqrt(1.0 / ivar)
            dy[~np.isfinite(dy)] = np.inf 
        else:
            dy = np.zeros_like(y)
            
        mask_in = data['mask'].data if 'mask' in data.columns else np.zeros_like(y, dtype=int)

        # --- CRITICAL FIX: CLEAN AND SORT ---
        # 1. Filter invalid wavelengths (<= 0) and NaNs
        valid_mask = (x > 1.0) & np.isfinite(x) & np.isfinite(y)
        
        x = x[valid_mask]
        y = y[valid_mask]
        dy = dy[valid_mask]
        mask_in = mask_in[valid_mask]
        
        # 2. Sort by Wavelength
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        y = y[sort_idx]
        dy = dy[sort_idx]
        mask_in = mask_in[sort_idx]
        
        # 3. Units
        x_unit = au.Angstrom
        if np.median(x) > 3000: 
            x = x / 10.0
            x_unit = au.nm

        y_unit = au.dimensionless_unscaled
        
        # 4. Build Structure
        x_col = DataColumnV2(x, x_unit, "Wavelength")
        xmin_col, xmax_col = _auto_limits(x, x_unit)
        aux_cols = {'quality': DataColumnV2(mask_in, au.dimensionless_unscaled, "PypeIt Mask")}

        with fits.open(file_path) as hdul:
            meta = dict(hdul[0].header)
        meta['ORIGIN'] = 'pypeit_coadd'

        return SpectrumDataV2(
            x=x_col, xmin=xmin_col, xmax=xmax_col,
            y=DataColumnV2(y, y_unit), 
            dy=DataColumnV2(dy, y_unit),
            aux_cols=aux_cols,
            meta=meta
        )
    except Exception as e:
        logging.error(f"Failed to load PypeIt Coadd: {e}")
        raise


@register_loader("pypeit_spec1d_multi")
def load_pypeit_spec1d_multi(file_path: str, **kwargs) -> SpectrumDataV2:
    """
    Loader for PypeIt 'spec1d' Multi-Extension FITS.
    Stitches echelle orders (extensions named 'ORDERxxxx') into one structure.
    """
    try:
        x_list, y_list, dy_list, mask_list, order_list = [], [], [], [], []
        
        with fits.open(file_path) as hdul:
            meta = dict(hdul[0].header)
            
            # Iterate over all HDUs
            for hdu in hdul:
                # Identify Order extensions
                if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)) and 'ORDER' in hdu.name:
                    
                    # Extract Order Number from name (e.g. "...ORDER0030")
                    try:
                        order_num = int(re.search(r'ORDER(\d+)', hdu.name).group(1))
                    except:
                        order_num = 0
                        
                    data = hdu.data
                    
                    # Grab Optimal Extraction columns (OPT_*)
                    if 'OPT_WAVE' not in data.names: continue
                    
                    wave = data['OPT_WAVE']
                    flux = data['OPT_FLAM'] # Flux density
                    
                    # Error (SIG preferred)
                    if 'OPT_FLAM_SIG' in data.names:
                        err = data['OPT_FLAM_SIG']
                    elif 'OPT_FLAM_IVAR' in data.names:
                        ivar = data['OPT_FLAM_IVAR']
                        with np.errstate(divide='ignore', invalid='ignore'):
                            err = np.sqrt(1.0 / ivar)
                        err[~np.isfinite(err)] = np.inf
                    else:
                        err = np.zeros_like(flux)
                        
                    mask = data['OPT_MASK'] if 'OPT_MASK' in data.names else np.zeros_like(flux, dtype=int)
                    
                    # Filter valid pixels within the order immediately
                    valid = (wave > 1.0) & np.isfinite(wave) & np.isfinite(flux)
                    
                    x_list.append(wave[valid])
                    y_list.append(flux[valid])
                    dy_list.append(err[valid])
                    mask_list.append(mask[valid])
                    order_list.append(np.full_like(wave[valid], order_num))

        if not x_list:
            raise ValueError("No valid 'ORDER' extensions found in PypeIt spec1d file.")

        # Concatenate all orders
        x = np.concatenate(x_list)
        y = np.concatenate(y_list)
        dy = np.concatenate(dy_list)
        q = np.concatenate(mask_list)
        o = np.concatenate(order_list)

        # Sort the stitched result by wavelength
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        y = y[sort_idx]
        dy = dy[sort_idx]
        q = q[sort_idx]
        o = o[sort_idx]

        # Units
        x_unit = au.Angstrom
        if np.median(x) > 3000:
            x = x / 10.0
            x_unit = au.nm
        
        y_unit = au.dimensionless_unscaled
        
        # Build structure
        x_col = DataColumnV2(x, x_unit, "Wavelength")
        xmin_col, xmax_col = _auto_limits(x, x_unit)
        
        aux_cols = {
            'quality': DataColumnV2(q, au.dimensionless_unscaled, "PypeIt Mask"),
            'order': DataColumnV2(o, au.dimensionless_unscaled, "Echelle Order")
        }

        meta['ORIGIN'] = 'pypeit_spec1d_multi'
        meta['plot_style'] = 'scatter' # Suggest scatter for echelle data

        return SpectrumDataV2(
            x=x_col, xmin=xmin_col, xmax=xmax_col,
            y=DataColumnV2(y, y_unit), 
            dy=DataColumnV2(dy, y_unit),
            aux_cols=aux_cols,
            meta=meta
        )

    except Exception as e:
        logging.error(f"Failed to load PypeIt spec1d Multi-Order: {e}")
        raise

@register_loader("pypeit_stack")
def load_pypeit_stack(file_path: str, **kwargs) -> SpectrumDataV2:
    """
    Loader for PypeIt 'OrderStack' files (Arrays stored in a single table row).
    """
    try:
        # Load the table, accessing the first row immediately
        t = Table.read(file_path, hdu='SPECTRUM')
        row = t[0] # The data is in the cells of the first row

        # 1. Wavelength
        x = row['wave_stack']
        x_unit = au.Angstrom
        if np.median(x) > 3000:
            x = x / 10.0
            x_unit = au.nm

        # 2. Flux
        y = row['flux_stack']
        y_unit = au.dimensionless_unscaled

        # 3. Error
        if 'sigma_stack' in t.columns:
             dy = row['sigma_stack']
        elif 'ivar_stack' in t.columns:
             ivar = row['ivar_stack']
             with np.errstate(divide='ignore', invalid='ignore'):
                 dy = np.sqrt(1.0 / ivar)
             dy[~np.isfinite(dy)] = np.inf
        else:
             dy = np.zeros_like(y)

        # 4. Aux
        aux_cols = {}
        if 'mask_stack' in t.columns:
            aux_cols['quality'] = DataColumnV2(
                row['mask_stack'], au.dimensionless_unscaled, "PypeIt Stack Mask"
            )

        # 5. Metadata
        with fits.open(file_path) as hdul:
            meta = dict(hdul[0].header)
        meta['ORIGIN'] = 'pypeit_stack'

        x_col = DataColumnV2(x, x_unit)
        xmin_col, xmax_col = _auto_limits(x, x_unit)

        return SpectrumDataV2(
            x=x_col, xmin=xmin_col, xmax=xmax_col,
            y=DataColumnV2(y, y_unit), 
            dy=DataColumnV2(dy, y_unit),
            aux_cols=aux_cols,
            meta=meta
        )
    except Exception as e:
        logging.error(f"Failed to load PypeIt Stack: {e}")
        raise

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

@register_loader("espresso_s2d_fits")
def load_espresso_s2d_fits(file_path: str, **kwargs) -> SpectrumDataV2:
    """
    Loader for ESPRESSO S2D (Order-by-Order) FITS files.
    Correctly handles pixel widths by computing them in 2D before flattening.
    """
    try:
        with fits.open(file_path) as hdul:
            header = hdul[0].header
            
            # Helper to safely get data or return None
            def get_data(key):
                if key in hdul: return hdul[key].data
                return None

            wave_data = get_data('WAVEDATA_VAC_BARY')
            if wave_data is None: wave_data = get_data('WAVEDATA_AIR') # Fallback
            
            sci_data = get_data('SCIDATA')
            err_data = get_data('ERRDATA')
            qual_data = get_data('QUALDATA')
            tell_data = get_data('TELLURIC_CORRECTION') 

            if wave_data is None or sci_data is None:
                raise ValueError("Essential extensions (WAVEDATA, SCIDATA) missing.")

            # --- CRITICAL FIX START ---
            # 1. Compute pixel widths and limits in 2D (per order)
            # This preserves the physical dispersion of the instrument.
            pixel_widths = np.gradient(wave_data, axis=1)
            
            # Divide Flux and Error by dlambda to get Flux Density
            with np.errstate(divide='ignore', invalid='ignore'):
                sci_data = sci_data / pixel_widths
                if err_data is not None:
                    err_data = err_data / pixel_widths
            
            # Calculate physical edges (2D)
            half_widths = np.abs(pixel_widths) / 2.0
            xmin_data = wave_data - half_widths
            xmax_data = wave_data + half_widths
            # --- CRITICAL FIX END ---

            # 2. Create 'Order' Index Array
            n_orders, n_pixels = wave_data.shape
            order_idx = np.repeat(np.arange(n_orders), n_pixels).reshape(n_orders, n_pixels)

            # 3. Flatten (Ravel) 2D arrays into 1D
            x = np.ravel(wave_data)/10
            y = np.ravel(sci_data)
            dy = np.ravel(err_data) if err_data is not None else np.zeros_like(y)
            q = np.ravel(qual_data) if qual_data is not None else np.zeros_like(y, dtype=int)
            o = np.ravel(order_idx)
            t = np.ravel(tell_data) if tell_data is not None else None
            
            # Flatten the pre-computed limits (convert to nm)
            xmin = np.ravel(xmin_data)/10
            xmax = np.ravel(xmax_data)/10

            # 4. Filter Valid Wavelengths
            valid_wave = (x > 0)
            x, y, dy, q, o = x[valid_wave], y[valid_wave], dy[valid_wave], q[valid_wave], o[valid_wave]
            xmin, xmax = xmin[valid_wave], xmax[valid_wave]
            if t is not None: t = t[valid_wave]

            # 5. Sort by Wavelength
            # We sort everything to keep V2 algorithms happy, but we carry 
            # the physically correct xmin/xmax with us.
            sort_idx = np.argsort(x)
            x = x[sort_idx]
            y = y[sort_idx]
            dy = dy[sort_idx]
            q = q[sort_idx]
            o = o[sort_idx]
            xmin = xmin[sort_idx]
            xmax = xmax[sort_idx]
            if t is not None: t = t[sort_idx]

            # 6. Apply Quality Cut
            good_quality = (q < 1)
            x, y, dy, q, o = x[good_quality], y[good_quality], dy[good_quality], q[good_quality], o[good_quality]
            xmin, xmax = xmin[good_quality], xmax[good_quality]
            if t is not None: t = t[good_quality]

            # 7. Units and Columns
            x_unit = au.nm
            y_unit = au.dimensionless_unscaled

            # DO NOT call _auto_limits(x). Use our pre-computed columns.
            xmin_col = DataColumnV2(xmin, x_unit)
            xmax_col = DataColumnV2(xmax, x_unit)

            aux_cols = {
                'quality': DataColumnV2(q, au.dimensionless_unscaled, "Quality Flag"),
                'order': DataColumnV2(o, au.dimensionless_unscaled, "Echelle Order Index")
            }
            if t is not None:
                aux_cols['telluric_corr'] = DataColumnV2(t, au.dimensionless_unscaled, "Telluric Correction")

            meta = dict(header)
            meta['ORIGIN'] = 'espresso_s2d_loader'
            meta['plot_style'] = 'scatter' 

            return SpectrumDataV2(
                x=DataColumnV2(x, x_unit), xmin=xmin_col, xmax=xmax_col,
                y=DataColumnV2(y, y_unit), dy=DataColumnV2(dy, y_unit),
                aux_cols=aux_cols,
                meta=meta
            )

    except Exception as e:
        logging.error(f"Failed to load ESPRESSO S2D: {e}")
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

@register_loader("generic_spectrum")
def load_generic_spectrum(file_path: str, **kwargs) -> SpectrumDataV2:
    """Fallback loader."""
    # ... (Keep your existing fallback logic here) ...
    raise NotImplementedError("Generic loader not fully implemented in this snippet.")