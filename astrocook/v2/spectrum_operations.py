import bisect
from astropy import units as au
from astropy.constants import c as c_light
from astropy.stats import sigma_clip
import logging
import numpy as np
from numpy.polynomial.polynomial import polyfit
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d, label
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter
from scipy.special import erf, erfinv
from typing import Dict, List, Optional, Tuple

from .atomic_data import STANDARD_MULTIPLETS, xem_d
from ..v1.functions import x_convert as x_convert_v1 
from ..v2.structures import DataColumnV2, SpectrumDataV2 

# We assume 'xem' (the Ly-alpha rest wavelength) is accessible here for the custom equivalency
# You must ensure this global is either imported or defined here (Technical Debt).
XEM_LYA = 121.567 * au.nm 
XEM_LYA_KM_S_ZERO = (XEM_LYA * (1 + 0.0)).to(au.km/au.s, equivalencies=au.equivalencies.doppler_optical(XEM_LYA)).value

def convert_axis_velocity(
    x_quantity: au.Quantity, 
    z_rf: float, # Renamed from zem
    target_unit: au.Unit
) -> au.Quantity:
    """
    Converts X-axis data between wavelength (nm) and velocity (km/s) space,
    using the V1 zero-point logic (Ly-alpha rest-frame relative to z_rf).
    """
    
    v_c = c_light.to(au.km/au.s).value
    
    # Define a custom conversion function
    # Lambda to Velocity (Forward)
    def lambda_to_v(lambda_val):
        # Use z_rf as the zero-point for the velocity conversion
        return np.log(lambda_val / (XEM_LYA.value * (1 + z_rf))) * v_c

    # Velocity to Lambda (Reverse)
    def v_to_lambda(v_val):
        return np.exp(v_val / v_c) * (XEM_LYA.value * (1 + z_rf))

    custom_equiv = [
        (au.nm, au.km/au.s, lambda_to_v, v_to_lambda),
        (au.Angstrom, au.km/au.s, lambda_to_v, v_to_lambda) # Add Angstrom
    ]
    
    if target_unit.is_equivalent(au.km/au.s) or target_unit.is_equivalent(au.nm):
        return x_quantity.to(target_unit, equivalencies=custom_equiv)
    else:
        return x_quantity.to(target_unit)

def convert_x_axis(spec_data: SpectrumDataV2, z_rf: float, xunit: au.Unit) -> Tuple[au.Quantity, au.Quantity, au.Quantity]:
    """
    Performs X-axis conversion based on V1 logic and returns new Quantity arrays.
    """
    
    x = spec_data.x.quantity
    xmin = spec_data.xmin.quantity
    xmax = spec_data.xmax.quantity
    
    if xunit.is_equivalent(au.km/au.s) or x.unit.is_equivalent(au.km/au.s):
        # Use the V2 dedicated velocity conversion logic
        x_new = convert_axis_velocity(x, z_rf, xunit)
        xmin_new = convert_axis_velocity(xmin, z_rf, xunit)
        xmax_new = convert_axis_velocity(xmax, z_rf, xunit)
        
    else:
        # Simple wavelength-to-wavelength
        x_new = x.to(xunit)
        xmin_new = xmin.to(xunit)
        xmax_new = xmax.to(xunit)
    
    return x_new, xmin_new, xmax_new

def convert_y_axis(spec_data: SpectrumDataV2, yunit: au.Unit) -> Tuple[DataColumnV2, DataColumnV2]:
    """
    Performs Y-axis unit conversion on flux and error data using native Astropy methods.
    """
    
    # Get Quantities from V2 Data
    y = spec_data.y.quantity
    dy = spec_data.dy.quantity

    # Apply standard Astropy conversion
    y_new = y.to(yunit) 
    dy_new = dy.to(yunit) 
    
    # Create new DataColumnV2 instances
    new_y = DataColumnV2(y_new.value, yunit)
    new_dy = DataColumnV2(dy_new.value, yunit)

    new_aux_cols = {}
    original_y_unit = spec_data.y.unit # The unit of the data before conversion

    # Iterate over all auxiliary columns
    for name, data_col in spec_data.aux_cols.items():
        
        # Check if the column is flux-like (shares the same unit as the main Y-axis)
        if data_col.unit.is_equivalent(original_y_unit):
            
            # Apply standard Astropy conversion
            converted_quantity = data_col.quantity.to(yunit)
            
            # Create a new DataColumnV2 entry
            new_aux_cols[name] = DataColumnV2(
                converted_quantity.value, 
                yunit, 
                description=data_col.description
            )
        else:
            # If the column unit is different (e.g., quality flags, resolution), copy it directly
            new_aux_cols[name] = data_col

    return new_y, new_dy, new_aux_cols

def smooth_spectrum(
    x: au.Quantity,
    y: np.ndarray,
    sigma_kms: float, # std in km/s
    z_rf: float = 0.0
) -> np.ndarray:
    """
    Smooths a flux array using a Gaussian filter with a given std in km/s.
    """
    
    # 1. Smooth the result (V1's _gauss_convolve logic)
    if sigma_kms > 0:
        # 2. Convert x-axis to km/s to get pixel size
        try:
            # We need the x-axis as a Quantity to convert it
            x_kms = convert_axis_velocity(x, z_rf=z_rf, target_unit=au.km/au.s).value
            avg_dv = np.mean(np.gradient(x_kms))
            
            # 3. Convert smoothing std (km/s) to pixels
            smooth_std_pix = sigma_kms / avg_dv
            
            # 4. Apply the filter
            y_final = gaussian_filter1d(y, smooth_std_pix)
        except Exception as e:
            logging.error(f"Failed to calculate smoothing in km/s: {e}. Returning unsmoothed array.")
            return y
    else:
        # No smoothing requested
        y_final = y
        
    return y_final

def rebin_spectrum(
        x_in: au.Quantity, xmin_in: au.Quantity, xmax_in: au.Quantity,
        spec_data: SpectrumDataV2, # Still pass the full data for y/dy/aux_cols
        xstart: Optional[au.Quantity], xend: Optional[au.Quantity], 
        dx: au.Quantity, calc_unit: au.Unit, # Pass the unit used for calculation
        kappa: Optional[float], filling: float
    ) -> Tuple[au.Quantity, au.Quantity, au.Quantity, au.Quantity, au.Quantity]:
    """
    Pure function to rebin spectrum data. MUST contain the logic from Spectrum._rebin (V1).
    Returns new arrays for x, xmin, xmax, y, dy.
    """
    # --- V1 Pre-Processing and Initialization ---
    
    # NOTE: In V1, the spectrum was converted in-place before rebinning. 
    # In V2, we work with the current units but create the output in the target dx unit.
    
    # Target unit is derived from dx, assuming its value represents the intended unit (e.g., km/s)
    xunit_target = dx.unit 
        
    # X-axis inputs (Already converted to calc_unit by SpectrumV2.rebin)
    x_value = x_in.value
    xmin_value = xmin_in.value
    xmax_value = xmax_in.value

    # Y-axis inputs (Extracted from the V2 Data Core)
    y_value = spec_data.y.quantity.value
    dy_value = spec_data.dy.quantity.value
    y_unit = spec_data.y.quantity.unit
    
    # Calculation unit is derived from the dx input
    dx_value = dx.value 

    # 1. Handling xstart
    if xstart is None:
        # Use the minimum of the input data (which is already in x_value, the calculation space)
        xstart_value = np.nanmin(x_value) 
    else:
        # If provided, convert the input Quantity (xstart) to the NumPy value in the calculation unit.
        # Note: We rely on the upstream call ensuring xstart is convertible to calc_unit.
        xstart_value = xstart.to(xunit_target, equivalencies=au.equivalencies.spectral()).value

    # 2. Handling xend
    if xend is None:
        # Use the maximum of the input data
        xend_value = np.nanmax(x_value)
    else:
        # Convert the input Quantity (xend) to the NumPy value in the calculation unit.
        xend_value = xend.to(xunit_target, equivalencies=au.equivalencies.spectral()).value

    # Create new bin centers (x_out) and boundaries
    x_out_value = np.arange(xstart_value, xend_value, dx_value)
    
    # Helper to create xmin/xmax for the new array (V1 logic simplified)
    mean_val = 0.5 * (x_out_value[1:] + x_out_value[:-1])
    xmin_out_value = np.append(x_out_value[0] - dx_value / 2, mean_val)
    xmax_out_value = np.append(mean_val, x_out_value[-1] + dx_value / 2)
    x_out_value = 0.5 * (xmin_out_value + xmax_out_value)

    # Initialize output arrays
    y_out_value = np.empty_like(x_out_value)
    dy_out_value = np.empty_like(x_out_value)
    
    # --- V1 Core Rebinning Loop (Weighted Average) ---
    
    for i, (m, M) in enumerate(zip(xmin_out_value, xmax_out_value)):
        # Find all input bins that overlap with the new output bin [m, M]
        im = bisect.bisect_left(xmax_value, m)
        iM = bisect.bisect_right(xmin_value, M)
        
        # Calculate the fractional overlap for each input bin
        frac = (np.minimum(M, xmax_value[im:iM]) - 
                np.maximum(m, xmin_value[im:iM])) / dx_value
        
        # Select data and error arrays for the overlapping region
        ysel = y_value[im:iM]
        dysel = dy_value[im:iM]
        
        # Filter for non-NaN input flux data
        not_nan_w = ~np.isnan(ysel)
        ysel = ysel[not_nan_w]
        dysel = dysel[not_nan_w]
        frac = frac[not_nan_w]
        
        # Optional Kappa-Sigma Clipping (V1 logic)
        if kappa is not None and len(ysel) > 0:
            clipped = sigma_clip(ysel, sigma=kappa, masked=True)
            # Find indices where data was not masked
            not_masked_w = ~clipped.mask
            
            if np.sum(not_masked_w) > 0:
                ysel = ysel[not_masked_w]
                dysel = dysel[not_masked_w]
                frac = frac[not_masked_w]

        # Final weighted average calculation
        if len(frac) > 0:
            # V1 logic: If dy=0 or is missing, use frac as weight (simple average)
            # If dy exists, use 1/dy^2 * frac as weight (weighted average)
            
            # Filter out zero errors (which lead to infinite weights)
            safe_error_w = (dysel != 0.0) & (~np.isnan(dysel))
            
            if np.any(safe_error_w):
                # Weighted average using variance (1/dy^2) and fractional overlap
                weights = frac[safe_error_w] / dysel[safe_error_w]**2
                y_out_value[i] = np.average(ysel[safe_error_w], weights=weights)
                
                # Propagate error (weighted variance calculation - V1 logic simplified)
                # V1 used a simplified error propagation that summed weights:
                sum_weights = np.nansum(weights)
                if sum_weights > 0:
                    dy_out_value[i] = np.sqrt(np.nansum(weights**2 * dysel[safe_error_w]**2)) / sum_weights
                else:
                    dy_out_value[i] = filling
            else:
                # Simple average weighted by fractional overlap
                y_out_value[i] = np.average(ysel, weights=frac)
                dy_out_value[i] = filling # Cannot calculate error without dy_in
        else:
            # No data in this bin
            y_out_value[i] = filling
            dy_out_value[i] = filling
            
    # --- Final V2 Output ---
    
    # Convert final NumPy arrays back to Quantities
    x_out = au.Quantity(x_out_value, xunit_target)
    xmin_out = au.Quantity(xmin_out_value, xunit_target)
    xmax_out = au.Quantity(xmax_out_value, xunit_target)
    y_out = au.Quantity(y_out_value, y_unit)
    dy_out = au.Quantity(dy_out_value, y_unit)
    
    # Return the new immutable data required by SpectrumV2
    return x_out, xmin_out, xmax_out, y_out, dy_out

def running_mean(y: np.ndarray, h: int) -> np.ndarray:
    """
    Computes a running mean of an array.
    """
    window_size = 2 * h + 1
    if window_size > len(y):
        logging.warning(f"Running mean window ({window_size}) is larger than array ({len(y)}), returning array mean.")
        return np.full_like(y, np.mean(y))
        
    # Use 'same' mode to keep array length, symmetric padding
    return np.convolve(y, np.ones(window_size) / window_size, mode='same')

def running_std(y: np.ndarray, h: int) -> np.ndarray:
    """
    Computes a running standard deviation (RMS with respect to the local mean).
    Uses the identity: var = <y^2> - <y>^2
    """
    window_size = 2 * h + 1
    if window_size > len(y):
        logging.warning(f"Window ({window_size}) larger than array ({len(y)}), returning overall std.")
        return np.full_like(y, np.nanstd(y))
        
    kernel = np.ones(window_size) / window_size
    
    # 1. Calculate running mean of y: <y>
    mean_y = np.convolve(y, kernel, mode='same')
    
    # 2. Calculate running mean of y squared: <y^2>
    mean_y_sq = np.convolve(y**2, kernel, mode='same')
    
    # 3. Calculate variance: <y^2> - <y>^2
    running_var = mean_y_sq - mean_y**2
    
    # 4. Clip small negative values due to floating-point errors
    running_var[running_var < 0] = 0
    
    return np.sqrt(running_var)

def find_absorbed_regions(
    x: au.Quantity,
    y: np.ndarray, 
    dy: np.ndarray, 
    z_em: float,
    smooth_len_lya: au.Quantity, 
    smooth_len_out: au.Quantity, 
    kappa: float, 
    template: bool,
    maxiter: int = 100
) -> np.ndarray:
    """
    Pure function to find absorbed regions, porting the V1 'clip_flux' logic.
    
    Returns:
        mask_abs (np.ndarray): A boolean mask where True means "unabsorbed".
    """
    
    # --- 1. Initial Masking ---
    mask_unabs = np.ones_like(y, dtype=bool)
    valid_data = ~np.isnan(y) & ~np.isnan(dy) & (dy > 0)
    mask_unabs[~valid_data] = False
    mask_unabs[y <= 3.0 * dy] = False # Low S/N
    
    if z_em == 0.0:
        logging.warning("z_em is 0.0. Cannot split at Ly-alpha. Using 'smooth_len_out' for all regions.")
        lya_prox_nm = -np.inf * au.nm
    else:
        # Calculate Ly-alpha proximity region (V1 logic)
        prox = 5000 * au.km/au.s
        lya_obs = XEM_LYA * (1 + z_em)
        lya_prox_nm = (lya_obs * (1 - prox / c_light)).to(au.nm)

    # --- 2. Handle Template (Stubbed) ---
    y_proc = y.copy()
    if template:
        # V1 logic: y_proc = y / template_interp
        logging.warning("Template correction is not yet implemented. Skipping.")
        # We would need to pass in the qso_composite data
        pass

    # --- 3. Calculate pixel size in km/s (dv) ---
    x_kms = convert_axis_velocity(x, z_rf=0.0, target_unit=au.km/au.s).value
    dv = np.gradient(x_kms)

    # --- 4. Define Lya-forest and outside regions ---
    # We must operate on the *full* array, not a slice
    in_lya_forest = (x < lya_prox_nm) & mask_unabs
    out_of_forest = (x >= lya_prox_nm) & mask_unabs
    
    intervals = [
        (in_lya_forest, smooth_len_lya),
        (out_of_forest, smooth_len_out)
    ]
    
    mask_clipped = np.zeros_like(y, dtype=bool) # Mask of *newly* clipped points

    # --- 5. Iterative Clipping Loop ---
    for i in range(maxiter):
        
        y_rm_full = np.full_like(y, np.nan)
        
        # --- 6. Process Lya and Outside regions separately ---
        for mask_interval, smooth_len in intervals:
            if np.sum(mask_interval) == 0:
                continue # Skip if this region is empty
                
            # Combine current interval mask with *not-yet-clipped* mask
            current_mask = mask_interval & ~mask_clipped
            if np.sum(current_mask) < 10: # Safety check
                continue
                
            y_interval = y_proc[current_mask]
            
            # Convert smoothing length (km/s) to pixel window
            avg_dv = np.mean(dv[current_mask])
            hwindow = int(smooth_len.to_value(au.km/au.s) / avg_dv / 8) # V1 logic: h = len/dv/8
            hwindow = max(1, hwindow) # Ensure window is at least 1
            
            if len(y_interval) < hwindow * 2 + 1:
                y_rm_interval = np.full_like(y_interval, np.mean(y_interval))
            else:
                y_rm_interval = running_mean(y_interval, h=hwindow)
            
            # Interpolate this mean back to the interval's indices
            indices_full = np.arange(len(y))
            indices_interval = indices_full[current_mask]
            y_rm_interp = np.interp(indices_full[mask_interval], indices_interval, y_rm_interval)
            
            # Store the result
            y_rm_full[mask_interval] = y_rm_interp

        # 7. Find new outliers across the *whole* spectrum
        is_outlier = (y_proc - y_rm_full) < (-kappa * dy)
        is_outlier[~mask_unabs] = False # Don't clip already-masked points
        
        if np.sum(is_outlier) == 0:
            logging.info(f"Clipping converged after {i+1} iterations.")
            break
        
        mask_clipped = mask_clipped | is_outlier # Add new outliers to the clipped mask
        mask_unabs = mask_unabs & ~is_outlier  # Update the master unabsorbed mask
        
        if i == maxiter - 1:
            logging.warning(f"Clipping did not converge after {maxiter} iterations.")

    return ~mask_unabs

def fit_continuum_interp(
    x: np.ndarray,
    y: np.ndarray,
    mask_abs: np.ndarray,
    fudge: float = 1.0,
    smooth_std_kms: float = 500.0, # std in km/s
    z_rf: float = 0.0
) -> np.ndarray:
    """
    Fits a continuum by interpolating unabsorbed points and smoothing.
    Based on V1 'clip_flux' logic.
    """
    
    mask_unabs = ~mask_abs

    # 1. Get unabsorbed points
    x_unabs = x[mask_unabs]
    y_unabs = y[mask_unabs]
    
    if len(x_unabs) < 2:
        raise ValueError("Cannot fit continuum: < 2 unabsorbed points found.")
        
    # 2. Interpolate unabsorbed points back to the full grid
    cont_interp = np.interp(x, x_unabs, y_unabs)
    
    # 3. Apply fudge factor
    cont_interp *= fudge
    
    # 4. Smooth the result (V1's _gauss_convolve)
    if smooth_std_kms > 0:
        # Convert x-axis to km/s to get pixel size
        x_kms = convert_axis_velocity(x * au.nm, z_rf=z_rf, target_unit=au.km/au.s).value
        avg_dv = np.mean(np.gradient(x_kms))
        
        # Convert smoothing std (km/s) to pixels
        smooth_std_pix = smooth_std_kms / avg_dv
        
        cont_final = gaussian_filter1d(cont_interp, smooth_std_pix)
    else:
        cont_final = cont_interp
        
    return cont_final

def fit_powerlaw_to_regions(
    x: np.ndarray,
    y: np.ndarray,
    dy: np.ndarray, # <-- ADD dy
    z_em: float,
    regions_str: str,
    kappa: float = 3.0, # <-- ADD kappa
    max_iter: int = 10  # Safety limit for iterations
) -> np.ndarray:
    """
    Fits a power-law continuum (linear in log-log space) to specified
    rest-frame regions.
    """
    
    # 1. Parse regions and select initial data
    all_x = []
    all_y = []
    all_dy = [] # We need dy for clipping
    
    try:
        for region in regions_str.split(','):
            if '-' not in region: continue
            rf_min, rf_max = map(float, region.split('-'))
            obs_min = rf_min * (1.0 + z_em)
            obs_max = rf_max * (1.0 + z_em)
            
            # Select valid, positive data within region
            indices = np.where((x >= obs_min) & (x <= obs_max) & (y > 0) & np.isfinite(y))[0]
            if len(indices) > 0:
                all_x.append(x[indices])
                all_y.append(y[indices])
                all_dy.append(dy[indices])
        
        if not all_x:
            raise ValueError("No valid data points found in specified regions.")
            
        x_sel = np.concatenate(all_x)
        y_sel = np.concatenate(all_y)
        dy_sel = np.concatenate(all_dy)
        
        # 2. Iterative Fit and Clip Loop
        mask = np.ones(len(x_sel), dtype=bool) # Start with all selected points
        intercept, slope = 0.0, 0.0

        for i in range(max_iter):
            if np.sum(mask) < 2:
                raise ValueError("Kappa-sigma clipping removed too many points to fit.")

            # A. Fit in log-log space using ONLY currently unmasked points
            x_fit = np.log10(x_sel[mask])
            y_fit = np.log10(y_sel[mask])
            intercept, slope = polyfit(x_fit, y_fit, 1)
            
            # B. Evaluate model on ALL selected points to check for outliers
            y_model_sel = (10**intercept) * (x_sel**slope)
            
            # C. Identify outliers: data is below model by more than kappa * sigma
            # We primarily care about clipping absorption (negative outliers)
            resid = y_sel - y_model_sel
            is_outlier = resid < (-kappa * dy_sel)
            
            # D. Update mask
            n_old = np.sum(mask)
            mask = mask & ~is_outlier # Keep only points that were good AND are not new outliers
            n_new = np.sum(mask)
            
            logging.debug(f"Power-law iter {i+1}: slope={slope:.3f}, clipped {n_old - n_new} points.")
            
            if n_new == n_old:
                break # Converged

        logging.info(f"Power-law fit converged: slope={slope:.3f}, intercept={intercept:.3f}, used {np.sum(mask)}/{len(x_sel)} points.")
        
        # 3. Generate final model over the FULL x-axis
        y_pl = (10**intercept) * (x**slope)
        return y_pl
        
    except Exception as e:
        logging.error(f"Failed to fit power-law: {e}")
        raise

def detect_regions(
    mask: np.ndarray,
    min_pix: int = 3
) -> Tuple[np.ndarray, int]:
    """
    Detects contiguous regions in a boolean mask (where True = region to detect).
    Returns an integer mask (0=background, 1, 2, 3...=regions) and the count.
    """
    # 1. Label contiguous regions
    # scipy.ndimage.label treats 0/False as background, non-zero as objects.
    region_map, num_regions = label(mask)
    
    # 2. Filter by minimum size
    if min_pix > 1 and num_regions > 0:
        # bincount gives the size of each region (index 0 is background)
        sizes = np.bincount(region_map.ravel())
        
        # Find labels that are too small
        # (Start from 1 because 0 is background)
        bad_labels = np.where((sizes < min_pix) & (np.arange(len(sizes)) > 0))[0]
        
        if len(bad_labels) > 0:
            # Set pixels of bad regions back to 0
            # np.isin is fast for checking multiple values
            mask_to_remove = np.isin(region_map, bad_labels)
            region_map[mask_to_remove] = 0
            
            # Re-label to ensure consecutive IDs (optional but nice)
            # (If we don't re-label, we might have gaps like 1, 3, 4...)
            # A quick way to re-label is just to run label() again on the cleaned mask
            region_map, num_regions = label(region_map > 0)

    return region_map, num_regions

def merge_regions_by_velocity(
    region_map: np.ndarray,
    x: au.Quantity,
    z_rf: float,
    merge_dv: float
) -> Tuple[np.ndarray, int]:
    """
    Merges contiguous absorption regions separated by less than a given
    velocity (dv) threshold.
    """
    if merge_dv <= 0:
        # No merging requested
        num_regions = len(np.unique(region_map)) - 1
        return region_map, num_regions

    try:
        # 1. Get velocity axis
        x_kms = convert_axis_velocity(x, z_rf=z_rf, target_unit=au.km/au.s).value
    except Exception as e:
        logging.error(f"Cannot merge regions: failed to convert x-axis to km/s: {e}")
        num_regions = len(np.unique(region_map)) - 1
        return region_map, num_regions # Return original

    # 2. Find the start/end *indices* of each region
    # Find where the region ID changes
    change_indices = np.where(np.diff(region_map) != 0)[0] + 1
    
    # Add start (0) and end (len-1)
    region_boundaries = np.unique(np.concatenate(([0], change_indices, [len(region_map)])))
    
    if len(region_boundaries) < 2:
        # No regions or only one region
        num_regions = len(np.unique(region_map)) - 1
        return region_map, num_regions

    # 3. Store region properties: {id: (v_min, v_max, [indices])}
    region_info = {}
    
    for i in range(len(region_boundaries) - 1):
        idx_start = region_boundaries[i]
        idx_end = region_boundaries[i+1] - 1 # Inclusive end index
        
        region_id = region_map[idx_start]
        if region_id == 0:
            continue # Skip continuum gaps

        # Get velocity range
        v_min = x_kms[idx_start]
        v_max = x_kms[idx_end]
        
        if region_id not in region_info:
            region_info[region_id] = {
                'v_min': v_min, 
                'v_max': v_max, 
                'indices': [idx_start, idx_end] # Store as min/max index
            }
        else:
            # This handles regions split by a single pixel
            info = region_info[region_id]
            info['v_min'] = min(info['v_min'], v_min)
            info['v_max'] = max(info['v_max'], v_max)
            info['indices'][0] = min(info['indices'][0], idx_start)
            info['indices'][1] = max(info['indices'][1], idx_end)

    if not region_info:
        logging.debug("merge_regions: No regions found to merge.")
        return region_map, 0

    # 4. Sort regions by starting velocity
    sorted_regions = sorted(region_info.items(), key=lambda item: item[1]['v_min'])
    
    # 5. Iterate and merge
    merged_map = np.zeros_like(region_map, dtype=int)
    if not sorted_regions:
        return merged_map, 0
        
    new_id_counter = 1
    
    # Start with the first region
    current_id, current_info = sorted_regions[0]
    current_v_max = current_info['v_max']
    current_idx_min = current_info['indices'][0]
    current_idx_max = current_info['indices'][1]
    
    for i in range(1, len(sorted_regions)):
        next_id, next_info = sorted_regions[i]
        
        # Calculate velocity gap
        v_gap = next_info['v_min'] - current_v_max
        
        if v_gap < merge_dv:
            # --- Merge ---
            # Extend the v_max and index_max of the current merged region
            current_v_max = max(current_v_max, next_info['v_max'])
            current_idx_max = max(current_idx_max, next_info['indices'][1])
        else:
            # --- Finalize previous region ---
            # Fill the merged_map with the new ID
            merged_map[current_idx_min:current_idx_max+1] = new_id_counter
            new_id_counter += 1
            
            # Start a new region
            current_v_max = next_info['v_max']
            current_idx_min = next_info['indices'][0]
            current_idx_max = next_info['indices'][1]

    # Finalize the last region
    merged_map[current_idx_min:current_idx_max+1] = new_id_counter
    
    # 6. Re-run 'detect_regions' to clean up any pixels in the map
    # that were part of the *original* mask but *outside* the new
    # min/max indices (e.g., gaps we just filled).
    final_map, final_num_regions = detect_regions(merged_map > 0, min_pix=1)

    logging.info(f"Merged {len(sorted_regions)} regions into {final_num_regions} (dv={merge_dv} km/s).")
    
    return final_map, final_num_regions

def compute_identification_signal(
    spectrum: 'SpectrumV2',
    series_name: str,
    z_grid: np.ndarray,
    modul: float
) -> np.ndarray:
    """
    Computes the "squashed significance product" signal for a given series.
    Replicates the core logic of V1's _abs_like.
    
    Args:
        spectrum (SpectrumV2): The V2 spectrum object.
        series_name (str): The name of the series (e.g., 'CIV', 'Ly_a').
        z_grid (np.ndarray): The common, regular redshift grid to interpolate onto.
        modul (float): The V1 modulation parameter.
        
    Returns:
        np.ndarray: The computed signal on the z_grid.
    """
    try:
        from ..v1.functions import trans_parse
    except ImportError:
        logging.error("V1 'trans_parse' not available. Cannot compute signal.")
        return np.full_like(z_grid, np.nan)

    if series_name not in STANDARD_MULTIPLETS:
        logging.error(f"Series '{series_name}' not in STANDARD_MULTIPLETS.")
        return np.full_like(z_grid, np.nan)
    transitions = STANDARD_MULTIPLETS[series_name]
    if not transitions:
        logging.warning(f"No transitions found for series '{series_name}'.")
        return np.full_like(z_grid, np.nan)

    # 1. Get spectral data
    x_nm = spectrum.x.to_value(au.nm)
    y = spectrum.y.value
    dy = spectrum.dy.value
    
    # Use continuum if present, otherwise assume 1.0
    cont = spectrum.cont.value if spectrum.cont is not None else np.ones_like(y)
    
    # 2. Calculate "squashed significance"
    safe_dy = dy.copy()
    safe_dy[safe_dy <= 0] = np.nan
    with np.errstate(invalid='ignore', divide='ignore'):
        significance = (cont - y) / safe_dy
    
    # V1 _abs_like uses erf(np.abs(...))
    squashed_sig = erf(np.abs(significance) / np.sqrt(2) / modul)
    squashed_sig = np.nan_to_num(squashed_sig, nan=0.0) # V1 logic treats NaNs as 0

    # 3. Interpolate each transition onto the common z_grid
    abs_int = []
    for t in transitions:
        if t not in xem_d:
            logging.warning(f"Transition '{t}' not in xem_d. Skipping.")
            continue
        
        # Convert native x-axis (nm) to z for this transition
        xem_nm = xem_d[t].to_value(au.nm)
        with np.errstate(invalid='ignore', divide='ignore'):
            z_native = (x_nm / xem_nm) - 1.0
        
        # Interpolate the squashed significance onto the z_grid
        # (Handle unsorted z_native, which can happen with e.g. echelle)
        sort_idx = np.argsort(z_native)
        interp_sig = np.interp(
            z_grid, 
            z_native[sort_idx], 
            squashed_sig[sort_idx],
            left=np.nan, right=np.nan # Use NaN for areas outside spectral range
        )
        abs_int.append(interp_sig)

    # 4. Compute the signal product
    if not abs_int:
        logging.error(f"No valid transitions found for '{series_name}'.")
        return np.full_like(z_grid, np.nan)
        
    with np.errstate(invalid='ignore'):
        signal_prod = np.nanprod(abs_int, axis=0)
        
    signal = erfinv(signal_prod) * np.sqrt(2) * modul
    signal[np.isinf(signal)] = -np.inf # Handle erfinv(0)

    return signal

def find_signal_peaks(
    signals_dict: Dict[str, np.ndarray],
    z_grid: np.ndarray,
    sigma: float,
    distance_pix: int,
    prominence: Optional[float]
) -> List[Tuple[str, float, float]]:
    """
    Finds and filters peaks in the identification signals.
    Replicates the core logic of V1's _systs_like.
    
    Args:
        signals_dict: Dict of {series_name: signal_array}
        z_grid: The common redshift grid.
        sigma: The minimum peak height (signal score).
        distance_pix: Minimum distance between peaks in *pixels* (indices).
        prominence: Minimum prominence of peaks (absolute signal value).
        
    Returns:
        List of (series_name, z_peak, score)
    """
    all_peaks_list = []

    smoothed_signals = {} # Store smoothed signals
    
    # 1. Smooth each signal individually
    for series_name, signal in signals_dict.items():
        if signal is None:
            continue
            
        # --- *** FIXED: Clean NaNs before smoothing *** ---
        # The input signal contains NaNs from interpolation edges.
        # Replace NaNs with 0.0 (neutral) before smoothing.
        signal_clean = np.nan_to_num(signal, nan=0.0)
        
        # V1 logic used a light smoothing
        smoothed_signal = savgol_filter(signal_clean, window_length=5, polyorder=1)
        smoothed_signals[series_name] = smoothed_signal
        # --- *** END FIX *** ---

    # 2. Create the "total signal" envelope from the *SMOOTHED* signals
    if not smoothed_signals:
        logging.warning("find_signal_peaks: No valid signals to process.")
        return []
        
    with np.errstate(invalid='ignore'): # Ignore warnings from all-NaN slices
        total_signal = np.nanmax(
            [s for s in smoothed_signals.values()], 
            axis=0
        )
    # Replace any remaining NaNs (if all signals were NaN) with -inf
    total_signal = np.nan_to_num(total_signal, nan=-np.inf)

    # 3. Find peaks for each *SMOOTHED* series
    for series_name, smoothed_signal in smoothed_signals.items():

        logging.debug(f"find_signal_peaks: Searching '{series_name}'...")
        logging.debug(f"  Signal stats (min/mean/max): {np.min(smoothed_signal):.2f} / {np.mean(smoothed_signal):.2f} / {np.max(smoothed_signal):.2f}")
        logging.debug(f"  Thresholds: sigma={sigma}, distance_pix={distance_pix}, prominence={prominence}")
        
        # 4. Find peaks
        peaks_idx, _ = find_peaks(
            smoothed_signal, # <<< Use the smoothed signal
            distance=distance_pix, 
            prominence=prominence
        )
        
        # --- *** ADDED DEBUGGING *** ---
        logging.debug(f" Found {len(peaks_idx)} raw peaks (before sigma/best-guess filter).")
        # --- *** END DEBUGGING *** ---
        
        if not peaks_idx.size:
            continue
            
        # 5. Filter by sigma (height threshold)
        peaks_idx = peaks_idx[smoothed_signal[peaks_idx] > sigma]
        
        # --- *** ADDED DEBUGGING *** ---
        logging.debug(f"  Found {len(peaks_idx)} peaks after 'sigma' filter.")
        # --- *** END DEBUGGING *** ---

        # 6. Best-Guess Filter
        for p_idx in peaks_idx:
            peak_score = smoothed_signal[p_idx]
            
            # Check if this peak is the highest signal at this redshift
            # --- *** FIXED: Comparing smoothed_score to smoothed_total_signal *** ---
            if np.isclose(peak_score, total_signal[p_idx]):
                z_peak = z_grid[p_idx]
                all_peaks_list.append((series_name, z_peak, peak_score))

    # 7. Sort final list by redshift
    all_peaks_list.sort(key=lambda x: x[1])
    
    return all_peaks_list