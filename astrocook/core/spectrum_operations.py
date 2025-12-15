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
from scipy.stats import binned_statistic
from typing import Dict, List, Optional, Tuple

from astrocook.core.atomic_data import ATOM_DATA, get_multiplet_velocity_lags, STANDARD_MULTIPLETS, xem_d
from astrocook.legacy.functions import x_convert as x_convert_v1 
from astrocook.core.structures import DataColumnV2, SpectrumDataV2 
from astrocook.gui.debug_utils import GLOBAL_PLOTTER

# We assume 'xem' (the Ly-alpha rest wavelength) is accessible here for the custom equivalency
# You must ensure this global is either imported or defined here (Technical Debt).
XEM_LYA = 121.567 * au.nm 
XEM_LYA_KM_S_ZERO = (XEM_LYA * (1 + 0.0)).to(au.km/au.s, equivalencies=au.equivalencies.doppler_optical(XEM_LYA)).value

def convert_axis_velocity(
    x_quantity: au.Quantity, 
    z_ref: float, # Renamed from zem
    target_unit: au.Unit
) -> au.Quantity:
    """
    Converts X-axis data between wavelength (nm) and velocity (km/s) space,
    using the V1 zero-point logic (Ly-alpha rest-frame relative to z_ref).
    """
    
    v_c = c_light.to(au.km/au.s).value
    
    # Define a custom conversion function
    # Lambda to Velocity (Forward)
    def lambda_to_v(lambda_val):
        # Use z_ref as the zero-point for the velocity conversion
        return np.log(lambda_val / (XEM_LYA.value * (1 + z_ref))) * v_c

    # Velocity to Lambda (Reverse)
    def v_to_lambda(v_val):
        return np.exp(v_val / v_c) * (XEM_LYA.value * (1 + z_ref))

    custom_equiv = [
        (au.nm, au.km/au.s, lambda_to_v, v_to_lambda),
        (au.Angstrom, au.km/au.s, lambda_to_v, v_to_lambda) # Add Angstrom
    ]
    
    if target_unit.is_equivalent(au.km/au.s) or target_unit.is_equivalent(au.nm):
        return x_quantity.to(target_unit, equivalencies=custom_equiv)
    else:
        return x_quantity.to(target_unit)

def convert_x_axis(spec_data: SpectrumDataV2, z_ref: float, xunit: au.Unit) -> Tuple[au.Quantity, au.Quantity, au.Quantity]:
    """
    Performs X-axis conversion based on V1 logic and returns new Quantity arrays.
    """
    
    x = spec_data.x.quantity
    xmin = spec_data.xmin.quantity
    xmax = spec_data.xmax.quantity
    
    if xunit.is_equivalent(au.km/au.s) or x.unit.is_equivalent(au.km/au.s):
        # Use the V2 dedicated velocity conversion logic
        x_new = convert_axis_velocity(x, z_ref, xunit)
        xmin_new = convert_axis_velocity(xmin, z_ref, xunit)
        xmax_new = convert_axis_velocity(xmax, z_ref, xunit)
        
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
    z_ref: float = 0.0
) -> np.ndarray:
    """
    Smooths a flux array using a Gaussian filter with a given std in km/s.
    """
    
    # 1. Smooth the result (V1's _gauss_convolve logic)
    if sigma_kms > 0:
        # 2. Convert x-axis to km/s to get pixel size
        try:
            # We need the x-axis as a Quantity to convert it
            x_kms = convert_axis_velocity(x, z_ref=z_ref, target_unit=au.km/au.s).value
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

def compute_flux_scaling(
    x_ref: np.ndarray, 
    y_ref: np.ndarray, 
    x_in: np.ndarray, 
    y_in: np.ndarray, 
    order: int = 0,
    smooth_window_kms: float = 1000.0  # The "Heavy" Smoothing Kernel
) -> np.ndarray:
    """
    Computes a flux scaling factor by comparing heavily smoothed versions of the spectra.
    
    Strategy:
    1. Resample both spectra onto a coarse velocity grid (~50 km/s).
    2. Apply a large Gaussian smooth (smooth_window_kms) to wash out local features.
    3. Compute the ratio of the smoothed profiles.
    4. Fit a polynomial (if order >= 0) or return the spline (if order < 0).
    """
    correction = np.ones_like(y_in)

    # 1. Identify Global Overlap
    x_min_ov = max(np.nanmin(x_ref), np.nanmin(x_in))
    x_max_ov = min(np.nanmax(x_ref), np.nanmax(x_in))
    
    if x_max_ov <= x_min_ov:
        return correction

    # 2. Define a Coarse Calculation Grid (Linear in Velocity / Log Lambda)
    # We use ~50 km/s pixels. This is high enough to trace the continuum 
    # but low enough to make smoothing instant.
    # c_kms = 299792.458
    # d_log_lam = 50.0 / c_kms
    # num_bins = int(np.log(x_max_ov / x_min_ov) / d_log_lam)
    
    # Simple Linear Approximation (sufficient for scaling)
    # Average x ~ 500nm -> 50km/s is ~0.08 nm.
    step_nm = (x_min_ov + x_max_ov) / 2.0 * (50.0 / 300000.0)
    num_bins = int((x_max_ov - x_min_ov) / step_nm)
    num_bins = max(num_bins, 10) # Safety floor

    bin_edges = np.linspace(x_min_ov, x_max_ov, num_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # 3. Bin the Spectra (Handling "Shredded" Pixels naturally)
    # We ignore NaNs and negative fluxes to get a clean continuum estimate.
    mask_ref = (x_ref >= x_min_ov) & (x_ref <= x_max_ov) & (y_ref > 0) & np.isfinite(y_ref)
    mask_in = (x_in >= x_min_ov) & (x_in <= x_max_ov) & (y_in > 0) & np.isfinite(y_in)
    
    if np.sum(mask_ref) < 10 or np.sum(mask_in) < 10:
        return correction

    # binned_statistic handles unsorted/overlapped arrays automatically
    ref_binned, _, _ = binned_statistic(x_ref[mask_ref], y_ref[mask_ref], statistic='median', bins=bin_edges)
    in_binned, _, _ = binned_statistic(x_in[mask_in], y_in[mask_in], statistic='median', bins=bin_edges)

    # 4. Handle Gaps (Interpolate over empty bins)
    def fill_nans(arr):
        mask = np.isfinite(arr)
        if not np.any(mask): return arr
        return np.interp(np.arange(len(arr)), np.arange(len(arr))[mask], arr[mask])

    ref_filled = fill_nans(ref_binned)
    in_filled = fill_nans(in_binned)

    # 5. Heavy Gaussian Smoothing
    # Convert window (km/s) to bin pixels.
    # We used 50 km/s bins, so sigma ~ window / 50.
    sigma_pix = smooth_window_kms / 50.0
    
    ref_smooth = gaussian_filter1d(ref_filled, sigma_pix)
    in_smooth = gaussian_filter1d(in_filled, sigma_pix)

    # 6. Compute Ratio of Smoothed Profiles
    # Avoid div by zero
    valid = (in_smooth > 0)
    if np.sum(valid) < 5:
        return correction

    ratio_grid = np.ones_like(ref_smooth)
    ratio_grid[valid] = ref_smooth[valid] / in_smooth[valid]
    
    # 7. Generate Correction Model
    if order < 0:
        # Spline Interpolation (Follows the smooth ratio exactly)
        interp = interp1d(bin_centers, ratio_grid, kind='linear', fill_value='extrapolate')
        return interp(x_in)
        
    elif order == 0:
        # Robust Scalar (Median of the smooth ratio)
        scalar = np.nanmedian(ratio_grid)
        return np.full_like(y_in, scalar)
        
    else:
        # Polynomial Fit to the smooth ratio
        coeffs = np.polyfit(bin_centers[valid], ratio_grid[valid], order)
        from numpy.polynomial.polynomial import polyval
        # Note: np.polyfit returns coeffs highest-to-lowest. 
        # numpy.polynomial.polynomial.polyval expects lowest-to-highest.
        # We use np.polyval (standard) for consistency with polyfit.
        return np.polyval(coeffs, x_in)
    

def rebin_spectrum(
        x_in: au.Quantity, xmin_in: au.Quantity, xmax_in: au.Quantity,
        spec_data: SpectrumDataV2,
        xstart: Optional[au.Quantity], xend: Optional[au.Quantity], 
        dx: au.Quantity, calc_unit: au.Unit,
        kappa: Optional[float], filling: float
    ) -> Tuple[au.Quantity, au.Quantity, au.Quantity, au.Quantity, au.Quantity]:
    """
    High-performance vectorized rebinning (Inverse/Drizzle method)
    with Iterative Local Sigma Clipping.
    """
    # 1. Prepare Units and Arrays
    xunit_target = dx.unit
    z_ref = spec_data.z_em if spec_data.z_em else 0.0

    def _to_target_val(q: au.Quantity):
        if q is None: return None
        return convert_axis_velocity(q, z_ref, xunit_target).value

    # Output Grid
    if xstart is None: xstart_val = np.nanmin(_to_target_val(x_in))
    else: xstart_val = _to_target_val(xstart)

    if xend is None: xend_val = np.nanmax(_to_target_val(x_in))
    else: xend_val = _to_target_val(xend)
        
    dx_val = dx.to_value(xunit_target)
    edges_out = np.arange(xstart_val, xend_val + dx_val, dx_val)
    x_out_val = 0.5 * (edges_out[:-1] + edges_out[1:])
    N_out = len(x_out_val)
    
    # Input Arrays
    xmin_in_val = _to_target_val(xmin_in)
    xmax_in_val = _to_target_val(xmax_in)
    x_in_val = 0.5 * (xmin_in_val + xmax_in_val) # Center for interpolation
    
    y_val = spec_data.y.quantity.to_value(spec_data.y.unit)
    dy_val = spec_data.dy.quantity.to_value(spec_data.y.unit)
    
    # Initial Mask (NaNs and bad pixels)
    valid_mask = np.isfinite(y_val) & (xmax_in_val > xmin_in_val) & (dy_val > 0)
    
    # --- HELPER: Core Drizzle Logic ---
    def _drizzle_pass(mask_in):
        # Filter inputs
        xmin_v = xmin_in_val[mask_in]
        xmax_v = xmax_in_val[mask_in]
        y_v = y_val[mask_in]
        dy_v = dy_val[mask_in]
        
        # Map Inputs to Outputs
        idx_start = np.searchsorted(edges_out, xmin_v, side='right') - 1
        idx_end = np.searchsorted(edges_out, xmax_v, side='left') - 1
        
        num = np.zeros(N_out)
        den = np.zeros(N_out)
        
        max_span = np.max(idx_end - idx_start) if len(idx_start) > 0 else 0
        
        for k in range(int(max_span) + 1):
            target_bin_indices = idx_start + k
            active = (target_bin_indices <= idx_end) & \
                     (target_bin_indices >= 0) & \
                     (target_bin_indices < N_out)
            
            if not np.any(active): continue
            
            valid_sub = np.where(active)[0]
            bin_idx = target_bin_indices[valid_sub]
            
            bin_edge_left = edges_out[bin_idx]
            bin_edge_right = edges_out[bin_idx+1]
            
            ov_min = np.maximum(xmin_v[valid_sub], bin_edge_left)
            ov_max = np.minimum(xmax_v[valid_sub], bin_edge_right)
            
            overlap_width = np.maximum(0.0, ov_max - ov_min)
            frac = overlap_width / dx_val
            
            sigma_sq = dy_v[valid_sub]**2
            weights = frac / sigma_sq
            
            np.add.at(num, bin_idx, y_v[valid_sub] * weights)
            np.add.at(den, bin_idx, weights)
            
        return num, den

    # --- PASS 1: Consensus ---
    num1, den1 = _drizzle_pass(valid_mask)
    
    # --- ITERATIVE LOCAL CLIPPING ---
    if kappa is not None:
        # 1. Compute Model (Pass 1 Result)
        valid_bins = den1 > 0
        y_model_grid = np.zeros(N_out)
        y_model_grid[valid_bins] = num1[valid_bins] / den1[valid_bins]
        
        # 2. Interpolate Model to Input Coordinates
        # Use nearest neighbor or linear to map output bin value to input pixel
        # (Linear is safer for sub-pixel accuracy)
        model_interp = interp1d(x_out_val, y_model_grid, kind='linear', 
                                bounds_error=False, fill_value=0.0)
        
        y_predicted = model_interp(x_in_val)
        
        # 3. Compute Residuals (Local deviation)
        residuals = np.abs(y_val - y_predicted)
        
        # 4. Clip
        # We assume 'dy_val' is the sigma. 
        # Points > kappa * sigma away from the consensus are rejected.
        is_outlier = residuals > (kappa * dy_val)
        
        # Update Mask (Keep original valid points that are NOT outliers)
        valid_mask &= ~is_outlier

    # --- PASS 2: Final Rebinning ---
    num_final, den_final = _drizzle_pass(valid_mask)

    # Finalize Output
    valid_out = den_final > 0
    y_out_val = np.full(N_out, filling)
    dy_out_val = np.full(N_out, filling)
    
    y_out_val[valid_out] = num_final[valid_out] / den_final[valid_out]
    dy_out_val[valid_out] = 1.0 / np.sqrt(den_final[valid_out])
    
    # Reconstruct Quantities
    x_out = au.Quantity(x_out_val, xunit_target)
    xmin_out = au.Quantity(edges_out[:-1], xunit_target)
    xmax_out = au.Quantity(edges_out[1:], xunit_target)
    y_out = au.Quantity(y_out_val, spec_data.y.unit)
    dy_out = au.Quantity(dy_out_val, spec_data.y.unit)
    
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
    Finds absorbed regions using iterative clipping with robust V1 parity logic.
    Fixes:
      1. Window size calculated from min(dv) instead of mean(dv) to preserve narrow features.
      2. Interpolation performed in wavelength space, not index space.
    """
    
    # --- 1. Initial Masking ---
    mask_unabs = np.ones_like(y, dtype=bool)
    valid_data = ~np.isnan(y) & ~np.isnan(dy) & (dy > 0)
    mask_unabs[~valid_data] = False
    
    # Exclude low S/N regions from being 'continuum anchors'
    with np.errstate(invalid='ignore'):
        mask_unabs[y <= 3.0 * dy] = False 
    
    # Define Ly-alpha proximity
    if z_em == 0.0:
        logging.warning("z_em is 0.0. Using 'smooth_len_out' for all regions.")
        lya_prox_nm = -np.inf * au.nm
    else:
        # We need the global XEM_LYA constant here. 
        # Ensure it is imported or defined as 121.567 * au.nm
        lya_obs = (121.567 * au.nm) * (1 + z_em) 
        prox_kms = 5000 * au.km/au.s
        c_kms = 299792.458 * au.km/au.s
        lya_prox_nm = lya_obs * (1 - prox_kms / c_kms)

    # --- 2. Handle Template (Placeholder) ---
    y_proc = y.copy()
    if template:
        pass # Template logic would go here

    # --- 3. Calculate pixel size (dv) & Setup X ---
    # We need dv in km/s for window size calculation
    x_kms = convert_axis_velocity(x, z_ref=0.0, target_unit=au.km/au.s).value
    dv = np.gradient(x_kms)
    x_val = x.value # Work with raw values for speed

    # --- 4. Define Regions ---
    # Convert mask to arrays for boolean indexing
    in_lya_forest = (x < lya_prox_nm)
    out_of_forest = (x >= lya_prox_nm)
    
    intervals = [
        (in_lya_forest, smooth_len_lya),
        (out_of_forest, smooth_len_out)
    ]
    
    mask_clipped = np.zeros_like(y, dtype=bool) 

    # --- 5. Iterative Clipping Loop ---
    for i in range(maxiter):
        
        y_rm_full = np.full_like(y, np.nan)
        
        for region_bool, smooth_len in intervals:
            # Combine geometric region with the current unabsorbed mask
            mask_interval = region_bool & mask_unabs
            
            # Points belonging to this interval AND not yet clipped locally
            current_mask = mask_interval & ~mask_clipped
            
            if np.sum(current_mask) < 10: 
                continue
                
            # V1 FIX: Use np.min(dv) for a conservative (smaller pixels -> larger window)
            # Use abs() because dv can be negative if spectrum is flipped
            dv_segment = np.abs(dv[current_mask])
            min_dv = np.min(dv_segment) if len(dv_segment) > 0 else 1.0
            if min_dv == 0: min_dv = 1.0 
            
            # V1 Logic: window = length / min_dv / 8
            # (The '/8' is a V1 heuristic for 'local' smoothing)
            smooth_val = smooth_len.to_value(au.km/au.s)
            hwindow = int(smooth_val / min_dv / 8)
            hwindow = max(3, hwindow) # Safety floor
            
            y_interval = y_proc[current_mask]
            
            if len(y_interval) < hwindow * 2 + 1:
                y_rm_interval = np.full_like(y_interval, np.mean(y_interval))
            else:
                y_rm_interval = running_mean(y_interval, h=hwindow)
            
            # V1 FIX: Interpolate using X coordinates, not indices
            # We interpolate the running mean (defined on 'current_mask') 
            # onto the full unabsorbed set ('mask_interval')
            x_target = x_val[mask_interval]      
            x_source = x_val[current_mask]       
            
            # np.interp requires sorted x. 
            if x_source[0] < x_source[-1]:
                y_rm_interp = np.interp(x_target, x_source, y_rm_interval)
            else:
                # Handle descending x-axis
                sort_idx = np.argsort(x_source)
                y_rm_interp = np.interp(x_target, x_source[sort_idx], y_rm_interval[sort_idx])
            
            y_rm_full[mask_interval] = y_rm_interp

        # 6. Find new outliers
        # Compare raw data against the smoothed model
        # V1 Logic: outlier if y < model - kappa * sigma
        # (We rely on NaNs in y_rm_full to skip invalid regions)
        with np.errstate(invalid='ignore'):
            is_outlier = (y_proc - y_rm_full) < (-kappa * dy)
        
        # Only mark as outlier if it was previously considered valid/unabsorbed
        is_outlier[~mask_unabs] = False 
        
        # Check convergence
        if np.sum(is_outlier) == 0:
            logging.info(f"Clipping converged after {i+1} iterations.")
            break
        
        mask_clipped = mask_clipped | is_outlier 
        mask_unabs = mask_unabs & ~is_outlier
        
        if i == maxiter - 1:
            logging.warning(f"Clipping did not converge after {maxiter} iterations.")

    return ~mask_unabs

def compute_auto_fudge(
    y: np.ndarray,
    mask_abs: np.ndarray,
    cont_model: np.ndarray
) -> float:
    """
    Calculates a scalar fudge factor to center the continuum residuals.
    Fixes the 'Low Bias' issue where sigma-clipped continuums sit too low.
    """
    mask_unabs = ~mask_abs
    if np.sum(mask_unabs) < 10: 
        return 1.0
    
    y_valid = y[mask_unabs]
    model_valid = cont_model[mask_unabs]
    
    # Calculate Residuals (Data - Model)
    # If model is low, residuals are positive.
    with np.errstate(invalid='ignore'):
        residuals = y_valid - model_valid
    
    # Use Median to be robust against non-Gaussian outliers
    avg_residual = np.nanmedian(residuals)
    avg_model = np.nanmedian(model_valid)
    
    if avg_model == 0 or np.isnan(avg_model): 
        return 1.0
    
    # V1 Logic: fudge = 1 + (residual / model)
    fudge_factor = 1.0 + (avg_residual / avg_model)
    
    return float(fudge_factor)

def interpolate_continuum_mask(
    x: np.ndarray,
    y: np.ndarray,
    mask_abs: np.ndarray
) -> np.ndarray:
    """
    Interpolates the unabsorbed points to create a raw continuum guess.
    """
    mask_unabs = ~mask_abs
    x_unabs = x[mask_unabs]
    y_unabs = y[mask_unabs]
    
    if len(x_unabs) < 2:
        raise ValueError("Cannot fit continuum: < 2 unabsorbed points found.")
        
    # Simple linear interpolation over the gaps
    cont_interp = np.interp(x, x_unabs, y_unabs)
    return cont_interp

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
    z_ref: float,
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
        x_kms = convert_axis_velocity(x, z_ref=z_ref, target_unit=au.km/au.s).value
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
        from astrocook.legacy.functions import trans_parse
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

def _find_kinematic_doublet_candidates(
    region_map: np.ndarray,
    x: au.Quantity,
    z_ref: float,
    multiplet_name: str,
    z_em: float
) -> List[Tuple[int, float]]:
    """
    Performs a kinematic region-overlap check to find doublet candidates.
    This is Step 2 & 3 of the new algorithm.
    """
    if not multiplet_name in STANDARD_MULTIPLETS:
        return []
    
    lines = STANDARD_MULTIPLETS[multiplet_name]
    if len(lines) < 2:
        return [] # Not a doublet

    try:
        # 1. Get velocity lags
        lags_kms = get_multiplet_velocity_lags(multiplet_name)
        primary_line, secondary_line = lines[0], lines[1]
        dv_kms = lags_kms[secondary_line] # e.g., ~498 km/s for CIV
        
        # 2. Get x-axis in velocity space
        x_kms = convert_axis_velocity(x, z_ref=z_ref, target_unit=au.km/au.s).value
        
        # 3. Create a "shifted" region map
        # Interpolate the region map onto an axis shifted by dv_kms
        # We use 'nearest' to preserve integer region IDs
        interp = interp1d(x_kms + dv_kms, region_map, kind='nearest', 
                          bounds_error=False, fill_value=0)
        shifted_map = interp(x_kms).astype(int)
        
        # 4. Find overlaps
        # Find where both original and shifted maps have a *non-zero* region ID
        overlap_mask = (region_map > 0) & (shifted_map > 0)
        if not np.any(overlap_mask):
            return []
            
        # 5. Get all (region_id, shifted_region_id) pairs
        overlapping_pairs = np.unique(
            np.stack((region_map[overlap_mask], shifted_map[overlap_mask])), 
            axis=1
        )
        
        # 6. For each unique pair, calculate the best-guess redshift
        candidates = []
        x_nm = x.to_value(au.nm)
        secondary_wave_nm = xem_d[lines[1]].to_value(au.nm)
        
        for (rid_1, rid_2) in overlapping_pairs.T:
            # Find the pixels corresponding to this specific overlap
            pair_mask = (region_map == rid_1) & (shifted_map == rid_2)
            
            # Calculate the geometric mean of the x-range of this overlap
            x_overlap_nm = x_nm[pair_mask]
            if len(x_overlap_nm) == 0: continue
            
            x_min = np.min(x_overlap_nm)
            x_max = np.max(x_overlap_nm)
            x_center = np.sqrt(x_min * x_max) # Geometric mean is good for log-space
            
            # Calculate the z_test for this overlap
            z_test = (x_center / secondary_wave_nm) - 1.0
            
            # Add (region_id_of_primary, z_test)
            candidates.append((rid_1, rid_2, z_test))

        return candidates
        
    except Exception as e:
        logging.error(f"Failed kinematic check for {multiplet_name}: {e}", exc_info=True)
        return []

def rate_doublet_candidate(
    spectrum: 'SpectrumV2',
    region_mask_1: np.ndarray,
    region_mask_2: np.ndarray,
    multiplet_name: str,
    z_test: float,
    debug_rating: bool = False
) -> float:
    """
    Rates a doublet candidate using the R^2 score on the 
    NORMALIZED FLUX profiles (Power-Law Model).
    
    Physics: tau_1 = gamma * tau_2  ==>  F_1 = F_2 ^ gamma
    
    NOTE: region_mask_1 is the BLUER region (e.g., CIV_1548)
          region_mask_2 is the REDDER region (e.g., CIV_1550)
    """
    try:
        # 1. Get lines and data
        lines = STANDARD_MULTIPLETS[multiplet_name]
        primary_line, secondary_line = lines[0], lines[1] 
        
        lags_kms = get_multiplet_velocity_lags(multiplet_name)
        dv_kms = lags_kms[secondary_line] 

        f_1 = ATOM_DATA[primary_line]['f']
        f_2 = ATOM_DATA[secondary_line]['f']
        lam_1 = ATOM_DATA[primary_line]['wave']
        lam_2 = ATOM_DATA[secondary_line]['wave']
        
        if f_2 == 0: return 0.0
        
        # Exponent for the flux relationship: tau_1/tau_2 = (f1*lam1)/(f2*lam2)
        gamma_exp = (f_1 * lam_1) / (f_2 * lam_2)

        # 2. Get Velocity Grid (relative to primary)
        primary_wave_obs = xem_d[primary_line] * (1.0 + z_test)
        c_kms = c_light.to(au.km/au.s)
        
        x_kms_test_frame = (spectrum.x - primary_wave_obs).to(
            au.km/au.s, equivalencies=au.equivalencies.doppler_optical(primary_wave_obs)
        ).value

        # 3. Get Global Normalized Flux
        cont = spectrum.cont.value if spectrum.cont is not None else np.ones_like(spectrum.y.value)
        
        # Handle bad continuum or zero division
        norm_flux = np.ones_like(spectrum.y.value)
        valid_cont = (cont > 0)
        norm_flux[valid_cont] = spectrum.y.value[valid_cont] / cont[valid_cont]
        
        # Clip unphysical values for stability (flux must be >= 0 for power law)
        norm_flux = np.clip(norm_flux, 0.0, 1.5)

        # Global Interpolator
        interp_flux = interp1d(x_kms_test_frame, norm_flux, kind='linear', bounds_error=False, fill_value=1.0)

        # 4. Define Intersection Window
        v_region_blue = x_kms_test_frame[region_mask_2]
        v_region_red = x_kms_test_frame[region_mask_1]
        if len(v_region_blue) < 2 or len(v_region_red) < 2: return 0.0

        v_red_shifted = v_region_red - dv_kms 

        v_start = max(np.min(v_region_blue), np.min(v_red_shifted))
        v_end = min(np.max(v_region_blue), np.max(v_red_shifted))
        
        padding_kms = 30.0
        v_start -= padding_kms
        v_end += padding_kms

        if v_start >= v_end: return 0.0

        # Select grid points
        mask_compare = (x_kms_test_frame >= v_start) & (x_kms_test_frame <= v_end)
        v_compare = x_kms_test_frame[mask_compare]
        if len(v_compare) < 3: return 0.0

        # 5. Construct Vectors (Absorption Signal: 1 - F)
        
        # Y_data: The Blue Flux (Observed)
        flux_blue_raw = interp_flux(v_compare)
        Y_data = 1.0 - flux_blue_raw

        # Y_model: The Red Flux transformed to look like Blue
        # We sample Red at (v + dv)
        flux_red_raw = interp_flux(v_compare + dv_kms)
        
        # Apply Physics: F_blue_model = (F_red)^(gamma)
        flux_blue_model = np.power(flux_red_raw, gamma_exp)
        
        # Convert to "Absorption Signal" for correlation
        Y_model = 1.0 - flux_blue_model

        # 6. Calculate R^2
        # Note: We do NOT mean-subtract here. We want to match the absolute zero-point (continuum).
        ss_res = np.sum((Y_data - Y_model)**2)
        ss_tot = np.sum((Y_data - np.mean(Y_data))**2)
        
        if ss_tot == 0: return 0.0
            
        r2_score = 1.0 - (ss_res / ss_tot)
        
        if debug_rating:
            plot_data = {
                'title': f"{multiplet_name} z={z_test:.4f} (R2={r2_score:.3f})",
                'v_compare': v_compare,
                'Y_data': Y_data,   # Blue Absorption (1-F)
                'Y_model': Y_model, # Modeled Blue Absorption
            }
            GLOBAL_PLOTTER.plot_requested.emit(plot_data)

        return max(0.0, r2_score)

    except Exception as e:
        logging.warning(f"Failed to rate {multiplet_name} at z={z_test}: {e}", exc_info=True)
        return 0.0
    
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
    
    *** MODIFIED to be a generic peak finder. ***
    *** It no longer smooths or runs a 'best-guess' filter. ***
    
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

    # 1. Find peaks for each series
    for series_name, signal in signals_dict.items():
        if signal is None:
            continue
        
        # --- *** MODIFIED: Clean NaNs but do NOT smooth *** ---
        signal_clean = np.nan_to_num(signal, nan=-np.inf)
        # --- *** END MODIFIED *** ---

        logging.debug(f"find_signal_peaks: Searching '{series_name}'...")
        logging.debug(f"  Signal stats (min/mean/max): {np.min(signal_clean):.2f} / {np.mean(signal_clean):.2f} / {np.max(signal_clean):.2f}")
        logging.debug(f"  Thresholds: sigma={sigma}, distance_pix={distance_pix}, prominence={prominence}")
        
        # 4. Find peaks
        peaks_idx, _ = find_peaks(
            signal_clean,
            distance=distance_pix, 
            prominence=prominence,
            height=sigma # Use 'height' for simple threshold
        )
        
        # --- *** ADDED DEBUGGING *** ---
        logging.debug(f" Found {len(peaks_idx)} raw peaks (after height/distance/prominence).")
        # --- *** END DEBUGGING *** ---
        
        if not peaks_idx.size:
            continue
            
        # 5. Filter by sigma (height threshold)
        # (This is now handled by the 'height' param in find_peaks)
        
        # 6. Best-Guess Filter (REMOVED)

        # 7. Add all found peaks
        for p_idx in peaks_idx:
            peak_score = signal_clean[p_idx]
            z_peak = z_grid[p_idx]
            all_peaks_list.append((series_name, z_peak, peak_score))

    # 8. Sort final list by redshift
    all_peaks_list.sort(key=lambda x: x[1])
    
    return all_peaks_list

