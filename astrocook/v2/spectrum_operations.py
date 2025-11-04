import bisect
from astropy import units as au
from astropy.constants import c as c_light
from astropy.stats import sigma_clip
import logging
import numpy as np
from typing import Optional, Tuple

from ..v1.functions import x_convert as x_convert_v1 
from ..v2.structures import DataColumnV2, SpectrumDataV2 

# We assume 'xem' (the Ly-alpha rest wavelength) is accessible here for the custom equivalency
# You must ensure this global is either imported or defined here (Technical Debt).
XEM_LYA = 121.567 * au.nm 

def convert_axis_velocity(
    x_quantity: au.Quantity, 
    zem: float, 
    target_unit: au.Unit
) -> au.Quantity:
    """
    Converts X-axis data between wavelength (nm) and velocity (km/s) space,
    using the V1 zero-point logic (Ly-alpha rest-frame relative to zem).
    
    This function handles both forward (Wavelength -> Velocity) and 
    reverse (Velocity -> Wavelength) conversion paths.
    """
    
    # 1. Define the V1 Custom Equivalency (using Ly-alpha and light speed)
    # The zero-point (xem) used for the conversion must be defined or accessed.
    # We use a standard value for the V1 utility dependency for now.
    
    # NOTE: The custom equivalency depends on a zero-point (xem). The V1 code 
    # often used Lya, or calculated xem from the input z. Since this conversion 
    # relies on a single fixed zero-point for the velocity scale, we'll use Lya for the example.
    
    # We must replicate the *effect* of your original custom equiv:
    v_c = c_light.to(au.km/au.s).value
    
    # Define a custom conversion function using Astropy's equivalence infrastructure
    # Since you used a fixed zero point for your conversion:
    
    # Lambda to Velocity (Forward)
    def lambda_to_v(lambda_val):
        return np.log(lambda_val / (XEM_LYA.value * (1 + zem))) * v_c

    # Velocity to Lambda (Reverse)
    def v_to_lambda(v_val):
        return np.exp(v_val / v_c) * (XEM_LYA.value * (1 + zem))

    custom_equiv = [
        (au.nm, au.km/au.s, lambda_to_v, v_to_lambda)
    ]
    
    # 2. Perform the conversion
    if target_unit.is_equivalent(au.km/au.s) or target_unit.is_equivalent(au.nm):
        return x_quantity.to(target_unit, equivalencies=custom_equiv)
    else:
        # For simple wavelength-to-wavelength conversion (e.g., Angstrom to nm)
        return x_quantity.to(target_unit)

def convert_x_axis(spec_data: SpectrumDataV2, zem: float, xunit: au.Unit) -> Tuple[au.Quantity, au.Quantity, au.Quantity]:
    """
    Performs X-axis conversion based on V1 logic and returns new Quantity arrays.
    """
    
    x = spec_data.x.quantity
    xmin = spec_data.xmin.quantity
    xmax = spec_data.xmax.quantity
    
    # Determine if conversion requires complex velocity math (V1 logic)
    if xunit.is_equivalent(au.km/au.s) or x.unit.is_equivalent(au.km/au.s):
        # Use the V2 dedicated velocity conversion logic (forward or reverse)
        x_new = convert_axis_velocity(x, zem, xunit)
        xmin_new = convert_axis_velocity(xmin, zem, xunit)
        xmax_new = convert_axis_velocity(xmax, zem, xunit)
        
    else:
        # Simple wavelength-to-wavelength or velocity-to-velocity conversion (standard Astropy)
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
    This is a standard numpy implementation of the 'h' parameter logic.
    """
    # Use 'same' mode to keep the array length, handling edges
    window_size = 2 * h + 1
    # Handle edge case where window is larger than array
    if window_size > len(y):
        logging.warning("Running mean window is larger than array, returning array mean.")
        return np.full_like(y, np.mean(y))
        
    return np.convolve(y, np.ones(window_size) / window_size, mode='same')

def find_unabsorbed_regions(
    y: np.ndarray, 
    dy: np.ndarray, 
    hwindow: int, 
    kappa: float, 
    maxiter: int = 100
) -> np.ndarray:
    """
    Pure function to find unabsorbed regions via kappa-sigma clipping.
    Based on the core logic of V1's 'clip_flux'.
    
    Returns:
        mask_unabs (np.ndarray): A boolean mask where True means "unabsorbed".
    """
    # Start with a mask where everything is unabsorbed
    mask_unabs = np.ones_like(y, dtype=bool)
    
    # Exclude regions with very low S/N, as they are unreliable
    # Also exclude NaN values from the start
    valid_data = ~np.isnan(y) & ~np.isnan(dy) & (dy > 0)
    mask_unabs[~valid_data] = False
    mask_unabs[y <= 3.0 * dy] = False # Low S/N
    
    for i in range(maxiter):
        # 1. Get the data from the *current* unabsorbed regions
        y_unabs = y[mask_unabs]
        
        if len(y_unabs) < hwindow * 2 + 1:
            logging.warning("Clipping stopped: not enough unabsorbed points to continue.")
            break # Not enough points to calculate a running mean
            
        # 2. Compute the running mean *only* on those points
        y_rm_unabs = running_mean(y_unabs, h=hwindow)
        
        # 3. Interpolate this mean back to the *full* array's size
        indices_full = np.arange(len(y))
        indices_unabs = indices_full[mask_unabs]
        y_rm_full = np.interp(indices_full, indices_unabs, y_rm_unabs)

        # 4. Find new outliers (absorbed regions)
        # An outlier is a point *below* the mean by kappa*sigma
        is_outlier = (y - y_rm_full) < (-kappa * dy)
        
        # 5. Update the mask
        new_mask_unabs = mask_unabs & ~is_outlier
        
        if np.all(new_mask_unabs == mask_unabs):
            logging.info(f"Clipping converged after {i+1} iterations.")
            break # No new outliers found, we are done
            
        mask_unabs = new_mask_unabs
        
        if i == maxiter - 1:
            logging.warning(f"Clipping did not converge after {maxiter} iterations.")
            
    return mask_unabs