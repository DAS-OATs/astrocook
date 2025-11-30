import numpy as np
import pytest
from astropy import units as au
from astrocook.v2.structures import SpectrumDataV2, DataColumnV2
from astrocook.v2.spectrum_operations import rebin_spectrum

# --- GOLDEN STANDARD CASE 1: Weighted Average and Error Propagation ---

# Input Data: 3 pixels, perfect overlap (simplified weighted average test)
def get_mock_input_case1():
    x_unit = au.nm
    y_unit = au.erg / (au.cm**2 * au.s * au.nm)
    
    # Input Pixels P1, P2, P3
    x_in_val = np.array([100.00, 100.05, 100.10])
    y_in_val = np.array([1.0, 0.5, 2.0])
    dy_in_val = np.array([0.10, 0.05, 0.20])
    
    # Bin limits (dx = 0.05 nm for all)
    xmin_in_val = x_in_val - 0.025
    xmax_in_val = x_in_val + 0.025
    
    spec_data = SpectrumDataV2(
        x=DataColumnV2(x_in_val, x_unit),
        xmin=DataColumnV2(xmin_in_val, x_unit),
        xmax=DataColumnV2(xmax_in_val, x_unit),
        y=DataColumnV2(y_in_val, y_unit),
        dy=DataColumnV2(dy_in_val, y_unit),
    )
    
    # Rebinning Parameters
    params = {
        'xstart': au.Quantity(100.025, au.nm), # Start just before P2
        'xend': au.Quantity(100.125, au.nm),   # End just after P3
        'dx': au.Quantity(0.05, au.nm),        # New bin size (same as input for simple overlap)
        'kappa': 100.0,                        # No clipping
        'filling': np.nan,
    }
    
    # Expected Output (Calculated Golden Standard for P2 + P3 overlap)
    # Output bin at 100.075 nm (P2 + P3 overlap)
    raw_expected_y = [0.5882352941176471]
    raw_expected_dy = [0.04850731464335431]
    
    expected_y = au.Quantity(raw_expected_y, y_unit)
    expected_dy = au.Quantity(raw_expected_dy, y_unit)    

    return spec_data, params, expected_y, expected_dy

# --- Test Function ---

def test_rebin_weighted_average(tolerance=1e-9):
    spec_data, params, expected_y, expected_dy = get_mock_input_case1()
    
    # 1. Extract X arrays from spec_data for positional use
    x_in = spec_data.x.quantity
    xmin_in = spec_data.xmin.quantity
    xmax_in = spec_data.xmax.quantity
    
    # 2. Define the calculation unit (which is the unit of the input X arrays here)
    calc_unit = x_in.unit 
    
    # Call the pure function with all required positional arguments
    x_out, xmin_out, xmax_out, y_out, dy_out = rebin_spectrum(
        # Positional Arguments:
        x_in, xmin_in, xmax_in, 
        spec_data, 
        # Keyword Arguments (from params dictionary):
        calc_unit=calc_unit, 
        **params # This includes xstart, xend, dx, kappa, filling
    )
    
    # Isolate the output bin we calculated (Index 1)
    y_test_value = y_out.value[1]
    dy_test_value = dy_out.value[1]

    # Extract the scalar expected value
    expected_y_value = expected_y.value[0]
    expected_dy_value = expected_dy.value[0]
    
    # Use np.isclose for scalar comparison (Astropy handles units internally)
    assert np.isclose(y_test_value, expected_y_value, atol=tolerance)
    assert np.isclose(dy_test_value, expected_dy_value, atol=tolerance)

# This first test ensures the core math for variance-weighted average is correct.