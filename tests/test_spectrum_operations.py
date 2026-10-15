import numpy as np
from astropy import units as au
from astrocook.core.structures import SpectrumDataV2, DataColumnV2
from astrocook.core.spectrum_operations import rebin_spectrum

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
        'xstart': au.Quantity(100.025, au.nm), # Left edge of P2
        'xend': au.Quantity(100.125, au.nm),   # Right edge of P3
        'dx': au.Quantity(0.1, au.nm),         # One output bin covering P2 and P3
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

# --- Test Functions ---

def _rebin_case1():
    spec_data, params, expected_y, expected_dy = get_mock_input_case1()
    x_out, xmin_out, xmax_out, y_out, dy_out, resol_out = rebin_spectrum(
        spec_data.x.quantity, spec_data.xmin.quantity, spec_data.xmax.quantity,
        spec_data,
        calc_unit=spec_data.x.quantity.unit,
        **params  # xstart, xend, dx, kappa, filling
    )
    # A single output bin, centred at 100.075 nm
    assert len(y_out) == 1
    assert np.isclose(x_out.value[0], 100.075)
    return y_out, dy_out, expected_y, expected_dy


def test_rebin_weighted_average(tolerance=1e-9):
    y_out, _, expected_y, _ = _rebin_case1()
    assert np.isclose(y_out.value[0], expected_y.value[0], atol=tolerance)


def test_rebin_weighted_average_error(tolerance=1e-9):
    _, dy_out, _, expected_dy = _rebin_case1()
    assert np.isclose(dy_out.value[0], expected_dy.value[0], atol=tolerance)
