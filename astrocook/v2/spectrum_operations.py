from ..v1.functions import x_convert as x_convert_v1 
from ..v2.structures import DataColumnV2, SpectrumDataV2 
from astropy import units as au
from typing import Optional, Tuple
import numpy as np

def convert_x_axis(spec_data: SpectrumDataV2, zem: float, xunit: au.Unit) -> Tuple[au.Quantity, au.Quantity, au.Quantity]:
    """
    Performs X-axis conversion based on V1 logic and returns new Quantity arrays.
    """
    
    # 1. Get Quantities from V2 Data
    x = spec_data.x.quantity
    xmin = spec_data.xmin.quantity
    xmax = spec_data.xmax.quantity
    
    # 2. Apply V1 Conversion Logic (assumes V1 x_convert utility works on Quantities)
    x_new = x_convert_v1(x, zem, xunit)
    xmin_new = x_convert_v1(xmin, zem, xunit)
    xmax_new = x_convert_v1(xmax, zem, xunit)
    
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

def rebin_spectrum(spec_data: SpectrumDataV2, xstart: Optional[au.Quantity], 
                   xend: Optional[au.Quantity], dx: au.Quantity, 
                   kappa: Optional[float], filling: float) -> Tuple[au.Quantity, au.Quantity, au.Quantity, au.Quantity, au.Quantity]:
    """
    Pure function to rebin spectrum data. MUST contain the logic from Spectrum._rebin (V1).
    Returns new arrays for x, xmin, xmax, y, dy.
    """
    # NOTE: This function currently represents technical debt. 
    # It must eventually contain all the complex logic from Spectrum._rebin (V1).
    # For now, we will use a mock to verify the V2 control flow.

    logging.warning("REBIN logic is using a V2 control flow mock. Implementation of V1 math needed.")
    
    # --- MOCK IMPLEMENTATION (To unblock V2 flow) ---
    x = spec_data.x.quantity
    y = spec_data.y.quantity
    dy = spec_data.dy.quantity
    
    # Return a dummy slice to prove immutability works
    return x[::2], x[::2], x[::2], y[::2], dy[::2]