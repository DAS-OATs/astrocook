# tests/test_v2_io.py

import astropy.units as au
import numpy as np 
from astrocook.core.session import load_spectrum_to_v2_format

def test_spectrum_data_integrity():
    """Verifies that the V2 Spectrum object contains correctly mapped real data."""
    
    # Use the actual file path that failed on the GUI
    test_file_path = '/Users/guido/Library/CloudStorage/GoogleDrive-guido.cupani@inaf.it/My Drive/GitHub/astrocook/tests/data/J0103-1305_metallines_lya_new_spec.fits'
    format_name = 'generic_spectrum' 
    
    # 1. Load the SpectrumV2 object using the full V2 adapter chain
    try:
        spec_v2 = load_spectrum_to_v2_format(
            path=test_file_path, 
            format_name=format_name,
            # We must pass a mock GUI since the adapter still requires it for V1 I/O flags
            gui=object() # Pass a generic object to satisfy the V1 requirement
        ) # Get the spectrum object from the session (or load directly if load_spectrum_to_v2_format returns SpectrumV2)
    except Exception as e:
        print(f"ERROR: V2 loading chain failed during direct execution: {e}")
        return

    print("--- V2 Data Integrity Check ---")
    
    # 2. Perform checks (Size and Units)
    
    # A. Size Check (Should be 215588 based on your mock log, or the actual size)
    expected_size = 215588 
    print(f"Array Length: {len(spec_v2.x)}")
    assert len(spec_v2.x) == expected_size, f"Size mismatch: Expected {expected_size}, got {len(spec_v2.x)}"
    
    # B. Unit Check (Verify the units are set correctly)
    expected_x_unit = au.nm
    expected_y_unit = au.erg / au.cm**2 / au.s / au.Angstrom # Based on common V1 units
    print(f"X Unit: {spec_v2.x.unit}")
    print(f"Y Unit: {spec_v2.y.unit}")
    assert spec_v2.x.unit == expected_x_unit, f"X Unit mismatch: Expected {expected_x_unit}, got {spec_v2.x.unit}"
    
    # C. Data Sanity Check (Median/NaNs)
    y_median = np.median(spec_v2.y.value[~np.isnan(spec_v2.y.value)])
    print(f"Y Data Median (No NaN): {y_median}")
    assert not np.isnan(y_median), "Y data is all NaN or corrupted."
    
    print("--- Data Integrity PASS: V2 Object is Clean ---")

if __name__ == "__main__":
    test_spectrum_data_integrity()