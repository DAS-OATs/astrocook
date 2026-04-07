import sys
from astropy.io import fits
import numpy as np

def inspect_fits(file_path):
    print(f"\n{'='*60}")
    print(f"Inspecting: {file_path}")
    print(f"{'='*60}")
    
    try:
        with fits.open(file_path) as hdul:
            hdul.info()
            
            for i, hdu in enumerate(hdul):
                print(f"\n--- HDU {i}: {hdu.name} ---")
                header = hdu.header
                
                # Check for Image Data (1D/2D arrays)
                if hdu.is_image and hdu.data is not None:
                    print(f"Type: Image Data")
                    print(f"Shape: {hdu.data.shape}")
                    # Check for WCS Keywords (if it's a 1D spectrum in an image)
                    wcs_keys = ['CRVAL1', 'CDELT1', 'CRPIX1', 'CUNIT1']
                    found_keys = {k: header.get(k) for k in wcs_keys if k in header}
                    if found_keys:
                        print(f"WCS Keywords: {found_keys}")
                        
                # Check for Table Data (Binary Tables)
                elif isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                    print(f"Type: Table Data")
                    print(f"Columns: {hdu.columns.names}")
                    # Print units if available
                    if hasattr(hdu.columns, 'units'):
                        print(f"Units: {hdu.columns.units}")
                    
                # specific check for PypeIt/Custom signatures
                if 'FLUX' in header: print("Note: 'FLUX' found in header keys.")
                
    except Exception as e:
        print(f"ERROR reading file: {e}")

if __name__ == "__main__":
    # Replace these with your actual paths
    files = [
        "~/Downloads/spec1d_XSHOO.2025-01-28T07_51_24.382-P159-02_XShooter_VIS_20250128T075124.381.fits", 
        "~/Downloads/PSOJ159_02_VIS_coadd.fits",
        "~/Downloads/PSOJ159_02_VIS_coadd_orderstackVIS.fits"
    ]
    
    for f in files:
        inspect_fits(f)