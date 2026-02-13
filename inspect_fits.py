
import os
from astropy.io import fits

file_path = "/Users/guido/Library/CloudStorage/GoogleDrive-guido.cupani@inaf.it/My Drive/teaching/ac/spectrum_BD+14_3079_LabAstro23.fits"

try:
    with fits.open(file_path) as hdul:
        print(f"Info:")
        hdul.info()
        
        h0 = hdul[0].header
        print(f"\nPrimary Header Keys:")
        keys = ['NAXIS', 'NAXIS1', 'NAXIS2', 'NAXIS3', 'CRVAL1', 'CDELT1', 'CD1_1', 'CRPIX1', 'CTYPE1', 'CUNIT1', 'BUNIT', 'TELESCOP', 'INSTRUME', 'OBJECT']
        for k in keys:
            if k in h0:
                print(f"{k} = {h0[k]}")
        
        if hdul[0].data is not None:
            print(f"\nData Shape: {hdul[0].data.shape}")
            print(f"Data Type: {hdul[0].data.dtype}")
            print(f"Min/Max: {hdul[0].data.min()}/{hdul[0].data.max()}")
            
except Exception as e:
    print(f"Error: {e}")
