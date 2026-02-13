
import logging
import numpy as np
from astrocook.io.loaders import detect_file_format, get_loader

file_path = "/Users/guido/Library/CloudStorage/GoogleDrive-guido.cupani@inaf.it/My Drive/teaching/ac/spectrum_BD+14_3079_LabAstro23.fits"

# Configure logging
logging.basicConfig(level=logging.INFO)

print(f"Testing file: {file_path}")

# 1. Test Detection
fmt = detect_file_format(file_path)
print(f"Detected format: {fmt}")

if fmt == "primary_image_fits":
    # 2. Test Loading
    loader = get_loader(fmt)
    try:
        data = loader(file_path)
        print("Load successful!")
        print(f"X range: {data.x.min().value:.2f} - {data.x.max().value:.2f} {data.x.unit}")
        print(f"Y mean: {np.mean(data.y.value):.2f}")
        print(f"Data length: {len(data.x)}")
    except Exception as e:
        print(f"Load failed: {e}")
else:
    print("Format detection failed (expected 'primary_image_fits')")
