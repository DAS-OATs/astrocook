
import sys
try:
    print("Importing loaders...", file=sys.stderr)
    from astrocook.io.loaders import detect_file_format, get_loader
    print("Import successful.", file=sys.stderr)

    file_path = "/Users/guido/Library/CloudStorage/GoogleDrive-guido.cupani@inaf.it/My Drive/teaching/ac/spectrum_BD+14_3079_LabAstro23.fits"
    
    print(f"Detecting format for {file_path}...", file=sys.stderr)
    fmt = detect_file_format(file_path)
    print(f"Detected format: {fmt}", file=sys.stderr)
    
    if fmt == "primary_image_fits":
        loader = get_loader(fmt)
        print("Loading data...", file=sys.stderr)
        data = loader(file_path)
        print("Success!", file=sys.stderr)
        print(f"Stats: X={len(data.x.value)} pts, Range={data.x.value.min():.2f}-{data.x.value.max():.2f}", file=sys.stderr)
        
except Exception as e:
    print(f"CRITICAL ERROR: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
