import os
import sys
import logging
import numpy as np
from scipy.ndimage import gaussian_filter1d
import astropy.units as au
from astropy.table import Table

# Ensure we can import astrocook
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from astrocook.core.session import SessionV2
from astrocook.core.structures import SystemListDataV2
from astrocook.core.atomic_data import ATOM_DATA

class HeadlessContext:
    def __init__(self):
        self._flags = []
        self.session_histories = [] 
    def _flags_cond(self, flag): return False
    def _flags_extr(self, flag): return None
    def _launch_recipe_dialog(self, *args):
        logging.warning("Recipe tried to launch GUI dialog in headless mode.")

def get_auto_redshift(session: SessionV2, series: str = 'Ly_a') -> float:
    """
    Scans the spectrum to find the redshift of the strongest absorption feature.
    """
    x = session.spec.x.value
    y = session.spec.y.value
    
    # Smooth to ignore noise spikes
    y_smooth = gaussian_filter1d(y, 2.0)
    
    # Avoid edges
    if len(y_smooth) > 20:
        valid_slice = slice(10, -10)
        sl_y = y_smooth[valid_slice]
        sl_x = x[valid_slice]
        min_idx = np.argmin(sl_y)
        lam_obs = sl_x[min_idx]
        flux_min = sl_y[min_idx]
    else:
        min_idx = np.argmin(y_smooth)
        lam_obs = x[min_idx]
        flux_min = y_smooth[min_idx]
        
    # Determine Rest Wavelength
    # Check alias first
    key = series
    
    if key in ATOM_DATA:
        lam_rest_ang = ATOM_DATA[key]['wave']
        lam_rest_nm = lam_rest_ang / 10.0
    else:
        # Fallback default for Ly_a if lookup fails
        lam_rest_nm = 121.567

    z_found = lam_obs / lam_rest_nm - 1.0
    
    logging.info(f"Auto-Detection ({series}): Deepest at {lam_obs:.2f} nm (Flux={flux_min:.2f}).")
    logging.info(f"  -> Inferred z = {z_found:.5f}")
    
    return z_found

def trim_to_velocity_window(session: SessionV2, transition: str, z_target: float, dv_kms: float = 500.0) -> SessionV2:
    """
    Uses the 'split' recipe to crop a session.
    """
    # Normalize Key
    key = transition

    if key not in ATOM_DATA:
        logging.warning(f"Trim: Unknown transition '{key}'. Skipping trim.")
        return session

    lam_rest_ang = ATOM_DATA[key]['wave']
    lam_rest_nm = lam_rest_ang / 10.0
    
    lam_obs = lam_rest_nm * (1 + z_target)
    c_kms = 299792.458
    d_lam = (dv_kms / c_kms) * lam_obs
    
    lam_min = lam_obs - d_lam
    lam_max = lam_obs + d_lam
    
    logging.info(f"Trimming {transition} to {lam_min:.2f}-{lam_max:.2f} nm (z={z_target:.4f})")
    return session.edit.split(f"(x > {lam_min:.4f}) & (x < {lam_max:.4f})")

def run_sightline_pipeline(data_dir: str, output_dir: str):
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
    context = HeadlessContext()
    
    files = {
        'HI': 'HI_spectra_los_203_snr10.txt',
        'OVI_1': 'OVI_doublet_los_203_snr10.txt',
        'OVI_2': 'OVI_doublet2_spectra_los_203_snr10.txt'
    }
    
    sessions = {}
    
    logging.info("--- 1. Loading Spectra ---")
    for key, fname in files.items():
        path = os.path.join(data_dir, fname)
        if not os.path.exists(path):
            logging.error(f"File not found: {path}")
            return
        sess = SessionV2.open_new(path, key, context, "ascii_resvel_header")
        if sess == 0: raise RuntimeError(f"Failed to load {fname}")
        sessions[key] = sess

    # --- 2. AUTO-DISCOVERY (Separate for HI and OVI) ---
    logging.info("--- 2. Auto-Detecting Redshifts ---")
    
    # A. Detect HI z
    z_hi = get_auto_redshift(sessions['HI'], 'Ly_a')
    
    # B. Detect OVI z (Independent guess)
    z_ovi = get_auto_redshift(sessions['OVI_1'], 'OVI_1031')
    
    # Check consistency (just for log)
    dz_kms = 299792.458 * abs(z_hi - z_ovi) / (1 + z_hi)
    logging.info(f"Velocity offset between HI and OVI guesses: {dz_kms:.1f} km/s")

    # --- 3. TRIMMING ---
    trim_window = 500.0 # km/s
    
    s_hi_trim = trim_to_velocity_window(sessions['HI'], 'Ly_a', z_hi, trim_window)
    # Note: Use z_ovi for OVI trimming to center the window correctly
    s_ovi1_trim = trim_to_velocity_window(sessions['OVI_1'], 'OVI_1031', z_ovi, trim_window)
    s_ovi2_trim = trim_to_velocity_window(sessions['OVI_2'], 'OVI_1037', z_ovi, trim_window)

    # --- 4. STITCHING ---
    logging.info("--- 4. Stitching Sessions ---")
    s_full = s_hi_trim.edit.stitch([s_ovi1_trim, s_ovi2_trim])
    logging.info(f"Stitched Spectrum: {len(s_full.spec.x.value)} valid pixels.")

    # --- 5. MODELING (Seed) ---
    logging.info("--- 5. Defining Seed System ---")
    
    # Use z_hi for Ly-a component
    s_full = s_full.absorbers.add_component(series='Ly_a', z=z_hi, logN=14.0, b=25.0)
    uuid_lya = s_full.systs.components[-1].uuid

    # Use z_ovi for OVI component
    s_full = s_full.absorbers.add_component(series='OVI', z=z_ovi, logN=13.5, b=15.0)
    uuid_ovi = s_full.systs.components[-1].uuid

    # Link OVI b to Ly_a b
    """
    logging.info(f"Linking OVI b-parameter to Ly-a...")
    s_full = s_full.absorbers.update_constraint(
        uuid=uuid_ovi, param='b', is_free=False, 
        expression=f"p['{uuid_lya}'].b * 0.25"
    )
    """
    
    # --- 6. OPTIMIZATION ---
    logging.info("--- 6. Running Iterative Optimization ---")
    
    # We assume the HI line is the complex one requiring multicomponents.
    # We optimize around the Ly-a UUID.
    s_final = s_full.absorbers.optimize_system(
        uuid=uuid_lya,          
        max_components=5,       # Increased to allow complex profiles
        threshold_sigma=2.0,    # Lowered threshold (more sensitive)
        aic_penalty=0.0,        # Zero penalty (greedy fit: keep if chi2 improves at all)
        z_window_kms=str(trim_window)
    )
    
    # --- 7. EXPORT ---
    logging.info("--- 7. Exporting Results ---")
    if not os.path.exists(output_dir): os.makedirs(output_dir)
        
    t_results = s_final.systs.t
    csv_path = os.path.join(output_dir, "results_los_203_auto.csv")
    t_results.write(csv_path, format='ascii.csv', overwrite=True)
    
    session_path = os.path.join(output_dir, "session_los_203_auto.acs2")
    s_final.save(session_path, initial_session=s_full)
    
    print("\n" + "="*40)
    print(f"Done! Saved to: {output_dir}")
    print(f"Detected z(HI):  {z_hi:.5f}")
    print(f"Detected z(OVI): {z_ovi:.5f}")
    print("Fitted Components:")
    t_results.pprint_all()
    print("="*40)

if __name__ == "__main__":
    DATA_DIR = "/Users/guido/Downloads/HI_OVI_spectra" 
    OUTPUT_DIR = "astrocook_output"
    run_sightline_pipeline(DATA_DIR, OUTPUT_DIR)