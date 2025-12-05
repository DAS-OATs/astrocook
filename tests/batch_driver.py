import argparse
import os
import sys
import logging
import subprocess 
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
    
    if len(x) == 0:
        return 0.0

    # Smooth to ignore noise spikes
    y_smooth = gaussian_filter1d(y, 2.0)
    
    # Avoid edges if possible
    if len(y_smooth) > 20:
        valid_slice = slice(10, -10)
        # Find index relative to the slice
        min_idx_rel = np.argmin(y_smooth[valid_slice]) 
        # Convert to absolute index
        min_idx = min_idx_rel + 10
    else:
        min_idx = np.argmin(y_smooth)
        
    lam_obs = x[min_idx]
    flux_min = y_smooth[min_idx]
    
    # Determine Rest Wavelength
    key = series
    
    if key in ATOM_DATA:
        lam_rest_ang = ATOM_DATA[key]['wave']
        # Loader converts Angstrom -> nm, so we must do the same
        lam_rest_nm = lam_rest_ang / 10.0
    else:
        lam_rest_nm = 121.567

    z_found = lam_obs / lam_rest_nm - 1.0
    
    logging.info(f"Auto-Detection ({series}): Deepest at {lam_obs:.2f} nm (Flux={flux_min:.2f}).")
    logging.info(f"  -> Inferred z = {z_found:.5f}")
    
    return z_found

def trim_to_velocity_window(session: SessionV2, transition: str, z_target: float, dv_kms: float = 500.0) -> SessionV2:
    """
    Uses the 'split' recipe to crop a session.
    """
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

def run_sightline_pipeline(data_dir: str, output_dir: str, los_id: str):
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
    context = HeadlessContext()
    
    files = {
        'HI': f'HI_spectra_los_{los_id}_snr10.txt',
        'OVI_1': f'OVI_doublet_los_{los_id}_snr10.txt',
        'OVI_2': f'OVI_doublet2_spectra_los_{los_id}_snr10.txt'
    }
    
    sessions = {}
    
    logging.info(f"--- Processing LOS {los_id} ---")
    for key, fname in files.items():
        path = os.path.join(data_dir, fname)
        if not os.path.exists(path):
            logging.error(f"File not found: {path}")
            return None
        sess = SessionV2.open_new(path, key, context, "ascii_resvel_header")
        if sess == 0: raise RuntimeError(f"Failed to load {fname}")
        sessions[key] = sess

    # --- 2. HIERARCHICAL DISCOVERY ---
    logging.info("--- Auto-Detecting Redshifts ---")
    
    # A. Anchor: Detect HI Redshift first
    z_hi = get_auto_redshift(sessions['HI'], 'Ly_a')
    
    # --- 3. TRIM SEGMENTS (Conditioned on z_hi) ---
    logging.info("--- Trimming ---")
    trim_window = 500.0 
    
    # Trim HI around its detected center
    s_hi_trim = trim_to_velocity_window(sessions['HI'], 'Ly_a', z_hi, trim_window)
    
    # Trim OVI segments around the HI center (This limits the search space!)
    s_ovi1_trim = trim_to_velocity_window(sessions['OVI_1'], 'OVI_1031', z_hi, trim_window)
    s_ovi2_trim = trim_to_velocity_window(sessions['OVI_2'], 'OVI_1037', z_hi, trim_window)

    # B. Refinement: Detect OVI Redshift inside the trimmed window
    # If the trim was empty (e.g. chunk didn't cover z_hi), fallback to z_hi
    if len(s_ovi1_trim.spec.x.value) > 0:
        z_ovi = get_auto_redshift(s_ovi1_trim, 'OVI_1031')
    else:
        logging.warning("OVI chunk seems empty after trimming. Falling back to z_hi.")
        z_ovi = z_hi
        
    dz_kms = 299792.458 * (z_ovi - z_hi) / (1 + z_hi)
    logging.info(f"Refined OVI z={z_ovi:.5f} (Offset from HI: {dz_kms:.1f} km/s)")

    # --- 4. STITCHING ---
    logging.info("--- Stitching Sessions ---")
    s_full = s_hi_trim.edit.stitch([s_ovi1_trim, s_ovi2_trim])

    # --- 5. MODELING ---
    logging.info("--- Modeling ---")
    # Place Ly-a at its detected z
    s_full = s_full.absorbers.add_component(series='Ly_a', z=z_hi, logN=14.0, b=25.0)
    uuid_lya = s_full.systs.components[-1].uuid
    
    # Place OVI at its REFINED z
    s_full = s_full.absorbers.add_component(series='OVI', z=z_ovi, logN=13.5, b=15.0)
    uuid_ovi = s_full.systs.components[-1].uuid

    #logging.info(f"Linking OVI b-parameter...")
    #s_full = s_full.absorbers.update_constraint(
    #    uuid=uuid_ovi, param='b', is_free=False, expression=f"p['{uuid_lya}'].b * 0.25"
    #)

    # --- 6. OPTIMIZATION ---
    logging.info("--- Optimization ---")
    s_final = s_full.absorbers.optimize_system(
        uuid=uuid_lya,          
        max_components=5,       
        threshold_sigma=2.0,    
        aic_penalty=0.0,        
        z_window_kms=str(trim_window),
        patience=2
    )
    
    if s_final == 0 or s_final is None:
        logging.error("Optimization failed.")
        return None

    # --- 7. EXPORT ---
    logging.info("--- Exporting ---")
    if not os.path.exists(output_dir): os.makedirs(output_dir)
        
    t_results = s_final.systs.t
    csv_path = os.path.join(output_dir, f"results_los_{los_id}.csv")
    t_results.write(csv_path, format='ascii.csv', overwrite=True)
    
    session_path = os.path.join(output_dir, f"session_los_{los_id}.acs2")
    s_final.save(session_path, initial_session=s_full)
    
    print(f"\nDone! Saved to: {output_dir}")
    print(f"z_HI: {z_hi:.5f}, z_OVI: {z_ovi:.5f}")
    print("Fitted Components:")
    t_results.pprint_all()
    print("="*40)
    
    return session_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Astrocook pipeline on a sightline.")
    parser.add_argument('los_id', type=str, help="ID of the Line of Sight (e.g., 203)")
    parser.add_argument('--data', type=str, default="/Users/guido/Downloads/HI_OVI_spectra", help="Data directory")
    parser.add_argument('--out', type=str, default="astrocook_output", help="Output directory")
    parser.add_argument('--inspect', action='store_true', help="Open the result in the GUI after processing.")
    
    args = parser.parse_args()
    
    saved_session = run_sightline_pipeline(args.data, args.out, args.los_id)
    
    if args.inspect and saved_session and os.path.exists(saved_session):
        logging.info(f"Launching Inspector for {saved_session}...")
        script_dir = os.path.dirname(os.path.abspath(__file__))
        gui_script = os.path.join(script_dir, "launch_pyside_app.py")
        if os.path.exists(gui_script):
            subprocess.run([sys.executable, gui_script, saved_session])
        else:
            logging.error(f"Could not find GUI launcher at {gui_script}")