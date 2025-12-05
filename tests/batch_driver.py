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
    
    if len(x) == 0: return 0.0

    y_smooth = gaussian_filter1d(y, 2.0)
    
    if len(y_smooth) > 20:
        valid_slice = slice(10, -10)
        min_idx_rel = np.argmin(y_smooth[valid_slice]) 
        min_idx = min_idx_rel + 10
    else:
        min_idx = np.argmin(y_smooth)
        
    lam_obs = x[min_idx]
    
    # Determine Rest Wavelength
    if series in ATOM_DATA:
        lam_rest_nm = ATOM_DATA[series]['wave'] / 10.0
    else:
        lam_rest_nm = 121.567

    z_found = lam_obs / lam_rest_nm - 1.0
    return z_found

def get_doublet_redshift(sess1: SessionV2, trans1: str, sess2: SessionV2, trans2: str) -> float:
    """
    Finds the redshift where BOTH transitions show simultaneous absorption.
    Calculates the 'Correlation Score' = Depth1 * Depth2.
    """
    # 1. Get Data
    x1, y1 = sess1.spec.x.value, sess1.spec.y.value
    x2, y2 = sess2.spec.x.value, sess2.spec.y.value
    
    if len(x1) == 0 or len(x2) == 0:
        raise ValueError("One of the doublet spectra is empty.")

    # 2. Get Rest Wavelengths (nm)
    lam1 = ATOM_DATA[trans1]['wave'] / 10.0
    lam2 = ATOM_DATA[trans2]['wave'] / 10.0

    # 3. Convert both to Redshift Space
    z1 = x1 / lam1 - 1.0
    z2 = x2 / lam2 - 1.0

    # 4. Create Common Grid (Overlap Region)
    z_min = max(z1.min(), z2.min())
    z_max = min(z1.max(), z2.max())
    
    if z_max <= z_min:
        raise ValueError(f"No spectral overlap between {trans1} and {trans2}.")

    # Use resolution of first spectrum
    dz = np.median(np.diff(z1))
    z_grid = np.arange(z_min, z_max, dz)
    
    # 5. Interpolate Fluxes
    y1_int = np.interp(z_grid, z1, y1)
    y2_int = np.interp(z_grid, z2, y2)
    
    # Smooth to suppress noise
    y1_s = gaussian_filter1d(y1_int, 2.0)
    y2_s = gaussian_filter1d(y2_int, 2.0)

    # 6. Compute Joint Signal (Product of Depths)
    # Depth = 1 - Flux. We clip at 0 (ignore emission)
    depth1 = np.maximum(0.0, 1.0 - y1_s)
    depth2 = np.maximum(0.0, 1.0 - y2_s)
    
    score = depth1 * depth2
    
    best_idx = np.argmax(score)
    z_best = z_grid[best_idx]
    
    logging.info(f"Doublet Scan: Max Coincidence at z={z_best:.5f} (Score={score[best_idx]:.3f})")
    return z_best

def trim_to_velocity_window(session: SessionV2, transition: str, z_target: float, dv_kms: float = 500.0) -> SessionV2:
    key = transition
    if key not in ATOM_DATA: return session

    lam_rest_nm = ATOM_DATA[key]['wave'] / 10.0
    lam_obs = lam_rest_nm * (1 + z_target)
    c_kms = 299792.458
    d_lam = (dv_kms / c_kms) * lam_obs
    
    lam_min = lam_obs - d_lam
    lam_max = lam_obs + d_lam
    
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

    # --- 2. HIERARCHICAL DISCOVERY (NEW STRATEGY) ---
    # A. Anchor: Detect OVI Redshift first (Doublet Coincidence)
    logging.info("--- Auto-Detecting OVI Doublet ---")
    try:
        z_ovi = get_doublet_redshift(sessions['OVI_1'], 'OVI_1031',
                                     sessions['OVI_2'], 'OVI_1037')
    except Exception as e:
        logging.error(f"Failed to detect doublet: {e}")
        return None
    
    # --- 3. TRIM SEGMENTS (Conditioned on z_ovi) ---
    logging.info("--- Trimming ---")
    trim_window = 500.0 
    
    # Trim EVERYTHING based on OVI redshift
    s_ovi1_trim = trim_to_velocity_window(sessions['OVI_1'], 'OVI_1031', z_ovi, trim_window)
    s_ovi2_trim = trim_to_velocity_window(sessions['OVI_2'], 'OVI_1037', z_ovi, trim_window)
    s_hi_trim   = trim_to_velocity_window(sessions['HI'],    'Ly_a',     z_ovi, trim_window)

    # B. Refinement: Detect HI Redshift inside the trimmed window
    # We look for the strongest Ly_a line that is CLOSE to the OVI system
    if len(s_hi_trim.spec.x.value) > 0:
        z_hi = get_auto_redshift(s_hi_trim, 'Ly_a')
        dz_kms = 299792.458 * (z_hi - z_ovi) / (1 + z_ovi)
        logging.info(f"Refined HI z={z_hi:.5f} (Offset from OVI: {dz_kms:.1f} km/s)")
    else:
        logging.warning("HI chunk seems empty after trimming. Falling back to z_ovi.")
        z_hi = z_ovi

    # --- 4. STITCHING ---
    logging.info("--- Stitching Sessions ---")
    s_full = s_hi_trim.edit.stitch([s_ovi1_trim, s_ovi2_trim])

    # --- 5. MODELING ---
    logging.info("--- Modeling ---")
    
    # 1. Place OVI (The Anchor)
    s_full = s_full.absorbers.add_component(series='OVI', z=z_ovi, logN=13.5, b=15.0)
    uuid_ovi = s_full.systs.components[-1].uuid

    # 2. Place Ly-a (The Target)
    s_full = s_full.absorbers.add_component(series='Ly_a', z=z_hi, logN=14.0, b=25.0)
    uuid_lya = s_full.systs.components[-1].uuid

    # --- 6. OPTIMIZATION ---
    logging.info("--- Optimization ---")
    
    # [STEP A] Optimize OVI Anchor
    # We now use optimize_system instead of just fit_component. 
    # This allows it to add components if the metal system is complex/clumpy.
    logging.info("Optimizing OVI Anchor...")
    s_full = s_full.absorbers.optimize_system(
        uuid=uuid_ovi,
        max_components=3,      # OVI is usually simpler, keep low
        threshold_sigma=2.5,
        z_window_kms=str(trim_window),
        patience=1
    )

    # [STEP B] Optimize Ly_a Structure
    # OVI is now optimized (and unfrozen, thanks to absorbers.py fix).
    # Since they are not grouped, OVI parameters will stay fixed during Ly_a fit.
    logging.info("Optimizing Ly_a Structure...")
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
    print(f"z_OVI (Anchor): {z_ovi:.5f}, z_HI: {z_hi:.5f}")
    
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