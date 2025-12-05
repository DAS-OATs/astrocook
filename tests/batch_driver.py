import argparse
import os
import sys
import logging
import subprocess 
import numpy as np
from scipy.ndimage import gaussian_filter1d
from astropy.table import Table

# Ensure we can import astrocook
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from astrocook.core.session import SessionV2
from astrocook.core.atomic_data import ATOM_DATA

# --- CONFIGURATION: Supported Doublets ---
SPECIES_CONFIG = {
    'OVI': {
        'lines': ['OVI_1031', 'OVI_1037'],
        'series': 'OVI',
        'logN_init': 13.5,
        'b_init': 15.0
    },
    'CIV': {
        'lines': ['CIV_1548', 'CIV_1550'],
        'series': 'CIV',
        'logN_init': 13.5,
        'b_init': 12.0
    },
    'MgII': {
        'lines': ['MgII_2796', 'MgII_2803'],
        'series': 'MgII',
        'logN_init': 13.0,
        'b_init': 8.0
    }
}

class HeadlessContext:
    def __init__(self):
        self._flags = []
    def _flags_cond(self, flag): return False
    def _flags_extr(self, flag): return None
    def _launch_recipe_dialog(self, *args):
        logging.warning("Recipe tried to launch GUI dialog in headless mode.")

# --- HELPER 1: Velocity Logic ---
def get_velocity_bounds(session: SessionV2, z_target: float, trans: str):
    """Calculates min/max velocity (km/s) covered by a session relative to z_target."""
    if trans not in ATOM_DATA: return -np.inf, np.inf
    
    x = session.spec.x.value
    if len(x) == 0: return 0.0, 0.0
    
    lam_rest_nm = ATOM_DATA[trans]['wave'] / 10.0
    lam_obs = lam_rest_nm * (1 + z_target)
    c_kms = 299792.458
    
    v = c_kms * (x - lam_obs) / lam_obs
    return np.min(v), np.max(v)

def find_common_velocity_range(sessions_map: dict, z_target: float, requested_width: float) -> tuple:
    """
    Finds the INTERSECTION of velocity coverage across all chunks.
    Ensures we don't try to fit regions where data is missing.
    """
    v_mins = []
    v_maxs = []
    
    # 1. Calculate physical bounds of available data
    for key, sess in sessions_map.items():
        # Heuristic: Key usually contains series name (e.g. 'CIV_1', 'HI')
        # We need to guess the transition to convert lambda -> velocity
        # For simplicity in this batch script, we look up the primary line.
        
        trans = 'Ly_a' # Default
        if 'HI' not in key:
            # It's a metal. Find which one.
            for species, conf in SPECIES_CONFIG.items():
                if species in key:
                    # key 'CIV_1' -> index 0, 'CIV_2' -> index 1
                    idx = 0 if '1' in key else 1
                    trans = conf['lines'][idx]
                    break
        
        v_min, v_max = get_velocity_bounds(sess, z_target, trans)
        v_mins.append(v_min)
        v_maxs.append(v_max)
        
    # 2. Intersection
    common_min = max(v_mins)
    common_max = min(v_maxs)
    
    # 3. Apply User Request (Trimming)
    # Center the request around 0
    req_min = -requested_width
    req_max = requested_width
    
    final_min = max(common_min, req_min)
    final_max = min(common_max, req_max)
    
    if final_max <= final_min:
        logging.warning(f"No common velocity overlap! Data: [{common_min:.0f}, {common_max:.0f}], Req: [{req_min}, {req_max}]")
        return -requested_width, requested_width # Fallback
        
    logging.info(f"Trimming Intersection: [{final_min:.1f}, {final_max:.1f}] km/s")
    return final_min, final_max

def trim_by_velocity(session: SessionV2, transition: str, z_target: float, v_min: float, v_max: float) -> SessionV2:
    if transition not in ATOM_DATA: return session
    lam_rest_nm = ATOM_DATA[transition]['wave'] / 10.0
    lam_obs = lam_rest_nm * (1 + z_target)
    c_kms = 299792.458
    
    lam_lower = lam_obs * (1 + v_min / c_kms)
    lam_upper = lam_obs * (1 + v_max / c_kms)
    
    return session.edit.split(f"(x > {lam_lower:.4f}) & (x < {lam_upper:.4f})")

# --- HELPER 2: Redshift Detection ---
def get_auto_redshift(session: SessionV2, series: str) -> float:
    x = session.spec.x.value
    y = session.spec.y.value
    if len(x) == 0: return 0.0
    y_smooth = gaussian_filter1d(y, 2.0)
    
    if len(y_smooth) > 20:
        valid_slice = slice(10, -10)
        min_idx = np.argmin(y_smooth[valid_slice]) + 10
    else:
        min_idx = np.argmin(y_smooth)
        
    lam_obs = x[min_idx]
    if series in ATOM_DATA:
        lam_rest_nm = ATOM_DATA[series]['wave'] / 10.0
    else:
        lam_rest_nm = 121.567
    return lam_obs / lam_rest_nm - 1.0

def get_doublet_redshift(sess1, trans1, sess2, trans2):
    x1, y1 = sess1.spec.x.value, sess1.spec.y.value
    x2, y2 = sess2.spec.x.value, sess2.spec.y.value
    if len(x1)==0 or len(x2)==0: raise ValueError("Empty spectrum")

    lam1 = ATOM_DATA[trans1]['wave'] / 10.0
    lam2 = ATOM_DATA[trans2]['wave'] / 10.0
    z1 = x1 / lam1 - 1.0
    z2 = x2 / lam2 - 1.0

    z_min, z_max = max(z1.min(), z2.min()), min(z1.max(), z2.max())
    if z_max <= z_min: raise ValueError("No redshift overlap")

    dz = np.median(np.diff(z1))
    z_grid = np.arange(z_min, z_max, dz)
    
    y1_int = np.interp(z_grid, z1, y1)
    y2_int = np.interp(z_grid, z2, y2)
    
    depth1 = np.maximum(0.0, 1.0 - gaussian_filter1d(y1_int, 2.0))
    depth2 = np.maximum(0.0, 1.0 - gaussian_filter1d(y2_int, 2.0))
    
    score = depth1 * depth2
    best_idx = np.argmax(score)
    return z_grid[best_idx]

# --- PIPELINE ---
def run_sightline_pipeline(data_dir: str, output_dir: str, los_id: str, 
                           species: str = 'OVI', trim_window: float = 500.0):
    
    # 1. Setup Configuration
    if species not in SPECIES_CONFIG:
        logging.error(f"Unknown species: {species}")
        return None
    
    conf = SPECIES_CONFIG[species]
    metal_key_1, metal_key_2 = conf['lines']
    metal_series = conf['series']
    
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
    context = HeadlessContext()
    
    # Dynamic Filenames based on Species
    files = {
        'HI': f'HI_spectra_los_{los_id}_snr10.txt',
        'METAL_1': f'{species}_doublet_los_{los_id}_snr10.txt',
        'METAL_2': f'{species}_doublet2_spectra_los_{los_id}_snr10.txt'
    }
    
    sessions = {}
    logging.info(f"--- Processing LOS {los_id} ({species}) ---")
    
    for key, fname in files.items():
        path = os.path.join(data_dir, fname)
        if not os.path.exists(path):
            logging.error(f"File not found: {path}")
            return None
        # Assuming Loader can handle generic columns
        sess = SessionV2.open_new(path, key, context, "ascii_resvel_header")
        if sess == 0: raise RuntimeError(f"Failed to load {fname}")
        sessions[key] = sess

    # --- 2. ANCHOR DETECTION ---
    logging.info(f"--- Auto-Detecting {species} Doublet ---")
    try:
        z_anchor = get_doublet_redshift(
            sessions['METAL_1'], metal_key_1,
            sessions['METAL_2'], metal_key_2
        )
        logging.info(f"Anchor found at z={z_anchor:.5f}")
    except Exception as e:
        logging.error(f"Failed to detect doublet: {e}")
        return None
    
    # --- 3. ROBUST TRIMMING ---
    logging.info("--- Calculating Trim Intersection ---")
    
    # Calculate the valid intersection of all chunks + user window
    v_min, v_max = find_common_velocity_range(sessions, z_anchor, trim_window)
    
    s_metal1_trim = trim_by_velocity(sessions['METAL_1'], metal_key_1, z_anchor, v_min, v_max)
    s_metal2_trim = trim_by_velocity(sessions['METAL_2'], metal_key_2, z_anchor, v_min, v_max)
    s_hi_trim     = trim_by_velocity(sessions['HI'],      'Ly_a',      z_anchor, v_min, v_max)

    # --- 4. HI REFINEMENT ---
    # Find strongest HI inside the (now strictly valid) trimmed window
    if len(s_hi_trim.spec.x.value) > 5:
        z_hi = get_auto_redshift(s_hi_trim, 'Ly_a')
        logging.info(f"Refined HI z={z_hi:.5f}")
    else:
        logging.warning("HI chunk empty after intersection trim. Using anchor z.")
        z_hi = z_anchor

    # --- 5. STITCH & MODEL ---
    s_full = s_hi_trim.edit.stitch([s_metal1_trim, s_metal2_trim])

    # Place Components
    # 1. Anchor (Metal)
    s_full = s_full.absorbers.add_component(
        series=metal_series, z=z_anchor, 
        logN=conf['logN_init'], b=conf['b_init']
    )
    uuid_anchor = s_full.systs.components[-1].uuid

    # 2. Target (HI)
    s_full = s_full.absorbers.add_component(
        series='Ly_a', z=z_hi, logN=14.0, b=25.0
    )
    uuid_lya = s_full.systs.components[-1].uuid

    # --- 6. OPTIMIZATION ---
    logging.info(f"--- Optimizing {species} ---")
    
    # Step A: Optimize Metal Anchor
    s_full = s_full.absorbers.optimize_system(
        uuid=uuid_anchor,
        max_components=3,
        threshold_sigma=2.5,
        z_window_kms=str(trim_window), # This now serves as the search window
        patience=1
    )

    # Step B: Optimize HI
    logging.info("--- Optimizing HI ---")
    s_final = s_full.absorbers.optimize_system(
        uuid=uuid_lya,          
        max_components=5,       
        threshold_sigma=2.0,    
        aic_penalty=0.0,        
        z_window_kms=str(trim_window),
        patience=2
    )

    if not s_final: return None

    # --- 7. EXPORT ---
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    csv_path = os.path.join(output_dir, f"results_los_{los_id}_{species}.csv")
    s_final.systs.t.write(csv_path, format='ascii.csv', overwrite=True)
    
    sess_path = os.path.join(output_dir, f"session_los_{los_id}_{species}.acs2")
    s_final.save(sess_path, initial_session=s_full)
    
    print(f"Done! {species} z={z_anchor:.5f}, HI z={z_hi:.5f}")
    return sess_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('los_id', type=str)
    parser.add_argument('--species', type=str, default='OVI', help="OVI, CIV, or MgII")
    parser.add_argument('--trim', type=float, default=500.0, help="Window half-width (km/s)")
    parser.add_argument('--data', type=str, default=".")
    parser.add_argument('--out', type=str, default="output")
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    
    saved = run_sightline_pipeline(args.data, args.out, args.los_id, args.species, args.trim)
    
    if args.inspect and saved and os.path.exists(saved):
        script = os.path.join(os.path.dirname(__file__), "launch_pyside_app.py")
        subprocess.run([sys.executable, script, saved])