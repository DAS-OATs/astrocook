import argparse
import os
import sys
import logging
import subprocess 
import numpy as np

# Ensure we can import astrocook
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from astrocook.core.session import SessionV2
from astrocook.core.atomic_data import ATOM_DATA

# --- CONFIGURATION ---
SPECIES_CONFIG = {
    'OVI': {
        'lines': ['OVI_1031', 'OVI_1037'],
        'series': 'OVI',
        'logN_init': 13.5, 'b_init': 15.0
    },
    'CIV': {
        'lines': ['CIV_1548', 'CIV_1550'],
        'series': 'CIV',
        'logN_init': 13.5, 'b_init': 12.0
    },
    'MgII': {
        'lines': ['MgII_2796', 'MgII_2803'],
        'series': 'MgII',
        'logN_init': 13.0, 'b_init': 8.0
    }
}

class HeadlessContext:
    def __init__(self): self._flags = []
    def _flags_cond(self, flag): return False
    def _flags_extr(self, flag): return None
    def _launch_recipe_dialog(self, *args): pass

def run_sightline_pipeline(data_dir: str, output_dir: str, los_id: str, 
                           species: str = 'OVI', trim_window: float = 500.0):
    
    # 1. Setup
    if species not in SPECIES_CONFIG:
        logging.error(f"Unknown species: {species}")
        return None
    
    conf = SPECIES_CONFIG[species]
    metal_key_1, metal_key_2 = conf['lines']
    
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
    context = HeadlessContext()
    
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
        sess = SessionV2.open_new(path, key, context, "ascii_resvel_header")
        if sess == 0: raise RuntimeError(f"Failed to load {fname}")
        sessions[key] = sess

    # --- 2. ANCHOR DETECTION (Using Core API) ---
    logging.info(f"--- Detecting {species} Anchor ---")
    # We call the Core method directly because we have Objects, not GUI strings
    z_anchor = sessions['METAL_1'].spec.detect_doublet_z(
        sessions['METAL_2'].spec, metal_key_1, metal_key_2
    )
    
    if z_anchor <= 0:
        logging.error("Could not find a valid doublet overlap.")
        return None
    logging.info(f"Anchor found at z={z_anchor:.5f}")
    
    # --- 3. ROBUST TRIMMING (Using Core API + Edit Recipe) ---
    logging.info("--- Trimming to Common Intersection ---")
    
    # Define the group of spectra to intersect
    group_specs = [sessions['METAL_1'].spec, sessions['METAL_2'].spec, sessions['HI'].spec]
    group_trans = [metal_key_1, metal_key_2, 'Ly_a']
    
    # Calculate the split expression for EACH session
    # (They share the velocity window, but have different wavelengths)
    trimmed_sessions = {}
    
    for key, sess in sessions.items():
        trans = metal_key_1 if 'METAL_1' in key else (metal_key_2 if 'METAL_2' in key else 'Ly_a')
        
        # Calculate intersection expression (Core API)
        expr = sess.spec.calc_intersection_expression(
            others_specs=[s for s in group_specs if s is not sess.spec],
            z_target=z_anchor,
            trans_self=trans,
            trans_others=[t for s, t in zip(group_specs, group_trans) if s is not sess.spec],
            window_kms=trim_window
        )
        
        if not expr:
            logging.warning(f"No overlap for {key}. Skipping LOS.")
            return None
            
        # Apply Split (Edit Recipe)
        trimmed_sessions[key] = sess.edit.split(expr)

    # --- 4. STITCHING (Using Edit Recipe) ---
    # The 'stitch' recipe now accepts a list of Session objects
    s_full = trimmed_sessions['HI'].edit.stitch([
        trimmed_sessions['METAL_1'], 
        trimmed_sessions['METAL_2']
    ])

    # --- 5. MODELING & OPTIMIZATION (Using Absorbers Recipe) ---
    logging.info("--- Modeling ---")
    
    # A. Place Anchor (Metal)
    s_full = s_full.absorbers.add_component(
        series=conf['series'], z=z_anchor, 
        logN=conf['logN_init'], b=conf['b_init']
    )
    uuid_anchor = s_full.systs.components[-1].uuid

    # B. Place HI (Seed at Anchor Z)
    s_full = s_full.absorbers.add_component(
        series='Ly_a', z=z_anchor, logN=14.0, b=25.0
    )
    uuid_lya = s_full.systs.components[-1].uuid

    # C. Optimize Metal
    logging.info(f"--- Optimizing {species} ---")
    s_full = s_full.absorbers.optimize_system(
        uuid=uuid_anchor,
        max_components=3,
        threshold_sigma=2.5,
        z_window_kms=str(trim_window),
        patience=1
    )

    # D. Optimize HI
    # optimize_system will find the strongest residual (the actual HI peak)
    # even if it's slightly offset from z_anchor.
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

    # --- 6. EXPORT ---
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    
    csv_name = f"results_los_{los_id}_{species}.csv"
    s_final.systs.t.write(os.path.join(output_dir, csv_name), format='ascii.csv', overwrite=True)
    
    acs_name = f"session_los_{los_id}_{species}.acs2"
    s_final.save(os.path.join(output_dir, acs_name), initial_session=s_full)
    
    print(f"Done! {species} z={z_anchor:.5f}")
    return os.path.join(output_dir, acs_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('los_id', type=str)
    parser.add_argument('--species', type=str, default='OVI')
    parser.add_argument('--trim', type=float, default=500.0)
    # Restore user default path
    parser.add_argument('--data', type=str, default="/Users/guido/Downloads/HI_OVI_spectra")
    parser.add_argument('--out', type=str, default="output")
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    
    saved = run_sightline_pipeline(args.data, args.out, args.los_id, args.species, args.trim)
    
    if args.inspect and saved and os.path.exists(saved):
        script = os.path.join(os.path.dirname(__file__), "launch_pyside_app.py")
        subprocess.run([sys.executable, script, saved])