import os
import sys
import logging
import numpy as np
from scipy.ndimage import gaussian_filter1d
import astropy.units as au

# GUI Imports for Debugging
from PySide6.QtWidgets import QApplication
from astrocook.gui.main_window import MainWindowV2
from astrocook.gui.debug_utils import GLOBAL_PLOTTER

# Ensure we can import astrocook
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from astrocook.core.session import SessionV2
from astrocook.core.atomic_data import ATOM_DATA

# --- 1. Helper Functions ---

def get_auto_redshift(session: SessionV2, series: str = 'Ly_a') -> float:
    """Scans the spectrum to find the redshift of the strongest absorption feature."""
    x = session.spec.x.value
    y = session.spec.y.value
    
    y_smooth = gaussian_filter1d(y, 2.0)
    
    if len(y_smooth) > 20:
        valid_slice = slice(10, -10)
        min_idx = np.argmin(y_smooth[valid_slice]) + 10
    else:
        min_idx = np.argmin(y_smooth)
        
    lam_obs = x[min_idx]
    flux_min = y_smooth[min_idx]
    
    key = series
    
    if key in ATOM_DATA:
        lam_rest_ang = ATOM_DATA[key]['wave']
        lam_rest_nm = lam_rest_ang / 10.0
    else:
        lam_rest_nm = 121.567

    z_found = lam_obs / lam_rest_nm - 1.0
    
    logging.info(f"Auto-Detection ({series}): Deepest at {lam_obs:.2f} nm (Flux={flux_min:.2f}).")
    logging.info(f"  -> Inferred z = {z_found:.5f}")
    return z_found

def trim_to_velocity_window(session: SessionV2, transition: str, z_target: float, dv_kms: float = 500.0) -> SessionV2:
    """Uses the 'split' recipe to crop a session."""
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

# --- 2. The Pipeline Logic ---

def run_sightline_pipeline(data_dir: str, output_dir: str):
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
    
    # --- GUI CONTEXT SETUP ---
    # We initialize the Qt App so signals work
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    
    # We create a hidden MainWindow to act as the "Brain"
    # It handles logging, signals, and recipe orchestration
    gui_context = MainWindowV2(initial_session=None, initial_log_object=None)
    
    # [CRITICAL] Connect the plotter signal
    # When absorbers.py emits a plot dict, MainWindow will catch it and show plt.show()
    GLOBAL_PLOTTER.plot_requested.connect(gui_context._on_debug_plot_requested)
    
    # -------------------------
    
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
        # Pass the REAL gui_context
        sess = SessionV2.open_new(path, key, gui_context, "ascii_resvel_header")
        if sess == 0: raise RuntimeError(f"Failed to load {fname}")
        sessions[key] = sess

    logging.info("--- 2. Auto-Detecting Redshifts ---")
    z_hi = get_auto_redshift(sessions['HI'], 'Ly_a')
    z_ovi = get_auto_redshift(sessions['OVI_1'], 'OVI_1031')
    
    dz_kms = 299792.458 * abs(z_hi - z_ovi) / (1 + z_hi)
    logging.info(f"Velocity offset between HI and OVI guesses: {dz_kms:.1f} km/s")

    # --- 3. TRIMMING ---
    trim_window = 500.0 
    s_hi_trim = trim_to_velocity_window(sessions['HI'], 'Ly_a', z_hi, trim_window)
    s_ovi1_trim = trim_to_velocity_window(sessions['OVI_1'], 'OVI_1031', z_ovi, trim_window)
    s_ovi2_trim = trim_to_velocity_window(sessions['OVI_2'], 'OVI_1037', z_ovi, trim_window)

    # --- 4. STITCHING ---
    logging.info("--- 4. Stitching Sessions ---")
    # Use the strings list or objects list? Objects are safer here.
    s_full = s_hi_trim.edit.stitch([s_ovi1_trim, s_ovi2_trim])
    logging.info(f"Stitched Spectrum: {len(s_full.spec.x.value)} valid pixels.")

    # --- 5. MODELING (Seed) ---
    logging.info("--- 5. Defining Seed System ---")
    s_full = s_full.absorbers.add_component(series='Ly_a', z=z_hi, logN=14.0, b=25.0)
    uuid_lya = s_full.systs.components[-1].uuid

    s_full = s_full.absorbers.add_component(series='OVI', z=z_ovi, logN=13.5, b=15.0)
    uuid_ovi = s_full.systs.components[-1].uuid

    # --- 6. OPTIMIZATION ---
    logging.info("--- 6. Running Iterative Optimization ---")
    logging.info("NOTE: Matplotlib windows will pop up. Close them to proceed to the next step.")
    
    s_final = s_full.absorbers.optimize_system(
        uuid=uuid_lya,          
        max_components=5,       
        threshold_sigma=2.0,    
        aic_penalty=0.0,        
        z_window_kms=str(trim_window),
        patience=2
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
    print("Fitted Components:")
    t_results.pprint_all()
    print("="*40)

if __name__ == "__main__":
    DATA_DIR = "/Users/guido/Downloads/HI_OVI_spectra" 
    OUTPUT_DIR = "astrocook_output"
    run_sightline_pipeline(DATA_DIR, OUTPUT_DIR)