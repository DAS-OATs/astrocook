import argparse
import os
import sys
import glob
import logging
import subprocess
import dataclasses
import numpy as np
import astropy.units as au
from scipy.ndimage import gaussian_filter1d

# Boilerplate to ensure we can import astrocook from parent directory
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import astrocook # Import package to locate run_astrocook.py
from astrocook.core.session import SessionV2
from astrocook.core.structures import DataColumnV2
from astrocook.core.atomic_data import STANDARD_MULTIPLETS, ATOM_DATA
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.core.system_list import SystemListV2

# --- CONFIGURATION ---
TARGET_ION = 'CII_1334'

# --- PATCH 1: ATOMIC DATA ---
if 'CII' in STANDARD_MULTIPLETS:
    del STANDARD_MULTIPLETS['CII']

# --- PATCH 2: OPTIMIZE HIERARCHY (Robust Version) ---
def patched_optimize_hierarchy(self, spec, uuid_seed, 
                               max_components=5, threshold_sigma=2.5,
                               aic_penalty=0.0, z_window_kms=100.0,
                               min_dv=10.0, group_depth=2, patience=2):
    
    current_systs = self
    target_comp = current_systs.get_component_by_uuid(uuid_seed)
    if not target_comp: return self

    series_name = target_comp.series
    lam_target_ang = ATOM_DATA[series_name]['wave']
    lam_target_nm = lam_target_ang / 10.0
    c_kms = 299792.458

    lam_obs_center = lam_target_nm * (1 + target_comp.z)
    dv_arr = c_kms * (spec.x.value - lam_obs_center) / lam_obs_center
    static_mask = (np.abs(dv_arr) < z_window_kms) & np.isfinite(spec.y.value) & (spec.dy.value > 0)
    
    if np.sum(static_mask) == 0: return self

    rejected_mask = np.zeros_like(spec.x.value, dtype=bool)

    def fit_and_evaluate(systs_obj, target_uuids):
        with systs_obj.fitting_context(target_uuids, group_depth=group_depth):
            fitter = VoigtFitterV2(spec, systs_obj)
            new_systs, _, res = fitter.fit(max_nfev=3000, z_window_kms=z_window_kms, verbose=0)
            
            if new_systs.constraint_model:
                new_systs.constraint_model.set_active_components(None)
            
            if not res: return float('inf'), systs_obj 

            eval_fitter = VoigtFitterV2(spec, new_systs)
            _, model_flux = eval_fitter.compute_model_flux()
            
            y_win = spec.y.value[static_mask]
            mod_win = model_flux[static_mask]
            dy_win = spec.dy.value[static_mask]
            
            chi2_static = np.sum(((y_win - mod_win) / dy_win)**2)
            k = len(res.x)
            aic = chi2_static + 2 * k
            return aic, new_systs

    best_aic, best_systs = fit_and_evaluate(current_systs, [uuid_seed])
    logging.info(f"Baseline AIC: {best_aic:.1f}")
    
    current_systs = best_systs
    trials_without_imp = 0
    
    for i in range(max_components):
        eval_fitter = VoigtFitterV2(spec, best_systs)
        _, model_phys = eval_fitter.compute_model_flux()
        
        with np.errstate(divide='ignore', invalid='ignore'):
            sigma_resid = (spec.y.value - model_phys) / spec.dy.value
        
        sigma_resid[~static_mask] = 0.0
        sigma_resid[rejected_mask] = 0.0 
        
        min_idx = np.argmin(sigma_resid)
        min_val = sigma_resid[min_idx]
        
        if min_val > -threshold_sigma:
            logging.info(f"Converged. Max resid {min_val:.1f}s > {threshold_sigma}")
            break
            
        z_new = spec.x.value[min_idx] / lam_target_nm - 1.0
        
        too_close = False
        for c in best_systs.components:
            if c.series == series_name:
                dz = abs(c.z - z_new)
                dv_sep = c_kms * dz / (1+z_new)
                if dv_sep < min_dv: 
                    too_close = True; break
        
        if too_close:
            sl = slice(max(0, min_idx-2), min(len(spec.x.value), min_idx+3))
            rejected_mask[sl] = True
            continue 

        logging.info(f"Iter {i+1}: Candidate at {z_new:.5f} (Sigma={min_val:.1f})")
        
        cand_systs = best_systs.add_component(series=series_name, z=z_new, logN=13.0, b=10.0)
        new_uuid = cand_systs.components[-1].uuid
        
        trial_aic, trial_systs = fit_and_evaluate(cand_systs, [uuid_seed, new_uuid])
        
        delta = trial_aic - best_aic
        if delta < -aic_penalty:
            logging.info(f"  -> Accepted. AIC {best_aic:.1f} -> {trial_aic:.1f} (Delta: {delta:.1f})")
            best_aic = trial_aic
            best_systs = trial_systs
            trials_without_imp = 0
        else:
            logging.info(f"  -> Rejected. AIC Delta {delta:.1f}")
            trials_without_imp += 1
            sl = slice(max(0, min_idx-2), min(len(spec.x.value), min_idx+3))
            rejected_mask[sl] = True
        
        if trials_without_imp >= patience:
            logging.info("Patience exhausted.")
            break
    
    return best_systs

SystemListV2.optimize_hierarchy = patched_optimize_hierarchy

# --- PATCH 3: FIX VOIGT FITTER EDGE CRASH ---
def patched_find_feature_limits(self, z_center: float, trans_name: str, z_window_kms: float = 150.0):
    """
    Patched version of _find_feature_limits that clamps indices to array bounds.
    """
    atom = ATOM_DATA.get(trans_name)
    if not atom: return None, None
    
    lambda_obs = atom['wave'] * (1 + z_center)
    c_kms = 299792.458
    
    dx_search = (z_window_kms / c_kms) * lambda_obs
    
    idx_center = np.searchsorted(self._x_ang, lambda_obs)
    # Safety: Ensure idx_center is within valid range
    idx_center = max(0, min(len(self._x_ang)-1, idx_center))

    if len(self._x_ang) > 10:
            sample_dx = np.mean(np.diff(self._x_ang[max(0, idx_center-10):min(len(self._x_ang), idx_center+10)]))
    else: sample_dx = 1.0
    
    px_width = int(dx_search / sample_dx) if sample_dx > 0 else 50
    
    start_idx = max(0, idx_center - px_width)
    end_idx = min(len(self._x_ang), idx_center + px_width)
    
    if end_idx - start_idx < 5: return None, None
    
    y_window = self._y_norm[start_idx:end_idx]
    dy_window = self._dy_norm[start_idx:end_idx]
    y_smooth = gaussian_filter1d(y_window, sigma=2.0)
    
    center_rel = idx_center - start_idx
    center_rel = max(0, min(len(y_window)-1, center_rel))

    # Scan Left
    left_boundary = 0
    lowest_point = y_smooth[center_rel]
    for i in range(center_rel, -1, -1):
        val = y_smooth[i]
        err = dy_window[i]
        if val < lowest_point: lowest_point = val
        is_local_max = (i > 0 and i < len(y_smooth)-1 and y_smooth[i] >= y_smooth[i-1] and y_smooth[i] >= y_smooth[i+1]) or (i == 0)
        if val > 0.9: left_boundary = i; break
        if is_local_max:
            if (val - lowest_point) > 4.0 * err: left_boundary = i; break
            
    # Scan Right
    right_boundary = len(y_window) - 1
    lowest_point = y_smooth[center_rel]
    for i in range(center_rel, len(y_window)):
        val = y_smooth[i]
        err = dy_window[i]
        if val < lowest_point: lowest_point = val
        is_local_max = (i > 0 and i < len(y_smooth)-1 and y_smooth[i] >= y_smooth[i-1] and y_smooth[i] >= y_smooth[i+1]) or (i == len(y_window)-1)
        if val > 0.9: right_boundary = i; break
        if is_local_max:
            if (val - lowest_point) > 4.0 * err: right_boundary = i; break

    # [FIX]: Clamp indices to valid array range
    idx_min = max(0, start_idx + left_boundary - 1)
    idx_max = min(len(self._x_ang) - 1, start_idx + right_boundary + 1) 
    
    if idx_max <= idx_min: return None, None
    return self._x_ang[idx_min], self._x_ang[idx_max]

# Apply the patch
VoigtFitterV2._find_feature_limits = patched_find_feature_limits

# --- HELPER FUNCTIONS ---
class HeadlessContext:
    def __init__(self): self._flags = []
    def _flags_cond(self, flag): return False
    def _flags_extr(self, flag): return None
    def _launch_recipe_dialog(self, *args): pass

def refine_z_seed(session, z_guess, window_kms):
    """Refines seed based on smoothed depth (robust peak finding)."""
    spec = session.spec
    x_ang = spec.x.to(au.Angstrom).value
    y = spec.y.value
    lam_rest = ATOM_DATA[TARGET_ION]['wave']
    lam_obs = lam_rest * (1 + z_guess)
    c_kms = 299792.458
    dv = c_kms * (x_ang - lam_obs) / lam_obs
    mask = (np.abs(dv) < window_kms)
    if np.sum(mask) < 5: return z_guess 
    
    depth = 1.0 - y[mask]
    depth_smoothed = gaussian_filter1d(depth, sigma=2.0)
    
    peak_idx = np.argmax(depth_smoothed)
    lam_peak = x_ang[mask][peak_idx]
    
    return (lam_peak / lam_rest) - 1.0

def merge_nearby_components(systs, velocity_threshold=50.0):
    """Merges components closer than velocity_threshold (km/s)."""
    from astrocook.core.system_list import SystemListV2
    components = systs.components
    if not components: return systs

    c_kms = 299792.458
    by_series = {}
    for c in components:
        if c.series not in by_series: by_series[c.series] = []
        by_series[c.series].append(c)
        
    new_components = []
    for series, comps in by_series.items():
        comps.sort(key=lambda x: x.z)
        merged = []
        if not comps: continue
        current = comps[0]
        
        for next_comp in comps[1:]:
            dz = next_comp.z - current.z
            dv = c_kms * dz / (1 + current.z)
            
            if dv < velocity_threshold:
                N1, N2 = 10**current.logN, 10**next_comp.logN
                N_tot = N1 + N2
                w1, w2 = N1/N_tot, N2/N_tot
                
                new_z = current.z * w1 + next_comp.z * w2
                new_b = current.b * w1 + next_comp.b * w2
                new_btur = current.btur * w1 + next_comp.btur * w2
                
                current = dataclasses.replace(current, z=new_z, logN=np.log10(N_tot), b=new_b, btur=new_btur)
            else:
                merged.append(current)
                current = next_comp
        merged.append(current)
        new_components.extend(merged)
            
    for i, c in enumerate(new_components):
        object.__setattr__(c, 'id', i+1)
            
    return SystemListV2(dataclasses.replace(systs._data, components=new_components))

def simulate_observation(session, snr, fwhm_kms, z_em_val):
    logging.info(f"Simulating Observation: SNR={snr}, FWHM={fwhm_kms} km/s")
    
    # 1. CONVOLUTION
    sigma_kms = fwhm_kms / 2.35482
    sess = session.spec.smooth(sigma_kms=sigma_kms)
    session = session.with_new_spectrum(sess)

    # 2. NOISE GENERATION
    sigma_noise = 1.0 / snr
    flux = session.spec.y.value
    noise_vector = np.random.normal(loc=0.0, scale=sigma_noise, size=len(flux))
    
    new_spec = session.spec.apply_expression('y', 'y + noise', {'noise': noise_vector})
    new_spec = new_spec.apply_expression('dy', f'{sigma_noise}')
    
    # 3. INJECT RESOLUTION METADATA
    c_kms = 299792.458
    resol_R = c_kms / fwhm_kms
    
    resol_array = np.full(len(flux), resol_R)
    col_resol = DataColumnV2(resol_array, au.dimensionless_unscaled, description="Resolving Power R")
    
    new_aux = new_spec._data.aux_cols.copy()
    new_aux['resol'] = col_resol
    new_spec_data = dataclasses.replace(new_spec._data, aux_cols=new_aux)
    
    # Update properties (z_em, resol)
    final_spec = new_spec.__class__(new_spec_data).with_properties(z_em=z_em_val, resol=resol_R)
    final_spec._data.meta['RESOL'] = resol_R
    
    return session.with_new_spectrum(final_spec)

def process_mock(file_path, output_dir, snr, fwhm_kms, do_merge):
    filename = os.path.basename(file_path)
    los_id = os.path.splitext(filename)[0]
    
    logging.info(f"--- Processing {los_id} ---")

    context = HeadlessContext()
    # Load raw file
    sess = SessionV2.open_new(file_path, "spec", context, "ascii_resvel_header")
    if sess is None: return None

    # --- DYNAMIC Z_SEED ---
    x_ang = sess.spec.x.to(au.Angstrom).value
    median_wave = np.nanmedian(x_ang)
    rest_wave = ATOM_DATA[TARGET_ION]['wave']
    z_seed_dynamic = (median_wave / rest_wave) - 1.0
    logging.info(f"Dynamic Z guess: {z_seed_dynamic:.4f} (Median wave: {median_wave:.1f} A)")

    # --- SIMULATION & FITTING ---
    sess = simulate_observation(sess, snr, fwhm_kms, z_em_val=z_seed_dynamic)

    z_best = refine_z_seed(sess, z_seed_dynamic, window_kms=600.0)
    logging.info(f"Seeding {TARGET_ION} at z={z_best:.5f}")
    sess = sess.absorbers.add_component(series=TARGET_ION, z=z_best, logN=13.5, b=10.0)
    seed_uuid = sess.systs.components[-1].uuid

    logging.info("Running AIC Optimization...")
    sess = sess.absorbers.optimize_system(
        uuid=seed_uuid,
        max_components=10,      
        threshold_sigma=2.0,    
        aic_penalty=1.0,        
        z_window_kms=2000.0,    
        min_dv=5.0,             
        patience=5              
    )

    if do_merge:
        logging.info("Merging components < 50 km/s...")
        new_systs = merge_nearby_components(sess.systs, velocity_threshold=50.0)
        sess = sess.with_new_system_list(new_systs)
    else:
        logging.info("Skipping merge (Raw Fit).")

    if not os.path.exists(output_dir): os.makedirs(output_dir)
    save_path = os.path.join(output_dir, f"{los_id}_snr{snr}.acs2")
    sess.save(save_path)
    
    csv_path = os.path.join(output_dir, f"{los_id}_snr{snr}_components.csv")
    sess.systs.t.write(csv_path, format='ascii.csv', overwrite=True)
    
    return save_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate and Fit CII Mocks")
    parser.add_argument('--input', type=str, required=True, help="Folder with raw mocks")
    parser.add_argument('--out', type=str, default="output_cii")
    parser.add_argument('--snr', type=float, default=20.0, help="Target SNR")
    parser.add_argument('--fwhm', type=float, default=30.0, help="Target FWHM (km/s)")
    parser.add_argument('--merge', action='store_true', help="Merge components < 50 km/s")
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

    files = glob.glob(os.path.join(args.input, "*.dat"))
    for f in files:
        try:
            saved = process_mock(f, args.out, args.snr, args.fwhm, args.merge)
            
            # --- INSPECT LOGIC ---
            if args.inspect and saved:
                abs_saved_path = os.path.abspath(saved)
                # Look for run_astrocook.py in the repository root
                launcher_path = os.path.join(os.path.dirname(os.path.dirname(astrocook.__file__)), "run_astrocook.py")
                
                if not os.path.exists(launcher_path):
                    logging.warning(f"Could not find run_astrocook.py at {launcher_path}. Check installation.")
                else:
                    logging.info(f"Launching Inspector: {launcher_path}")
                    subprocess.run([sys.executable, launcher_path, abs_saved_path])
                
        except Exception as e:
            logging.error(f"Error on {f}: {e}")