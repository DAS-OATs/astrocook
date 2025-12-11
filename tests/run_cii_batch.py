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

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from astrocook.core.session import SessionV2
from astrocook.core.structures import DataColumnV2
from astrocook.core.atomic_data import STANDARD_MULTIPLETS, ATOM_DATA
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.core.system_list import SystemListV2

# --- CONFIGURATION ---
TARGET_ION = 'CII_1334'  
DEFAULT_Z = 6.01         
RESOLUTION = 7000.0     

# --- PATCH 1: ATOMIC DATA ---
if 'CII' in STANDARD_MULTIPLETS:
    del STANDARD_MULTIPLETS['CII']

# --- PATCH 2: OPTIMIZE HIERARCHY ---
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

# --- HELPER FUNCTIONS ---
class HeadlessContext:
    def __init__(self): self._flags = []
    def _flags_cond(self, flag): return False
    def _flags_extr(self, flag): return None
    def _launch_recipe_dialog(self, *args): pass

def refine_z_seed(session, z_guess, window_kms):
    spec = session.spec
    x_ang = spec.x.to(au.Angstrom).value
    y = spec.y.value
    lam_rest = ATOM_DATA[TARGET_ION]['wave']
    lam_obs = lam_rest * (1 + z_guess)
    c_kms = 299792.458
    dv = c_kms * (x_ang - lam_obs) / lam_obs
    mask = (np.abs(dv) < window_kms)
    if np.sum(mask) < 5: return z_guess 
    
    y_smooth = gaussian_filter1d(y[mask], sigma=2.0)
    min_idx = np.argmin(y_smooth)
    lam_min = x_ang[mask][min_idx]
    
    return (lam_min / lam_rest) - 1.0

def process_mock(file_path, output_dir, z_seed=DEFAULT_Z, scale_err=1.0):
    filename = os.path.basename(file_path)
    los_id = os.path.splitext(filename)[0]
    logging.info(f"--- Processing {los_id} ---")

    context = HeadlessContext()
    sess = SessionV2.open_new(file_path, "spec", context, "ascii_resvel_header")
    if sess is None: return None

    # 1. SCALING ERRORS (Crucial Step)
    if scale_err != 1.0:
        logging.info(f"Rescaling errors by factor {scale_err}")
        # Apply expression returns a new SpectrumV2
        new_spec_scaled = sess.spec.apply_expression('dy', f'dy * {scale_err}')
        sess = sess.with_new_spectrum(new_spec_scaled)

    # 2. Inject Resolution
    logging.info(f"Injecting Resolution R={RESOLUTION}")
    new_aux = sess.spec._data.aux_cols.copy()
    new_aux['resol'] = DataColumnV2(np.full(len(sess.spec.x), RESOLUTION), au.dimensionless_unscaled)
    new_spec = sess.spec.__class__(dataclasses.replace(sess.spec._data, aux_cols=new_aux)).with_properties(z_em=z_seed, resol=RESOLUTION)
    new_spec._data.meta['RESOL'] = RESOLUTION
    sess = sess.with_new_spectrum(new_spec)

    z_best = refine_z_seed(sess, z_seed, window_kms=600.0)
    sess = sess.absorbers.add_component(series=TARGET_ION, z=z_best, logN=13.5, b=10.0)
    seed_uuid = sess.systs.components[-1].uuid

    logging.info("Running AIC Optimization...")
    # Restored to reasonable scientific defaults
    sess = sess.absorbers.optimize_system(
        uuid=seed_uuid,
        max_components=10,      
        threshold_sigma=1.0,    
        aic_penalty=0.0,        
        z_window_kms=1200.0,    
        min_dv=5.0,             
        patience=5              
    )

    logging.info("Merging components < 50 km/s...")
    # NOTE: merge_nearby_components must be defined in your script (as in previous versions)
    # I am re-including it here for completeness
    
    # Re-inserting function for standalone safety:
    def merge_nearby_components_local(systs, velocity_threshold=50.0):
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
        for i, c in enumerate(new_components): object.__setattr__(c, 'id', i+1)
        return SystemListV2(dataclasses.replace(systs._data, components=new_components))

    new_systs = merge_nearby_components_local(sess.systs, velocity_threshold=50.0)
    sess = sess.with_new_system_list(new_systs)

    if not os.path.exists(output_dir): os.makedirs(output_dir)
    save_path = os.path.join(output_dir, f"{los_id}.acs2")
    sess.save(save_path)
    csv_path = os.path.join(output_dir, f"{los_id}_components.csv")
    sess.systs.t.write(csv_path, format='ascii.csv', overwrite=True)
    
    return save_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True)
    parser.add_argument('--out', type=str, default="output_cii")
    parser.add_argument('--scale_err', type=float, default=1.0, help="Scale factor for error column")
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

    files = glob.glob(os.path.join(args.input, "*.dat"))
    for f in files:
        try:
            saved = process_mock(f, args.out, scale_err=args.scale_err)
            if args.inspect and saved:
                subprocess.run([sys.executable, "launch_pyside_app.py", saved])
        except Exception as e:
            logging.error(f"Error on {f}: {e}")