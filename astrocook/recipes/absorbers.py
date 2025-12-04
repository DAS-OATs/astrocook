import astropy.units as au
import json
import logging
import numpy as np
from typing import TYPE_CHECKING, Optional

from astrocook.core.atomic_data import STANDARD_MULTIPLETS, ATOM_DATA
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import HistoryLogV2
try:
    from astrocook.legacy.functions import trans_parse
except ImportError:
    # logging.error("V1 'trans_parse' not available. Identification recipe will fail.")
    trans_parse = lambda s: []
from astrocook.legacy.message import msg_param_fail

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2

ABSORBERS_RECIPES_SCHEMAS = {
    "identify_lines": {
        "brief": "Identify likely absorption systems.",
        "details": "Finds absorption regions, computes correlation signals for multiplets, and identifies the most likely candidates.",
        "params": [
            {"name": "multiplets", "type": str, "default": "CIV,SiIV,MgII,Ly_ab", "doc": "Multiplets to search for (comma-separated)."},
            {"name": "mask_col", "type": str, "default": "abs_mask", "doc": "Mask column (True=Unabsorbed)"},
            {"name": "min_pix_region", "type": int, "default": 3, "doc": "Minimum width in pixels to count as a region"},
            {"name": "merge_dv", "type": float, "default": 10.0, "doc": "Max velocity (km/s) to merge adjacent regions."},
            {"name": "score_threshold", "type": float, "default": 0.5, "doc": "Min R^2 score (0 to 1) to accept a doublet candidate."},
            {"name": "bypass_scoring", "type": bool, "default": False, "doc": "Accept all kinematic candidates, ignoring R^2 score."},
            {"name": "debug_rating", "type": bool, "default": False, "doc": "Show debug plots for R^2 ratings."},
            {"name": "auto_populate", "type": bool, "default": True, "doc": "Automatically add identified components to the system list?"}
        ],
        "url": "absorbers_cb.html#identify_lines"
    },
    "populate_from_identification": {
        "brief": "Populate components from identification.",
        "details": "Populates the system list using the regions and candidates found by 'identify_lines'.",
        "params": [
            {"name": "region_id_col", "type": str, "default": "abs_ids", "doc": "Column containing region IDs."},
            {"name": "series_map_json", "type": str, "default": "", "doc": "JSON mapping Region ID -> Series."},
            {"name": "z_map_json", "type": str, "default": "", "doc": "JSON mapping Region ID -> Redshift."}
        ],
        "url": "absorbers_cb.html#populate_from_identification",
        "gui_hidden": True
    },
    "add_component": {
        "brief": "Add a single component.",
        "details": "Adds a new manual Voigt component at a specific redshift.",
        "params": [
            {"name": "series", "type": str, "default": "Ly_a", "doc": "Transition or Multiplet name (e.g. 'CIV')"},
            {"name": "z", "type": float, "default": 0.0, "doc": "Redshift"},
            {"name": "logN", "type": float, "default": 13.5, "doc": "Column Density (log)"},
            {"name": "b", "type": float, "default": 10.0, "doc": "Doppler parameter (km/s)"},
            {"name": "btur", "type": float, "default": 0.0, "doc": "Turbulent broadening (km/s)"},
            {"name": "z_window_kms", "type": float, "default": 20.0, "doc": "Max shift for initial fit (km/s)"}
        ],
        "url": "absorbers_cb.html#add_component" 
    },
    "update_component": {
        "brief": "Update component.",
        "details": "Modify parameters of an existing component.",
        "params": [
            {"name": "uuid", "type": str, "default": "", "doc": "Component UUID", "gui_hidden": True},
            {"name": "z", "type": float, "default": None, "doc": "Redshift"},
            {"name": "logN", "type": float, "default": None, "doc": "logN"},
            {"name": "b", "type": float, "default": None, "doc": "b (km/s)"},
            {"name": "series", "type": str, "default": None, "doc": "Series"}
        ],
        "url": "absorbers_cb.html#update_component",
        "gui_hidden": True
    },
    "delete_component": {
        "brief": "Delete component.",
        "details": "Remove a component from the list.",
        "params": [
            {"name": "uuid", "type": str, "default": "", "doc": "Component UUID", "gui_hidden": True}
        ],
        "url": "absorbers_cb.html#delete_component",
        "gui_hidden": True
    },
    "fit_component": {
        "brief": "Fit component.",
        "details": "Fits the selected component (and its neighbors) using the Voigt engine.",
        "params": [
            {"name": "uuid", "type": str, "default": "", "doc": "Target Component UUID", "gui_hidden": True},
            {"name": "max_nfev", "type": int, "default": 200, "doc": "Max iterations"},
            {"name": "z_window_kms", "type": float, "default": 20.0, "doc": "Max velocity shift allowed (km/s)"}
        ],
        "url": "absorbers_cb.html#fit_component",
        "gui_hidden": True
    },
    "update_constraint": {
        "brief": "Update constraint.",
        "details": "Freeze, unfreeze, or link component parameters.",
        "params": [
            {"name": "uuid", "type": str, "default": "", "doc": "Component UUID", "gui_hidden": True},
            {"name": "param", "type": str, "default": "", "doc": "Parameter name (z, logN, b)"},
            {"name": "is_free", "type": bool, "default": None, "doc": "Is the parameter free to vary?"},
            {"name": "expression", "type": str, "default": None, "doc": "Math expression for linking"},
            {"name": "target_uuid", "type": str, "default": None, "doc": "UUID of the source component (for linking)"} 
        ],
        "url": "absorbers_cb.html#update_constraint",
        "gui_hidden": True
    },
    # --- NEW SCHEMA ---
    "optimize_system": {
        "brief": "Optimize system (Iterative Fit).",
        "details": "Iteratively adds components to a system to minimize residuals using Akaike Information Criterion (AIC).",
        "params": [
            {"name": "uuid", "type": str, "default": "", "doc": "Target Component UUID (seed)", "gui_hidden": True},
            {"name": "max_components", "type": int, "default": 3, "doc": "Maximum number of extra components to try adding."},
            {"name": "threshold_sigma", "type": float, "default": 3.0, "doc": "Residual threshold (sigma) required to trigger a new component."},
            {"name": "aic_penalty", "type": float, "default": 5.0, "doc": "Minimum AIC improvement required to accept a new component."},
            {"name": "z_window_kms", "type": float, "default": 100.0, "doc": "Velocity window (km/s) around the system to analyze."}
        ],
        "url": "absorbers_cb.html#optimize_system"
    }
}

class RecipeAbsorbersV2:
    """
    Recipes for identification, fitting, and management of absorption lines.
    Accessed via ``session.absorbers``.
    """
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'abs'

    def identify_lines(self, 
                       multiplets: str = "CIV,SiIV,MgII,Ly_ab",
                       mask_col: str = 'abs_mask', 
                       min_pix_region: str = '3',
                       merge_dv: str = '10.0',
                       score_threshold: str = '0.5',
                       bypass_scoring: str = 'False',
                       debug_rating: str = 'False',
                       auto_populate: str = 'True') -> 'SessionV2':
        try:
            min_pix_i = int(min_pix_region)
            merge_dv_f = float(merge_dv)
            score_thresh_f = float(score_threshold)
            bypass_scoring_b = bypass_scoring.lower() == 'true'
            debug_rating_b = debug_rating.lower() == 'true'
            auto_populate_b = str(auto_populate).lower() == 'true'
            multiplet_list = [m.strip() for m in multiplets.split(',') if m.strip()]
        except ValueError:
            logging.error(msg_param_fail); return 0
        
        try:
            logging.debug("identify_lines recipe: Calling spec.identify_lines...")
            new_spec_v2 = self._session.spec.identify_lines(
                multiplet_list=multiplet_list,
                mask_col=mask_col,
                min_pix_region=min_pix_i,
                merge_dv=merge_dv_f,
                score_threshold=score_thresh_f,
                bypass_scoring=bypass_scoring_b,
                debug_rating=debug_rating_b
            )
            
            series_map_json = new_spec_v2.meta.get('series_map_json')
            z_map_json = new_spec_v2.meta.get('z_map_json')

            if not series_map_json or not z_map_json:
                logging.warning("identify_lines: Identification ran but produced no component maps.")
                return self._session.with_new_spectrum(new_spec_v2)

            try:
                series_map = json.loads(series_map_json)
                z_map = json.loads(z_map_json)
                series_map = {int(k): v for k, v in series_map.items()}
                z_map = {int(k): v for k, v in z_map.items()}
            except Exception as e:
                logging.error(f"identify_lines: Failed to deserialize component maps: {e}")
                return self._session.with_new_spectrum(new_spec_v2)

            if not auto_populate_b:
                return self._session.with_new_spectrum(new_spec_v2)

            temp_session = self._session.with_new_spectrum(new_spec_v2)
            if temp_session.systs is None:
                from astrocook.core.structures import SystemListDataV2
                from astrocook.core.system_list import SystemListV2
                temp_session = temp_session.with_new_system_list(SystemListV2(SystemListDataV2()))

            temp_recipe = RecipeAbsorbersV2(temp_session)
            final_session = temp_recipe.populate_from_identification(
                region_id_col='abs_ids',
                series_map_json=series_map_json,
                z_map_json=z_map_json
            )
            
            return final_session
        except Exception as e:
            logging.error(f"Failed during identify_lines: {e}", exc_info=True)
            return 0
    
    def populate_from_identification(self, region_id_col: str = 'abs_ids',
                       series_map_json: str = None, z_map_json: str = None) -> 'SessionV2':
        if not series_map_json or not z_map_json:
            logging.error("populate_from_identification: Missing maps.")
            return 0
        try:
            series_map = json.loads(series_map_json)
            z_map = json.loads(z_map_json)
            series_map = {int(k): v for k, v in series_map.items()}
            z_map = {int(k): v for k, v in z_map.items()}
        except Exception as e:
            logging.error(f"populate_from_identification: Failed to deserialize maps: {e}")
            return 0

        try:
            new_systs_v2 = self._session.systs.add_components_from_regions(
                spec=self._session.spec,
                region_id_col=region_id_col,
                series_map=series_map,
                z_map=z_map
            )
            return self._session.with_new_system_list(new_systs_v2)
        except Exception as e:
            logging.error(f"Failed during populate_from_identification: {e}", exc_info=True)
            return 0
        
    def add_component(self, series: str = 'Ly_a', z: str = '0.0', 
                      logN: str = '13.5', b: str = '10.0', btur: str = '0.0',
                      z_window_kms: str = '20.0') -> 'SessionV2':
        try:
            z_f = float(z); logN_f = float(logN); b_f = float(b); btur_f = float(btur)
            z_win_f = float(z_window_kms)
        except ValueError: return 0
        
        try:
            if self._session.systs is None:
                from astrocook.core.structures import SystemListDataV2
                from astrocook.core.system_list import SystemListV2
                empty_systs = SystemListV2(SystemListDataV2())
                new_systs = empty_systs.add_component(series, z_f, logN_f, b_f, btur_f)
            else:
                new_systs = self._session.systs.add_component(series, z_f, logN_f, b_f, btur_f)
            
            added_comp = new_systs.components[-1]
            temp_session = self._session.with_new_system_list(new_systs)
            temp_recipe = RecipeAbsorbersV2(temp_session)
            
            logging.info(f"Auto-fitting new component {added_comp.series}...")
            final_session = temp_recipe.fit_component(added_comp.uuid, z_window_kms=str(z_win_f))
            return final_session
        except Exception as e:
            logging.error(f"Failed add_component: {e}"); return 0
        
    def update_component(self, uuid: str, z: str = 'None', logN: str = 'None', 
                         b: str = 'None', series: str = 'None') -> 'SessionV2':
        try:
            changes = {}
            if z != 'None': changes['z'] = float(z)
            if logN != 'None': changes['logN'] = float(logN)
            if b != 'None': changes['b'] = float(b)
            if series != 'None': changes['series'] = str(series)
            
            new_systs = self._session.systs.update_component(uuid, **changes)
            return self._session.with_new_system_list(new_systs)
        except Exception as e:
            logging.error(f"Failed update_component: {e}"); return 0

    def delete_component(self, uuid: str) -> 'SessionV2':
        try:
            new_systs = self._session.systs.delete_component(uuid)
            fitter = VoigtFitterV2(self._session.spec, new_systs)
            _, model_flux = fitter.compute_model_flux()
            
            from astrocook.core.structures import DataColumnV2
            from copy import deepcopy
            import dataclasses
            
            new_aux_cols = deepcopy(self._session.spec._data.aux_cols)
            new_aux_cols['model'] = DataColumnV2(
                values=model_flux,
                unit=self._session.spec.y.unit,
                description="Voigt Fit Model"
            )
            
            new_spec_data = dataclasses.replace(self._session.spec._data, aux_cols=new_aux_cols)
            new_spec = self._session.spec.__class__(new_spec_data) 
            session_with_spec = self._session.with_new_spectrum(new_spec)
            return session_with_spec.with_new_system_list(new_systs)

        except Exception as e:
            logging.error(f"Failed delete_component: {e}", exc_info=True)
            return 0

    def fit_component(self, uuid: str, max_nfev: str = '2000', z_window_kms: str = '20.0') -> 'SessionV2':
        try:
            max_nfev_i = int(max_nfev)
            z_win_f = float(z_window_kms)
            if not self._session.systs: return 0
            
            self._session.systs.constraint_model.set_active_components([uuid])
            fitter = VoigtFitterV2(self._session.spec, self._session.systs)
            
            new_systs, model_flux, res = fitter.fit(max_nfev=max_nfev_i, z_window_kms=z_win_f)
            
            if not res or not res.success:
                logging.warning(f"Fit failed: {res.message if res else 'No result'}")
            
            from astrocook.core.structures import DataColumnV2
            from copy import deepcopy
            import dataclasses
            
            new_aux_cols = deepcopy(self._session.spec._data.aux_cols)
            new_aux_cols['model'] = DataColumnV2(
                values=model_flux,
                unit=self._session.spec.y.unit,
                description="Voigt Fit Model"
            )
            new_spec_data = dataclasses.replace(self._session.spec._data, aux_cols=new_aux_cols)
            new_spec = self._session.spec.__class__(new_spec_data) 
            
            session_with_spec = self._session.with_new_spectrum(new_spec)
            return session_with_spec.with_new_system_list(new_systs)
            
        except Exception as e:
            logging.error(f"Failed fit_component: {e}", exc_info=True); return 0
        
    def update_constraint(self, uuid: str, param: str, is_free: bool = None, 
                          expression: str = None, target_uuid: str = None) -> 'SessionV2':
        try:
            new_systs = self._session.systs.update_constraint(
                uuid, param, is_free, expression, target_uuid
            )
            
            if expression and not is_free:
                import numpy as np
                class CompProxy:
                    def __init__(self, c): self.c = c
                    @property
                    def z(self): return self.c.z
                    @property
                    def logN(self): return self.c.logN
                    @property
                    def b(self): return self.c.b
                    @property
                    def btur(self): return self.c.btur

                comp_map = {c.uuid: CompProxy(c) for c in new_systs.components}
                eval_globals = {'p': comp_map, 'np': np, 'sqrt': np.sqrt, 'log': np.log10}
                
                try:
                    new_val = float(eval(expression, eval_globals))
                    logging.info(f"Link established. Auto-updating {param} to {new_val:.4f}")
                    new_systs = new_systs.update_component(uuid, **{param: new_val})
                except Exception as e:
                    logging.warning(f"Could not auto-update linked value: {e}")

            return self._session.with_new_system_list(new_systs)
            
        except Exception as e:
            logging.error(f"Failed update_constraint: {e}", exc_info=True)
            return 0

    # --- NEW METHOD ---
    def optimize_system(self, uuid: str, max_components: str = '3', 
                        threshold_sigma: str = '3.0', aic_penalty: str = '5.0',
                        z_window_kms: str = '100.0') -> 'SessionV2':
        """
        Iteratively improves a fit by adding components where residuals are high (Negative Flux).
        
        Logic:
        1. Performs an initial fit of the target system.
        2. Computes residuals (Flux - Model).
        3. Scans for the strongest absorption dip (negative residual) in the vicinity.
        4. If a dip exceeds threshold_sigma, tentatively adds a new component there.
        5. Re-fits and calculates the Akaike Information Criterion (AIC).
        6. If AIC improves (decreases) by at least aic_penalty, keeps the new component.
        7. Repeats until max_components is reached or no improvement is found.

        Parameters
        ----------
        uuid : str
            UUID of the seed component (defines the velocity window center and series).
        max_components : str (int)
            Maximum number of additional components to try adding.
        threshold_sigma : str (float)
            Significance level (e.g. 3.0 sigma) of the residual dip required to trigger an addition.
        aic_penalty : str (float)
            The AIC penalty term. A new model is accepted only if AIC_new < AIC_old - penalty.
        z_window_kms : str (float)
            The velocity window size (km/s) around the seed component to search for residuals.

        Returns
        -------
        SessionV2
            The optimized session (or the original if no improvement found).
        """
        try:
            max_comp_i = int(max_components)
            thresh_f = float(threshold_sigma)
            aic_pen_f = float(aic_penalty)
            z_win_f = float(z_window_kms)
        except ValueError:
            logging.error(msg_param_fail); return 0
            
        logging.info(f"Optimizing system (max {max_comp_i} extra comps)...")
        
        # Start with the current state as the best known state
        current_session = self._session
        
        # 1. Baseline Fit
        # Ensure we start with a fresh fit of the current model to get accurate Chi2/AIC baseline
        # We use the large window to ensure the chi2 calculation covers the whole complex
        current_session = current_session.absorbers.fit_component(uuid, z_window_kms=str(z_win_f))
        if current_session == 0: return 0

        # Retrieve the seed component to know series/redshift
        target_comp = current_session.systs.get_component_by_uuid(uuid)
        if not target_comp: 
            logging.error(f"Optimize: Target component {uuid} not found.")
            return 0
        
        # Calculate Baseline AIC
        # We need to peek inside the last fit result. 
        # fit_component doesn't return the 'res' object, only the session.
        # So we must manually re-run the fitter class here to get metrics.
        # This is cheap because the parameters are already optimized.
        
        # Optimization Loop
        for i in range(max_comp_i):
            # A. Metric Calculation (Manual Fitter Invocation)
            active_uuids = [uuid] # The Fitter expands this group automatically
            systs = current_session.systs
            systs.constraint_model.set_active_components(active_uuids)
            fitter = VoigtFitterV2(current_session.spec, systs)
            
            # Run a quick fit (should converge instantly) to get stats
            _, model_flux, res = fitter.fit(max_nfev=100, z_window_kms=z_win_f, verbose=0)
            
            if not res or not res.success:
                logging.warning("Optimization: Baseline fit check failed. Aborting.")
                break

            resid_vector = res.fun # Weighted residuals (y - model) / dy from fitter
            chi2 = np.sum(resid_vector**2)
            k = len(res.x) # Number of free parameters
            aic_current = chi2 + 2 * k
            
            # B. Find Residual Peak
            # We need physical residuals to find the dip location
            # Mask out pixels far from target to avoid fitting noise elsewhere
            x_full = current_session.spec.x.value
            y_full = current_session.spec.y.value
            dy_full = current_session.spec.dy.value
            
            # Recompute full model (physical)
            _, model_full_phys = fitter.compute_model_flux()
            
            # Sigma Residuals (y - model) / err
            # Absorption is y < model => residual is NEGATIVE
            sigma_resid = (y_full - model_full_phys) / dy_full
            sigma_resid[~np.isfinite(sigma_resid)] = 0.0
            
            # Create Window Mask
            c_kms = 299792.458
            lam_target = 0
            if target_comp.series in ATOM_DATA:
                lam_target = ATOM_DATA[target_comp.series]['wave'] * (1 + target_comp.z)
            elif target_comp.series in STANDARD_MULTIPLETS:
                # Use first line of multiplet as center reference
                primary = STANDARD_MULTIPLETS[target_comp.series][0]
                lam_target = ATOM_DATA[primary]['wave'] * (1 + target_comp.z)
            
            if lam_target == 0: break
                
            dv = c_kms * (x_full - lam_target) / lam_target
            window_mask = np.abs(dv) < z_win_f
            
            masked_resid = sigma_resid.copy()
            masked_resid[~window_mask] = 0.0
            
            # Find strongest absorption feature (minimum value)
            min_idx = np.argmin(masked_resid)
            min_val = masked_resid[min_idx]
            
            logging.info(f"Iter {i+1}: Current AIC={aic_current:.1f}. Max Resid Sigma={min_val:.1f} (at idx {min_idx})")
            
            # Threshold Check (Absorption dips are negative)
            if min_val > -thresh_f:
                logging.info("Residuals within threshold. Optimization complete.")
                break
                
            # C. Trial: Add Component
            z_new = x_full[min_idx] / (lam_target / (1 + target_comp.z)) - 1.0
            
            logging.info(f"Trial: Adding component at z={z_new:.5f} (series={target_comp.series})")
            
            # Create TRIAL session (branching)
            # Default guesses: logN=13.0 (weak), b=15.0
            trial_session = current_session.absorbers.add_component(
                series=target_comp.series, z=z_new, logN=13.0, b=15.0, z_window_kms=str(z_win_f)
            )
            
            if trial_session == 0: 
                logging.warning("Failed to add trial component.")
                break

            # D. Evaluate Trial
            # We must re-run the fitter manually on the trial session to get the exact AIC
            # (add_component runs a fit, but we need the raw numbers to compare apples-to-apples)
            trial_systs = trial_session.systs
            trial_systs.constraint_model.set_active_components([uuid]) # The group includes the new one now
            trial_fitter = VoigtFitterV2(trial_session.spec, trial_systs)
            
            _, _, res_trial = trial_fitter.fit(max_nfev=500, z_window_kms=z_win_f, verbose=0)
            
            if not res_trial or not res_trial.success:
                logging.warning("Trial fit failed/diverged.")
                break
                
            chi2_trial = np.sum(res_trial.fun**2)
            k_trial = len(res_trial.x)
            aic_trial = chi2_trial + 2 * k_trial
            
            logging.info(f"Trial AIC={aic_trial:.1f} vs Current {aic_current:.1f}")
            
            # E. Decision
            if aic_trial < aic_current - aic_pen_f:
                logging.info("Accepting new component.")
                # Update current_session to the trial result
                current_session = trial_session
                # Update the seed uuid? No, keep the original anchor.
            else:
                logging.info("Rejecting new component (insufficient improvement).")
                break
                
        return current_session