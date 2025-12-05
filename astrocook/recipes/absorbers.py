import astropy.units as au
import json
import logging
import numpy as np
from typing import TYPE_CHECKING, Optional

from astrocook.core.atomic_data import STANDARD_MULTIPLETS, ATOM_DATA
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import HistoryLogV2
from astrocook.gui.debug_utils import GLOBAL_PLOTTER

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
            {"name": "z_window_kms", "type": float, "default": 20.0, "doc": "Max velocity shift allowed (km/s)"},
            {"name": "group_depth", "type": int, "default": 2, "doc": "Grouping depth (1=Neighbors, 2=FoF)"} 
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
    "optimize_system": {
        "brief": "Optimize system (Iterative Fit).",
        "details": "Iteratively adds components to a system to minimize residuals using Akaike Information Criterion (AIC).",
        "params": [
            {"name": "uuid", "type": str, "default": "", "doc": "Target Component UUID (seed)", "gui_hidden": True},
            {"name": "max_components", "type": int, "default": 3, "doc": "Maximum number of extra components to try adding."},
            {"name": "threshold_sigma", "type": float, "default": 3.0, "doc": "Residual threshold (sigma) required to trigger a new component."},
            {"name": "aic_penalty", "type": float, "default": 5.0, "doc": "Minimum AIC improvement required to accept a new component."},
            {"name": "z_window_kms", "type": float, "default": 100.0, "doc": "Velocity window (km/s) around the system to analyze."},
            {"name": "group_depth", "type": int, "default": 2, "doc": "Grouping depth for finding neighbors."},
            {"name": "patience", "type": int, "default": 2, "doc": "Trials allowed without AIC improvement."}
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

    def fit_component(self, uuid: str, max_nfev: str = '2000', z_window_kms: str = '20.0', group_depth: str = '2') -> 'SessionV2':
        try:
            max_nfev_i = int(max_nfev)
            z_win_f = float(z_window_kms)
            depth_i = int(group_depth)
            
            if not self._session.systs: return 0
            
            # [CRITICAL] 1. Capture Original State
            # We must restore the input session's state later to avoid side effects in the GUI or subsequent scripts.
            original_active = self._session.systs.constraint_model._active_uuids
            
            try:
                # 2. Apply Transient Fit Mask
                # This restricts the fit to the specific group, but DOES NOT change the permanent 'is_free' flags in the data.
                self._session.systs.constraint_model.set_active_components([uuid], group_depth=depth_i)
                
                # 3. Run Fit
                fitter = VoigtFitterV2(self._session.spec, self._session.systs)
                new_systs, model_flux, res = fitter.fit(max_nfev=max_nfev_i, z_window_kms=z_win_f)
                
                if not res or not res.success:
                    logging.warning(f"Fit failed: {res.message if res else 'No result'}")

            finally:
                # [CRITICAL] 4. Restore Input Session State
                # Whether the fit succeeds or fails, we unmask the input session.
                self._session.systs.constraint_model.set_active_components(original_active)

            # 5. sanitize Output Session
            # The new system list is created fresh from data, so it should be clean.
            # But we explicitly force it to None to be 100% sure it doesn't carry over internal state.
            if new_systs.constraint_model:
                new_systs.constraint_model.set_active_components(None)

            # 6. Construct Result
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
    def optimize_system(self, uuid: str, max_components: str = '10', 
                        threshold_sigma: str = '2.5', aic_penalty: str = '0.0',
                        z_window_kms: str = '100.0', group_depth: str = '2',
                        patience: str = '2') -> 'SessionV2':
        """
        Full implementation of the iterative optimization logic.
        """
        try:
            max_comp_i = int(max_components)
            thresh_f = float(threshold_sigma)
            aic_pen_f = float(aic_penalty)
            z_win_f = float(z_window_kms)
            depth_i = int(group_depth)
            patience_i = int(patience)
        except ValueError:
            logging.error(msg_param_fail); return 0
            
        logging.info(f"Optimizing system (max {max_comp_i} comps, patience {patience_i})...")
        
        # 1. Initial Baseline Fit
        # fit_component now handles the cleanup, so 'current_session' is clean.
        current_session = self.fit_component(
            uuid, max_nfev='2000', z_window_kms=str(z_win_f), group_depth=str(depth_i)
        )
        if current_session == 0: return 0

        target_comp = current_session.systs.get_component_by_uuid(uuid)
        if not target_comp: return 0
        
        # Helper: Calculate AIC for the CURRENT state
        def get_aic(session):
            # We must apply the mask to calculate Chi2/AIC correctly for the group
            original_active = session.systs.constraint_model._active_uuids
            try:
                session.systs.constraint_model.set_active_components([uuid], group_depth=depth_i)
                fitter = VoigtFitterV2(session.spec, session.systs)
                
                # Quick fit (0 iters if just evaluating, or small # to refine)
                # Here we just want the Chi2 of the CURRENT parameters
                _, _, res = fitter.fit(max_nfev=10, z_window_kms=z_win_f, verbose=0)
                
                if not res: return float('inf'), None, 0.0, 0
                
                chi2 = np.sum(res.fun**2)
                k = len(res.x) # Number of free parameters
                aic = chi2 + 2 * k
                return aic, fitter, chi2, k
            finally:
                # Always clean up
                session.systs.constraint_model.set_active_components(original_active)

        best_aic, best_fitter, best_chi2, best_k = get_aic(current_session)
        
        if best_fitter is None:
            logging.warning("Baseline AIC check failed.")
            return current_session

        best_session = current_session
        logging.info(f"Baseline: AIC={best_aic:.1f}, Chi2={best_chi2:.1f}, k={best_k}")
        
        trials_without_improvement = 0
        
        for i in range(max_comp_i):
            # Recalculate AIC/Model for current session to find residuals
            current_aic, current_fitter, _, _ = get_aic(current_session)
            if current_fitter is None: break 

            x_full = current_session.spec.x.value
            y_full = current_session.spec.y.value
            dy_full = current_session.spec.dy.value
            
            # Compute full model to find residuals
            _, model_full_phys = current_fitter.compute_model_flux()
            
            sigma_resid = (y_full - model_full_phys) / dy_full
            sigma_resid[~np.isfinite(sigma_resid)] = 0.0
            
            # Define Search Window
            c_kms = 299792.458
            series_name = target_comp.series
            if series_name in ATOM_DATA:
                lam_target_ang = ATOM_DATA[series_name]['wave']
            elif series_name in STANDARD_MULTIPLETS:
                primary = STANDARD_MULTIPLETS[series_name][0]
                lam_target_ang = ATOM_DATA[primary]['wave']
            else:
                logging.warning(f"Unknown series {series_name}, stopping.")
                break
                
            lam_target_nm = lam_target_ang / 10.0
            lam_obs_nm = lam_target_nm * (1 + target_comp.z)
            
            dv = c_kms * (x_full - lam_obs_nm) / lam_obs_nm
            window_mask = np.abs(dv) < z_win_f
            
            # Mask residuals outside window
            masked_resid = sigma_resid.copy()
            masked_resid[~window_mask] = np.inf 
            
            # Find deepest dip
            min_idx = np.argmin(masked_resid)
            min_val = masked_resid[min_idx]
            
            logging.info(f"Iter {i+1}: Max Resid Sigma={min_val:.1f}")

            # Threshold Check
            if min_val > -thresh_f and trials_without_improvement == 0:
                logging.info("Residuals within threshold.")
                break

            # Add Candidate Component
            z_new = x_full[min_idx] / lam_target_nm - 1.0
            
            # add_component automatically runs fit_component on the new UUID
            # This returns a clean session with the new component optimized
            trial_session = current_session.absorbers.add_component(
                series=target_comp.series, z=z_new, logN=13.0, b=15.0, z_window_kms=str(z_win_f)
            )
            
            if trial_session == 0: break

            # Evaluate Trial
            trial_aic, _, t_chi2, t_k = get_aic(trial_session)
            logging.info(f"  -> Trial AIC: {trial_aic:.1f} (Chi2: {t_chi2:.1f}, k={t_k}) vs Best: {best_aic:.1f}")

            # Always update current to explore further
            current_session = trial_session

            # Decision
            if trial_aic < best_aic - aic_pen_f:
                logging.info("  -> New Best Found! Resetting patience.")
                best_aic = trial_aic
                best_session = trial_session
                trials_without_improvement = 0
            else:
                trials_without_improvement += 1
                logging.info(f"  -> No improvement. Patience: {trials_without_improvement}/{patience_i}")
                
            if trials_without_improvement >= patience_i:
                logging.info("Patience exhausted. Rolling back.")
                break
                
        return best_session