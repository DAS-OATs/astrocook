import astropy.units as au
from copy import deepcopy
import json
import logging
import numpy as np
from typing import TYPE_CHECKING, Optional

from astrocook.core.atomic_data import STANDARD_MULTIPLETS, ATOM_DATA
from astrocook.fitting.voigt_fitter import VoigtFitterV2
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
            {"name": "series", "type": str, "default": None, "doc": "Series"},
            # [CHANGE] Added resol parameter
            {"name": "resol", "type": float, "default": None, "doc": "Resolution (R or FWHM)"}
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
    "detect_anchor": {
        "brief": "Detect doublet anchor redshift.",
        "details": "Finds the redshift of maximum overlap between this session and a sibling session (doublet).",
        "params": [
            {"name": "sibling_name", "type": str, "default": "", "doc": "Name of the sibling session."},
            {"name": "trans_self", "type": str, "default": "OVI_1031", "doc": "Transition in this session."},
            {"name": "trans_sibling", "type": str, "default": "OVI_1037", "doc": "Transition in sibling session."}
        ],
        "url": "absorbers_cb.html#detect_anchor"
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
            {"name": "min_dv", "type": float, "default": 10.0, "doc": "Min separation (km/s). Lower for de-blending."},
            {"name": "group_depth", "type": int, "default": 2, "doc": "Grouping depth for finding neighbors."},
            {"name": "patience", "type": int, "default": 2, "doc": "Trials allowed without AIC improvement."}
        ],
        "url": "absorbers_cb.html#optimize_system"
    },
    "refit_all": {
        "brief": "Refit all systems.",
        "details": "Iteratively fits all disjoint groups of components (islands) in the spectrum.",
        "params": [
            {"name": "max_nfev", "type": int, "default": 200, "doc": "Max iterations per group"},
            {"name": "z_window_kms", "type": float, "default": 20.0, "doc": "Fit window around lines (km/s)"},
            {"name": "group_depth", "type": int, "default": 2, "doc": "Grouping depth (1=Neighbors, 2=FoF)"}
        ],
        "url": "absorbers_cb.html#refit_all"
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
            multiplet_list_arr = [m.strip() for m in multiplets.split(',') if m.strip()]
        except ValueError:
            logging.error(msg_param_fail); return 0
        
        try:
            logging.debug("identify_lines recipe: Calling spec.identify_lines...")
            new_spec_v2 = self._session.spec.identify_lines(
                multiplet_list=multiplet_list_arr,
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
                logging.warning("identify_lines: No candidates found.")
                return self._session.with_new_spectrum(new_spec_v2)

            if not auto_populate_b:
                return self._session.with_new_spectrum(new_spec_v2)

            # --- Auto-Populate Logic ---
            temp_session = self._session.with_new_spectrum(new_spec_v2)
            
            # Ensure SystemList exists
            if temp_session.systs is None:
                from astrocook.core.structures import SystemListDataV2
                from astrocook.core.system_list import SystemListV2
                temp_session = temp_session.with_new_system_list(SystemListV2(SystemListDataV2()))

            # Call the NEW populate recipe
            temp_recipe = RecipeAbsorbersV2(temp_session)
            final_session = temp_recipe.populate_from_identification(
                region_id_col='abs_ids', # Unused but kept for API compat
                series_map_json=series_map_json,
                z_map_json=z_map_json
            )
            
            return final_session
        except Exception as e:
            logging.error(f"Failed during identify_lines: {e}", exc_info=True)
            return 0
    
    def populate_from_identification(self, region_id_col: str = 'abs_ids',
                       series_map_json: str = None, z_map_json: str = None) -> 'SessionV2':
        """
        Populates the system list from identification maps.
        Uses the upgraded add_component to create linked multiplets and fits them immediately.
        """
        from astrocook.fitting.voigt_fitter import VoigtFitterV2
        from astrocook.core.atomic_data import STANDARD_MULTIPLETS

        if not series_map_json or not z_map_json:
            logging.error("populate_from_identification: Missing maps.")
            return 0
        
        try:
            series_map = json.loads(series_map_json)
            z_map = json.loads(z_map_json)
            # Ensure keys are integers
            series_map = {int(k): v for k, v in series_map.items()}
            z_map = {int(k): v for k, v in z_map.items()}
        except Exception as e:
            logging.error(f"populate_from_identification: Failed to deserialize maps: {e}")
            return 0

        running_systs = self._session.systs
        spec = self._session.spec

        # 1. Trigger Continuum if needed
        if spec.norm is None:
            logging.info("Flux not normalized. Triggering continuum estimate...")
            from astrocook.recipes.continuum import RecipeContinuumV2
            cont_sess = RecipeContinuumV2(self._session).estimate_auto()
            spec = cont_sess.spec 
        
        logging.info(f"Populating {len(series_map)} systems...")

        # 2. Iterate and Build
        for sys_id in sorted(series_map.keys()):
            series = series_map[sys_id]
            z_val = z_map.get(sys_id, 0.0)
            
            # A. Add Single Component (e.g. 'CIV')
            running_systs = running_systs.add_component(series, z_val, logN=13.5, b=15.0)
            
            # B. Immediate Local Fit
            # Even though it's one component, VoigtFitter will fit all lines in the doublet
            if running_systs.components:
                target_uuid = running_systs.components[-1].uuid
                
                try:
                    with running_systs.fitting_context([target_uuid], group_depth=0):
                        fitter = VoigtFitterV2(spec, running_systs)
                        new_systs, _, res = fitter.fit(max_nfev=100, z_window_kms=20.0, verbose=0)
                        
                        if res and res.success:
                            running_systs = new_systs
                except Exception as e:
                    logging.warning(f"Fit failed for {series}: {e}")

        # 3. Cleanup and Model Generation
        if running_systs.constraint_model:
            running_systs.constraint_model.set_active_components(None)

        final_fitter = VoigtFitterV2(spec, running_systs)
        _, model_flux = final_fitter.compute_model_flux()
        
        # 4. Return Result
        from astrocook.core.structures import DataColumnV2
        import dataclasses
        from copy import deepcopy
        
        new_aux_cols = deepcopy(spec._data.aux_cols)
        new_aux_cols['model'] = DataColumnV2(
            values=model_flux,
            unit=spec.y.unit,
            description="Voigt Fit Model"
        )
        new_spec_data = dataclasses.replace(spec._data, aux_cols=new_aux_cols)
        new_spec = spec.__class__(new_spec_data) 
        
        return self._session.with_new_spectrum(new_spec).with_new_system_list(running_systs)
        
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
                         b: str = 'None', series: str = 'None', resol: str = 'None') -> 'SessionV2':
        try:
            changes = {}
            if z != 'None': changes['z'] = float(z)
            if logN != 'None': changes['logN'] = float(logN)
            if b != 'None': changes['b'] = float(b)
            if series != 'None': changes['series'] = str(series)
            if resol != 'None': changes['resol'] = float(resol) # Handle resolution update
            
            # Pass dictionary unpacking to SystemListV2.update_component
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

    def detect_anchor(self, sibling_name: str, trans_self: str, trans_sibling: str) -> float:
        try:
            sibling_session = None
            if hasattr(self._session, '_gui'):
                for hist in self._session._gui.session_histories:
                    if hist.display_name == sibling_name:
                        sibling_session = hist.current_state
                        break
            
            if not sibling_session:
                logging.error(f"Sibling session '{sibling_name}' not found.")
                return 0.0
            
            # --- CALL CORE LOGIC ---
            z_found = self._session.spec.detect_doublet_z(
                sibling_session.spec, trans_self, trans_sibling
            )
            logging.info(f"Detected anchor at z={z_found:.5f}")
            return z_found
            
        except Exception as e:
            logging.error(f"detect_anchor failed: {e}")
            return 0.0
        
    # Let's implement the logic we wrote in batch_driver, tailored for the Recipe API.
    # Since Recipes usually operate on 'self._session', we might need to pass the second session explicitly.
    
    def detect_doublet_anchor(self, sibling_session: 'SessionV2', 
                              trans_primary: str, trans_secondary: str) -> float:
        from scipy.ndimage import gaussian_filter1d
        
        # 1. Get Data
        s1 = self._session
        s2 = sibling_session
        x1, y1 = s1.spec.x.value, s1.spec.y.value
        x2, y2 = s2.spec.x.value, s2.spec.y.value
        
        if len(x1) == 0 or len(x2) == 0: return 0.0

        # 2. Physics
        if trans_primary not in ATOM_DATA or trans_secondary not in ATOM_DATA:
            logging.error("Unknown transitions"); return 0.0

        lam1 = ATOM_DATA[trans_primary]['wave'] / 10.0
        lam2 = ATOM_DATA[trans_secondary]['wave'] / 10.0

        # 3. Correlation
        z1 = x1 / lam1 - 1.0
        z2 = x2 / lam2 - 1.0
        
        z_min, z_max = max(z1.min(), z2.min()), min(z1.max(), z2.max())
        if z_max <= z_min: return 0.0
        
        dz = np.median(np.diff(z1))
        z_grid = np.arange(z_min, z_max, dz)
        
        y1_int = np.interp(z_grid, z1, y1)
        y2_int = np.interp(z_grid, z2, y2)
        
        depth1 = np.maximum(0.0, 1.0 - gaussian_filter1d(y1_int, 2.0))
        depth2 = np.maximum(0.0, 1.0 - gaussian_filter1d(y2_int, 2.0))
        
        score = depth1 * depth2
        best_idx = np.argmax(score)
        
        logging.info(f"Anchor detected at z={z_grid[best_idx]:.5f}")
        return float(z_grid[best_idx])

    def fit_component(self, uuid: str, max_nfev: str = '2000', z_window_kms: str = '20.0', group_depth: str = '2') -> 'SessionV2':
        try:
            max_nfev_i = int(max_nfev)
            z_win_f = float(z_window_kms)
            depth_i = int(group_depth)
            
            if not self._session.systs: return 0

            current_session = self._session

            # [CHANGE] Trigger Continuum Auto-Estimation if needed
            if current_session.spec.norm is None:
                logging.info("Flux not normalized and no continuum. Triggering estimate_auto...")
                from astrocook.recipes.continuum import RecipeContinuumV2
                cont_rec = RecipeContinuumV2(current_session)
                # This returns a NEW session with continuum
                current_session = cont_rec.estimate_auto()
            
            # Use the possibly-updated session for fitting
            with current_session.systs.fitting_context([uuid], group_depth=depth_i):
                fitter = VoigtFitterV2(current_session.spec, current_session.systs)
                new_systs, model_flux, res = fitter.fit(max_nfev=max_nfev_i, z_window_kms=z_win_f)
            
            if new_systs.constraint_model:
                new_systs.constraint_model.set_active_components(None)

            from astrocook.core.structures import DataColumnV2
            from copy import deepcopy
            import dataclasses
            
            new_aux_cols = deepcopy(current_session.spec._data.aux_cols)
            new_aux_cols['model'] = DataColumnV2(
                values=model_flux,
                unit=current_session.spec.y.unit,
                description="Voigt Fit Model"
            )
            new_spec_data = dataclasses.replace(current_session.spec._data, aux_cols=new_aux_cols)
            new_spec = current_session.spec.__class__(new_spec_data) 
            
            session_with_spec = current_session.with_new_spectrum(new_spec)
            return session_with_spec.with_new_system_list(new_systs)
            
        except Exception as e:
            logging.error(f"Failed fit_component: {e}", exc_info=True); return 0
            
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

    def optimize_system(self, uuid: str, max_components: str = '5', 
                        threshold_sigma: str = '2.5', aic_penalty: str = '0.0',
                        z_window_kms: str = '100.0', min_dv: str = '10.0',
                        group_depth: str = '2', patience: str = '2') -> 'SessionV2':
        try:
            # 1. Parse Args
            max_c = int(max_components); thresh = float(threshold_sigma)
            aic_p = float(aic_penalty); z_win = float(z_window_kms)
            min_dv_f = float(min_dv)
            depth = int(group_depth); pat = int(patience)

            if not self._session.systs: return 0

            current_session = self._session

            # [CHANGE] Trigger Continuum Auto-Estimation if needed
            if current_session.spec.norm is None:
                logging.info("Flux not normalized and no continuum. Triggering estimate_auto...")
                from astrocook.recipes.continuum import RecipeContinuumV2
                cont_rec = RecipeContinuumV2(current_session)
                current_session = cont_rec.estimate_auto()

            # Pass the (possibly updated) session to optimize_hierarchy
            new_systs = current_session.systs.optimize_hierarchy(
                spec=current_session.spec,
                uuid_seed=uuid,
                max_components=max_c,
                threshold_sigma=thresh,
                aic_penalty=aic_p,
                z_window_kms=z_win,
                min_dv=min_dv_f,
                group_depth=depth,
                patience=pat
            )
            
            if new_systs.constraint_model:
                new_systs.constraint_model.set_active_components(None)

            from astrocook.fitting.voigt_fitter import VoigtFitterV2
            # Use current_session.spec (which might have new continuum)
            fitter = VoigtFitterV2(current_session.spec, new_systs)
            _, model_flux = fitter.compute_model_flux()
            
            from astrocook.core.structures import DataColumnV2
            from copy import deepcopy
            import dataclasses
            
            new_aux_cols = deepcopy(current_session.spec._data.aux_cols)
            new_aux_cols['model'] = DataColumnV2(
                values=model_flux,
                unit=current_session.spec.y.unit,
                description="Voigt Fit Model"
            )
            new_spec_data = dataclasses.replace(current_session.spec._data, aux_cols=new_aux_cols)
            new_spec = current_session.spec.__class__(new_spec_data) 
            
            return current_session.with_new_spectrum(new_spec).with_new_system_list(new_systs)
            
        except Exception as e:
             logging.error(f"optimize_system failed: {e}", exc_info=True)
             return 0
        

    def refit_all(self, max_nfev: str = '200', z_window_kms: str = '20.0', group_depth: str = '2') -> 'SessionV2':
        """
        Refits all components in the session by processing disjoint groups (islands) sequentially.
        """
        try:
            # 1. Parse Arguments
            max_nfev_i = int(max_nfev)
            z_win_f = float(z_window_kms)
            depth_i = int(group_depth)

            if not self._session.systs or not self._session.systs.components:
                logging.warning("refit_all: No components to fit.")
                return 0

            # 2. Trigger Continuum Auto-Estimate if needed (Crucial for Batch/Empty starts)
            current_session = self._session
            if current_session.spec.norm is None:
                logging.info("Flux not normalized. Triggering continuum estimate...")
                from astrocook.recipes.continuum import RecipeContinuumV2
                current_session = RecipeContinuumV2(current_session).estimate_auto()

            # 3. Identify Islands (Graph-based grouping)
            # We clone the list logic to find groups without modifying the object yet
            all_uuids = set(c.uuid for c in current_session.systs.components)
            processed_uuids = set()
            islands = []

            logging.info(f"Refit All: Scanning {len(all_uuids)} components for disjoint islands...")

            while all_uuids:
                seed = next(iter(all_uuids))
                # Use the current system list to find connections
                group = current_session.systs.get_connected_group([seed], max_depth=depth_i)
                
                # Verify group is not empty (sanity check)
                if not group: 
                    all_uuids.remove(seed)
                    continue

                islands.append(list(group))
                
                # Mark as processed
                processed_uuids.update(group)
                all_uuids.difference_update(group)

            logging.info(f"Refit All: Found {len(islands)} disjoint islands. Starting sequential fit...")

            # 4. Sequential Fit Loop
            # We update 'running_systs' iteratively. 
            running_systs = current_session.systs

            for i, island_uuids in enumerate(islands):
                logging.info(f"  -> Fitting Island {i+1}/{len(islands)} ({len(island_uuids)} components)...")
                
                # Create a Fitter for the CURRENT state of the system list
                fitter = VoigtFitterV2(current_session.spec, running_systs)
                
                # Apply the Mask for this island
                with running_systs.fitting_context(island_uuids, group_depth=0): # Depth 0 because we already grouped them!
                    # Fit
                    new_systs, _, res = fitter.fit(max_nfev=max_nfev_i, z_window_kms=z_win_f, verbose=0)
                    
                    # Update our running state
                    running_systs = new_systs

            # 5. Final Model Generation
            # Clean up any active masks
            if running_systs.constraint_model:
                running_systs.constraint_model.set_active_components(None)
            
            # Compute one final global model trace
            final_fitter = VoigtFitterV2(current_session.spec, running_systs)
            _, model_flux = final_fitter.compute_model_flux()

            # 6. Pack Result
            from astrocook.core.structures import DataColumnV2
            import dataclasses
            
            new_aux_cols = deepcopy(current_session.spec._data.aux_cols)
            new_aux_cols['model'] = DataColumnV2(
                values=model_flux,
                unit=current_session.spec.y.unit,
                description="Voigt Fit Model"
            )
            
            new_spec_data = dataclasses.replace(current_session.spec._data, aux_cols=new_aux_cols)
            new_spec = current_session.spec.__class__(new_spec_data) 
            
            logging.info("Refit All: Complete.")
            return current_session.with_new_spectrum(new_spec).with_new_system_list(running_systs)

        except Exception as e:
            logging.error(f"refit_all failed: {e}", exc_info=True)
            return 0