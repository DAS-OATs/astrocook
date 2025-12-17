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
    "doublets_auto": {
        "brief": "Auto-detect and fit doublets.",
        "details": "Automated pipeline: Identifies kinematic doublet candidates, populates them as systems, and refits the spectrum.",
        "params": [
            {"name": "multiplets", "type": str, "default": "CIV,SiIV,MgII,Ly_ab", "doc": "Multiplets to search for."},
            {"name": "score_threshold", "type": float, "default": 0.5, "doc": "Min R^2 score (0-1) to accept a candidate."},
            {"name": "merge_dv", "type": float, "default": 10.0, "doc": "Max velocity (km/s) to merge adjacent regions."},
            {"name": "min_pix_region", "type": int, "default": 3, "doc": "Min pixels for a region."},
            {"name": "mask_col", "type": str, "default": "abs_mask", "doc": "Mask column."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "doc": "Window for fitting (km/s)."},
        ],
        "url": "absorbers_cb.html#doublets_auto"
    },
    "lya_auto": {
        "brief": "Auto-fit Lyman-alpha forest.",
        "details": "Pipeline optimized for the forest: Detects Ly-a lines, handles saturation, and filters narrow interlopers.",
        "params": [
            {"name": "score_threshold", "type": float, "default": 0.1, "doc": "Lower threshold for single lines."},
            {"name": "min_b", "type": float, "default": 10.0, "doc": "Minimum b parameter (km/s) to keep."},
            {"name": "z_start", "type": float, "default": 0.0, "doc": "Start redshift (0=start of spectrum)."},
            {"name": "z_end", "type": float, "default": 10.0, "doc": "End redshift (10=end of spectrum)."},
            {"name": "z_window_kms", "type": float, "default": 100.0, "doc": "Fit window (km/s)."},
        ],
        "url": "absorbers_cb.html#lya_auto"
    },
    "forest_auto": {
        "brief": "Auto-fit full Lyman series.",
        "details": "Pipeline that models the entire Lyman series (Ly-alpha to Ly-limit) simultaneously.",
        "params": [
            {"name": "score_threshold", "type": float, "default": 0.1, "doc": "Threshold for detection."},
            {"name": "min_b", "type": float, "default": 10.0, "doc": "Minimum b parameter (km/s)."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "doc": "Fit window (km/s)."},
        ],
        "url": "absorbers_cb.html#forest_auto"
    },
    "metals_auto": {
        "brief": "Check for associated metals.",
        "details": "Scans existing systems and attempts to add other common metal ions (e.g. NV, SiIV) at the same redshift.",
        "params": [
            {"name": "ions", "type": str, "default": "CIV,SiIV,NV,OVI,CII,MgII", "doc": "Comma-separated list of ions to check."},
            {"name": "aic_penalty", "type": float, "default": 5.0, "doc": "AIC improvement required to keep a new ion."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "doc": "Fit window (km/s)."}
        ],
        "url": "absorbers_cb.html#metals_auto"
    },
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
            {"name": "group_depth", "type": int, "default": 2, "doc": "Grouping depth (1=Neighbors, 2=FoF)"},
            {"name": "max_group_size", "type": int, "default": 20, "doc": "Max components per fit group (chunking)"} # <--- NEW
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

    def doublets_auto(self, multiplets: str = "CIV,SiIV,MgII,Ly_ab",
                      score_threshold: str = "0.5",
                      merge_dv: str = "10.0",
                      min_pix_region: str = "3",
                      mask_col: str = "abs_mask",
                      z_window_kms: str = "20.0") -> 'SessionV2':
        """
        Orchestrates the identification, population, optimization, and fitting of doublets.
        Includes a 'Smart Skip' to avoid optimizing simple systems that fit well.
        """
        logging.info("Starting 'doublets_auto' pipeline...")
        
        # 1. IDENTIFY
        session_step1 = self.identify_lines(
            multiplets=multiplets, mask_col=mask_col, min_pix_region=min_pix_region,
            merge_dv=merge_dv, score_threshold=score_threshold,
            bypass_scoring="False", debug_rating="False"
        )
        
        spec_meta = session_step1.spec.meta
        series_map_json = spec_meta.get('series_map_json')
        z_map_json = spec_meta.get('z_map_json')
        
        if not series_map_json:
            logging.warning("doublets_auto: No candidates found.")
            return session_step1

        # 2. POPULATE (and Initial Fit)
        rec_step2 = RecipeAbsorbersV2(session_step1)
        session_step2 = rec_step2.populate_from_identification(
            region_id_col='abs_ids',
            series_map_json=series_map_json,
            z_map_json=z_map_json
        )

        # 4. REFIT ALL
        rec_step3 = RecipeAbsorbersV2(session_step2)
        final_session = rec_step3.refit_all(
            max_nfev='200',
            z_window_kms=z_window_kms,
            group_depth='2',
            max_group_size='10'
        )
        
        logging.info("doublets_auto pipeline complete.")
        return final_session
    
    def lya_auto(self, score_threshold: str = "0.1", min_b: str = "10.0",
                 z_start: str = "0.0", z_end: str = "10.0",
                 z_window_kms: str = "100.0") -> 'SessionV2':
        """
        Pipeline for the Lyman-alpha forest.
        """
        # Reuse doublets_auto logic but targeting 'Ly_a' specifically
        # and adding a post-processing cleanup step.
        
        logging.info("Starting 'lya_auto' pipeline...")
        
        # 1. Run Standard Pipeline with 'Ly_a' target
        # We assume the user wants to populate 'Ly_a' single lines.
        # (For Ly-series matching, use forest_auto)
        final_session = self.doublets_auto(
            multiplets="Ly_a", # Force Ly_a only
            score_threshold=score_threshold,
            z_window_kms=z_window_kms,
            mask_col='abs_mask' # Standard mask
        )
        
        # 2. Cleanup: Filter by b-parameter (Interlopers)
        # We need to access the SystemList directly
        try:
            min_b_val = float(min_b)
            systs = final_session.systs
            
            # Simple filter: Keep components with b >= min_b OR not Ly_a (don't delete metals!)
            valid_comps = []
            deleted_count = 0
            
            for c in systs.components:
                is_lya = (c.series == 'Ly_a')
                if is_lya and c.b < min_b_val:
                    deleted_count += 1
                    continue # Drop it
                valid_comps.append(c)
            
            if deleted_count > 0:
                logging.info(f"lya_auto: Removed {deleted_count} narrow interlopers (b < {min_b_val} km/s).")
                
                # Update System List
                import dataclasses
                from astrocook.core.system_list import SystemListV2
                new_data = dataclasses.replace(systs._data, components=valid_comps)
                new_systs = SystemListV2(new_data)
                
                # Refit Global Model (fastest way to update plot)
                final_session = final_session.with_new_system_list(new_systs)
                
                # We should probably refit once more to be safe, but let's just update model
                from astrocook.fitting.voigt_fitter import VoigtFitterV2
                fitter = VoigtFitterV2(final_session.spec, new_systs)
                _, model_flux = fitter.compute_model_flux()
                
                # Update Spectrum
                from astrocook.core.structures import DataColumnV2
                from copy import deepcopy
                new_aux = deepcopy(final_session.spec._data.aux_cols)
                new_aux['model'] = DataColumnV2(model_flux, final_session.spec.y.unit, "Voigt Fit Model")
                new_spec_d = dataclasses.replace(final_session.spec._data, aux_cols=new_aux)
                final_session = final_session.with_new_spectrum(final_session.spec.__class__(new_spec_d))

                # C. [CRITICAL FIX] Refit the survivors!
                # The remaining lines need to re-adjust now that their neighbors are gone.
                logging.info("lya_auto: Refitting survivors to fill gaps...")
                rec_polish = RecipeAbsorbersV2(final_session)
                final_session = rec_polish.refit_all(
                    max_nfev='200',
                    z_window_kms=z_window_kms,
                    group_depth='1',
                    max_group_size='25'
                )

        except Exception as e:
            logging.warning(f"lya_auto cleanup failed: {e}")

        return final_session

    def forest_auto(self, score_threshold: str = "0.1", min_b: str = "10.0",
                    z_window_kms: str = "20.0") -> 'SessionV2':
        """
        Pipeline for the full Lyman series. 
        Identifies 'Ly_a' candidates but populates them as 'Lyman' systems.
        """
        logging.info("Starting 'forest_auto' pipeline...")
        
        # 1. IDENTIFY 'Ly_a' lines first (we use Ly_a as the finder)
        session_step1 = self.identify_lines(
            multiplets="Ly_a", 
            score_threshold=score_threshold
        )
        
        # 2. MODIFY MAPS: Change 'Ly_a' to 'Lyman'
        # This is the trick: We found Ly_a, but we tell populate to create the full series.
        import json
        spec_meta = session_step1.spec.meta
        series_json = spec_meta.get('series_map_json')
        
        if series_json:
            series_map = json.loads(series_json)
            # Swap 'Ly_a' -> 'Lyman'
            for k, v in series_map.items():
                if v == 'Ly_a':
                    series_map[k] = 'Lyman'
            
            # Save back to meta (this is a bit hacky but efficient)
            spec_meta['series_map_json'] = json.dumps(series_map)
            # Note: SpectrumV2 is immutable, so we technically need a new one, 
            # but since we are just passing metadata to the next recipe step, 
            # we can cheat by passing the map directly to populate_from_identification 
            # or creating a temp session. Let's pass arguments explicitly.
            
            new_series_json = json.dumps(series_map)
        else:
            new_series_json = None

        # 3. POPULATE & FIT (using 'Lyman' series)
        rec_step2 = RecipeAbsorbersV2(session_step1)
        session_step2 = rec_step2.populate_from_identification(
            series_map_json=new_series_json,
            z_map_json=spec_meta.get('z_map_json')
        )
        
        # 4. Cleanup & Refit (Reuse lya_auto logic for filtering?)
        # For now, just refit all.
        rec_step3 = RecipeAbsorbersV2(session_step2)
        return rec_step3.refit_all(z_window_kms=z_window_kms)

    def metals_auto(self, ions: str = "CIV,SiIV,NV,OVI,CII,MgII", 
                    aic_penalty: str = "5.0", z_window_kms: str = "20.0") -> 'SessionV2':
        """
        Scans existing systems and tries to add other metal ions at the same redshift.
        """
        try:
            ion_list = [i.strip() for i in ions.split(',')]
            aic_pen = float(aic_penalty)
            z_win = float(z_window_kms)
        except ValueError: return 0
        
        if not self._session.systs or not self._session.systs.components:
            logging.warning("metals_auto: No existing systems to check.")
            return 0

        logging.info(f"metals_auto: Checking {len(ion_list)} ions for existing systems...")
        
        current_session = self._session
        # Get unique redshifts of existing systems (grouped by proximity)
        # Actually, simpler: Iterate over every component, try to add missing ions at that Z.
        
        # Optimization: Group components by Z to avoid checking the same Z 10 times
        from collections import defaultdict
        systems_by_z = defaultdict(list)
        for c in current_session.systs.components:
            z_round = round(c.z, 4)
            systems_by_z[z_round].append(c)
            
        # We need a Fitter for AIC calculation
        from astrocook.fitting.voigt_fitter import VoigtFitterV2
        
        # Base AIC
        fitter = VoigtFitterV2(current_session.spec, current_session.systs)
        # We need a way to calculate global AIC or local AIC. 
        # Since we are adding lines, global AIC is fine if we are consistent.
        # But global fit is slow. Let's do local checks.
        
        # Actually, let's use the 'optimize_hierarchy' logic but focused on SPECIES not COMPONENTS.
        # Since we don't have a specific 'add_species' core method yet, we will script it here.
        
        running_systs = current_session.systs
        
        count_added = 0
        
        for z_val, comps in systems_by_z.items():
            existing_series = {c.series for c in comps}
            
            # Check each target ion
            for ion in ion_list:
                # 1. Skip if already present
                if ion in existing_series: continue
                
                # 2. Skip if physically unlikely (e.g. Low Ion in High Ion system)? 
                # For now, check everything.
                
                # 3. Trial Add
                logging.debug(f"  Checking {ion} at z={z_val}...")
                
                # Use a lightweight trial system list
                # We need to copy running_systs to trial_systs
                # This is expensive. 
                # BETTER: Add it, fit locally, check AIC, remove if bad.
                
                # Save state (constraints are complex, so deepcopy logic needed or careful management)
                # Let's try: Add -> Fit(Island) -> Calc AIC -> Decide
                
                # A. Add
                trial_systs = running_systs.add_component(ion, z_val, logN=13.0, b=10.0)
                new_uuid = trial_systs.components[-1].uuid
                
                # B. Fit (Local Island)
                # We need to fit the NEW ion + the EXISTING system at this Z.
                # Find connected group for the new UUID
                group = trial_systs.get_connected_group([new_uuid])
                
                try:
                    # Calculate AIC BEFORE (on the group without the new line)
                    # This is tricky because the new line changes the model.
                    # Standard approach: Fit with new line. If LogN drops to 0 or errors are huge, reject.
                    # Or use strict AIC.
                    
                    with trial_systs.fitting_context(list(group), group_depth=0):
                        fitter = VoigtFitterV2(current_session.spec, trial_systs)
                        fitted_systs, _, res = fitter.fit(max_nfev=100, z_window_kms=z_win, verbose=0)
                        
                        # C. Evaluate
                        # Check significance: Is the fitted LogN significant? Is the line visible?
                        # Heuristic: If LogN < 11.5 or b > 100 or b < 1, it's probably noise.
                        
                        # Find the fitted component
                        fitted_comp = fitted_systs.get_component_by_uuid(new_uuid)
                        if not fitted_comp: continue # Should not happen
                        
                        keep = True
                        if fitted_comp.logN < 12.0: keep = False
                        if fitted_comp.b > 50.0: keep = False
                        if fitted_comp.b < 2.0: keep = False
                        
                        # Formal AIC check would require comparing Chi2 of 'running_systs' vs 'fitted_systs' 
                        # on the SAME pixel mask.
                        
                        if keep:
                            logging.info(f"  -> Found {ion} at z={z_val} (logN={fitted_comp.logN:.2f})")
                            running_systs = fitted_systs
                            count_added += 1
                            existing_series.add(ion) # Don't add twice
                except Exception as e:
                    logging.debug(f"  Fit failed: {e}")

        if count_added > 0:
            logging.info(f"metals_auto: Added {count_added} new systems.")
            # Final Polish
            return self.refit_all()
        else:
            logging.info("metals_auto: No new metals found.")
            return current_session

    def identify_lines(self, 
                       multiplets: str = "CIV,SiIV,MgII,Ly_ab",
                       mask_col: str = 'abs_mask', 
                       min_pix_region: str = '3',
                       merge_dv: str = '10.0',
                       score_threshold: str = '0.5',
                       bypass_scoring: str = 'False',
                       debug_rating: str = 'False') -> 'SessionV2':
        try:
            min_pix_i = int(min_pix_region)
            merge_dv_f = float(merge_dv)
            score_thresh_f = float(score_threshold)
            bypass_scoring_b = bypass_scoring.lower() == 'true'
            debug_rating_b = debug_rating.lower() == 'true'
            multiplet_list = [m.strip() for m in multiplets.split(',') if m.strip()]
        except ValueError:
            logging.error(msg_param_fail); return 0
        
        try:
            logging.debug("identify_lines: Calling spec.identify_lines...")
            # This relies on SpectrumV2.identify_lines saving SYSTEM names to JSON
            new_spec_v2 = self._session.spec.identify_lines(
                multiplet_list=multiplet_list,
                mask_col=mask_col,
                min_pix_region=min_pix_i,
                merge_dv=merge_dv_f,
                score_threshold=score_thresh_f,
                bypass_scoring=bypass_scoring_b,
                debug_rating=debug_rating_b
            )
            
            return self._session.with_new_spectrum(new_spec_v2)
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
            #if running_systs.components:
            #    target_uuid = running_systs.components[-1].uuid
            #    
            #    try:
            #        with running_systs.fitting_context([target_uuid], group_depth=0):
            #            fitter = VoigtFitterV2(spec, running_systs)
            #            new_systs, _, res = fitter.fit(max_nfev=100, z_window_kms=20.0, verbose=0)
            #            
            #            if res and res.success:
            #                running_systs = new_systs
            #    except Exception as e:
            #        logging.warning(f"Fit failed for {series}: {e}")

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
        

    def refit_all(self, max_nfev: str = '100', z_window_kms: str = '20.0', 
                  group_depth: str = '0', max_group_size: str = '20') -> 'SessionV2':
        """
        Refits components. 
        If group_depth=0 (Recommended), it fits components one-by-one (Iterative),
        which is robust against de-blending degeneracies.
        """
        try:
            max_nfev_i = int(max_nfev)
            z_win_f = float(z_window_kms)
            depth_i = int(group_depth)
            max_size_i = int(max_group_size)

            if not self._session.systs or not self._session.systs.components:
                logging.warning("refit_all: No components to fit.")
                return 0

            # Trigger Continuum if needed
            current_session = self._session
            if current_session.spec.norm is None:
                from astrocook.recipes.continuum import RecipeContinuumV2
                current_session = RecipeContinuumV2(current_session).estimate_auto()

            # 1. Define Work Units
            # If depth=0, we FORCE one-by-one fitting.
            if depth_i == 0:
                logging.info(f"Refit All: Strategy = Iterative Single-Component (Robust).")
                # Sort by redshift to sweep through the spectrum logically
                sorted_comps = sorted(current_session.systs.components, key=lambda c: c.z)
                work_units = [[c.uuid] for c in sorted_comps]
            else:
                # Standard Friends-of-Friends (Only use if you trust the groups are non-degenerate)
                logging.info(f"Refit All: Strategy = Connected Groups (depth={depth_i}).")
                all_uuids = set(c.uuid for c in current_session.systs.components)
                work_units = []
                
                comp_map = {c.uuid: c for c in current_session.systs.components}

                while all_uuids:
                    seed = next(iter(all_uuids))
                    group = current_session.systs.get_connected_group([seed], max_depth=depth_i)
                    
                    if not group: 
                        all_uuids.remove(seed); continue
                    
                    group_list = list(group)
                    
                    # If group is acceptable, add it
                    if len(group_list) <= max_size_i:
                        work_units.append(group_list)
                    else:
                        # GROUP IS TOO BIG: SPLIT IT
                        # 1. Sort by redshift
                        sorted_group = sorted([comp_map[u] for u in group_list], key=lambda c: c.z)
                        sorted_uuids = [c.uuid for c in sorted_group]
                        
                        # 2. Chunk it
                        # We overlap chunks slightly? No, standard partitioning is safer for stability.
                        # The "Context Window" in fit loop handles the boundaries (fixing neighbors).
                        for k in range(0, len(sorted_uuids), max_size_i):
                            chunk = sorted_uuids[k : k + max_size_i]
                            work_units.append(chunk)
                        
                        logging.info(f"  -> Split massive group ({len(group_list)} comps) into {len(sorted_uuids)//max_size_i + 1} chunks.")

                    all_uuids.difference_update(group)

            logging.info(f"Refit All: Processing {len(work_units)} units (Batch Mode)...")

            # 2. Sequential Local Fit Loop
            reference_systs = current_session.systs
            updated_components_map = {} 

            import dataclasses
            from astrocook.core.system_list import SystemListV2

            for i, target_uuids in enumerate(work_units):
                if not target_uuids: continue
                
                # A. Define Local Context Window
                # We look up the targets in the REFERENCE list
                target_comps = [c for c in reference_systs.components if c.uuid in target_uuids]
                if not target_comps: continue
                
                min_z = min(c.z for c in target_comps)
                max_z = max(c.z for c in target_comps)
                
                # Context Buffer: +/- 1000 km/s is enough for wings
                c_kms = 299792.458
                dz_buffer = (1 + max_z) * (1000.0 / c_kms)
                z_start = min_z - dz_buffer
                z_end = max_z + dz_buffer
                
                # B. Build Local System List
                # This includes neighbors, but they will be FIXED (parameters not in target_uuids)
                local_comps = [c for c in reference_systs.components if z_start <= c.z <= z_end]
                local_systs_data = dataclasses.replace(reference_systs._data, components=local_comps)
                local_systs = SystemListV2(local_systs_data)
                
                # C. Fit Locally
                try:
                    fitter = VoigtFitterV2(current_session.spec, local_systs)
                    
                    # group_depth=0 ensures ONLY 'target_uuids' are free. 
                    # Neighbors in 'local_comps' act as fixed background.
                    with local_systs.fitting_context(target_uuids, group_depth=0):
                        new_local_systs, _, res = fitter.fit(max_nfev=max_nfev_i, z_window_kms=z_win_f, verbose=0)
                    
                    # D. Collect Results
                    for fitted_comp in new_local_systs.components:
                        if fitted_comp.uuid in target_uuids:
                            updated_components_map[fitted_comp.uuid] = fitted_comp
                            
                except Exception as e:
                    logging.warning(f"Fit failed for unit {i}: {e}")

            # 3. Batch Merge
            logging.info(f"Refit All: Merging {len(updated_components_map)} fitted components...")
            
            final_comps = []
            for c in reference_systs.components:
                if c.uuid in updated_components_map:
                    final_comps.append(updated_components_map[c.uuid])
                else:
                    final_comps.append(c)
            
            new_systs_data = dataclasses.replace(reference_systs._data, components=final_comps)
            final_systs = SystemListV2(new_systs_data)

            # 4. Final Global Model
            if final_systs.constraint_model:
                final_systs.constraint_model.set_active_components(None)
            
            final_fitter = VoigtFitterV2(current_session.spec, final_systs)
            _, model_flux = final_fitter.compute_model_flux()

            from astrocook.core.structures import DataColumnV2
            from copy import deepcopy
            
            new_aux_cols = deepcopy(current_session.spec._data.aux_cols)
            new_aux_cols['model'] = DataColumnV2(
                values=model_flux,
                unit=current_session.spec.y.unit,
                description="Voigt Fit Model"
            )
            
            new_spec_data = dataclasses.replace(current_session.spec._data, aux_cols=new_aux_cols)
            new_spec = current_session.spec.__class__(new_spec_data) 
            
            logging.info("Refit All: Complete.")
            return current_session.with_new_spectrum(new_spec).with_new_system_list(final_systs)

        except Exception as e:
            logging.error(f"refit_all failed: {e}", exc_info=True)
            return 0

