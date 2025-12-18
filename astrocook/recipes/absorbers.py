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
    trans_parse = lambda s: []
from astrocook.legacy.message import msg_param_fail

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2

ABSORBERS_RECIPES_SCHEMAS = {
    # --- 1. THE MASTER PIPELINE ---
    "components_auto": {
        "brief": "Auto-detect and fit components.",
        "details": "General pipeline: Identifies candidates (doublets or singles), populates them as systems, and refits the spectrum.",
        "params": [
            {"name": "multiplets", "type": str, "default": "CIV,SiIV,MgII,Ly_ab", "doc": "Multiplets to search for."},
            {"name": "score_threshold", "type": float, "default": 0.5, "doc": "Min R^2 score (0-1) to accept a candidate."},
            {"name": "merge_dv", "type": float, "default": 10.0, "doc": "Max velocity (km/s) to merge adjacent regions."},
            {"name": "min_pix_region", "type": int, "default": 3, "doc": "Min pixels for a region."},
            {"name": "mask_col", "type": str, "default": "abs_mask", "doc": "Mask column."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "doc": "Window for fitting (km/s)."},
            # New generic region limiter
            {"name": "region_limit", "type": str, "default": "None", "doc": "Limit search to 'red_side', 'lya_forest', or 'None'."}
        ],
        "url": "absorbers_cb.html#components_auto"
    },
    # --- 2. WRAPPERS ---
    "doublets_auto": {
        "brief": "Auto-fit metal doublets (Red Side).",
        "details": "Wrapper for components_auto targeted at the red side (z < z_em). Finds CIV, SiIV, MgII.",
        "params": [
            {"name": "multiplets", "type": str, "default": "CIV,SiIV,MgII", "doc": "Multiplets to search for."},
            {"name": "score_threshold", "type": float, "default": 0.5, "doc": "Min score."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "doc": "Fit window (km/s)."}
        ],
        "url": "absorbers_cb.html#doublets_auto"
    },
    "lya_auto": {
        "brief": "Auto-fit Ly-alpha forest.",
        "details": "Wrapper for components_auto targeted at the forest. Includes interloper filtering.",
        "params": [
            {"name": "score_threshold", "type": float, "default": 0.1, "doc": "Lower threshold for single lines."},
            {"name": "min_b", "type": float, "default": 10.0, "doc": "Minimum b parameter (km/s) to keep."},
            {"name": "z_window_kms", "type": float, "default": 50.0, "doc": "Fit window (km/s)."}
        ],
        "url": "absorbers_cb.html#lya_auto"
    },
    "forest_auto": {
        "brief": "Auto-fit full Lyman series.",
        "details": "Pipeline that models the entire Lyman series (Ly-alpha to Ly-limit) simultaneously.",
        "params": [
            {"name": "score_threshold", "type": float, "default": 0.1, "doc": "Threshold for detection."},
            {"name": "min_b", "type": float, "default": 10.0, "doc": "Minimum b parameter (km/s)."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "doc": "Fit window (km/s)."}
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
            {"name": "max_group_size", "type": int, "default": 20, "doc": "Max components per fit group (chunking)"},
            {"name": "region_limit", "type": str, "default": "None", "doc": "Limit fitting to 'red_side' or 'lya_forest'."}
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

    def _is_stop_requested(self):
        """
        Robustly checks if a stop has been requested.
        Looks at the local session flag AND the active GUI worker.
        """
        # 1. Check direct flag (works if we are on the root session)
        if getattr(self._session, '_stop_flag', False):
            return True
            
        # 2. Check via GUI worker (works for derived sessions in pipelines)
        # The _gui reference is preserved across session copies.
        if hasattr(self._session, '_gui') and self._session._gui:
            gui = self._session._gui
            # Access the worker stored in MainWindow
            if hasattr(gui, '_current_worker') and gui._current_worker:
                # Check the worker's internal flag
                if getattr(gui._current_worker, '_is_stopped', False):
                    return True
                    
        return False

    def components_auto(self, multiplets: str = "CIV,SiIV,MgII,Ly_ab",
                      score_threshold: str = "0.5",
                      merge_dv: str = "10.0",
                      min_pix_region: str = "3",
                      mask_col: str = "abs_mask",
                      z_window_kms: str = "20.0",
                      region_limit: str = "None") -> 'SessionV2':
        """
        The Master Pipeline: Identifies, Populates, and Fits components.
        Supports masking by region ('red_side' or 'lya_forest').
        """
        logging.info(">> PROGRESS: 0")
        logging.info(f"Initializing component search ({region_limit})...")
        
        # --- 1. PREPARE REGION MASK ---
        # We create a temporary session with a restricted mask to force identify_lines
        # to ignore specific regions.
        
        session_current = self._session
        
        # --- 0. FORCE MASK REGENERATION ---
        # Always generate a fresh mask to ensure it matches current z_em/flux state.
        # This prevents "stale mask" bugs when running recipes sequentially.
        logging.info(f"Generating fresh '{mask_col}'...")
        try:
            from astrocook.recipes.continuum import RecipeContinuumV2
            rec_cont = RecipeContinuumV2(session_current)
            session_current = rec_cont.find_absorbed(
                smooth_len_lya="5000.0", 
                smooth_len_out="400.0",
                kappa="2.0",
                template="False"
            )
        except Exception as e:
            logging.error(f"Mask generation failed: {e}")
            return 0
        
        # --- 0.5 SMART CLEAN: Remove OLD components in this region ---
        # This ensures re-execution is idempotent (no duplicates).
        if region_limit != "None" and session_current.systs:
            try:
                spec = session_current.spec
                obs_min, obs_max = spec.get_region_bounds(region_limit)
                
                from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS
                
                # 1. Identify Survivors
                surviving_comps = []
                deleted_count = 0
                
                for c in session_current.systs.components:
                    # Resolve Wavelength
                    lam_rest = None
                    if c.series in ATOM_DATA:
                        lam_rest = ATOM_DATA[c.series]['wave']
                    elif c.series in STANDARD_MULTIPLETS:
                         # Use primary line for location check
                         prim = STANDARD_MULTIPLETS[c.series][0]
                         if prim in ATOM_DATA:
                             lam_rest = ATOM_DATA[prim]['wave']
                    
                    should_delete = False
                    if lam_rest:
                        lam_obs = lam_rest * (1 + c.z) * au.Angstrom
                        # Check overlap
                        if obs_min <= lam_obs <= obs_max:
                            should_delete = True
                    
                    if should_delete:
                        deleted_count += 1
                    else:
                        surviving_comps.append(c)
                
                # 2. Bulk Delete
                if deleted_count > 0:
                    logging.info(f"Clean Start: Removed {deleted_count} existing components in '{region_limit}'.")
                    
                    import dataclasses
                    from astrocook.core.system_list import SystemListV2
                    
                    # Create new system list with only survivors
                    new_data = dataclasses.replace(session_current.systs._data, components=surviving_comps)
                    new_systs = SystemListV2(new_data)
                    session_current = session_current.with_new_system_list(new_systs)
                    
            except Exception as e:
                logging.warning(f"Smart clean failed: {e}")
        
        # --- 1. PREPARE REGION MASK (Using SSOT) ---
        target_mask_col = mask_col
        if region_limit != "None":
            
            try:
                spec = session_current.spec
                # 1. Get Bounds (e.g. Ly-beta to Ly-alpha)
                obs_min, obs_max = spec.get_region_bounds(region_limit)
                logging.debug(f"Applying {region_limit} mask: {obs_min:.1f} - {obs_max:.1f}")

                # 2. Get Data Arrays
                # We work in Angstroms for the mask logic (internal convention)
                x_ang = spec.x.to(au.Angstrom).value
                min_ang = obs_min.to(au.Angstrom).value
                max_ang = obs_max.to(au.Angstrom).value
                
                # 3. Modify Mask
                # Copy original mask
                mask_vals = spec.get_column(mask_col).value.astype(bool).copy()
                
                # STRICT MASKING: False if OUTSIDE bounds
                mask_vals[x_ang < min_ang] = False
                mask_vals[x_ang > max_ang] = False
                
                # 4. Inject
                new_spec = session_current.spec.update_column(
                    target_mask_col, mask_vals, unit=au.dimensionless_unscaled
                )
                session_current = session_current.with_new_spectrum(new_spec)
                
            except ValueError as e:
                # Handle missing z_em or bad region name gracefully
                logging.warning(f"Region masking skipped: {e}")
                # Proceed without masking (fallback to full spectrum)
        
        # --- 2. IDENTIFY (using the masked session) ---
        logging.info("Step 1/3: Scanning spectrum for candidate pairs...")
        
        # Capture the OLD abs_ids before they are overwritten
        old_abs_ids = None
        if session_current.spec.has_aux_column('abs_ids'):
            old_abs_ids = session_current.spec.get_column('abs_ids').value.copy()
            
        # Create a temp recipe wrapper for the masked session
        rec_temp = RecipeAbsorbersV2(session_current)
        
        session_step1 = rec_temp.identify_lines(
            multiplets=multiplets, 
            mask_col=target_mask_col, # Use the restricted mask
            min_pix_region=min_pix_region,
            merge_dv=merge_dv, 
            score_threshold=score_threshold,
            bypass_scoring="False", 
            debug_rating="False"
        )
        
        logging.info(">> PROGRESS: 30")
        if self._is_stop_requested(): return self._session

        # --- 2.5 MERGE IDENTIFICATION IDs (Additive Visualization) ---
        if old_abs_ids is not None and session_step1.spec.has_aux_column('abs_ids'):
            try:
                new_ids = session_step1.spec.get_column('abs_ids').value
                # If new_ids has a value (detected region), keep it.
                # If it's 0 (empty), restore the old ID (from previous recipe).
                merged_ids = np.where(new_ids != 0, new_ids, old_abs_ids)
                
                # Inject back
                spec_merged = session_step1.spec.update_column(
                    'abs_ids', merged_ids, unit=au.dimensionless_unscaled
                )
                session_step1 = session_step1.with_new_spectrum(spec_merged)
                logging.debug("Merged new identifications with existing ones.")
            except Exception as e:
                logging.warning(f"Failed to merge abs_ids: {e}")
        
        # Check results
        spec_meta = session_step1.spec.meta
        series_map_json = spec_meta.get('series_map_json')
        z_map_json = spec_meta.get('z_map_json')
        
        if not series_map_json:
            logging.warning("components_auto: No candidates found in region.")
            # Return original session (undoing the temp mask creation)
            return self._session

        # --- 3. POPULATE ---
        logging.info("Step 2/3: Creating component models...")
        rec_step2 = RecipeAbsorbersV2(session_step1)
        session_step2 = rec_step2.populate_from_identification(
            region_id_col='abs_ids',
            series_map_json=series_map_json,
            z_map_json=z_map_json
        )

        logging.info(">> PROGRESS: 50")
        if self._is_stop_requested(): return self._session
        
        # --- 4. REFIT ALL ---
        logging.info("Step 3/3: Re-fitting all systems...")
        rec_step3 = RecipeAbsorbersV2(session_step2)
        final_session = rec_step3.refit_all(
            max_nfev='200',
            z_window_kms=z_window_kms,
            group_depth='1',
            max_group_size='25',
            region_limit=region_limit
        )
        
        # Handle the case where refit_all returns 0 (Empty/Clean)
        # Instead of crashing, we return the session from Step 2 (populated but empty/unfitted)
        if final_session == 0:
            logging.info("Refit skipped (no components to fit). Returning clean state.")
            return session_step2

        # Check for stop request on the valid final session

        if self._is_stop_requested(): return self._session
        
        logging.info("Component pipeline complete.")
        return final_session

    def doublets_auto(self, multiplets: str = "CIV,SiIV,MgII",
                      score_threshold: str = "0.5",
                      z_window_kms: str = "20.0") -> 'SessionV2':
        """
        Wrapper: Auto-detect metal doublets on the Red Side.
        """
        return self.components_auto(
            multiplets=multiplets,
            score_threshold=score_threshold,
            z_window_kms=z_window_kms,
            region_limit="red_side"  # <--- FORCE RED SIDE
        )

    def lya_auto(self, score_threshold: str = "0.1", min_b: str = "10.0",
                 z_window_kms: str = "50.0") -> 'SessionV2':
        """
        Wrapper: Auto-detect Ly-alpha in the Forest.
        """
        logging.info(">> PROGRESS: 0")
        logging.info("Initializing Lyman-alpha forest analysis...")
        
        # 1. Run Pipeline (FORCE region='lya_forest')
        final_session = self.components_auto(
            multiplets="Ly_a",
            score_threshold=score_threshold,
            z_window_kms=z_window_kms,
            mask_col='abs_mask',
            region_limit='lya_forest' # <--- FORCE FOREST
        )
        if self._is_stop_requested(): return self._session
        
        # 2. Cleanup Phase (Interlopers)
        logging.info(">> PROGRESS: 0")
        logging.info("Analysing line shapes to remove interlopers...")
        
        try:
            min_b_val = float(min_b)
            systs = final_session.systs
            
            valid_comps = []
            deleted_count = 0
            
            for c in systs.components:
                is_lya = (c.series == 'Ly_a')
                if is_lya and c.b < min_b_val:
                    deleted_count += 1
                    continue
                valid_comps.append(c)
            
            #if deleted_count > 0:
            #    logging.info(f"Rejected {deleted_count} lines (too narrow).")
            #    
            #    import dataclasses
            #    from astrocook.core.system_list import SystemListV2
            #    new_data = dataclasses.replace(systs._data, components=valid_comps)
            #    new_systs = SystemListV2(new_data)
            #    final_session = final_session.with_new_system_list(new_systs)
            #    
            #    # 3. Final Refit
            #    logging.info("Refining fit for confirmed systems...")
            #    rec_polish = RecipeAbsorbersV2(final_session)
            #    final_session = rec_polish.refit_all(
            #        max_nfev='200',
            #        z_window_kms=z_window_kms,
            #        group_depth='1',
            #        max_group_size='25'
            #    )
            #else:
            #    logging.info("No interlopers found. Clean.")
            #    logging.info(">> PROGRESS: 100")

        except Exception as e:
            logging.warning(f"lya_auto cleanup failed: {e}")

        return final_session

    def forest_auto(self, score_threshold: str = "0.1", min_b: str = "10.0",
                    z_window_kms: str = "20.0") -> 'SessionV2':
        """
        Pipeline for the full Lyman series. 
        """
        logging.info("Starting 'forest_auto' pipeline...")
        
        # 1. IDENTIFY 'Ly_a' lines
        # (We can use components_auto logic, but we need to intercept the map before populate)
        # So we manually call steps 1 & 2.
        
        # Force forest region
        session_id_only = self.identify_lines(
            multiplets="Ly_a", 
            score_threshold=score_threshold,
            # We can't easily pass region mask to identify_lines directly via public API without
            # duplicating the logic in components_auto. 
            # Ideally, identify_lines should support region_limit, but for now, 
            # let's trust the user or just run on whole spectrum (since forest is dominant).
            # OR better: Use components_auto but hacking the 'multiplets' is hard.
        )
        
        # 2. MODIFY MAPS: Swap 'Ly_a' -> 'Lyman'
        import json
        spec_meta = session_id_only.spec.meta
        series_json = spec_meta.get('series_map_json')
        
        if series_json:
            series_map = json.loads(series_json)
            for k, v in series_map.items():
                if v == 'Ly_a': series_map[k] = 'Lyman'
            new_series_json = json.dumps(series_map)
        else:
            new_series_json = None

        # 3. POPULATE & FIT
        rec_step2 = RecipeAbsorbersV2(session_id_only)
        session_step2 = rec_step2.populate_from_identification(
            series_map_json=new_series_json,
            z_map_json=spec_meta.get('z_map_json')
        )
        
        # 4. Refit All
        rec_step3 = RecipeAbsorbersV2(session_step2)
        return rec_step3.refit_all(z_window_kms=z_window_kms, group_depth='1')

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
        # Optimization: Group components by Z
        from collections import defaultdict
        systems_by_z = defaultdict(list)
        for c in current_session.systs.components:
            z_round = round(c.z, 4)
            systems_by_z[z_round].append(c)
            
        from astrocook.fitting.voigt_fitter import VoigtFitterV2
        
        running_systs = current_session.systs
        count_added = 0
        
        for z_val, comps in systems_by_z.items():
            existing_series = {c.series for c in comps}
            
            for ion in ion_list:
                if ion in existing_series: continue
                
                logging.debug(f"  Checking {ion} at z={z_val}...")
                
                # A. Add
                trial_systs = running_systs.add_component(ion, z_val, logN=13.0, b=10.0)
                new_uuid = trial_systs.components[-1].uuid
                
                # B. Fit (Local Island)
                group = trial_systs.get_connected_group([new_uuid])
                
                try:
                    with trial_systs.fitting_context(list(group), group_depth=0):
                        fitter = VoigtFitterV2(current_session.spec, trial_systs)
                        fitted_systs, _, res = fitter.fit(max_nfev=100, z_window_kms=z_win, verbose=0)
                        
                        # C. Evaluate
                        fitted_comp = fitted_systs.get_component_by_uuid(new_uuid)
                        if not fitted_comp: continue 
                        
                        keep = True
                        if fitted_comp.logN < 12.0: keep = False
                        if fitted_comp.b > 50.0: keep = False
                        if fitted_comp.b < 2.0: keep = False
                        
                        if keep:
                            logging.info(f"  -> Found {ion} at z={z_val} (logN={fitted_comp.logN:.2f})")
                            running_systs = fitted_systs
                            count_added += 1
                            existing_series.add(ion) 
                except Exception as e:
                    logging.debug(f"  Fit failed: {e}")

        if count_added > 0:
            logging.info(f"metals_auto: Added {count_added} new systems.")
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
        """
        from astrocook.fitting.voigt_fitter import VoigtFitterV2

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
            
            z_found = self._session.spec.detect_doublet_z(
                sibling_session.spec, trans_self, trans_sibling
            )
            logging.info(f"Detected anchor at z={z_found:.5f}")
            return z_found
            
        except Exception as e:
            logging.error(f"detect_anchor failed: {e}")
            return 0.0

    def fit_component(self, uuid: str, max_nfev: str = '2000', z_window_kms: str = '20.0', group_depth: str = '2') -> 'SessionV2':
        try:
            max_nfev_i = int(max_nfev)
            z_win_f = float(z_window_kms)
            depth_i = int(group_depth)
            
            if not self._session.systs: return 0

            current_session = self._session

            if current_session.spec.norm is None:
                logging.info("Flux not normalized and no continuum. Triggering estimate_auto...")
                from astrocook.recipes.continuum import RecipeContinuumV2
                cont_rec = RecipeContinuumV2(current_session)
                current_session = cont_rec.estimate_auto()
            
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
            max_c = int(max_components); thresh = float(threshold_sigma)
            aic_p = float(aic_penalty); z_win = float(z_window_kms)
            min_dv_f = float(min_dv)
            depth = int(group_depth); pat = int(patience)

            if not self._session.systs: return 0

            current_session = self._session

            if current_session.spec.norm is None:
                logging.info("Flux not normalized and no continuum. Triggering estimate_auto...")
                from astrocook.recipes.continuum import RecipeContinuumV2
                cont_rec = RecipeContinuumV2(current_session)
                current_session = cont_rec.estimate_auto()

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
                  group_depth: str = '0', max_group_size: str = '20',
                  region_limit: str = "None") -> 'SessionV2':
        """
        Refits components using Chunking.
        """
        try:
            max_nfev_i = int(max_nfev)
            z_win_f = float(z_window_kms)
            depth_i = int(group_depth)
            max_size_i = int(max_group_size)

            if not self._session.systs or not self._session.systs.components:
                logging.warning("refit_all: No components to fit.")
                return 0

            current_session = self._session
            if current_session.spec.norm is None:
                from astrocook.recipes.continuum import RecipeContinuumV2
                current_session = RecipeContinuumV2(current_session).estimate_auto()

            # --- 0. DETERMINE FILTER BOUNDS ---
            obs_min_limit = 0.0 * au.nm
            obs_max_limit = 1e6 * au.nm
            use_region_filter = False

            if region_limit != "None":
                try:
                    obs_min_limit, obs_max_limit = current_session.spec.get_region_bounds(region_limit)
                    use_region_filter = True
                    logging.info(f"Refit restricted to {region_limit} ({obs_min_limit:.1f} - {obs_max_limit:.1f}).")
                except ValueError:
                    logging.warning(f"Could not apply region limit '{region_limit}'. Fitting all.")

            # --- 1. IDENTIFY VALID COMPONENTS (Filtering Step) ---
            from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS
            
            uuids_in_region = set()
            
            for c in current_session.systs.components:
                is_valid = True
                if use_region_filter:
                    # Resolve primary rest wavelength
                    series = c.series
                    lam_rest = None
                    if series in ATOM_DATA:
                        lam_rest = ATOM_DATA[series]['wave']
                    elif series in STANDARD_MULTIPLETS:
                        prim = STANDARD_MULTIPLETS[series][0]
                        if prim in ATOM_DATA:
                            lam_rest = ATOM_DATA[prim]['wave']
                    
                    if lam_rest:
                        lam_obs_ang = lam_rest * (1 + c.z)
                        lam_obs_q = lam_obs_ang * au.Angstrom
                        if not (obs_min_limit <= lam_obs_q <= obs_max_limit):
                            is_valid = False
                
                if is_valid:
                    uuids_in_region.add(c.uuid)

            if not uuids_in_region:
                logging.info("No components found in the specified region to fit.")
                return current_session

            # --- 2. BUILD WORK UNITS (Grouping Step) ---
            work_units = []

            # 1. Define Work Units
            if depth_i == 0:
                logging.debug(f"Refit All: Strategy = Iterative Single-Component (Robust).")
                # Filter the list, then sort by Z
                valid_comps = [c for c in current_session.systs.components if c.uuid in uuids_in_region]
                sorted_comps = sorted(valid_comps, key=lambda c: c.z)
                work_units = [[c.uuid] for c in sorted_comps]
            else:
                logging.debug(f"Refit All: Strategy = Connected Groups (depth={depth_i}).")
                # Initialize 'todo' set with only valid UUIDs
                # (We only start groups from valid seeds, but neighbors outside region can still be pulled in)
                todo_uuids = uuids_in_region.copy()
                comp_map = {c.uuid: c for c in current_session.systs.components}

                while todo_uuids:
                    seed = next(iter(todo_uuids))
                    
                    # Find group (physics-based)
                    group = current_session.systs.get_connected_group([seed], max_depth=depth_i)
                    
                    if not group: 
                        todo_uuids.remove(seed); continue
                    
                    group_list = list(group)
                    
                    # Chunking Logic
                    if len(group_list) <= max_size_i:
                        work_units.append(group_list)
                    else:
                        sorted_group = sorted([comp_map[u] for u in group_list], key=lambda c: c.z)
                        sorted_uuids = [c.uuid for c in sorted_group]
                        
                        for k in range(0, len(sorted_uuids), max_size_i):
                            chunk = sorted_uuids[k : k + max_size_i]
                            work_units.append(chunk)
                        
                        logging.debug(f"  -> Split massive group ({len(group_list)} comps) into {len(sorted_uuids)//max_size_i + 1} chunks.")

                    # Remove found items from todo list
                    todo_uuids.difference_update(group)

            total_units = len(work_units)
            logging.info(f"Preparing to fit {total_units} absorption systems...")
            logging.info(">> PROGRESS: 0")

            # 2. Sequential Local Fit Loop
            reference_systs = current_session.systs
            updated_components_map = {} 

            import dataclasses
            from astrocook.core.system_list import SystemListV2

            for i, target_uuids in enumerate(work_units):

                if self._is_stop_requested():
                    logging.warning("Refit All: Stop requested. Aborting.")
                    return self._session  # Return original state (revert)
            
                pct = int(100 * (i + 1) / total_units)
                logging.info(f">> PROGRESS: {pct}")
                
                if target_uuids:
                    peek_comp = next((c for c in reference_systs.components if c.uuid == target_uuids[0]), None)
                    if peek_comp:
                        series_str = peek_comp.series
                        z_str = f"{peek_comp.z:.4f}"
                        logging.info(f"Fitting System {i+1}/{total_units}: {series_str} at z={z_str}")
                    else:
                        logging.info(f"Fitting Group {i+1}/{total_units}...")

                if not target_uuids: continue
                
                # A. Define Local Context Window
                target_comps = [c for c in reference_systs.components if c.uuid in target_uuids]
                if not target_comps: continue
                
                min_z = min(c.z for c in target_comps)
                max_z = max(c.z for c in target_comps)
                
                c_kms = 299792.458
                dz_buffer = (1 + max_z) * (1000.0 / c_kms)
                z_start = min_z - dz_buffer
                z_end = max_z + dz_buffer
                
                local_comps = [c for c in reference_systs.components if z_start <= c.z <= z_end]
                local_systs_data = dataclasses.replace(reference_systs._data, components=local_comps)
                local_systs = SystemListV2(local_systs_data)
                
                try:
                    fitter = VoigtFitterV2(current_session.spec, local_systs)
                    
                    with local_systs.fitting_context(target_uuids, group_depth=0):
                        new_local_systs, _, res = fitter.fit(max_nfev=max_nfev_i, z_window_kms=z_win_f, verbose=0)
                    
                    for fitted_comp in new_local_systs.components:
                        if fitted_comp.uuid in target_uuids:
                            updated_components_map[fitted_comp.uuid] = fitted_comp
                            
                except Exception as e:
                    logging.warning(f"Fit failed for unit {i}: {e}")

            # 3. Batch Merge
            logging.info(">> PROGRESS: 100")
            logging.debug(f"Refit All: Merging {len(updated_components_map)} fitted components...")
            
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