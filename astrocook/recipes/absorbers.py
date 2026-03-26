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
            {"name": "z_window_kms", "type": float, "default": 100.0, "doc": "Fit window (km/s)."}
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
            {"name": "btur", "type": float, "default": None, "doc": "btur (km/s)"},
            {"name": "series", "type": str, "default": None, "doc": "Series"},
            {"name": "resol", "type": float, "default": None, "doc": "Resolution (R or FWHM)"}
        ],
        "url": "absorbers_cb.html#update_component",
        "gui_hidden": True
    },
    "update_components": {
        "brief": "Batch update components.",
        "params": [
            {"name": "uuids", "type": list, "default": None, "gui_hidden": True},
            {"name": "z", "type": float, "default": None},
            {"name": "logN", "type": float, "default": None},
            {"name": "b", "type": float, "default": None},
            {"name": "btur", "type": float, "default": None},
            {"name": "resol", "type": float, "default": None}
        ]
    },
    "delete_component": {
        "brief": "Delete component(s).",
        "params": [
            {"name": "uuid", "type": str, "default": None, "gui_hidden": True},
            {"name": "uuids", "type": list, "default": None, "gui_hidden": True}
        ]
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
    "update_constraints_column": {
        "brief": "Bulk update constraints.",
        "details": "Freeze or unfreeze a specific parameter for ALL components.",
        "params": [
            {"name": "param", "type": str, "default": "", "doc": "Parameter name (z, logN, b, btur)"},
            {"name": "is_free", "type": bool, "default": True, "doc": "Set to True (free) or False (frozen)."}
        ],
        "url": "absorbers_cb.html#update_constraints_column",
        "gui_hidden": True
    },
    "add_linked_system": {
        "brief": "Add linked system.",
        "details": "Add multiple components (e.g. CIV, SiIV) and link their redshifts.",
        "params": [
            {"name": "series_list", "type": str, "default": "CIV,SiIV", "doc": "Comma-separated list of series."},
            {"name": "z", "type": float, "default": 0.0, "doc": "Redshift"},
            {"name": "logN", "type": float, "default": 13.5, "doc": "Column Density"},
            {"name": "b", "type": float, "default": 10.0, "doc": "Doppler parameter"},
            {"name": "z_window_kms", "type": float, "default": 20.0, "doc": "Fit window"}
        ],
        "url": "absorbers_cb.html#add_linked_system",
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
            {"name": "z_window_kms", "type": float, "default": 100.0, "doc": "Fit window around lines (km/s)"},
            {"name": "group_depth", "type": int, "default": 2, "doc": "Grouping depth (1=Neighbors, 2=FoF)"},
            {"name": "max_group_size", "type": int, "default": 20, "doc": "Max components per fit group (chunking)"},
            {"name": "region_limit", "type": str, "default": "None", "doc": "Limit fitting to 'red_side' or 'lya_forest'."},
            {"name": "selected_uuids", "type": list, "default": None, "gui_hidden": True}
        ],
        "url": "absorbers_cb.html#refit_all"
    },
    "clean_negligible": {
        "brief": "Remove unphysical/negligible components.",
        "details": "Removes components that are too weak, too broad, or unphysically narrow.",
        "params": [
            {"name": "min_logN", "type": float, "default": 11.5, "doc": "Minimum column density to keep."},
            {"name": "max_b", "type": float, "default": 80.0, "doc": "Maximum Doppler parameter (km/s)."},
            {"name": "min_b", "type": float, "default": 3.0, "doc": "Minimum Doppler parameter (km/s)."},
            {"name": "combined_check", "type": bool, "default": True, "doc": "If True, only removes broad lines if they are ALSO weak (b > max_b AND logN < 13.0)."}
        ],
        "url": "absorbers_cb.html#clean_negligible"
    },
}

class RecipeAbsorbersV2:
    """
    Recipes for identification, fitting, and management of absorption lines.

    These methods are accessed via the ``session.absorbers`` attribute and delegate
    logic to :class:`~astrocook.core.system_list.SystemListV2`,
    :class:`~astrocook.core.spectrum.SpectrumV2`, and the fitting engine.
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
        Auto-detect and fit components.

        The Master Pipeline: Identifies candidates (doublets or singles), populates
        them as systems, and refits the spectrum. Supports masking by region.
        Delegates to :meth:`identify_lines`, :meth:`populate_from_identification`,
        and :meth:`refit_all`.

        Parameters
        ----------
        multiplets : str, optional
            Comma-separated list of multiplets to search for. Defaults to ``"CIV,SiIV,MgII,Ly_ab"``.
        score_threshold : str, optional
            Min R^2 score (0-1) to accept a candidate. Defaults to ``"0.5"``.
        merge_dv : str, optional
            Max velocity (km/s) to merge adjacent regions. Defaults to ``"10.0"``.
        min_pix_region : str, optional
            Minimum width in pixels to count as a region. Defaults to ``"3"``.
        mask_col : str, optional
            Name of the mask column (True=Absorbed). Defaults to ``"abs_mask"``.
        z_window_kms : str, optional
            Velocity window for fitting (km/s). Defaults to ``"20.0"``.
        region_limit : str, optional
            Limit search to ``'red_side'``, ``'lya_forest'``, or ``'None'``.
            Defaults to ``"None"``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with detected systems,
            or 0 on failure.
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
        
        # --- 0.5 SMART CLEAN (Refactored) ---
        if region_limit != "None" and session_current.systs:
            try:
                spec = session_current.spec
                # Get bounds as Quantities
                obs_min, obs_max = spec.get_region_bounds(region_limit)
                
                # Delegate logic to SystemListV2
                new_systs = session_current.systs.remove_in_region(obs_min, obs_max)
                session_current = session_current.with_new_system_list(new_systs)
                    
            except Exception as e:
                logging.warning(f"Smart clean failed: {e}")
        
        # --- 1. PREPARE REGION MASK ---
        target_mask_col = mask_col
        if region_limit != "None":
            try:
                obs_min, obs_max = session_current.spec.get_region_bounds(region_limit)
                # Public API Call (Logic moved to Core)
                new_spec = session_current.spec.mask_region(mask_col, obs_min, obs_max)
                session_current = session_current.with_new_spectrum(new_spec)
            except ValueError as e:
                logging.warning(f"Region masking skipped: {e}")
        
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

        # --- 2.5 MERGE IDENTIFICATION IDs ---
        if old_abs_ids is not None and session_step1.spec.has_aux_column('abs_ids'):
            try:
                # 1. Get the newly generated IDs
                new_ids = session_step1.spec.get_column('abs_ids').value
                
                # 2. Shift new IDs to avoid collisions
                import numpy as np
                max_old = np.max(old_abs_ids) if len(old_abs_ids) > 0 else 0
                shifted_new = np.zeros_like(new_ids)
                mask_new = new_ids > 0
                shifted_new[mask_new] = new_ids[mask_new] + max_old
                
                # 3. Combine: Keep OLD, fill gaps with NEW
                final_ids = old_abs_ids.copy()
                fill_mask = (shifted_new > 0) & (final_ids == 0)
                final_ids[fill_mask] = shifted_new[fill_mask]
                
                # 4. Save
                final_spec = session_step1.spec.update_column('abs_ids', final_ids)
                session_step1 = session_step1.with_new_spectrum(final_spec)
                
            except Exception as e:
                logging.warning(f"Failed to merge abs_ids: {e}")

        spec_meta = session_step1.spec.meta
        series_map_json = spec_meta.get('series_map_json')
        z_map_json = spec_meta.get('z_map_json')
        
        if not series_map_json:
            logging.warning("components_auto: No candidates found.")
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
        Auto-fit metal doublets (Red Side).

        Wrapper for :meth:`components_auto` targeted at the red side (z < z_em).
        Finds CIV, SiIV, MgII.

        Parameters
        ----------
        multiplets : str, optional
            Multiplets to search for. Defaults to ``"CIV,SiIV,MgII"``.
        score_threshold : str, optional
            Min R^2 score (0-1). Defaults to ``"0.5"``.
        z_window_kms : str, optional
            Fit window (km/s). Defaults to ``"20.0"``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with detected doublets.
        """
        return self.components_auto(
            multiplets=multiplets,
            score_threshold=score_threshold,
            z_window_kms=z_window_kms,
            region_limit="red_side"  # <--- FORCE RED SIDE
        )

    def lya_auto(self, score_threshold: str = "0.1", min_b: str = "10.0",
                 z_window_kms: str = "100.0") -> 'SessionV2':
        """
        Auto-fit Ly-alpha forest.

        Wrapper for :meth:`components_auto` targeted at the forest. Includes
        filtering of interlopers based on the b-parameter.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.filter_by_criteria`.

        Parameters
        ----------
        score_threshold : str, optional
            Lower threshold for single lines (R^2). Defaults to ``"0.1"``.
        min_b : str, optional
            Minimum b parameter (km/s) to keep. Defaults to ``"10.0"``.
        z_window_kms : str, optional
            Fit window (km/s). Defaults to ``"100.0"``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with the forest model.
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
            
            # Use new Core method
            # Logic: Keep if (NOT Ly_a) OR (Ly_a AND b >= min_b)
            condition = lambda c: (c.series != 'Ly_a') or (c.b >= min_b_val)
            
            if final_session.systs:
                new_systs = final_session.systs.filter_by_criteria(condition)
                
                # Check if components were actually removed
                if len(new_systs.components) < len(final_session.systs.components):
                     logging.info("Refining fit after removing narrow lines...")
                     
                     # Update session with filtered list
                     final_session = final_session.with_new_system_list(new_systs)
                     
                     # Refit
                     rec_polish = RecipeAbsorbersV2(final_session)
                     final_session = rec_polish.refit_all(
                        max_nfev='200',
                        z_window_kms=z_window_kms,
                        group_depth='1',
                        max_group_size='25'
                     )
                else:
                     logging.info("No interlopers found. Clean.")

        except Exception as e:
            logging.warning(f"lya_auto cleanup failed: {e}", exc_info=True)

        return final_session

    def forest_auto(self, score_threshold: str = "0.1", min_b: str = "10.0",
                    z_window_kms: str = "20.0") -> 'SessionV2':
        """
        Auto-fit full Lyman series.

        Pipeline that models the entire Lyman series (Ly-alpha to Ly-limit)
        simultaneously. Uses 'Ly_a' detection but maps it to 'Lyman' series.

        Parameters
        ----------
        score_threshold : str, optional
            Threshold for detection. Defaults to ``"0.1"``.
        min_b : str, optional
            Minimum b parameter (km/s). Defaults to ``"10.0"``.
        z_window_kms : str, optional
            Fit window (km/s). Defaults to ``"20.0"``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with Lyman series systems.
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
        Check for associated metals.

        Scans existing systems and attempts to add other common metal ions
        (e.g., NV, SiIV) at the same redshift, keeping them if fit quality holds.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.add_component`.

        Parameters
        ----------
        ions : str, optional
            Comma-separated list of ions to check. Defaults to ``"CIV,SiIV,NV,OVI,CII,MgII"``.
        aic_penalty : str, optional
            AIC improvement required to keep a new ion (Not currently used in logic).
            Defaults to ``"5.0"``.
        z_window_kms : str, optional
            Fit window (km/s) for local island fitting. Defaults to ``"20.0"``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with added metals,
            or 0 on failure.
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
                       min_pix_region: str = '1',       
                       merge_dv: str = '2.0',           
                       score_threshold: str = '0.05',   
                       bypass_scoring: str = 'True',    
                       debug_rating: str = 'False') -> 'SessionV2':
        """
        Identify likely absorption systems.

        Finds absorption regions, computes correlation signals for multiplets, and
        identifies the most likely candidates using a scoring system.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.identify_lines`.

        Parameters
        ----------
        multiplets : str, optional
            Multiplets to search for (comma-separated). Defaults to ``"CIV,SiIV,MgII,Ly_ab"``.
        mask_col : str, optional
            Mask column (True=Unabsorbed, used to find absorptions). Defaults to ``"abs_mask"``.
        min_pix_region : str, optional
            Minimum width in pixels to count as a region. Defaults to ``"3"``.
        merge_dv : str, optional
            Max velocity (km/s) to merge adjacent regions. Defaults to ``"10.0"``.
        score_threshold : str, optional
            Min R^2 score (0 to 1) to accept a doublet candidate. Defaults to ``"0.5"``.
        bypass_scoring : str, optional
            If ``'True'``, accept all kinematic candidates ignoring score.
            Defaults to ``"False"``.
        debug_rating : str, optional
            If ``'True'``, show debug plots. Defaults to ``"False"``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with identification maps
            in metadata, or 0 on failure.
        """
        try:
            min_pix_i = int(min_pix_region)
            merge_dv_f = float(merge_dv)
            score_thresh_f = float(score_threshold)
            bypass_scoring_b = bypass_scoring.lower() == 'true'
            debug_rating_b = debug_rating.lower() == 'true'
            multiplet_list = [m.strip() for m in multiplets.split(',') if m.strip()]
            
            logging.info(f"Scanning with high sensitivity: dv={merge_dv_f}, min_pix={min_pix_i}")

            # Call the engine ONLY with the arguments it explicitly accepts
            new_spec_v2 = self._session.spec.identify_lines(
                multiplet_list=multiplet_list,
                mask_col=mask_col,
                min_pix_region=min_pix_i,
                merge_dv=merge_dv_f,
                score_threshold=score_thresh_f,
                bypass_scoring=bypass_scoring_b,
                debug_rating=debug_rating_b
            )
            
            if new_spec_v2 is None:
                logging.warning("identify_lines returned None. Returning current session.")
                return self._session
                
            return self._session.with_new_spectrum(new_spec_v2)

        except Exception as e:
            logging.error(f"Failed during identify_lines: {e}", exc_info=True)
            return self._session
    
    def populate_from_identification(self, region_id_col: str = 'abs_ids',
                       series_map_json: str = None, z_map_json: str = None) -> 'SessionV2':
        """
        Populate components from identification.

        Populates the system list using the regions and candidates found by
        :meth:`identify_lines`. Delegates to :meth:`astrocook.core.system_list.SystemListV2.add_component`
        and :meth:`astrocook.core.spectrum.SpectrumV2.update_model`.

        Parameters
        ----------
        region_id_col : str, optional
            Column containing region IDs. Defaults to ``"abs_ids"``.
        series_map_json : str, optional
            JSON string mapping Region ID -> Series Name.
        z_map_json : str, optional
            JSON string mapping Region ID -> Redshift.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with populated systems
            and updated model, or 0 on failure.
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

        # 3. Cleanup and Model Generation (Refactored)
        if running_systs.constraint_model:
            running_systs.constraint_model.set_active_components(None)

        final_fitter = VoigtFitterV2(spec, running_systs)
        _, model_flux = final_fitter.compute_model_flux()
        
        # Delegate to SpectrumV2
        new_spec = spec.update_model(model_flux)
        
        # Remove temporary JSON maps to prevent FITS header overflow
        # These keys are no longer needed once components are populated.
        if hasattr(new_spec, 'remove_meta_keys'):
            new_spec = new_spec.remove_meta_keys([
                'series_map_json', 
                'z_map_json', 
                'logN_map_json', 
                'b_map_json'
            ])

        return self._session.with_new_spectrum(new_spec).with_new_system_list(running_systs)
        
    def add_component(self, series: str = 'Ly_a', z: str = '0.0', 
                      logN: str = '13.5', b: str = '10.0', btur: str = '0.0',
                      z_window_kms: str = '20.0') -> 'SessionV2':
        """
        Add a single component.

        Adds a new manual Voigt component at a specific redshift and auto-fits it.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.add_component`.

        Parameters
        ----------
        series : str, optional
            Transition or Multiplet name (e.g. ``'CIV'``). Defaults to ``"Ly_a"``.
        z : str, optional
            Redshift. Defaults to ``"0.0"``.
        logN : str, optional
            Column Density (log cm^-2). Defaults to ``"13.5"``.
        b : str, optional
            Doppler parameter (km/s). Defaults to ``"10.0"``.
        btur : str, optional
            Turbulent broadening (km/s). Defaults to ``"0.0"``.
        z_window_kms : str, optional
            Max shift for initial fit (km/s). Defaults to ``"20.0"``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with the added component,
            or 0 on failure.
        """
        try:
            z_f = float(z); logN_f = float(logN); b_f = float(b); btur_f = float(btur)
            z_win_f = float(z_window_kms)
        except ValueError: return 0

        if ',' in series and series not in STANDARD_MULTIPLETS:
            # 1. Parse the string "FeII_2586,FeII_2600" -> ["FeII_2586", "FeII_2600"]
            members = [s.strip() for s in series.split(',') if s.strip()]
            
            # 2. Register globally so the Fitter can find it
            STANDARD_MULTIPLETS[series] = members
            logging.info(f"Registered custom multiplet: '{series}' -> {members}")
        
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
                         b: str = 'None', btur: str = 'None', series: str = 'None', resol: str = 'None') -> 'SessionV2':
        """
        Update component parameters.

        Modifies the physical parameters of an existing component.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.update_component`.

        Parameters
        ----------
        uuid : str
            Component UUID.
        z : str, optional
            Redshift. Defaults to ``'None'``.
        logN : str, optional
            Column Density (log). Defaults to ``'None'``.
        b : str, optional
            Doppler parameter (km/s). Defaults to ``'None'``.
        btur : str, optional
            Turbulent broadening (km/s). Defaults to ``'None'``.
        series : str, optional
            Series name. Defaults to ``'None'``.
        resol : str, optional
            Resolution (R or FWHM). Defaults to ``'None'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with updated component,
            or 0 on failure.
        """
        try:
            changes = {}
            if z != 'None': changes['z'] = float(z)
            if logN != 'None': changes['logN'] = float(logN)
            if b != 'None': changes['b'] = float(b)
            if btur != 'None': changes['btur'] = float(btur)
            if series != 'None': changes['series'] = str(series)
            if resol != 'None': changes['resol'] = float(resol) # Handle resolution update
            
            new_systs = self._session.systs.update_component(uuid, **changes)
            return self._session.with_new_system_list(new_systs)
        except Exception as e:
            logging.error(f"Failed update_component: {e}"); return 0
        
    def update_components(self, uuids: list, z: str = 'None', logN: str = 'None', 
                          b: str = 'None', btur: str = 'None',
                          series: str = 'None', resol: str = 'None') -> 'SessionV2':
        """
        Batch update component parameters.

        Modifies the physical parameters of multiple existing components simultaneously.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.update_components`.

        Parameters
        ----------
        uuids : list
            List of component UUIDs to update.
        z : str, optional
            Redshift. Defaults to ``'None'``.
        logN : str, optional
            Column Density (log). Defaults to ``'None'``.
        b : str, optional
            Doppler parameter (km/s). Defaults to ``'None'``.
        btur : str, optional
            Turbulent broadening (km/s). Defaults to ``'None'``.
        series : str, optional
            Series name. Defaults to ``'None'``.
        resol : str, optional
            Resolution (R or FWHM). Defaults to ``'None'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with updated components,
            or 0 on failure.
        """
        try:
            # Parse values
            changes = {}
            if z != 'None': changes['z'] = float(z)
            if logN != 'None': changes['logN'] = float(logN)
            if b != 'None': changes['b'] = float(b)
            if btur != 'None': changes['btur'] = float(btur)
            if series != 'None': changes['series'] = str(series)
            if resol != 'None': changes['resol'] = float(resol)
            
            if not changes: return self._session

            # Delegate to the new bulk method
            new_systs = self._session.systs.update_components(uuids, **changes)
            return self._session.with_new_system_list(new_systs)
            
        except Exception as e:
            logging.error(f"Failed update_components: {e}", exc_info=True)
            return 0

    def delete_component(self, uuid: str = None, uuids: list = None) -> 'SessionV2':
        """
        Delete component(s).

        Removes specified components and updates the model.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.delete_component`
        and :meth:`astrocook.core.spectrum.SpectrumV2.update_model`.

        Parameters
        ----------
        uuid : str, optional
            UUID of a single component to delete.
        uuids : list, optional
            List of UUIDs to delete.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with components removed.
        """
        from astrocook.fitting.voigt_fitter import VoigtFitterV2
        
        try:
            # 1. Delegate to SystemListV2 (which we updated to handle lists)
            # We pass both arguments; SystemListV2 handles the set logic.
            new_systs = self._session.systs.delete_component(uuid=uuid, uuids=uuids)
            
            # 2. Recompute Model (Only once!)
            fitter = VoigtFitterV2(self._session.spec, new_systs)
            _, model_flux = fitter.compute_model_flux()
            
            # 3. Update Spectrum
            new_spec = self._session.spec.update_model(model_flux)
            
            return self._session.with_new_spectrum(new_spec).with_new_system_list(new_systs)
            
        except Exception as e:
            logging.error(f"Failed delete_component: {e}", exc_info=True)
            return self._session

    def detect_anchor(self, sibling_name: str, trans_self: str, trans_sibling: str) -> float:
        """
        Detect doublet anchor redshift.

        Finds the redshift of maximum overlap between this session and a sibling
        session (useful for doublet identification across orders/exposures).
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.detect_doublet_z`.

        Parameters
        ----------
        sibling_name : str
            Name of the sibling session.
        trans_self : str
            Transition in this session (e.g., ``'OVI_1031'``).
        trans_sibling : str
            Transition in sibling session (e.g., ``'OVI_1037'``).

        Returns
        -------
        float
            The detected redshift, or 0.0 on failure.
        """
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
        """
        Fit component.

        Fits the selected component (and its connected neighbors) using the Voigt engine.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.fitting_context`
        and :class:`astrocook.fitting.voigt_fitter.VoigtFitterV2`.

        Parameters
        ----------
        uuid : str
            Target Component UUID.
        max_nfev : str, optional
            Max iterations for the optimizer. Defaults to ``"2000"``.
        z_window_kms : str, optional
            Max velocity shift allowed during fit (km/s). Defaults to ``"20.0"``.
        group_depth : str, optional
            Grouping depth (1=Neighbors, 2=FoF). Defaults to ``"2"``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with updated parameters
            and model, or 0 on failure.
        """
        try:
            max_nfev_i = int(max_nfev)
            z_win_f = float(z_window_kms)
            depth_i = int(group_depth)
            
            if not self._session.systs: return 0

            current_session = self._session

            # Define the callback using the recipe's stop check
            def check_stop(): return self._is_stop_requested()

            try:
                if current_session.spec.norm is None:
                    logging.info("Flux not normalized and no continuum. Triggering estimate_auto...")
                    from astrocook.recipes.continuum import RecipeContinuumV2
                    cont_rec = RecipeContinuumV2(current_session)
                    current_session = cont_rec.estimate_auto()

                with current_session.systs.fitting_context([uuid], group_depth=depth_i):
                    fitter = VoigtFitterV2(current_session.spec, current_session.systs, check_stop=check_stop)
                    new_systs, model_flux, res = fitter.fit(max_nfev=max_nfev_i, z_window_kms=z_win_f)
            
            except RuntimeError as e:
                if "stopped" in str(e):
                    logging.info("Fit stopped by user. Reverting to pre-fit state.")
                    return self._session # Return the session as it was before fitting started
                raise e # Re-raise if it's a real error
            
            if new_systs.constraint_model:
                new_systs.constraint_model.set_active_components(None)

            # Delegate to SpectrumV2
            new_spec = current_session.spec.update_model(model_flux)
            
            return current_session.with_new_spectrum(new_spec).with_new_system_list(new_systs)
            
        except Exception as e:
            logging.error(f"Failed fit_component: {e}", exc_info=True); return 0
        
    def update_constraint(self, uuid: str, param: str, is_free: bool = None, 
                          expression: str = None, target_uuid: str = None) -> 'SessionV2':
        """
        Update constraint.

        Freeze, unfreeze, or link component parameters.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.update_constraint`.

        Parameters
        ----------
        uuid : str
            Component UUID.
        param : str
            Parameter name (``z``, ``logN``, ``b``).
        is_free : bool, optional
            Is the parameter free to vary?
        expression : str, optional
            Math expression for linking (e.g., ``"p['uuid'].z"``).
        target_uuid : str, optional
            UUID of the source component (for linking tracking).

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with updated constraints,
            or 0 on failure.
        """
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

    def update_constraints_column(self, param: str, is_free: bool) -> 'SessionV2':
        """
        Bulk update constraints.

        Updates the constraint status (free/frozen) for a specific parameter
        across ALL components in the session.

        Parameters
        ----------
        param : str
            Parameter name (``z``, ``logN``, ``b``, ``btur``).
        is_free : bool
            Set to ``True`` (free) or ``False`` (frozen).

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with updated constraints,
            or 0 on failure.
        """
        if not self._session.systs: return 0
        
        current_systs = self._session.systs
        # Iterate over all components and apply the update_constraint logic
        # Since update_constraint returns a NEW system list, we chain the updates.
        
        # Optimization: We can modify the map directly if we access the underlying dict,
        # but let's use the API loop to be safe, modifying the wrapper in place.
        
        running_systs = current_systs
        
        for c in current_systs.components:
            # We explicitly pass None for expression/target to clear links if setting to free,
            # or preserve them if just freezing? 
            # Usually "Freeze All" implies breaking links or just fixing values.
            # The `update_constraint` logic handles "is_free=True/False" by updating the flag.
            
            # Note: If is_free=True, update_constraint clears links.
            # If is_free=False, it keeps links if they exist, or just freezes value.
            running_systs = running_systs.update_constraint(
                c.uuid, param, is_free=is_free
            )
            
        return self._session.with_new_system_list(running_systs)

    def add_linked_system(self, series_list: str, z: float, logN: float = 13.5, 
                          b: float = 10.0, z_window_kms: float = 20.0) -> 'SessionV2':
        """
        Add linked system.

        Add multiple components (e.g. CIV, SiIV) and link their redshifts to the
        first component in the list.

        Parameters
        ----------
        series_list : str
            Comma-separated list of series (e.g., ``"CIV,SiIV"``).
        z : float
            Redshift.
        logN : float, optional
            Column Density. Defaults to ``13.5``.
        b : float, optional
            Doppler parameter. Defaults to ``10.0``.
        z_window_kms : float, optional
            Fit window. Defaults to ``20.0``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with linked components
            fitted, or 0 on failure.
        """
        series_names = [s.strip() for s in series_list.split(',') if s.strip()]
        if not series_names: return 0
        
        current_session = self._session
        
        # 1. Add Primary (The Anchor)
        primary_series = series_names[0]
        logging.info(f"Adding Primary: {primary_series} at z={z:.5f}")
        
        # Use existing add_component recipe logic to handle creation
        # (This handles empty list creation etc.)
        # But we need the UUID, so we call the SystemList API directly
        
        if current_session.systs is None:
            from astrocook.core.structures import SystemListDataV2
            from astrocook.core.system_list import SystemListV2
            new_systs = SystemListV2(SystemListDataV2())
        else:
            new_systs = current_session.systs
            
        # Add Primary
        new_systs = new_systs.add_component(primary_series, z, logN, b)
        primary_uuid = new_systs.components[-1].uuid
        
        # 2. Add Secondaries (Linked)
        secondary_uuids = []
        
        for s_name in series_names[1:]:
            logging.info(f"Adding Linked: {s_name}")
            new_systs = new_systs.add_component(s_name, z, logN, b)
            sec_uuid = new_systs.components[-1].uuid
            secondary_uuids.append(sec_uuid)
            
            # Link Z to Primary
            expr = f"p['{primary_uuid}'].z"
            new_systs = new_systs.update_constraint(
                sec_uuid, 'z', is_free=False, expression=expr, target_uuid=primary_uuid
            )
            
        # 3. Create Temp Session & Fit
        # We fit the Primary (which pulls in the linked secondaries via grouping)
        temp_session = current_session.with_new_system_list(new_systs)
        temp_recipe = RecipeAbsorbersV2(temp_session)
        
        logging.info(f"Fitting linked system group...")
        # Fit centered on Primary
        final_session = temp_recipe.fit_component(primary_uuid, z_window_kms=str(z_window_kms))
        
        return final_session

    def optimize_system(self, uuid: str, max_components: str = '5', 
                        threshold_sigma: str = '2.0', aic_penalty: str = '0.0',
                        z_window_kms: str = '100.0', min_dv: str = '5.0',
                        group_depth: str = '2', patience: str = '2') -> 'SessionV2':
        """
        Optimize system by iteratively analyzing fit residuals.

        Scans the normalized residuals of the current fit and iteratively adds 
        new components to unexplained absorption features (flux dips). The process 
        stops when no residual exceeds the specified sigma threshold.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.optimize_hierarchy`.

        Parameters
        ----------
        uuid : str
            Target Component UUID (seed). Used to anchor the initial search.
        max_components : str, optional
            Maximum number of new components to add in this region. Defaults to ``"5"``.
        threshold_sigma : str, optional
            Residual threshold (sigma) required to trigger a new component addition. 
            Defaults to ``"2.0"``.
        aic_penalty : str, optional
            Min AIC improvement to accept new component. Defaults to ``"0.0"``.
        z_window_kms : str, optional
            Velocity window (km/s) to analyze. Defaults to ``"100.0"``.
        min_dv : str, optional
            Min separation (km/s). Defaults to ``"5.0"``.
        group_depth : str, optional
            Grouping depth. Defaults to ``"2"``.
        patience : str, optional
            Trials allowed without AIC improvement. Defaults to ``"2"``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with optimized systems,
            or 0 on failure.
        """
        try:
            max_c = int(max_components); thresh = float(threshold_sigma)
            z_win = float(z_window_kms); depth = int(group_depth)

            if not self._session.systs: return 0
            current_session = self._session

            # Ensure continuum exists to accurately calculate residuals
            if current_session.spec.norm is None:
                from astrocook.recipes.continuum import RecipeContinuumV2
                current_session = RecipeContinuumV2(current_session).estimate_auto()

            running_session = current_session
            
            # Iteratively search for unexplained absorption features
            for i in range(max_c):
                from astrocook.fitting.voigt_fitter import VoigtFitterV2
                fitter = VoigtFitterV2(running_session.spec, running_session.systs)
                _, model_flux = fitter.compute_model_flux()
                
                # Calculate normalized residuals: (Data - Model) / Error
                residuals = (running_session.spec.flux - model_flux) / running_session.spec.sig
                
                # Break the loop if no significant residual is found
                # Using 2.0 sigma as a balance between noise and weak blended lines
                if np.max(np.abs(residuals)) < thresh:
                    logging.info(f"Optimization: No residuals > {thresh} sigma. Stopping iteration.")
                    break
                
                # Identify the peak of the worst residual to use as a new seed
                worst_idx = np.argmax(np.abs(residuals))
                peak_w = running_session.spec.wave[worst_idx]
                logging.info(f"Optimization: Significant residual detected near wave={peak_w:.2f}. Fitting new component.")
                
                # Attempt to fit a new component at the location of the residual peak
                # min_dv=5.0 allows the logN=13 line to be placed close to larger neighbors
                new_systs = running_session.systs.optimize_hierarchy(
                    spec=running_session.spec,
                    uuid_seed=uuid, 
                    max_components=1,
                    threshold_sigma=thresh,
                    z_window_kms=z_win,
                    group_depth=depth,
                    min_dv=5.0
                )
                
                running_session = running_session.with_new_system_list(new_systs)

            # Recompute the final global model after all additions
            fitter = VoigtFitterV2(running_session.spec, running_session.systs)
            _, final_model = fitter.compute_model_flux()
            new_spec = running_session.spec.update_model(final_model)
            
            return running_session.with_new_spectrum(new_spec)
                        
        except Exception as e:
             logging.error(f"Failed optimize_system: {e}", exc_info=True)
             return 0
        

    def refit_all(self, max_nfev: str = '100', z_window_kms: str = '100.0', 
                  group_depth: str = '0', max_group_size: str = '20',
                  region_limit: str = "None",
                  selected_uuids: list = None) -> 'SessionV2':
        """
        Refit all systems.

        Iteratively fits all disjoint groups of components (islands) in the spectrum.
        Uses chunking for large groups.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.get_connected_group`
        and :class:`astrocook.fitting.voigt_fitter.VoigtFitterV2`.

        Parameters
        ----------
        max_nfev : str, optional
            Max iterations per group. Defaults to ``"100"``.
        z_window_kms : str, optional
            Fit window around lines (km/s). Defaults to ``"100.0"``.
        group_depth : str, optional
            Grouping depth (0=Iterative Single, 1=Neighbors). Defaults to ``"0"``.
        max_group_size : str, optional
            Max components per fit group (chunking). Defaults to ``"20"``.
        region_limit : str, optional
            Limit fitting to ``'red_side'`` or ``'lya_forest'``. Defaults to ``"None"``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with all systems refitted,
            or 0 on failure.
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

            # --- START NEW CONSOLIDATED LOGIC FOR BETA.4 ---
            # If selected_uuids is provided (e.g., from 'Refit Selected' in UI),
            # we fit the entire connected group once instead of looping.
            if selected_uuids:
                logging.info(f"Refitting {len(selected_uuids)} component{'s' if len(selected_uuids) > 1 else ''}...")
                
                # Use our new spatial-index-powered grouping logic
                combined_active_set = current_session.systs.get_connected_group(
                    selected_uuids, max_depth=depth_i
                )

                try:
                    # Execute a single fit pass for the entire 'Super-Group'
                    with current_session.systs.fitting_context(list(combined_active_set), group_depth=0):
                        fitter = VoigtFitterV2(current_session.spec, current_session.systs, 
                                              check_stop=self._is_stop_requested)
                        new_systs, model_flux, res = fitter.fit(
                            max_nfev=max_nfev_i, z_window_kms=z_win_f
                        )
                    
                    # Update model and return immediately
                    new_spec = current_session.spec.update_model(model_flux)
                    return current_session.with_new_spectrum(new_spec).with_new_system_list(new_systs)
                except Exception as e:
                    logging.error(f"Consolidated fit failed: {e}")
                    return 0
            # --- END NEW CONSOLIDATED LOGIC ---

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
            
            # Prepare selection set
            target_set = set(selected_uuids) if selected_uuids else None

            for c in current_session.systs.components:
                # Check selection first
                if target_set is not None and c.uuid not in target_set:
                    continue

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

            # Delegate to SpectrumV2
            new_spec = current_session.spec.update_model(model_flux)
            
            return current_session.with_new_spectrum(new_spec).with_new_system_list(final_systs)

        except Exception as e:
            logging.error(f"refit_all failed: {e}", exc_info=True)
            return 0
        
    def clean_negligible(self, min_logN: str = "11.5", max_b: str = "80.0", 
                         min_b: str = "3.0", combined_check: str = "True") -> 'SessionV2':
        """
        Clean negligible components using an adaptive Signal-to-Noise threshold.

        Removes unphysical or negligible components from the system list. Features an 
        adaptive threshold mechanism: if the median S/N of the spectrum is poor, the 
        column density threshold is relaxed to protect faint lines buried in noise.

        Parameters
        ----------
        min_logN : str, optional
            Minimum column density (log cm^-2) to keep. Defaults to ``"11.5"``.
        max_b : str, optional
            Maximum Doppler parameter (km/s). Defaults to ``"80.0"``.
        min_b : str, optional
            Minimum Doppler parameter (km/s). Defaults to ``"3.0"``.
        combined_check : str, optional
            If ``"True"``, broad components are only removed if they are also weak.
            Defaults to ``"True"``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with flagged components removed.
        """
        try:
            lim_N = float(min_logN)
            lim_b_max = float(max_b)
            lim_b_min = float(min_b)
            use_combined = (combined_check.lower() == "true")
            
            # Evaluate median Signal-to-Noise Ratio for adaptive thresholding
            snr_array = self._session.spec.flux / self._session.spec.sig
            median_snr = np.nanmedian(snr_array)
            
            # Relax the minimum logN requirement if the spectrum is highly noisy
            effective_logN = lim_N
            if median_snr < 10.0:
                effective_logN = lim_N - 0.3
                logging.info(f"Low median S/N detected ({median_snr:.1f}). Adjusting min_logN threshold to {effective_logN:.2f}")

        except ValueError:
            logging.error("Invalid parameters provided to clean_negligible."); return 0

        if not self._session.systs: return 0
        
        to_remove = []
        for c in self._session.systs.components:
            # 1. Flag unphysically narrow lines (noise spikes)
            if c.b < lim_b_min:
                logging.info(f"Flagged {c.series} (z={c.z:.4f}): Below minimum b-value (b={c.b:.1f})")
                to_remove.append(c.uuid)
                continue

            # 2. Flag weak lines (using the S/N adjusted threshold)
            if c.logN < effective_logN:
                logging.info(f"Flagged {c.series} (z={c.z:.4f}): Below minimum logN (logN={c.logN:.2f})")
                to_remove.append(c.uuid)
                continue

            # 3. Flag overly broad lines (continuum drifts / ghosts)
            if c.b > lim_b_max:
                if use_combined:
                    if c.logN < 13.0: 
                        logging.info(f"Flagged {c.series} (z={c.z:.4f}): Broad and weak (b={c.b:.1f}, logN={c.logN:.2f})")
                        to_remove.append(c.uuid)
                else:
                    logging.info(f"Flagged {c.series} (z={c.z:.4f}): Above maximum b-value (b={c.b:.1f})")
                    to_remove.append(c.uuid)

        if to_remove:
            logging.info(f"Cleaning: Removing {len(to_remove)} flagged outlier components...")
            return self.delete_component(uuids=to_remove)
        else:
            logging.info("Cleaning: No outliers found.")
            return self._session