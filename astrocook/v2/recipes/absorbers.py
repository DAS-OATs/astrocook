import astropy.units as au
import json
import logging
from typing import TYPE_CHECKING, Optional

from ...v1.vars import xem_d
from ..atomic_data import STANDARD_MULTIPLETS
from ..spectrum import SpectrumV2
from ..structures import HistoryLogV2
try:
    from ...v1.functions import trans_parse
except ImportError:
    logging.error("V1 'trans_parse' not available. Identification recipe will fail.")
    trans_parse = lambda s: []
from ...v1.message import msg_param_fail

if TYPE_CHECKING:
    from ..session import SessionV2

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
            {"name": "debug_rating", "type": bool, "default": False, "doc": "Show debug plots for R^2 ratings."}
        ],
        "url": "absorbers_cb.html#identify_lines"
    },
    "add_component": {
        "brief": "Add a single component.",
        "details": "Adds a new manual Voigt component at a specific redshift.",
        "params": [
            {"name": "series", "type": str, "default": "Ly_a", "doc": "Transition or Multiplet name (e.g. 'CIV')"},
            {"name": "z", "type": float, "default": 0.0, "doc": "Redshift"},
            {"name": "logN", "type": float, "default": 13.5, "doc": "Column Density (log)"},
            {"name": "b", "type": float, "default": 10.0, "doc": "Doppler parameter (km/s)"},
            {"name": "btur", "type": float, "default": 0.0, "doc": "Turbulent broadening (km/s)"}
        ],
        # We hide this from the menu because it's usually triggered by mouse, 
        # but keeping it in the schema allows scripting/logging.
        "url": "absorbers_cb.html#add_component" 
    },
}

class RecipeAbsorbersV2:
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
                       debug_rating: str = 'False') -> 'SessionV2':
        """
        API: Orchestrator recipe to identify absorption lines.
        (Thin wrapper for SpectrumV2.identify_lines)
        """
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
            # 1. Run the identification, which returns a new spectrum
            #    with the identification maps in its metadata.
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
            
            # 2. Get the maps from the new spectrum's metadata
            series_map_json = new_spec_v2.meta.get('series_map_json')
            z_map_json = new_spec_v2.meta.get('z_map_json')

            if not series_map_json or not z_map_json:
                logging.warning("identify_lines: Identification ran but produced no component maps. Only spectrum was updated.")
                return self._session.with_new_spectrum(new_spec_v2)

            logging.debug("identify_lines recipe: Calling systs.add_components_from_regions...")
            
            # 3. Deserialize maps (they are JSON strings)
            try:
                series_map = json.loads(series_map_json)
                z_map = json.loads(z_map_json)
                # Convert string keys back to int
                series_map = {int(k): v for k, v in series_map.items()}
                z_map = {int(k): v for k, v in z_map.items()}
            except Exception as e:
                logging.error(f"identify_lines: Failed to deserialize component maps: {e}")
                # Return the spectrum-only update
                return self._session.with_new_spectrum(new_spec_v2)

            # 4. Call the SystemList API to add components
            #    This uses the *original* session's system list as the base.
            new_systs_v2 = self._session.systs.add_components_from_regions(
                spec=new_spec_v2, # Pass the *new* spectrum for its region map
                region_id_col='abs_ids',
                series_map=series_map,
                z_map=z_map
            )

            # 5. Return a *new* session with *both* the new spectrum and new system list
            # We must use the base .with_new_spectrum() and .with_new_system_list()
            # to create the final, fully updated SessionV2 state.
            
            # Start from the original session
            session_with_new_spec = self._session.with_new_spectrum(new_spec_v2)
            final_session = session_with_new_spec.with_new_system_list(new_systs_v2)
            
            return final_session
        except Exception as e:
            logging.error(f"Failed during identify_lines: {e}", exc_info=True)
            return 0
    
    def add_components(self,
                       region_id_col: str = 'abs_ids',
                       series_map_json: str = None,
                       z_map_json: str = None) -> 'SessionV2':
        """
        API: Populates the SystemList with components found by identify_lines.
        (This is called internally by the identify_lines orchestrator)
        """
        if not series_map_json or not z_map_json:
            logging.error("add_components: Missing series_map or z_map.")
            return 0
            
        try:
            # Deserialize the maps
            series_map = json.loads(series_map_json)
            z_map = json.loads(z_map_json)
            # Convert string keys back to int
            series_map = {int(k): v for k, v in series_map.items()}
            z_map = {int(k): v for k, v in z_map.items()}

        except Exception as e:
            logging.error(f"add_components: Failed to deserialize maps: {e}")
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
            logging.error(f"Failed during add_components: {e}", exc_info=True)
            return 0
        
    def add_component(self, series: str = 'Ly_a', z: str = '0.0', 
                      logN: str = '13.5', b: str = '10.0', btur: str = '0.0') -> 'SessionV2':
        """
        API: Recipe to add a single component.
        """
        try:
            z_f = float(z)
            logN_f = float(logN)
            b_f = float(b)
            btur_f = float(btur)
        except ValueError:
            logging.error(msg_param_fail)
            return 0
            
        try:
            # Handle case where session has no system list yet
            if self._session.systs is None:
                # We need to bootstrap an empty system list. 
                # Ideally SessionV2 should handle this, but for now we can handle it here 
                # or assume the session creation initialized an empty one.
                # If 'systs' is None, we can't call add_component on it.
                # Let's check `session.py`. If it allows None, we must create one.
                from ..structures import SystemListDataV2
                from ..system_list import SystemListV2
                empty_systs = SystemListV2(SystemListDataV2())
                new_systs = empty_systs.add_component(series, z_f, logN_f, b_f, btur_f)
            else:
                # Call the API layer
                new_systs = self._session.systs.add_component(series, z_f, logN_f, b_f, btur_f)
            
            return self._session.with_new_system_list(new_systs)
            
        except Exception as e:
            logging.error(f"Failed during add_component: {e}", exc_info=True)
            return 0