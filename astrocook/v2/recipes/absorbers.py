import astropy.units as au
import json
import logging
from typing import TYPE_CHECKING, Optional

from ...v1.vars import xem_d
from ..atomic_data import STANDARD_MULTIPLETS
from ..fitting.voigt_fitter import VoigtFitterV2
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
            {"name": "btur", "type": float, "default": 0.0, "doc": "Turbulent broadening (km/s)"}
        ],
        # We hide this from the menu because it's usually triggered by mouse, 
        # but keeping it in the schema allows scripting/logging.
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
        "gui_hidden": True # Usually triggered by table edit
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
            {"name": "max_nfev", "type": int, "default": 2000, "doc": "Max iterations"}
        ],
        "url": "absorbers_cb.html#fit_component",
        "gui_hidden": True
    }
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
                       debug_rating: str = 'False',
                       auto_populate: str = 'True') -> 'SessionV2':
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
            auto_populate_b = str(auto_populate).lower() == 'true'
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

            if not auto_populate_b:
                logging.info("identify_lines: Auto-populate disabled. Returning spectrum with suggestions.")
                return self._session.with_new_spectrum(new_spec_v2)

            # 3. LINKAGE: Call populate_from_identification
            temp_session = self._session.with_new_spectrum(new_spec_v2)
            
            # Ensure systs is initialized
            if temp_session.systs is None:
                from ..structures import SystemListDataV2
                from ..system_list import SystemListV2
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
    
    def populate_from_identification(self,
                       region_id_col: str = 'abs_ids',
                       series_map_json: str = None,
                       z_map_json: str = None) -> 'SessionV2':
        """
        API: Populates the SystemList using results from identification.
        """
        if not series_map_json or not z_map_json:
            logging.error("populate_from_identification: Missing maps.")
            return 0
            
        try:
            # Deserialize the maps
            series_map = json.loads(series_map_json)
            z_map = json.loads(z_map_json)
            # Convert string keys back to int
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
                      logN: str = '13.5', b: str = '10.0', btur: str = '0.0') -> 'SessionV2':
        """
        API: Recipe to add a single component (Manual).
        Automatically fits the component after addition.
        """
        try:
            z_f = float(z); logN_f = float(logN); b_f = float(b); btur_f = float(btur)
        except ValueError:
            logging.error(msg_param_fail); return 0
            
        try:
            # 1. Add the component
            if self._session.systs is None:
                from ..structures import SystemListDataV2
                from ..system_list import SystemListV2
                empty_systs = SystemListV2(SystemListDataV2())
                new_systs = empty_systs.add_component(series, z_f, logN_f, b_f, btur_f)
            else:
                new_systs = self._session.systs.add_component(series, z_f, logN_f, b_f, btur_f)
            
            # Get the added component's UUID (it's the last one)
            added_comp = new_systs.components[-1]
            
            # Update session with the unfitted component first (so fit_component can find it)
            temp_session = self._session.with_new_system_list(new_systs)
            
            # 2. Run Fit immediately
            # We create a new recipe instance on the temp session to call fit
            temp_recipe = RecipeAbsorbersV2(temp_session)
            
            logging.info(f"Auto-fitting new component {added_comp.series}...")
            final_session = temp_recipe.fit_component(added_comp.uuid)
            
            return final_session
            
        except Exception as e:
            logging.error(f"Failed during add_component: {e}", exc_info=True)
            return 0
        
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
            return self._session.with_new_system_list(new_systs)
        except Exception as e:
            logging.error(f"Failed delete_component: {e}"); return 0

    def fit_component(self, uuid: str, max_nfev: str = '2000') -> 'SessionV2':
        try:
            max_nfev_i = int(max_nfev)
            if not self._session.systs: return 0
            
            # 1. Configure Fluid Group
            self._session.systs.constraint_model.set_active_components([uuid])
            
            # 2. Initialize Fitter
            fitter = VoigtFitterV2(self._session.spec, self._session.systs)
            
            # 3. Run Fit (Unpacking model_flux now)
            new_systs, model_flux, res = fitter.fit(max_nfev=max_nfev_i)
            
            if not res or not res.success:
                logging.warning(f"Fit failed: {res.message if res else 'No result'}")
            
            # 4. Update Spectrum Model
            # We need to create a DataColumnV2 for the model and update the spectrum
            from ..structures import DataColumnV2
            from copy import deepcopy
            import dataclasses
            
            new_aux_cols = deepcopy(self._session.spec._data.aux_cols)
            new_aux_cols['model'] = DataColumnV2(
                values=model_flux,
                unit=self._session.spec.y.unit,
                description="Voigt Fit Model"
            )
            
            new_spec_data = dataclasses.replace(self._session.spec._data, aux_cols=new_aux_cols)
            new_spec = self._session.spec.__class__(new_spec_data) # SpectrumV2(new_data)
            
            # 5. Return session with BOTH updated
            session_with_spec = self._session.with_new_spectrum(new_spec)
            return session_with_spec.with_new_system_list(new_systs)
            
        except Exception as e:
            logging.error(f"Failed fit_component: {e}", exc_info=True); return 0