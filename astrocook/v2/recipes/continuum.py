import logging
from typing import TYPE_CHECKING, Optional
import astropy.units as au

from ..structures import HistoryLogV2
from ...v1.message import msg_param_fail

if TYPE_CHECKING:
    from ..session import SessionV2

# --- Schemas for the Continuum Menu ---
CONTINUUM_RECIPES_SCHEMAS = {
    "estimate_auto": {
        "brief": "Auto-estimate continuum (single-click).",
        "details": "A single-click recipe that runs 'find_unabsorbed' and then 'fit_continuum' in sequence, based on the V1 'clip_flux' algorithm.",
        "params": [
            {"name": "smooth_len_lya", "type": float, "default": 5000.0, "doc": "Smoothing length in Ly-a forest (km/s)"},
            {"name": "smooth_len_out", "type": float, "default": 400.0, "doc": "Smoothing length outside Ly-a forest (km/s)"},
            {"name": "kappa", "type": float, "default": 2.0, "doc": "Sigma threshold for clipping"},
            {"name": "fudge", "type": float, "default": 1.0, "doc": "Continuum fudge factor"},
            {"name": "smooth_std", "type": float, "default": 500.0, "doc": "Final Gaussian smoothing std (km/s)"},
            {"name": "template", "type": bool, "default": False, "doc": "Use QSO template (NOT IMPLEMENTED)"},
        ],
        "url": "continuum_cb.html#estimate_auto"
    },
    "find_unabsorbed": {
        "brief": "Find unabsorbed regions (V1 logic).",
        "details": "Run the V1 'clip_flux' kappa-sigma algorithm (with Ly-a split) to create a 'mask_unabs' column.",
        "params": [
            {"name": "smooth_len_lya", "type": float, "default": 5000.0, "doc": "Smoothing length in Ly-a forest (km/s)"},
            {"name": "smooth_len_out", "type": float, "default": 400.0, "doc": "Smoothing length outside Ly-a forest (km/s)"},
            {"name": "kappa", "type": float, "default": 2.0, "doc": "Sigma threshold for clipping"},
            {"name": "template", "type": bool, "default": False, "doc": "Use QSO template (NOT IMPLEMENTED)"},
        ],
        "url": "continuum_cb.html#find_unabsorbed"
    },
    "fit_continuum": {
        "brief": "Fit continuum to mask (V1 logic).",
        "details": "Interpolate unabsorbed regions from a mask, apply fudge factor, and smooth the result to create 'cont'.",
        "params": [
            {"name": "fudge", "type": float, "default": 1.0, "doc": "Continuum fudge factor"},
            {"name": "smooth_std", "type": float, "default": 500.0, "doc": "Final Gaussian smoothing std (km/s)"},
            {"name": "mask_col", "type": str, "default": "mask_unabs", "doc": "Mask column to use"},
        ],
        "url": "continuum_cb.html#fit_continuum"
    }
}


class RecipeContinuumV2:
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cont' # V1 tag for this menu

    def find_unabsorbed(self, smooth_len_lya: str = '5000.0', smooth_len_out: str = '400.0', 
                        kappa: str = '2.0', template: str = 'False') -> 'SessionV2':
        """
        API: Finds unabsorbed regions using V1 'clip_flux' logic.
        """
        try:
            smooth_len_lya_q = float(smooth_len_lya) * au.km/au.s
            smooth_len_out_q = float(smooth_len_out) * au.km/au.s
            kappa_f = float(kappa)
            template_b = str(template) == 'True'
            z_em_f = self._session.spec._data.z_em # Read z_em from session
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        if z_em_f == 0.0:
            logging.warning("z_em is 0.0. Run 'Edit > Set Properties' before continuum fitting.")
            # We don't stop, as the operation function handles z_em=0.0
            
        try:
            new_spec_v2 = self._session.spec.find_unabsorbed(
                smooth_len_lya=smooth_len_lya_q,
                smooth_len_out=smooth_len_out_q,
                kappa=kappa_f,
                template=template_b
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during find_unabsorbed: {e}", exc_info=True)
            return 0

    def fit_continuum(self, fudge: str = '1.0', smooth_std: str = '500.0', 
                       mask_col: str = 'mask_unabs') -> 'SessionV2':
        """
        API: Fits a continuum to a mask using V1 'interp-and-smooth' logic.
        """
        try:
            fudge_f = float(fudge)
            smooth_std_f = float(smooth_std) # This is in km/s
        except ValueError:
            logging.error(msg_param_fail)
            return 0
        
        try:
            new_spec_v2 = self._session.spec.fit_continuum(
                fudge=fudge_f,
                smooth_std_kms=smooth_std_f,
                mask_col=mask_col
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during fit_continuum: {e}", exc_info=True)
            return 0
            
    def estimate_auto(self, smooth_len_lya: str = '5000.0', smooth_len_out: str = '400.0', 
                      kappa: str = '2.0', fudge: str = '1.0', 
                      smooth_std: str = '500.0', template: str = 'False') -> 'SessionV2':
        """
        API: Single-click recipe to estimate continuum.
        """
        try:
            smooth_len_lya_q = float(smooth_len_lya) * au.km/au.s
            smooth_len_out_q = float(smooth_len_out) * au.km/au.s
            kappa_f = float(kappa)
            template_b = str(template) == 'True'
            fudge_f = float(fudge)
            smooth_std_f = float(smooth_std) # This is in km/s
            z_em_f = self._session.spec._data.z_em # Read z_em from session
        except ValueError:
            logging.error(msg_param_fail)
            return 0
            
        if z_em_f == 0.0:
            logging.warning("z_em is 0.0. Run 'Edit > Set Properties' first.")

        try:
            logging.info("Auto-continuum: Finding unabsorbed regions...")
            # 1. Call the first API method
            spec_with_mask = self._session.spec.find_unabsorbed(
                smooth_len_lya=smooth_len_lya_q,
                smooth_len_out=smooth_len_out_q,
                kappa=kappa_f,
                template=template_b
            )
            
            logging.info("Auto-continuum: Fitting continuum to mask...")
            # 2. Call the second API method *on the result of the first*
            spec_with_cont = spec_with_mask.fit_continuum(
                fudge=fudge_f,
                smooth_std_kms=smooth_std_f,
                mask_col='mask_unabs' # Use the mask we just created
            )
            
            # 3. Return the final new state
            return self._session.with_new_spectrum(spec_with_cont)
        except Exception as e:
            logging.error(f"Failed during estimate_auto: {e}", exc_info=True)
            return 0