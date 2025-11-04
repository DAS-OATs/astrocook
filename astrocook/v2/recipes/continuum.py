import logging
from typing import TYPE_CHECKING, Optional

from ..structures import HistoryLogV2

if TYPE_CHECKING:
    from ..session import SessionV2

# --- Schemas for the Continuum Menu ---
CONTINUUM_RECIPES_SCHEMAS = {
    "find_unabsorbed": {
        "brief": "Find unabsorbed regions (kappa-sigma).",
        "details": "Run a kappa-sigma clipping algorithm to identify unabsorbed regions and create a 'mask_unabs' column.",
        "params": [
            {"name": "hwindow", "type": int, "default": 100, "doc": "Half-window size for running mean (pixels)"},
            {"name": "kappa", "type": float, "default": 2.0, "doc": "Sigma threshold for clipping"},
        ],
        "url": "continuum_cb.html#find_unabsorbed"
    },
    "fit_polynomial": {
        "brief": "Fit polynomial to mask.",
        "details": "Fit a polynomial to the regions marked 'True' in the 'mask_unabs' column and save it as 'cont'.",
        "params": [
            {"name": "deg", "type": int, "default": 5, "doc": "Degree of the polynomial"},
            {"name": "mask_col", "type": str, "default": "mask_unabs", "doc": "Mask column to use"},
        ],
        "url": "continuum_cb.html#fit_polynomial"
    },
    "estimate_auto": {
        "brief": "Auto-estimate continuum (single-click).",
        "details": "A single-click recipe that runs 'find_unabsorbed' and then 'fit_polynomial' in sequence.",
        "params": [
            {"name": "hwindow", "type": int, "default": 100, "doc": "Half-window size for running mean (pixels)"},
            {"name": "kappa", "type": float, "default": 2.0, "doc": "Sigma threshold for clipping"},
            {"name": "deg", "type": int, "default": 5, "doc": "Degree of the polynomial"},
        ],
        "url": "continuum_cb.html#estimate_auto"
    }
}


class RecipeContinuumV2:
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cont' # V1 tag for this menu

    def find_unabsorbed(self, hwindow: str = '100', kappa: str = '2.0') -> 'SessionV2':
        """
        API: Finds unabsorbed regions.
        """
        try:
            hwindow_i = int(hwindow)
            kappa_f = float(kappa)
        except ValueError:
            logging.error("Invalid parameters: hwindow must be int, kappa must be float.")
            return 0

        try:
            # 1. Call the API method
            new_spec_v2 = self._session.spec.find_unabsorbed(
                hwindow=hwindow_i,
                kappa=kappa_f
            )
            # 2. Return new state
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during find_unabsorbed: {e}", exc_info=True)
            return 0

    def fit_polynomial(self, deg: str = '5', mask_col: str = 'mask_unabs') -> 'SessionV2':
        """
        API: Fits a polynomial to a mask.
        """
        try:
            deg_i = int(deg)
        except ValueError:
            logging.error("Invalid parameter: deg must be an integer.")
            return 0
        
        try:
            # 1. Call the API method
            new_spec_v2 = self._session.spec.fit_polynomial(
                deg=deg_i,
                mask_col=mask_col
            )
            # 2. Return new state
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during fit_polynomial: {e}", exc_info=True)
            return 0
            
    def estimate_auto(self, hwindow: str = '100', kappa: str = '2.0', deg: str = '5') -> 'SessionV2':
        """
        API: Single-click recipe to estimate continuum.
        """
        try:
            hwindow_i = int(hwindow)
            kappa_f = float(kappa)
            deg_i = int(deg)
        except ValueError:
            logging.error("Invalid parameters for auto_estimate.")
            return 0

        try:
            logging.info("Auto-continuum: Finding unabsorbed regions...")
            # 1. Call the first API method
            spec_with_mask = self._session.spec.find_unabsorbed(
                hwindow=hwindow_i,
                kappa=kappa_f
            )
            
            logging.info("Auto-continuum: Fitting polynomial to mask...")
            # 2. Call the second API method *on the result of the first*
            spec_with_cont = spec_with_mask.fit_polynomial(
                deg=deg_i,
                mask_col='mask_unabs' # Use the mask we just created
            )
            
            # 3. Return the final new state
            return self._session.with_new_spectrum(spec_with_cont)
        except Exception as e:
            logging.error(f"Failed during estimate_auto: {e}", exc_info=True)
            return 0