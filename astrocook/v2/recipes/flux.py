from astropy import units as au
import logging
import numpy as np
from typing import TYPE_CHECKING, Optional, Union

from ..spectrum import DataColumnV2, SpectrumDataV2, SpectrumV2
from ..spectrum_operations import rebin_spectrum
from ...v1.message import msg_param_fail

if TYPE_CHECKING:
    from ..session import SessionV2 

FLUX_RECIPES_SCHEMAS = {
    "rebin": {
        "brief": "Rebin spectrum.",
        "details": "Apply a new binning to a spectrum, with a constant bin size.",
        "params": [
            # CRITICAL: Use string "None" to represent the Python None object
            {"name": "xstart", "type": str, "default": "None", "doc": "Start wavelength (nm, None for auto)"},
            {"name": "xend", "type": str, "default": "None", "doc": "End wavelength (nm, None for auto)"},
            
            # CRITICAL: dx needs to be a number (10.0) but the V1 flow expects km/s unit conversion later.
            # We use the raw value here, and the adapter handles the unit.
            {"name": "dx", "type": float, "default": 10.0, "doc": "Step in x"}, 
            
            # Use string "km/s" for the unit input field
            {"name": "xunit", "type": str, "default": "km/s", "doc": "Unit of wavelength or velocity"},
            
            # Use string "None"
            {"name": "kappa", "type": str, "default": "None", "doc": "Sigma clipping threshold (None for off)"},
            
            # CRITICAL: V1 uses 'norm' as a boolean flag. Keep it as string 'False' for dialog logic.
            {"name": "norm", "type": str, "default": "False", "doc": "Return normalized spectrum, if continuum exists"}, 

            # CRITICAL: Use string representation 'nan'
            {"name": "filling", "type": str, "default": "nan", "doc": "Value to fill region without data"}, 
        ],
        "url": "edit_cb.html#rebin"
    }
}


class RecipeFluxV2:
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cb'

    def rebin(self, xstart: str = 'None', xend: str = 'None', dx: str = '10.0', xunit: str = 'km/s',
              kappa: str = '5.0', filling: str = 'nan', norm: str = 'False') -> 'SessionV2':
        """
        API: Rebins the spectrum, logs the action, and returns a NEW SpectrumV2 instance.
        """

        # 2. Log the action *before* executing
        # We log to the session's log instance.
        # append_full also adds the '_refresh' command
        try:
            # self._tag is 'cb' from the __init__
            self._session.log.append_full(self._tag, 'rebin', params)
            logging.debug(f"Logged recipe: {self._tag}.rebin")
        except Exception as e:
            # Don't stop the recipe, but log the failure
            logging.error(f"Failed to log rebin action: {e}")

        try:
            # Type Casting based on string input
            xstart_q = au.Quantity(float(xstart), au.nm) if xstart != 'None' else None
            xend_q = au.Quantity(float(xend), au.nm) if xend != 'None' else None
            
            # V1 logic from cookbook_flux.py uses the unit attached to xunit for dx
            dx_unit = au.Unit(xunit) 
            dx_q = au.Quantity(float(dx), dx_unit) 

            kappa_f = float(kappa) if kappa != 'None' else None
            
            norm_b = str(norm) == 'True' # Handle boolean flag
            filling_f = float(filling) if filling != 'nan' else np.nan
        except ValueError:
            # Replicates V1 error handling path to avoid crash in the GUI context
            logging.error(msg_param_fail)

        # 1. Execute the immutable V2 operation
        new_spec_v2 = self._session.spec.rebin(
            xstart=xstart_q, xend=xend_q, dx=dx_q, kappa=kappa_f, filling=filling_f
        )
        
        # 2. Return a NEW SessionV2 instance (the new state)
        return self._session.with_new_spectrum(new_spec_v2)