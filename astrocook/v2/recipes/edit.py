from astropy import units as au
from typing import TYPE_CHECKING, Optional, Union

from ...v1.message import msg_param_fail

if TYPE_CHECKING:
    from ..session import SessionV2

EDIT_RECIPES_SCHEMAS = {
    "x_convert": {
        "brief": "Convert X axis.",
        "details": "Convert the x axis units and zero point to wavelength or velocity units.",
        "params": [
            {"name": "zem", "type": float, "default": 0.0, "doc": "Emission redshift"},
            {"name": "xunit", "type": str, "default": "km/s", "doc": "Unit of wavelength or velocity"},
        ],
        "url": "edit_cb.html#convert-x-axis"
    },
    "y_convert": {
        "brief": "Convert Y axis units.",
        "details": "Converts the flux and error arrays to the specified unit.",
        "params": [
            {"name": "yunit", "type": str, "default": "erg/(nm s cm^2)", "doc": "Unit of flux density"},
        ],
        "url": "edit_cb.html#convert-y-axis"
    }
}

class RecipeEditV2:
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cb'

    def x_convert(self, zem: str = '0.0', xunit: str = 'km/s') -> Optional['SessionV2']:
        """
        Intercepts the V1 x_convert call, performs validation, and returns a NEW SessionV2.
        """
        # V1-style validation (simplified)
        try:
            zem_float = float(zem)
            xunit_obj = au.Unit(xunit)
        except ValueError:
            # Replicates V1 error handling path to avoid crash in the GUI context
            logging.error(msg_param_fail)
            return 0 # Returning 0 prevents crashing and triggers V1 logging

        # 1. Get the new SpectrumV2 instance (the immutable operation)
        new_spec_v2 = self._session.spec.x_convert(zem=zem_float, xunit=xunit_obj)
        
        # 2. Return a NEW SessionV2 instance with the updated spectrum
        return self._session.with_new_spectrum(new_spec_v2) 
    
    def y_convert(self, yunit: str = 'erg/(nm s cm^2)') -> Optional['SessionV2']:
        """
        [V1-COMPATIBLE SIGNATURE]
        Intercepts the V1 call, performs validation, and returns a NEW SessionV2.
        """
        # C. Parameter Validation and Type Casting
        try:
            yunit_obj = au.Unit(yunit)
        except ValueError:
            logging.error(msg_param_fail)
            return 0 

        # 1. Execute the immutable V2 operation
        new_spec_v2 = self._session.spec.y_convert(yunit=yunit_obj)
        
        # 2. Return a NEW SessionV2 instance (the new state)
        return self._session.with_new_spectrum(new_spec_v2)

