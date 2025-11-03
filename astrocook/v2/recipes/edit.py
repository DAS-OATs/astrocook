from astropy import units as au
import logging
from typing import TYPE_CHECKING, Optional, Union

from ..structures import HistoryLogV2  # <<< *** ADD THIS IMPORT ***
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
    },
    "apply_expression": {
        "brief": "Apply expression to columns.",
        "details": "Apply a NumPy-style expression. Use column names (x, y, cont...) as variables. E.g., 'y / cont', 'y * 2.0', or 'log10(y)'.",
        "params": [
            {"name": "target_col", "type": str, "default": "y", "doc": "Target column (to be overwritten or created)"},
            {"name": "expression", "type": str, "default": "y / cont", "doc": "Expression to evaluate (e.g., 'y / 2.0')"},
        ],
        "url": "edit_cb.html#apply_expression" # Placeholder URL
    },
    "mask_expression": {
        "brief": "Mask column by expression.",
        "details": "Mask a target column by setting values to NaN where the expression is True. E.g., '(x < 300) | (x > 400)', or 'dy <= 0'.",
        "params": [
            {"name": "target_col", "type": str, "default": "y", "doc": "Target column to mask (e.g., y)"},
            {"name": "expression", "type": str, "default": "dy <= 0", "doc": "Boolean expression (e.g., '(x > 300) & (x < 400)')"},
        ],
        "url": "edit_cb.html#mask_expression" # Placeholder URL
    },
    "split": {
        "brief": "Split spectrum by range.",
        "details": "Extract a portion of the spectrum into a new, separate session.",
        "params": [
            {"name": "col_check", "type": str, "default": "x", "doc": "Column to check for range (e.g., x)"},
            {"name": "min_val", "type": float, "default": 400.0, "doc": "Minimum value to keep (in column's units)"},
            {"name": "max_val", "type": float, "default": 500.0, "doc": "Maximum value to keep (in column's units)"},
        ],
        "url": "edit_cb.html#split" # Placeholder URL
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
    
    def apply_expression(self, target_col: str, expression: str) -> 'SessionV2':
        """
        API: Applies a numerical expression to the spectrum data.
        """
        if not expression:
            logging.error("Expression cannot be empty.")
            return 0
        try:
            new_spec_v2 = self._session.spec.apply_expression(
                target_col=target_col,
                expression=expression
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during apply_expression: {e}", exc_info=True)
            return 0 # V1-style failure
            
    def mask_expression(self, target_col: str, expression: str) -> 'SessionV2':
        """
        API: Masks a column based on a boolean expression.
        """
        if not expression:
            logging.error("Expression cannot be empty.")
            return 0
        try:
            new_spec_v2 = self._session.spec.mask_expression(
                target_col=target_col,
                expression=expression
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during mask_expression: {e}", exc_info=True)
            return 0

    def split(self, col_check: str, min_val: str, max_val: str) -> 'SessionV2':
        """
        API: Splits the spectrum, returning a new spectrum of the sub-range.
        """
        try:
            # 1. Type casting
            min_val_f = float(min_val)
            max_val_f = float(max_val)
            
            # 2. Execute the immutable V2 operation
            new_spec_v2 = self._session.spec.split(
                col_check=col_check,
                min_val=min_val_f,
                max_val=max_val_f
            )
            
            # 3. Return a NEW SessionV2 instance
            # The GUI will handle this as a "branching" action
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during split: {e}", exc_info=True)
            return 0

