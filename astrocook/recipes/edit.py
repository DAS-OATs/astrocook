from astropy import units as au
import logging
import numpy as np
import re
from typing import TYPE_CHECKING, Dict, Optional, Union

from astrocook.core.structures import HistoryLogV2  # <<< *** ADD THIS IMPORT ***
from astrocook.legacy.message import msg_param_fail

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2

EDIT_RECIPES_SCHEMAS = {
    "set_properties": {
        "brief": "Set session properties.",
        "details": "Set core session properties like emission redshift (z_em) and rest-frame redshift (z_rf).",
        "params": [
            {"name": "z_em", "type": float, "default": "_current_", "doc": "Emission Redshift (z_em)"},
            {"name": "resol", "type": float, "default": "_current_", "doc": "Resolving Power R (e.g. 50000). Set to 0 if 'resol' column exists."}
        ],
        "url": "edit_cb.html#set_properties"
    },
    "x_convert": {
        "brief": "Convert X axis.",
        "details": "Convert the x axis units and zero point to wavelength or velocity units.",
        "params": [
            {"name": "z_rf", "type": float, "default": "_current_", "doc": "Emission redshift"},
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
    "smooth_column": {
        "brief": "Smooth a single column.",
        "details": "Apply a Gaussian filter to a single data column (e.g., 'dy' or 'snr').",
        "params": [
            {"name": "target_col", "type": str, "default": "dy", "doc": "Target column to smooth"},
            {"name": "sigma_kms", "type": float, "default": 100.0, "doc": "Standard deviation for Gaussian kernel (km/s)"},
            {"name": "renorm_model", "type": bool, "default": False, "doc": "If smoothing 'cont', also re-normalize 'model'?", "gui_hidden": True},
        ],
        "url": "edit_cb.html#smooth_column" # Placeholder URL
    },
    "split": {
        "brief": "Split spectrum by expression.",
        "details": "Extract a portion of the spectrum into a new session using a boolean expression. E.g., '(x > 400) & (x < 500)'.",
        "params": [
            {"name": "expression", "type": str, "default": "(x > 400) & (x < 500)", "doc": "Boolean expression to select data"}
        ],
        "url": "edit_cb.html#split"
    },
    "extract_preset": {
        "brief": "Extract region by preset.",
        "details": "Extract a specific standard region (like the Ly-alpha forest) into a new session, based on the current emission redshift (z_em).",
        "params": [
            # We use a string for now. Options: 'lya_forest', 'all_ly_forest', 'red_side'
            {"name": "region", "type": str, "default": "lya_forest", "doc": "Preset region ('lya_forest', 'all_ly_forest', 'red_side')"}
        ],
        "url": "edit_cb.html#extract_preset"
    }
}

SESSION_VAR_REGEX = re.compile(r'([a-zA-Z0-9_]+)\s*\.\s*([a-zA-Z0-9_]+)')

class RecipeEditV2:
    """
    Recipes for general structure editing and arithmetic operations.
    Accessed via ``session.edit``.
    """
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cb'

    def _resample_and_replace_match(self, match_obj: re.Match, 
                                  alias_map: Dict[str, str], 
                                  extra_vars: Dict[str, np.ndarray]) -> str:
        """
        Internal helper for regex substitution in multi-session expressions.
        
        This function is called for every pattern match like 's1.y'. It finds
        the referenced session, resamples its data onto the current session's
        grid, adds the data to ``extra_vars``, and returns a safe variable name.
        """
        alias = match_obj.group(1)
        col_name = match_obj.group(2)
        
        if alias not in alias_map:
            return match_obj.group(0)

        safe_var_name = f"{alias}_{col_name}" # e.g., "s1_y"
        
        if safe_var_name in extra_vars:
            return safe_var_name

        try:
            # 1. Find the target session
            full_name = alias_map[alias]
            target_hist = None
            gui = self._session._gui
            if not (gui and hasattr(gui, 'session_histories')):
                raise RuntimeError("GUI context not found.")
                
            for hist in gui.session_histories:
                if hist.display_name == full_name:
                    target_hist = hist
                    break
            
            if not target_hist:
                raise ValueError(f"Could not find session '{full_name}' for alias '{alias}'.")
            
            target_spec = target_hist.current_state.spec
            current_spec = self._session.spec

            if len(target_spec.x) > len(current_spec.x):
                logging.warning(f"Oversampling: Target '{alias}' grid ({len(target_spec.x)} pts) "
                                f"is finer than current grid ({len(current_spec.x)} pts).")
            
            # 2. Get the *current* spec's grid (in nm)
            current_x_grid_nm = current_spec.x.to_value(au.nm)
            
            # 3. Call the lightweight API method
            logging.debug(f"Auto-resampling {alias}.{col_name} onto current grid...")
            resampled_array = target_spec.get_resampled_column(
                col_name,
                current_x_grid_nm
            )
            
            # --- *** THIS IS THE FIX *** ---
            # 4. Check if the column was non-numerical
            if resampled_array is None:
                # Raise a clear error that the user will see in the dialog
                raise ValueError(f"Column '{col_name}' in session '{alias}' is non-numerical and cannot be used in expression.")
            # --- *** END FIX *** ---
                
            # 5. Add the new array to our extra_vars dict
            extra_vars[safe_var_name] = resampled_array
            
            return safe_var_name 
            
        except Exception as e:
            # This will now catch our new ValueError
            logging.error(f"Failed to resample {alias}.{col_name}: {e}", exc_info=True)
            # Re-raise the error so the recipe can catch it
            raise e
    
    def _prepare_expression_contexts(self, expression: str, alias_map: Dict[str, str]) -> (str, Dict[str, np.ndarray]):
        """
        Parses an expression string to handle multi-session variable references.

        It scans for patterns like ``s1.y``, resamples the data from session ``s1``
        onto the current grid, and rewrites the expression to use local variable names.

        Parameters
        ----------
        expression : str
            The user-provided expression (e.g. ``"y - s1.y"``).
        alias_map : dict
            Mapping of aliases to full session names (e.g. ``{'s1': 'QSO_1423'}``).

        Returns
        -------
        final_expression : str
            Rewritten expression (e.g. ``"y - s1_y"``).
        extra_vars : dict
            Dictionary of resampled arrays to be passed to ``numexpr``.
        """
        if not alias_map:
            return expression, {} # No aliases, nothing to do

        extra_vars = {}

        try:
            replacer = lambda m: self._resample_and_replace_match(m, alias_map, extra_vars)
            final_expression = SESSION_VAR_REGEX.sub(replacer, expression)
        except Exception as e:
            # Propagate the clear error message (e.g., "Cannot resample...")
            raise e
        
        # --- *** END FIX *** ---

        return final_expression, extra_vars
    
    def set_properties(self, z_em: str = '0.0', resol: str = '0.0') -> Optional['SessionV2']:
        """
        Sets the core physical properties of the session.

        Parameters
        ----------
        z_em : str (float)
            Emission redshift of the source.
        resol : str (float)
            Resolving power :math:`R = \lambda / \Delta\lambda`.
            If 0.0, the recipe will try to read it from metadata.

        Returns
        -------
        SessionV2
            A new session instance with updated metadata.
        """
        try:
            z_em_f = float(z_em)
            resol_f = float(resol)
        except ValueError:
            logging.error(msg_param_fail)
            return 0
            
        try:
            # Check if values actually changed
            # [FIX] Check against stored R
            current_resol = self._session.spec._data.resol
            if current_resol == 0.0:
                 current_resol = self._session.spec.meta.get('resol', 0.0)

            if z_em_f == self._session.spec._data.z_em and resol_f == current_resol:
                raise ValueError("No changes were made to the properties.")
            
            # [FIX] Pass 'resol' (R) to with_properties
            new_spec = self._session.spec.with_properties(z_em=z_em_f, resol=resol_f)
            return self._session.with_new_spectrum(new_spec)
        except Exception as e:
            raise e

    def x_convert(self, z_rf: str = '0.0', xunit: str = 'km/s') -> Optional['SessionV2']:
        """
        Converts the X-axis (wavelength/velocity) to a new unit.

        Parameters
        ----------
        z_rf : str (float)
            Rest-frame redshift used for velocity conversions.
        xunit : str
            Target unit (e.g. ``'nm'``, ``'Angstrom'``, ``'km/s'``).

        Returns
        -------
        SessionV2
            A new session with the converted X-axis.
        """
        # V1-style validation (simplified)
        try:
            z_rf_float = float(z_rf)
            xunit_obj = au.Unit(xunit)
        except ValueError:
            # Replicates V1 error handling path to avoid crash in the GUI context
            logging.error(msg_param_fail)
            return 0 # Returning 0 prevents crashing and triggers V1 logging

        # 1. Get the new SpectrumV2 instance (the immutable operation)
        new_spec_v2 = self._session.spec.x_convert(z_rf=z_rf_float, xunit=xunit_obj)
        
        # 2. Return a NEW SessionV2 instance with the updated spectrum
        return self._session.with_new_spectrum(new_spec_v2) 
    
    def y_convert(self, yunit: str = 'erg/(nm s cm^2)') -> Optional['SessionV2']:
        """
        Converts the Y-axis (flux) to a new unit.

        Parameters
        ----------
        yunit : str
            Target flux unit.

        Returns
        -------
        SessionV2
            A new session with the converted Y-axis.
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
    
    def apply_expression(self, target_col: str, expression: str, alias_map: Dict[str, str] = None) -> 'SessionV2':
        """
        Applies a numerical expression to columns (arithmetic).

        Supports basic math (``+ - * /``) and NumPy functions (``log``, ``sqrt``, ``min``, etc.).
        Also supports **multi-session operations** via aliases (e.g., ``y - s1.y``).

        Parameters
        ----------
        target_col : str
            The column to create or overwrite (e.g. ``'y'``, ``'ratio'``).
        expression : str
            The mathematical expression string (e.g. ``"y / cont"``).
        alias_map : dict, optional
            A mapping of aliases to session names for cross-session math.
            (Provided automatically by the GUI dialog).

        Returns
        -------
        SessionV2
            A new session with the computed column.
        """
        expression = expression.strip()
        if not expression:
            logging.error("Expression cannot be empty.")
            return 0
        try:
            # 1. Prepare multi-session vars
            final_expression, extra_vars = self._prepare_expression_contexts(expression, alias_map)
            logging.debug(f"Original expression: '{expression}'")
            logging.debug(f"Final expression: '{final_expression}'")
            logging.debug(f"Extra vars: {list(extra_vars.keys())}")

            # 2. Call the API with the prepared data
            new_spec_v2 = self._session.spec.apply_expression(
                target_col=target_col,
                expression=final_expression,
                extra_vars=extra_vars
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during apply_expression: {e}", exc_info=True)
            return 0 
            
    def mask_expression(self, target_col: str, expression: str, alias_map: Dict[str, str] = None) -> 'SessionV2':
        """
        Masks a column (sets values to NaN) based on a boolean expression.

        Parameters
        ----------
        target_col : str
            The column to apply the mask to.
        expression : str
            A boolean expression (e.g. ``"x > 500"``). Where True, the data is masked.
        alias_map : dict, optional
            Mapping for multi-session variables.

        Returns
        -------
        SessionV2
            A new session with the masked data.
        """
        expression = expression.strip()
        if not expression:
            logging.error("Expression cannot be empty.")
            return 0
        try:
            # 1. Prepare multi-session vars
            final_expression, extra_vars = self._prepare_expression_contexts(expression, alias_map)
            logging.debug(f"Original expression: '{expression}'")
            logging.debug(f"Final expression: '{final_expression}'")
            logging.debug(f"Extra vars: {list(extra_vars.keys())}")
            
            # 2. Call the API with the prepared data
            new_spec_v2 = self._session.spec.mask_expression(
                target_col=target_col,
                expression=final_expression,
                extra_vars=extra_vars
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during mask_expression: {e}", exc_info=True)
            return 0

    def smooth_column(self, target_col: str = 'dy', sigma_kms: str = '100.0', renorm_model: str = 'False') -> 'SessionV2':
        """
        Applies Gaussian smoothing to a single column.

        Parameters
        ----------
        target_col : str
            The name of the column to smooth.
        sigma_kms : str (float)
            Standard deviation of the Gaussian kernel in km/s.
        renorm_model : str (bool)
            If 'True' and ``target_col`` is 'cont', re-normalize the 'model' column.

        Returns
        -------
        SessionV2
            A new session with the smoothed column.
        """
        try:
            sigma_kms_f = float(sigma_kms)
            # --- ADD THIS LINE ---
            renorm_model_b = str(renorm_model) == 'True'
        except ValueError:
            logging.error(msg_param_fail)
            return 0
            
        try:
            # 1. Call the immutable V2 operation
            new_spec_v2 = self._session.spec.smooth_column(
                target_col=target_col,
                sigma_kms=sigma_kms_f,
                # --- ADD THIS ARGUMENT ---
                renorm_model=renorm_model_b
            )
            
            # 2. Return a NEW SessionV2 instance
            return self._session.with_new_spectrum(new_spec_v2)
            
        except Exception as e:
            logging.error(f"Failed during smooth_column: {e}", exc_info=True)
            return 0

    def split(self, expression: str, alias_map: Dict[str, str] = None) -> 'SessionV2':
        """
        Splits the spectrum, creating a new session with a subset of the data.

        Parameters
        ----------
        expression : str
            A boolean expression defining the rows to KEEP (e.g. ``"x > 121.6"``).
        alias_map : dict, optional
            Mapping for multi-session variables.

        Returns
        -------
        SessionV2
            A new session containing only the selected data points.
        """
        expression = expression.strip()
        if not expression:
            logging.error("Expression cannot be empty.")
            return 0
            
        # 1. Prepare multi-session vars (re-using our existing helper!)
        final_expression, extra_vars = self._prepare_expression_contexts(expression, alias_map)
        
        # 2. Call the immutable V2 operation
        new_spec_v2 = self._session.spec.split(
            expression=final_expression,
            extra_vars=extra_vars
        )
        
        # 3. Return a NEW SessionV2 instance (GUI handles branching)
        return self._session.with_new_spectrum(new_spec_v2)
    
    def extract_preset(self, region: str = 'lya_forest') -> 'SessionV2':
        """
        Extracts a standard spectral region based on emission redshift (z_em).

        Useful for quickly isolating the Lyman-alpha forest.

        Parameters
        ----------
        region : str
            Preset identifier:
            * ``'lya_forest'``: Between Ly-beta and Ly-alpha emission.
            * ``'all_ly_forest'``: From Lyman limit to Ly-alpha emission.
            * ``'red_side'``: Everything redward of Ly-alpha emission.

        Returns
        -------
        SessionV2
            A new session containing the extracted region.
        """
        z_em = self._session.spec._data.z_em
        if z_em == 0.0:
             logging.error("Cannot extract preset region: z_em is 0.0. Run 'Set Properties' first.")
             return 0

        # 1. Define Presets (Rest-Frame Wavelengths in nm)
        # These are standard approximations.
        # Ly-alpha is ~121.567 nm, Ly-beta is ~102.57 nm, Lyman limit is ~91.18 nm
        presets = {
            'lya_forest': (102.57, 121.567),   # Between Ly-beta and Ly-alpha emissions (with buffer)
            'all_ly_forest': (91.18, 121.567),  # Full series down to the limit
            'red_side': (121.567, 10000.0)    # Everything redward of Ly-alpha emission
        }
        
        if region not in presets:
            logging.error(f"Unknown preset region '{region}'. Available: {list(presets.keys())}")
            return 0
            
        rf_min, rf_max = presets[region]
        
        # 2. Convert to Observed Frame
        obs_min = rf_min * (1.0 + z_em)
        obs_max = rf_max * (1.0 + z_em)
        
        # 3. Build the split expression
        expression = f"(x > {obs_min:.4f}) & (x < {obs_max:.4f})"
        logging.info(f"Extracting '{region}': {expression} (z_em={z_em:.4f})")
        
        # 4. Re-use the existing split logic
        # We can call our own split method directly!
        return self.split(expression=expression)