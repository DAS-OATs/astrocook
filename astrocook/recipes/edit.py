from astropy import units as au
from copy import deepcopy
import dataclasses
import logging
import numpy as np
import re
from typing import TYPE_CHECKING, Any, Dict, Optional, Union, List, Tuple

from astrocook.core.structures import HistoryLogV2, DataColumnV2, SpectrumDataV2, SystemListDataV2  # <<< *** ADD THIS IMPORT ***
from astrocook.core.system_list import SystemListV2
from astrocook.legacy.message import msg_param_fail

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2

EDIT_RECIPES_SCHEMAS = {
    "set_properties": {
        "brief": "Set session properties.",
        "details": "Set core session properties like name, object, redshift, and resolution.",
        "params": [
            {"name": "name", "type": str, "default": "_current_", "doc": "Session Name (display label)"},
            {"name": "object", "type": str, "default": "_current_", "doc": "Object Name (header)"},
            {"name": "z_em", "type": str, "default": "_current_", "doc": "Emission Redshift (z_em)"},
            {"name": "resol", "type": str, "default": "_current_", "doc": "Resolving Power R (e.g. 50000) OR Pixel FWHM (e.g. '3px')"}, 
            {"name": "overwrite_resol_col", "type": bool, "default": True, "gui_hidden": True} # <--- [NEW] Add this hidden flag
        ],
        "url": "edit_cb.html#set_properties"
    },
    "x_convert": {
        "brief": "Convert X axis.",
        "details": "Convert the x axis units and zero point to wavelength or velocity units.",
        "params": [
            {"name": "z_ref", "type": float, "default": "_current_", "doc": "Reference redshift for velocity conversion (defaults to z_em)"},
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
        "details": "Apply a NumPy-style expression. Use column names (λ, F, cont...) as variables. E.g., 'F / cont', 'F * 2.0', or 'log10(F)'.",
        "params": [
            {"name": "target_col", "type": str, "default": "F", "doc": "Target column (to be overwritten or created)"},
            {"name": "expression", "type": str, "default": "F / dF", "doc": "Expression to evaluate (e.g., 'F / 2.0')"},
        ],
        "url": "edit_cb.html#apply_expression" # Placeholder URL
    },
    "mask_expression": {
        "brief": "Mask column by expression.",
        "details": "Mask a target column by setting values to NaN where the expression is True. E.g., '(λ < 300) | (λ > 400)', or 'dF <= 0'.",
        "params": [
            {"name": "target_col", "type": str, "default": "F", "doc": "Target column to mask (e.g., F)"},
            {"name": "expression", "type": str, "default": "dF <= 0", "doc": "Boolean expression (e.g., '(λ > 300) & (λ < 400)')"},
        ],
        "url": "edit_cb.html#mask_expression" # Placeholder URL
    },
    "smooth_column": {
        "brief": "Smooth a single column.",
        "details": "Apply a Gaussian filter to a single data column (e.g., 'dF' or 'snr').",
        "params": [
            {"name": "target_col", "type": str, "default": "dF", "doc": "Target column to smooth"},
            {"name": "sigma_kms", "type": float, "default": 100.0, "doc": "Standard deviation for Gaussian kernel (km/s)"},
            {"name": "renorm_model", "type": bool, "default": False, "doc": "If smoothing 'cont', also re-normalize 'model'?", "gui_hidden": True},
        ],
        "url": "edit_cb.html#smooth_column" # Placeholder URL
    },
    "import_systems": {
        "brief": "Import systems from session.",
        "details": "Copy absorption systems from another open session into this one.",
        "params": [
            {"name": "source_session", "type": str, "default": "None", "doc": "Name of the session to import from."},
            {"name": "append", "type": bool, "default": True, "doc": "If True, adds to existing. If False, replaces them."},
            {"name": "refit", "type": bool, "default": False, "doc": "Refit the systems to the new data immediately?"}
        ],
        "url": "absorbers_cb.html#import_systems"
    },
    "delete": {
        "brief": "Delete elements.",
        "details": "Remove columns from the spectrum or clear the line list.",
        "params": [
            {"name": "targets", "type": str, "default": "", "doc": "Comma-separated list of items to delete (e.g. 'cont, lines')"}
        ],
        "url": "edit_cb.html#delete"
    },
    "split": {
        "brief": "Split spectrum by expression.",
        "details": "Extract a portion of the spectrum into a new session using a boolean expression. E.g., '(λ > 400) & (λ < 500)'.",
        "params": [
            {"name": "expression", "type": str, "default": "(λ > 400) & (λ < 500)", "doc": "Boolean expression to select data"}
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
    },
    "trim_common": {
        "brief": "Trim to common velocity coverage.",
        "details": "Trims the current session to the velocity range shared by a list of other sessions (intersection).",
        "params": [
            {"name": "others_names", "type": str, "default": "", "doc": "Comma-separated names of other sessions"},
            {"name": "z_target", "type": float, "default": "_current_", "doc": "Target redshift"},
            {"name": "trans_self", "type": str, "default": "Ly_a", "doc": "Transition for current session"},
            {"name": "trans_others", "type": str, "default": "", "doc": "Comma-separated transitions for other sessions"},
            {"name": "window_kms", "type": float, "default": 500.0, "doc": "Symmetric window width (km/s)"}
        ],
        "url": "edit_cb.html#trim_common"
    },
    "stitch": {
        "brief": "Merge multiple sessions",
        "params": [],
        "gui_hidden": True  # Usually called via script
    },
    "coadd": {
        "brief": "Co-add multiple sessions.",
        "details": "Stitches multiple sessions together and rebins them onto a single common grid in one pass (Inverse Variance Weighting).",
        "params": [
            {"name": "session_names", "type": str, "default": "", "doc": "Comma-separated names of sessions to co-add"},
            {"name": "xstart", "type": str, "default": "None", "doc": "Start wavelength (nm, None for auto)"},
            {"name": "xend", "type": str, "default": "None", "doc": "End wavelength (nm, None for auto)"},
            {"name": "dx", "type": float, "default": 10.0, "doc": "Step size for the final grid"},
            {"name": "xunit", "type": str, "default": "km/s", "doc": "Unit for step size (km/s or nm)"},
            {"name": "kappa", "type": str, "default": "5.0", "doc": "Sigma clipping threshold (None for off)"},
            {"name": "equalize_order", "type": int, "default": 0, "doc": "Polynomial order for flux scaling (-1=Off, 0=Scalar, 1=Linear, etc.)"},
        ],
        "url": "edit_cb.html#coadd"
    }
}

SESSION_VAR_REGEX = re.compile(r'([a-zA-Z0-9_]+)\s*\.\s*([a-zA-Z0-9_]+)')

class RecipeEditV2:
    """
    Recipes for general structure editing and arithmetic operations.

    These methods are accessed via the ``session.edit`` attribute and delegate
    heavy logic to the Core structures :class:`~astrocook.core.spectrum.SpectrumV2`
    and :class:`~astrocook.core.system_list.SystemListV2`.
    """
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cb'

    def _translate_ui_names(self, name_str: str) -> str:
        """
        Translates scientific symbols back to internal keys.
        Supports: λ -> x, F -> y, dF -> dy.
        Original V1 names (x, y, dy) remain untouched for compatibility.
        """
        if not isinstance(name_str, str):
            return name_str

        # Mapping: UI Display -> Internal Key
        mapping = {
            'λ': 'x',
            'λmin': 'xmin',
            'λmax': 'xmax',
            'F': 'y',
            'dF': 'dy'
        }
        
        processed_str = name_str
        # Sort keys by length (descending)
        sorted_keys = sorted(mapping.keys(), key=len, reverse=True)
        
        for ui_name in sorted_keys:
            internal_name = mapping[ui_name]
            # Use regex word boundaries to ensure we only replace exact matches
            processed_str = re.sub(r'\b' + re.escape(ui_name) + r'\b', internal_name, processed_str)
        
        return processed_str

    def _prepare_expression_contexts(self, expression: str, alias_map: Dict[str, str]) -> Tuple[str, Dict[str, np.ndarray]]:
        """
        Internal helper to resolve session aliases and delegate parsing to the Core.

        This prepares the context for multi-session arithmetic (e.g., subtracting
        one spectrum from another) by resolving GUI session names to actual
        Spectrum objects.

        Parameters
        ----------
        expression : str
            The raw expression string provided by the user (e.g. ``"F - s1.F"``).
        alias_map : dict
            A dictionary mapping aliases (e.g. ``"s1"``) to full session names.

        Returns
        -------
        tuple
            A tuple containing:
            1. **final_expression** (*str*): The processed expression compatible with ``numexpr``.
            2. **extra_vars** (*dict*): A dictionary of numpy arrays representing the
               resampled data from other sessions.

        See Also
        --------
        :meth:`astrocook.core.spectrum.SpectrumV2.prepare_expression_context`
        """
        if not alias_map:
            return expression, {}

        # 1. Resolve Aliases to Objects (GUI Logic)
        resolved_specs = {}
        
        # Access the GUI context
        gui = getattr(self._session, '_gui', None)
        if not gui or not hasattr(gui, 'session_histories'):
            # If running in script without GUI, we can't resolve aliases easily 
            # unless we passed the objects directly. For now, assume GUI context or fail gracefully.
            if alias_map:
                logging.warning("GUI context missing: Cannot resolve multi-session aliases.")
            return expression, {}

        for alias, full_name in alias_map.items():
            found = False
            for hist in gui.session_histories:
                if hist.display_name == full_name:
                    resolved_specs[alias] = hist.current_state.spec
                    found = True
                    break
            if not found:
                logging.warning(f"Could not find session '{full_name}' for alias '{alias}'.")

        # 2. Delegate to Core (Business Logic)
        return self._session.spec.prepare_expression_context(expression, resolved_specs)
    
    def _parse_resolution_string(self, resol_str: str) -> float:
        """
        Parses user input for resolution.

        Accepts either a raw floating point number (R) or a pixel width string
        (e.g., "3px") which is converted to R based on the current spectral sampling.

        Parameters
        ----------
        resol_str : str
            Input string, either a raw number (R) or a pixel width (e.g., "3px").

        Returns
        -------
        float
            The resolving power R. Returns 0.0 if parsing fails or spectrum is invalid.
        """
        clean_str = resol_str.strip().lower()
        
        # Check for pixel syntax
        if 'px' in clean_str or 'pixel' in clean_str:
            # Extract the number part
            match = re.match(r"([\d\.]+)", clean_str)
            if not match:
                raise ValueError(f"Could not parse pixel width from '{resol_str}'")
            
            pix_fwhm = float(match.group(1))
            logging.info(f"Estimating Resolution (R) assuming {pix_fwhm:.1f} pixel FWHM...")
            
            # 1. Get raw x array
            x_vals = self._session.spec.x.value
            
            if len(x_vals) > 1:
                # 2. Calculate Median Sampling (dx)
                dx = np.nanmedian(np.diff(x_vals))
                
                # 3. Calculate Median Wavelength (lambda)
                mid_x = np.nanmedian(x_vals)
                
                if dx > 0 and mid_x > 0:
                    # 4. Calculate R = lambda / (N_pix * dx)
                    fwhm_wave = pix_fwhm * dx
                    estimated_R = mid_x / fwhm_wave
                    
                    logging.warning(f"  -> Estimated R ~ {estimated_R:.0f} (at {mid_x:.1f}). "
                                    f"This is an approximation assuming constant sampling.")
                    return estimated_R
                else:
                    logging.warning("Cannot estimate resolution: invalid grid spacing or wavelength.")
                    return 0.0
            else:
                logging.warning("Cannot estimate resolution: spectrum too short.")
                return 0.0
                
        # Fallback: Try standard float conversion
        try:
            return float(clean_str)
        except ValueError:
            raise ValueError(f"Invalid resolution format: '{resol_str}'. Use a number (R) or 'Npx'.")
        

    def set_properties(self, name: str = '_current_', object: str = '_current_', 
                       z_em: str = '_current_', resol: str = '_current_',
                       overwrite_resol_col: bool = True) -> Optional['SessionV2']:
        """
        Set session properties.

        Set core session properties like name, object, redshift, and resolution.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.with_properties`.

        Parameters
        ----------
        name : str, optional
            The session Name (display label). Defaults to ``'_current_'`` (no change).
        object : str, optional
            The object Name (header target). Defaults to ``'_current_'`` (no change).
        z_em : str, optional
            Emission Redshift. Defaults to ``'_current_'`` (no change).
        resol : str, optional
            Resolving Power R (e.g., ``50000``) OR Pixel FWHM (e.g., ``'3px'``).
            Defaults to ``'_current_'`` (no change).
        overwrite_resol_col : bool, optional
            If ``True``, updates the 'resol' column in the spectrum data with the new resolution value. Defaults to ``True``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with updated properties,
            or 0 on failure.
        """
        # 1. Parse numerical inputs (z_em, resol)
        try:
            if str(z_em) == '_current_':
                z_em_f = self._session.spec._data.z_em
            else:
                z_em_f = float(z_em)
                
            if str(resol) == '_current_':
                # Use stored resolution (R) directly
                resol_f = self._session.spec._data.resol
            else:
                resol_f = self._parse_resolution_string(str(resol))
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        # 2. Handle text inputs
        current_meta = self._session.spec.meta
        
        # Determine new values, falling back to current if '_current_' is passed
        new_name = name if name != '_current_' else self._session.name
        new_object = object if object != '_current_' else current_meta.get('OBJECT', '')
        
        try:
            # 3. Apply changes via SpectrumV2.with_properties
            new_spec = self._session.spec.with_properties(
                z_em=z_em_f, 
                resol=resol_f,
                object_name=new_object,
                overwrite_resol_col=overwrite_resol_col 
            )
            
            # 4. Handle Session Name (Lives on the Session wrapper)
            new_session = self._session.with_new_spectrum(new_spec)
            if new_name:
                new_session.name = new_name
            
            return new_session
            
        except Exception as e:
            raise e

    def x_convert(self, z_ref: str = '0.0', xunit: str = 'km/s') -> Optional['SessionV2']:
        """
        Convert X axis.

        Convert the x axis units and zero point to wavelength or velocity units.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.x_convert`.

        Parameters
        ----------
        z_ref : str, optional
            Reference redshift for velocity conversion. Defaults to ``'0.0'``
            (can use ``'_current_'`` to use the session's z_em).
        xunit : str, optional
            Unit of wavelength or velocity (e.g., ``'nm'``, ``'km/s'``).
            Defaults to ``'km/s'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with the converted X-axis,
            or 0 on failure.
        """
        # Handle default logic here
        if z_ref == '_current_':
            z_ref_float = self._session.spec._data.z_em
        else:
            try:
                z_ref_float = float(z_ref)
            except ValueError:
                logging.error(msg_param_fail)
                return 0
            
        try:
            xunit_obj = au.Unit(xunit)
        except ValueError:
            # Replicates V1 error handling path to avoid crash in the GUI context
            logging.error(msg_param_fail)
            return 0 # Returning 0 prevents crashing and triggers V1 logging

        # 1. Get the new SpectrumV2 instance (the immutable operation)
        new_spec_v2 = self._session.spec.x_convert(z_ref=z_ref_float, xunit=xunit_obj)
        
        # 2. Return a NEW SessionV2 instance with the updated spectrum
        return self._session.with_new_spectrum(new_spec_v2) 
    
    def y_convert(self, yunit: str = 'erg/(nm s cm^2)') -> Optional['SessionV2']:
        """
        Convert Y axis units.

        Converts the flux and error arrays to the specified unit.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.y_convert`.

        Parameters
        ----------
        yunit : str, optional
            Target unit of flux density. Defaults to ``'erg/(nm s cm^2)'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with the converted Y-axis,
            or 0 on failure.
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
        Apply expression to columns.

        Apply a NumPy-style expression. Use column names (λ, F, cont...) as variables.
        Supports multi-session arithmetic via aliases.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.apply_expression`.

        Parameters
        ----------
        target_col : str
            Target column (to be overwritten or created).
        expression : str
            Expression to evaluate (e.g., ``'F / dF'``, ``'F * 2.0'``).
        alias_map : dict, optional
            Mapping of aliases to session names (e.g. ``{'s1': 'Session 1'}``).

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with the computed column,
            or 0 on failure.
        """
        # Translate UI symbols back to internal keys
        target_col = self._translate_ui_names(target_col.strip())
        expression = self._translate_ui_names(expression.strip())
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
        Mask column by expression.

        Mask a target column by setting values to NaN where the expression is True.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.mask_expression`.

        Parameters
        ----------
        target_col : str
            Target column to mask (e.g., ``'F'``).
        expression : str
            Boolean expression (e.g., ``'(λ < 300) | (λ > 400)'``, ``'dF <= 0'``).
        alias_map : dict, optional
            Mapping of aliases to session names.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with the masked column,
            or 0 on failure.
        """
        # Translate UI symbols back to internal keys
        target_col = self._translate_ui_names(target_col.strip())
        expression = self._translate_ui_names(expression.strip())
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

    def smooth_column(self, target_col: str = 'dF', sigma_kms: str = '100.0', renorm_model: str = 'False') -> 'SessionV2':
        """
        Smooth a single column.

        Apply a Gaussian filter to a single data column.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.smooth_column`.

        Parameters
        ----------
        target_col : str, optional
            Target column to smooth. Defaults to ``'dF'``.
        sigma_kms : str, optional
            Standard deviation for Gaussian kernel (km/s). Defaults to ``'100.0'``.
        renorm_model : str, optional
            If ``'True'`` and smoothing ``'cont'``, the ``'model'`` column is
            re-normalized. Defaults to ``'False'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with the smoothed column,
            or 0 on failure.
        """
        # Translate UI symbols back to internal keys
        target_col = self._translate_ui_names(target_col.strip())
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
        
    def import_systems(self, source_session: str, append: bool = True, refit: bool = False) -> 'SessionV2':
        """
        Import systems from session.

        Copy absorption systems from another open session into this one, filtering
        those that fall outside the current spectral range.
        Delegates filtering to :meth:`astrocook.core.system_list.SystemListV2.filter_by_range`.

        Parameters
        ----------
        source_session : str
            Name of the session to import from.
        append : bool, optional
            If ``True``, adds to existing systems. If ``False``, replaces them.
            Defaults to ``True``.
        refit : bool, optional
            If ``True``, runs a refit on the imported systems. Defaults to ``False``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with imported systems.
        """

        do_append = str(append) == 'True'
        do_refit = str(refit) == 'True'

        import numpy as np
        from astrocook.core.structures import SystemListDataV2
        from astrocook.core.system_list import SystemListV2

        # 1. Find the source session
        target_hist = None
        if hasattr(self._session, '_gui'):
            for hist in self._session._gui.session_histories:
                if hist.display_name == source_session:
                    target_hist = hist
                    break
        
        if not target_hist:
            logging.error(f"Session '{source_session}' not found.")
            return self._session

        source_systs = target_hist.current_state.systs
        if not source_systs or not source_systs.components:
            logging.warning(f"Source session '{source_session}' has no systems.")
            return self._session

        # 2. Prepare current list
        if do_append and self._session.systs:
            base_systs = self._session.systs
        else:
            base_systs = SystemListV2(SystemListDataV2())

        # 3. Filter Source Components (Delegated to Core)
        if self._session.spec and len(self._session.spec.x) > 0:
            try:
                xmin = np.min(self._session.spec.x.value)
                xmax = np.max(self._session.spec.x.value)
                
                filtered_source_systs = source_systs.filter_by_range(xmin, xmax)
                
                skipped = len(source_systs.components) - len(filtered_source_systs.components)
                if skipped > 0:
                    logging.info(f"Skipped {skipped} systems outside spectral range.")
                    
                source_systs = filtered_source_systs 
                
            except Exception as e:
                logging.warning(f"Could not filter systems by range: {e}")

        if not source_systs.components:
            logging.warning("No systems fell within the target spectral range.")
            return self._session

        # 4. Merge
        logging.info(f"Importing {len(source_systs.components)} systems from '{source_session}'...")
        new_systs = base_systs.merge(source_systs, copy_uuids=False)
        
        # Create intermediate session with imported (but unfitted) systems
        intermediate_session = self._session.with_new_system_list(new_systs)

        # 5. Optional Refit
        if do_refit:
            logging.info("Refitting imported systems...")
            try:
                # Local import to avoid circular dependency
                from astrocook.recipes.absorbers import RecipeAbsorbersV2
                
                # Wrap the intermediate session in the Absorbers recipe
                rec_abs = RecipeAbsorbersV2(intermediate_session)
                
                # Run refit_all
                # You can tweak z_window_kms or max_nfev here if needed, 
                # or expose them as advanced params in the schema.
                refitted_session = rec_abs.refit_all(z_window_kms=20.0, group_depth=1)
                
                if refitted_session != 0:
                    return refitted_session
                else:
                    logging.warning("Refit failed. Returning imported systems without fitting.")
                    return intermediate_session
            except Exception as e:
                logging.error(f"Error during import auto-refit: {e}")
                return intermediate_session

        # Return unfitted result
        return intermediate_session
        
    def delete(self, targets: str) -> 'SessionV2':
        """
        Delete elements.

        Remove columns from the spectrum or clear the line list.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.remove_columns`.

        Parameters
        ----------
        targets : str
            Comma-separated list of items to delete (e.g. ``'cont, lines'``).
            Use ``'lines'`` or ``'systems'`` to clear the component list.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with the elements removed.
        """
        targets = self._translate_ui_names(targets)
        target_list = [t.strip() for t in targets.split(',') if t.strip()]
        if not target_list:
            return self._session

        # 1. Handle System List
        new_systs = self._session.systs
        if 'lines' in target_list or 'systems' in target_list:
            # Re-initialize empty
            new_systs = SystemListV2(SystemListDataV2())
            logging.info("Marked system list for deletion.")
            # Remove keywords from list to avoid trying to delete them as columns
            target_list = [t for t in target_list if t not in ['lines', 'systems']]

        # 2. Handle Columns via Core API
        if target_list:
            new_spec = self._session.spec.remove_columns(target_list)
        else:
            new_spec = self._session.spec
            
        return self._session.with_new_spectrum(new_spec).with_new_system_list(new_systs)
    
    def split(self, expression: str, alias_map: Dict[str, str] = None) -> 'SessionV2':
        """
        Split spectrum by expression.

        Extract a portion of the spectrum into a new session using a boolean expression.
        Only rows where the expression is True are kept.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.split`.

        Parameters
        ----------
        expression : str
            Boolean expression to select data (e.g., ``'(λ > 400) & (λ < 500)'``).
        alias_map : dict, optional
            Mapping for multi-session variables.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` containing the subset of data,
            or 0 on failure.
        """
        expression = self._translate_ui_names(expression.strip())
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

        new_systs = self._session.systs
        if new_spec_v2 and len(new_spec_v2.x) > 0 and new_systs:
            import numpy as np
            import astropy.units as au
            
            # Safely get min/max in nanometers, regardless of current plot units
            xmin_nm = np.min(new_spec_v2.x.to(au.nm).value)
            xmax_nm = np.max(new_spec_v2.x.to(au.nm).value)
            
            original_count = len(new_systs.components)
            new_systs = new_systs.filter_by_range(xmin_nm, xmax_nm)
            filtered_count = len(new_systs.components)
            
            if original_count != filtered_count:
                logging.info(f"Split: Filtered out {original_count - filtered_count} systems outside the new range.")
        
        # 4. Return a NEW SessionV2 instance (GUI handles branching)
        return self._session.with_new_spectrum(new_spec_v2).with_new_system_list(new_systs)
    
    def extract_preset(self, region: str = 'lya_forest') -> Union['SessionV2', List['SessionV2']]:
        """
        Extract region by preset.

        Extract a specific standard region (like the Ly-alpha forest) into a new
        session, based on the current emission redshift (z_em).
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.get_region_bounds`
        and :meth:`split`.

        Parameters
        ----------
        region : str, optional
            Preset region identifier. Options: ``'lya_forest'``, ``'all_ly_forest'``,
            ``'red_side'``. Can be a comma-separated list. Defaults to ``'lya_forest'``.

        Returns
        -------
        SessionV2 or List[SessionV2] or int
            One or more new :class:`~astrocook.core.session.SessionV2` instances,
            or 0 on failure.
        """
        # Parse inputs (split by comma and strip whitespace)
        region_list = [r.strip() for r in region.split(',') if r.strip()]
        
        results = []
        
        for reg in region_list:
            try:
                # Calculate bounds
                obs_min, obs_max = self._session.spec.get_region_bounds(reg)
                
                # Convert to native units for expression
                spec_unit = self._session.spec.x.unit
                v_min = obs_min.to(spec_unit).value
                v_max = obs_max.to(spec_unit).value
                
                # Build expression
                expression = f"(λ > {v_min:.4f}) & (λ < {v_max:.4f})"
                logging.info(f"Extracting '{reg}': {expression} (z_em={self._session.spec._data.z_em:.4f})")
                
                # Create the split session
                new_session = self.split(expression=expression)
                
                # (Optional) Store the intended suffix name in the session metadata 
                # so the GUI could potentially use it later (though GUI currently auto-names branches)
                # new_session.name += f"_{reg}" 
                
                results.append(new_session)
                
            except ValueError as e:
                logging.error(f"Skipping region '{reg}': {e}")
            except Exception as e:
                logging.error(f"Failed to extract '{reg}': {e}", exc_info=True)

        if not results:
            return 0 # Return failure code if nothing worked

        # Return single object if one, or list if multiple
        if len(results) == 1:
            return results[0]
        
        return results
    
    def trim_common(self, others_names: str, z_target: str, trans_self: str, 
                    trans_others: str, window_kms: str = '500.0') -> 'SessionV2':
        """
        Trim to common velocity coverage.

        Trims the current session to the velocity range shared by a list of other
        sessions (intersection).
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.calc_intersection_expression`.

        Parameters
        ----------
        others_names : str
            Comma-separated names of other sessions.
        z_target : str
            Target redshift (float parsed from string).
        trans_self : str
            Transition for current session (e.g. ``'Ly_a'``).
        trans_others : str
            Comma-separated transitions for other sessions.
        window_kms : str, optional
            Symmetric window width in km/s. Defaults to ``'500.0'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` trimmed to the common coverage,
            or 0 on failure.
        """
        try:
            # 1. Parse Inputs
            z_f = float(z_target)
            win_f = float(window_kms)
            others_list = [n.strip() for n in others_names.split(',') if n.strip()]
            trans_list_others = [t.strip() for t in trans_others.split(',') if t.strip()]
            
            # 2. Resolve Session Objects (GUI Logic)
            resolved_others_specs = []
            if hasattr(self._session, '_gui'):
                for name in others_list:
                    found = False
                    for hist in self._session._gui.session_histories:
                        if hist.display_name == name:
                            resolved_others_specs.append(hist.current_state.spec)
                            found = True
                            break
                    if not found:
                        logging.warning(f"trim_common: Session '{name}' not found in GUI history.")
            
            if len(resolved_others_specs) != len(trans_list_others):
                logging.error(f"Mismatch: Found {len(resolved_others_specs)} sessions but {len(trans_list_others)} transitions provided.")
                return 0
                
            # 3. Delegate Calculation to Core (SSOT)
            # The core calculates the physics and returns the boolean expression string
            split_expression = self._session.spec.calc_intersection_expression(
                others_specs=resolved_others_specs,
                z_target=z_f,
                trans_self=trans_self,
                trans_others=trans_list_others,
                window_kms=win_f
            )
            
            if not split_expression:
                logging.warning("trim_common: No overlap found or invalid inputs.")
                return self._session

            logging.info(f"Applying intersection trim: {split_expression}")
            
            # 4. Execute Split
            return self.split(split_expression)

        except Exception as e:
            logging.error(f"trim_common failed: {e}", exc_info=True)
            return 0
    
    def stitch(self, other_sessions: List[Union['SessionV2', str]], sort: bool = True) -> 'SessionV2':
        """
        Merge multiple sessions.

        Concatenates the current session with a list of other sessions.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.stitch`.

        Parameters
        ----------
        other_sessions : list
            A list of :class:`~astrocook.core.session.SessionV2` objects OR strings
            (names of sessions loaded in the GUI).
        sort : bool, optional
            If ``True``, sorts the final spectrum by wavelength. Defaults to ``True``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` containing the merged data,
            or 0 on failure.
        """
        # Delayed import to avoid circular dependency
        from astrocook.core.session import SessionV2 

        if not other_sessions:
            logging.warning("No sessions provided to stitch.")
            return self._session

        # 1. Resolve Inputs (GUI/Script Logic)
        resolved_specs = []
        
        for item in other_sessions:
            if isinstance(item, SessionV2):
                resolved_specs.append(item.spec)
            elif isinstance(item, str):
                # GUI Lookup
                found = None
                if hasattr(self._session, '_gui') and hasattr(self._session._gui, 'session_histories'):
                    for hist in self._session._gui.session_histories:
                        if hist.display_name == item:
                            found = hist.current_state
                            break
                if found:
                    resolved_specs.append(found.spec)
                else:
                    logging.warning(f"Stitch: Could not find loaded session named '{item}'. Skipping.")
            else:
                logging.warning(f"Stitch: Invalid item type {type(item)}. Skipping.")

        if not resolved_specs:
            logging.error("No valid sessions found to stitch.")
            return self._session

        try:
            # 2. Delegate Data Manipulation to Core
            new_spec = self._session.spec.stitch(resolved_specs, sort=sort)

            # Remove NaNs which break the 'step' plot (causing "scatter" look)
            # We use the public split API to keep valid data only.
            try:
                # Check for NaNs in X or Y
                x_val = new_spec.x.value
                y_val = new_spec.y.value
                if not (np.all(np.isfinite(x_val)) and np.all(np.isfinite(y_val))):
                    logging.info("Stitch: Cleaning up NaNs to fix plot display...")
                    # Filter: Keep rows where X and Y are finite
                    # Note: We rely on the core's split expression evaluator
                    new_spec = new_spec.split("isfinite(x) & isfinite(y)")
            except Exception as clean_err:
                logging.warning(f"Stitch cleanup warning: {clean_err}")
            
            # 3. Return New Session
            return self._session.with_new_spectrum(new_spec)
            
        except Exception as e:
            logging.error(f"Stitch failed: {e}", exc_info=True)
            return 0
        
    def coadd(self, session_names: str, xstart: str = 'None', xend: str = 'None', 
              dx: str = '0.01', xunit: str = 'nm', kappa: str = '5.0', 
              equalize_order: str = '0') -> Dict[str, Any]:
        """
        Co-add multiple sessions.

        Stitches multiple sessions together and rebins them onto a single common
        grid in one pass (Inverse Variance Weighting).
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.coadd`.

        Parameters
        ----------
        session_names : str
            Comma-separated names of sessions to co-add.
        xstart : str, optional
            Start wavelength. ``'None'`` for auto. Defaults to ``'None'``.
        xend : str, optional
            End wavelength. ``'None'`` for auto. Defaults to ``'None'``.
        dx : str, optional
            Step size for the final grid. Defaults to ``'0.01'``.
        xunit : str, optional
            Unit for the step size (e.g., ``'nm'`` or ``'km/s'``). Defaults to ``'nm'``.
        kappa : str, optional
            Sigma clipping threshold. ``'None'`` for off. Defaults to ``'5.0'``.
        equalize_order : str, optional
            Polynomial order for flux scaling (-1=Off, 0=Scalar, 1=Linear).
            Defaults to ``'0'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` containing the co-added spectrum,
            or 0 on failure.
        """
        
        names_list = [n.strip() for n in session_names.split(',') if n.strip()]
        if not names_list:
            logging.warning("Coadd: No sessions selected.")
            return 0

        # 1. Resolve Session Objects
        # We assume the first session in the list is the "primary" (self is implied context, 
        # but we need to find the specific objects named in the list).
        
        resolved_specs = []
        if hasattr(self._session, '_gui'):
            for name in names_list:
                for hist in self._session._gui.session_histories:
                    if hist.display_name == name:
                        resolved_specs.append(hist.current_state.spec)
                        break
        
        if not resolved_specs:
            logging.error("Coadd: No valid sessions found.")
            return 0

        # 2. Parse Parameters
        try:
            eq_order = int(equalize_order)
            xstart_q = au.Quantity(float(xstart), au.nm) if xstart != 'None' else None
            xend_q = au.Quantity(float(xend), au.nm) if xend != 'None' else None
            dx_unit = au.Unit(xunit) 
            dx_q = au.Quantity(float(dx), dx_unit) 
            kappa_f = float(kappa) if kappa != 'None' else None

            # 3. Delegate to Core
            # We use the first spectrum in the list as the 'base' to call the method on.
            base_spec = resolved_specs[0]
            others = resolved_specs[1:]
            
            logging.info(f"Co-adding {len(resolved_specs)} sessions onto {dx_q} grid...")
            
            new_spec = base_spec.coadd(
                others=others,
                xstart=xstart_q,
                xend=xend_q,
                dx=dx_q,
                kappa=kappa_f,
                equalize_order=eq_order
            )
            
            # 4. Return New Session
            return self._session.with_new_spectrum(new_spec)

        except Exception as e:
            logging.error(f"Coadd failed: {e}", exc_info=True)
            return 0