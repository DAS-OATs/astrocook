from astropy import units as au
from copy import deepcopy
import dataclasses
import logging
import numpy as np
from typing import TYPE_CHECKING, Optional, Union

from astrocook.core.photometry import calculate_synthetic_ab_mag, STANDARD_FILTERS
from astrocook.core.spectrum import DataColumnV2, SpectrumDataV2, SpectrumV2
from astrocook.core.spectrum_operations import rebin_spectrum
from astrocook.core.structures import HistoryLogV2
from astrocook.legacy.message import msg_param_fail

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2 

FLUX_RECIPES_SCHEMAS = {
    "calculate_running_std": {
        "brief": "Calculate running standard deviation.",
        "details": "Calculate a running Root-Mean-Square on a column (e.g., 'y') to estimate local error.",
        "params": [
            {"name": "input_col", "type": str, "default": "y", "doc": "Column to use for StdDev calculation (e.g., y)"},
            {"name": "output_col", "type": str, "default": "running_std", "doc": "Name of the new column to create (e.g., running_std)"},
            {"name": "window_pix", "type": int, "default": 21, "doc": "Total window size in pixels (should be odd)"}
        ],
        "url": "edit_cb.html#calculate_running_std" # Placeholder URL
    },
    "smooth": {
        "brief": "Smooth spectrum.",
        "details": "Apply a Gaussian filter to the spectrum flux (y-axis). The error (dy-axis) is not changed.",
        "params": [
            {"name": "sigma_kms", "type": float, "default": 100.0, "doc": "Standard deviation for Gaussian kernel (km/s)"},
        ],
        "url": "edit_cb.html#smooth" # Placeholder URL
    },
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
    },
    "resample": {
        "brief": "Resample on another grid.",
        "details": "Resample this spectrum onto the exact wavelength grid of another open session. This is the first step for multi-session arithmetics.",
        "params": [
            {"name": "target_session", "type": str, "default": "None", "doc": "Name of the session to use as the grid"}
        ],
        "url": "edit_cb.html#resample" # Placeholder URL
    },
    "calibrate_from_magnitudes": {
        "brief": "Flux calibrate (multi-band).",
        "details": "Warp the spectrum SED to match multiple photometric points. Enter as 'Filter=Mag' pairs separated by commas.",
        "params": [
            {"name": "magnitudes", "type": str, "default": "SDSS_r=17.0", "doc": "List of 'Filter=Mag' pairs (e.g., 'SDSS_g=17.5, SDSS_r=17.0')"}
        ],
        "url": "edit_cb.html#calibrate"
    }
}


class RecipeFluxV2:
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cb'

    def calculate_running_std(self, 
                              input_col: str = 'y', 
                              output_col: str = 'running_std', 
                              window_pix: str = '21') -> 'SessionV2':
        """
        API: Calculates a running StdDev on a column.
        """
        try:
            window_pix_i = int(window_pix)
            if window_pix_i < 3:
                window_pix_i = 3
            if window_pix_i % 2 == 0: # Ensure it's odd
                window_pix_i += 1
                
        except ValueError:
            logging.error(msg_param_fail)
            return 0
            
        try:
            # 1. Call the immutable V2 operation
            new_spec_v2 = self._session.spec.calculate_running_std(
                input_col=input_col,
                output_col=output_col,
                window_pix=window_pix_i
            )
            
            # 2. Return a NEW SessionV2 instance
            return self._session.with_new_spectrum(new_spec_v2)
            
        except Exception as e:
            logging.error(f"Failed during calculate_running_std: {e}", exc_info=True)
            return 0

    def smooth(self, sigma_kms: str = '100.0') -> 'SessionV2':
        """
        API: Applies Gaussian smoothing to the spectrum.
        """
        try:
            sigma_kms_f = float(sigma_kms)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        try:
            new_spec_v2 = self._session.spec.smooth(
                sigma_kms=sigma_kms_f
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during smooth: {e}", exc_info=True)
            return 0

    def rebin(self, xstart: str = 'None', xend: str = 'None', dx: str = '10.0', xunit: str = 'km/s',
              kappa: str = '5.0', filling: str = 'nan', norm: str = 'False') -> 'SessionV2':
        """
        API: Rebins the spectrum, logs the action, and returns a NEW SpectrumV2 instance.
        """

        # 4. Perform Type Casting
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

        # 5. Execute the immutable V2 operation
        new_spec_v2 = self._session.spec.rebin(
            xstart=xstart_q, xend=xend_q, dx=dx_q, kappa=kappa_f, filling=filling_f
        )
        
        # 6. Return a NEW SessionV2 instance (the new state)
        return self._session.with_new_spectrum(new_spec_v2)
    
    def resample(self, target_session: str = 'None') -> 'SessionV2':
        """
        API: Resamples the spectrum onto the grid of a target session.
        """
        if target_session == 'None' or target_session == '':
            logging.error("No target session was selected for resampling.")
            return 0

        try:
            # 1. Find the target session object from the main window
            target_hist = None
            if hasattr(self._session, '_gui') and hasattr(self._session._gui, 'session_histories'):
                for hist in self._session._gui.session_histories:
                    if hist.display_name == target_session:
                        target_hist = hist
                        break
            
            if target_hist is None:
                raise ValueError(f"Could not find an open session named '{target_session}'.")

            target_spec = target_hist.current_state.spec
            if target_spec is None:
                 raise ValueError(f"Target session '{target_session}' has no spectrum data.")

            # 2. Execute the immutable V2 operation
            new_spec_v2 = self._session.spec.resample_on_grid(
                target_grid_spec=target_spec
            )
            
            # 3. Return a NEW SessionV2 instance
            return self._session.with_new_spectrum(new_spec_v2)
            
        except Exception as e:
            logging.error(f"Failed during resample: {e}", exc_info=True)
            return 0
        
    def calibrate_from_magnitudes(self, magnitudes: str = "SDSS_r=17.0") -> 'SessionV2':
        """
        API: Rescales (and warps) the spectrum to match target magnitudes.
        """
        try:
            # The API layer now handles all the complexity
            new_spec_v2 = self._session.spec.calibrate(magnitudes)
            return self._session.with_new_spectrum(new_spec_v2)
            
        except Exception as e:
            # Standard error logging for the worker
            logging.error(f"Failed during calibrate_from_magnitudes: {e}", exc_info=True)
            return 0