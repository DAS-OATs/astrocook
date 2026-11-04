from astropy import units as au
import logging
import numpy as np
from typing import TYPE_CHECKING, Optional

from astrocook.recipes.messages import msg_param_fail, msg_session_not_found
from astrocook.recipes.params import typed_recipe

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2 

logger = logging.getLogger(__name__)

FLUX_RECIPES_SCHEMAS = {
    "calculate_running_std": {
        "brief": "Calculate running standard deviation.",
        "details": "Calculate a running Root-Mean-Square on a column (e.g., 'y') to estimate local error.",
        "params": [
            {"name": "input_col", "type": str, "default": "y", "doc": "Column to use for StdDev calculation (e.g., y)"},
            {"name": "output_col", "type": str, "default": "running_std", "doc": "Name of the new column to create (e.g., running_std)"},
            {"name": "window_pix", "type": int, "default": 21, "gt": 0, "doc": "Total window size in pixels (should be odd)"}
        ],
        "url": "edit_cb.html#calculate_running_std" # Placeholder URL
    },
    "smooth": {
        "brief": "Smooth spectrum.",
        "details": "Apply a Gaussian filter to the spectrum flux (y-axis). The error (dy-axis) is not changed.",
        "params": [
            {"name": "sigma_kms", "type": float, "default": 100.0, "gt": 0, "doc": "Standard deviation for Gaussian kernel (km/s)"},
        ],
        "url": "edit_cb.html#smooth" # Placeholder URL
    },
    "rebin": {
        "brief": "Rebin spectrum.",
        "details": "Apply a new binning to a spectrum, with a constant bin size.",
        "params": [
            # None (written 'None' in the dialog and the log): the data range
            {"name": "xstart", "type": float, "default": None, "optional": True, "doc": "Start wavelength (nm, None for auto)"},
            {"name": "xend", "type": float, "default": None, "optional": True, "doc": "End wavelength (nm, None for auto)"},
            
            # CRITICAL: dx needs to be a number (10.0) but the V1 flow expects km/s unit conversion later.
            # We use the raw value here, and the adapter handles the unit.
            {"name": "dx", "type": float, "default": 10.0, "gt": 0, "doc": "Step in x"}, 
            
            # Use string "km/s" for the unit input field
            {"name": "xunit", "type": str, "default": "km/s", "doc": "Unit of wavelength or velocity"},
            
            # None (written 'None'): no clipping
            {"name": "kappa", "type": float, "default": 5.0, "optional": True, "gt": 0, "doc": "Sigma clipping threshold (None for off)"},
            
            # 'norm' (V1: return the normalized spectrum) is not offered: it
            # was never applied in V2 (findings register F37, #560)

            # NaN (written 'nan') is the default
            {"name": "filling", "type": float, "default": float("nan"), "allow_nan": True, "doc": "Value to fill region without data"},
        ],
        "url": "edit_cb.html#rebin"
    },
    "equalize": {
        "brief": "Equalize flux to reference.",
        "details": "Scale the flux of this session to match another 'reference' session.",
        "params": [
            {"name": "reference_session", "type": str, "default": "None", "doc": "Name of the reference session"},
            {"name": "order", "type": int, "default": 0, "doc": "Polynomial order (-1=spline, 0=scalar)"}
        ]
    },
    "resample": {
        "brief": "Resample on another grid.",
        "branching": True,
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
    """
    Recipes for manipulating spectral flux and grid sampling.

    These methods are accessed via the ``session.flux`` attribute and delegate
    logic to :class:`~astrocook.core.spectrum.SpectrumV2`.
    """
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cb'

    @typed_recipe
    def calculate_running_std(self,
                              input_col: str = 'y',
                              output_col: str = 'running_std',
                              window_pix: int = 21) -> 'SessionV2':
        """
        Calculate running standard deviation.

        Calculates a running Root-Mean-Square (RMS) on a specific column to estimate
        local error or variability.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.calculate_running_std`.

        Parameters
        ----------
        input_col : str, optional
            The column to use for StdDev calculation (e.g., ``'y'``). Defaults to ``'y'``.
        output_col : str, optional
            Name of the new column to create (e.g., ``'running_std'``).
            Defaults to ``'running_std'``.
        window_pix : int, optional
            Total window size in pixels (should be odd: an even one is
            increased by 1, and one below 3 is 3). Defaults to ``21``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with the calculated column added,
            or 0 on failure.
        """
        window_pix_i = max(window_pix, 3)
        if window_pix_i % 2 == 0: # Ensure it's odd
            window_pix_i += 1

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
            logger.error(f"Could not compute the running standard deviation: {e}", exc_info=True)
            return 0

    @typed_recipe
    def smooth(self, sigma_kms: float = 100.0) -> 'SessionV2':
        """
        Smooth spectrum.

        Applies a Gaussian filter to the spectrum flux (y-axis), error (dy-axis),
        and auxiliary columns.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.smooth`.

        Parameters
        ----------
        sigma_kms : float, optional
            Standard deviation for Gaussian kernel in km/s. Defaults to ``100.0``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with smoothed columns,
            or 0 on failure.
        """
        try:
            new_spec_v2 = self._session.spec.smooth(
                sigma_kms=sigma_kms
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logger.error(f"Could not smooth the spectrum: {e}", exc_info=True)
            return 0

    @typed_recipe
    def rebin(self, xstart: Optional[float] = None, xend: Optional[float] = None,
              dx: float = 10.0, xunit: str = 'km/s', kappa: Optional[float] = 5.0,
              filling: float = np.nan, norm: str = 'False') -> 'SessionV2':
        """
        Rebin spectrum.

        Applies a new binning to a spectrum with a constant bin size.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.rebin`.

        Parameters
        ----------
        xstart : float or None, optional
            Start wavelength (nm). None for auto (data minimum).
            Defaults to None.
        xend : float or None, optional
            End wavelength (nm). None for auto (data maximum).
            Defaults to None.
        dx : float, optional
            Step size in x. Defaults to ``10.0``.
        xunit : str, optional
            Unit of wavelength or velocity for the step size (e.g., ``'km/s'``, ``'nm'``).
            Defaults to ``'km/s'``.
        kappa : float or None, optional
            Sigma clipping threshold for binning. None for off. Defaults to ``5.0``.
        filling : float, optional
            Value to fill regions without data. Defaults to NaN.
        norm : str, optional
            Ignored, kept so that scripts and logs written for V1 still run:
            the rebinned spectrum is never normalized, and ``'True'`` only
            logs a warning. To see the normalized flux, use the normalized
            view of the plot, or ``apply_expression`` with ``y/cont``.
            Defaults to ``'False'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with the rebinned spectrum,
            or 0 on failure.
        """

        xstart_q = au.Quantity(xstart, au.nm) if xstart is not None else None
        xend_q = au.Quantity(xend, au.nm) if xend is not None else None
        try:
            # V1 logic from cookbook_flux.py uses the unit attached to xunit for dx
            dx_q = au.Quantity(dx, au.Unit(xunit))
        except (ValueError, TypeError):
            # Return 0 like the other recipes (findings register F23)
            logger.error(f"{msg_param_fail} xunit: {xunit!r} is not a unit.")
            return 0

        # The V1 normalization was never applied in V2 and is not offered
        # any more (F37, #560): say so instead of ignoring it silently
        if str(norm) == 'True':
            logger.warning("The parameter 'norm' is ignored: the rebinned spectrum is "
                           "not normalized (use the normalized view of the "
                           "plot, or apply_expression with y/cont).")

        # 5. Execute the immutable V2 operation
        new_spec_v2 = self._session.spec.rebin(
            xstart=xstart_q, xend=xend_q, dx=dx_q, kappa=kappa, filling=filling
        )
        
        # 6. Return a NEW SessionV2 instance (the new state)
        return self._session.with_new_spectrum(new_spec_v2)
    
    @typed_recipe
    def equalize(self, reference_session: str = 'None', order: int = 0) -> 'SessionV2':
        """
        Equalize flux to a reference session.

        Scales the flux of the current session to match the flux level of a 
        reference session using a polynomial scaling model.

        Parameters
        ----------
        reference_session : str, optional
            The name of the reference session. Defaults to ``'None'``.
        order : int, optional
            The polynomial order for the scaling model (``0`` = scalar, ``-1`` = spline).
            Defaults to ``0``.

        Returns
        -------
        SessionV2
            A new session with the scaled flux. Returns ``0`` (int) on failure.
        """
        if reference_session in ['None', '']:
            logger.error("No reference session is selected.")
            return 0
        try:
            # Find the reference session object via the GUI bridge
            ref_hist = next((h for h in self._session._gui.session_histories 
                            if h.display_name == reference_session), None)
            
            if not ref_hist:
                raise ValueError(msg_session_not_found.format(name=reference_session))

            new_spec = self._session.spec.equalize_to_reference(
                ref_hist.current_state.spec, order=order
            )
            return self._session.with_new_spectrum(new_spec)
        except Exception as e:
            logger.error(f"Could not equalize the flux: {e}", exc_info=True)
            return 0
    
    @typed_recipe
    def resample(self, target_session: str = 'None') -> 'SessionV2':
        """
        Resample on another grid.

        Resamples this spectrum onto the exact wavelength grid of another open session.
        This is often a prerequisite for multi-session arithmetic.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.resample_on_grid`.

        Parameters
        ----------
        target_session : str
            Name of the session to use as the grid source.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` interpolated to the target grid,
            or 0 on failure.
        """
        if target_session == 'None' or target_session == '':
            logger.error("No target session is selected.")
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
            logger.error(f"Could not resample the spectrum: {e}", exc_info=True)
            return 0
        
    @typed_recipe
    def calibrate_from_magnitudes(self, magnitudes: str = "SDSS_r=17.0") -> 'SessionV2':
        """
        Flux calibrate (multi-band).

        Warps the spectrum SED to match multiple photometric points.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.calibrate`.

        Parameters
        ----------
        magnitudes : str, optional
            List of 'Filter=Mag' pairs separated by commas (e.g., ``"SDSS_g=17.5, SDSS_r=17.0"``).
            Defaults to ``"SDSS_r=17.0"``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with calibrated flux,
            or 0 on failure.
        """
        try:
            # The API layer now handles all the complexity
            new_spec_v2 = self._session.spec.calibrate(magnitudes)
            return self._session.with_new_spectrum(new_spec_v2)
            
        except Exception as e:
            # Standard error logging for the worker
            logger.error(f"Could not calibrate the flux: {e}", exc_info=True)
            return 0