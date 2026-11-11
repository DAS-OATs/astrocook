from astropy import units as au
import logging
import numpy as np
from typing import TYPE_CHECKING, Optional

from astrocook.recipes.errors import InvalidParameter, RecipeFailed
from astrocook.recipes.messages import msg_session_not_found
from astrocook.recipes._spectrum import spec_of
from astrocook.recipes.params import typed_recipe

if TYPE_CHECKING:
    from astrocook.session.session import SessionV2 

logger = logging.getLogger(__name__)

# --- The text of the methods (phase 6-ter, PR 5c: astrocook.methods_text) ---
from astrocook.methods_text import fmt
def _methods_rebin(p):
    from astrocook.recipes.edit import _clipping, _grid
    filling = "" if fmt(p.raw('filling')) == "nan" else f"; the pixels without data were set to {p['filling']}"
    return (f"The spectrum was rebinned onto {_grid(p)}, each new pixel being the "
            f"inverse-variance weighted mean of the old ones{_clipping(p)}{filling}{p.menu}")


def _methods_equalize(p):
    order = p.raw('order')
    how = ("a constant factor" if order == 0 else "a spline" if order == -1
           else f"a polynomial of order {order}")
    return (f"The flux was scaled to match that of the session {p['reference_session']} "
            f"where they overlap, by {how}{p.menu}")


FLUX_RECIPES_SCHEMAS = {
    "calculate_running_std": {
        "methods": "The standard deviation of the column {input_col} was computed in a running window of {window_pix} pixels and written in the column {output_col}{menu}",
        "brief": "Calculate running standard deviation.",
        "details": "Computes the standard deviation of a column in a running window, e.g. to estimate the local noise of the flux; written in a new column.",
        "params": [
            {"name": "input_col", "type": str, "default": "y", "doc": "Column of which the standard deviation is computed, e.g. 'y' (the flux).", "kind": "column"},
            {"name": "output_col", "type": str, "default": "running_std", "doc": "Name of the column to write."},
            {"name": "window_pix", "type": int, "default": 21, "gt": 0, "doc": "Width of the window (pixels; an odd number)."}
        ],
        "url": "edit_cb.html#calculate_running_std" # Placeholder URL
    },
    "smooth": {
        "methods": "The flux, its error and the columns like the flux were smoothed with a Gaussian kernel of standard deviation {sigma_kms} km/s{menu}",
        "brief": "Smooth spectrum.",
        "details": "Apply a Gaussian filter to the flux, its error and the columns like the flux (continuum, model, telluric transmission).",
        "params": [
            {"name": "sigma_kms", "type": float, "default": 100.0, "gt": 0, "doc": "Standard deviation of the Gaussian kernel (km/s)."},
        ],
        "url": "edit_cb.html#smooth" # Placeholder URL
    },
    "rebin": {
        "methods": _methods_rebin,
        "brief": "Rebin spectrum.",
        "details": "Rebins the spectrum onto a new grid with a constant step, in velocity (uniform in log wavelength) or in wavelength; the flux of each new pixel is the inverse-variance weighted mean of the old ones.",
        "params": [
            # None (written 'None' in the dialog and the log): the data range
            {"name": "xstart", "type": float, "default": None, "optional": True, "doc": "Start of the new grid (nm); None: the start of the data."},
            {"name": "xend", "type": float, "default": None, "optional": True, "doc": "End of the new grid (nm); None: the end of the data."},
            
            # CRITICAL: dx needs to be a number (10.0) but the V1 flow expects km/s unit conversion later.
            # We use the raw value here, and the adapter handles the unit.
            {"name": "dx", "type": float, "default": 10.0, "gt": 0, "doc": "Step of the new grid, in the unit of xunit."}, 
            
            # Use string "km/s" for the unit input field
            {"name": "xunit", "type": str, "default": "km/s", "doc": "Unit of the step: a velocity (km/s) or a wavelength (nm)."},
            
            # None (written 'None'): no clipping
            {"name": "kappa", "type": float, "default": 5.0, "optional": True, "gt": 0, "doc": "Threshold (in sigma) of the clipping of outliers; None: no clipping."},
            
            # 'norm' (V1: return the normalized spectrum) is not offered: it
            # was never applied in V2 (findings register F37, #560)

            # NaN (written 'nan') is the default
            {"name": "filling", "type": float, "default": float("nan"), "allow_nan": True, "doc": "Value of the new pixels without data (NaN by default)."},
        ],
        "url": "edit_cb.html#rebin"
    },
    "equalize": {
        "methods": _methods_equalize,
        "brief": "Equalize flux to reference.",
        "details": "Scales the flux of this session to match that of another session where they overlap.",
        "params": [
            {"name": "reference_session", "type": str, "default": "None", "doc": "Name of the open session whose flux is matched.", "kind": "session"},
            {"name": "order", "type": int, "default": 0, "doc": "Model of the scaling: 0 a constant factor, n > 0 a polynomial of order n, -1 a spline."}
        ]
    },
    "resample": {
        "methods": "The spectrum was interpolated onto the wavelength grid of the session {target_session}, into a new session{menu}",
        "brief": "Resample on another grid.",
        "branching": True,
        "details": "Interpolates this spectrum onto the wavelength grid of another open session, into a new session; the first step to combine two sessions pixel by pixel.",
        "params": [
            {"name": "target_session", "type": str, "default": "None", "doc": "Name of the open session whose grid is used.", "kind": "session"}
        ],
        "url": "edit_cb.html#resample" # Placeholder URL
    },
    "calibrate_from_magnitudes": {
        "methods": "The flux was scaled, smoothly with wavelength, to match the magnitudes {magnitudes}{menu}",
        "brief": "Flux calibrate (multi-band).",
        "details": "Scales the flux, smoothly with wavelength, to match the magnitudes of the object in one or more photometric filters.",
        "params": [
            {"name": "magnitudes", "type": str, "default": "SDSS_r=17.0", "doc": "Comma-separated 'filter=magnitude' pairs, e.g. 'SDSS_g=17.5, SDSS_r=17.0'."}
        ],
        "url": "edit_cb.html#calibrate"
    }
}


def _find_session(session: 'SessionV2', name: str):
    """The history of the open session called ``name`` (RecipeFailed if
    there is none, or no GUI to find it in)."""
    gui = getattr(session, '_gui', None)
    for hist in getattr(gui, 'session_histories', None) or []:
        if hist.display_name == name:
            return hist
    raise RecipeFailed(msg_session_not_found.format(name=name))


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
        SessionV2
            A new :class:`~astrocook.session.session.SessionV2` with the calculated column added.
        """
        window_pix_i = max(window_pix, 3)
        if window_pix_i % 2 == 0: # Ensure it's odd
            window_pix_i += 1

        new_spec_v2 = spec_of(self._session).calculate_running_std(
            input_col=input_col,
            output_col=output_col,
            window_pix=window_pix_i
        )
        return self._session.with_new_spectrum(new_spec_v2)

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
        SessionV2
            A new :class:`~astrocook.session.session.SessionV2` with smoothed columns.
        """
        new_spec_v2 = spec_of(self._session).smooth(sigma_kms=sigma_kms)
        return self._session.with_new_spectrum(new_spec_v2)

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
        SessionV2
            A new :class:`~astrocook.session.session.SessionV2` with the rebinned spectrum.
        """

        xstart_q = au.Quantity(xstart, au.nm) if xstart is not None else None
        xend_q = au.Quantity(xend, au.nm) if xend is not None else None
        try:
            # V1 logic from cookbook_flux.py uses the unit attached to xunit for dx
            dx_q = au.Quantity(dx, au.Unit(xunit))
        except (ValueError, TypeError):
            # (Findings register F23: it crashed the GUI)
            raise InvalidParameter(f"xunit: {xunit!r} is not a unit") from None

        # The V1 normalization was never applied in V2 and is not offered
        # any more (F37, #560): say so instead of ignoring it silently
        if str(norm) == 'True':
            logger.warning("The parameter 'norm' is ignored: the rebinned spectrum is "
                           "not normalized (use the normalized view of the "
                           "plot, or apply_expression with y/cont).")

        # 5. Execute the immutable V2 operation
        new_spec_v2 = spec_of(self._session).rebin(
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
            A new session with the scaled flux. 
        """
        if reference_session in ['None', '']:
            raise RecipeFailed("No reference session is selected.")
        ref_hist = _find_session(self._session, reference_session)
        new_spec = spec_of(self._session).equalize_to_reference(
            ref_hist.current_state.spec, order=order
        )
        return self._session.with_new_spectrum(new_spec)
    
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
        SessionV2
            A new :class:`~astrocook.session.session.SessionV2` interpolated to the target grid.
        """
        if target_session == 'None' or target_session == '':
            raise RecipeFailed("No target session is selected.")
        target_spec = _find_session(self._session, target_session).current_state.spec
        if target_spec is None:
            raise RecipeFailed(f"The session '{target_session}' has no spectrum.")
        new_spec_v2 = spec_of(self._session).resample_on_grid(target_grid_spec=target_spec)
        return self._session.with_new_spectrum(new_spec_v2)
        
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
        SessionV2
            A new :class:`~astrocook.session.session.SessionV2` with calibrated flux.
        """
        new_spec_v2 = spec_of(self._session).calibrate(magnitudes)
        return self._session.with_new_spectrum(new_spec_v2)