from astrocook._renamed import renamed
from astropy import units as au
import logging
import numpy as np
from typing import TYPE_CHECKING, Optional, Union

from astrocook.core.column_names import CORE_COLUMNS as CORE_COLUMN_NAMES, internal_column, is_core_column
from astrocook.core.photometry import STANDARD_FILTERS
from astrocook.recipes.errors import InvalidParameter, RecipeFailed
from astrocook.recipes._combine import check_axes, flux_relation
from astrocook.recipes._sessions import is_session, open_session
from astrocook.recipes._spectrum import spec_of
from astrocook.recipes.params import typed_recipe

if TYPE_CHECKING:
    from astrocook.session.session import Session 

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


def _window_used(window_pix) -> int:
    """The window that :meth:`RecipeFlux.calculate_running_std` uses: an
    odd number of pixels, at least 3."""
    window = max(int(window_pix), 3)
    return window + 1 if window % 2 == 0 else window


def _methods_running_std(p):
    # The window used, not the one passed (phase 6-quater, PR 7: an even
    # window passed was written as it was)
    try:
        window = str(_window_used(p.raw('window_pix')))
    except (TypeError, ValueError):
        window = p['window_pix']  # a log written by hand: as written
    return (f"The standard deviation of the column {p['input_col']} was computed in a "
            f"running window of {window} pixels and written "
            f"in the column {p['output_col']}{p.menu}")


#: The columns of every spectrum: not a name for a column written by a recipe
#: (phase 6-quater, PR 7: the result went to a column that nothing shows),
#: by the names of the data; those of the GUI (λ, F, dF, ...) too since PR 8
CORE_COLUMNS = CORE_COLUMN_NAMES


FLUX_RECIPES_SCHEMAS = {
    "calculate_running_std": {
        "methods": _methods_running_std,
        "brief": "Calculate running standard deviation.",
        "url": "tutorials/flux_operations.html#calculating-running-statistics",
        "menu": "Flux > Calculate Running &StdDev",
        "details": "Computes the standard deviation of a column in a running window, e.g. to estimate the local noise of the flux; written in a new column.",
        "params": [
            {"name": "input_col", "type": str, "default": "F", "doc": "Column of which the standard deviation is computed, e.g. 'F' (the flux).", "kind": "column"},
            {"name": "output_col", "type": str, "default": "running_std", "doc": "Name of the column to write (a new one, or one to replace; not λ, λmin, λmax, F, dF)."},
            {"name": "window_pix", "type": int, "default": 21, "gt": 0, "doc": "Width of the window (pixels): an odd number; an even one is increased by 1, and one below 3 is 3."}
        ],
    },
    "smooth": {
        "methods": "The flux, its error and the columns like the flux were smoothed with a Gaussian kernel of standard deviation {sigma_kms} km/s{menu}",
        "brief": "Smooth spectrum.",
        "url": "tutorials/flux_operations.html#smoothing",
        "menu": "Flux > &Smooth Spectrum",
        "details": "Apply a Gaussian filter to the flux, its error and the columns like the flux (continuum, model, telluric transmission).",
        "params": [
            {"name": "sigma_kms", "type": float, "default": 100.0, "gt": 0, "doc": "Standard deviation of the Gaussian kernel (km/s)."},
        ],
    },
    "rebin": {
        "methods": _methods_rebin,
        "brief": "Rebin spectrum.",
        "url": "tutorials/flux_operations.html#rebinning",
        "menu": "Flux > &Rebin Spectrum",
        "details": "Rebins the spectrum onto a new grid with a constant step, in velocity (uniform in log wavelength) or in wavelength; the flux of each new pixel is the inverse-variance weighted mean of the old ones.",
        "params": [
            # None (written 'None' in the dialog and the log): the data range
            {"name": "xstart", "type": float, "default": None, "optional": True, "doc": "Start of the new grid (nm); None: the start of the data."},
            {"name": "xend", "type": float, "default": None, "optional": True, "doc": "End of the new grid (nm), above xstart; None: the end of the data."},
            
            # CRITICAL: dx needs to be a number (10.0) but the V1 flow expects km/s unit conversion later.
            # We use the raw value here, and the adapter handles the unit.
            {"name": "dx", "type": float, "default": 10.0, "gt": 0, "doc": "Step of the new grid, in the unit of xunit."}, 
            
            # Use string "km/s" for the unit input field
            {"name": "xunit", "type": str, "default": "km/s", "doc": "Unit of the step: a velocity (km/s, a grid uniform in log wavelength) or a wavelength (nm)."},
            
            # None (written 'None'): no clipping
            {"name": "kappa", "type": float, "default": 5.0, "optional": True, "gt": 0, "doc": "Threshold (in sigma) of the clipping of outliers; None: no clipping."},
            
            # 'norm' (V1: return the normalized spectrum) is not offered: it
            # was never applied in V2 (findings register F37, #560); the
            # method still takes it, for the scripts and logs of V1

            # NaN (written 'nan') is the default
            {"name": "filling", "type": float, "default": float("nan"), "allow_nan": True, "doc": "Value of the new pixels without data (NaN by default)."},
        ],
    },
    "equalize": {
        "methods": _methods_equalize,
        "brief": "Equalize flux to reference.",
        "url": "tutorials/flux_operations.html#equalizing",
        "menu": "Flux > &Equalize to Reference",
        "details": "Scales the flux of this session to match that of another session where they overlap, compared in the units of this session; a flux in arbitrary units takes the unit of the reference once scaled, a normalized flux is not scaled to one that is not.",
        "params": [
            # Required (phase 6-quater, PR 7: the default was the string 'None',
            # which still means "no session" and fails as before)
            {"name": "reference_session", "type": str, "default": "", "doc": "Name of the open session whose flux is matched.", "kind": "session"},
            # -1 is the lowest (phase 6-quater, PR 7: -2, -3... were a spline
            # too, without a word). Not the -1 of coadd, which means no
            # scaling (findings register F123)
            {"name": "order", "type": int, "default": 0, "min": -1, "doc": "Model of the scaling: 0 a constant factor, n > 0 a polynomial of order n, -1 a spline."}
        ]
    },
    "resample": {
        "methods": "The spectrum was interpolated onto the wavelength grid of the session {target_session}, into a new session{menu}",
        "brief": "Resample on another grid.",
        "url": "tutorials/flux_operations.html#resampling",
        "menu": "Flux > Re&sample on Grid",
        "branching": True,
        "details": "Interpolates this spectrum onto the wavelength grid of another open session, into a new session; the first step to combine two sessions pixel by pixel.",
        "params": [
            # Required, as reference_session of equalize (phase 6-quater, PR 7)
            {"name": "target_session", "type": str, "default": "", "doc": "Name of the open session whose grid is used.", "kind": "session"}
        ],
    },
    "calibrate_from_magnitudes": {
        "methods": "The flux was scaled, smoothly with wavelength, to match the magnitudes {magnitudes}{menu}",
        "brief": "Flux calibrate (multi-band).",
        "url": "tutorials/flux_operations.html#flux-calibration",
        "menu": "Flux > Flux &Calibrate",
        "details": "Scales the flux, smoothly with wavelength, to match the magnitudes of the object in one or more photometric filters.",
        "params": [
            {"name": "magnitudes", "type": str, "default": "SDSS_r=17.0", "doc": "Comma-separated 'filter=magnitude' pairs (AB magnitudes), e.g. 'SDSS_g=17.5, SDSS_r=17.0'. Filters: " + ", ".join(STANDARD_FILTERS) + "."}
        ],
    }
}


def _check_magnitudes(magnitudes: str) -> None:
    """:class:`InvalidParameter` if ``magnitudes`` is not a list of
    ``filter=magnitude`` pairs of known filters; the same reading as
    :func:`astrocook.core.photometry.generate_calibration_curve`."""
    fmt_error = (f"magnitudes: {magnitudes!r} is not a list of 'filter=magnitude' "
                 f"pairs, e.g. 'SDSS_g=17.5, SDSS_r=17.0'")
    for pair in magnitudes.split(','):
        parts = pair.split('=')
        if len(parts) != 2:
            raise InvalidParameter(fmt_error)
        name, mag = parts[0].strip(), parts[1].strip()
        try:
            float(mag)
        except ValueError:
            raise InvalidParameter(fmt_error) from None
        if name not in STANDARD_FILTERS:
            raise InvalidParameter(f"magnitudes: the filter {name!r} is not known; the "
                                   f"filters are {', '.join(STANDARD_FILTERS)}")


def _find_session(session: 'Session', name, param: str) -> 'Session':
    """The open session called ``name``, or ``name`` itself in a script (the
    rule of :mod:`astrocook.recipes._sessions`, phase 6-quater, PR 9);
    RecipeFailed if there is none."""
    return open_session(session, name, param)


class RecipeFlux:
    """
    Recipes for manipulating spectral flux and grid sampling.

    These methods are accessed via the ``session.flux`` attribute and delegate
    logic to :class:`~astrocook.core.spectrum.Spectrum`.
    """
    def __init__(self, session_v2: 'Session'):
        self._session = session_v2
        self._tag = 'cb'

    @typed_recipe
    def calculate_running_std(self,
                              input_col: str = 'F',
                              output_col: str = 'running_std',
                              window_pix: int = 21) -> 'Session':
        """
        Calculate running standard deviation.

        Calculates a running Root-Mean-Square (RMS) on a specific column to estimate
        local error or variability.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.calculate_running_std`.

        Parameters
        ----------
        input_col : str, optional
            The column to use for StdDev calculation, by the names of the GUI
            (``'F'``, ``'dF'``, ``'cont'``, ...; ``'y'``, ``'dy'`` are
            accepted). Defaults to ``'F'`` (``'y'`` until phase 6-quater,
            PR 8: the same column).
        output_col : str, optional
            Name of the column to write (e.g., ``'running_std'``): a new
            column, or an auxiliary one to replace; not empty, nor one of
            the columns of every spectrum (``λ``, ``λmin``, ``λmax``, ``F``,
            ``dF``, or their names in the data, ``x``, ..., ``dy``).
            Defaults to ``'running_std'``.
        window_pix : int, optional
            Total window size in pixels (should be odd: an even one is
            increased by 1, and one below 3 is 3). Defaults to ``21``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with the calculated column added.

        Raises
        ------
        InvalidParameter
            If ``output_col`` is empty or the name of a column of every
            spectrum (phase 6-quater, PR 7: the result went to a hidden
            column, findings register F122; PR 8: also by the names of the
            GUI, which made a column that no expression could name).
        """
        if not output_col.strip() or is_core_column(output_col.strip()):
            raise InvalidParameter(f"output_col: {output_col!r} is not a name for a new "
                                   f"column (nor one of λ, λmin, λmax, F, dF)")
        window_pix_i = _window_used(window_pix)

        new_spec_v2 = spec_of(self._session).calculate_running_std(
            # The names of the GUI to those of the data (F120, PR 8: 'F',
            # offered by the dialog, was not found)
            input_col=internal_column(input_col.strip()),
            output_col=output_col,
            window_pix=window_pix_i
        )
        return self._session.with_new_spectrum(new_spec_v2)

    @typed_recipe
    def smooth(self, sigma_kms: float = 100.0) -> 'Session':
        """
        Smooth spectrum.

        Applies a Gaussian filter to the spectrum flux (y-axis), error (dy-axis),
        and auxiliary columns.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.smooth`.

        Parameters
        ----------
        sigma_kms : float, optional
            Standard deviation for Gaussian kernel in km/s. Defaults to ``100.0``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with smoothed columns.
        """
        new_spec_v2 = spec_of(self._session).smooth(sigma_kms=sigma_kms)
        return self._session.with_new_spectrum(new_spec_v2)

    @typed_recipe
    def rebin(self, xstart: Optional[float] = None, xend: Optional[float] = None,
              dx: float = 10.0, xunit: str = 'km/s', kappa: Optional[float] = 5.0,
              filling: float = np.nan, norm: Union[bool, str] = False) -> 'Session':
        """
        Rebin spectrum.

        Applies a new binning to a spectrum with a constant bin size.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.rebin`.

        Parameters
        ----------
        xstart : float or None, optional
            Start wavelength (nm). None for auto (data minimum).
            Defaults to None.
        xend : float or None, optional
            End wavelength (nm), above ``xstart``. None for auto (data
            maximum). Defaults to None.
        dx : float, optional
            Step size, in the unit of ``xunit``. Defaults to ``10.0``.
        xunit : str, optional
            Unit of the step: a velocity (``'km/s'``, a grid uniform in log
            wavelength) or a wavelength (``'nm'``). Defaults to ``'km/s'``.
        kappa : float or None, optional
            Sigma clipping threshold for binning. None for off. Defaults to ``5.0``.
        filling : float, optional
            Value to fill regions without data. Defaults to NaN.
        norm : bool or str, optional
            Ignored, kept so that scripts and logs written for V1 still run:
            the rebinned spectrum is never normalized, and a true value
            (``True``, ``'True'``, ``'true'``, ``1``) only logs a warning. To
            see the normalized flux, use the normalized view of the plot, or
            ``apply_expression`` with ``y/cont``. Defaults to ``False``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with the rebinned spectrum.

        Raises
        ------
        InvalidParameter
            If ``xunit`` is not a unit, or ``xstart`` is not below ``xend``
            (phase 6-quater, PR 7: an empty spectrum, findings register
            F124).
        """
        if xstart is not None and xend is not None and not xstart < xend:
            raise InvalidParameter(f"xstart: {xstart} is not below xend ({xend})")

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
        if str(norm).strip().lower() in ('true', '1'):
            logger.warning("The parameter 'norm' is ignored: the rebinned spectrum is "
                           "not normalized (use the normalized view of the "
                           "plot, or apply_expression with y/cont).")

        # 5. Execute the immutable V2 operation
        new_spec_v2 = spec_of(self._session).rebin(
            xstart=xstart_q, xend=xend_q, dx=dx_q, kappa=kappa, filling=filling
        )
        
        # 6. Return a NEW Session instance (the new state)
        return self._session.with_new_spectrum(new_spec_v2)
    
    @typed_recipe
    def equalize(self, reference_session: str = '', order: int = 0) -> 'Session':
        """
        Equalize flux to a reference session.

        Scales the flux of the current session to match the flux level of a 
        reference session using a polynomial scaling model.

        Parameters
        ----------
        reference_session : str
            The name of the open session whose flux is matched (required;
            ``''`` or ``'None'``, the default of the old logs, fail saying
            that no session is selected).
        order : int, optional
            The model of the scaling: ``0`` a constant factor, ``n > 0`` a
            polynomial of order ``n``, ``-1`` a spline (the lowest value).
            Defaults to ``0``.

        Returns
        -------
        Session
            A new session with the scaled flux, in the unit of its flux (an
            arbitrary flux scaled to a physical one takes the unit of the
            reference).

        Raises
        ------
        RecipeFailed
            If the reference is not an open session, its axis cannot be
            compared with this one (a velocity against a wavelength), or a
            flux is normalized and the other is not (phase 6-quater, PR 9,
            F138: the values were compared as they were, so an axis in Å
            against one in nm was not scaled, and a flux per Å against one
            per nm came out off by their ratio).
        """
        if reference_session in ['None', '']:
            raise RecipeFailed("No reference session is selected.")
        reference = _find_session(self._session, reference_session, 'reference_session')
        spec, ref = spec_of(self._session), spec_of(reference)
        ref_name = reference.name if is_session(reference_session) else reference_session
        check_axes(spec.x.unit, ref.x.unit, self._session.name, ref_name)
        how = flux_relation(spec.y.unit, ref.y.unit, self._session.name, ref_name, scaled=True)
        new_spec = spec.equalize_to_reference(ref, order=order)
        if how == 'declare':
            new_spec = new_spec.with_flux_unit(ref.y.unit)
        return self._session.with_new_spectrum(new_spec)
    
    @typed_recipe
    def resample(self, target_session: str = '') -> 'Session':
        """
        Resample on another grid.

        Resamples this spectrum onto the exact wavelength grid of another open session.
        This is often a prerequisite for multi-session arithmetic.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.resample_on_grid`.

        Parameters
        ----------
        target_session : str
            Name of the open session whose grid is used (required; ``''``
            or ``'None'``, the default of the old logs, fail saying that no
            session is selected).

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` interpolated to the target grid.
        """
        if target_session in ['None', '']:
            raise RecipeFailed("No target session is selected.")
        target_spec = _find_session(self._session, target_session, 'target_session').spec
        if target_spec is None:
            raise RecipeFailed(f"The session '{target_session}' has no spectrum.")
        new_spec_v2 = spec_of(self._session).resample_on_grid(target_grid_spec=target_spec)
        return self._session.with_new_spectrum(new_spec_v2)
        
    @typed_recipe
    def calibrate_from_magnitudes(self, magnitudes: str = "SDSS_r=17.0") -> 'Session':
        """
        Flux calibrate (multi-band).

        Warps the spectrum SED to match multiple photometric points.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.calibrate`.

        Parameters
        ----------
        magnitudes : str, optional
            List of 'Filter=Mag' pairs (AB magnitudes) separated by commas
            (e.g., ``"SDSS_g=17.5, SDSS_r=17.0"``), the filters among those
            of :data:`astrocook.core.photometry.STANDARD_FILTERS`. Defaults
            to ``"SDSS_r=17.0"``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with calibrated flux.

        Raises
        ------
        InvalidParameter
            If ``magnitudes`` is not a list of pairs, or names a filter that
            is not known (phase 6-quater, PR 7: a RecipeFailed from the
            core).
        """
        _check_magnitudes(magnitudes)
        new_spec_v2 = spec_of(self._session).calibrate(magnitudes)
        return self._session.with_new_spectrum(new_spec_v2)

# Old names (phase 6-quater, PR 6): see astrocook/_renamed.py
__getattr__ = renamed(__name__, {"RecipeFluxV2": "RecipeFlux"})
