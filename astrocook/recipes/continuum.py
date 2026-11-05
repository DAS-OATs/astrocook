import astropy.units as au
import logging
from typing import TYPE_CHECKING

from astrocook.recipes.errors import InvalidParameter, MissingPrerequisite
from astrocook.recipes.messages import msg_z_em_unset
from astrocook.recipes.params import typed_recipe

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2

logger = logging.getLogger(__name__)

# --- Schemas for the Continuum Menu ---
CONTINUUM_RECIPES_SCHEMAS = {
    "estimate_auto": {
        "brief": "Auto-estimate continuum (single-click).",
        "provides": ("continuum",),
        "requires": ("z_em",),
        "details": "Finds the absorbed pixels (as Find absorbed regions) and fits the continuum to the others (as Fit continuum to mask), in one step; writes 'abs_mask' and 'cont'.",
        "params": [
            {"name": "smooth_len_lya", "type": float, "default": 5000.0, "gt": 0, "doc": "Smoothing length in the Lyman-alpha forest (km/s): long, as most pixels there are absorbed."},
            {"name": "smooth_len_out", "type": float, "default": 400.0, "gt": 0, "doc": "Smoothing length outside the forest (km/s)."},
            {"name": "kappa", "type": float, "default": 2.0, "gt": 0, "doc": "Threshold (in sigma) below which a pixel is taken as absorbed."},
            {"name": "fudge", "type": str, "default": "auto", "doc": "Factor that multiplies the continuum: 'auto' to centre it on the unabsorbed pixels, or a number (e.g. '1.0')."}, 
            {"name": "smooth_fwhm", "type": float, "default": 800.0, "gt": 0, "doc": "FWHM of the Gaussian smoothing of the final continuum (km/s)."},
            {"name": "template", "type": bool, "default": False, "doc": "Trace the Lyman-alpha emission with a quasar template (needs z_em)."},
            {"name": "renorm_model", "type": bool, "default": True, "doc": "Rescale the model to the new continuum (asked by the GUI).", "gui_hidden": True, "agents": False},
        ],
        "url": "continuum_cb.html#estimate_auto"
    },
    "find_absorbed": {
        "brief": "Find absorbed regions (V1 logic).",
        "requires": ("z_em",),
        "details": "Finds the absorbed pixels with an iterative kappa-sigma clipping of the smoothed flux (the V1 'clip_flux' algorithm), with a longer smoothing in the Lyman-alpha forest; writes the 'abs_mask' column (True for the absorbed pixels).",
        "params": [
            {"name": "smooth_len_lya", "type": float, "default": 5000.0, "gt": 0, "doc": "Smoothing length in the Lyman-alpha forest (km/s): long, as most pixels there are absorbed."},
            {"name": "smooth_len_out", "type": float, "default": 400.0, "gt": 0, "doc": "Smoothing length outside the forest (km/s)."},
            {"name": "kappa", "type": float, "default": 2.0, "gt": 0, "doc": "Threshold (in sigma) below which a pixel is taken as absorbed."},
            {"name": "template", "type": bool, "default": False, "doc": "Not implemented: ignored."},
        ],
        "url": "continuum_cb.html#find_absorbed"
    },
    "fit_continuum": {
        "brief": "Fit continuum to mask (V1 logic).",
        "provides": ("continuum",),
        "details": "Fits the continuum to the pixels that a mask marks as unabsorbed: interpolates them, smooths the result and optionally centres it; writes 'cont'.",
        "params": [
            # CHANGED: type is now str, default is "auto"
            {"name": "fudge", "type": str, "default": "auto", "doc": "Factor that multiplies the continuum: 'auto' to centre it on the unabsorbed pixels, or a number (e.g. '1.0')."}, 
            {"name": "smooth_fwhm", "type": float, "default": 800.0, "gt": 0, "doc": "FWHM of the Gaussian smoothing of the continuum (km/s)."},
            {"name": "mask_col", "type": str, "default": "abs_mask", "doc": "Column of the mask (True for the absorbed pixels), e.g. 'abs_mask' of Find absorbed regions."},
            {"name": "renorm_model", "type": bool, "default": True, "doc": "Rescale the model to the new continuum (asked by the GUI).", "gui_hidden": True, "agents": False},
        ],
        "url": "continuum_cb.html#fit_continuum"
    },
    "fit_powerlaw": {
        "brief": "Fit power-law continuum.",
        "requires": ("z_em",),
        "details": "Fits a power law (a line in log flux vs log wavelength) to rest-frame regions free of emission and absorption, with clipping of outliers; writes 'cont_pl' (copy it to 'cont' with Apply expression to use it).",
        "params": [
            {"name": "regions", "type": str, "default": "128.0-129.0, 131.5-132.5, 134.5-136.0", "doc": "Comma-separated rest-frame wavelength ranges (nm), e.g. '128.0-129.0, 131.5-132.5'."},
            {"name": "kappa", "type": float, "default": 3.0, "gt": 0, "doc": "Threshold (in sigma) of the clipping of outliers in the regions."},
        ],
        "url": "continuum_cb.html#fit_powerlaw" # Placeholder URL
    },
    "update_from_knots": {
        "brief": "Update continuum from manual knots.",
        "provides": ("continuum",),
        "details": "Interpolates a spline through knots placed on the plot and writes it as the continuum.",
        "params": [
            {"name": "knots_x", "type": list, "default": [], "doc": "Wavelengths of the knots.", "gui_hidden": True},
            {"name": "knots_y", "type": list, "default": [], "doc": "Flux of the continuum at the knots.", "gui_hidden": True},
            {"name": "renorm_model", "type": bool, "default": True, "doc": "Rescale the model to the new continuum (asked by the GUI).", "gui_hidden": True, "agents": False},
        ],
        # No URL needed for internal interactive recipes, but you can add one if you have docs
    },
}


FWHM_TO_SIGMA = 1.0 / 2.354820045  # 1 / (2 sqrt(2 ln 2))


def _fudge(fudge):
    """The fudge factor: 'auto' (any case) or a number given as a string."""
    if str(fudge).strip().lower() == 'auto':
        return 'auto'
    try:
        return float(fudge)
    except ValueError:
        raise InvalidParameter(f"fudge: {fudge!r} is not 'auto' or a number") from None


def _smoothing_sigma(smooth_fwhm, smooth_std):
    """Sigma (km/s) of the final continuum smoothing.

    The width is given as a FWHM (findings register F16, decided with
    Guido; default 800 km/s, within 0.5% of the V1 continuum on J0103-1305).
    ``smooth_std`` is the deprecated V2 parameter (a sigma, for old scripts
    and logs) and takes precedence when given.
    """
    if smooth_std not in (None, 'None'):
        return float(smooth_std)
    return float(smooth_fwhm) * FWHM_TO_SIGMA


class RecipeContinuumV2:
    """
    Recipes for estimating and fitting the quasar continuum.

    These methods are accessed via the ``session.continuum`` attribute and delegate
    heavy logic to :class:`~astrocook.core.spectrum.SpectrumV2`.
    """
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cont' # V1 tag for this menu

    @typed_recipe
    def estimate_auto(self, smooth_len_lya: float = 5000.0, smooth_len_out: float = 400.0,
                      kappa: float = 2.0, fudge: str = 'auto',
                      smooth_fwhm: float = 800.0, template: bool = False, renorm_model: bool = True,
                      smooth_std: str = None) -> 'SessionV2':
        """
        Auto-estimate continuum (single-click).

        A single-click recipe that runs ``find_absorbed`` and then ``fit_continuum``
        in sequence, based on the V1 'clip_flux' algorithm.

        Parameters
        ----------
        smooth_len_lya : float, optional
            Smoothing length (km/s) for initial guess in the Ly-alpha forest.
            Defaults to ``5000.0``.
        smooth_len_out : float, optional
            Smoothing length (km/s) for initial guess outside the forest.
            Defaults to ``400.0``.
        kappa : float, optional
            Sigma threshold for finding absorption. Defaults to ``2.0``.
        fudge : str, optional
            A multiplicative factor. Can be a float (e.g., ``'1.0'``) or ``'auto'``.
            If ``'auto'``, the factor is calculated to center the continuum residuals.
            Defaults to ``'auto'``, as the GUI (the same result as ``'1.0'``,
            see :meth:`fit_continuum`).
        smooth_fwhm : float, optional
            FWHM (km/s) of the final Gaussian smoothing. Defaults to ``800.0``.
        smooth_std : str, optional
            Deprecated: the standard deviation (km/s) used before; overrides
            ``smooth_fwhm`` when given. Defaults to ``None``.
        template : bool, optional
            If ``True``, use a QSO template. Defaults to ``False``.
        renorm_model : bool, optional
            If ``True``, re-normalize the 'model' column to match the new continuum.
            Defaults to ``True``.

        The parameters can also be given as strings (``'5000.0'``,
        ``'True'``), as the dialogs and the session log give them: they are
        converted to their types before the recipe runs.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with updated continuum
            and absorption mask.

        Raises
        ------
        RecipeError
            If a parameter cannot be used, or the continuum cannot be
            estimated (:mod:`astrocook.recipes.errors`).
        """
        smooth_len_lya_q = smooth_len_lya * au.km/au.s
        smooth_len_out_q = smooth_len_out * au.km/au.s
        z_em_f = self._session.spec._data.z_em # Read z_em from session
        fudge_f = _fudge(fudge)
        smooth_std_f = _smoothing_sigma(smooth_fwhm, smooth_std)  # sigma, km/s

        if z_em_f == 0.0:
            logger.warning(msg_z_em_unset + " The continuum of the Lyman-alpha forest "
                           "is estimated as outside it.")

        logger.info("Finding the absorbed regions.")
        spec_with_mask = self._session.spec.find_absorbed(
            smooth_len_lya=smooth_len_lya_q,
            smooth_len_out=smooth_len_out_q,
            kappa=kappa,
            template=template
        )

        logger.info("Fitting the continuum to the unabsorbed pixels.")
        spec_with_cont = spec_with_mask.fit_continuum(
            fudge=fudge_f,
            smooth_std_kms=smooth_std_f,
            mask_col='abs_mask',
            renorm_model=renorm_model,
            template=template
        )
        return self._session.with_new_spectrum(spec_with_cont)

    @typed_recipe
    def find_absorbed(self, smooth_len_lya: float = 5000.0, smooth_len_out: float = 400.0,
                      kappa: float = 2.0, template: bool = False) -> 'SessionV2':
        """
        Find absorbed regions (V1 logic).

        Runs the V1 'clip_flux' kappa-sigma algorithm (with Ly-alpha forest split)
        to create a boolean mask (``abs_mask``) where True indicates absorption.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.find_absorbed`.

        Parameters
        ----------
        smooth_len_lya : float, optional
            Smoothing length (km/s) for initial guess in the Ly-alpha forest.
            Defaults to ``5000.0``.
        smooth_len_out : float, optional
            Smoothing length (km/s) for initial guess outside the forest.
            Defaults to ``400.0``.
        kappa : float, optional
            Sigma threshold for clipping. Flux below ``mean - kappa * sigma`` is
            considered absorbed. Defaults to ``2.0``.
        template : bool, optional
            If ``True``, use a QSO template (Not Implemented). Defaults to ``False``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with the ``abs_mask``
            column added/updated.

        Raises
        ------
        RecipeError
            If the recipe cannot run (:mod:`astrocook.recipes.errors`).
        """
        smooth_len_lya_q = smooth_len_lya * au.km/au.s
        smooth_len_out_q = smooth_len_out * au.km/au.s
        z_em_f = self._session.spec._data.z_em # Read z_em from session

        if z_em_f == 0.0:
            logger.warning(msg_z_em_unset + " The absorbed regions of the Lyman-alpha "
                           "forest are found as outside it.")
            # We don't stop, as the operation function handles z_em=0.0
            
        new_spec_v2 = self._session.spec.find_absorbed(
            smooth_len_lya=smooth_len_lya_q,
            smooth_len_out=smooth_len_out_q,
            kappa=kappa,
            template=template
        )
        return self._session.with_new_spectrum(new_spec_v2)

    @typed_recipe
    def fit_continuum(self, fudge: str = 'auto', smooth_fwhm: float = 800.0,
                      mask_col: str = 'abs_mask', renorm_model: bool = True, template: str = 'False',
                      smooth_std: str = None) -> 'SessionV2':
        """
        Fit continuum to mask (V1 logic).

        Interpolates unabsorbed regions from a mask, optionally applies an auto-fudge
        factor to center residuals, and smooths the result.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.fit_continuum`.

        Parameters
        ----------
        fudge : str, optional
            Fudge factor. Can be a float (e.g., ``'1.0'``) or ``'auto'``.
            Defaults to ``'auto'``.
        smooth_fwhm : float, optional
            FWHM (km/s) of the final Gaussian smoothing. Defaults to ``800.0``.
        smooth_std : str, optional
            Deprecated: the standard deviation (km/s) used before; overrides
            ``smooth_fwhm`` when given. Defaults to ``None``.
        mask_col : str, optional
            Name of the mask column to use (True=Absorbed). Defaults to ``'abs_mask'``.
        renorm_model : bool, optional
            If ``True``, re-normalize the 'model' column. Defaults to ``True``.
        template : str, optional
            If ``'True'``, use a QSO template. Defaults to ``'False'``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with the ``cont``
            column added/updated.

        Raises
        ------
        RecipeError
            If the recipe cannot run (:mod:`astrocook.recipes.errors`).
        """
        smooth_std_f = _smoothing_sigma(smooth_fwhm, smooth_std)
        template_b = str(template) == 'True'
        fudge_arg = _fudge(fudge)
        new_spec_v2 = self._session.spec.fit_continuum(
            fudge=fudge_arg,
            smooth_std_kms=smooth_std_f,
            mask_col=mask_col,
            renorm_model=renorm_model,
            template=template_b
        )
        return self._session.with_new_spectrum(new_spec_v2)
        
    @typed_recipe
    def fit_powerlaw(self, regions: str = '128.0-129.0, 131.5-132.5, 134.5-136.0',
                     kappa: float = 3.0) -> 'SessionV2':
        """
        Fit power-law continuum.

        Fits a power-law (linear in log-log space) to user-defined, unabsorbed
        rest-frame regions. Creates the ``cont_pl`` column.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.fit_powerlaw`.

        Parameters
        ----------
        regions : str, optional
            Comma-separated rest-frame wavelength regions (nm).
            Example: ``'128.0-129.0, 131.5-132.5'``.
            Defaults to ``'128.0-129.0, 131.5-132.5, 134.5-136.0'``.
        kappa : float, optional
            Sigma threshold for clipping outliers within regions. Defaults to ``3.0``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with the ``cont_pl``
            column added/updated.

        Raises
        ------
        RecipeError
            If the recipe cannot run (:mod:`astrocook.recipes.errors`).
        """
        if self._session.spec._data.z_em == 0.0:
            raise MissingPrerequisite(msg_z_em_unset + " The power law is fitted in "
                                      "rest-frame regions.")
        new_spec_v2 = self._session.spec.fit_powerlaw(
            regions_str=regions,
            kappa=kappa
        )
        return self._session.with_new_spectrum(new_spec_v2)
    
    @typed_recipe
    def update_from_knots(self, knots_x: list, knots_y: list, renorm_model: bool = True) -> 'SessionV2':
        """
        Update continuum from manual knots.

        Interpolates a Cubic Spline through manually placed knots and updates the
        continuum.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.update_continuum_from_knots`.

        Parameters
        ----------
        knots_x : list
            List of x coordinates (wavelength/velocity) for the knots.
        knots_y : list
            List of y coordinates (flux) for the knots.
        renorm_model : bool, optional
            If ``True``, re-normalize the 'model' column to match the new continuum.
            Defaults to ``True``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with the updated continuum.
        """
        # 1. Delegate Business Logic to Core
        new_spec = self._session.spec.update_continuum_from_knots(
            knots_x, knots_y, renorm_model
        )

        # 2. Logging (Specific to the Recipe layer)
        if self._session.log_manager:
            log_params = {
                'knots_x': list(knots_x), # Ensure list
                'knots_y': list(knots_y),
                'renorm_model': renorm_model
            }
            self._session.log_manager.add_entry("update_from_knots", log_params)

        # 3. Return new session
        return self._session.with_new_spectrum(new_spec)