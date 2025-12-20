import astropy.units as au
import logging
import numpy as np
from scipy.interpolate import CubicSpline
from typing import TYPE_CHECKING, Optional

from astrocook.core.structures import HistoryLogV2, LogEntryV2
from astrocook.legacy.message import msg_param_fail

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2

# --- Schemas for the Continuum Menu ---
CONTINUUM_RECIPES_SCHEMAS = {
    "estimate_auto": {
        "brief": "Auto-estimate continuum (single-click).",
        "details": "A single-click recipe that runs 'find_absorbed' and then 'fit_continuum' in sequence, based on the V1 'clip_flux' algorithm.",
        "params": [
            {"name": "smooth_len_lya", "type": float, "default": 5000.0, "doc": "Smoothing length in Ly-a forest (km/s)"},
            {"name": "smooth_len_out", "type": float, "default": 400.0, "doc": "Smoothing length outside Ly-a forest (km/s)"},
            {"name": "kappa", "type": float, "default": 2.0, "doc": "Sigma threshold for clipping"},
            {"name": "fudge", "type": str, "default": "auto", "doc": "Fudge factor. Use 'auto' or a float (e.g., '1.0')."}, 
            {"name": "smooth_std", "type": float, "default": 500.0, "doc": "Final Gaussian smoothing std (km/s)"},
            {"name": "template", "type": bool, "default": False, "doc": "Use QSO template (NOT IMPLEMENTED)"},
            {"name": "renorm_model", "type": bool, "default": True, "doc": "Also re-normalize 'model'?", "gui_hidden": True},
        ],
        "url": "continuum_cb.html#estimate_auto"
    },
    "find_absorbed": {
        "brief": "Find absorbed regions (V1 logic).",
        "details": "Run the V1 'clip_flux' kappa-sigma algorithm (with Ly-a split) to create a 'mask_unabs' column.",
        "params": [
            {"name": "smooth_len_lya", "type": float, "default": 5000.0, "doc": "Smoothing length in Ly-a forest (km/s)"},
            {"name": "smooth_len_out", "type": float, "default": 400.0, "doc": "Smoothing length outside Ly-a forest (km/s)"},
            {"name": "kappa", "type": float, "default": 2.0, "doc": "Sigma threshold for clipping"},
            {"name": "template", "type": bool, "default": False, "doc": "Use QSO template (NOT IMPLEMENTED)"},
        ],
        "url": "continuum_cb.html#find_absorbed"
    },
    "fit_continuum": {
        "brief": "Fit continuum to mask (V1 logic).",
        "details": "Interpolate unabsorbed regions from a mask, optionally auto-center (fudge='auto'), and smooth result.",
        "params": [
            # CHANGED: type is now str, default is "auto"
            {"name": "fudge", "type": str, "default": "auto", "doc": "Fudge factor. Use 'auto' or a float (e.g., '1.0')."}, 
            {"name": "smooth_std", "type": float, "default": 500.0, "doc": "Final Gaussian smoothing std (km/s)"},
            {"name": "mask_col", "type": str, "default": "abs_mask", "doc": "Mask column to use (True=Absorbed)"},
            {"name": "renorm_model", "type": bool, "default": True, "doc": "Also re-normalize 'model'?", "gui_hidden": True},
        ],
        "url": "continuum_cb.html#fit_continuum"
    },
    "fit_powerlaw": {
        "brief": "Fit power-law continuum.",
        "details": "Fit a power-law (linear in log-log space) to user-defined, unabsorbed rest-frame regions. Creates 'cont_pl' column.",
        "params": [
            {"name": "regions", "type": str, "default": "128.0-129.0, 131.5-132.5, 134.5-136.0", "doc": "Comma-separated rest-frame regions (nm), e.g., '128-129,131-132'"},
            {"name": "kappa", "type": float, "default": 3.0, "doc": "Sigma threshold for clipping outliers"},
        ],
        "url": "continuum_cb.html#fit_powerlaw" # Placeholder URL
    },
    "update_from_knots": {
        "brief": "Update continuum from manual knots.",
        "details": "Interpolate a spline through manually placed knots and update the continuum.",
        "params": [
            {"name": "knots_x", "type": list, "default": [], "doc": "List of x coordinates", "gui_hidden": True},
            {"name": "knots_y", "type": list, "default": [], "doc": "List of y coordinates", "gui_hidden": True},
            {"name": "renorm_model", "type": bool, "default": True, "doc": "Also re-normalize 'model'?", "gui_hidden": True},
        ],
        # No URL needed for internal interactive recipes, but you can add one if you have docs
    },
}


class RecipeContinuumV2:
    """
    Recipes for estimating and fitting the quasar continuum.

    These methods are accessed via the ``session.continuum`` attribute and delegate
    heavy logic to :class:`~astrocook.core.spectrum.SpectrumV2`.
    """
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cont' # V1 tag for this menu

    def estimate_auto(self, smooth_len_lya: str = '5000.0', smooth_len_out: str = '400.0', 
                      kappa: str = '2.0', fudge: str = '1.0', 
                      smooth_std: str = '500.0', template: str = 'False', renorm_model: str = 'True') -> 'SessionV2':
        """
        Auto-estimate continuum (single-click).

        A single-click recipe that runs ``find_absorbed`` and then ``fit_continuum``
        in sequence, based on the V1 'clip_flux' algorithm.

        Parameters
        ----------
        smooth_len_lya : str, optional
            Smoothing length (km/s) for initial guess in the Ly-alpha forest.
            Defaults to ``'5000.0'``.
        smooth_len_out : str, optional
            Smoothing length (km/s) for initial guess outside the forest.
            Defaults to ``'400.0'``.
        kappa : str, optional
            Sigma threshold for finding absorption. Defaults to ``'2.0'``.
        fudge : str, optional
            A multiplicative factor. Can be a float (e.g., ``'1.0'``) or ``'auto'``.
            If ``'auto'``, the factor is calculated to center the continuum residuals.
            Defaults to ``'1.0'``.
        smooth_std : str, optional
            Standard deviation (km/s) for the final Gaussian smoothing.
            Defaults to ``'500.0'``.
        template : str, optional
            If ``'True'``, use a QSO template (Not Implemented). Defaults to ``'False'``.
        renorm_model : str, optional
            If ``'True'``, re-normalize the 'model' column to match the new continuum.
            Defaults to ``'True'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with updated continuum
            and absorption mask, or 0 on failure.
        """
        try:
            smooth_len_lya_q = float(smooth_len_lya) * au.km/au.s
            smooth_len_out_q = float(smooth_len_out) * au.km/au.s
            kappa_f = float(kappa)
            template_b = str(template) == 'True'
            # Allow 'auto' or parse float
            if fudge.lower() == 'auto':
                fudge_f = 'auto'
            else:
                fudge_f = float(fudge)
            smooth_std_f = float(smooth_std) # This is in km/s
            z_em_f = self._session.spec._data.z_em # Read z_em from session
            renorm_model_b = str(renorm_model) == 'True'
        except ValueError:
            logging.error(msg_param_fail)
            return 0
            
        if z_em_f == 0.0:
            logging.warning("z_em is 0.0. Run 'Edit > Set Properties' first.")

        try:
            logging.info("Auto-continuum: Finding absorbed regions...")
            # 1. Call the first API method
            spec_with_mask = self._session.spec.find_absorbed(
                smooth_len_lya=smooth_len_lya_q,
                smooth_len_out=smooth_len_out_q,
                kappa=kappa_f,
                template=template_b
            )
            
            logging.info("Auto-continuum: Fitting continuum to mask...")
            spec_with_cont = spec_with_mask.fit_continuum(
                fudge=fudge_f,
                smooth_std_kms=smooth_std_f,
                mask_col='abs_mask',
                renorm_model=renorm_model_b
            )
            
            # 3. Return the final new state
            return self._session.with_new_spectrum(spec_with_cont)
        except Exception as e:
            logging.error(f"Failed during estimate_auto: {e}", exc_info=True)
            return 0

    def find_absorbed(self, smooth_len_lya: str = '5000.0', smooth_len_out: str = '400.0', 
                        kappa: str = '2.0', template: str = 'False') -> 'SessionV2':
        """
        Find absorbed regions (V1 logic).

        Runs the V1 'clip_flux' kappa-sigma algorithm (with Ly-alpha forest split)
        to create a boolean mask (``abs_mask``) where True indicates absorption.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.find_absorbed`.

        Parameters
        ----------
        smooth_len_lya : str, optional
            Smoothing length (km/s) for initial guess in the Ly-alpha forest.
            Defaults to ``'5000.0'``.
        smooth_len_out : str, optional
            Smoothing length (km/s) for initial guess outside the forest.
            Defaults to ``'400.0'``.
        kappa : str, optional
            Sigma threshold for clipping. Flux below ``mean - kappa * sigma`` is
            considered absorbed. Defaults to ``'2.0'``.
        template : str, optional
            If ``'True'``, use a QSO template (Not Implemented). Defaults to ``'False'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with the ``abs_mask``
            column added/updated, or 0 on failure.
        """
        try:
            smooth_len_lya_q = float(smooth_len_lya) * au.km/au.s
            smooth_len_out_q = float(smooth_len_out) * au.km/au.s
            kappa_f = float(kappa)
            template_b = str(template) == 'True'
            z_em_f = self._session.spec._data.z_em # Read z_em from session
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        if z_em_f == 0.0:
            logging.warning("z_em is 0.0. Run 'Edit > Set Properties' before continuum fitting.")
            # We don't stop, as the operation function handles z_em=0.0
            
        try:
            new_spec_v2 = self._session.spec.find_absorbed(
                smooth_len_lya=smooth_len_lya_q,
                smooth_len_out=smooth_len_out_q,
                kappa=kappa_f,
                template=template_b
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during find_absorbed: {e}", exc_info=True)
            return 0

    def fit_continuum(self, fudge: str = 'auto', smooth_std: str = '500.0', 
                      mask_col: str = 'abs_mask', renorm_model: str = 'True') -> 'SessionV2':
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
        smooth_std : str, optional
            Final Gaussian smoothing standard deviation (km/s). Defaults to ``'500.0'``.
        mask_col : str, optional
            Name of the mask column to use (True=Absorbed). Defaults to ``'abs_mask'``.
        renorm_model : str, optional
            If ``'True'``, re-normalize the 'model' column. Defaults to ``'True'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with the ``cont``
            column added/updated, or 0 on failure.
        """
        try:
            smooth_std_f = float(smooth_std)
            renorm_model_b = str(renorm_model) == 'True'
            
            # Allow 'auto' or parse float
            if fudge.lower() == 'auto':
                fudge_arg = 'auto'
            else:
                fudge_arg = float(fudge)
                
        except ValueError:
            logging.error(msg_param_fail)
            return 0
        
        try:
            new_spec_v2 = self._session.spec.fit_continuum(
                fudge=fudge_arg,
                smooth_std_kms=smooth_std_f,
                mask_col=mask_col,
                renorm_model=renorm_model_b 
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during fit_continuum: {e}", exc_info=True)
            return 0
        
    def fit_powerlaw(self, regions: str = '128.0-129.0, 131.5-132.5, 134.5-136.0',
                     kappa: str = '3.0') -> 'SessionV2':
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
        kappa : str, optional
            Sigma threshold for clipping outliers within regions. Defaults to ``'3.0'``.

        Returns
        -------
        SessionV2 or int
            A new :class:`~astrocook.core.session.SessionV2` with the ``cont_pl``
            column added/updated, or 0 on failure.
        """
        if self._session.spec._data.z_em == 0.0:
            logging.error("Cannot fit power-law: Emission Redshift (z_em) is 0.0. "
                          "Run 'Edit > Set Properties' first.")
            return 0
            
        try:
            kappa_f = float(kappa) # Parse kappa

            new_spec_v2 = self._session.spec.fit_powerlaw(
                regions_str=regions,
                kappa=kappa_f
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during fit_powerlaw: {e}", exc_info=True)
            return 0
    
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