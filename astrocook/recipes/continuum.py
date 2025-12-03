import logging
from typing import TYPE_CHECKING, Optional
import astropy.units as au

from astrocook.core.structures import HistoryLogV2
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
            {"name": "fudge", "type": float, "default": 1.0, "doc": "Continuum fudge factor"},
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
        "details": "Interpolate unabsorbed regions from a mask, apply fudge factor, and smooth the result to create 'cont'.",
        "params": [
            {"name": "fudge", "type": float, "default": 1.0, "doc": "Continuum fudge factor"},
            {"name": "smooth_std", "type": float, "default": 500.0, "doc": "Final Gaussian smoothing std (km/s)"},
            {"name": "mask_col", "type": str, "default": "mask_unabs", "doc": "Mask column to use"},
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
    }
}


class RecipeContinuumV2:
    """
    Recipes for estimating and fitting the quasar continuum.
    Accessed via ``session.continuum``.
    """
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cont' # V1 tag for this menu

    def find_absorbed(self, smooth_len_lya: str = '5000.0', smooth_len_out: str = '400.0', 
                        kappa: str = '2.0', template: str = 'False') -> 'SessionV2':
        """
        Identify absorbed regions using an iterative sigma-clipping algorithm.

        This recipe creates a boolean mask (``abs_mask``) where True indicates
        pixels that are statistically likely to be absorption features rather
        than continuum.

        Parameters
        ----------
        smooth_len_lya : str (float)
            Smoothing length (km/s) used for the initial continuum guess in the 
            Lyman-alpha forest region (blue side of Ly-alpha emission).
        smooth_len_out : str (float)
            Smoothing length (km/s) used redward of the Ly-alpha emission.
        kappa : str (float)
            Sigma threshold. Flux below ``mean - kappa * sigma`` is considered absorbed.
        template : str (bool)
            If 'True', use a PCA template for the initial guess (Not Implemented).

        Returns
        -------
        SessionV2
            A new session with the ``abs_mask`` column added/updated.
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

    def fit_continuum(self, fudge: str = '1.0', smooth_std: str = '500.0', 
                      mask_col: str = 'abs_mask', renorm_model: str = 'True') -> 'SessionV2':
        """
        Fit a continuum to the unabsorbed regions of the spectrum.

        This recipe interpolates across the regions identified by ``mask_col``
        and applies a Gaussian smoothing to produce the final continuum.

        Parameters
        ----------
        fudge : str (float)
            A multiplicative factor to adjust the continuum level manually.
        smooth_std : str (float)
            Standard deviation (km/s) of the Gaussian kernel used to smooth 
            the interpolated spline.
        mask_col : str
            Name of the mask column to use (default: ``'abs_mask'``).
        renorm_model : str (bool)
            If 'True', any existing absorption model (``model`` column) is 
            re-normalized to match the new continuum level.

        Returns
        -------
        SessionV2
            A new session with the ``cont`` column added/updated.
        """
        try:
            fudge_f = float(fudge)
            smooth_std_f = float(smooth_std) # This is in km/s
            # --- ADD THIS LINE ---
            renorm_model_b = str(renorm_model) == 'True'
        except ValueError:
            logging.error(msg_param_fail)
            return 0
        
        try:
            new_spec_v2 = self._session.spec.fit_continuum(
                fudge=fudge_f,
                smooth_std_kms=smooth_std_f,
                mask_col=mask_col,
                # --- ADD THIS ARGUMENT (must be added to spectrum.py) ---
                renorm_model=renorm_model_b 
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during fit_continuum: {e}", exc_info=True)
            return 0
            
    def estimate_auto(self, smooth_len_lya: str = '5000.0', smooth_len_out: str = '400.0', 
                      kappa: str = '2.0', fudge: str = '1.0', 
                      smooth_std: str = '500.0', template: str = 'False', renorm_model: str = 'True') -> 'SessionV2':
        """
        Automatically estimate the continuum in a single step.

        This recipe chains ``find_absorbed`` and ``fit_continuum`` together.

        Parameters
        ----------
        smooth_len_lya : str (float)
            Smoothing length (km/s) for initial guess (Ly-alpha forest).
        smooth_len_out : str (float)
            Smoothing length (km/s) for initial guess (redward).
        kappa : str (float)
            Sigma threshold for finding absorption.
        fudge : str (float)
            Manual adjustment factor for the final continuum.
        smooth_std : str (float)
            Final smoothing width (km/s).
        template : str (bool)
            Use PCA template (Not Implemented).
        renorm_model : str (bool)
            Re-normalize existing model?

        Returns
        -------
        SessionV2
            A new session with both ``abs_mask`` and ``cont`` columns updated.
        """
        try:
            smooth_len_lya_q = float(smooth_len_lya) * au.km/au.s
            smooth_len_out_q = float(smooth_len_out) * au.km/au.s
            kappa_f = float(kappa)
            template_b = str(template) == 'True'
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
        
    def fit_powerlaw(self, regions: str = '128.0-129.0, 131.5-132.5, 134.5-136.0',
                     kappa: str = '3.0') -> 'SessionV2':
        """
        Fit a power-law continuum to specified rest-frame spectral windows.

        This method fits a linear model in log-log space (flux vs wavelength)
        to the data within the given regions, excluding outliers via sigma-clipping.

        Parameters
        ----------
        regions : str
            Comma-separated list of wavelength ranges in nm (rest-frame).
            Example: ``'128.0-129.0, 131.5-132.5'``.
        kappa : str (float)
            Sigma threshold for clipping outliers within the fitting regions.

        Returns
        -------
        SessionV2
            A new session with the ``cont_pl`` column added/updated.
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