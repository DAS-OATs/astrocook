import logging
import numpy as np
import astropy.units as au
from typing import TYPE_CHECKING, Dict, Optional, Union

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2
from astrocook.core.spectrum import SpectrumV2


FEATURES_RECIPES_SCHEMAS = {
    "compute_ew": {
        "brief": "Compute Equivalent Width.",
        "details": "Compute the Equivalent Width (EW) of a feature within a given wavelength range.",
        "params": [
            {"name": "xmin", "type": float, "default": 0.0, "doc": "Start wavelength (nm)"},
            {"name": "xmax", "type": float, "default": 0.0, "doc": "End wavelength (nm)"},
            {"name": "z_start", "type": float, "default": 0.0, "doc": "Redshift for rest-frame conversion"},
            {"name": "use_continuum_col", "type": str, "default": "False", "doc": "Use 'cont' column if available (True/False)"}
        ],
        "url": "features.html#compute_ew" # Placeholder
    }
}

class RecipeFeaturesV2:
    """
    Recipes for extracting spectral features (Equivalent Widths, etc.).
    """
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'feat' 

    def compute_ew(self, xmin: float, xmax: float, z_start: float = 0.0, use_continuum_col: bool = False) -> Dict[str, Union[float, str]]:
        """
        Compute the Equivalent Width (EW) of a feature within a given wavelength range.

        Parameters
        ----------
        xmin : float
            Start wavelength (in nm).
        xmax : float
            End wavelength (in nm).
        z_start : float, optional
            Redshift for converting observed W to rest-frame W. Defaults to 0.0.
        use_continuum_col : bool, optional
            If True, uses the 'cont' column (if available) as the continuum level.
            If False (default), estimates a linear continuum between the edges of the selection.

        Returns
        -------
        dict
            A dictionary containing:
            - 'ew': Equivalent Width (in the same unit as x, usually nm).
            - 'ew_err': Error on EW.
            - 'centroid': Centroid wavelength.
            - 'flux': Total integrated flux (sum(1 - F/Fc)).
            - 'continuum': The continuum value used (mean or interpolated).
        """
        spec = self._session.spec
        if not spec:
            logging.error("No spectrum loaded.")
            return {}

        # 1. Extract Data in Range
        # Assuming x is in nm. 
        # TODO: Handle units robustly if x is not nm.
        
        # Get data arrays
        x = spec.x.value
        y = spec.y.value
        dy = spec.dy.value
        
        mask = (x >= xmin) & (x <= xmax)
        if np.sum(mask) < 2:
            logging.warning("Not enough points in selection for EW calculation.")
            return {}

        x_sel = x[mask]
        y_sel = y[mask]
        dy_sel = dy[mask]
        
        # 2. Determine Continuum
        if use_continuum_col and spec.has_aux_column('cont'):
            cont_sel = spec.get_column('cont').value[mask]
        else:
            # Linear Interpolation between edges
            # We take a small window around the edges to estimate the continuum level
            # Or just use the first and last point of the selection
            y_start = y_sel[0]
            y_end = y_sel[-1]
            x_start = x_sel[0]
            x_end = x_sel[-1]
            
            # Linear model: y = mx + q
            slope = (y_end - y_start) / (x_end - x_start)
            intercept = y_start - slope * x_start
            
            cont_sel = slope * x_sel + intercept

        # 3. Compute EW
        # Formula: EW = Integral(1 - F_lambda / F_cont) d_lambda
        # Discrete: Sum[(1 - y_i / c_i) * dx_i]
        
        # Determine dx (pixel widths)
        # For the integration, we can use trapezoidal rule or simple rectangle
        # Using simple rectangle with gradient-based dx
        dx_sel = np.gradient(x_sel)
        
        # Avoid division by zero
        valid_cont = (cont_sel != 0)
        
        integrand = 1.0 - (y_sel / cont_sel)
        # Apply mask for valid continuum
        integrand[~valid_cont] = 0.0
        
        ew = np.nansum(integrand * dx_sel)
        
        # 4. Compute Error on EW
        # Propagate error from dy
        # Variance(EW) = Sum [ (dy_i * dx_i / c_i)^2 ] + Continuum Error (neglected for manual linear cont)
        # Ignoring continuum error for the linear interpolation case for simplicity in V1
        
        ew_var = np.nansum( ( (dy_sel * dx_sel) / cont_sel )**2 )
        ew_err = np.sqrt(ew_var)

        # 5. Compute Centroid
        # Centroid = Sum(lambda * (1 - F/Fc)) / Sum(1 - F/Fc) = Sum(lambda * integrand) / Sum(integrand)
        # Use abs(integrand) for centroid weighting? Standard is usually just depths.
        # But if it's an emission line, 1-F/Fc is negative.
        # Let's use the flux deficit/excess directly as weights.
        
        flux_diff = (cont_sel - y_sel) # Positive for absorption
        centroid = np.nansum(x_sel * flux_diff) / np.nansum(flux_diff)
        
        return {
            'ew': ew,
            'ew_err': ew_err,
            'centroid': centroid,
            'xmin': xmin,
            'xmax': xmax,
            'cont_source': 'column' if use_continuum_col else 'linear_interp'
        }
