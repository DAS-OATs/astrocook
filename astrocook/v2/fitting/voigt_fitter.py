import astropy.constants as const
import astropy.units as au
from astropy.convolution import convolve, Gaussian1DKernel
import dataclasses
import logging
import numpy as np
from scipy.optimize import least_squares
from scipy.special import wofz
from typing import List, Tuple, Dict, Any, Optional

from ..atomic_data import ATOM_DATA, STANDARD_MULTIPLETS
from ..spectrum import SpectrumV2
from ..structures import ComponentDataV2, SystemListDataV2
from ..system_list import SystemListV2
from .voigt_model import VoigtModelConstraintV2

class VoigtFitterV2:
    """
    The Physics Engine for Astrocook V2.
    
    Responsibilities:
    1. Translate SystemListV2 + Constraints -> Flux Model.
    2. Handle "Intrinsic" constraints (Multiplet expansion) natively.
    3. Convolve with Instrument Resolution (Gaussian approximation).
    4. Minimize Chi^2 using SciPy.
    """
    
    def __init__(self, spectrum: SpectrumV2, system_list: SystemListV2):
        self._spectrum = spectrum
        self._system_list = system_list
        self._constraints = system_list.constraint_model
        
        # 1. Prepare Data for Fitting (Performance Optimization)
        # We assume the user has already sliced/chunked the spectrum if they only want to fit a region.
        # We work in Angstroms for physics calculations.
        self._x_ang = self._spectrum.x.to(au.Angstrom).value
        
        # We fit Normalized Flux: Flux / Continuum
        if self._spectrum.cont is not None:
            self._cont = self._spectrum.cont.value
            # Avoid division by zero
            self._cont[self._cont == 0] = np.nan 
        else:
            logging.warning("No continuum found. Assuming flux is already normalized.")
            self._cont = np.ones_like(self._spectrum.y.value)

        # Normalize Flux and Error
        self._y_norm = self._spectrum.y.value / self._cont
        self._dy_norm = self._spectrum.dy.value / self._cont
        
        # Mask NaNs or bad pixels
        self._mask = np.isfinite(self._y_norm) & np.isfinite(self._dy_norm) & (self._dy_norm > 0)
        
        # 2. Prepare Resolution (Instrumental broadening)
        # Check for 'resol' column (Resolution R = lambda / d_lambda)
        self._resol_array = None
        if self._spectrum.has_aux_column('resol'):
            resol_col = self._spectrum.get_column('resol')
            if resol_col is not None:
                vals = resol_col.value
                # Use median resolution for convolution kernel to avoid varying kernel cost 
                # (Optimization for typical UVES/HIRES data where R is roughly constant in a chunk)
                # TODO: Implement variable kernel convolution if R varies significantly.
                self._resol_array = vals
                self._median_resol = np.nanmedian(vals)
                if np.isnan(self._median_resol) or self._median_resol <= 0:
                    self._median_resol = None
            else:
                self._median_resol = None
        else:
            self._median_resol = None
            
        logging.info(f"VoigtFitter initialized. Data points: {np.sum(self._mask)}. "
                     f"Resolution: {self._median_resol if self._median_resol else 'None'}")

    def _voigt_optical_depth(self, wave_grid_ang: np.ndarray, 
                             lambda_0: float, f_val: float, gamma: float,
                             z: float, N_col: float, b_kms: float) -> np.ndarray:
        """
        Calculates the optical depth tau(lambda) for a single transition.
        Formula: tau(x) ~ N * f * H(a, x)
        """
        # 1. Central wavelength observed
        lambda_c = lambda_0 * (1.0 + z)
        
        # 2. Doppler width in Angstroms: sigma_lambda = (b/c) * lambda_c
        c_kms = 2.99792458e5
        # Ensure b is not zero to prevent division errors (thermal floor)
        b_safe = np.maximum(b_kms, 0.1) 
        dop_width_ang = (b_safe / c_kms) * lambda_c
        
        # 3. Normalized frequency coordinate x
        # x = (lambda - lambda_c) / dop_width
        x = (wave_grid_ang - lambda_c) / dop_width_ang
        
        # 4. Damping parameter a (Lorentzian width / Doppler width)
        # gamma is damping constant (s^-1)
        # c in Ang/s = 2.998e18
        c_ang_s = 2.99792458e18
        # a = (gamma * lambda_c^2) / (4 * pi * c * dop_width_ang)
        # Note: Physics convention often defines 'a' at rest frame or observed frame.
        # Standard approx: a = (gamma * lambda_0) / (4 * pi * b_vel) ? 
        # Let's use the explicit derivation:
        a = (gamma * lambda_c**2) / (4.0 * np.pi * c_ang_s * dop_width_ang)
        
        # 5. Voigt Profile H(a, x)
        z_complex = x + 1j * a
        H_ax = wofz(z_complex).real
        
        # 6. Optical Depth Pre-factor
        # Constant K = sqrt(pi) * e^2 / (m_e * c_cgs) * (1e-13 for Angstroms/km/s units handling)
        # Using the standard numerical constant for N in cm^-2, lambda in Ang, b in km/s:
        # tau = 1.497e-2 * N * f * lambda_0(Ang) / b(km/s) * H(a,x)
        # (Note: lambda_0 here is Rest Wavelength, consistent with VPFIT/V1)
        
        prefactor = 1.4974e-15
        
        tau = prefactor * N_col * f_val * lambda_0 / b_safe * H_ax
        
        return tau

    def _compute_model(self, p_free: np.ndarray) -> np.ndarray:
        """
        The Core Model Generator.
        1. Maps Optimization Vector (p_free) -> Physical Parameters (p_full).
        2. Iterates Components ("Intrinsic" Multiplet Logic).
        3. Sums Optical Depths.
        4. Convolves with Resolution.
        5. Returns Normalized Flux Model.
        """
        
        # 1. Expand Constraints (Extrinsic Links)
        # The constraint model handles the user-defined links here.
        p_full = self._constraints.map_p_free_to_full(p_free)
        
        # Initialize Tau
        tau_total = np.zeros_like(self._x_ang)
        
        # 2. Iterate Components
        # We stride through p_full in blocks of 4 (z, logN, b, btur)
        components = self._system_list.components
        num_params_per_comp = 4
        
        for i, comp in enumerate(components):
            idx = i * num_params_per_comp
            
            # Extract parameters for this component
            z = p_full[idx]
            logN = p_full[idx + 1]
            b_val = p_full[idx + 2]
            btur = p_full[idx + 3]
            
            # Combine thermal and turbulent broadening
            # b_eff^2 = b_therm^2 + b_tur^2
            b_eff = np.sqrt(b_val**2 + btur**2)
            N_col = 10**logN
            
            # --- INTRINSIC CONSTRAINT LOGIC (Multiplets) ---
            # Determine which lines this component generates.
            # If comp.series is a key in STANDARD_MULTIPLETS (e.g. 'CIV'), 
            # we generate ALL lines in that multiplet using the SAME z, N, b.
            
            transitions_to_compute = []
            if comp.series in STANDARD_MULTIPLETS:
                transitions_to_compute = STANDARD_MULTIPLETS[comp.series]
            elif comp.series in ATOM_DATA:
                transitions_to_compute = [comp.series]
            else:
                # Fallback: maybe the series string is a comma-separated list? (Legacy V1)
                if ',' in comp.series:
                    transitions_to_compute = comp.series.split(',')
                else:
                    # Log warning only if verbose? 
                    continue

            # Sum contributions from all lines in the multiplet
            for trans_name in transitions_to_compute:
                atom = ATOM_DATA.get(trans_name)
                if not atom:
                    continue
                
                # Check if line is within spectral range (optimization)
                obs_cent = atom['wave'] * (1 + z)
                if obs_cent < self._x_ang[0] - 50 or obs_cent > self._x_ang[-1] + 50:
                    continue

                tau_comp = self._voigt_optical_depth(
                    self._x_ang,
                    atom['wave'],
                    atom['f'],
                    atom['gamma'],
                    z,
                    N_col,
                    b_eff
                )
                tau_total += tau_comp

        # 3. Compute Flux
        flux_model = np.exp(-tau_total)
        
        # 4. Instrument Convolution
        if self._median_resol:
            # FWHM_obs = lambda / R
            # sigma_obs = FWHM / 2.355
            # We approximate this with a constant kernel over the chunk if R is constant.
            
            mid_wave = np.median(self._x_ang)
            fwhm_ang = mid_wave / self._median_resol
            sigma_ang = fwhm_ang / 2.35482
            
            # Convert sigma to pixels
            dx_avg = np.mean(np.diff(self._x_ang))
            sigma_pix = sigma_ang / dx_avg
            
            if sigma_pix > 0.5:
                kernel = Gaussian1DKernel(stddev=sigma_pix)
                # handle boundary to avoid edge effects
                flux_model = convolve(flux_model, kernel, boundary='extend')

        return flux_model

    def _residual_function(self, p_free: np.ndarray) -> np.ndarray:
        """
        Cost function for Least Squares.
        Returns weighted residuals: (Data - Model) / Error
        """
        model = self._compute_model(p_free)
        
        # Calculate residual on normalized flux
        resid = (self._y_norm - model) / self._dy_norm
        
        # Apply mask (remove NaNs, bad pixels)
        return resid[self._mask]

    def fit(self, max_nfev: int = 5000, method: str = 'trf', verbose: int = 1) -> Tuple[SystemListV2, Any]:
        logging.info("Starting Voigt Fit...")
        
        p0 = self._constraints.p_free_vector
        
        if len(p0) == 0:
            logging.warning("No free parameters to fit.")
            return self._system_list, None

        # Get Bounds from Constraints
        bounds = self._constraints.get_bounds()

        # Run Optimization with Scaling
        # x_scale='jac' automatically rescales parameters so they have comparable influence
        res = least_squares(
            self._residual_function, 
            p0, 
            bounds=bounds, 
            method=method,
            loss='linear',
            max_nfev=max_nfev,
            x_scale='jac',  # <--- CRITICAL FIX FOR CONVERGENCE
            verbose=verbose
        )
        
        if not res.success:
            logging.warning(f"Fit warning: {res.message}")
            
        # Map results back
        p_fitted_full = self._constraints.map_p_free_to_full(res.x)
        
        new_components = []
        old_components = self._system_list.components
        num_params = 4
        
        for i, old_comp in enumerate(old_components):
            idx = i * num_params
            
            new_comp = dataclasses.replace(
                old_comp,
                z=p_fitted_full[idx],
                logN=p_fitted_full[idx + 1],
                b=p_fitted_full[idx + 2],
                btur=p_fitted_full[idx + 3]
            )
            new_components.append(new_comp)
            
        new_data_core = dataclasses.replace(
            self._system_list._data,
            components=new_components
        )
        
        new_system_list = SystemListV2(new_data_core)
        
        logging.info(f"Fit complete. Final Chi2: {np.sum(res.fun**2):.2f}")
        
        return new_system_list, res    
    
    def compute_model_flux(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Public API to get the current model flux (scaled to continuum) for plotting.
        Does NOT fit, just evaluates.
        """
        p_current_free = self._constraints.p_free_vector
        y_model_norm = self._compute_model(p_current_free)
        
        # Scale back to physical flux
        y_model_phys = y_model_norm * self._cont
        
        return self._x_ang, y_model_phys