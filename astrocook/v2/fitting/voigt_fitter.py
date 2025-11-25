import astropy.constants as const
import astropy.units as au
from scipy.ndimage import gaussian_filter1d # <<< FASTER CONVOLUTION
import dataclasses
import logging
import numpy as np
from scipy.optimize import least_squares
from scipy.special import wofz
from scipy.linalg import inv
from typing import List, Tuple, Dict, Any, Optional

from ..atomic_data import ATOM_DATA, STANDARD_MULTIPLETS
from ..spectrum import SpectrumV2
from ..structures import ComponentDataV2, SystemListDataV2
from ..system_list import SystemListV2
from .voigt_model import VoigtModelConstraintV2

class VoigtFitterV2:
    def __init__(self, spectrum: SpectrumV2, system_list: SystemListV2):
        self._spectrum = spectrum
        self._system_list = system_list
        self._constraints = system_list.constraint_model
        
        self._x_ang = self._spectrum.x.to(au.Angstrom).value
        
        if self._spectrum.cont is not None:
            self._cont = self._spectrum.cont.value
            self._cont[self._cont == 0] = np.nan 
        else:
            logging.warning("No continuum found. Assuming flux is already normalized.")
            self._cont = np.ones_like(self._spectrum.y.value)

        self._y_norm = self._spectrum.y.value / self._cont
        
        # Clip dy to a minimum floor to prevent Inf Chi2
        raw_dy = self._spectrum.dy.value
        min_err = 1e-9 * np.nanmedian(self._cont) 
        safe_dy = np.maximum(raw_dy, min_err)
        
        self._dy_norm = safe_dy / self._cont
        
        self._mask = np.isfinite(self._y_norm) & np.isfinite(self._dy_norm) & (self._dy_norm > 0)
        
        # Prepare Resolution (Speed Optimization)
        # We now calculate sigma_pix once during init
        self._resol_sigma_pix = None
        
        if self._spectrum.has_aux_column('resol'):
            resol_col = self._spectrum.get_column('resol')
            if resol_col is not None:
                vals = resol_col.value
                median_resol = np.nanmedian(vals)
                if median_resol > 0:
                    # FWHM_ang = lambda / R
                    # sigma_ang = FWHM / 2.355
                    # sigma_pix = sigma_ang / dx
                    
                    mid_wave = np.median(self._x_ang)
                    fwhm_ang = mid_wave / median_resol
                    sigma_ang = fwhm_ang / 2.35482
                    
                    dx_avg = np.mean(np.diff(self._x_ang))
                    self._resol_sigma_pix = sigma_ang / dx_avg
                    
                    logging.info(f"VoigtFitter: R={median_resol:.0f} -> sigma_pix={self._resol_sigma_pix:.3f}")

    def _voigt_optical_depth(self, wave_grid_ang: np.ndarray, 
                             lambda_0: float, f_val: float, gamma: float,
                             z: float, N_col: float, b_kms: float) -> np.ndarray:
        lambda_c = lambda_0 * (1.0 + z)
        c_kms = 2.99792458e5
        b_safe = np.maximum(b_kms, 0.1) 
        dop_width_ang = (b_safe / c_kms) * lambda_c
        dop_width_ang = np.maximum(dop_width_ang, 1e-10)

        x = (wave_grid_ang - lambda_c) / dop_width_ang
        c_ang_s = 2.99792458e18
        a = (gamma * lambda_c**2) / (4.0 * np.pi * c_ang_s * dop_width_ang)
        
        z_complex = x + 1j * a
        H_ax = wofz(z_complex).real
        
        prefactor = 1.4974e-15
        tau = prefactor * N_col * f_val * lambda_0 / b_safe * H_ax
        return tau

    def _compute_model(self, p_free: np.ndarray) -> np.ndarray:
        p_full = self._constraints.map_p_free_to_full(p_free)
        tau_total = np.zeros_like(self._x_ang)
        
        components = self._system_list.components
        num_params_per_comp = 4
        
        for i, comp in enumerate(components):
            idx = i * num_params_per_comp
            z = p_full[idx]
            logN = p_full[idx + 1]
            b_val = p_full[idx + 2]
            btur = p_full[idx + 3]
            
            b_eff = np.sqrt(b_val**2 + btur**2)
            N_col = 10**logN
            
            transitions_to_compute = []
            if comp.series in STANDARD_MULTIPLETS:
                transitions_to_compute = STANDARD_MULTIPLETS[comp.series]
            elif comp.series in ATOM_DATA:
                transitions_to_compute = [comp.series]
            elif ',' in comp.series:
                transitions_to_compute = comp.series.split(',')
            else:
                continue

            for trans_name in transitions_to_compute:
                trans_name = trans_name.strip()
                atom = ATOM_DATA.get(trans_name)
                if not atom: continue
                
                obs_cent = atom['wave'] * (1 + z)
                if obs_cent < self._x_ang[0] - 100 or obs_cent > self._x_ang[-1] + 100:
                    continue

                tau_comp = self._voigt_optical_depth(
                    self._x_ang, atom['wave'], atom['f'], atom['gamma'],
                    z, N_col, b_eff
                )
                tau_total += tau_comp

        flux_model = np.exp(-tau_total)
        
        # Optimization: Use fast Gaussian filter instead of Astropy convolve
        if self._resol_sigma_pix:
            flux_model = gaussian_filter1d(flux_model, self._resol_sigma_pix)

        return flux_model

    def _residual_function(self, p_free: np.ndarray) -> np.ndarray:
        model = self._compute_model(p_free)
        resid = (self._y_norm - model) / self._dy_norm
        resid[~np.isfinite(resid)] = 0.0
        return resid[self._mask]

    def fit(self, max_nfev: int = 2000, method: str = 'trf', verbose: int = 1) -> Tuple[SystemListV2, np.ndarray, Any]:
        logging.info("Starting Voigt Fit...")
        
        p0 = self._constraints.p_free_vector
        bounds = self._constraints.get_bounds()
        
        if len(p0) == 0:
            logging.warning("No free parameters to fit.")
            model = self._compute_model(p0) * self._cont
            return self._system_list, model, None

        res = least_squares(
            self._residual_function, p0, bounds=bounds, 
            method=method, loss='linear', max_nfev=max_nfev, 
            x_scale='jac', verbose=verbose
        )
        
        p_free_errors = np.zeros_like(res.x)
        try:
            n_data = np.sum(self._mask)
            n_param = len(p0)
            dof = max(1, n_data - n_param)
            chi2 = np.sum(res.fun**2)
            red_chi2 = chi2 / dof
            
            J = res.jac
            # Unscaled errors first
            cov = inv(J.T @ J) 
            p_free_errors = np.sqrt(np.diag(cov))
            # Optional: Multiply by sqrt(red_chi2) if needed, 
            # but unscaled is often preferred for Voigt fitting.
        except Exception as e:
            logging.warning(f"Could not calculate parameter errors: {e}")
            p_free_errors[:] = np.nan

        p_fitted_full = self._constraints.map_p_free_to_full(res.x)
        is_free_mask = self._constraints._param_map['is_free']
        p_full_errors = np.zeros_like(p_fitted_full)
        p_full_errors[is_free_mask] = p_free_errors
        
        new_components = []
        old_components = self._system_list.components
        num_params = 4
        
        for i, old_comp in enumerate(old_components):
            idx = i * num_params
            
            dz_new = p_full_errors[idx] if is_free_mask[idx] else old_comp.dz
            dlogN_new = p_full_errors[idx+1] if is_free_mask[idx+1] else old_comp.dlogN
            db_new = p_full_errors[idx+2] if is_free_mask[idx+2] else old_comp.db
            
            new_comp = dataclasses.replace(
                old_comp,
                z=p_fitted_full[idx],
                logN=p_fitted_full[idx + 1],
                b=p_fitted_full[idx + 2],
                btur=p_fitted_full[idx + 3],
                dz=dz_new,
                dlogN=dlogN_new,
                db=db_new,
            )
            new_components.append(new_comp)
            
        new_data_core = dataclasses.replace(self._system_list._data, components=new_components)
        new_system_list = SystemListV2(new_data_core)
        
        final_model_norm = self._compute_model(res.x)
        final_model_flux = final_model_norm * self._cont

        logging.info(f"Fit complete. Chi2: {chi2:.2f} (Red: {red_chi2:.2f})")
        return new_system_list, final_model_flux, res

    def compute_model_flux(self) -> Tuple[np.ndarray, np.ndarray]:
        p_current_free = self._constraints.p_free_vector
        y_model_norm = self._compute_model(p_current_free)
        y_model_phys = y_model_norm * self._cont
        return self._x_ang, y_model_phys