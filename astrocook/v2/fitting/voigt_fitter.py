import astropy.constants as const
import astropy.units as au
from scipy.ndimage import gaussian_filter1d
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
        
        raw_dy = self._spectrum.dy.value
        min_err = 1e-9 * np.nanmedian(self._cont) 
        safe_dy = np.maximum(raw_dy, min_err)
        self._dy_norm = safe_dy / self._cont
        
        self._valid_mask = np.isfinite(self._y_norm) & np.isfinite(self._dy_norm) & (self._dy_norm > 0)
        
        # --- 1. Dynamic Fit Mask Construction ---
        self._fit_mask = np.zeros_like(self._valid_mask, dtype=bool)
        
        c_kms = 299792.458
        threshold = 0.9
        
        active_uuids = self._constraints._active_uuids
        has_active = False
        
        for comp in self._system_list.components:
            if active_uuids is not None and comp.uuid not in active_uuids:
                continue
            has_active = True
            
            # Identify transitions
            if comp.series in STANDARD_MULTIPLETS:
                trans_list = STANDARD_MULTIPLETS[comp.series]
            elif comp.series in ATOM_DATA:
                trans_list = [comp.series]
            else: continue
            
            for trans in trans_list:
                if trans not in ATOM_DATA: continue
                lam_0 = ATOM_DATA[trans]['wave']
                lam_obs = lam_0 * (1.0 + comp.z)
                
                # Find nearest pixel
                idx_center = np.searchsorted(self._x_ang, lam_obs)
                if idx_center <= 0 or idx_center >= len(self._x_ang) - 1: continue
                
                start_idx = idx_center
                end_idx = idx_center
                
                while start_idx > 0 and self._y_norm[start_idx] < threshold:
                    start_idx -= 1
                
                while end_idx < len(self._y_norm) - 1 and self._y_norm[end_idx] < threshold:
                    end_idx += 1
                
                width_pix = end_idx - start_idx
                if width_pix < 3:
                    dw_fallback = (50.0 / c_kms) * lam_obs
                    mask_fallback = (self._x_ang > lam_obs - dw_fallback) & (self._x_ang < lam_obs + dw_fallback)
                    self._fit_mask |= mask_fallback
                else:
                    start_idx = max(0, start_idx - 1)
                    end_idx = min(len(self._x_ang), end_idx + 1)
                    self._fit_mask[start_idx:end_idx] = True

        if not has_active or np.sum(self._fit_mask) == 0:
            self._fit_mask = self._valid_mask
        else:
            self._fit_mask &= self._valid_mask
            
        logging.info(f"VoigtFitter: Fit mask defined on {np.sum(self._fit_mask)} pixels.")

        self._resol_sigma_pix = None
        if self._spectrum.has_aux_column('resol'):
            resol_col = self._spectrum.get_column('resol')
            if resol_col is not None:
                vals = resol_col.value
                median_resol = np.nanmedian(vals)
                if median_resol > 0:
                    mid_wave = np.median(self._x_ang)
                    fwhm_ang = mid_wave / median_resol
                    sigma_ang = fwhm_ang / 2.35482
                    dx_avg = np.mean(np.diff(self._x_ang))
                    self._resol_sigma_pix = sigma_ang / dx_avg
                    logging.info(f"VoigtFitter: R={median_resol:.0f} -> sigma_pix={self._resol_sigma_pix:.3f}")
        
        elif self._spectrum.meta.get('resol_fwhm', 0.0) > 0:
            resol_fwhm_kms = self._spectrum.meta['resol_fwhm']
            mid_wave = np.median(self._x_ang)
            fwhm_ang = (resol_fwhm_kms / c_kms) * mid_wave
            sigma_ang = fwhm_ang / 2.35482
            dx_avg = np.mean(np.diff(self._x_ang))
            self._resol_sigma_pix = sigma_ang / dx_avg
            logging.info(f"VoigtFitter: Meta FWHM={resol_fwhm_kms} km/s -> sigma_pix={self._resol_sigma_pix:.3f}")

        # --- CACHING OPTIMIZATION ---
        # We will pre-calculate tau for frozen components in fit()
        self._cached_tau_static = None

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

    def _compute_tau_component(self, comp: ComponentDataV2, z: float, logN: float, b: float, btur: float) -> np.ndarray:
        """Helper to compute tau for a single component."""
        tau_comp = np.zeros_like(self._x_ang)
        b_eff = np.sqrt(b**2 + btur**2)
        N_col = 10**logN
        
        transitions_to_compute = []
        if comp.series in STANDARD_MULTIPLETS:
            transitions_to_compute = STANDARD_MULTIPLETS[comp.series]
        elif comp.series in ATOM_DATA:
            transitions_to_compute = [comp.series]
        elif ',' in comp.series:
            transitions_to_compute = comp.series.split(',')
        else: return tau_comp

        for trans_name in transitions_to_compute:
            trans_name = trans_name.strip()
            atom = ATOM_DATA.get(trans_name)
            if not atom: continue
            
            obs_cent = atom['wave'] * (1 + z)
            if obs_cent < self._x_ang[0] - 100 or obs_cent > self._x_ang[-1] + 100:
                continue

            tau_comp += self._voigt_optical_depth(
                self._x_ang, atom['wave'], atom['f'], atom['gamma'],
                z, N_col, b_eff
            )
        return tau_comp

    def _compute_model(self, p_free: np.ndarray) -> np.ndarray:
        p_full = self._constraints.map_p_free_to_full(p_free)
        
        # Initialize with static background if available, else 0
        if self._cached_tau_static is not None:
            tau_total = self._cached_tau_static.copy()
        else:
            tau_total = np.zeros_like(self._x_ang)
            
        components = self._system_list.components
        active_uuids = self._constraints._active_uuids
        num_params_per_comp = 4
        
        for i, comp in enumerate(components):
            # Optimization: Skip components that were pre-calculated
            # If we have a cache, we ONLY process Active components here
            if self._cached_tau_static is not None:
                if active_uuids is not None and comp.uuid not in active_uuids:
                    continue
            
            idx = i * num_params_per_comp
            z = p_full[idx]
            logN = p_full[idx + 1]
            b_val = p_full[idx + 2]
            btur = p_full[idx + 3]
            
            tau_total += self._compute_tau_component(comp, z, logN, b_val, btur)

        flux_model = np.exp(-tau_total)
        if self._resol_sigma_pix:
            flux_model = gaussian_filter1d(flux_model, self._resol_sigma_pix)
        return flux_model

    def _residual_function(self, p_free: np.ndarray) -> np.ndarray:
        model = self._compute_model(p_free)
        resid = (self._y_norm - model) / self._dy_norm
        resid[~np.isfinite(resid)] = 0.0
        return resid[self._fit_mask]

    def _find_feature_limits(self, z_center: float, trans_name: str, z_window_kms: float = 100.0) -> Tuple[Optional[float], Optional[float]]:
        atom = ATOM_DATA.get(trans_name)
        if not atom: return None, None
        
        lambda_obs = atom['wave'] * (1 + z_center)
        c_kms = 299792.458
        
        dx_search = (z_window_kms / c_kms) * lambda_obs
        mask_search = (self._x_ang > lambda_obs - dx_search) & (self._x_ang < lambda_obs + dx_search)
        
        if np.sum(mask_search) < 3: return None, None
        
        x_search = self._x_ang[mask_search]
        y_search = self._y_norm[mask_search]
        
        min_idx = np.argmin(y_search)
        if y_search[min_idx] > 0.9: return None, None
        
        start_idx = min_idx
        end_idx = min_idx
        threshold = 0.9
        
        while start_idx > 0 and y_search[start_idx] < threshold:
            start_idx -= 1
        while end_idx < len(y_search) - 1 and y_search[end_idx] < threshold:
            end_idx += 1
            
        if end_idx - start_idx < 2: return None, None
        
        return x_search[max(0, start_idx - 1)], x_search[min(len(y_search)-1, end_idx + 1)]

    def _smart_guess(self, p0: np.ndarray, z_window_kms: float) -> np.ndarray:
        p_refined = p0.copy()
        p_full = self._constraints.map_p_free_to_full(p0)
        is_free_mask = self._constraints._param_map['is_free']
        components = self._system_list.components
        num_params = 4
        free_idx_map = np.cumsum(is_free_mask) - 1
        c_kms = 299792.458
        
        for i, comp in enumerate(components):
            base_idx = i * num_params
            if not is_free_mask[base_idx]: continue
            
            current_logN = p_full[base_idx + 1]
            current_b = p_full[base_idx + 2]
            if not (np.isclose(current_logN, 13.5) and np.isclose(current_b, 10.0)):
                continue

            if comp.series in STANDARD_MULTIPLETS: primary = STANDARD_MULTIPLETS[comp.series][0]
            elif comp.series in ATOM_DATA: primary = comp.series
            else: continue
            
            lam_min, lam_max = self._find_feature_limits(comp.z, primary, z_window_kms=z_window_kms)
            if lam_min is None: continue
            
            mask_feat = (self._x_ang >= lam_min) & (self._x_ang <= lam_max)
            x_feat = self._x_ang[mask_feat]
            y_feat = self._y_norm[mask_feat]
            
            atom = ATOM_DATA[primary]
            lambda_obs = atom['wave'] * (1 + comp.z)
            
            A = 1.0 - y_feat
            A = np.clip(A, 0.0, 1.0)
            if np.sum(A) == 0: continue

            centroid = np.average(x_feat, weights=A)
            new_z = (centroid / atom['wave']) - 1.0
            
            variance = np.average((x_feat - centroid)**2, weights=A)
            sigma_ang = np.sqrt(variance)
            b_guess = (sigma_ang / centroid) * c_kms * 1.414
            b_guess = np.clip(b_guess, 5.0, 60.0)
            
            peak_flux = np.min(y_feat)
            safe_peak = max(0.01, min(0.98, peak_flux))
            tau0 = -np.log(safe_peak)
            N_guess = (tau0 * b_guess) / (1.4974e-15 * atom['f'] * atom['wave'])
            logN_guess = np.clip(np.log10(N_guess), 12.0, 17.0)
            
            idx_z = free_idx_map[base_idx]
            idx_logN = free_idx_map[base_idx + 1]
            idx_b = free_idx_map[base_idx + 2]
            
            p_refined[idx_z] = new_z
            p_refined[idx_logN] = logN_guess
            p_refined[idx_b] = b_guess
            
            logging.info(f"  -> Smart Guess: z={new_z:.5f}, logN={logN_guess:.2f}, b={b_guess:.1f}")
            
        return p_refined

    def fit(self, max_nfev: int = 2000, method: str = 'trf', 
            z_window_kms: float = 20.0,
            verbose: int = 1) -> Tuple[SystemListV2, np.ndarray, Any]:
        
        logging.info(f"Starting Voigt Fit (z_window={z_window_kms} km/s)...")
        p0 = self._constraints.p_free_vector
        
        # --- 1. PRE-CALCULATE FROZEN BACKGROUND ---
        # If we have active components, calculate the tau contribution 
        # of everything ELSE once.
        active_uuids = self._constraints._active_uuids
        if active_uuids is not None:
            logging.info("Pre-calculating frozen background...")
            self._cached_tau_static = np.zeros_like(self._x_ang)
            for comp in self._system_list.components:
                if comp.uuid not in active_uuids:
                    # This component is frozen. Add it to static.
                    self._cached_tau_static += self._compute_tau_component(
                        comp, comp.z, comp.logN, comp.b, comp.btur
                    )
        else:
            self._cached_tau_static = None
        # -------------------------------------------
        
        # --- CLAMP BOUNDS ---
        lower_bounds, upper_bounds = self._constraints.get_bounds()
        p_full = self._constraints.map_p_free_to_full(p0)
        is_free_mask = self._constraints._param_map['is_free']
        free_idx_map = np.cumsum(is_free_mask) - 1
        components = self._system_list.components
        num_params = 4
        c_kms = 299792.458
        
        for i, comp in enumerate(components):
            base_idx = i * num_params
            if not is_free_mask[base_idx]: continue 
            
            if comp.series in STANDARD_MULTIPLETS: primary = STANDARD_MULTIPLETS[comp.series][0]
            elif comp.series in ATOM_DATA: primary = comp.series
            else: continue
            
            dz_window = (1.0 + comp.z) * (z_window_kms / c_kms)
            z_min_win = comp.z - dz_window
            z_max_win = comp.z + dz_window
            
            lam_min, lam_max = self._find_feature_limits(comp.z, primary, z_window_kms=z_window_kms)
            
            if lam_min is not None:
                atom = ATOM_DATA[primary]
                z_min_feat = (lam_min / atom['wave']) - 1.0
                z_max_feat = (lam_max / atom['wave']) - 1.0
                z_final_min = max(z_min_win, z_min_feat)
                z_final_max = min(z_max_win, z_max_feat)
            else:
                z_final_min = z_min_win
                z_final_max = z_max_win

            idx_z = free_idx_map[base_idx]
            current_z = p0[idx_z]
            
            lower_bounds[idx_z] = min(z_final_min, current_z - 1e-6)
            upper_bounds[idx_z] = max(z_final_max, current_z + 1e-6)
            
            logging.debug(f"  -> Clamped {comp.series}: [{lower_bounds[idx_z]:.5f}, {upper_bounds[idx_z]:.5f}]")
        
        bounds = (lower_bounds, upper_bounds)
        
        # Smart Guess
        p0 = self._smart_guess(p0, z_window_kms)

        res = least_squares(
            self._residual_function, p0, bounds=bounds, 
            method=method, loss='linear', max_nfev=max_nfev, 
            x_scale='jac', verbose=verbose
        )
        
        # Reset cache for final model calculation (which might involve everything)
        # Actually, we can keep it, but let's clear it to be safe for the final output
        self._cached_tau_static = None 
        
        p_free_errors = np.zeros_like(res.x)
        chi2 = 0; red_chi2 = 0
        try:
            n_data = np.sum(self._fit_mask) 
            n_param = len(p0)
            dof = max(1, n_data - n_param)
            residuals = res.fun 
            chi2 = np.sum(residuals**2)
            red_chi2 = chi2 / dof
            J = res.jac
            cov = inv(J.T @ J) 
            p_free_errors = np.sqrt(np.diag(cov))
        except Exception as e:
            logging.warning(f"Could not calculate parameter errors: {e}")
            p_free_errors[:] = np.nan

        p_fitted_full = self._constraints.map_p_free_to_full(res.x)
        p_full_errors = np.zeros_like(p_fitted_full)
        p_full_errors[is_free_mask] = p_free_errors
        
        new_components = []
        old_components = self._system_list.components
        
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
        
        # Generate FULL model for the final output
        final_model_norm = self._compute_model(res.x)
        final_model_flux = final_model_norm * self._cont

        logging.info(f"Fit complete. Chi2: {chi2:.2f} (Red: {red_chi2:.2f}) over {n_data} pixels.")
        return new_system_list, final_model_flux, res

    def compute_model_flux(self) -> Tuple[np.ndarray, np.ndarray]:
        p_current_free = self._constraints.p_free_vector
        # No cache used here, calculate full model
        self._cached_tau_static = None 
        y_model_norm = self._compute_model(p_current_free)
        y_model_phys = y_model_norm * self._cont
        return self._x_ang, y_model_phys