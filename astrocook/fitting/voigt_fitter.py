import astropy.constants as const
import astropy.units as au
from scipy.ndimage import gaussian_filter1d
from scipy.signal import argrelextrema
import dataclasses
import logging
import numpy as np
from scipy.optimize import least_squares
from scipy.special import wofz
from scipy.linalg import pinv
from typing import List, Tuple, Dict, Any, Optional

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import ComponentDataV2, SystemListDataV2
from astrocook.core.system_list import SystemListV2
from .voigt_model import VoigtModelConstraintV2

class VoigtFitterV2:
    """
    Engine for fitting Voigt profiles to spectral data.

    This class handles:
    1.  Resolution logic (determining R/FWHM from metadata).
    2.  Dynamic masking (fitting only relevant pixels).
    3.  Optimization (using ``scipy.optimize.least_squares``).
    4.  Error estimation (covariance matrix inversion).

    Parameters
    ----------
    spectrum : SpectrumV2
        The spectrum data to fit.
    system_list : SystemListV2
        The list of components to model.
    """
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
        
        # Clip dy
        raw_dy = self._spectrum.dy.value
        min_err = 1e-9 * np.nanmedian(self._cont) 
        safe_dy = np.maximum(raw_dy, min_err)
        self._dy_norm = safe_dy / self._cont
        
        self._valid_mask = np.isfinite(self._y_norm) & np.isfinite(self._dy_norm) & (self._dy_norm > 0)
        
        # Placeholders
        self._fit_mask = None
        self._sl = slice(None)
        self._x_calc = None
        self._y_norm_calc = None
        self._dy_norm_calc = None
        self._fit_mask_calc = None
        self._cached_tau_static = None

        # --- RESOLUTION LOGIC ---
        # We determine R and sigma_pix here.
        # This allows us to use the exact same logic for fitting and saving.
        self._resol_val, self._resol_sigma_pix = self._determine_resolution()
        
        if self._resol_sigma_pix:
            logging.info(f"VoigtFitter: Resolution R={self._resol_val:.0f} (sigma_pix={self._resol_sigma_pix:.3f})")
        else:
            logging.info("VoigtFitter: No resolution found. Fitting unconvolved profiles.")

    def _determine_resolution(self) -> Tuple[float, Optional[float]]:
        """
        Robustly finds resolution from Columns, Attributes, or Metadata.
        Returns: (R_value, sigma_pix)
        """
        # 1. Try Column ('resol')
        if self._spectrum.has_aux_column('resol'):
            col = self._spectrum.get_column('resol')
            if col is not None:
                val = np.nanmedian(col.value)
                if np.isfinite(val) and val > 0:
                    return val, self._calc_sigma_pix_from_R(val)

        # 2. Try Attribute (spec.resol) - covers SpectrumDataV2 field
        # Check wrapper property first, then data core
        val = getattr(self._spectrum, 'resol', 0.0)
        if (not np.isfinite(val) or val <= 0) and hasattr(self._spectrum, '_data'):
            val = getattr(self._spectrum._data, 'resol', 0.0)
        
        if np.isfinite(val) and val > 0:
            return val, self._calc_sigma_pix_from_R(val)

        # 3. Try Metadata ('resol' or 'resol_fwhm')
        meta = self._spectrum.meta
        
        # 3a. 'resol' (R)
        if 'resol' in meta:
            val = float(meta['resol'])
            if np.isfinite(val) and val > 0:
                return val, self._calc_sigma_pix_from_R(val)
        
        # 3b. 'resol_fwhm' (Velocity km/s)
        if 'resol_fwhm' in meta:
            val_kms = float(meta['resol_fwhm'])
            if np.isfinite(val_kms) and val_kms > 0:
                # Convert FWHM km/s to R approx for return value
                # R = c / FWHM
                c_kms = 299792.458
                r_approx = c_kms / val_kms
                return r_approx, self._calc_sigma_pix_from_fwhm_kms(val_kms)

        return 0.0, None

    def _calc_sigma_pix_from_R(self, R: float) -> float:
        """Helper: Convert R (lambda/d_lambda) to Gaussian Sigma in pixels."""
        mid_wave = np.median(self._x_ang)
        fwhm_ang = mid_wave / R
        sigma_ang = fwhm_ang / 2.35482
        dx_avg = np.mean(np.diff(self._x_ang))
        return sigma_ang / dx_avg

    def _calc_sigma_pix_from_fwhm_kms(self, fwhm_kms: float) -> float:
        """Helper: Convert FWHM (km/s) to Gaussian Sigma in pixels."""
        mid_wave = np.median(self._x_ang)
        c_kms = 299792.458
        fwhm_ang = (fwhm_kms / c_kms) * mid_wave
        sigma_ang = fwhm_ang / 2.35482
        dx_avg = np.mean(np.diff(self._x_ang))
        return sigma_ang / dx_avg

    def _voigt_optical_depth(self, wave_grid_ang: np.ndarray, 
                             lambda_0: float, f_val: float, gamma: float,
                             z: float, N_col: float, b_kms: float) -> np.ndarray:
        """Computes the optical depth tau(lambda) for a single Voigt line."""
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
        tau_comp = np.zeros_like(self._x_calc)
        b_eff = np.sqrt(b**2 + btur**2)
        N_col = 10**logN
        
        if comp.series in STANDARD_MULTIPLETS: trans_list = STANDARD_MULTIPLETS[comp.series]
        elif comp.series in ATOM_DATA: trans_list = [comp.series]
        elif ',' in comp.series: trans_list = comp.series.split(',')
        else: return tau_comp

        for trans_name in trans_list:
            trans_name = trans_name.strip()
            atom = ATOM_DATA.get(trans_name)
            if not atom: continue
            
            obs_cent = atom['wave'] * (1 + z)
            if obs_cent < self._x_calc[0] - 100 or obs_cent > self._x_calc[-1] + 100:
                continue

            tau_comp += self._voigt_optical_depth(
                self._x_calc, atom['wave'], atom['f'], atom['gamma'],
                z, N_col, b_eff
            )
        return tau_comp

    def _compute_model(self, p_free: np.ndarray) -> np.ndarray:
        """
        Computes the normalized flux model for the current free parameters.
        Includes convolution with resolution if applicable.
        """
        p_full = self._constraints.map_p_free_to_full(p_free)
        
        if self._cached_tau_static is not None:
            tau_total = self._cached_tau_static.copy()
        else:
            tau_total = np.zeros_like(self._x_calc)
            
        components = self._system_list.components
        active_uuids = self._constraints._active_uuids
        num_params_per_comp = 4
        
        for i, comp in enumerate(components):
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
        resid = (self._y_norm_calc - model) / self._dy_norm_calc
        resid[~np.isfinite(resid)] = 0.0
        final_resid = resid[self._fit_mask_calc]
        
        p_full = self._constraints.map_p_free_to_full(p_free)
        is_free_mask = self._constraints._param_map['is_free']
        num_params = 4
        penalty_list = []
        
        for i in range(len(p_full) // num_params):
            idx_b = i * num_params + 2
            if is_free_mask[idx_b]:
                b_val = p_full[idx_b]
                if b_val > 50.0:
                    penalty = 5.0 * (b_val - 50.0)
                    penalty_list.append(penalty)
                else:
                    penalty_list.append(0.0)
        
        if penalty_list:
            final_resid = np.concatenate((final_resid, penalty_list))

        return final_resid

    def _find_feature_limits(self, z_center: float, trans_name: str, z_window_kms: float = 150.0) -> Tuple[Optional[float], Optional[float]]:
        atom = ATOM_DATA.get(trans_name)
        if not atom: return None, None
        
        lambda_obs = atom['wave'] * (1 + z_center)
        c_kms = 299792.458
        
        dx_search = (z_window_kms / c_kms) * lambda_obs
        
        idx_center = np.searchsorted(self._x_ang, lambda_obs)
        
        if len(self._x_ang) > 10:
             sample_dx = np.mean(np.diff(self._x_ang[max(0, idx_center-10):min(len(self._x_ang), idx_center+10)]))
        else: sample_dx = 1.0
        
        px_width = int(dx_search / sample_dx) if sample_dx > 0 else 50
        
        start_idx = max(0, idx_center - px_width)
        end_idx = min(len(self._x_ang), idx_center + px_width)
        
        if end_idx - start_idx < 5: return None, None
        
        y_window = self._y_norm[start_idx:end_idx]
        dy_window = self._dy_norm[start_idx:end_idx]
        y_smooth = gaussian_filter1d(y_window, sigma=2.0)
        
        center_rel = idx_center - start_idx
        center_rel = max(0, min(len(y_window)-1, center_rel))

        # Scan Left
        left_boundary = 0
        lowest_point = y_smooth[center_rel]
        for i in range(center_rel, -1, -1):
            val = y_smooth[i]
            err = dy_window[i]
            if val < lowest_point: lowest_point = val
            is_local_max = (i > 0 and i < len(y_smooth)-1 and y_smooth[i] >= y_smooth[i-1] and y_smooth[i] >= y_smooth[i+1]) or (i == 0)
            if val > 0.9: left_boundary = i; break
            if is_local_max:
                if (val - lowest_point) > 4.0 * err: left_boundary = i; break
                
        # Scan Right
        right_boundary = len(y_window) - 1
        lowest_point = y_smooth[center_rel]
        for i in range(center_rel, len(y_window)):
            val = y_smooth[i]
            err = dy_window[i]
            if val < lowest_point: lowest_point = val
            is_local_max = (i > 0 and i < len(y_smooth)-1 and y_smooth[i] >= y_smooth[i-1] and y_smooth[i] >= y_smooth[i+1]) or (i == len(y_window)-1)
            if val > 0.9: right_boundary = i; break
            if is_local_max:
                if (val - lowest_point) > 4.0 * err: right_boundary = i; break

        idx_min = max(0, start_idx + left_boundary - 1)
        idx_max = min(len(self._x_ang), start_idx + right_boundary + 1)
        
        if idx_max <= idx_min: return None, None
        return self._x_ang[idx_min], self._x_ang[idx_max]

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
            if not (np.isclose(current_logN, 13.5) and np.isclose(current_b, 10.0)): continue

            if comp.series in STANDARD_MULTIPLETS: primary = STANDARD_MULTIPLETS[comp.series][0]
            elif comp.series in ATOM_DATA: primary = comp.series
            else: continue
            
            lam_min, lam_max = self._find_feature_limits(comp.z, primary, z_window_kms=z_window_kms)
            if lam_min is None: continue
            
            mask_feat = (self._x_ang >= lam_min) & (self._x_ang <= lam_max)
            x_feat = self._x_ang[mask_feat]
            y_feat = self._y_norm[mask_feat]
            atom = ATOM_DATA[primary]
            
            A = np.clip(1.0 - y_feat, 0.0, 1.0)
            if np.sum(A) == 0: continue

            centroid = np.average(x_feat, weights=A)
            new_z = (centroid / atom['wave']) - 1.0
            variance = np.average((x_feat - centroid)**2, weights=A)
            sigma_ang = np.sqrt(variance)
            
            b_guess = np.clip((sigma_ang / centroid) * c_kms * 1.414, 5.0, 25.0) 
            peak_flux = np.min(y_feat)
            safe_peak = max(0.01, min(0.98, peak_flux))
            tau0 = -np.log(safe_peak)
            N_guess = (tau0 * b_guess) / (1.4974e-15 * atom['f'] * atom['wave'])
            logN_guess = np.clip(np.log10(N_guess), 12.0, 15.0)
            
            idx_z = free_idx_map[base_idx]
            idx_logN = free_idx_map[base_idx + 1]
            idx_b = free_idx_map[base_idx + 2]
            p_refined[idx_z] = new_z
            p_refined[idx_logN] = logN_guess
            p_refined[idx_b] = b_guess
            logging.info(f"  -> Smart Guess: z={new_z:.5f}, logN={logN_guess:.2f}, b={b_guess:.1f}")
        return p_refined

    def fit(self, max_nfev: int = 2000, method: str = 'trf', 
            z_window_kms: float = 20.0, verbose: int = 1) -> Tuple[SystemListV2, np.ndarray, Any]:
        """
        Executes the optimization using scipy.optimize.least_squares.

        This method:
        1.  Determines the **Fit Mask**: Only pixels within ``z_window_kms`` of the
            target lines are used for chi-squared calculation.
        2.  Refines the initial guess using ``_smart_guess``.
        3.  Runs the optimizer.
        4.  Computes parameter errors from the covariance matrix.
        5.  Updates the components in the system list.

        Parameters
        ----------
        max_nfev : int
            Maximum number of function evaluations allowed.
        method : str
            Optimization algorithm ('trf', 'dogbox', 'lm').
        z_window_kms : float
            Velocity window (km/s) around line centers to include in the fit mask.
        verbose : int
            Verbosity level for the optimizer.

        Returns
        -------
        new_system_list : SystemListV2
            A new list containing the best-fit component parameters.
        final_model_flux : np.ndarray
            The flux array of the best-fit model (convolved and multiplied by continuum).
        res : OptimizeResult
            The raw result object from ``scipy.optimize``.
        """
        logging.info(f"Starting Voigt Fit (Window={z_window_kms} km/s)...")
        p0 = self._constraints.p_free_vector
        
        # 1. DYNAMIC FIT MASK
        self._fit_mask = np.zeros_like(self._valid_mask, dtype=bool)
        c_kms = 299792.458
        active_uuids = self._constraints._active_uuids
        
        def get_trans_list(series):
            if series in STANDARD_MULTIPLETS: return STANDARD_MULTIPLETS[series]
            if series in ATOM_DATA: return [series]
            return []

        active_components = [c for c in self._system_list.components if c.uuid in active_uuids] if active_uuids else self._system_list.components

        if not active_components:
            self._fit_mask = self._valid_mask
        else:
            for comp in active_components:
                for trans in get_trans_list(comp.series):
                    lam_min, lam_max = self._find_feature_limits(comp.z, trans, z_window_kms=z_window_kms)
                    if lam_min is not None:
                        mask_feat = (self._x_ang >= lam_min) & (self._x_ang <= lam_max)
                        self._fit_mask |= mask_feat
                    else:
                        lam_0 = ATOM_DATA.get(trans, {}).get('wave')
                        if lam_0:
                            lam_obs = lam_0 * (1.0 + comp.z)
                            dw = (z_window_kms / c_kms) * lam_obs
                            self._fit_mask |= (self._x_ang > lam_obs - dw) & (self._x_ang < lam_obs + dw)
            self._fit_mask &= self._valid_mask
        
        if np.any(self._fit_mask):
            indices = np.where(self._fit_mask)[0]
            self._sl = slice(max(0, np.min(indices) - 50), min(len(self._x_ang), np.max(indices) + 51))
        else:
            self._sl = slice(None)

        self._x_calc = self._x_ang[self._sl]
        self._y_norm_calc = self._y_norm[self._sl]
        self._dy_norm_calc = self._dy_norm[self._sl]
        self._fit_mask_calc = self._fit_mask[self._sl]
        logging.info(f"Fit Mask: {np.sum(self._fit_mask)} pixels (Slice: {len(self._x_calc)}).")

        # 2. PRE-CALCULATE BACKGROUND
        if active_uuids is not None:
            self._cached_tau_static = np.zeros_like(self._x_calc)
            for comp in self._system_list.components:
                if comp.uuid not in active_uuids:
                    self._cached_tau_static += self._compute_tau_component(comp, comp.z, comp.logN, comp.b, comp.btur)
        else:
            self._cached_tau_static = None

        # Bounds Logic
        lower_bounds, upper_bounds = self._constraints.get_bounds()
        p_full = self._constraints.map_p_free_to_full(p0)
        is_free_mask = self._constraints._param_map['is_free']
        free_idx_map = np.cumsum(is_free_mask) - 1
        
        for i, comp in enumerate(self._system_list.components):
            base_idx = i * 4
            if not is_free_mask[base_idx]: continue 
            
            if comp.series in STANDARD_MULTIPLETS: primary = STANDARD_MULTIPLETS[comp.series][0]
            elif comp.series in ATOM_DATA: primary = comp.series
            else: continue
            
            dz_window = (1.0 + comp.z) * (z_window_kms / c_kms)
            z_min_win, z_max_win = comp.z - dz_window, comp.z + dz_window
            lam_min, lam_max = self._find_feature_limits(comp.z, primary, z_window_kms=z_window_kms)
            
            z_final_min, z_final_max = z_min_win, z_max_win
            if lam_min is not None:
                atom = ATOM_DATA[primary]
                z_final_min = max(z_min_win, (lam_min / atom['wave']) - 1.0)
                z_final_max = min(z_max_win, (lam_max / atom['wave']) - 1.0)

            idx_z = free_idx_map[base_idx]
            cur = p0[idx_z]
            lower_bounds[idx_z] = min(z_final_min, cur - 1e-6)
            upper_bounds[idx_z] = max(z_final_max, cur + 1e-6)
        
        bounds = (lower_bounds, upper_bounds)
        p0 = self._smart_guess(p0, z_window_kms)

        # Fit
        res = least_squares(self._residual_function, p0, bounds=bounds, method=method, loss='linear', max_nfev=max_nfev, x_scale='jac', verbose=verbose)
        
        # Statistics
        p_free_errors = np.zeros_like(res.x)
        chi2 = 0; red_chi2 = 0
        try:
            residuals = res.fun 
            chi2 = np.sum(residuals**2)
            n_data = np.sum(self._fit_mask_calc)
            dof = max(1, n_data - len(p0))
            red_chi2 = chi2 / dof
            J = res.jac
            cov = pinv(J.T @ J) 
            p_free_errors = np.sqrt(np.diag(cov))
        except Exception as e:
            logging.warning(f"Error calc failed: {e}")
            p_free_errors[:] = np.nan

        p_fitted_full = self._constraints.map_p_free_to_full(res.x)
        p_full_errors = np.zeros_like(p_fitted_full)
        p_full_errors[is_free_mask] = p_free_errors
        
        # [FIX] Use the exact resolution calculated in __init__
        # This ensures the metadata saved matches the fit physics.
        fit_resol_val = self._resol_val

        new_components = []
        old_components = self._system_list.components
        active_uuids = self._constraints._active_uuids

        for i, old_comp in enumerate(old_components):
            idx = i * 4
            
            # 1. Update Parameters (always safe: fixed params just get their old value back)
            dz_new = p_full_errors[idx] if is_free_mask[idx] else old_comp.dz
            dlogN_new = p_full_errors[idx+1] if is_free_mask[idx+1] else old_comp.dlogN
            db_new = p_full_errors[idx+2] if is_free_mask[idx+2] else old_comp.db
            
            # 2. Determine Metadata (Chi2 / Resol)
            # Only apply new stats if the component was part of the optimization group.
            # If active_uuids is None, it means "Fit All", so everyone gets updated.
            if active_uuids is None or old_comp.uuid in active_uuids:
                comp_chi2 = red_chi2
                comp_resol = fit_resol_val
            else:
                # Keep existing metadata for untouched components
                comp_chi2 = old_comp.chi2
                comp_resol = old_comp.resol

            new_comp = dataclasses.replace(old_comp,
                z=p_fitted_full[idx], 
                logN=p_fitted_full[idx+1], 
                b=p_fitted_full[idx+2], 
                btur=p_fitted_full[idx+3],
                dz=dz_new, 
                dlogN=dlogN_new, 
                db=db_new,
                chi2=comp_chi2,   # [FIX] Conditional Update
                resol=comp_resol  # [FIX] Conditional Update
            )
            new_components.append(new_comp)
            
        new_data_core = dataclasses.replace(self._system_list._data, components=new_components)
        new_system_list = SystemListV2(new_data_core)
        
        self._cached_tau_static = None
        temp_x = self._x_calc
        self._x_calc = self._x_ang
        final_model_norm = self._compute_model(res.x)
        self._x_calc = temp_x
        final_model_flux = final_model_norm * self._cont

        logging.info(f"Fit complete. Chi2: {chi2:.2f} (Red: {red_chi2:.2f}) over {n_data} pixels.")
        return new_system_list, final_model_flux, res

    def compute_model_flux(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the model flux for the *current* parameters without running a fit.

        Returns
        -------
        x_ang : np.ndarray
            The wavelength grid in Angstroms.
        y_model : np.ndarray
            The model flux (physical units) over the full spectral range.
        """
        p_current_free = self._constraints.p_free_vector
        temp_x = self._x_calc
        self._x_calc = self._x_ang
        self._cached_tau_static = None
        y_model_norm = self._compute_model(p_current_free)
        self._x_calc = temp_x 
        y_model_phys = y_model_norm * self._cont
        return self._x_ang, y_model_phys