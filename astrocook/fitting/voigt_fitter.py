import astropy.constants as const
import astropy.units as au
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares
from scipy.special import wofz
from scipy.linalg import pinv
import dataclasses
import logging
import numpy as np
from typing import List, Tuple, Dict, Any, Optional

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import ComponentDataV2
from astrocook.core.system_list import SystemListV2
from .voigt_model import VoigtModelConstraintV2

def convolve_flux(flux: np.ndarray, x_ang: np.ndarray, resol: Any, resol_unit: str = 'R') -> np.ndarray:
    r"""
    Convolves a flux array with a Gaussian kernel defined by resolution.

    Handles both constant resolution (scalar) and variable resolution (array)
    across the wavelength grid.

    Parameters
    ----------
    flux : np.ndarray
        The input flux array (either normalized 0-1 or physical units).
    x_ang : np.ndarray
        The wavelength grid in Angstroms. Must match ``flux`` length.
    resol : float or np.ndarray
        Resolution value(s).
        - If float: Constant resolution across the spectrum.
        - If array: Pixel-by-pixel resolution (must match ``flux`` length).
    resol_unit : str, optional
        Unit of the resolution value:
        - ``'R'``: Resolving Power (:math:`R = \lambda / \Delta\lambda`).
        - ``'km/s'``: FWHM velocity.
        Defaults to ``'R'``.

    Returns
    -------
    np.ndarray
        The convolved flux array with the same shape as input.
    """
    if len(flux) == 0: return flux
    if len(flux) != len(x_ang):
        raise ValueError(f"Shape mismatch: flux {len(flux)} vs wave {len(x_ang)}")

    # Estimate sampling (dx)
    if len(x_ang) > 1:
        dx = np.nanmedian(np.diff(x_ang))
    else:
        dx = 1.0
    
    if dx <= 0: return flux

    c_kms = 299792.458

    def get_sigma_pix(res_val, wav):
        if resol_unit == 'km/s':
            # FWHM_lambda = (v / c) * lambda
            # Sigma = FWHM / 2.355
            sigma_ang = (res_val / c_kms) * wav / 2.35482
        else:
            # FWHM_lambda = lambda / R
            sigma_ang = (wav / res_val) / 2.35482
        return sigma_ang / dx

    # 1. Variable Resolution (Array)
    if isinstance(resol, (np.ndarray, list)) and len(resol) == len(flux):
        resol_arr = np.asarray(resol)
        final_convolved = np.zeros_like(flux)
        
        if resol_unit == 'R':
            rounded_resol = np.round(resol_arr / 100.0) * 100.0
        else:
            rounded_resol = np.round(resol_arr * 2.0) / 2.0
            
        unique_resols = np.unique(rounded_resol)
        
        for r_val in unique_resols:
            if r_val <= 0: continue
            
            mask_r = np.isclose(resol_arr, r_val, atol=1e-3)
            
            # Use median wavelength of the masked region for kernel width approximation
            local_wave = np.median(x_ang[mask_r])
            sigma_pix = get_sigma_pix(r_val, local_wave)
            
            if sigma_pix > 0.1:
                # Convolve full array to handle edges, then apply mask
                temp_conv = gaussian_filter1d(flux, sigma_pix)
                final_convolved[mask_r] = temp_conv[mask_r]
            else:
                final_convolved[mask_r] = flux[mask_r]

        # Fill gaps (where resol <= 0)
        mask_no_res = (resol_arr <= 0)
        if np.any(mask_no_res):
            final_convolved[mask_no_res] = flux[mask_no_res]
            
        return final_convolved

    # 2. Constant Resolution (Scalar)
    elif isinstance(resol, (int, float)) and resol > 0:
        local_wave = np.median(x_ang)
        sigma_pix = get_sigma_pix(resol, local_wave)
        if sigma_pix > 0.1:
            return gaussian_filter1d(flux, sigma_pix)
            
    return flux

class VoigtFitterV2:
    """
    Engine for fitting Voigt profiles to spectral data.

    This class orchestrates the optimization process. It handles:

    1.  **Resolution Logic**: determining effective resolution (R or FWHM) from columns or metadata.
    2.  **Dynamic Masking**: identifying relevant pixels for the fit to speed up computation.
    3.  **Optimization**: utilizing ``scipy.optimize.least_squares`` to minimize residuals.
    4.  **Vectorization**: calculating optical depths for multiple lines simultaneously.

    Parameters
    ----------
    spectrum : SpectrumV2
        The spectral data to fit. Must have a valid normalization (continuum) or be pre-normalized.
    system_list : SystemListV2
        The list of absorption components to model.
    """
    def __init__(self, spectrum: SpectrumV2, system_list: SystemListV2, check_stop=None):
        self._spectrum = spectrum
        self._system_list = system_list
        self._constraints = system_list.constraint_model
        self._check_stop = check_stop

        # 1. Prepare Data Arrays
        # Convert to Angstroms for internal Voigt calc
        self._x_ang = self._spectrum.x.to(au.Angstrom).value
        
        norm_q = self._spectrum.norm
        
        if norm_q is None:
            # Raise explicit error so Recipe can catch it and trigger estimate_auto
            raise ValueError("Normalization vector is missing (Flux is physical and no continuum found).")
            
        self._norm = norm_q.value
        # Mask invalid normalization (zeros or negative values)
        self._norm[self._norm <= 0] = np.nan 

        # Normalize Flux and Error
        self._y_norm = self._spectrum.y.value / self._norm
        
        raw_dy = self._spectrum.dy.value
        min_err = 1e-9 * np.nanmedian(self._norm) 
        safe_dy = np.maximum(raw_dy, min_err)
        self._dy_norm = safe_dy / self._norm
        
        self._valid_mask = np.isfinite(self._y_norm) & np.isfinite(self._dy_norm) & (self._dy_norm > 0)
        
        # 2. Resolution Setup
        self._resol_column = None
        self._use_variable_resolution = False
        self._variable_resol_unit = 'R' # Default
        
        if self._spectrum.has_aux_column('resol'):
            self._resol_column = self._spectrum.get_column('resol').value
            # Check if column has valid data (not just NaNs/zeros)
            if np.any(np.isfinite(self._resol_column) & (self._resol_column > 0)):
                self._use_variable_resolution = True
                
                # Heuristic for unit
                med_res = np.nanmedian(self._resol_column)
                if med_res > 500.0:
                    self._variable_resol_unit = 'R'
                    logging.debug(f"VoigtFitter: Variable resolution detected. Median={med_res:.0f} -> Assuming 'R'.")
                else:
                    self._variable_resol_unit = 'km/s'
                    logging.debug(f"VoigtFitter: Variable resolution detected. Median={med_res:.1f} -> Assuming 'km/s'.")

        # Fallback to global metadata
        self._global_resol_val, self._global_sigma_pix = self._determine_resolution()
        
        # Determine global unit based on value magnitude (heuristic)
        self._global_resol_unit = 'R'
        if self._global_resol_val and self._global_resol_val < 500.0:
             self._global_resol_unit = 'km/s'

        if not self._use_variable_resolution:
            if self._global_sigma_pix:
                logging.debug(f"VoigtFitter: Global R={self._global_resol_val:.0f} ({self._global_resol_unit})")
            else:
                logging.debug("VoigtFitter: No resolution found. Fitting unconvolved profiles.")

        # Placeholders
        self._fit_mask = None
        self._sl = slice(None)
        self._x_calc = None
        self._y_norm_calc = None
        self._dy_norm_calc = None
        self._fit_mask_calc = None
        self._resol_calc = None
        
        # Cache is now a dict: { (resol_val, unit_str): tau_array }
        self._cached_tau_groups = {}
        
        # [OPTIMIZATION] Cache list of active components for the inner loop
        # Stores tuples: (index_in_system_list, component_object)
        self._active_fit_components: Optional[List[Tuple[int, ComponentDataV2]]] = None


    def _determine_resolution(self) -> Tuple[float, Optional[float]]:
        """
        Robustly finds resolution from Columns, Attributes, or Metadata.
        Returns: (R_value, sigma_pix)
        """
        # 1. Try Column ('resol') - Handled in __init__, but fallback logic here
        # (We already checked column in init, so here we check global props)

        # 2. Try Attribute (spec.resol)
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

    def _calc_sigma_pix_from_fwhm_kms(self, fwhm_kms: float, local_wave: float = None) -> float:
        """
        Convert FWHM (km/s) to Gaussian Sigma in pixels.
        If local_wave is provided, calculates for that specific wavelength.
        """
        wave = local_wave if local_wave is not None else np.median(self._x_ang)
        c_kms = 299792.458
        fwhm_ang = (fwhm_kms / c_kms) * wave
        sigma_ang = fwhm_ang / 2.35482
        
        # Estimate local sampling (dx)
        if len(self._x_ang) > 1:
            if local_wave is not None:
                idx = np.searchsorted(self._x_ang, local_wave)
                idx = np.clip(idx, 0, len(self._x_ang)-2)
                dx = self._x_ang[idx+1] - self._x_ang[idx]
            else:
                dx = np.mean(np.diff(self._x_ang))
        else:
            dx = 1.0
            
        return sigma_ang / dx

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
            
            # --- CRITICAL OPTIMIZATION Check ---
            # Even with pre-filtering, keep this for safety in the inner loop
            # Check if line center is vaguely in window
            obs_cent = atom['wave'] * (1 + z)
            if obs_cent < self._x_calc[0] - 100 or obs_cent > self._x_calc[-1] + 100:
                continue

            tau_comp += self._voigt_optical_depth(
                self._x_calc, atom['wave'], atom['f'], atom['gamma'],
                z, N_col, b_eff
            )
        return tau_comp
    
    def _get_component_resolution(self, comp: ComponentDataV2) -> Tuple[float, str]:
        """
        Determines the resolution (val, unit) for a specific component.
        Priority: Variable Column > Component Attribute > Global Default
        """
        # 1. Variable Column (Overrides everything if active)
        if self._use_variable_resolution:
            return -1.0, 'column' # Sentinel value for "Use Column"

        # 2. Component Attribute
        if comp.resol is not None and comp.resol > 0:
            # Heuristic for unit
            unit = 'R' if comp.resol > 500.0 else 'km/s'
            return float(comp.resol), unit
            
        # 3. Global Default
        if self._global_resol_val > 0:
            return self._global_resol_val, self._global_resol_unit
            
        return 0.0, 'none'

    def _compute_model(self, p_free: np.ndarray) -> np.ndarray:
        """
        Computes the normalized flux model.
        Uses vectorization if available (during fit), otherwise falls back to 
        standard iteration (during plotting).
        """
        p_full = self._constraints.map_p_free_to_full(p_free)
        
        # 1. Load Background (Tau computed from static components)
        tau_total_dict = {k: v.copy() for k, v in self._cached_tau_groups.items()}
        
        # 2. Add Active/Dynamic Components
        # CASE A: Vectorized (Fastest - Used inside fit())
        if hasattr(self, '_active_lines_data') and self._active_lines_data:
            for res_key, data in self._active_lines_data.items():
                # A. Extract parameters using fancy indexing
                z_arr = p_full[data['idx_z']]
                logN_arr = p_full[data['idx_N']]
                b_arr = p_full[data['idx_b']]
                btur_arr = p_full[data['idx_btur']]
                
                # B. Calculate Effective b
                b_eff_arr = np.sqrt(b_arr**2 + btur_arr**2)
                N_arr = 10**logN_arr
                
                # C. Compute Tau Vectorized
                tau_active = self._voigt_optical_depth_vectorized(
                    self._x_calc, 
                    data['lambda'], data['f'], data['gamma'],
                    z_arr, N_arr, b_eff_arr
                )
                
                # D. Add to accumulator
                if res_key not in tau_total_dict:
                    tau_total_dict[res_key] = tau_active
                else:
                    tau_total_dict[res_key] += tau_active

        # CASE B: Standard Loop (Fallback - Used for plotting / compute_model_flux)
        # This handles cases where we aren't using the pre-compiled vector cache.
        else:
            num_params_per_comp = 4
            
            # If active_uuids is None (default for plotting), iterate ALL components
            active_uuids = self._constraints._active_uuids
            
            for i, comp in enumerate(self._system_list.components):
                # Skip if we are in a partial-fit mode and this component isn't active
                if active_uuids is not None and comp.uuid not in active_uuids:
                    continue

                # Identify Resolution Group
                res_key = self._get_component_resolution(comp)
                if res_key not in tau_total_dict:
                    tau_total_dict[res_key] = np.zeros_like(self._x_calc)
                
                # Compute parameters
                idx = i * num_params_per_comp
                z = p_full[idx]
                logN = p_full[idx + 1]
                b_val = p_full[idx + 2]
                btur = p_full[idx + 3]
                
                # Add to group Tau
                tau_total_dict[res_key] += self._compute_tau_component(comp, z, logN, b_val, btur)

        # 3. Convolve and Combine
        final_model = np.ones_like(self._x_calc)
        
        for (res_val, res_unit), tau_arr in tau_total_dict.items():
            flux_group = np.exp(-tau_arr)
            
            if res_unit == 'column' and self._resol_calc is not None:
                 flux_group = convolve_flux(flux_group, self._x_calc, self._resol_calc, resol_unit=self._variable_resol_unit)
            elif res_val > 0:
                 flux_group = convolve_flux(flux_group, self._x_calc, res_val, resol_unit=res_unit)
            
            final_model *= flux_group
            
        return final_model

    def _residual_function(self, p_free: np.ndarray) -> np.ndarray:
        # Check for stop request
        if self._check_stop and self._check_stop():
            raise RuntimeError("Optimization stopped by user.")
        
        model = self._compute_model(p_free)
        resid = (self._y_norm_calc - model) / self._dy_norm_calc
        
        # Handle NaNs/Infs
        resid[~np.isfinite(resid)] = 0.0
        
        final_resid = resid[self._fit_mask_calc]
        
        # Optional: Soft penalty for b > 50 km/s to prevent runaways
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
        idx_max = min(len(self._x_ang) - 1, start_idx + right_boundary + 1)
        
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
            p_refined[idx_z] = new_z
            
            if is_free_mask[base_idx + 1]:
                idx_logN = free_idx_map[base_idx + 1]
                p_refined[idx_logN] = logN_guess
                
            if is_free_mask[base_idx + 2]:
                idx_b = free_idx_map[base_idx + 2]
                p_refined[idx_b] = b_guess
            logging.info(f"Smart Guess: z={new_z:.5f}, logN={logN_guess:.2f}, b={b_guess:.1f}")
        return p_refined

    def _voigt_optical_depth_vectorized(self, wave_grid_ang: np.ndarray, 
                                      lambda_0: np.ndarray, f_val: np.ndarray, 
                                      gamma: np.ndarray, z: np.ndarray, 
                                      N_col: np.ndarray, b_kms: np.ndarray) -> np.ndarray:
        """
        Vectorized Voigt profile calculation.
        Computes profiles for M lines on N pixels simultaneously.
        
        Input arrays (lambda_0, z, etc.) must be shape (M,).
        wave_grid_ang must be shape (N,).
        Returns summed tau array of shape (N,).
        """
        # Broadcast dimensions: (Lines, Pixels)
        # lambda_0: (M, 1)
        # wave_grid: (1, N)
        
        lambda_0 = lambda_0[:, np.newaxis]
        f_val = f_val[:, np.newaxis]
        gamma = gamma[:, np.newaxis]
        z = z[:, np.newaxis]
        N_col = N_col[:, np.newaxis]
        b_safe = np.maximum(b_kms, 0.1)[:, np.newaxis]
        
        wave_grid = wave_grid_ang[np.newaxis, :]
        
        # 1. Physics
        lambda_c = lambda_0 * (1.0 + z)
        c_kms = 2.99792458e5
        dop_width_ang = (b_safe / c_kms) * lambda_c
        dop_width_ang = np.maximum(dop_width_ang, 1e-10)

        # 2. Coordinate Transformation
        # x shape: (M, N)
        x = (wave_grid - lambda_c) / dop_width_ang
        
        c_ang_s = 2.99792458e18
        a = (gamma * lambda_c**2) / (4.0 * np.pi * c_ang_s * dop_width_ang)
        
        # 3. Faddeeva (Voigt) - Expensive Step
        # Computes entire matrix in C-speed
        z_complex = x + 1j * a
        H_ax = wofz(z_complex).real
        
        # 4. Optical Depth
        prefactor = 1.4974e-15
        tau_matrix = prefactor * N_col * f_val * lambda_0 / b_safe * H_ax
        
        # 5. Collapse lines (Sum over axis 0)
        return np.sum(tau_matrix, axis=0)

    def fit(self, max_nfev: int = 2000, method: str = 'trf', 
            z_window_kms: float = 20.0, verbose: int = 0) -> Tuple[SystemListV2, np.ndarray, Any]:
        """
        Execute the Voigt profile optimization.

        Performs a least-squares fit of the active components to the spectrum.
        Includes automatic spatial pre-filtering (fitting only pixels near the lines)
        and vectorized optical depth calculation.

        Parameters
        ----------
        max_nfev : int, optional
            Maximum number of function evaluations for the optimizer. Defaults to ``2000``.
        method : str, optional
            Optimization method (e.g., ``'trf'``, ``'lm'``). Defaults to ``'trf'``.
        z_window_kms : float, optional
            Velocity window (km/s) around component centers to include in the fit mask.
            Defaults to ``20.0``.
        verbose : int, optional
            Verbosity level for ``scipy.optimize.least_squares``. Defaults to ``0``.

        Returns
        -------
        tuple
            A tuple containing:
            
            1. **new_system_list** (:class:`~astrocook.core.system_list.SystemListV2`):
               A new system list containing the optimized parameter values and errors.
            2. **final_model_flux** (*np.ndarray*):
               The full model flux array (physical units), evaluated on the entire grid.
            3. **res** (*OptimizeResult*):
               The raw result object from ``scipy.optimize``.
        """
        logging.debug(f"Starting Voigt Fit (Window={z_window_kms} km/s)...")
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
        
        # Optimization: Slice the arrays to only the relevant fitting regions
        if np.any(self._fit_mask):
            indices = np.where(self._fit_mask)[0]
            self._sl = slice(max(0, np.min(indices) - 50), min(len(self._x_ang), np.max(indices) + 51))
        else:
            self._sl = slice(None)

        # Slice Data
        self._x_calc = self._x_ang[self._sl]
        self._y_norm_calc = self._y_norm[self._sl]
        self._dy_norm_calc = self._dy_norm[self._sl]
        self._fit_mask_calc = self._fit_mask[self._sl]
        
        # Slice Resolution (Crucial for Stitching)
        if self._use_variable_resolution and self._resol_column is not None:
            self._resol_calc = self._resol_column[self._sl]
        else:
            self._resol_calc = None

        logging.debug(f"Fit Mask: {np.sum(self._fit_mask)} pixels (Slice: {len(self._x_calc)}). Variable Resol: {self._use_variable_resolution}")

        # --- OPTIMIZATION 1: PRE-CALCULATE BACKGROUND ON RELEVANT COMPONENTS ONLY ---
        self._cached_tau_groups = {}
        
        # Get window bounds for spatial filtering (with buffer)
        if len(self._x_calc) > 0:
            x_min_win = self._x_calc[0]
            x_max_win = self._x_calc[-1]
        else:
            x_min_win, x_max_win = -np.inf, np.inf
        buffer_ang = 100.0
        
        relevant_static_comps = []
        if active_uuids is not None:
            # We filter the 1000+ static components to just the ~5-10 nearby ones
            for comp in self._system_list.components:
                if comp.uuid in active_uuids: continue
                
                is_relevant = False
                t_list = get_trans_list(comp.series)
                for t in t_list:
                    atom = ATOM_DATA.get(t)
                    if not atom: continue
                    obs_w = atom['wave'] * (1 + comp.z)
                    # Simple bounds check (much faster than computing profile)
                    if (x_min_win - buffer_ang) < obs_w < (x_max_win + buffer_ang):
                        is_relevant = True
                        break
                
                if is_relevant:
                    relevant_static_comps.append(comp)
        
        # Now compute background only for the survivors
        for comp in relevant_static_comps:
            # Determine Group
            res_key = self._get_component_resolution(comp)
            
            if res_key not in self._cached_tau_groups:
                self._cached_tau_groups[res_key] = np.zeros_like(self._x_calc)
            
            # Add Tau
            self._cached_tau_groups[res_key] += self._compute_tau_component(comp, comp.z, comp.logN, comp.b, comp.btur)

        # --- OPTIMIZATION 2: CACHE ACTIVE LIST FOR INNER LOOP ---
        self._active_fit_components = []
        if active_uuids is not None:
            for i, comp in enumerate(self._system_list.components):
                if comp.uuid in active_uuids:
                    self._active_fit_components.append((i, comp))
        else:
            # If fitting everything, we don't need a special cache
            self._active_fit_components = None

        # --- OPTIMIZATION 3: PREPARE VECTORIZED LINE DATA ---
        # We flatten components -> lines for the active set
        self._active_lines_data = {} # Key: res_key, Value: dict of arrays
        
        # Determine iterator (either the cached active list or full list)
        if self._active_fit_components is not None:
            iterator = self._active_fit_components
        else:
            iterator = enumerate(self._system_list.components)

        grouped_lines = {}
        
        for idx, comp in iterator:
            res_key = self._get_component_resolution(comp)
            if res_key not in grouped_lines:
                grouped_lines[res_key] = {
                    'lambda': [], 'f': [], 'gamma': [], 
                    'p_idx_z': [], 'p_idx_logN': [], 'p_idx_b': [], 'p_idx_btur': [] 
                }
            
            # Get transitions
            if comp.series in STANDARD_MULTIPLETS: t_list = STANDARD_MULTIPLETS[comp.series]
            elif comp.series in ATOM_DATA: t_list = [comp.series]
            else: continue
            
            # Store Mapping to p_full indices
            base_idx = idx * 4
            
            for t in t_list:
                atom = ATOM_DATA.get(t)
                if not atom: continue
                
                # Append Atomic Data
                g = grouped_lines[res_key]
                g['lambda'].append(atom['wave'])
                g['f'].append(atom['f'])
                g['gamma'].append(atom['gamma'])
                
                # Append Parameter Indices (So we can pull values from p_full rapidly)
                g['p_idx_z'].append(base_idx)
                g['p_idx_logN'].append(base_idx + 1)
                g['p_idx_b'].append(base_idx + 2)
                g['p_idx_btur'].append(base_idx + 3)

        # Convert to Numpy Arrays and store
        for rk, data in grouped_lines.items():
            if not data['lambda']: continue
            self._active_lines_data[rk] = {
                'lambda': np.array(data['lambda']),
                'f':      np.array(data['f']),
                'gamma':  np.array(data['gamma']),
                'idx_z':  np.array(data['p_idx_z']),
                'idx_N':  np.array(data['p_idx_logN']),
                'idx_b':  np.array(data['p_idx_b']),
                'idx_btur': np.array(data['p_idx_btur'])
            }

        # --- Bounds Logic ---
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

        # Clamp p0 to bounds to prevent "Initial guess outside bounds" error
        # Because smart_guess might have proposed values that violate strict constraints.
        p0 = np.clip(p0, lower_bounds, upper_bounds)

        # Fit
        res = least_squares(self._residual_function, p0, bounds=bounds, method=method, loss='linear', max_nfev=max_nfev, x_scale='jac', verbose=verbose)
        
        # --- CLEANUP (Reset Optimized Caches) ---
        self._active_fit_components = None
        self._active_lines_data = None

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
        
        # We try to pick the most representative resolution used in this fit.
        fit_resol_val = 0.0
        if self._use_variable_resolution and self._resol_calc is not None:
             fit_resol_val = np.nanmedian(self._resol_calc)
             if self._variable_resol_unit == 'km/s': fit_resol_val = 299792.458 / fit_resol_val # Store as R? Or keep native? Let's keep native value.
        elif self._global_resol_val > 0:
             fit_resol_val = self._global_resol_val

        new_components = []
        old_components = self._system_list.components
        active_uuids = self._constraints._active_uuids

        for i, old_comp in enumerate(old_components):
            idx = i * 4
            
            dz_new = p_full_errors[idx] if is_free_mask[idx] else old_comp.dz
            dlogN_new = p_full_errors[idx+1] if is_free_mask[idx+1] else old_comp.dlogN
            db_new = p_full_errors[idx+2] if is_free_mask[idx+2] else old_comp.db
            
            # Metadata Update
            comp_chi2 = old_comp.chi2
            comp_resol = old_comp.resol
            
            if active_uuids is None or old_comp.uuid in active_uuids:
                comp_chi2 = red_chi2
                
                # Logic: If we used the component's OWN resolution, preserve it.
                # If we used the global/column resolution, update it to reflect the fit conditions.
                used_res_val, used_res_unit = self._get_component_resolution(old_comp)
                
                if used_res_unit == 'column':
                    # If variable, we usually store the median R of the fit
                    comp_resol = fit_resol_val
                elif used_res_val > 0:
                    # If we used a specific scalar (Global or Attribute), store that.
                    comp_resol = used_res_val
            
            new_comp = dataclasses.replace(old_comp,
                z=p_fitted_full[idx], 
                logN=p_fitted_full[idx+1], 
                b=p_fitted_full[idx+2], 
                btur=p_fitted_full[idx+3],
                dz=dz_new, 
                dlogN=dlogN_new, 
                db=db_new,
                chi2=comp_chi2,
                resol=comp_resol
            )
            new_components.append(new_comp)
            
        new_data_core = dataclasses.replace(self._system_list._data, components=new_components)
        new_system_list = SystemListV2(new_data_core)
        
        # 4. Final Flux Model (Un-sliced)
        self._cached_tau_groups = {} # Clear cache

        # [FIX] Temporarily disable active_uuids so _compute_model calculates EVERYTHING
        # We need the full model (static + active) for the green line.
        stored_active_uuids = self._constraints._active_uuids
        self._constraints._active_uuids = None
        self._active_fit_components = None 
        self._active_lines_data = None # Ensure we iterate everything
        
        # Restore full arrays
        temp_x = self._x_calc
        temp_resol = self._resol_calc
        self._x_calc = self._x_ang
        if self._use_variable_resolution:
             self._resol_calc = self._resol_column
        else:
             self._resol_calc = None
             
        # This will now loop over ALL components because active_uuids is None
        final_model_norm = self._compute_model(res.x)
        
        # Restore state (good practice)
        self._constraints._active_uuids = stored_active_uuids
        self._x_calc = temp_x
        self._resol_calc = temp_resol
        
        final_model_flux = final_model_norm * self._norm

        logging.debug(f"Fit complete. Chi2: {chi2:.2f} (Red: {red_chi2:.2f}) over {n_data} pixels.")
        return new_system_list, final_model_flux, res

    def compute_model_flux(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the model flux for the current parameters.

        Calculates the Voigt profile model based on the current state of the
        system list, without performing any fitting/optimization.

        Returns
        -------
        tuple
            A tuple containing:
            1. **wavelength** (*np.ndarray*): The wavelength grid in Angstroms.
            2. **model_flux** (*np.ndarray*): The computed model flux (physical units).
        """
        p_current_free = self._constraints.p_free_vector
        temp_x = self._x_calc
        temp_resol = self._resol_calc
        
        # Switch to full spectrum context
        self._x_calc = self._x_ang
        if self._use_variable_resolution:
             self._resol_calc = self._resol_column
        else:
             self._resol_calc = None

        # Reset caches to force Case B (Standard Loop) in _compute_model
        self._cached_tau_groups = {}
        self._active_lines_data = None 
        
        y_model_norm = self._compute_model(p_current_free)
        
        # Restore context
        self._x_calc = temp_x 
        self._resol_calc = temp_resol
        
        y_model_phys = y_model_norm * self._norm
        return self._x_ang, y_model_phys