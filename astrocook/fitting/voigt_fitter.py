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
from astrocook.fitting.lsf import PixelLSF
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import ComponentDataV2
from astrocook.core.system_list import SystemListV2

#: A line is evaluated only where its optical depth exceeds this (#526)
TAU_CUT = 1e-6


def _transitions(series: str) -> List[str]:
    """Transitions of a series in the atomic data (multiplet, single line
    or comma-separated list)."""
    if series in STANDARD_MULTIPLETS:
        names = STANDARD_MULTIPLETS[series]
    elif series in ATOM_DATA:
        names = [series]
    elif ',' in series:
        names = series.split(',')
    else:
        names = []
    return [n.strip() for n in names if ATOM_DATA.get(n.strip())]


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
            if not r_val > 0: continue  # also skips NaN
            
            # Select by the ROUNDED value: selecting by the raw value left the
            # pixels with a non-round resolution at 0 (findings register F10)
            mask_r = rounded_resol == r_val
            
            # Use median wavelength of the masked region for kernel width approximation
            local_wave = np.median(x_ang[mask_r])
            sigma_pix = get_sigma_pix(r_val, local_wave)
            
            if sigma_pix > 0.1:
                # Convolve full array to handle edges, then apply mask
                temp_conv = gaussian_filter1d(flux, sigma_pix)
                final_convolved[mask_r] = temp_conv[mask_r]
            else:
                final_convolved[mask_r] = flux[mask_r]

        # Fill gaps (where resol <= 0 or NaN)
        mask_no_res = ~(rounded_resol > 0)
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
        # Pixel edges, for the exact convolution with the LSF (#526)
        self._xmin_ang = self._spectrum.xmin.to(au.Angstrom).value
        self._xmax_ang = self._spectrum.xmax.to(au.Angstrom).value
        self._x_tau = None
        self._lsf = None
        
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
        
        # 2. Resolution Setup: the 'resol' column of the spectrum, as R
        # (#526; before, five sources and a guess of the unit by value). A
        # flat column takes the scalar path, as a global value did before
        self._resol_column = None
        self._use_variable_resolution = False
        self._variable_resol_unit = 'R'
        self._global_resol_val, self._global_sigma_pix = 0.0, None
        self._global_resol_unit = 'R'
        r_pix = self._spectrum.resolution()
        finite = r_pix[np.isfinite(r_pix)]
        if finite.size == r_pix.size and finite.size and np.all(finite == finite[0]):
            self._global_resol_val = float(finite[0])
            self._global_sigma_pix = self._calc_sigma_pix_from_R(self._global_resol_val)
            logging.debug(f"VoigtFitter: Global R={self._global_resol_val:.0f}")
        elif finite.size:
            self._resol_column = r_pix
            self._use_variable_resolution = True
            logging.debug(f"VoigtFitter: Variable resolution, median R={np.median(finite):.0f}")
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

    def _compute_tau_component(self, comp: ComponentDataV2, z: float, logN: float, b: float, btur: float) -> np.ndarray:
        atoms = [ATOM_DATA[t] for t in _transitions(comp.series)]
        if not atoms:
            return np.zeros_like(self._x_tau)
        n = len(atoms)
        return self._voigt_optical_depth_vectorized(
            self._x_tau,
            np.array([at['wave'] for at in atoms]), np.array([at['f'] for at in atoms]),
            np.array([at['gamma'] for at in atoms]), np.full(n, float(z)),
            np.full(n, 10**logN), np.full(n, np.sqrt(b**2 + btur**2)))

    def _get_component_resolution(self, comp: ComponentDataV2) -> Tuple[float, str]:
        """
        Determines the resolution (val, unit) for a specific component.
        Priority: Variable Column > Component Attribute > Global Default
        """
        # 1. Component Attribute (Highest Priority: User manual override),
        # always R (#526; legacy FWHM values are converted at load)
        if comp.resol is not None and comp.resol > 0:
            return float(comp.resol), 'R'

        # 2. Variable Column (Instrument baseline)
        if self._use_variable_resolution:
            return -1.0, 'column' # Sentinel value for "Use Column"
            
        # 3. Global Default (Fallback)
        if self._global_resol_val > 0:
            return self._global_resol_val, self._global_resol_unit
            
        return 0.0, 'none'

    def _set_lsf(self, sl: slice, cache: bool = True) -> None:
        """Fine grid and LSF for the pixels ``sl`` (#526): the optical
        depth is computed on ``self._x_tau``, the model on the pixels.
        Without ``cache``, for a model computed once (see ``PixelLSF``)."""
        r_values = []
        if self._use_variable_resolution and self._resol_column is not None:
            r_values.extend(self._resol_column[sl][np.isfinite(self._resol_column[sl])])
        if self._global_resol_val > 0:
            r_values.append(self._global_resol_val)
        r_values.extend(c.resol for c in self._system_list.components
                        if c.resol is not None and c.resol > 0)
        r_min = float(np.min(r_values)) if r_values else None
        r_max = float(np.max(r_values)) if r_values else None
        # The narrowest line in the pixels sets the fine step, at the
        # current b: a fit that narrows a line below it loses some accuracy
        # (~4e-4 at half the b), recovered by the final model, which is
        # computed again with the fitted values
        xmin, xmax = self._xmin_ang[sl], self._xmax_ang[sl]
        lo, hi = (xmin[0], xmax[-1]) if len(xmin) else (np.inf, -np.inf)
        b_values = []
        for comp in self._system_list.components:
            for t in _transitions(comp.series):
                if lo <= ATOM_DATA[t]['wave'] * (1 + comp.z) <= hi:
                    b_values.append(np.hypot(comp.b, comp.btur or 0.0))
                    break
        b_min = float(np.min(b_values)) if b_values else None
        self._lsf = PixelLSF(self._xmin_ang[sl], self._xmax_ang[sl], r_min=r_min, r_max=r_max,
                             b_min=b_min, cache=cache)
        self._x_tau = self._lsf.x_fine

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
                    self._x_tau, 
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
                    tau_total_dict[res_key] = np.zeros_like(self._x_tau)
                
                # Compute parameters
                idx = i * num_params_per_comp
                z = p_full[idx]
                logN = p_full[idx + 1]
                b_val = p_full[idx + 2]
                btur = p_full[idx + 3]
                
                # Add to group Tau
                tau_total_dict[res_key] += self._compute_tau_component(comp, z, logN, b_val, btur)

        # 3. Convolve and Combine: each pixel gets the fine flux convolved
        # with its own LSF and averaged over the pixel (#526)
        final_model = np.ones_like(self._x_calc)
        
        for (res_val, res_unit), tau_arr in tau_total_dict.items():
            flux_fine = np.exp(-tau_arr)
            
            if res_unit == 'column' and self._resol_calc is not None:
                r, key = self._resol_calc, 'column'
            elif res_val > 0:
                r, key = res_val, res_val
            else:
                r, key = None, 'none'
            final_model *= self._lsf.apply(flux_fine, r, key=key)
            
        return final_model

    @staticmethod
    def _scaled_errors(residuals: np.ndarray, jac: np.ndarray) -> Tuple[np.ndarray, float]:
        """Parameter errors sqrt(diag((J^T J)^-1)) scaled by sqrt(reduced chi2).

        Scaled as in V1 (lmfit ``scale_covar=True``; findings register F29,
        decided with Guido): right also when the error spectrum is off by a
        constant factor. The reduced chi2 is computed only on the pixels that
        constrain the fitted parameters (rows of the Jacobian that are not
        negligible), so that unrelated absorption or misfit elsewhere in the
        fit window does not inflate the errors (issue #534). It is returned
        and stored in each fitted component (``chi2``).
        """
        n_par = jac.shape[1]
        row_norm = np.sqrt(np.sum(jac ** 2, axis=1))
        used = row_norm > 1e-3 * np.max(row_norm) if row_norm.size else row_norm.astype(bool)
        dof = max(1, int(np.sum(used)) - n_par)
        red_chi2 = float(np.sum(residuals[used] ** 2) / dof)
        cov = pinv(jac.T @ jac)
        # A negative diagonal element (ill-conditioned fit) gives a NaN
        # error: reported by fit() instead of a bare RuntimeWarning (F20)
        with np.errstate(invalid='ignore'):
            errors = np.sqrt(np.diag(cov)) * np.sqrt(red_chi2)
        return errors, red_chi2

    def evaluate_uncertainties(self, z_window_kms: float = 20.0) -> SystemListV2:
        """
        Recompute the uncertainties of the active components at their current
        values, without fitting.

        Used at the end of the automatic pipelines: ``refit_all`` fits each
        group while its neighbours still have their initial values, so the
        chi2 of that fit is not the chi2 of the final model (issue #534).
        One numerical Jacobian is computed at the current parameters.

        Returns
        -------
        SystemListV2
            A new system list with ``dz``, ``dlogN``, ``db``, ``dbtur`` (of the
            free parameters) and ``chi2`` of the active components updated.
        """
        lower, upper = self.prepare_fit_context(z_window_kms=z_window_kms)
        p0 = np.clip(self._constraints.p_free_vector, lower, upper)
        r0 = self._residual_function(p0)
        jac = np.empty((len(r0), len(p0)))
        for k in range(len(p0)):
            h = 1e-6 * max(1.0, abs(p0[k]))
            step = h if p0[k] + h <= upper[k] else -h
            p = p0.copy(); p[k] += step
            jac[:, k] = (self._residual_function(p) - r0) / step
        self._active_fit_components = None
        self._active_lines_data = None
        try:
            errors, red_chi2 = self._scaled_errors(r0, jac)
        except Exception as e:
            logging.warning(f"Error calc failed: {e}")
            return self._system_list
        is_free_mask = self._constraints._param_map['is_free']
        full_errors = np.full(len(is_free_mask), np.nan)
        full_errors[is_free_mask] = errors
        active = self._constraints._active_uuids
        comps = []
        for i, c in enumerate(self._system_list.components):
            if active is not None and c.uuid not in active:
                comps.append(c); continue
            idx = i * 4
            new = {'chi2': red_chi2}
            for name, j in (('dz', idx), ('dlogN', idx + 1), ('db', idx + 2),
                            ('dbtur', idx + 3)):
                if is_free_mask[j]:
                    new[name] = full_errors[j]
            comps.append(dataclasses.replace(c, **new))
        return SystemListV2(dataclasses.replace(self._system_list._data, components=comps))

    def _cap_upper(self, bounds: Tuple[np.ndarray, np.ndarray], offset: int,
                   cap: float) -> Tuple[np.ndarray, np.ndarray]:
        """Lower the upper bound of a free parameter of every component
        (``offset`` in z, logN, b, btur) to ``cap``."""
        lower, upper = bounds[0].copy(), bounds[1].copy()
        is_free_mask = self._constraints._param_map['is_free']
        free_idx_map = np.cumsum(is_free_mask) - 1
        for idx in range(offset, len(is_free_mask), 4):
            if is_free_mask[idx]:
                j = free_idx_map[idx]
                upper[j] = max(min(upper[j], cap), lower[j] + 1e-6)
        return lower, upper

    def _cap_b(self, bounds: Tuple[np.ndarray, np.ndarray],
               b_max: float) -> Tuple[np.ndarray, np.ndarray]:
        """Lower the upper bound of the free b parameters to ``b_max``."""
        return self._cap_upper(bounds, 2, b_max)

    def _residual_function(self, p_free: np.ndarray) -> np.ndarray:
        # Check for stop request
        if self._check_stop and self._check_stop():
            raise RuntimeError("Optimization stopped by user.")
        
        model = self._compute_model(p_free)
        resid = (self._y_norm_calc - model) / self._dy_norm_calc
        
        # Handle NaNs/Infs
        resid[~np.isfinite(resid)] = 0.0
        
        final_resid = resid[self._fit_mask_calc]
        
        return final_resid
    
    def _find_feature_limits(self, z_center: float, trans_name: str, 
                             z_window_kms: float = 150.0, logN: float = 14.0) -> Tuple[Optional[float], Optional[float]]:
        atom = ATOM_DATA.get(trans_name)
        if not atom: return None, None
        
        lambda_obs = atom['wave'] * (1 + z_center)
        c_kms = 299792.458
        
        # --- DLA LOGIC ---
        is_dla = logN > 19.0
        
        # 1. Setup Windows
        # If DLA, force a massive search window (at least 3000 km/s) to catch wings
        if is_dla:
            eff_window_kms = max(z_window_kms, 4000.0)
            continuum_threshold = 0.98 # Must fully recover
        else:
            eff_window_kms = z_window_kms
            continuum_threshold = 0.95

        # 1. Find Center Index
        idx_center = np.searchsorted(self._x_ang, lambda_obs)
        idx_center = np.clip(idx_center, 0, len(self._x_ang)-1)

        # 2. Determine Search Range (Indices)
        if len(self._x_ang) > 10:
             sample_dx = np.mean(np.diff(self._x_ang[max(0, idx_center-10):min(len(self._x_ang), idx_center+10)]))
        else: sample_dx = 1.0
        
        dx_search = (eff_window_kms / c_kms) * lambda_obs
        px_width = int(dx_search / sample_dx) if sample_dx > 0 else 50
        
        # Absolute safety bounds
        max_px_width = int((15000.0 / c_kms) * lambda_obs / sample_dx)
        
        # 3. Prepare Smoothed Flux
        safe_start = max(0, idx_center - max_px_width)
        safe_end = min(len(self._x_ang), idx_center + max_px_width)
        if safe_end - safe_start < 5: return None, None
        
        y_view = self._y_norm[safe_start:safe_end]
        dy_view = self._dy_norm[safe_start:safe_end]
        y_smooth = gaussian_filter1d(y_view, sigma=2.0)
        center_rel = idx_center - safe_start
        center_rel = max(0, min(len(y_smooth)-1, center_rel))

        # --- Helper: Scan Logic ---
        def scan_direction(direction):
            lowest_point = y_smooth[center_rel]
            step = direction
            curr = center_rel
            
            while 0 <= curr < len(y_smooth):
                val = y_smooth[curr]
                err = dy_view[curr]
                
                if val < lowest_point: lowest_point = val
                
                # Check Local Max (Separation between blended lines)
                is_local_max = False
                if 0 < curr < len(y_smooth)-1:
                    if y_smooth[curr] >= y_smooth[curr-1] and y_smooth[curr] >= y_smooth[curr+1]:
                        is_local_max = True
                
                # CRITERIA 1: Continuum Reached
                if val > continuum_threshold: 
                    return curr
                
                # CRITERIA 2: Separation (Saddle Point)
                # [FIX] If DLA, IGNORE local maxima (blends). Bulldoze through to find the wing.
                # If Normal line, stop at peaks to separate blends.
                if is_local_max and not is_dla:
                    rise = val - lowest_point
                    is_significant = rise > 4.0 * err
                    is_above_mud = val > 0.1 
                    
                    if is_significant and is_above_mud:
                        return curr

                # CRITERIA 3: Distance Limit
                dist_pix = abs(curr - center_rel)
                
                # If DLA, we respect the massive px_width set above (4000 km/s).
                # If Normal, we respect the user's window, but allow expansion if still deep.
                if dist_pix > px_width:
                    if is_dla:
                        return curr # Hard stop at 4000 km/s
                    else:
                        if val > 0.7: return curr # Normal expansion stop

                curr += step
            return curr

        # 4. Execute Scan
        left_rel = scan_direction(-1)
        right_rel = scan_direction(1)
        
        idx_min = max(0, safe_start + left_rel - 1)
        idx_max = min(len(self._x_ang) - 1, safe_start + right_rel + 1)
        
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
        Computes profiles for M lines on N samples simultaneously.
        
        Input arrays (lambda_0, z, etc.) must be shape (M,).
        wave_grid_ang must be shape (N,).
        Returns summed tau array of shape (N,).
        """
        # Each line only where its optical depth exceeds TAU_CUT (#526: the
        # grid is finer than the pixels): the samples of all lines are
        # gathered in one flat array, and summed back onto the grid
        c_kms = 2.99792458e5
        c_ang_s = 2.99792458e18
        prefactor = 1.4974e-15
        b_safe = np.maximum(b_kms, 0.1)
        lambda_c = lambda_0 * (1.0 + z)
        dop_width_ang = np.maximum((b_safe / c_kms) * lambda_c, 1e-10)
        # Damping parameter in the rest frame, a = Gamma lambda_0 / (4 pi b): it
        # was computed with the observed wavelength, (1 + z) too large (F61)
        a = (gamma * lambda_0) / (4.0 * np.pi * c_ang_s * (b_safe / c_kms))
        tau0 = prefactor * N_col * f_val * lambda_0 / b_safe

        # Gaussian core and Lorentzian wings, in Doppler widths
        with np.errstate(divide='ignore', invalid='ignore'):
            x_cut = np.maximum(np.sqrt(np.log(np.maximum(tau0 / TAU_CUT, np.e))),
                               np.sqrt(tau0 * a / (np.sqrt(np.pi) * TAU_CUT))) + 3.0
        x_cut = np.nan_to_num(x_cut, nan=0.0, posinf=0.0)
        lo = np.searchsorted(wave_grid_ang, lambda_c - x_cut * dop_width_ang, side='left')
        hi = np.searchsorted(wave_grid_ang, lambda_c + x_cut * dop_width_ang, side='right')
        lens = np.maximum(hi - lo, 0)
        total = int(lens.sum())
        if total == 0:
            return np.zeros_like(wave_grid_ang)
        line = np.repeat(np.arange(len(lens)), lens)
        idx = lo[line] + np.arange(total) - np.repeat(np.cumsum(lens) - lens, lens)

        # Faddeeva (Voigt): the expensive step
        x = (wave_grid_ang[idx] - lambda_c[line]) / dop_width_ang[line]
        H_ax = wofz(x + 1j * a[line]).real
        return np.bincount(idx, weights=tau0[line] * H_ax, minlength=len(wave_grid_ang))

    def prepare_fit_context(self, z_window_kms: float = 20.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Prepares the internal state for fitting (masks, background, vectorized data).
        
        Returns
        -------
        tuple
            (lower_bounds, upper_bounds) for the free parameters.
        """
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
                    lam_min, lam_max = self._find_feature_limits(
                        comp.z, trans, 
                        z_window_kms=z_window_kms, 
                        logN=comp.logN
                    )
                    
                    use_fallback = True
                    if lam_min is not None and lam_max is not None:
                        width_v = c_kms * (lam_max - lam_min) / lam_max
                        if width_v > 30.0:
                            use_fallback = False
                            mask_feat = (self._x_ang >= lam_min) & (self._x_ang <= lam_max)
                            self._fit_mask |= mask_feat
    
                    if use_fallback:
                        lam_0 = ATOM_DATA.get(trans, {}).get('wave')
                        if lam_0:
                            lam_obs = lam_0 * (1.0 + comp.z)
                            safe_window = max(z_window_kms, 150.0)
                            dw = (safe_window / c_kms) * lam_obs
                            self._fit_mask |= (self._x_ang > lam_obs - dw) & (self._x_ang < lam_obs + dw)
            self._fit_mask &= self._valid_mask
        
        if np.any(self._fit_mask):
            indices = np.where(self._fit_mask)[0]
            self._sl = slice(max(0, np.min(indices) - 50), min(len(self._x_ang), np.max(indices) + 51))
        else:
            logging.warning("Fit mask is empty! Falling back to full valid mask.")
            self._sl = slice(None)
            self._fit_mask = self._valid_mask.copy()

        self._x_calc = self._x_ang[self._sl]
        self._y_norm_calc = self._y_norm[self._sl]
        self._dy_norm_calc = self._dy_norm[self._sl]
        self._fit_mask_calc = self._fit_mask[self._sl]
        
        #logging.info(f"Fitting context: {len(self._x_calc)} pixels in slice, {np.sum(self._fit_mask_calc)} pixels in fit mask.")
        if np.sum(self._fit_mask_calc) == 0:
            logging.debug("FIT MASK IS EMPTY!")
            # Debug: what are we searching for?
            for idx, comp in enumerate(active_components):
                logging.debug(f"  Comp {idx}: z={comp.z:.5f}, series={comp.series}")
        
        if self._use_variable_resolution and self._resol_column is not None:
            self._resol_calc = self._resol_column[self._sl]
        else:
            self._resol_calc = None
        self._set_lsf(self._sl)

        # Background
        self._cached_tau_groups = {}
        if len(self._x_calc) > 0:
            x_min_win, x_max_win = self._x_calc[0], self._x_calc[-1]
        else:
            x_min_win, x_max_win = -np.inf, np.inf
        buffer_ang = 100.0
        
        relevant_static_comps = []
        if active_uuids is not None:
            for comp in self._system_list.components:
                if comp.uuid in active_uuids: continue
                is_relevant = False
                for t in get_trans_list(comp.series):
                    atom = ATOM_DATA.get(t)
                    if not atom: continue
                    obs_w = atom['wave'] * (1 + comp.z)
                    if (x_min_win - buffer_ang) < obs_w < (x_max_win + buffer_ang):
                        is_relevant = True; break
                if is_relevant: relevant_static_comps.append(comp)
        
        for comp in relevant_static_comps:
            res_key = self._get_component_resolution(comp)
            if res_key not in self._cached_tau_groups:
                self._cached_tau_groups[res_key] = np.zeros_like(self._x_tau)
            self._cached_tau_groups[res_key] += self._compute_tau_component(comp, comp.z, comp.logN, comp.b, comp.btur)

        # Active Data
        self._active_fit_components = []
        if active_uuids is not None:
            for i, comp in enumerate(self._system_list.components):
                if comp.uuid in active_uuids: self._active_fit_components.append((i, comp))
        else:
            self._active_fit_components = None

        self._active_lines_data = {}
        iterator = self._active_fit_components if self._active_fit_components is not None else enumerate(self._system_list.components)
        grouped_lines = {}
        
        for idx, comp in iterator:
            res_key = self._get_component_resolution(comp)
            if res_key not in grouped_lines:
                grouped_lines[res_key] = {'lambda': [], 'f': [], 'gamma': [], 'idx_z': [], 'idx_N': [], 'idx_b': [], 'idx_btur': []}
            t_list = get_trans_list(comp.series)
            base_idx = idx * 4
            for t in t_list:
                atom = ATOM_DATA.get(t)
                if not atom: continue
                g = grouped_lines[res_key]
                g['lambda'].append(atom['wave']); g['f'].append(atom['f']); g['gamma'].append(atom['gamma'])
                g['idx_z'].append(base_idx); g['idx_N'].append(base_idx + 1); g['idx_b'].append(base_idx + 2); g['idx_btur'].append(base_idx + 3)

        for rk, data in grouped_lines.items():
            if not data['lambda']: continue
            self._active_lines_data[rk] = {k: np.array(v) for k, v in data.items()}

        # Bounds
        lower_bounds, upper_bounds = self._constraints.get_bounds()
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
            
            idx_z = free_idx_map[base_idx]
            
            # Narrow the physical bounds based only on the requested z_window_kms
            z_final_min = max(lower_bounds[idx_z], z_min_win)
            z_final_max = min(upper_bounds[idx_z], z_max_win)

            cur = self._constraints.p_free_vector[idx_z]
            lower_bounds[idx_z] = min(z_final_min, cur - 1e-6)
            upper_bounds[idx_z] = max(z_final_max, cur + 1e-6)
        
        return lower_bounds, upper_bounds


    def get_fit_context_data(self) -> Dict[str, Any]:
        """
        Retrieve all internal data required for external likelihood calculations.

        Extracts the sliced arrays, fit masks, and the model computation function
        used by the frequentist fitter, allowing external engines (like Bayesian
        samplers) to evaluate the exact same physical model. Must be called
        after :meth:`prepare_fit_context`.

        Returns
        -------
        dict
            A dictionary containing:
            - ``'x'``: The wavelength grid array.
            - ``'y'``: The normalized flux array.
            - ``'dy'``: The normalized flux error array.
            - ``'mask'``: The boolean fit mask array.
            - ``'compute_model'``: A callable function that takes a free parameter vector and returns the normalized model array.
        """
        return {
            'x': self._x_calc,
            'y': self._y_norm_calc,
            'dy': self._dy_norm_calc,
            'mask': self._fit_mask_calc,
            'compute_model': self._compute_model
        }

    def fit(self, max_nfev: int = 2000, method: str = 'trf', 
            z_window_kms: float = 20.0, verbose: int = 0,
            b_max: Optional[float] = None,
            logN_max: Optional[float] = None) -> Tuple[SystemListV2, np.ndarray, Any]:
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
        b_max : float, optional
            Upper bound (km/s) on the free b parameters, tighter than the
            physical limit of the constraint model (300 km/s). Used by the
            automatic pipelines to keep the components from broadening into
            their neighbours; components that end at this bound may be broader
            and can be refitted without it. Defaults to ``None`` (no extra
            bound).
        logN_max : float, optional
            Upper bound on the free logN parameters, tighter than the
            physical limit of the constraint model (23). Used by the
            automatic pipelines (19): in saturated troughs a single
            component can otherwise fill the trough with its damping wings
            and run away to logN 20-23 (findings register F62). Defaults to
            ``None`` (no extra bound).

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
        p0 = self._constraints.p_free_vector
        #logging.info(f"Initial p0: {p0}")
        # 1. Refine initial guess using centroiding (Smart Guess)
        # p0 = self._smart_guess(p0, z_window_kms)
        #print(f"DEBUG: p0 after skipping smart_guess: {p0}")

        # 2. Setup the fitting context (pixels, background, vectors) based on refined guess
        # We need to update the system list temporarily to ensure prepare_fit_context 
        # uses the refined positions for mask calculation.
        p_full_refined = self._constraints.map_p_free_to_full(p0)
        # This is a bit hacky but prevents the mask from being off-center
        from dataclasses import replace
        old_systs = self._system_list
        new_comps = []
        for i, comp in enumerate(old_systs.components):
            idx = i * 4
            new_comps.append(replace(comp, z=p_full_refined[idx], logN=p_full_refined[idx+1], b=p_full_refined[idx+2], btur=p_full_refined[idx+3]))
        self._system_list = SystemListV2(replace(old_systs._data, components=new_comps))

        bounds = self.prepare_fit_context(z_window_kms=z_window_kms)
        if b_max is not None:
            bounds = self._cap_b(bounds, b_max)
        if logN_max is not None:
            bounds = self._cap_upper(bounds, 1, logN_max)

        # Clamp p0 to bounds
        p0 = np.clip(p0, bounds[0], bounds[1])

        # Fit
        res = least_squares(self._residual_function, p0, bounds=bounds, method=method, loss='linear', max_nfev=max_nfev, x_scale='jac', verbose=verbose)
        
        # Restore system list (it will be updated properly at the end)
        self._system_list = old_systs
        
        # --- CLEANUP (Reset Optimized Caches) ---
        self._active_fit_components = None
        self._active_lines_data = None

        # Statistics
        p_free_errors = np.zeros_like(res.x)
        is_free_mask = self._constraints._param_map['is_free']
        red_chi2 = 0
        try:
            p_free_errors, red_chi2 = self._scaled_errors(res.fun, res.jac)
            logging.debug(f"Errors calculated: {p_free_errors}")
        except Exception as e:
            logging.warning(f"Error calc failed: {e}")
            p_free_errors = np.zeros_like(res.x)
            p_free_errors[:] = np.nan

        p_fitted_full = self._constraints.map_p_free_to_full(res.x)
        p_full_errors = np.zeros_like(p_fitted_full)
        p_full_errors[is_free_mask] = p_free_errors
        
        new_components = []
        old_components = self._system_list.components
        active_uuids = None
        if hasattr(self, '_constraints'):
            # The constraint model stores the active set in _active_uuids
            if hasattr(self._constraints, '_active_uuids'):
                active_uuids = self._constraints._active_uuids

        for i, old_comp in enumerate(old_components):
            idx = i * 4

            # Spectator Check
            # A component is a spectator if a specific active set exists 
            # AND this component is NOT in it.
            is_spectator = (active_uuids is not None) and (old_comp.uuid not in active_uuids)
            
            # Helper to assign error
            def resolve_error(param_idx, old_err):
                if is_free_mask[param_idx]:
                    return p_full_errors[param_idx] # New calculated error
                elif is_spectator:
                    return old_err # PRESERVE old error for spectators
                else:
                    return 0.0 # ZERO error for user-frozen params
                
            dz_new = resolve_error(idx, old_comp.dz)
            dlogN_new = resolve_error(idx+1, old_comp.dlogN)
            db_new = resolve_error(idx+2, old_comp.db)
            # Stored like the other errors (before, dbtur kept its value from
            # before the fit: findings register F36). btur is frozen unless
            # the user frees it, so its error is usually 0
            dbtur_new = resolve_error(idx+3, old_comp.dbtur)
            
            # Metadata Update. The resolution of the component is left as it
            # is: it is the user's override, None to follow the spectrum.
            # Before #526 the fitter stamped the resolution it used, so the
            # fitted components ignored later changes of the spectrum
            comp_chi2 = old_comp.chi2
            comp_resol = old_comp.resol
            
            if active_uuids is None or old_comp.uuid in active_uuids:
                comp_chi2 = red_chi2
            
            new_comp = dataclasses.replace(old_comp,
                z=p_fitted_full[idx], 
                logN=p_fitted_full[idx+1], 
                b=p_fitted_full[idx+2], 
                btur=p_fitted_full[idx+3],
                dz=dz_new, 
                dlogN=dlogN_new, 
                db=db_new,
                dbtur=dbtur_new,
                chi2=comp_chi2,
                resol=comp_resol
            )
            new_components.append(new_comp)
            
        # Report the components whose uncertainties could not be computed
        # (findings register F20): the errors stay NaN, which the tables show
        nan_errors = [c for c, old in zip(new_components, old_components)
                      if (active_uuids is None or old.uuid in active_uuids)
                      and any(e is not None and np.isnan(e)
                              for e in (c.dz, c.dlogN, c.db, c.dbtur))]
        if nan_errors:
            listed = ", ".join(f"{c.series} at z = {c.z:.5f}" for c in nan_errors[:10])
            more = f" and {len(nan_errors) - 10} more" if len(nan_errors) > 10 else ""
            logging.warning(
                f"Uncertainties could not be computed (ill-conditioned fit) for "
                f"{len(nan_errors)} component(s): {listed}{more}. "
                f"Their errors are NaN; consider refitting or freezing parameters.")

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
        saved = (self._x_calc, self._resol_calc, self._x_tau, self._lsf)
        self._x_calc = self._x_ang
        if self._use_variable_resolution:
             self._resol_calc = self._resol_column
        else:
             self._resol_calc = None
        self._set_lsf(slice(None), cache=False)
             
        # This will now loop over ALL components because active_uuids is None
        final_model_norm = self._compute_model(res.x)
        
        # Restore state (good practice)
        self._constraints._active_uuids = stored_active_uuids
        self._x_calc, self._resol_calc, self._x_tau, self._lsf = saved
        
        final_model_flux = final_model_norm * self._norm

        logging.debug(f"Fit complete. Reduced chi2 (constraining pixels): {red_chi2:.2f}.")
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
        saved = (self._x_calc, self._resol_calc, self._x_tau, self._lsf)
        
        # Switch to full spectrum context
        self._x_calc = self._x_ang
        if self._use_variable_resolution:
             self._resol_calc = self._resol_column
        else:
             self._resol_calc = None
        self._set_lsf(slice(None), cache=False)

        # Reset caches to force Case B (Standard Loop) in _compute_model
        self._cached_tau_groups = {}
        self._active_lines_data = None 
        
        y_model_norm = self._compute_model(p_current_free)
        
        # Restore context
        self._x_calc, self._resol_calc, self._x_tau, self._lsf = saved
        
        y_model_phys = y_model_norm * self._norm
        return self._x_ang, y_model_phys