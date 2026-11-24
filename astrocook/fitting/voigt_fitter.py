from astrocook._renamed import renamed
import astropy.units as au
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares
from scipy.special import wofz
from scipy.linalg import pinv
from scipy import sparse
import dataclasses
import logging
import numpy as np
from typing import List, Literal, Tuple, Dict, Any, Optional, Union, overload

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS, is_multiplet, series_lines
from astrocook.fitting.lsf import PixelLSF
from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import ComponentData
from astrocook.core.system_list import SystemList

logger = logging.getLogger(__name__)

#: A line is evaluated only where its optical depth exceeds this (#526)
TAU_CUT = 1e-6

#: Flags that a fit leaves on the components it fits (#517), in
#: ComponentData.flags: a free parameter ended at a bound of the fit
#: (within AT_BOUND_RTOL of the width of its range), or the uncertainties
#: could not be computed. The upper bounds of the fit are those of the
#: constraint model (logN 23, b 300 km/s) or the tighter max_logN and max_b
#: of the automatic pipelines (19, 50 km/s); the z bound is the window of
#: the fit. A parameter at a bound has an unreliable value and error.
COMPONENT_FLAGS = ('z_at_window', 'logN_at_lower_bound', 'logN_at_upper_bound',
                   'b_at_lower_bound', 'b_at_upper_bound', 'nan_errors')
AT_BOUND_RTOL = 1e-3


def chi2_of(residuals: np.ndarray) -> float:
    """The chi2 of normalized residuals, the non-finite ones left out."""
    return float(np.sum(residuals[np.isfinite(residuals)] ** 2))


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


def _transitions_of(series: str) -> List[str]:
    """The lines of a series as the fit takes them: those of a multiplet,
    or the series itself if it is one line (a list with commas: none)."""
    if is_multiplet(series):
        return series_lines(series)
    if series in ATOM_DATA:
        return [series]
    return []


def _line_terms(lambda_0: np.ndarray, f_val: np.ndarray, gamma: np.ndarray,
                z: np.ndarray, N_col: np.ndarray, b_kms: np.ndarray
                ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """The observed centre (Angstrom), the Doppler width (Angstrom), the
    damping parameter, the central optical depth of each line, and the
    distance from the centre, in Doppler widths, beyond which its optical
    depth is below ``TAU_CUT``."""
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
    return lambda_c, dop_width_ang, a, tau0, x_cut


def line_reach(lambda_0: np.ndarray, f_val: np.ndarray, gamma: np.ndarray,
               z: np.ndarray, N_col: np.ndarray, b_kms: np.ndarray
               ) -> Tuple[np.ndarray, np.ndarray]:
    """The observed centre of each line and the half-width (both in
    Angstrom) beyond which :func:`voigt_optical_depth` leaves it out (its
    optical depth below ``TAU_CUT``): the reach of a line on the pixels,
    before the convolution."""
    lambda_c, dop_width_ang, _, _, x_cut = _line_terms(lambda_0, f_val, gamma, z, N_col, b_kms)
    return lambda_c, x_cut * dop_width_ang


def _line_samples(wave_grid_ang: np.ndarray, lambda_0: np.ndarray, f_val: np.ndarray,
                  gamma: np.ndarray, z: np.ndarray, N_col: np.ndarray, b_kms: np.ndarray):
    """The samples of the grid where each line is evaluated (its optical
    depth above ``TAU_CUT``, #526): the index of each sample in the grid, the
    line it belongs to, its distance from the centre in Doppler widths, the
    Faddeeva function there, and the terms of the lines
    (:func:`_line_terms`). ``None`` if no sample."""
    lambda_c, dop_width_ang, a, tau0, x_cut = _line_terms(
        lambda_0, f_val, gamma, z, N_col, b_kms)
    lo = np.searchsorted(wave_grid_ang, lambda_c - x_cut * dop_width_ang, side='left')
    hi = np.searchsorted(wave_grid_ang, lambda_c + x_cut * dop_width_ang, side='right')
    lens = np.maximum(hi - lo, 0)
    total = int(lens.sum())
    if total == 0:
        return None
    line = np.repeat(np.arange(len(lens)), lens)
    idx = lo[line] + np.arange(total) - np.repeat(np.cumsum(lens) - lens, lens)

    # Faddeeva (Voigt): the expensive step
    x = (wave_grid_ang[idx] - lambda_c[line]) / dop_width_ang[line]
    w = wofz(x + 1j * a[line])
    return idx, line, x, w, (lambda_c, dop_width_ang, a, tau0)


def voigt_optical_depth(wave_grid_ang: np.ndarray,
                        lambda_0: np.ndarray, f_val: np.ndarray,
                        gamma: np.ndarray, z: np.ndarray,
                        N_col: np.ndarray, b_kms: np.ndarray) -> np.ndarray:
    """
    Optical depth of M Voigt lines on N samples, summed: the physics of
    the fitter, used also by the profiles of the System Inspector.

    Input arrays (lambda_0, z, etc.) must be shape (M,).
    wave_grid_ang must be shape (N,).
    Returns summed tau array of shape (N,).
    """
    # Each line only where its optical depth exceeds TAU_CUT (#526: the
    # grid is finer than the pixels): the samples of all lines are
    # gathered in one flat array, and summed back onto the grid
    samples = _line_samples(wave_grid_ang, lambda_0, f_val, gamma, z, N_col, b_kms)
    if samples is None:
        return np.zeros_like(wave_grid_ang)
    idx, line, _, w, (_, _, _, tau0) = samples
    return np.bincount(idx, weights=tau0[line] * w.real, minlength=len(wave_grid_ang))


def voigt_optical_depth_derivatives(wave_grid_ang: np.ndarray,
                                    lambda_0: np.ndarray, f_val: np.ndarray,
                                    gamma: np.ndarray, z: np.ndarray,
                                    N_col: np.ndarray, b_kms: np.ndarray):
    """The optical depth of :func:`voigt_optical_depth` and its analytic
    derivatives (phase 7, PR 7c.4d-bis, F171), sample by sample: the index
    of each sample in the grid, the line it belongs to, and the derivatives
    of the optical depth of that line there with respect to its z, its
    log10 N and its b (the b given, ``b_kms``). ``tau`` is summed on the
    grid, as :func:`voigt_optical_depth` (the same operations: identical).

    With tau = tau0 H(x, a), H the real part of the Faddeeva function
    w(x + i a), x = (lambda - lambda_c) / Delta_D, Delta_D = b lambda_c / c,
    a = Gamma lambda_0 / (4 pi b), tau0 proportional to N / b:
    w' = -2 (x + i a) w + 2 i / sqrt(pi), dH/dx = Re w', dH/da = -Im w';
    d tau / d log10 N = ln(10) tau; d tau / dz = tau0 dH/dx dx/dz with
    dx/dz = -(x + c / b) / (1 + z); d tau / db = -(tau0 / b)
    (H + x dH/dx + a dH/da). The cut of each line at ``TAU_CUT`` moves with
    its parameters: its derivative (the samples entering or leaving the
    cut, below 1e-6 in optical depth) is left out.

    Returns
    -------
    tuple or None
        ``(tau, idx, line, dtau_dz, dtau_dlogN, dtau_db)``, ``None`` if no
        line has samples on the grid.
    """
    c_kms = 2.99792458e5
    samples = _line_samples(wave_grid_ang, lambda_0, f_val, gamma, z, N_col, b_kms)
    if samples is None:
        return None
    idx, line, x, w, (_, _, a, tau0) = samples
    a_s = a[line]
    t0 = tau0[line]
    tau_s = t0 * w.real
    tau = np.bincount(idx, weights=tau_s, minlength=len(wave_grid_ang))
    # w' = -2 (x + i a) w + 2 i / sqrt(pi)
    dh_dx = -2.0 * (x * w.real - a_s * w.imag)
    dh_da = 2.0 * (x * w.imag + a_s * w.real) - 2.0 / np.sqrt(np.pi)
    b_safe = np.maximum(b_kms, 0.1)[line]
    dtau_dz = t0 * dh_dx * (-(x + c_kms / b_safe) / (1.0 + z[line]))
    dtau_dlogN = np.log(10.0) * tau_s
    dtau_db = -(t0 / b_safe) * (w.real + x * dh_dx + a_s * dh_da)
    # Below the smallest b of the model (0.1 km/s) the optical depth does
    # not follow b
    dtau_db = np.where(b_kms[line] > 0.1, dtau_db, 0.0)
    return tau, idx, line, dtau_dz, dtau_dlogN, dtau_db


class NoPixelsToFit(ValueError):
    """No valid pixel is around the components to fit: their lines fall
    outside the spectrum, or on pixels that cannot be used (F111, phase
    6-ter, PR 7-quater). A ``ValueError``: the recipes report it as a fit
    that could not be done."""


class VoigtFitter:
    """
    Engine for fitting Voigt profiles to spectral data.

    This class orchestrates the optimization process. It handles:

    1.  **Resolution Logic**: determining effective resolution (R or FWHM) from columns or metadata.
    2.  **Dynamic Masking**: identifying relevant pixels for the fit to speed up computation.
    3.  **Optimization**: utilizing ``scipy.optimize.least_squares`` to minimize residuals.
    4.  **Vectorization**: calculating optical depths for multiple lines simultaneously.

    Parameters
    ----------
    spectrum : Spectrum
        The spectral data to fit. Must have a valid normalization (continuum) or be pre-normalized.
    system_list : SystemList
        The list of absorption components to model.
    """
    def __init__(self, spectrum: Spectrum, system_list: SystemList, check_stop=None):
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
        self._x_tau: Any = None
        self._lsf: Any = None
        
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
        self._resol_column: Any = None
        self._use_variable_resolution = False
        self._global_resol_val = 0.0
        self._global_resol_unit = 'R'
        r_pix = self._spectrum.resolution()
        finite = r_pix[np.isfinite(r_pix)]
        if finite.size == r_pix.size and finite.size and np.all(finite == finite[0]):
            self._global_resol_val = float(finite[0])
            logger.debug(f"VoigtFitter: Global R={self._global_resol_val:.0f}")
        elif finite.size:
            self._resol_column = r_pix
            self._use_variable_resolution = True
            logger.debug(f"VoigtFitter: Variable resolution, median R={np.median(finite):.0f}")
        else:
            logger.debug("VoigtFitter: No resolution found. Fitting unconvolved profiles.")

        # Placeholders
        self._fit_mask: Any = None
        self._sl = slice(None)
        self._x_calc: Any = None
        self._y_norm_calc: Any = None
        self._dy_norm_calc: Any = None
        self._fit_mask_calc: Any = None
        self._resol_calc: Any = None
        
        # Cache is now a dict: { (resol_val, unit_str): tau_array }
        self._cached_tau_groups: Dict[Any, Any] = {}
        
        # [OPTIMIZATION] Cache list of active components for the inner loop
        # Stores tuples: (index_in_system_list, component_object)
        self._active_fit_components: Optional[List[Tuple[int, ComponentData]]] = None
        self._active_lines_data: Any = None
        self._link_input_cache: Optional[List[int]] = None
        #: The model of the whole spectrum only on some pixels (F167): their
        #: indices and the fine samples (start, stop) they take
        self._model_pixels: Optional[Tuple[np.ndarray, int, int]] = None


    def _compute_tau_component(self, comp: ComponentData, z: float, logN: float, b: float, btur: float) -> np.ndarray:
        atoms = [ATOM_DATA[t] for t in _transitions(comp.series)]
        if not atoms:
            return np.zeros_like(self._x_tau)
        n = len(atoms)
        return self._voigt_optical_depth_vectorized(
            self._x_tau,
            np.array([at['wave'] for at in atoms]), np.array([at['f'] for at in atoms]),
            np.array([at['gamma'] for at in atoms]), np.full(n, float(z)),
            np.full(n, 10**logN), np.full(n, np.sqrt(b**2 + btur**2)))

    def _get_component_resolution(self, comp: ComponentData) -> Tuple[float, str]:
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
        r_values: List[Any] = []
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
        
        for res_key, tau_arr in tau_total_dict.items():
            final_model *= self._convolved(res_key, tau_arr)
            
        return final_model

    def _free_to_full(self, p_free: np.ndarray) -> sparse.csr_matrix:
        """The derivatives of the full parameter vector with respect to the
        free one (``map_p_free_to_full``), sparse (full x free): 1 for a free
        parameter, the gradient of its expression for a linked one (central
        differences on the free parameters the links read, as
        ``_linked_errors``: exact for the usual links, a copy or a shift),
        nothing for a frozen one."""
        is_free = self._constraints._param_map['is_free']
        n_full, n_free = len(is_free), len(p_free)
        rows = list(np.flatnonzero(is_free))
        cols = list(range(n_free))
        vals = [1.0] * n_free
        links = self._constraints._compiled_links
        if links and n_free:
            linked = np.array([idx for idx, _ in links])
            for k in self._link_inputs(p_free):
                h = 1e-6 * max(1.0, abs(p_free[k]))
                p_plus = p_free.copy(); p_plus[k] += h
                p_minus = p_free.copy(); p_minus[k] -= h
                grad = (self._constraints.map_p_free_to_full(p_plus)[linked]
                        - self._constraints.map_p_free_to_full(p_minus)[linked]) / (2 * h)
                nz = np.flatnonzero(grad)
                rows.extend(linked[nz]); cols.extend([k] * len(nz)); vals.extend(grad[nz])
        return sparse.csr_matrix((vals, (rows, cols)), shape=(n_full, n_free))

    def _link_inputs(self, p_free: np.ndarray) -> List[int]:
        """The free parameters that the links may read: those whose change
        moves a linked parameter (found once per fit, by moving each free
        parameter in turn; all of them if any link is not a function of the
        free parameters alone)."""
        if self._link_input_cache is not None:
            return self._link_input_cache
        links = self._constraints._compiled_links
        linked = np.array([idx for idx, _ in links])
        base = self._constraints.map_p_free_to_full(p_free)[linked]
        inputs = []
        for k in range(len(p_free)):
            p = p_free.copy(); p[k] += 1e-3 * max(1.0, abs(p_free[k]))
            if np.any(self._constraints.map_p_free_to_full(p)[linked] != base):
                inputs.append(k)
        self._link_input_cache = inputs
        return inputs

    def _jacobian(self, p_free: np.ndarray) -> np.ndarray:
        """The analytic Jacobian of ``_residual_function`` (phase 7, PR
        7c.4d-bis, F171), in place of finite differences (one model per free
        parameter): the derivatives of the optical depth of each line
        (:func:`voigt_optical_depth_derivatives`), times the fine flux,
        convolved with the LSF as the flux (``PixelLSF.matrix``), the
        resolution groups by the product rule, through the links
        (``_free_to_full``). Equal to the finite differences within their
        own error (tests)."""
        if self._check_stop and self._check_stop():
            raise RuntimeError("Optimization stopped by user.")
        p_full = self._constraints.map_p_free_to_full(p_free)
        n_full = len(p_full)
        d_free = self._free_to_full(p_free)
        groups = []
        for res_key in set(self._cached_tau_groups) | set(self._active_lines_data or {}):
            tau = self._cached_tau_groups.get(res_key)
            tau = np.zeros_like(self._x_tau) if tau is None else tau.copy()
            data = (self._active_lines_data or {}).get(res_key)
            d_fine = None
            if data is not None:
                btur = p_full[data['idx_btur']]
                b = p_full[data['idx_b']]
                b_eff = np.sqrt(b**2 + btur**2)
                out = voigt_optical_depth_derivatives(
                    self._x_tau, data['lambda'], data['f'], data['gamma'],
                    p_full[data['idx_z']], 10**p_full[data['idx_N']], b_eff)
                if out is not None:
                    tau_act, idx, line, d_z, d_logn, d_beff = out
                    tau += tau_act
                    with np.errstate(invalid='ignore', divide='ignore'):
                        db_db = np.where(b_eff > 0, b / b_eff, 1.0)[line]
                        db_dbtur = np.where(b_eff > 0, btur / b_eff, 0.0)[line]
                    cols = np.concatenate([data['idx_z'][line], data['idx_N'][line],
                                           data['idx_b'][line], data['idx_btur'][line]])
                    rows = np.tile(idx, 4)
                    vals = np.concatenate([d_z, d_logn, d_beff * db_db, d_beff * db_dbtur])
                    d_fine = (rows, cols, vals)
            d_pix = None
            if d_fine is not None:
                rows, cols, vals = d_fine
                g = sparse.csr_matrix((-np.exp(-tau[rows]) * vals, (rows, cols)),
                                      shape=(len(self._x_tau), n_full))
                d_pix = self._lsf_matrix(res_key) @ (g @ d_free)
            groups.append((res_key, tau, d_pix))
        # The model is the product of the convolved flux of the groups
        jac = np.zeros((len(self._x_calc), len(p_free)))
        for i, (_, _, d_pix) in enumerate(groups):
            if d_pix is None:
                continue
            others = np.ones(len(self._x_calc))
            for j, (res_key, tau, _) in enumerate(groups):
                if j != i:
                    others *= self._convolved(res_key, tau)
            jac += d_pix.multiply(others[:, None]).toarray()
        # The residuals are (data - model) / error, 0 where not finite
        mask = self._fit_mask_calc
        resid_ok = np.isfinite(self._y_norm_calc) & np.isfinite(self._dy_norm_calc)
        jac = -jac / self._dy_norm_calc[:, None]
        jac[~resid_ok] = 0.0
        return jac[mask]

    def _lsf_matrix(self, res_key: Tuple[float, str]) -> sparse.csr_matrix:
        """The LSF of ``_convolved`` as a matrix (``PixelLSF.matrix``)."""
        res_val, res_unit = res_key
        if res_unit == 'column' and self._resol_calc is not None:
            return self._lsf.matrix(self._resol_calc, key='column')
        if res_val > 0:
            return self._lsf.matrix(res_val, key=res_val)
        return self._lsf.matrix(None, key='none')

    @staticmethod
    def _scaled_errors(residuals: np.ndarray, jac: np.ndarray
                       ) -> Tuple[np.ndarray, float, np.ndarray]:
        """Parameter errors sqrt(diag((J^T J)^-1)) scaled by sqrt(reduced chi2).

        Scaled as in V1 (lmfit ``scale_covar=True``; findings register F29,
        decided with Guido): right also when the error spectrum is off by a
        constant factor. The reduced chi2 is computed only on the pixels that
        constrain the fitted parameters (rows of the Jacobian that are not
        negligible), so that unrelated absorption or misfit elsewhere in the
        fit window does not inflate the errors (issue #534). It is returned
        and stored in each fitted component (``chi2``), with the covariance
        matrix scaled in the same way (for the errors of the linked
        parameters, ``_linked_errors``).
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
        return errors, red_chi2, cov * red_chi2

    def _linked_errors(self, p_free: np.ndarray, cov: np.ndarray) -> Dict[int, float]:
        """Errors of the linked parameters, propagated from the free ones.

        A linked parameter is a function of the free parameters through the
        expression of its link (``map_p_free_to_full``): its error is
        sqrt(g C g^T), with C the covariance of the free parameters and g the
        gradient of the expression, computed numerically (central
        differences; exact for the usual links, a copy or a shift of the
        target). Before, a linked parameter was given error 0, which the
        ASCII list of the systems exported (issue #630). A parameter linked
        to targets that are not free in this fit (frozen, or outside the
        fitted group) gets 0, as a frozen parameter.

        Returns
        -------
        dict
            Index in the full parameter vector -> error.
        """
        linked = [idx for idx, _ in self._constraints._compiled_links]
        if not linked or len(p_free) == 0:
            return {}
        grad = np.empty((len(linked), len(p_free)))
        for k in range(len(p_free)):
            h = 1e-6 * max(1.0, abs(p_free[k]))
            p_plus = p_free.copy(); p_plus[k] += h
            p_minus = p_free.copy(); p_minus[k] -= h
            grad[:, k] = (self._constraints.map_p_free_to_full(p_plus)[linked]
                          - self._constraints.map_p_free_to_full(p_minus)[linked]) / (2 * h)
        with np.errstate(invalid='ignore'):
            errors = np.sqrt(np.einsum('ij,jk,ik->i', grad, cov, grad))
        return dict(zip(linked, errors, strict=True))

    def evaluate_uncertainties(self, z_window_kms: float = 20.0) -> SystemList:
        """
        Recompute the uncertainties of the active components at their current
        values, without fitting.

        Used at the end of the automatic pipelines: ``refit_all`` fits each
        group while its neighbours still have their initial values, so the
        chi2 of that fit is not the chi2 of the final model (issue #534).
        One Jacobian (analytic) is computed at the current parameters.

        Returns
        -------
        SystemList
            A new system list with ``dz``, ``dlogN``, ``db``, ``dbtur`` (of the
            free parameters) and ``chi2`` of the active components updated.
        """
        lower, upper = self.prepare_fit_context(z_window_kms=z_window_kms)
        p0 = np.clip(self._constraints.p_free_vector, lower, upper)
        r0 = self._residual_function(p0)
        if self._active_lines_data:
            # The analytic Jacobian, as the fit (phase 7, PR 7c.4d-bis):
            # before, forward differences with a step of 1e-6 of the
            # parameter (in z about 0.3 km/s at z = 3), which moved the
            # errors by up to about 1% from those of the fit
            jac = self._jacobian(p0)
        else:
            jac = np.empty((len(r0), len(p0)))
            for k in range(len(p0)):
                h = 1e-6 * max(1.0, abs(p0[k]))
                step = h if p0[k] + h <= upper[k] else -h
                p = p0.copy(); p[k] += step
                jac[:, k] = (self._residual_function(p) - r0) / step
        self._active_fit_components = None
        self._active_lines_data = None
        try:
            errors, red_chi2, cov = self._scaled_errors(r0, jac)
            linked_errors = self._linked_errors(p0, cov)
        except (ValueError, ArithmeticError, KeyError, TypeError) as e:
            # ValueError covers np.linalg.LinAlgError; KeyError and
            # TypeError come from the expression of a link
            logger.warning(f"The uncertainties could not be computed ({e}): the "
                           "components keep their previous uncertainties.")
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
            new: Dict[str, Any] = {'chi2': red_chi2}
            for name, j in (('dz', idx), ('dlogN', idx + 1), ('db', idx + 2),
                            ('dbtur', idx + 3)):
                if is_free_mask[j]:
                    new[name] = full_errors[j]
                elif j in linked_errors:
                    new[name] = linked_errors[j]
            comps.append(dataclasses.replace(c, **new))
        return SystemList(dataclasses.replace(self._system_list._data, components=comps))

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

    def _bound_flags(self, idx: int, p_free: np.ndarray,
                     bounds: Tuple[np.ndarray, np.ndarray]) -> Tuple[str, ...]:
        """Flags of the component whose parameters start at ``idx`` (in
        the full vector): its free z, logN and b at a bound (#517)."""
        is_free_mask = self._constraints._param_map['is_free']
        free_idx_map = np.cumsum(is_free_mask) - 1
        flags = []
        for k, (low, high) in enumerate((('z_at_window', 'z_at_window'),
                                         ('logN_at_lower_bound', 'logN_at_upper_bound'),
                                         ('b_at_lower_bound', 'b_at_upper_bound'))):
            if not is_free_mask[idx + k]:
                continue
            j = free_idx_map[idx + k]
            lo, hi, value = bounds[0][j], bounds[1][j], p_free[j]
            if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
                continue
            tol = AT_BOUND_RTOL * (hi - lo)
            if value - lo <= tol:
                flags.append(low)
            elif hi - value <= tol:
                flags.append(high)
        return tuple(dict.fromkeys(flags))

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
        idx_center: Any = np.searchsorted(self._x_ang, lambda_obs)
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

            if is_multiplet(comp.series): primary = series_lines(comp.series)[0]
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
            logger.debug(f"Initial guess: z = {new_z:.5f}, logN = {logN_guess:.2f}, "
                         f"b = {b_guess:.1f} km/s.")
        return p_refined

    def _voigt_optical_depth_vectorized(self, wave_grid_ang: np.ndarray,
                                        lambda_0: np.ndarray, f_val: np.ndarray,
                                        gamma: np.ndarray, z: np.ndarray,
                                        N_col: np.ndarray, b_kms: np.ndarray) -> np.ndarray:
        """See :func:`voigt_optical_depth`."""
        return voigt_optical_depth(wave_grid_ang, lambda_0, f_val, gamma, z, N_col, b_kms)

    def feature_mask(self, components: List[ComponentData],
                     z_window_kms: float = 20.0) -> np.ndarray:
        """The pixels a fit of ``components`` uses (those of
        :meth:`prepare_fit_context`): for each of their lines, the feature
        around it in the data (:meth:`_find_feature_limits`, out to the
        continuum or to a saddle between blended lines, the wings of a damped
        line included), or a fixed window of at least 150 km/s when the
        feature is narrower than 30 km/s; only the valid pixels. All the
        valid pixels without components."""
        c_kms = 299792.458
        if not components:
            return self._valid_mask
        mask = np.zeros_like(self._valid_mask, dtype=bool)
        for comp in components:
            for trans in _transitions_of(comp.series):
                lam_min, lam_max = self._find_feature_limits(
                    comp.z, trans, z_window_kms=z_window_kms, logN=comp.logN)

                # A window narrower than 30 km/s (a very narrow line, or
                # limits not found) is replaced by a fixed one of at
                # least 150 km/s, so that the fit of a very narrow
                # component still has pixels (its purpose as Guido
                # recalls it, November 2026, F111)
                use_fallback = True
                if lam_min is not None and lam_max is not None:
                    width_v = c_kms * (lam_max - lam_min) / lam_max
                    if width_v > 30.0:
                        use_fallback = False
                        mask |= (self._x_ang >= lam_min) & (self._x_ang <= lam_max)

                if use_fallback:
                    lam_0 = ATOM_DATA.get(trans, {}).get('wave')
                    if lam_0:
                        lam_obs = lam_0 * (1.0 + comp.z)
                        safe_window = max(z_window_kms, 150.0)
                        dw = (safe_window / c_kms) * lam_obs
                        mask |= (self._x_ang > lam_obs - dw) & (self._x_ang < lam_obs + dw)
        return mask & self._valid_mask

    def prepare_fit_context(self, z_window_kms: float = 20.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Prepares the internal state for fitting (masks, background, vectorized data).
        
        Returns
        -------
        tuple
            (lower_bounds, upper_bounds) for the free parameters.
        """
        c_kms = 299792.458
        active_uuids = self._constraints._active_uuids
        
        def get_trans_list(series):
            if is_multiplet(series): return series_lines(series)
            if series in ATOM_DATA: return [series]
            return []

        active_components = [c for c in self._system_list.components if c.uuid in active_uuids] if active_uuids else self._system_list.components

        self._fit_mask = self.feature_mask(active_components, z_window_kms)
        self._link_input_cache = None
        
        if np.any(self._fit_mask):
            indices = np.where(self._fit_mask)[0]
            self._sl = slice(max(0, np.min(indices) - 50), min(len(self._x_ang), np.max(indices) + 51))
        else:
            # It fitted the whole spectrum, with a WARNING: the components
            # (usually lines out of the spectrum) stayed at their first
            # values, and the fit took hours on a long spectrum (F111;
            # Guido's choice, phase 6-ter, PR 7-quater)
            listed = ", ".join(f"{c.series} at z = {c.z:.4f}" for c in active_components[:5])
            if len(active_components) > 5:
                listed += f" and {len(active_components) - 5} more"
            raise NoPixelsToFit(
                f"No valid pixels around the components to fit ({listed}): their lines "
                f"fall outside the spectrum ({self._x_ang.min() / 10:.1f}-"
                f"{self._x_ang.max() / 10:.1f} nm) or on pixels that cannot be used.")

        self._x_calc = self._x_ang[self._sl]
        self._y_norm_calc = self._y_norm[self._sl]
        self._dy_norm_calc = self._dy_norm[self._sl]
        self._fit_mask_calc = self._fit_mask[self._sl]
        

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
        grouped_lines: Dict[Any, Any] = {}
        
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
            if is_multiplet(comp.series): primary = series_lines(comp.series)[0]
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

    @overload
    def fit(self, max_nfev: int = ..., method: str = ..., z_window_kms: float = ...,
            verbose: int = ..., b_max: Optional[float] = ..., logN_max: Optional[float] = ...,
            jac: str = ..., model: Literal[True, 'fit_pixels'] = ...
            ) -> Tuple[SystemList, np.ndarray, Any]: ...

    @overload
    def fit(self, max_nfev: int = ..., method: str = ..., z_window_kms: float = ...,
            verbose: int = ..., b_max: Optional[float] = ..., logN_max: Optional[float] = ...,
            jac: str = ..., *, model: Literal[False]) -> Tuple[SystemList, None, Any]: ...

    def fit(self, max_nfev: int = 2000, method: str = 'trf',
            z_window_kms: float = 20.0, verbose: int = 0,
            b_max: Optional[float] = None,
            logN_max: Optional[float] = None, jac: str = 'analytic',
            model: Union[bool, str] = True) -> Tuple[SystemList, Optional[np.ndarray], Any]:
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
        jac : str, optional
            The Jacobian of the optimization: ``'analytic'`` (the default,
            phase 7, PR 7c.4d-bis), or a finite-difference scheme of
            ``scipy.optimize.least_squares`` (``'2-point'``, as before).
        model : bool, optional
            Return the model of the whole spectrum with all the components
            (the default). ``False`` for a caller that only needs the fitted
            list: the model is not computed and ``None`` is returned in its
            place (F167: on a long spectrum with many components it cost
            most of a short fit). ``'fit_pixels'``: the model only on the
            pixels of the fit (those of :meth:`chi2_on_fit_pixels`), NaN
            elsewhere, the same numbers there.

        Returns
        -------
        tuple
            A tuple containing:
            
            1. **new_system_list** (:class:`~astrocook.core.system_list.SystemList`):
               A new system list containing the optimized parameter values and errors.
            2. **final_model_flux** (*np.ndarray* or ``None``):
               The full model flux array (physical units), evaluated on the entire grid
               (``None`` with ``model=False``).
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
        self._system_list = SystemList(replace(old_systs._data, components=new_comps))

        bounds = self.prepare_fit_context(z_window_kms=z_window_kms)
        if b_max is not None:
            bounds = self._cap_b(bounds, b_max)
        if logN_max is not None:
            bounds = self._cap_upper(bounds, 1, logN_max)

        # Clamp p0 to bounds
        p0 = np.clip(p0, bounds[0], bounds[1])

        # Fit
        # The analytic Jacobian (phase 7, PR 7c.4d-bis, F171): before, finite
        # differences, one model per free parameter at each iteration
        jac_fn: Any = self._jacobian if (jac == 'analytic' and self._active_lines_data) else '2-point'
        res = least_squares(self._residual_function, p0, jac=jac_fn, bounds=bounds, method=method,
                            loss='linear', max_nfev=max_nfev, x_scale='jac', verbose=verbose)
        
        # Restore system list (it will be updated properly at the end)
        self._system_list = old_systs
        
        # --- CLEANUP (Reset Optimized Caches) ---
        self._active_fit_components = None
        self._active_lines_data = None

        # Statistics
        p_free_errors = np.zeros_like(res.x)
        is_free_mask = self._constraints._param_map['is_free']
        red_chi2: float = 0
        linked_errors = {}
        try:
            p_free_errors, red_chi2, cov = self._scaled_errors(res.fun, res.jac)
            linked_errors = self._linked_errors(res.x, cov)
            logger.debug(f"Errors calculated: {p_free_errors}")
        except (ValueError, ArithmeticError, KeyError, TypeError) as e:
            # The components with NaN errors are listed in a warning below
            logger.debug(f"The uncertainties could not be computed: {e}")
            p_free_errors = np.zeros_like(res.x)
            p_free_errors[:] = np.nan
            linked_errors = {idx: np.nan for idx, _ in self._constraints._compiled_links}

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
            def resolve_error(param_idx, old_err, is_spectator=is_spectator):
                if is_free_mask[param_idx]:
                    return p_full_errors[param_idx] # New calculated error
                elif is_spectator:
                    return old_err # PRESERVE old error for spectators
                elif param_idx in linked_errors:
                    return linked_errors[param_idx] # Propagated through the link (#630)
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
            
            comp_flags = old_comp.flags
            if active_uuids is None or old_comp.uuid in active_uuids:
                comp_flags = self._bound_flags(idx, res.x, bounds)
                if any(e is not None and np.isnan(e)
                       for e in (dz_new, dlogN_new, db_new, dbtur_new)):
                    comp_flags += ('nan_errors',)

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
                resol=comp_resol,
                flags=comp_flags
            )
            new_components.append(new_comp)
            
        # Report the components whose uncertainties could not be computed
        # (findings register F20): the errors stay NaN, which the tables show
        nan_errors = [c for c, old in zip(new_components, old_components, strict=True)
                      if (active_uuids is None or old.uuid in active_uuids)
                      and any(e is not None and np.isnan(e)
                              for e in (c.dz, c.dlogN, c.db, c.dbtur))]
        if nan_errors:
            listed = ", ".join(f"{c.series} at z = {c.z:.5f}" for c in nan_errors[:10])
            more = f" and {len(nan_errors) - 10} more" if len(nan_errors) > 10 else ""
            logger.warning(
                f"Uncertainties could not be computed (ill-conditioned fit) for "
                f"{len(nan_errors)} component(s): {listed}{more}. "
                f"Their errors are NaN; consider refitting or freezing parameters.")

        new_data_core = dataclasses.replace(self._system_list._data, components=new_components)
        new_system_list = SystemList(new_data_core)
        
        # 4. Final Flux Model (Un-sliced)
        self._cached_tau_groups = {} # Clear cache
        if not model:
            logger.debug(f"Fit complete. Reduced chi2 (constraining pixels): {red_chi2:.2f}.")
            return new_system_list, None, res

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
        if isinstance(model, str):
            self._restrict_model(self._fit_mask)

        # This will now loop over ALL components because active_uuids is None
        try:
            final_model_norm = self._compute_model(res.x)
        finally:
            # Restore state (good practice)
            self._model_pixels = None
            self._constraints._active_uuids = stored_active_uuids
            self._x_calc, self._resol_calc, self._x_tau, self._lsf = saved
        
        final_model_flux = final_model_norm * self._norm

        logger.debug(f"Fit complete. Reduced chi2 (constraining pixels): {red_chi2:.2f}.")
        return new_system_list, final_model_flux, res

    def normalized_residuals(self, model_flux: np.ndarray, mask: np.ndarray) -> np.ndarray:
        """(flux - model) / error on the pixels ``mask``, for a model flux in
        the units of :meth:`fit` and :meth:`compute_model_flux` (on the whole
        spectrum). The error is the one of the fits (bounded from below, a
        continuum <= 0 gives NaN). The one computation of the residuals of a
        model on a fixed set of pixels: :meth:`chi2_on_fit_pixels`
        (``metals_auto``, F111) and the evaluation of a model (F154, phase 7,
        PR 7a.5)."""
        return (self._y_norm[mask] - model_flux[mask] / self._norm[mask]) / self._dy_norm[mask]

    @property
    def fit_mask(self) -> np.ndarray:
        """The pixels of the last fit (on the whole spectrum)."""
        return self._fit_mask

    @property
    def valid_mask(self) -> np.ndarray:
        """The pixels with a finite normalized flux and a positive error: the
        only ones a fit or an evaluation can use."""
        return self._valid_mask

    def chi2_on_fit_pixels(self, model_flux: np.ndarray) -> Tuple[float, int]:
        """The chi2 of a model flux (in the units of :meth:`fit` and
        :meth:`compute_model_flux`, on the whole spectrum) on the pixels of
        the last fit, and their number: two models compared on the same
        pixels (``metals_auto``, F111)."""
        mask = self._fit_mask
        return chi2_of(self.normalized_residuals(model_flux, mask)), int(np.sum(mask))

    def component_profiles(self, uuids: List[str], sl: slice) -> Dict[str, np.ndarray]:
        """The normalized flux of each component of ``uuids`` alone, at the
        current parameters, convolved with its resolution and averaged over
        each pixel of ``sl`` (the profile I^a of Webb et al. 2021, eq. 6, as
        the pixels see it: their model, eq. 7, is the product of these
        profiles, which is exact only before the convolution)."""
        p_full = self._constraints.map_p_free_to_full(self._constraints.p_free_vector)
        saved = (self._x_calc, self._resol_calc, self._x_tau, self._lsf)
        self._x_calc = self._x_ang[sl]
        self._resol_calc = self._resol_column[sl] if self._use_variable_resolution else None
        self._set_lsf(sl, cache=False)
        index = {c.uuid: i for i, c in enumerate(self._system_list.components)}
        profiles = {}
        try:
            for uuid in uuids:
                i = index[uuid]
                comp = self._system_list.components[i]
                z, logN, b_val, btur = p_full[4 * i:4 * i + 4]
                tau = self._compute_tau_component(comp, z, logN, b_val, btur)
                profiles[uuid] = self._convolved(self._get_component_resolution(comp), tau)
        finally:
            self._x_calc, self._resol_calc, self._x_tau, self._lsf = saved
        return profiles

    def _convolved(self, res_key: Tuple[float, str], tau_arr: np.ndarray) -> np.ndarray:
        """The flux exp(-tau) of the fine grid convolved with the LSF of the
        resolution ``res_key`` and averaged over the pixels (#526)."""
        res_val, res_unit = res_key
        flux_fine = np.exp(-tau_arr)
        pix = None
        if self._model_pixels is not None:
            pix, j0, j1 = self._model_pixels
            flux_fine = np.ones(len(self._lsf.x_fine))
            flux_fine[j0:j1] = np.exp(-tau_arr)
        r: Any
        key: Any
        if res_unit == 'column' and self._resol_calc is not None:
            r, key = self._resol_calc, 'column'
        elif res_val > 0:
            r, key = res_val, res_val
        else:
            r, key = None, 'none'
        return self._lsf.apply(flux_fine, r, key=key, pixels=pix)

    def _restrict_model(self, pixels: np.ndarray) -> None:
        """Compute the model of the whole spectrum (set up by
        ``_set_lsf(slice(None), cache=False)``) only on the pixels
        ``pixels`` (a mask), the others NaN (F167): the optical depth only
        on the fine samples they take, the same samples of the whole grid,
        so the same numbers (up to the order of the sums)."""
        pix = np.flatnonzero(pixels)
        j0, j1 = self._lsf.fine_range(pix)
        self._x_tau = self._lsf.x_fine[j0:j1]
        self._model_pixels = (pix, j0, j1)

    def compute_model_flux(self, pixels: Optional[np.ndarray] = None
                           ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the model flux for the current parameters.

        Calculates the Voigt profile model based on the current state of the
        system list, without performing any fitting/optimization. With
        ``pixels`` (a boolean mask on the spectrum), only on those pixels,
        NaN elsewhere: the same numbers on them, for a caller that needs no
        more (F167, ``metals_auto``).

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
        if pixels is not None:
            self._restrict_model(pixels)

        # Reset caches to force Case B (Standard Loop) in _compute_model
        self._cached_tau_groups = {}
        self._active_lines_data = None 
        
        try:
            y_model_norm = self._compute_model(p_current_free)
        finally:
            # Restore context
            self._model_pixels = None
            self._x_calc, self._resol_calc, self._x_tau, self._lsf = saved
        
        y_model_phys = y_model_norm * self._norm
        return self._x_ang, y_model_phys

# Old names (phase 6-quater, PR 6): see astrocook/_renamed.py
__getattr__ = renamed(__name__, {"VoigtFitterV2": "VoigtFitter"})
