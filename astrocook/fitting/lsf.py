"""Exact convolution of a model with the line spread function (#526).

The model flux F is computed on a fine grid, uniform in velocity, and each
pixel i of the spectrum gets the mean over its own edges [v_lo_i, v_hi_i]
of F convolved with a Gaussian of the resolution *of that pixel*::

    model_i = (1 / (v_hi_i - v_lo_i)) * integral over [v_lo_i, v_hi_i] of (F * G_i)(v) dv
            = integral of F(v) K_i(v) dv / integral of K_i(v) dv,
    K_i(v) = Phi((v_hi_i - v) / sigma_i) - Phi((v_lo_i - v) / sigma_i),

with Phi the Gaussian cumulative distribution and
sigma_i = c / R_i / (2 sqrt(2 ln 2)). On the fine samples F_j (step 2h)
the integral is a sum with weights w_ij, normalised over j:

- when the kernel is resolved (sigma_i at least one step), w_ij = K_i(v_j):
  sampling a smooth kernel converges very fast (1e-8 of the continuum
  against the analytic result for Gaussian lines);
- without a resolution, or with sigma_i below one step, the kernel has
  sharp edges and its samples would shift each pixel by up to h, depending
  on how they fall within it (3e-3 of the continuum on a line with
  b = 12 km/s). The flux is then taken constant over each cell
  [v_j - h, v_j + h] and the weight is the exact integral of K_i over it,

      w_ij = r(v_hi_i - v_j + h) - r(v_hi_i - v_j - h)
             - r(v_lo_i - v_j + h) + r(v_lo_i - v_j - h),

  with r(t) = sigma_i G(t / sigma_i), G(u) = u Phi(u) + phi(u); without a
  resolution, the overlap of the cell with the pixel. The error is then
  that of a midpoint rule, ~h^2/24 F'' (up to 3e-3 on a line with
  b = 1 km/s).

This replaces the convolution on the pixel grid with one kernel width for
the whole window (taken at its median wavelength) and no integration over
the pixel (findings register F10, F9).
"""
from typing import Dict, Optional, Tuple, Union

import numpy as np
from scipy.special import ndtr

C_KMS = 299792.458
FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

#: Fine step (#526, F9): half the instrumental sigma, so that the kernel
#: is resolved, and a third of the b of the narrowest line in the pixels
#: (within 1e-5 of the continuum against a mock with 60 sub-pixels, also
#: on saturated lines; with b / 2, up to 4e-4), never below half the
#: smallest b the fitter allows (1 km/s, ``voigt_model``; 2e-4 there) nor
#: above the smallest pixel. Without a resolution, also OVERSAMPLE samples
#: per pixel. Sampling at the pixels alone (before #526) was off by up to
#: 6e-3 at R 45000 and 6e-2 at R 8900
SIGMA_SAMPLES = 2.0
B_SAMPLES = 3.0
MIN_STEP_KMS = 0.5
OVERSAMPLE = 5
#: Pixels per block when the weights are not cached (whole spectra): the
#: weights take n x K values, with K up to a few hundred samples per pixel
#: (low resolution, narrow lines), ~10 MB per array with 2000 pixels
CHUNK = 2000
#: The kernel is cut at this many sigma
KERNEL_SIGMAS = 5.0


def _ramp(t: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """Integral of the Gaussian cumulative distribution:
    sigma * G(t / sigma), G(u) = u Phi(u) + phi(u); max(t, 0) for
    sigma = 0. The weight of a fine cell [v - h, v + h] for a pixel
    [v_lo, v_hi] is a combination of four of these (see ``_weights``)."""
    with np.errstate(invalid='ignore', divide='ignore'):
        u = t / sigma
        smooth = sigma * (u * ndtr(u) + np.exp(-0.5 * u**2) / np.sqrt(2.0 * np.pi))
    return np.where(sigma > 0, smooth, np.maximum(t, 0.0))


class PixelLSF:
    """Fine grid and exact LSF convolution for a set of pixels.

    Parameters
    ----------
    xmin_ang, xmax_ang : np.ndarray
        Edges of the pixels (Angstrom), increasing.
    r_min, r_max : float, optional
        Smallest and largest resolving power to be used: the fine grid
        extends beyond the outer pixels by ``KERNEL_SIGMAS`` sigma of
        ``r_min``, and its step follows the sigma of ``r_max``.
    b_min : float, optional
        Doppler parameter (km/s) of the narrowest line in the pixels.
    cache : bool
        Keep the weights of each resolution, for the many calls of a fit.
        Without, the weights are computed at each call only for the pixels
        that have absorption within the reach of their kernel (the others
        are 1), for a model of the whole spectrum computed once.
    """

    def __init__(self, xmin_ang: np.ndarray, xmax_ang: np.ndarray,
                 r_min: Optional[float] = None, r_max: Optional[float] = None,
                 b_min: Optional[float] = None, cache: bool = True):
        xmin = np.asarray(xmin_ang, dtype=float)
        xmax = np.asarray(xmax_ang, dtype=float)
        self.n_pix = len(xmin)
        self._lam_ref = float(np.median(np.sqrt(xmin * xmax))) if self.n_pix else 1.0
        self.v_lo = C_KMS * np.log(xmin / self._lam_ref)
        self.v_hi = C_KMS * np.log(xmax / self._lam_ref)
        widths = self.v_hi - self.v_lo
        good = widths[np.isfinite(widths) & (widths > 0)]
        min_width = float(np.min(good)) if good.size else 1.0
        steps = [min_width]
        if r_max and r_max > 0:
            steps.append(C_KMS / r_max * FWHM_TO_SIGMA / SIGMA_SAMPLES)
        else:
            steps.append(min_width / OVERSAMPLE)
        if b_min is not None and b_min > 0:
            steps.append(max(b_min / B_SAMPLES, MIN_STEP_KMS))
        self.dv = min(steps)
        sigma_max = (C_KMS / r_min * FWHM_TO_SIGMA) if (r_min and r_min > 0) else 0.0
        margin = KERNEL_SIGMAS * sigma_max + self.dv
        if self.n_pix:
            v_start = np.nanmin(self.v_lo) - margin
            v_end = np.nanmax(self.v_hi) + margin
        else:
            v_start = v_end = 0.0
        n_fine = int(np.ceil((v_end - v_start) / self.dv)) + 1
        self.v_fine = v_start + self.dv * np.arange(n_fine)
        #: Fine grid in Angstrom, on which the optical depth is computed
        self.x_fine = self._lam_ref * np.exp(self.v_fine / C_KMS)
        self._use_cache = cache
        self._cache: Dict[Union[float, str], Tuple[np.ndarray, np.ndarray]] = {}

    def _weights(self, sigma: np.ndarray, sl=slice(None)
                 ) -> Tuple[np.ndarray, np.ndarray]:
        """Indices (n, K) into the fine grid and normalised weights, for
        the pixels ``sl`` (slice or indices)."""
        v_lo, v_hi = self.v_lo[sl], self.v_hi[sl]
        h = 0.5 * self.dv
        reach = KERNEL_SIGMAS * np.nan_to_num(sigma, nan=0.0) + h
        j_lo = np.searchsorted(self.v_fine, v_lo - reach, side='left')
        j_hi = np.searchsorted(self.v_fine, v_hi + reach, side='right')
        k = max(1, int(np.max(j_hi - j_lo)) if len(v_lo) else 1)
        idx = j_lo[:, None] + np.arange(k)[None, :]
        inside = idx < j_hi[:, None]
        idx = np.minimum(idx, len(self.v_fine) - 1)
        v = self.v_fine[idx]
        s = np.where(np.isfinite(sigma) & (sigma > 0), sigma, 0.0)[:, None]
        resolved = s[:, 0] >= self.dv
        w = np.empty_like(v)
        if np.any(resolved):
            r, sr = resolved, s[resolved]
            w[r] = (ndtr((v_hi[r, None] - v[r]) / sr) - ndtr((v_lo[r, None] - v[r]) / sr))
        if not np.all(resolved):
            c, sc = ~resolved, s[~resolved]
            vc, lo_c, hi_c = v[c], v_lo[c, None], v_hi[c, None]
            w[c] = (_ramp(hi_c - vc + h, sc) - _ramp(hi_c - vc - h, sc)
                    - _ramp(lo_c - vc + h, sc) + _ramp(lo_c - vc - h, sc))
        w = np.where(inside, np.maximum(w, 0.0), 0.0)
        norm = w.sum(axis=1, keepdims=True)
        # A pixel outside the fine grid (it should not happen) takes the
        # nearest sample
        empty = norm[:, 0] <= 0
        if np.any(empty):
            centre = 0.5 * (v_lo[empty] + v_hi[empty])
            j = np.clip(np.searchsorted(self.v_fine, centre), 0, len(self.v_fine) - 1)
            idx[empty, 0] = j
            w[empty] = 0.0
            w[empty, 0] = 1.0
            norm[empty] = 1.0
        return idx, w / norm

    def apply(self, flux_fine: np.ndarray, resol: Union[None, float, np.ndarray],
              key=None) -> np.ndarray:
        """Flux of each pixel: the fine flux convolved with the LSF of the
        pixel (resolving power ``resol``, scalar or one per pixel; None or
        NaN for no LSF) and averaged over the pixel. The weights are cached
        under ``key`` (the resolution itself if scalar or None), unless the
        pixels are more than ``CHUNK``."""
        if resol is None:
            r = np.full(self.n_pix, np.nan)
            key = 'none' if key is None else key
        elif np.ndim(resol) == 0:
            r = np.full(self.n_pix, float(resol))
            key = float(resol) if key is None else key
        else:
            r = np.asarray(resol, dtype=float)
        with np.errstate(invalid='ignore', divide='ignore'):
            sigma = np.where(np.isfinite(r) & (r > 0), C_KMS / r * FWHM_TO_SIGMA, np.nan)
        if not self._use_cache:
            return self._apply_absorbed(flux_fine, sigma)
        if self.n_pix > CHUNK:
            out = np.empty(self.n_pix)
            for start in range(0, self.n_pix, CHUNK):
                sl = slice(start, min(start + CHUNK, self.n_pix))
                idx, w = self._weights(sigma[sl], sl)
                out[sl] = np.sum(flux_fine[idx] * w, axis=1)
            return out
        if key is None or key not in self._cache:
            weights = self._weights(sigma)
            if key is None:
                idx, w = weights
                return np.sum(flux_fine[idx] * w, axis=1)
            self._cache[key] = weights
        idx, w = self._cache[key]
        return np.sum(flux_fine[idx] * w, axis=1)

    def _apply_absorbed(self, flux_fine: np.ndarray, sigma: np.ndarray) -> np.ndarray:
        """As ``apply``, without cache: only the pixels with absorption
        (fine flux below 1) within the reach of their kernel are computed;
        the others are exactly 1."""
        out = np.ones(self.n_pix)
        absorbed = np.concatenate([[0], np.cumsum(flux_fine < 1.0)])
        if absorbed[-1] == 0:
            return out
        reach = KERNEL_SIGMAS * np.nan_to_num(sigma, nan=0.0) + 0.5 * self.dv
        j_lo = np.searchsorted(self.v_fine, self.v_lo - reach, side='left')
        j_hi = np.searchsorted(self.v_fine, self.v_hi + reach, side='right')
        pix = np.flatnonzero(absorbed[j_hi] > absorbed[j_lo])
        for start in range(0, len(pix), CHUNK):
            p = pix[start:start + CHUNK]
            idx, w = self._weights(sigma[p], p)
            out[p] = np.sum(flux_fine[idx] * w, axis=1)
        return out
