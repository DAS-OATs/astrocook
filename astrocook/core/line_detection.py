"""Detection of the absorption lines by matched filters (phase 7, PR 7c.4b).

The direction of F57 and of section 5 of ``design/identification.md``:
the lines are searched on the flux normalised by the continuum with a bank
of filters matched to the expected line widths, instead of the pixels that
the clipping of the continuum rejects (``abs_mask``, which stays the mask
of the continuum).

For each pixel and each width the significance of a Gaussian absorption
centred there is

    S = sum_i w_i (1 - F_i) / sigma_i^2 / sqrt(sum_i w_i^2 / sigma_i^2),

with ``w`` the filter, F the normalised flux and sigma its error: the
amplitude of the filter fitted to the data by least squares, over its
error. On pure noise S is a standard normal variable at each pixel
(correlated between neighbours). The filter is a Gaussian in velocity
whose variance is that of the instrumental profile plus that of a line of
Doppler parameter b, b^2 / 2 (an optically thin line); three widths by
default, b = 5, 15, 40 km/s.

A *line* is a pixel where S > :data:`LINE_SIGMA` for any width; its
*region* the run of pixels around it where S > :data:`EXTENT_SIGMA` for
any width (a hysteresis threshold). The multiplets are also searched on
their joint significance at their lags: for the strongest line (highest f
lambda) at a pixel, the filters of the other lines are scaled by their f
lambda over that of the strongest (the ratio of their optical depths for
a thin line) and placed at their lags,

    J = sum_k r_k N_k / sqrt(sum_k r_k^2 D_k),

with N and D the numerator and the squared denominator of S of each
line. A multiplet where J > :data:`LINE_SIGMA`, each of its lines above
:data:`EXTENT_SIGMA` alone and the lines consistent with the ratio
(:func:`multiplet_seeds`) gives a line at each of its lines, also where no
line reaches the threshold alone (a doublet at 3 and 1.5 sigma has J of
about 3.35); their regions are found as those of the other lines.

Not used by the pipelines yet (PR 7c.4c will): no recipe changes its
numbers.
"""
import dataclasses
import logging
from typing import Dict, Iterable, Optional, Tuple

import astropy.units as au
import numpy as np
from scipy.ndimage import label, maximum_filter1d

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS

logger = logging.getLogger(__name__)

C_KMS = 299792.458

#: The Doppler parameters of the filters (km/s): a narrow metal line, a
#: typical one or a forest line, a broad one or a blend (section 5 of
#: design/identification.md, Guido, 18 November 2026)
DEFAULT_B_KMS = (5.0, 15.0, 40.0)
#: Significance of a line, for any width, and of a multiplet at its lags
#: (Guido, 18 November 2026: 3 sigma per line and on the joint
#: significance of a doublet)
LINE_SIGMA = 3.0
#: Significance that bounds the region of a line: the pixels around it
#: where the absorption is still above the noise
EXTENT_SIGMA = 1.0
#: The multiplets searched on their joint significance by default: the
#: anchors of the mock benchmark (the identification of PR 7c.4c chooses
#: its own)
DEFAULT_MULTIPLETS = ("CIV", "SiIV", "MgII")
#: Largest chi2 of the lines of a multiplet about the ratio of their f
#: lambda (:func:`multiplet_seeds`): 3 sigma for the one degree of freedom
#: of a doublet, whatever the threshold of J
CONSISTENCY_CHI2 = 9.0
#: Half width of a filter, in its standard deviations (the Gaussian is
#: below 1.1% of its peak beyond it)
KERNEL_SIGMAS = 3.0
#: Largest number of filter values held at once (pixels times filter
#: length), to bound the memory on long spectra
_CHUNK_VALUES = 4_000_000

_FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


@dataclasses.dataclass(frozen=True)
class LineDetection:
    """The result of :func:`matched_filter_detection`.

    Attributes
    ----------
    b_kms : tuple of float
        The Doppler parameters of the filters.
    significance : np.ndarray
        S of each filter at each pixel, shape (len(b_kms), n); NaN where
        the filter holds no valid pixel.
    joint : dict
        Multiplet -> its joint significance J, shape (len(b_kms), n), with
        its strongest line at each pixel (NaN where a line falls outside
        the spectrum, or beyond ``z_max``).
    lines : np.ndarray
        True at the pixels of a line: S > ``line_sigma`` for any width, or
        a line of a multiplet with J > ``line_sigma``.
    mask : np.ndarray
        True in the regions of the lines.
    regions : np.ndarray
        The regions, numbered from 1 (0 outside).
    """
    b_kms: Tuple[float, ...]
    significance: np.ndarray
    joint: Dict[str, np.ndarray]
    lines: np.ndarray
    mask: np.ndarray
    regions: np.ndarray

    @property
    def n_regions(self) -> int:
        """The number of regions."""
        return int(self.regions.max(initial=0))

    @property
    def max_significance(self) -> np.ndarray:
        """S of the most significant width at each pixel (NaN where no
        filter holds a valid pixel)."""
        return _nanmax(self.significance)


def _nanmax(values: np.ndarray) -> np.ndarray:
    """Maximum over the first axis, NaN where all the values are NaN
    (without the warning of numpy)."""
    finite = np.isfinite(values)
    out = np.where(finite, values, -np.inf).max(axis=0)
    return np.where(finite.any(axis=0), out, np.nan)


def filter_sigma_kms(resolving_power: np.ndarray, b_kms: float) -> np.ndarray:
    """Standard deviation (km/s) of the filter of a line of Doppler
    parameter ``b_kms`` at each pixel: the instrumental one (FWHM c / R)
    and b / sqrt(2) in quadrature."""
    instr = C_KMS / np.asarray(resolving_power, dtype=float) * _FWHM_TO_SIGMA
    return np.sqrt(instr ** 2 + b_kms ** 2 / 2.0)


def filter_sums(v_kms: np.ndarray, depth: np.ndarray, ivar: np.ndarray,
                sigma_kms: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """The sums of the matched filter centred at each pixel.

    Parameters
    ----------
    v_kms : np.ndarray
        The pixels in velocity (c ln lambda, increasing).
    depth : np.ndarray
        1 - F (0 where the pixel is not valid).
    ivar : np.ndarray
        1 / sigma^2 of 1 - F (0 where the pixel is not valid).
    sigma_kms : np.ndarray
        Standard deviation of the filter centred at each pixel.

    Returns
    -------
    (np.ndarray, np.ndarray)
        N = sum w depth ivar and D = sum w^2 ivar; S = N / sqrt(D).
    """
    n = len(v_kms)
    dv = np.abs(np.gradient(v_kms)) if n > 1 else np.ones(1)
    dv = np.where(dv > 0, dv, np.nan)
    with np.errstate(invalid="ignore"):
        reach = np.nanmax(KERNEL_SIGMAS * sigma_kms / dv) if n else 0.0
    k = max(1, int(np.ceil(reach))) if np.isfinite(reach) else 1
    pad = np.zeros(k)
    v_pad = np.concatenate((np.full(k, v_kms[0]), v_kms, np.full(k, v_kms[-1])))
    a_pad = np.concatenate((pad, depth * ivar, pad))
    i_pad = np.concatenate((pad, ivar, pad))
    view_v = np.lib.stride_tricks.sliding_window_view(v_pad, 2 * k + 1)
    view_a = np.lib.stride_tricks.sliding_window_view(a_pad, 2 * k + 1)
    view_i = np.lib.stride_tricks.sliding_window_view(i_pad, 2 * k + 1)
    num = np.empty(n)
    den = np.empty(n)
    step = max(1, _CHUNK_VALUES // (2 * k + 1))
    for lo in range(0, n, step):
        hi = min(n, lo + step)
        u = (view_v[lo:hi] - v_kms[lo:hi, None]) / sigma_kms[lo:hi, None]
        w = np.exp(-0.5 * u ** 2)
        w[np.abs(u) > KERNEL_SIGMAS] = 0.0
        num[lo:hi] = np.sum(w * view_a[lo:hi], axis=1)
        den[lo:hi] = np.sum(w ** 2 * view_i[lo:hi], axis=1)
    return num, den


def _significance(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / np.sqrt(den), np.nan)


def multiplet_lines(multiplet: str) -> Tuple[np.ndarray, np.ndarray]:
    """Lags (km/s, from the strongest line) and f lambda over that of the
    strongest line, of the lines of a standard multiplet."""
    lines = STANDARD_MULTIPLETS[multiplet]
    flam = np.array([ATOM_DATA[t]["f"] * ATOM_DATA[t]["wave"] for t in lines])
    wave = np.array([ATOM_DATA[t]["wave"] for t in lines], dtype=float)
    ref = int(np.argmax(flam))
    return C_KMS * np.log(wave / wave[ref]), flam / flam[ref]


def joint_significance(v_kms: np.ndarray, num: np.ndarray, den: np.ndarray,
                       lags_kms: np.ndarray, ratios: np.ndarray
                       ) -> Tuple[np.ndarray, np.ndarray]:
    """J of a multiplet with its strongest line at each pixel, and S of
    each of its lines, from the sums of one filter (:func:`filter_sums`)
    taken at the lags of the lines (linear interpolation in velocity); NaN
    where a line falls outside the spectrum.

    The lines are taken as independent: their filters must not overlap
    (true for the standard doublets, lags of 260 km/s and more, with
    filters of at most about 100 km/s).

    Returns
    -------
    (np.ndarray, np.ndarray)
        J, shape (n,), and S of the lines, shape (len(lags_kms), n).
    """
    tot_num = np.zeros(len(v_kms))
    tot_den = np.zeros(len(v_kms))
    s_lines = []
    for lag, r in zip(lags_kms, ratios, strict=True):
        at = v_kms + lag
        nu = np.interp(at, v_kms, num, left=np.nan, right=np.nan)
        de = np.interp(at, v_kms, den, left=np.nan, right=np.nan)
        tot_num += r * nu
        tot_den += r ** 2 * de
        s_lines.append(_significance(nu, de))
    return _significance(tot_num, tot_den), np.array(s_lines)


def multiplet_seeds(joint: np.ndarray, s_lines: np.ndarray, line_sigma: float,
                    extent_sigma: float) -> np.ndarray:
    """The pixels of the strongest line of a multiplet found by its joint
    significance (one filter): J > ``line_sigma``, every line above
    ``extent_sigma`` alone, and the lines consistent with the ratio of
    their f lambda, chi2 = sum S_k^2 - J^2 <= :data:`CONSISTENCY_CHI2` (the
    residuals of the amplitude fitted to all the lines; one degree of
    freedom less than the lines).

    The last two conditions keep a strong line from making a multiplet
    with the noise at its lags (a ghost at each lag of each multiplet
    around each strong line: on the mocks of the benchmark the whole
    forest and a quarter of the red side were in a region without them).
    A saturated multiplet, whose ratio is not that of f lambda, fails the
    consistency, but its lines are found alone."""
    with np.errstate(invalid="ignore"):
        chi2 = np.sum(s_lines ** 2, axis=0) - joint ** 2
        return ((joint > line_sigma) & np.all(s_lines > extent_sigma, axis=0)
                & (chi2 <= CONSISTENCY_CHI2))


def _runs_with_seeds(extent: np.ndarray, seeds: np.ndarray) -> np.ndarray:
    """The runs of ``extent`` that hold a pixel of ``seeds`` (the seeds are
    also in their own run)."""
    runs, _ = label(extent | seeds)
    keep = np.unique(runs[seeds & (runs > 0)])
    return np.isin(runs, keep) & (runs > 0)


def _project(v_kms: np.ndarray, seeds: np.ndarray, lag: float) -> np.ndarray:
    """The pixels nearest to the seeds moved by ``lag`` (those that fall
    inside the spectrum)."""
    out = np.zeros(len(v_kms), dtype=bool)
    at = v_kms[seeds] + lag
    at = at[(at >= v_kms[0]) & (at <= v_kms[-1])]
    if len(at) == 0:
        return out
    hi = np.clip(np.searchsorted(v_kms, at), 1, len(v_kms) - 1)
    nearer_lo = np.abs(at - v_kms[hi - 1]) <= np.abs(v_kms[hi] - at)
    out[np.where(nearer_lo, hi - 1, hi)] = True
    return out


def _multiplet_lines(v_kms: np.ndarray, sums, peak_half, name: str,
                     line_sigma: float, extent_sigma: float,
                     beyond: Optional[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    """J of a multiplet for each filter (NaN at the pixels ``beyond``),
    and the pixels of its lines where :func:`multiplet_seeds` finds it at
    a peak of J (its maximum within ``peak_half`` pixels, about one
    standard deviation of the filter: in the wing of a strong line J
    rises towards the line, and the wing with the noise at a lag is not a
    multiplet)."""
    lags, ratios = multiplet_lines(name)
    joints = []
    seeds = np.zeros(len(v_kms), dtype=bool)
    for (nu, de), half in zip(sums, peak_half, strict=True):
        j, s_lines = joint_significance(v_kms, nu, de, lags, ratios)
        if beyond is not None:
            j[beyond] = np.nan
        joints.append(j)
        j_fill = np.where(np.isfinite(j), j, -np.inf)
        peak = j_fill >= maximum_filter1d(j_fill, 2 * half + 1, mode="nearest")
        seeds |= multiplet_seeds(j, s_lines, line_sigma, extent_sigma) & peak
    lines = np.zeros(len(v_kms), dtype=bool)
    for lag in lags:
        lines |= _project(v_kms, seeds, lag)
    return np.array(joints).reshape(len(sums), len(v_kms)), lines


def matched_filter_detection(
        x_nm: np.ndarray, flux: np.ndarray, err: np.ndarray,
        resolving_power, b_kms: Iterable[float] = DEFAULT_B_KMS,
        line_sigma: float = LINE_SIGMA, extent_sigma: float = EXTENT_SIGMA,
        multiplets: Iterable[str] = DEFAULT_MULTIPLETS,
        z_max: Optional[float] = None,
        multiplet_sigma: Optional[float] = None) -> LineDetection:
    """The lines and their regions by matched filters (see the module).

    Parameters
    ----------
    x_nm : np.ndarray
        The wavelengths of the pixels (nm, increasing).
    flux, err : np.ndarray
        The flux normalised by the continuum and its error; the pixels with
        a value or an error not finite, or an error not positive, are left
        out of the sums.
    resolving_power : float or np.ndarray
        R of each pixel (or one for all).
    b_kms : iterable of float
        The Doppler parameters of the filters.
    line_sigma, extent_sigma : float
        The thresholds of a line and of its region.
    multiplets : iterable of str
        Standard multiplets (``STANDARD_MULTIPLETS``) searched on their
        joint significance (:func:`multiplet_seeds`); empty for none.
    z_max : float, optional
        No multiplet with its strongest line beyond this redshift (e.g.
        z_em and the margin of the associated absorbers, F55).
    multiplet_sigma : float, optional
        The threshold of J, if not ``line_sigma`` (each multiplet searched
        is a new set of trials: on pure noise the three doublets of
        :data:`DEFAULT_MULTIPLETS` at 3 sigma add several times the
        regions of the lines alone).

    Returns
    -------
    LineDetection
    """
    x_nm = np.asarray(x_nm, dtype=float)
    flux = np.asarray(flux, dtype=float)
    err = np.asarray(err, dtype=float)
    n = len(x_nm)
    r = np.broadcast_to(np.asarray(resolving_power, dtype=float), (n,))
    if not np.all(np.isfinite(r) & (r > 0)):
        raise ValueError("The matched filters need the resolution of every pixel "
                         "(the 'resol' column).")
    valid = np.isfinite(flux) & np.isfinite(err) & (err > 0)
    ivar = np.where(valid, 1.0 / np.where(valid, err, 1.0) ** 2, 0.0)
    depth = np.where(valid, 1.0 - np.where(valid, flux, 1.0), 0.0)
    v = C_KMS * np.log(x_nm)
    b_kms = tuple(float(b) for b in b_kms)
    sigmas = [filter_sigma_kms(r, b) for b in b_kms]
    sums = [filter_sums(v, depth, ivar, sig) for sig in sigmas]
    dv = np.abs(np.gradient(v)) if n > 1 else np.ones(n)
    peak_half = [max(1, int(round(float(np.median(sig / dv))))) if n else 1
                 for sig in sigmas]
    significance = np.array([_significance(nu, de) for nu, de in sums]).reshape(len(b_kms), n)
    best = _nanmax(significance)
    with np.errstate(invalid="ignore"):
        lines = best > line_sigma
    joint = {}
    for name in multiplets:
        beyond = None
        if z_max is not None:
            _, ratios = multiplet_lines(name)
            ref = STANDARD_MULTIPLETS[name][int(np.argmax(ratios))]
            beyond = x_nm / (ATOM_DATA[ref]["wave"] / 10.0) - 1.0 > z_max
        joint[name], m_lines = _multiplet_lines(
            v, sums, peak_half, name,
            line_sigma if multiplet_sigma is None else multiplet_sigma,
            extent_sigma, beyond)
        lines |= m_lines
    with np.errstate(invalid="ignore"):
        mask = _runs_with_seeds(best > extent_sigma, lines)
    regions, _ = label(mask)
    return LineDetection(b_kms=b_kms, significance=significance, joint=joint,
                         lines=lines, mask=mask, regions=regions)


def detect_lines(spec, **kwargs) -> LineDetection:
    """:func:`matched_filter_detection` on a spectrum: its flux and error
    over its normalization (``norm``: the continuum, times the telluric
    model if any), its resolution (``resol``). ``kwargs`` go to
    :func:`matched_filter_detection`. ``ValueError`` if the spectrum has
    no continuum or no resolution."""
    from astrocook.core.spectrum_lines import _normalized_flux, _wavelength_nm
    x_nm = _wavelength_nm(spec)
    flux = _normalized_flux(spec, "the spectrum")
    norm = np.asarray(au.Quantity(spec.norm).to_value(spec.y.unit), dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        err = np.asarray(spec.dy.to_value(spec.y.unit), dtype=float) / norm
    return matched_filter_detection(x_nm, flux, err, spec.resolution(), **kwargs)
