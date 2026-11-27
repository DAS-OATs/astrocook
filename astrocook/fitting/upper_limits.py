"""The computation of the upper limits on the column of a species not detected
(phase 7, PR 7d.7; findings register F156; ``ROADMAP.md``, phase 7, "Upper
limits").

Two methods:

- ``profile`` (the default): the profile likelihood on fixed pixels. A
  component of the species is added at the redshift, with b fixed and its
  log N fixed at each step; the components with a line on the pixels are
  refitted at each step, from their values in the list (so each step is a
  function of log N alone, not of the order of the steps); the chi2 is
  measured on fixed pixels. The limit is the log N above the minimum of the
  chi2 where it rises by ``sigma``^2: blends are accounted for, as the
  refitted components absorb what they can.
- ``ew``: the quick estimate from the noise: the error of the equivalent
  width on the pixels within one FWHM of each line, the linear part of the
  curve of growth, the lines combined by inverse variance, the limit
  ``sigma`` times the error of the column. It ignores what absorbs there:
  where the pixels absorb it is a sensitivity, not a limit (the recipe
  warns).

The pixels of a line are those within +-(3 FWHM + 5 b) of it, the window
of the evaluation of the identification (``fitting/evaluation_window.py``,
decided with Guido on 18 November 2026): the core of the line and the
continuum on both sides. For ``profile`` the pixels that the fit of the
refitted components takes from the data are added, computed once.

A limit describes the data and the model on its pixels when it was computed;
:func:`drop_stale` removes it when they change (Guido, 27 November 2026):
the session calls it whenever its spectrum or its system list changes.
"""
from __future__ import annotations

import dataclasses
import hashlib
import logging
from datetime import datetime, timezone
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import astropy.units as au
import numpy as np
from scipy.optimize import brentq, minimize_scalar

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.spectrum import Spectrum
from astrocook.core.system_list import SystemList
from astrocook.core.upper_limits import UpperLimit, limits_of, with_limits
from astrocook.fitting.component_flags import FWHM_OVER_B, b_instrumental
from astrocook.fitting.evaluation_window import WINDOW_B, WINDOW_FWHM
from astrocook.fitting.fit_record import current_recipe
from astrocook.fitting.lsf import C_KMS, FWHM_TO_SIGMA, KERNEL_SIGMAS
from astrocook.fitting.voigt_fitter import VoigtFitter, _transitions, chi2_of, line_reach

logger = logging.getLogger(__name__)

#: The linear part of the curve of growth: N = COG_LINEAR W_r / (f lambda^2),
#: N in cm^-2, W_r and lambda in Angstrom (m_e c^2 / (pi e^2) = 1.13e20)
COG_LINEAR = 1.13e20

#: The step of the scan of log N (dex). The chi2 grows about as N^2 in the
#: linear part of the curve of growth, so from the start (below the limit)
#: a few steps bracket it; the limit is then found to ``LOGN_XTOL``
LOGN_STEP = 0.25
#: The precision of the limit and of the minimum (dex)
LOGN_XTOL = 0.005
#: The scan starts this far below the estimate from the equivalent width
#: (the reference of no absorption is ``LOGN_FLOOR`` below it): at 2 dex
#: below the noise the line changes the chi2 by 1e-4 of sigma^2
LOGN_START = 1.0
LOGN_FLOOR = 2.0
#: The scan gives up above this: no limit below a damped column
LOGN_MAX = 22.0

#: The evaluations of the model in each refit of the components on the
#: pixels (as ``refit_all``)
MAX_NFEV = 200
#: The window of the redshift of the refitted components (km/s; as
#: ``fit_component``)
Z_WINDOW_KMS = 20.0


@dataclasses.dataclass(frozen=True)
class _Line:
    name: str
    wave: float
    """ Observed, Angstrom. """
    f: float
    rest: float
    """ Rest, Angstrom. """


def _lines_in(spec: Spectrum, series: str, z: float) -> Tuple[List[_Line], List[str]]:
    """The lines of ``series`` at ``z`` with a pixel in the spectrum, and the
    names of the others."""
    x = spec.x.to_value(au.Angstrom)
    lines, skipped = [], []
    for t in _transitions(series):
        at = ATOM_DATA[t]
        w = at['wave'] * (1.0 + z)
        if x[0] <= w <= x[-1]:
            lines.append(_Line(t, float(w), float(at['f']), float(at['wave'])))
        else:
            skipped.append(t)
    return lines, skipped


def _resol_at(spec: Spectrum, wave_ang: float) -> float:
    x = spec.x.to_value(au.Angstrom)
    r = spec.resolution()
    i = int(np.clip(np.searchsorted(x, wave_ang), 0, len(x) - 1))
    return float(r[i])


def _window(x_ang: np.ndarray, centre: float, half_kms: float) -> np.ndarray:
    lo, hi = centre * np.exp(-half_kms / C_KMS), centre * np.exp(half_kms / C_KMS)
    return (x_ang >= lo) & (x_ang <= hi)


def _ranges(x_ang: np.ndarray, mask: np.ndarray) -> Tuple[Tuple[float, float], ...]:
    """The runs of ``mask`` as ranges of the centres of their pixels (nm)."""
    idx = np.flatnonzero(mask)
    if not idx.size:
        return ()
    breaks = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate([[idx[0]], idx[breaks + 1]])
    ends = np.concatenate([idx[breaks], [idx[-1]]])
    return tuple((float(x_ang[a]) / 10.0, float(x_ang[b]) / 10.0)
                 for a, b in zip(starts, ends, strict=True))


def _in_ranges(x_ang: np.ndarray, ranges_nm: Sequence[Tuple[float, float]]) -> np.ndarray:
    """The pixels whose centres are in ``ranges_nm`` (a range given by the
    centres of its first and last pixels holds them, to rounding)."""
    mask = np.zeros(len(x_ang), dtype=bool)
    for lo, hi in ranges_nm:
        mask |= (x_ang >= lo * 10.0 * (1 - 1e-12)) & (x_ang <= hi * 10.0 * (1 + 1e-12))
    return mask


def _widened(ranges_nm: Sequence[Tuple[float, float]], pad_kms: float
             ) -> Tuple[Tuple[float, float], ...]:
    f = np.exp(pad_kms / C_KMS)
    return tuple((lo / f, hi * f) for lo, hi in ranges_nm)


def components_reaching(systs: SystemList, reach_nm: Sequence[Tuple[float, float]]) -> List[str]:
    """The uuids of the components with a line whose reach (its optical
    depth above the cut of the model, :func:`line_reach`) falls in
    ``reach_nm``: those that change the model on the pixels."""
    rows = [(c.uuid, ATOM_DATA[t], c) for c in systs.components for t in _transitions(c.series)]
    if not rows or not reach_nm:
        return []
    lam_c, half = line_reach(
        np.array([a['wave'] for _, a, _ in rows]), np.array([a['f'] for _, a, _ in rows]),
        np.array([a['gamma'] for _, a, _ in rows]), np.array([c.z for _, _, c in rows]),
        np.array([10.0 ** c.logN for _, _, c in rows]),
        np.array([np.hypot(c.b, c.btur or 0.0) for _, _, c in rows]))
    lo, hi = (lam_c - half) / 10.0, (lam_c + half) / 10.0
    hit = np.zeros(len(rows), dtype=bool)
    for a, b in reach_nm:
        hit |= (hi >= a) & (lo <= b)
    seen: Dict[str, None] = {}
    for (u, _, _), h in zip(rows, hit, strict=True):
        if h:
            seen[u] = None
    return list(seen)


def components_fingerprint(systs: SystemList, uuids: Sequence[str]) -> str:
    """A fingerprint of the parameters and constraints of the components
    ``uuids``: it changes when any of them changes, or one is added or
    removed."""
    by_uuid = {c.uuid: c for c in systs.components}
    cons = systs.parsed_constraints or {}
    parts = []
    for u in sorted(uuids):
        c = by_uuid.get(u)
        if c is None:
            continue
        links = []
        for p in ("z", "logN", "b", "btur"):
            k = cons.get((u, p))
            # No entry is the default: free, btur frozen (an archive read
            # back writes the defaults out)
            links.append((p != 'btur', None, None) if k is None else
                         (bool(getattr(k, 'is_free', True)), getattr(k, 'expression', None),
                          getattr(k, 'target_uuid', None)))
        parts.append(repr((u, c.series, float(c.z), float(c.logN), float(c.b),
                           float(c.btur or 0.0), c.resol, links)))
    return hashlib.sha1("\n".join(parts).encode()).hexdigest()


def data_fingerprint(spec: Spectrum, ranges_nm: Sequence[Tuple[float, float]]) -> str:
    """A fingerprint of the data on the pixels of ``ranges_nm``: wavelength,
    flux, error, normalization (continuum and tellurics) and resolution."""
    x = spec.x.to_value(au.Angstrom)
    sel = _in_ranges(x, ranges_nm)
    norm = spec.norm
    h = hashlib.sha1()
    for arr in (x, spec.y.value, spec.dy.value,
                None if norm is None else norm.value, spec.resolution()):
        if arr is None:
            h.update(b"none")
            continue
        h.update(np.ascontiguousarray(np.asarray(arr, dtype=float)[sel]).tobytes())
    return h.hexdigest()


def _valid(spec: Spectrum) -> np.ndarray:
    """The pixels with a finite normalized flux and a positive error (those
    of :attr:`VoigtFitter.valid_mask`)."""
    norm_q = spec.norm
    if norm_q is None:
        raise ValueError("the spectrum has no continuum")
    norm = np.asarray(norm_q.value, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        y = spec.y.value / norm
        dy = spec.dy.value / norm
    return np.isfinite(y) & np.isfinite(dy) & (dy > 0) & (norm > 0)


def _ew_estimate(spec: Spectrum, lines: Sequence[_Line], z: float, b: float,
                 b_inst: float) -> Tuple[float, float, np.ndarray, List[str]]:
    """The column measured from the equivalent width, its error (cm^-2),
    the pixels used and the lines with valid pixels: for each line the
    pixels within one FWHM (of b and b_inst in quadrature) of its centre, the
    lines combined by inverse variance."""
    x = spec.x.to_value(au.Angstrom)
    dl = spec.xmax.to_value(au.Angstrom) - spec.xmin.to_value(au.Angstrom)
    norm_q = spec.norm
    if norm_q is None:
        raise ValueError("the spectrum has no continuum")
    norm = np.asarray(norm_q.value, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        y = spec.y.value / norm
        dy = spec.dy.value / norm
    valid = np.isfinite(y) & np.isfinite(dy) & (dy > 0)
    fwhm = FWHM_OVER_B * float(np.hypot(b, b_inst))
    used = np.zeros(len(x), dtype=bool)
    weights, values, names = [], [], []
    for line in lines:
        sel = _window(x, line.wave, fwhm) & valid
        if not np.any(sel):
            continue
        coef = COG_LINEAR / (line.f * line.rest ** 2) / (1.0 + z)
        w = float(np.sum((1.0 - y[sel]) * dl[sel])) * coef
        s = float(np.sqrt(np.sum((dy[sel] * dl[sel]) ** 2))) * coef
        if not (np.isfinite(s) and s > 0):
            continue
        used |= sel
        weights.append(1.0 / s ** 2)
        values.append(w)
        names.append(line.name)
    if not weights:
        return float('nan'), float('nan'), used, names
    wsum = float(np.sum(weights))
    n_meas = float(np.dot(weights, values) / wsum)
    return n_meas, float(1.0 / np.sqrt(wsum)), used, names


def upper_limit(spec: Spectrum, systs: SystemList, series: str, z: float,
                b: Optional[float] = None, sigma: float = 3.0, method: str = "profile",
                ref_uuid: Optional[str] = None,
                check_stop: Optional[Callable[[], bool]] = None) -> UpperLimit:
    """The upper limit on the column of ``series`` at ``z`` (module
    docstring), with b fixed (None: the instrumental b at the strongest
    line).

    Raises
    ------
    ValueError
        If no line of the series falls in the spectrum with valid pixels,
        if b is not given and the resolution is not known, or if the scan
        finds no limit below ``LOGN_MAX``.
    """
    if method not in ("profile", "ew"):
        raise ValueError(f"unknown method {method!r}")
    lines, skipped = _lines_in(spec, series, z)
    if not lines:
        raise ValueError(f"no line of {series} at z = {z:.5f} falls in the spectrum")
    strongest = max(lines, key=lambda ln: ln.f * ln.rest)
    b_inst = b_instrumental(_resol_at(spec, strongest.wave))
    if b is None:
        if b_inst <= 0:
            raise ValueError(f"the resolution is not known at {strongest.name}: give b")
        b = b_inst

    n_meas, s_n, ew_pixels, ew_lines = _ew_estimate(spec, lines, z, b, b_inst)
    if not ew_lines:
        raise ValueError(f"no line of {series} at z = {z:.5f} has valid pixels")
    skipped += [ln.name for ln in lines if ln.name not in ew_lines]
    lines = [ln for ln in lines if ln.name in ew_lines]
    logN_ew = float(np.log10(sigma * s_n))
    x = spec.x.to_value(au.Angstrom)
    resol = spec.resolution()
    finite_r = resol[np.isfinite(resol) & (resol > 0)]
    pad_kms = (KERNEL_SIGMAS * C_KMS / float(np.min(finite_r)) * FWHM_TO_SIGMA
               if finite_r.size else 0.0)
    # What both methods record; the numbers are set below
    base = UpperLimit(series=series, z=float(z), b=float(b), logN_lim=float('nan'),
                      sigma=float(sigma), method=method,
                      transitions=tuple(ln.name for ln in lines), skipped=tuple(skipped),
                      b_inst=float(b_inst), ref_uuid=ref_uuid, recipe=current_recipe(),
                      time=datetime.now(timezone.utc).isoformat(timespec='seconds'))

    if method == "ew":
        ranges = _ranges(x, ew_pixels)
        return dataclasses.replace(
            base, logN_lim=logN_ew, ranges_nm=ranges, n_pix=int(np.sum(ew_pixels)),
            logN_best=float(np.log10(n_meas)) if n_meas > 0 else None,
            significance=float(max(n_meas, 0.0) / s_n), data=data_fingerprint(spec, ranges))

    # --- profile ---
    fitter0 = VoigtFitter(spec, systs) if systs.components else None
    valid = _valid(spec)
    window = np.zeros(len(x), dtype=bool)
    for ln in lines:
        r = _resol_at(spec, ln.wave)
        fwhm_kms = C_KMS / r if np.isfinite(r) and r > 0 else 0.0
        window |= _window(x, ln.wave, WINDOW_FWHM * fwhm_kms + WINDOW_B * b)
    window &= valid
    reach = _widened(_ranges(x, window), pad_kms)
    refitted = components_reaching(systs, reach)
    pixels = window.copy()
    if refitted and fitter0 is not None:
        comps = [c for c in systs.components if c.uuid in set(refitted)]
        pixels |= fitter0.feature_mask(comps, Z_WINDOW_KMS)
        pixels &= valid
    ranges = _ranges(x, pixels)
    # The components on the pixels now, with those of the fit added
    reach = _widened(ranges, pad_kms)
    refitted = components_reaching(systs, reach)

    cache: Dict[float, float] = {}

    def chi2_at(logN: float) -> float:
        key = round(float(logN), 6)
        if key in cache:
            return cache[key]
        if check_stop is not None and check_stop():
            raise InterruptedError("the computation of the upper limit was cancelled")
        trial = systs.add_component(series, float(z), logN=float(logN), b=float(b), btur=0.0)
        new_uuid = trial.components[-1].uuid
        for p in ("z", "logN", "b"):
            trial = trial.update_constraint(new_uuid, p, is_free=False)
        if refitted:
            with trial.fitting_context([new_uuid] + list(refitted), group_depth=0, spec=spec):
                fitter = VoigtFitter(spec, trial)
                trial, _, _ = fitter.fit(max_nfev=MAX_NFEV, z_window_kms=Z_WINDOW_KMS,
                                         model=False)
            if trial.constraint_model:
                trial.constraint_model.set_active_components(None)
        ev = VoigtFitter(spec, trial)
        _, model = ev.compute_model_flux(pixels=pixels)
        value = chi2_of(ev.normalized_residuals(model, pixels))
        cache[key] = value
        return value

    target = sigma ** 2
    start = logN_ew - LOGN_START
    xs = [logN_ew - LOGN_FLOOR]
    chi = [chi2_at(xs[0])]
    xv = start
    while True:
        xs.append(xv)
        chi.append(chi2_at(xv))
        if chi[-1] > min(chi) + target:
            break
        xv += LOGN_STEP
        if xv > LOGN_MAX:
            raise ValueError(f"the chi2 does not rise by {target:g} below log N {LOGN_MAX:g}: "
                             "no limit")
    i = int(np.argmin(chi))
    chi_min = chi[i]
    logN_best: Optional[float] = None
    if i > 0:
        res = minimize_scalar(chi2_at, bounds=(xs[i - 1], xs[i + 1]), method='bounded',
                              options={'xatol': LOGN_XTOL})
        if float(res.fun) < chi_min:
            chi_min, logN_best = float(res.fun), float(res.x)
        else:
            logN_best = xs[i]
    # The bracket: from the last point of the scan below the threshold, at
    # or above the minimum (the minimum itself if none), to the first above
    lower = logN_best if logN_best is not None else xs[0]
    for xk, ck in zip(xs[i:-1], chi[i:-1], strict=True):
        if ck <= chi_min + target and xk > lower:
            lower = xk
    lim = float(brentq(lambda v: chi2_at(v) - chi_min - target, lower, xs[-1],
                       xtol=LOGN_XTOL))
    return dataclasses.replace(
        base, logN_lim=lim, ranges_nm=ranges, n_pix=int(np.sum(pixels)), logN_best=logN_best,
        significance=float(np.sqrt(max(chi[0] - chi_min, 0.0))), refitted=tuple(refitted),
        reach_nm=reach, components=components_fingerprint(systs, refitted),
        data=data_fingerprint(spec, ranges))


def stale_reason(spec: Optional[Spectrum], systs: SystemList, lim: UpperLimit) -> Optional[str]:
    """Why ``lim`` no longer describes the data and the model, or None."""
    if lim.ref_uuid is not None and systs.get_component_by_uuid(lim.ref_uuid) is None:
        return "the component of its redshift was deleted"
    if spec is None or data_fingerprint(spec, lim.ranges_nm) != lim.data:
        return "the data on its pixels changed"
    if lim.method == "profile":
        now = components_reaching(systs, lim.reach_nm)
        if components_fingerprint(systs, now) != lim.components:
            return "the components on its pixels changed"
    return None


def drop_stale(spec: Optional[Spectrum], systs: SystemList
               ) -> Tuple[SystemList, List[Tuple[UpperLimit, str]]]:
    """The system list without the limits that no longer describe the data
    and the model (:func:`stale_reason`), and those removed with the
    reason."""
    limits = limits_of(systs)
    if not limits:
        return systs, []
    kept, removed = [], []
    for lim in limits:
        why = stale_reason(spec, systs, lim)
        if why is None:
            kept.append(lim)
        else:
            removed.append((lim, why))
    if not removed:
        return systs, []
    return with_limits(systs, tuple(kept)), removed
