"""The flags of the components whose values are not measured (phase 7, PR
7d.2; findings register F52, F116).

The fit leaves on each component it fits the flags of #517 (a parameter at
a bound, ``nan_errors``: ``voigt_fitter.COMPONENT_FLAGS``); this module adds
those that say a value is there but is not a measurement:

- ``unresolved``: b is not measured, the line is narrower than the
  resolution;
- ``saturated``: log N is not measured, the weakest line of the series in
  the fit is saturated;
- ``large_errors``: the errors have no meaning (F52);
- ``telluric_overlap``: a line of the series falls on a telluric
  absorption, so its errors are not to be trusted.

They describe the last fit, as those of #517, and are annotations only: no
number of the fit and no decision of a pipeline depends on them (design:
``ROADMAP.md``, phase 7, "The flags of 7d.2", agreed with Guido on 26
November 2026). :func:`measurement_flags` is a pure function of the
component, the atomic data, the resolution and the telluric column, called
by the fitter on the components it fits.

The thresholds are heuristics; what supports each is written next to it,
and the measurement on the mocks is in the PR of 7d.2.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np
from scipy.special import wofz

from astrocook.core.atomic_data import ATOM_DATA, series_lines
from astrocook.core.structures import ComponentData

C_KMS = 299792.458

#: The flags of this module, in the order they are added to
#: ``ComponentData.flags`` (after those of #517)
MEASUREMENT_FLAGS = ('unresolved', 'saturated', 'large_errors', 'telluric_overlap')

#: FWHM / b of a Gaussian profile, 2 sqrt(ln 2): the b of the instrumental
#: profile of resolving power R is c / (1.665 R) (4.0 km/s at R = 45000)
FWHM_OVER_B = 2.0 * np.sqrt(np.log(2.0))

#: ``unresolved``: b below this fraction of the instrumental b. A line with
#: b < b_inst / 2 has an observed width within 12 per cent of the
#: instrumental one (sqrt(1 + 1/4) = 1.12), so its b is set by a few per
#: cent of the observed width, below what a fit at the usual S/N resolves
#: (F116: b = 1.4-2.8 km/s at R = 45000, b_inst = 4.0 km/s, judged
#: unreliable by both agents of the trials of PR 7 of phase 6-ter)
UNRESOLVED_RATIO = 0.5

#: ``saturated``: the central optical depth of the weakest line of the
#: series above this. tau_0 = 3 is a transmission of 5 per cent at the
#: centre of the intrinsic profile: from there the core is flat and log N
#: moves the equivalent width only through the wings (the flat part of the
#: curve of growth)
SATURATED_TAU = 3.0

#: A damped H I (log N at or above this) is not flagged ``saturated``: its
#: Lorentzian wings constrain log N
DAMPED_LOGN = 19.0

#: ``large_errors`` (F52): the error of log N above this (dex), the rule of
#: the benchmark since phase 4 (``_flags`` of ``validate_mocks.py``: a
#: factor of 10 in N is no measurement). The other two: db > b, and the
#: error of z in velocity above max(b, b_inst), the width of the line
LARGE_DLOGN = 1.0

#: ``telluric_overlap``: the telluric depth 1 - T above this within the
#: reach of a line. 5 per cent is about the noise of a pixel at S/N 20, the
#: lowest of the real spectra of the benchmark: a shallower telluric line
#: is below the noise of every spectrum of the sample (the threshold of the
#: plan of 7d.2; its counts on the mocks with tellurics are in the PR)
TELLURIC_DEPTH = 0.05

#: ``telluric_overlap``: the reach of a line, in units of max(b, b_inst)
#: around its centre: two Doppler widths hold the core of a Gaussian line
#: down to exp(-4), 2 per cent of its central depth
TELLURIC_REACH = 2.0

#: Prefactor of the central optical depth, pi e^2 / (m_e c) in the units of
#: the fitter (N in cm^-2, lambda in Angstrom, b in km/s): the same of
#: ``voigt_fitter._line_terms``
_TAU_PREFACTOR = 1.4974e-15


@dataclass(frozen=True)
class FitPixels:
    """What the flags read of the spectrum of a fit, on its pixels: the
    centres (Angstrom), the resolving power R of each (NaN where it is not
    known), the pixels of the fit (None: all) and the telluric
    transmission (None: no ``telluric`` column)."""
    wave: np.ndarray
    resol: np.ndarray
    fit_mask: Optional[np.ndarray] = None
    telluric: Optional[np.ndarray] = None


@dataclass(frozen=True)
class FlagQuantities:
    """The quantities of the criteria of a component, so that a flag can be
    checked (the panel of 7d.9 will show them): the instrumental b and R
    at the strongest line, the central optical depth of the weakest line,
    the largest telluric depth within the reach of its lines (None without
    the column)."""
    b_inst: float
    resol: float
    tau0_weakest: float
    telluric_depth: Optional[float]
    strongest: str = ""
    """The name of the strongest line (that of R and b_inst)."""
    weakest: str = ""
    """The name of the weakest line in the fit (that of tau_0)."""


def spectrum_pixels(spectrum) -> FitPixels:
    """What the flags read of a whole spectrum (:class:`FitPixels`, no
    pixels of a fit): for the quantities of the criteria shown by the
    Details panel of the System Inspector (7d.9), as the fitter takes them."""
    import astropy.units as au
    tell = spectrum.get_column('telluric')
    return FitPixels(wave=spectrum.x.to(au.Angstrom).value,
                     resol=np.asarray(spectrum.resolution(), dtype=float),
                     telluric=None if tell is None else np.asarray(tell.value, dtype=float))


def b_instrumental(resol: float) -> float:
    """The b (km/s) of a Gaussian instrumental profile of resolving power
    ``resol``; 0 if it is not known."""
    if resol is None or not np.isfinite(resol) or resol <= 0:
        return 0.0
    return float(C_KMS / (FWHM_OVER_B * resol))


def _lines(series: str) -> List[dict]:
    return [ATOM_DATA[t] for t in series_lines(series) if t in ATOM_DATA]


def central_tau(atom: dict, logN: float, b_kms: float) -> float:
    """The optical depth at the centre of the intrinsic Voigt profile of a
    line (before the convolution with the instrumental profile)."""
    b_safe = max(float(b_kms), 0.1)  # as the fitter
    a = atom['gamma'] * atom['wave'] * 1e-13 / (4.0 * np.pi * b_safe)
    tau0 = _TAU_PREFACTOR * 10.0 ** logN * atom['f'] * atom['wave'] / b_safe
    return float(tau0 * wofz(1j * a).real)


def _window(pix: FitPixels, centre: float, half_kms: float) -> np.ndarray:
    """The indices of the pixels within ``half_kms`` of ``centre``; the
    nearest one if none."""
    half = centre * half_kms / C_KMS
    lo, hi = np.searchsorted(pix.wave, [centre - half, centre + half])
    if hi > lo:
        return np.arange(lo, hi)
    i = int(np.clip(np.searchsorted(pix.wave, centre), 0, len(pix.wave) - 1))
    if i > 0 and abs(pix.wave[i - 1] - centre) < abs(pix.wave[i] - centre):
        i -= 1
    return np.array([i])


def _resol_at(comp: ComponentData, pix: FitPixels, centre: float) -> float:
    """R at a wavelength, as the fit takes it: the component's own, else the
    spectrum's at the nearest pixel."""
    if comp.resol is not None and comp.resol > 0:
        return float(comp.resol)
    if not len(pix.wave):
        return float('nan')
    return float(pix.resol[_window(pix, centre, 0.0)[0]])


def flag_quantities(comp: ComponentData, pix: FitPixels
                    ) -> Tuple[Optional[FlagQuantities], List[dict]]:
    """The quantities of the criteria of ``comp`` and the lines of its
    series in the fit (those with a pixel of the fit within their reach;
    all those in the spectrum if none). ``(None, [])`` if no line of the
    series is known or falls in the spectrum."""
    lines = _lines(comp.series)
    if not lines or not len(pix.wave):
        return None, []
    lo, hi = pix.wave[0], pix.wave[-1]
    lines = [at for at in lines if lo <= at['wave'] * (1 + comp.z) <= hi]
    if not lines:
        return None, []
    b_tot = float(np.hypot(comp.b, comp.btur or 0.0))

    # R and b_inst at the strongest line (the largest f lambda)
    strongest = max(lines, key=lambda at: at['f'] * at['wave'])
    resol = _resol_at(comp, pix, strongest['wave'] * (1 + comp.z))
    b_inst = b_instrumental(resol)
    reach = TELLURIC_REACH * max(b_tot, b_inst)

    if pix.fit_mask is not None:
        fitted = [at for at in lines
                  if np.any(pix.fit_mask[_window(pix, at['wave'] * (1 + comp.z), reach)])]
        lines = fitted or lines
    weakest = min(lines, key=lambda at: at['f'] * at['wave'])
    tau0 = central_tau(weakest, comp.logN, b_tot)
    names = {id(ATOM_DATA[t]): t for t in series_lines(comp.series) if t in ATOM_DATA}

    depth: Optional[float] = None
    if pix.telluric is not None:
        depths = [1.0 - pix.telluric[_window(pix, at['wave'] * (1 + comp.z), reach)]
                  for at in lines]
        finite = np.concatenate(depths)
        finite = finite[np.isfinite(finite)]
        depth = float(finite.max()) if finite.size else None
    return FlagQuantities(b_inst=b_inst, resol=resol, tau0_weakest=tau0,
                          telluric_depth=depth, strongest=names.get(id(strongest), ""),
                          weakest=names.get(id(weakest), "")), lines


def large_error_params(comp: ComponentData, b_inst: float = 0.0) -> Tuple[str, ...]:
    """The parameters whose errors have no meaning (F52): ``'logN'`` if dlogN
    > 1, ``'b'`` if db > b, ``'z'`` if the error of z in velocity exceeds
    max(b, b_inst). Empty unless the three errors are finite (otherwise the
    flag is ``nan_errors``)."""
    dz, dlogN, db = comp.dz, comp.dlogN, comp.db
    if dz is None or dlogN is None or db is None or not np.all(np.isfinite([dz, dlogN, db])):
        return ()
    out = []
    if dlogN > LARGE_DLOGN:
        out.append('logN')
    if db > comp.b:
        out.append('b')
    if C_KMS * dz / (1.0 + comp.z) > max(comp.b, b_inst):
        out.append('z')
    return tuple(out)


def _is_hydrogen(series: str) -> bool:
    lines = _lines(series)
    return bool(lines) and all(at.get('mass', 99.0) < 2.0 for at in lines)


def measurement_flags(comp: ComponentData, pix: FitPixels) -> Tuple[str, ...]:
    """The flags of :data:`MEASUREMENT_FLAGS` of a fitted component, in their
    order: ``unresolved`` (b < 0.5 b_inst, R at the strongest line of the
    series), ``saturated`` (the central optical depth of the weakest line
    of the series in the fit above 3, on the intrinsic profile; not for a
    damped H I), ``large_errors`` (:func:`large_error_params`),
    ``telluric_overlap`` (a line of the series in the fit where the
    ``telluric`` column has 1 - T > 0.05 within 2 max(b, b_inst) of its
    centre; never without the column)."""
    q, _ = flag_quantities(comp, pix)
    b_inst = q.b_inst if q is not None else 0.0
    flags = []
    if q is not None:
        b_tot = float(np.hypot(comp.b, comp.btur or 0.0))
        if b_inst > 0 and b_tot < UNRESOLVED_RATIO * b_inst:
            flags.append('unresolved')
        damped = _is_hydrogen(comp.series) and comp.logN >= DAMPED_LOGN
        if q.tau0_weakest > SATURATED_TAU and not damped:
            flags.append('saturated')
    if large_error_params(comp, b_inst):
        flags.append('large_errors')
    if q is not None and q.telluric_depth is not None and q.telluric_depth > TELLURIC_DEPTH:
        flags.append('telluric_overlap')
    return tuple(flags)


def with_measurement_flags(flags: Sequence[str], new: Sequence[str]) -> Tuple[str, ...]:
    """``flags`` with those of :data:`MEASUREMENT_FLAGS` replaced by ``new``
    (the others kept, in their order)."""
    return tuple(f for f in flags if f not in MEASUREMENT_FLAGS) + tuple(new)
