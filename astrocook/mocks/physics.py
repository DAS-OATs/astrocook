"""Voigt optical depth and instrumental convolution for mocks.

Written independently of :mod:`astrocook.fitting.voigt_fitter`, on purpose.
"""
from typing import List, Union

import numpy as np
from astropy import constants as const
from astropy import units as au
from scipy.special import wofz

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS

C_KMS = const.c.to_value(au.km / au.s)
# pi e^2 / (m_e c) in cm^2 s^-1 (Gaussian units)
_PI_E2_MEC = (np.pi * const.e.esu**2 / (const.m_e * const.c)).cgs.value
FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

#: The H I Lyman series of the atomic data table, Ly_a to Ly_15 (the
#: ``Lyman`` multiplet of ``STANDARD_MULTIPLETS`` stops at Ly_8), as a
#: series string for :class:`~astrocook.mocks.MockComponent`
LYMAN_SERIES = ",".join(k for k in ATOM_DATA
                        if k.startswith("Ly_") and ATOM_DATA[k]["f"] > 0)
#: H I ionization edge (Angstrom), as in the atomic data table
LYMAN_LIMIT_ANG = ATOM_DATA["Ly_lim"]["wave"]
#: H I photoionization cross section at the edge (cm^2)
SIGMA_LYMAN_LIMIT = 6.30e-18


def transitions(series: str) -> List[str]:
    """Transitions of a series (single line, standard multiplet or list)."""
    if series in STANDARD_MULTIPLETS:
        return list(STANDARD_MULTIPLETS[series])
    names = [s.strip() for s in series.split(",") if s.strip()]
    missing = [n for n in names if n not in ATOM_DATA]
    if missing:
        raise ValueError(f"Unknown transitions: {missing}")
    return names


def optical_depth(lam_ang: np.ndarray, series: str, z: float, logN: float,
                  b: float, btur: float = 0.0) -> np.ndarray:
    """Optical depth of one component at observed wavelengths (Angstrom).

    tau = (pi e^2 / m_e c) f lambda_0 N / (sqrt(pi) b) H(a, u), with
    u = (lambda_rest - lambda_0) / Delta lambda_D,
    Delta lambda_D = (b / c) lambda_0, a = Gamma lambda_0 / (4 pi b) and
    H(a, u) = Re[w(u + i a)].
    """
    b_kms = np.sqrt(b**2 + btur**2)
    b_cms = b_kms * 1e5
    lam_rest = np.asarray(lam_ang, dtype=float) / (1.0 + z)
    n_col = 10.0**logN
    tau = np.zeros_like(lam_rest)
    for name in transitions(series):
        atom = ATOM_DATA[name]
        lam0_ang = atom["wave"]
        lam0_cm = lam0_ang * 1e-8
        u = (lam_rest - lam0_ang) / (b_kms / C_KMS * lam0_ang)
        a = atom["gamma"] * lam0_cm / (4.0 * np.pi * b_cms)
        h = wofz(u + 1j * a).real
        tau += _PI_E2_MEC * atom["f"] * lam0_cm * n_col / (np.sqrt(np.pi) * b_cms) * h
    return tau


def is_hydrogen(series: str) -> bool:
    """True if the series has H I Lyman lines."""
    return any(name.startswith("Ly_") for name in transitions(series))


def lyman_limit_optical_depth(lam_ang: np.ndarray, z: float,
                              logN: float) -> np.ndarray:
    """Optical depth of the H I continuum below the Lyman limit.

    tau = N sigma_LL (lambda_rest / lambda_LL)^3 for lambda_rest <=
    lambda_LL, 0 above. The (nu_LL / nu)^3 law is the usual approximation
    of the hydrogenic cross section (the exact one falls slightly more
    slowly, ~nu^-2.75, near the edge). The Lyman lines above Ly_15 are
    not modelled: between Ly_15 (915.3 A) and the edge they would merge
    into the edge of a strong absorber.
    """
    lam_rest = np.asarray(lam_ang, dtype=float) / (1.0 + z)
    ratio = lam_rest / LYMAN_LIMIT_ANG
    return np.where(ratio <= 1.0, 10.0**logN * SIGMA_LYMAN_LIMIT * ratio**3, 0.0)


def convolve_velocity(flux: np.ndarray, dv_kms: float,
                      resol: Union[float, np.ndarray, None]) -> np.ndarray:
    """Gaussian instrumental convolution on a grid uniform in velocity.

    The kernel at each pixel has FWHM = c / R(pixel), so a varying
    resolution is handled exactly (no grouping of pixels). The flux is
    extended with its edge values to avoid edge effects.
    """
    if resol is None:
        return np.array(flux, dtype=float)
    flux = np.asarray(flux, dtype=float)
    n = len(flux)
    sigma_pix = np.broadcast_to(
        C_KMS / np.asarray(resol, dtype=float) * FWHM_TO_SIGMA / dv_kms, (n,))
    half = int(np.ceil(5.0 * sigma_pix.max()))
    padded = np.concatenate([np.full(half, flux[0]), flux, np.full(half, flux[-1])])
    num = np.zeros(n)
    den = np.zeros(n)
    for k in range(-half, half + 1):
        w = np.exp(-0.5 * (k / sigma_pix) ** 2)
        num += w * padded[half + k: half + k + n]
        den += w
    return num / den
