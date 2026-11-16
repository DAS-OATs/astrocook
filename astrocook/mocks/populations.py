"""Random populations of absorbers for mocks.

The default parameters are standard functional forms with indicative
values, to be calibrated against the literature before scientific use:

- Lyman-alpha forest: column densities from f(N) ~ N^-beta (beta = 1.5),
  Doppler parameters log-normal around 25 km/s, line density
  dN/dz = A ((1 + z) / 4)^gamma. The default A = 1700 (12 < logN < 17,
  gamma = 2.5) gives a mean transmission within ~0.03 of
  exp(-tau_eff), tau_eff = 0.0018 (1 + z)^3.92 (Faucher-Giguere et al.
  2008), at z = 3 and 3.8.
- Metal systems: a few components within a small velocity range,
  optionally with an H I component at the system redshift.
- Associated absorbers: metal systems within a few thousand km/s of
  z_em, beyond it too; and a proximity zone of the forest
  (:func:`lya_forest`, ``proximity_kms``), phase 7, PR 7a.1.
- Strong H I absorbers (Lyman-limit systems, sub-DLAs, DLAs): column
  densities from f(N) ~ N^-beta (beta = 1.7) within a range, uniform
  Doppler parameters, uniform redshifts, all the Lyman lines of the
  atomic data table (:data:`~astrocook.mocks.physics.LYMAN_SERIES`).
  Combine with ``MockTruth(lyman_limit=True)`` for the absorption below
  the Lyman limit.
"""
from typing import List, Optional, Sequence, Tuple

import numpy as np

from astrocook.mocks.physics import C_KMS, LYMAN_SERIES
from astrocook.mocks.truth import MockComponent

#: Column density ranges (log10 cm^-2) of the strong H I absorbers
LLS_LOGN = (17.2, 19.0)
SUBDLA_LOGN = (19.0, 20.3)
DLA_LOGN = (20.3, 22.0)


def _power_law(rng, n, x_min, x_max, slope):
    """Samples from p(x) ~ x^-slope in [x_min, x_max] (slope != 1)."""
    e = 1.0 - slope
    u = rng.random(n)
    return (x_min**e + u * (x_max**e - x_min**e)) ** (1.0 / e)


def _redshifts(rng, n, z_min, z_max, gamma):
    """Samples from p(z) ~ (1 + z)^gamma in [z_min, z_max]."""
    return _power_law(rng, n, 1.0 + z_min, 1.0 + z_max, -gamma) - 1.0


def expected_forest_lines(z_min: float, z_max: float, dndz_at_3: float = 1700.0,
                          gamma: float = 2.5) -> float:
    """Expected number of forest lines between z_min and z_max."""
    a = dndz_at_3 / 4.0**gamma
    return a / (gamma + 1.0) * ((1 + z_max) ** (gamma + 1) - (1 + z_min) ** (gamma + 1))


def mean_transmission_fg08(z: float) -> float:
    """Mean Lyman-alpha forest transmission, exp(-0.0018 (1 + z)^3.92)."""
    return float(np.exp(-0.0018 * (1.0 + z) ** 3.92))


def lya_forest(rng: np.random.Generator, z_min: float, z_max: float,
               logN_range: Tuple[float, float] = (12.0, 17.0), beta: float = 1.5,
               b_median: float = 25.0, b_sigma_dex: float = 0.15,
               b_min: float = 10.0, dndz_at_3: float = 1700.0,
               gamma: float = 2.5, series: str = "Ly_a",
               strong_series: Optional[str] = None,
               strong_logN: float = 15.0, z_em: Optional[float] = None,
               proximity_kms: Optional[float] = None) -> List[MockComponent]:
    """Random Lyman-alpha forest between ``z_min`` and ``z_max``.

    The number of lines is Poisson-distributed around
    :func:`expected_forest_lines`; ``dndz_at_3`` is the line density at
    z = 3 for column densities in ``logN_range``. ``series`` is the series
    of every line: ``'Ly_a'`` (the default), ``'Ly_ab'``, or
    :data:`~astrocook.mocks.physics.LYMAN_SERIES` for the whole series.
    With ``strong_series`` the lines with logN >= ``strong_logN`` get that
    series instead, e.g. the whole Lyman series only where the higher
    lines matter (a few percent of the lines with the defaults, so the
    mock stays fast); the draws are the same.

    With ``z_em`` and ``proximity_kms`` the forest has a proximity zone
    (phase 7, PR 7a.1): the column of each line is divided by 1 + omega,
    the ratio of the ionizing flux of the quasar to the background, with
    omega = (``proximity_kms`` / dv)^2 at the velocity dv from ``z_em``
    (the 1/r^2 of the flux, with the distance proportional to dv near the
    quasar), so that ``proximity_kms`` is where the two are equal. With
    f(N) ~ N^-beta this reduces the density of the lines above any column
    by (1 + omega)^(1 - beta). The draws are the same as without it, and
    the series of a line follows its reduced column.
    """
    n = rng.poisson(expected_forest_lines(z_min, z_max, dndz_at_3, gamma))
    z = _redshifts(rng, n, z_min, z_max, gamma)
    logN = np.log10(_power_law(rng, n, 10**logN_range[0], 10**logN_range[1], beta))
    b = np.maximum(b_median * 10 ** (b_sigma_dex * rng.standard_normal(n)), b_min)
    if proximity_kms is not None:
        if z_em is None:
            raise ValueError("proximity_kms needs z_em")
        logN = logN - np.log10(1.0 + proximity_omega(z, z_em, proximity_kms))
    order = np.argsort(z)
    def series_of(n):
        if strong_series is not None and n >= strong_logN:
            return strong_series
        return series
    return [MockComponent(series_of(logN[i]), float(z[i]), float(logN[i]),
                          float(b[i])) for i in order]


def proximity_omega(z, z_em: float, proximity_kms: float) -> np.ndarray:
    """Ratio of the ionizing flux of the quasar to the background at
    redshift ``z`` (see :func:`lya_forest`); infinite at and beyond
    ``z_em``."""
    dv = (z_em - np.asarray(z, dtype=float)) / (1.0 + z_em) * C_KMS
    with np.errstate(divide="ignore"):
        return np.where(dv > 0, (proximity_kms / np.where(dv > 0, dv, 1.0)) ** 2,
                        np.inf)


def metal_systems(rng: np.random.Generator, n_systems: int, z_min: float,
                  z_max: float, species: Sequence[str] = ("CIV", "SiIV", "MgII"),
                  max_components: int = 4, dv_kms: float = 50.0,
                  logN_range: Tuple[float, float] = (12.0, 14.5),
                  b_range: Tuple[float, float] = (4.0, 15.0),
                  hi_logN_range: Optional[Tuple[float, float]] = None,
                  hi_b_range: Tuple[float, float] = (15.0, 35.0)
                  ) -> List[MockComponent]:
    """Random metal systems, each with 1 to ``max_components`` components
    of one species, spread within +/- ``dv_kms`` of the system redshift.

    With ``hi_logN_range`` each system also has one H I component (all the
    Lyman lines) at the system redshift, with logN uniform in that range
    and b uniform in ``hi_b_range``; without it the systems are the same
    as before the option existed (same random draws).
    """
    comps = []
    for _ in range(n_systems):
        series = str(rng.choice(list(species)))
        z_sys = rng.uniform(z_min, z_max)
        for _ in range(rng.integers(1, max_components + 1)):
            dv = rng.uniform(-dv_kms, dv_kms)
            comps.append(MockComponent(
                series, float(z_sys + dv / C_KMS * (1 + z_sys)),
                float(rng.uniform(*logN_range)), float(rng.uniform(*b_range))))
        if hi_logN_range is not None:
            comps.append(MockComponent(
                LYMAN_SERIES, float(z_sys), float(rng.uniform(*hi_logN_range)),
                float(rng.uniform(*hi_b_range))))
    return sorted(comps, key=lambda c: c.z)


def strong_hi_systems(rng: np.random.Generator, n_systems: int, z_min: float,
                      z_max: float, logN_range: Tuple[float, float] = DLA_LOGN,
                      beta: float = 1.7,
                      b_range: Tuple[float, float] = (15.0, 35.0)
                      ) -> List[MockComponent]:
    """Random strong H I absorbers, one component each, with all the
    Lyman lines: Lyman-limit systems, sub-DLAs or DLAs depending on
    ``logN_range`` (:data:`LLS_LOGN`, :data:`SUBDLA_LOGN`,
    :data:`DLA_LOGN`)."""
    z = rng.uniform(z_min, z_max, n_systems)
    logN = np.log10(_power_law(rng, n_systems, 10**logN_range[0],
                               10**logN_range[1], beta))
    b = rng.uniform(*b_range, n_systems)
    order = np.argsort(z)
    return [MockComponent(LYMAN_SERIES, float(z[i]), float(logN[i]), float(b[i]))
            for i in order]


def associated_absorbers(rng: np.random.Generator, n_systems: int, z_em: float,
                         v_range_kms: Tuple[float, float] = (-2000.0, 3000.0),
                         species: Sequence[str] = ("CIV", "NV", "SiIV"),
                         max_components: int = 3, dv_kms: float = 50.0,
                         logN_range: Tuple[float, float] = (13.0, 15.0),
                         b_range: Tuple[float, float] = (4.0, 15.0),
                         hi_logN_range: Optional[Tuple[float, float]] = (14.0, 17.0)
                         ) -> List[MockComponent]:
    """Random associated absorbers (phase 7, PR 7a.1): metal systems (each
    with its H I) at velocities from ``z_em`` in ``v_range_kms``, positive
    beyond z_em (z > z_em, as for infall or a z_em underestimated), i.e.
    redshifts uniform in z_em + v / c (1 + z_em). The margin of the
    identification beyond z_em is 3000 km/s (F55, ``z_em_margin_kms``): a
    range beyond it gives systems that the identification should reject.
    The systems are those of :func:`metal_systems`."""
    z_lo, z_hi = (z_em + v / C_KMS * (1.0 + z_em) for v in v_range_kms)
    return metal_systems(rng, n_systems, z_lo, z_hi, species=species,
                         max_components=max_components, dv_kms=dv_kms,
                         logN_range=logN_range, b_range=b_range,
                         hi_logN_range=hi_logN_range)
