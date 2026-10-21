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
               gamma: float = 2.5, series: str = "Ly_a") -> List[MockComponent]:
    """Random Lyman-alpha forest between ``z_min`` and ``z_max``.

    The number of lines is Poisson-distributed around
    :func:`expected_forest_lines`; ``dndz_at_3`` is the line density at
    z = 3 for column densities in ``logN_range``. ``series`` is the series
    of every line: ``'Ly_a'`` (the default), ``'Ly_ab'``, or
    :data:`~astrocook.mocks.physics.LYMAN_SERIES` for the whole series.
    """
    n = rng.poisson(expected_forest_lines(z_min, z_max, dndz_at_3, gamma))
    z = _redshifts(rng, n, z_min, z_max, gamma)
    logN = np.log10(_power_law(rng, n, 10**logN_range[0], 10**logN_range[1], beta))
    b = np.maximum(b_median * 10 ** (b_sigma_dex * rng.standard_normal(n)), b_min)
    order = np.argsort(z)
    return [MockComponent(series, float(z[i]), float(logN[i]), float(b[i]))
            for i in order]


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
