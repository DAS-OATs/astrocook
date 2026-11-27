"""The Laplace evidence of a hypothesis, with priors of population per
species (phase 7, PR 7c.4d-quater, an experiment: an option, no default
changed; ``design/identification.md`` section 3.6).

The information criteria (AIC, BIC, the SpIC of Webb et al. 2021) count the
parameters of a hypothesis, not how plausible their values are: a line of
the forest fitted as CII 1334 with b = 6 km/s, or as a Ly-a with b = 6 km/s,
has the same chi2 and the same k (F174). The evidence of a hypothesis H,

    Z = integral L(theta) pi(theta | H) d theta,

weighs the likelihood by the prior of its parameters: the expected number
of absorbers of the species (the incidence, per unit redshift), the
distribution of their columns and of their b. Its Laplace approximation, to
second order about the point of the short fit theta_0::

    ln Z = ln L_0 + ln pi_0 + (1/2) g^T H^-1 g + (k/2) ln 2 pi - (1/2) ln det H

with H = J^T J + B, J the Jacobian of the normalized residuals (the
analytic one of the fit, PR 7c.4d-bis: one Jacobian, so the evidence costs
about three models), B the curvature of -ln pi (its concave part: a convex
prior, as the log-uniform b, adds no curvature), g the gradient of ln pi at
theta_0 (the term is the step from the maximum of the likelihood to that of
the posterior). The likelihood is taken at its maximum, the end of the
short fit, and its own gradient left out: with it, a fit stopped before
convergence (300 evaluations) gained the step that the Gauss-Newton model
predicts, which overestimates it (a CIV fitted on two unrelated Ly-a lines,
b at 87 km/s, gained 2160 in ln Z, PR 7c.4d-quater). ln L_0 is
-chi2/2 on the pixels of the comparison (the constant of the errors is the
same for every hypothesis on the same pixels, and left out). The criterion
reported is ``evidence`` = -2 ln Z, in the units of the criteria (lower is
better; the weights exp(-Delta/2) of the comparison are then the posterior
probabilities of the hypotheses, the prior odds included).

The prior of a hypothesis: its components are independent absorbers of a
Poisson process, each with the density

    pi(z, log N, b) = rho_s(z) * phi_s(log N) * g_s(b),

rho_s the expected number of components of the species s per unit z above
a reference column (the incidence), phi_s the distribution of log N (a power
law f(N) ~ N^-beta, normalized above the reference column and continued
below it: rho_s phi_s is the differential incidence d^2N / dz dlogN, which
does not depend on the reference), g_s the distribution of b (normalized).
The factor exp(-expected number) of the Poisson process is the same for
every hypothesis of a window, and left out. A component of a species at
the redshift of a known system (the sweep, section 3.5 of the design): rho
the probability that the system has the species over the redshift range of
the profile of the anchor (uniform in z there).

The populations (:data:`POPULATIONS`) come from the literature, never from
the generator of the mocks (``astrocook.mocks.populations``): the mocks are
built from similar literature, so they are an optimistic test of these
priors; the twins and the development spectra are the real one. The values
were collected in a session without access to arXiv (PR 7c.4d-quater): the
numbers marked "to verify" are as quoted by later papers and reviews.

Not part of the evidence: the clustering of the components of a system (a
second CIV component next to a first is far more likely than the incidence
says: the evidence penalizes a system of many components), the bounds of the
fit (the Gaussian approximation ignores them), the evolution of the b with z
and with the column.
"""
import dataclasses
import math
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from astrocook.core.system_list import SystemList
from astrocook.fitting.lsf import C_KMS

LN10 = math.log(10.0)


@dataclasses.dataclass(frozen=True)
class DopplerPrior:
    """The distribution of b (km/s) of a population, normalized:
    ``"gauss"`` (``a`` the mean, ``s`` the dispersion), ``"lognormal"``
    (``a`` the median, ``s`` the dispersion of ln b) or ``"loguniform"``
    (between ``a`` and ``s``)."""
    kind: str
    a: float
    s: float

    def terms(self, b: float) -> Tuple[float, float, float]:
        """ln g(b), its derivative and the curvature of -ln g (its concave
        part only, >= 0)."""
        b = max(float(b), 1e-3)
        if self.kind == "gauss":
            lp = -0.5 * ((b - self.a) / self.s) ** 2 - math.log(self.s * math.sqrt(2 * math.pi))
            return lp, -(b - self.a) / self.s ** 2, 1.0 / self.s ** 2
        if self.kind == "lognormal":
            u = (math.log(b) - math.log(self.a)) / self.s
            lp = -math.log(b * self.s * math.sqrt(2 * math.pi)) - 0.5 * u * u
            d1 = -(1.0 + u / self.s) / b
            d2 = (1.0 + (u - 1.0 / self.s) / self.s) / b ** 2
            return lp, d1, max(0.0, -d2)
        if self.kind == "loguniform":
            lp = -math.log(b * math.log(self.s / self.a))
            return lp, -1.0 / b, 0.0
        raise ValueError(f"unknown prior of b '{self.kind}'")


@dataclasses.dataclass(frozen=True)
class Population:
    """The prior of the components of a species: the incidence above
    ``logn_ref`` (components per unit z at ``z_ref``, evolving as
    ((1 + z) / (1 + z_ref))^``gamma``), the slope ``beta`` of the column
    density distribution f(N) ~ N^-beta (None: log N uniform over
    ``logn_range``), the distribution of b, and the sources."""
    name: str
    rate: float
    z_ref: float
    gamma: float
    logn_ref: float
    beta: Optional[float]
    b: DopplerPrior
    source: str
    logn_range: Tuple[float, float] = (11.0, 15.0)

    def log_rate(self, z: float) -> float:
        """ln of the incidence per unit z above ``logn_ref`` at ``z``."""
        return math.log(self.rate) + self.gamma * math.log((1.0 + z) / (1.0 + self.z_ref))

    def logn_terms(self, logn: float) -> Tuple[float, float, float]:
        """ln phi(log N), its derivative and the curvature of -ln phi."""
        if self.beta is None:
            lo, hi = self.logn_range
            return -math.log(hi - lo), 0.0, 0.0
        slope = (self.beta - 1.0) * LN10
        return math.log(slope) - slope * (logn - self.logn_ref), -slope, 0.0


#: The b of the forest (Hu et al. 1995, AJ 110, 1526: at z ~ 3 a Gaussian of
#: mean 28 km/s and dispersion 10 km/s above a cutoff of about 20 km/s, as
#: summarized by Rauch 1998, ARA&A 36, 267). The cutoff is not applied: a
#: hard one makes every H I narrower than 20 km/s impossible, which the noise
#: and the blends of a real fit do not allow (the narrow lines of the forest
#: are mostly metals, which is what the prior should say, but not with
#: certainty); the Gaussian continued below it gives a line of 6 km/s about
#: 1/10 of the density of the peak
B_FOREST = DopplerPrior("gauss", 28.0, 10.0)
#: The b of the metal components (Rauch, Sargent, Womble & Barlow 1996,
#: ApJ 467, L5: CIV components with median total b of about 10 km/s, fewer
#: than 8% below 6 km/s above log N 12.5; mean thermal b 7.2 km/s). A
#: lognormal of median 10 km/s and dispersion 0.4 in ln b (10% below 6
#: km/s). The same for the low ions (narrower in the damped systems: the
#: prior is broad for them), to verify against Boksenberg & Sargent (2015,
#: ApJS 218, 7) for the other ions
B_METAL = DopplerPrior("lognormal", 10.0, 0.4)
#: The b of a line of unknown species: log-uniform over the bounds of b of
#: the short fit
B_ANY = DopplerPrior("loguniform", 1.0, 100.0)

#: The populations of the species, by the species of
#: ``hypotheses_sweep.SWEEP_SPECIES`` (``"HI"`` the forest). Sources:
#:
#: - H I: the incidence dn/dz = 6.1 (1 + z)^2.47 for 13.64 < log N < 17 at
#:   1.5 < z < 4 (Kim, Carswell, Cristiani, D'Odorico & Giallongo 2002,
#:   MNRAS 335, 555; to verify: as quoted with Kim, Cristiani & D'Odorico
#:   2001, A&A 373, 757, who give gamma 2.19 for log N 13.64-16); the slope
#:   of f(N), beta = 1.5 below log N 14 (Kim et al. 2001: 1.4-1.5; Hu et
#:   al. 1995: 1.46); b :data:`B_FOREST`.
#: - CIV: about 10 components per unit z above log N 12.5 at z 2-3 (the
#:   order of the systems of D'Odorico, Calura, Cristiani & Viel 2010,
#:   MNRAS 401, 2715, 5-10 per unit z, design section 3.4, with about two
#:   components per system: 1023 components in 508 systems there), no
#:   evolution over 1.5 < z < 4 at this precision; beta 1.8 (to verify:
#:   Songaila 2001, ApJ 561, L153; Ellison et al. 2000 find 1.44); b
#:   :data:`B_METAL`.
#: - SiIV: a third of the CIV (to verify: the SiIV/CIV of the systems of
#:   Boksenberg & Sargent 2015), the slope and b of the CIV.
#: - MgII: about 3 systems per unit z above W 0.02 A at z ~ 1 (the weak
#:   systems of Churchill, Rigby, Charlton & Vogt 1999, ApJS 120, 51, 1.74
#:   per unit z for 0.02 < W < 0.3 A, and the strong ones of Nestor,
#:   Turnshek & Rao 2005, ApJ 628, 637, about 0.8 above 0.3 A), about 5
#:   components above log N 12; the slope and b of the CIV (to verify).
#: - The other metals (NV, OVI, AlIII, CII, CIII, SiII, SiIII, OI, AlII,
#:   FeII, MgI and the ions of the damped systems): no incidence of their
#:   own here, they enter only at the redshift of a known system (the
#:   sweep, :func:`at_system`); as free lines, the incidence of the SiIV
#:   (a placeholder, not from the literature).
#: - ``"unidentified"``: a line of unknown species, its incidence measured on
#:   the spectrum itself (the regions per unit z of the placeholder Ly-a in
#:   their zone, :func:`unidentified`), log N uniform over 11-15, b
#:   :data:`B_ANY`.
POPULATIONS: Dict[str, Population] = {
    "HI": Population("HI", 6.1 * 4.0 ** 2.47, 3.0, 2.47, 13.64, 1.5, B_FOREST,
                     "Kim et al. 2002; Kim et al. 2001; Hu et al. 1995"),
    "CIV": Population("CIV", 10.0, 2.5, 0.0, 12.5, 1.8, B_METAL,
                      "D'Odorico et al. 2010; Songaila 2001; Rauch et al. 1996"),
    "SiIV": Population("SiIV", 3.0, 2.5, 0.0, 12.5, 1.8, B_METAL,
                       "a third of the CIV; Rauch et al. 1996"),
    "MgII": Population("MgII", 5.0, 1.0, 0.0, 12.0, 1.8, B_METAL,
                       "Churchill et al. 1999; Nestor et al. 2005; Rauch et al. 1996"),
}
#: The probability that a known system has a species (the sweep): not from
#: the literature, a choice of equal odds (present or not), the same for
#: every species; the H I of a metal system is always there (1)
P_AT_SYSTEM = 0.5
#: The slope of the columns of a species at a known system, above
#: :data:`LOGN_AT_SYSTEM` (that of the CIV)
BETA_AT_SYSTEM = 1.8
LOGN_AT_SYSTEM = 12.0
#: The weak curvature added to every parameter, against a Jacobian that
#: does not constrain it (a component without influence would have an
#: infinite evidence): a Gaussian of these widths (z in km/s, log N, b in
#: km/s), far wider than any measured uncertainty
WEAK_WIDTHS = (200.0, 3.0, 50.0)


def species_of(series: str) -> str:
    """The species of a series (``"Ly_a"``, ``"Ly_ab"`` -> ``"HI"``,
    ``"OI_1039,OI_1302"`` -> ``"OI"``, ``"CIV"`` -> ``"CIV"``)."""
    if series.startswith("Ly"):
        return "HI"
    return series.split(",")[0].split("_")[0]


@dataclasses.dataclass(frozen=True)
class ComponentPrior:
    """The prior of one component: ln of its density in z (the incidence
    per unit z, or the probability at a known system over its range in z),
    and the distributions of its log N and b."""
    log_rate: float
    population: Population

    def terms(self, logn: float, b: float) -> Tuple[float, np.ndarray, np.ndarray]:
        """ln pi, its gradient and the curvature of -ln pi in (z, log N, b)."""
        ln_n, dn, cn = self.population.logn_terms(logn)
        ln_b, db, cb = self.population.b.terms(b)
        return (self.log_rate + ln_n + ln_b, np.array([0.0, dn, db]),
                np.array([0.0, cn, cb]))


def free_line(species: str, z: float) -> ComponentPrior:
    """The prior of a free line of ``species`` (its incidence at ``z``); a
    species without a population of its own takes that of the SiIV."""
    pop = POPULATIONS.get(species, POPULATIONS["SiIV"])
    return ComponentPrior(pop.log_rate(z), pop)


def unidentified(rate_per_z: float) -> ComponentPrior:
    """The prior of an unidentified line, ``rate_per_z`` lines per unit z
    (measured on the spectrum)."""
    pop = Population("unidentified", max(rate_per_z, 1e-3), 0.0, 0.0, 11.0, None, B_ANY,
                     "the regions of the spectrum")
    return ComponentPrior(pop.log_rate(0.0), pop)


def at_system(species: str, z_range: Tuple[float, float],
              probability: Optional[float] = None) -> ComponentPrior:
    """The prior of a component of ``species`` at a known system: the
    probability that the system has it (:data:`P_AT_SYSTEM`, 1 for the
    H I) over the range ``z_range`` (uniform in z), the columns a power law
    (the H I that of the forest, the metals :data:`BETA_AT_SYSTEM` above
    :data:`LOGN_AT_SYSTEM`), b that of the population (the H I of the
    forest, the metals :data:`B_METAL`)."""
    if probability is None:
        probability = 1.0 if species == "HI" else P_AT_SYSTEM
    width = max(z_range[1] - z_range[0], 1e-7)
    if species == "HI":
        pop = POPULATIONS["HI"]
    else:
        pop = Population(species, 1.0, 0.0, 0.0, LOGN_AT_SYSTEM, BETA_AT_SYSTEM, B_METAL,
                         "at a known system")
    return ComponentPrior(math.log(probability) - math.log(width), pop)


@dataclasses.dataclass(frozen=True)
class Evidence:
    """The terms of the Laplace evidence of a hypothesis."""
    ln_z: float
    chi2: float
    ln_prior: float
    """ln pi at the point of the fit."""
    occam: float
    """(k/2) ln 2 pi - (1/2) ln det H: the volume of the posterior."""
    shift: float
    """(1/2) g^T H^-1 g: the step to the maximum of the posterior."""
    k: int

    @property
    def criterion(self) -> float:
        """-2 ln Z, in the units of the information criteria."""
        return -2.0 * self.ln_z


def laplace(chi2: float, jac: np.ndarray, resid: np.ndarray, ln_prior: float,
            grad_prior: np.ndarray, curv_prior: np.ndarray) -> Evidence:
    """The Laplace evidence from the chi2, the Jacobian ``jac`` of the
    normalized residuals (pixels x k; ``resid`` kept for the record), and ln
    pi with its gradient
    and the curvature of -ln pi (diagonal, k)."""
    k = jac.shape[1]
    if k == 0:
        return Evidence(-0.5 * chi2 + ln_prior, chi2, ln_prior, 0.0, 0.0, 0)
    h = jac.T @ jac + np.diag(curv_prior)
    g = np.asarray(grad_prior, dtype=float)
    sign, logdet = np.linalg.slogdet(h)
    if sign <= 0:
        raise np.linalg.LinAlgError("the curvature of the posterior is not positive")
    shift = 0.5 * float(g @ np.linalg.solve(h, g))
    occam = 0.5 * k * math.log(2 * math.pi) - 0.5 * float(logdet)
    return Evidence(-0.5 * chi2 + ln_prior + shift + occam, chi2, ln_prior, occam, shift, k)


def window_evidence(win, systs: SystemList, priors: Mapping[str, ComponentPrior]) -> Evidence:
    """The Laplace evidence of the hypothesis ``systs`` on the pixels of the
    :class:`~astrocook.fitting.evaluation_window.EvaluationWindow` ``win``
    (the fixed components of the window frozen), ``priors`` by the uuid of
    each component of ``systs`` (all of them). The chi2 is that of
    ``win.evaluate``.

    Raises
    ------
    KeyError
        If a component has no prior.
    numpy.linalg.LinAlgError
        If the curvature of the posterior is not positive.
    """
    from astrocook.fitting.voigt_fitter import VoigtFitter
    own = [c.uuid for c in systs.components]
    model = win._with_fixed(systs)
    n = len(win.index)
    if not own:
        ev = win.evaluate(systs, min_influence=0.0)
        return laplace(ev.chi2, np.zeros((n, 0)), np.zeros(n), 0.0, np.zeros(0), np.zeros(0))
    with model.fitting_context(own, group_depth=0):
        fitter = VoigtFitter(win.spectrum, model)
        fitter.prepare_fit_context(pixels=np.ones(n, dtype=bool))
        cm = fitter._constraints
        p = np.asarray(cm.p_free_vector, dtype=float)
        resid = fitter._residual_function(p)
        jac = fitter._jacobian(p)
        p_full = cm.map_p_free_to_full(p)
        is_free = np.asarray(cm._param_map['is_free'], dtype=bool)
        fitter._active_fit_components = None
        fitter._active_lines_data = None
    free_index = np.cumsum(is_free) - 1
    k = len(p)
    ln_prior = 0.0
    grad = np.zeros(k)
    curv = np.zeros(k)
    for i, comp in enumerate(model.components):
        if comp.uuid not in priors:
            if any(is_free[4 * i:4 * i + 4]):
                raise KeyError(f"no prior for the component {comp.uuid} ({comp.series})")
            continue
        z, logn, b = (float(p_full[4 * i + j]) for j in range(3))
        lp, gp, cp = priors[comp.uuid].terms(logn, b)
        ln_prior += lp
        weak = np.array([(1.0 + z) * WEAK_WIDTHS[0] / C_KMS, WEAK_WIDTHS[1], WEAK_WIDTHS[2]])
        for j in range(3):
            if is_free[4 * i + j]:
                m = free_index[4 * i + j]
                grad[m] += gp[j]
                curv[m] += cp[j] + 1.0 / weak[j] ** 2
    return laplace(float(np.sum(resid ** 2)), jac, resid, ln_prior, grad, curv)


def region_rates(centres_z: Sequence[float], zones: Sequence[str],
                 coverage: Mapping[str, Tuple[float, float]]) -> Dict[str, float]:
    """The regions per unit z (of the placeholder Ly-a) in each zone: the
    incidence of an unidentified line measured on the spectrum (F174
    direction 2). ``coverage`` the range of z of each zone."""
    out: Dict[str, float] = {}
    for zone, (lo, hi) in coverage.items():
        n = sum(1 for z, zn in zip(centres_z, zones, strict=True) if zn == zone)
        out[zone] = max(n, 1) / max(hi - lo, 1e-6)
    return out


def describe() -> List[str]:
    """The populations, one line each (for the logs and the docs)."""
    return [f"{p.name}: {p.rate:g} per unit z above log N {p.logn_ref:g} at z {p.z_ref:g}, "
            f"gamma {p.gamma:g}, beta {p.beta}, b {p.b.kind}({p.b.a:g}, {p.b.s:g}); {p.source}"
            for p in POPULATIONS.values()]
