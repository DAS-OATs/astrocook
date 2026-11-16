"""The evaluation of a model on a fixed set of pixels, without fitting
(findings register F154, phase 7, PR 7a.5).

The chi2 of a model and its information criteria (AIC, AICc, BIC and the
SpIC of Webb et al. 2021) on pixels chosen by the user, so that two
alternative models can be compared on the same pixels: the basis of the
comparison between alternative identifications (phase 7c) and of the
evaluation of an imported model (7a.3).

The model is the whole system list; a selection of its components decides
only which of them are counted: their free parameters (k) and their terms of
the SpIC. Of the selection, only the components that influence the pixels
are counted (:data:`MIN_INFLUENCE`).

SpIC (Webb, Lee, Carswell & Milakovic 2021, MNRAS 501, 2268,
doi:10.1093/mnras/staa3551, eq. 3, 5, 6)::

    SpIC = chi2 + sum_a [ 2 f k_a R_a / (R_a - k_a - 1) + (1 - f) k_a ln R_a ]
    R_a  = sum_i (1 - I^a_i) / sigma_i

over the velocity components a, with k_a their free parameters
(sum_a k_a = k), I^a_i the normalized profile of a alone at pixel i and
sigma_i the error of the normalized flux. R_a replaces the number of pixels
N of AICc and BIC: the pixels a component influences, weighted by its depth
and by the S/N. f = 1 is the analogue of AICc (SpIC_A), f = 0 of BIC
(SpIC_B), f = 1/2 their hybrid (SpIC_H), which the paper finds the best
all-rounder (section 5.3) and is ``spic`` here. R_a is bounded from below by
k_a + 2, where the first term would diverge (section 2.4).

A velocity component of the paper is one redshift shared by several
species (k_a = 6 for three species there): here, the components whose
redshifts are linked by the constraints (``z`` of one an expression of the
``z`` of another), one component alone when its redshift is not linked.
Not the groups of the fits (:meth:`SystemList.get_connected_group`), which
also join the blended components.
"""
import dataclasses
import logging
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from astrocook.core.spectrum import Spectrum
from astrocook.core.system_list import SystemList
from astrocook.fitting.voigt_fitter import VoigtFitter, chi2_of

logger = logging.getLogger(__name__)

#: The weight f of the AICc-like term of the SpIC (Webb et al. 2021, eq. 3),
#: for each SpIC reported: SpIC_H (``spic``, the one the paper recommends),
#: SpIC_A and SpIC_B
SPIC_WEIGHTS = {"spic": 0.5, "spic_a": 1.0, "spic_b": 0.0}

#: The influence of a component on the pixels, sum_i ((1 - I_i) / sigma_i)^2,
#: above which it is counted (its free parameters in k, its terms in the
#: SpIC). For a weak, isolated component the sum is about the increase of the
#: chi2 if it were removed: 1 means that the pixels see it at the level of
#: 1 sigma. Without it, the components far from the pixels counted too (on a
#: window of 5 nm of a mock with 140 components, k = 420 on 1286 pixels and a
#: reduced chi2 of 1.47 for a chi2/n of 0.99); a rule on the centres of the
#: lines missed the wings that do absorb in the pixels, so that a component
#: at the edge of a window improved the chi2 for free (Guido, 16 November
#: 2026: this criterion, threshold 1, variable)
MIN_INFLUENCE = 1.0

#: |residual| above which a pixel counts in ``frac_beyond_3sigma``
OUTLIER_SIGMA = 3.0

_UUID_IN_EXPRESSION = re.compile(r"""p\[\s*['"]([^'"]+)['"]\s*\]""")


@dataclasses.dataclass(frozen=True)
class VelocityComponent:
    """The selected components of one velocity component of the SpIC, their
    free parameters (k_a) and their R_a (before its lower bound)."""
    uuids: List[str]
    k: int
    r: float

    def as_dict(self) -> Dict[str, Any]:
        return {"uuids": list(self.uuids), "k": self.k, "R": float(self.r)}


@dataclasses.dataclass(frozen=True)
class ModelEvaluation:
    """The measure of a model on a fixed set of pixels."""
    chi2: float
    n_pixels: int
    k: int
    dof: int
    chi2_reduced: Optional[float]
    aic: float
    aicc: Optional[float]
    bic: Optional[float]
    spic: float
    spic_a: float
    spic_b: float
    resid_mean: Optional[float]
    frac_beyond_3sigma: Optional[float]
    components: List[str]
    #: The candidates left out for their influence below the threshold
    n_not_counted: int
    velocity_components: List[VelocityComponent]
    #: Normalized residuals on the whole spectrum, NaN out of the pixels
    residuals: np.ndarray = dataclasses.field(repr=False, compare=False)

    def as_dict(self) -> Dict[str, Any]:
        """The numbers, for the ``extra`` of the result and the agents (the
        residuals left out)."""
        out = {f.name: getattr(self, f.name) for f in dataclasses.fields(self)
               if f.name not in ("residuals", "velocity_components")}
        out["velocity_components"] = [v.as_dict() for v in self.velocity_components]
        return out


def information_criteria(chi2: float, n: int, k: int
                         ) -> Tuple[float, Optional[float], Optional[float]]:
    """(AIC, AICc, BIC): AIC = chi2 + 2k, AICc = AIC + 2k(k+1)/(n-k-1) (None if n - k - 1 <= 0)
    and BIC = chi2 + k ln n (None without pixels); chi2 stands for
    -2 ln L, whose constant does not matter between models on the same
    pixels (as in ``metals_auto``)."""
    aic = chi2 + 2 * k
    aicc = aic + 2 * k * (k + 1) / (n - k - 1) if n - k - 1 > 0 else None
    bic = chi2 + k * np.log(n) if n > 0 else None
    return (float(aic), None if aicc is None else float(aicc),
            None if bic is None else float(bic))


def spic(chi2: float, terms: Sequence[VelocityComponent], f: float) -> float:
    """The SpIC (Webb et al. 2021, eq. 3) for the weight ``f``; R_a bounded
    from below by k_a + 2 (section 2.4), a velocity component with no free
    parameter adds nothing."""
    total = chi2
    for t in terms:
        if t.k == 0:
            continue
        r = max(t.r, t.k + 2.0)
        total += 2 * f * t.k * r / (r - t.k - 1) + (1 - f) * t.k * np.log(r)
    return float(total)


def free_parameters(fitter: VoigtFitter) -> Dict[str, int]:
    """The free parameters of each component (frozen and linked ones do not
    count), from the constraint model of ``fitter`` with all components
    active."""
    is_free = fitter._constraints._param_map['is_free'].reshape(-1, 4)
    return {c.uuid: int(np.sum(row)) for c, row in
            zip(fitter._system_list.components, is_free, strict=True)}


def velocity_components(systs: SystemList, constraints: Dict[str, Dict[str, Any]]
                        ) -> List[List[str]]:
    """The components joined by the links of their redshifts (an expression
    of ``z`` that names another component), in the order of the list."""
    order = [c.uuid for c in systs.components]
    parent = {u: u for u in order}

    def root(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u

    for u in order:
        z = constraints.get(u, {}).get('z')
        expression = getattr(z, 'expression', None) if z is not None else None
        if not expression:
            continue
        for target in _UUID_IN_EXPRESSION.findall(expression):
            if target in parent:
                parent[root(u)] = root(target)
    groups: Dict[str, List[str]] = {}
    for u in order:
        groups.setdefault(root(u), []).append(u)
    return list(groups.values())


def evaluate_model(spec: Spectrum, systs: SystemList, mask: np.ndarray,
                   uuids: Optional[Sequence[str]] = None,
                   min_influence: float = MIN_INFLUENCE) -> ModelEvaluation:
    """Measure the model of ``systs`` (all its components) on the pixels
    ``mask`` (those without a valid flux or error left out), counting the
    components of ``uuids`` (None: all) whose influence on the pixels is at
    least ``min_influence`` (0: all of them; see :data:`MIN_INFLUENCE`).

    Raises
    ------
    ValueError
        If no valid pixel is in ``mask``.
    """
    fitter = VoigtFitter(spec, systs)
    pixels = np.asarray(mask, dtype=bool) & fitter.valid_mask
    if not np.any(pixels):
        raise ValueError("no valid pixel to evaluate the model on")
    if systs.components:
        _, model_flux = fitter.compute_model_flux()
    else:
        norm = spec.norm
        assert norm is not None  # VoigtFitter refuses a spectrum without it
        model_flux = np.asarray(norm.value, dtype=float)
    resid = fitter.normalized_residuals(model_flux, pixels)
    chi2 = chi2_of(resid)
    finite = resid[np.isfinite(resid)]
    n = int(np.sum(pixels))

    candidates = [c.uuid for c in systs.components] if uuids is None else list(uuids)
    k_of = free_parameters(fitter) if systs.components else {}
    counted: List[str] = []
    terms: List[VelocityComponent] = []
    if candidates:
        idx = np.flatnonzero(pixels)
        sl = slice(int(idx[0]), int(idx[-1]) + 1)
        profiles = fitter.component_profiles(candidates, sl)
        in_sl = pixels[sl]
        sigma = fitter._dy_norm[sl][in_sl]
        depth = {u: (1.0 - prof[in_sl]) / sigma for u, prof in profiles.items()}
        counted = [u for u in candidates
                   if min_influence <= 0 or float(np.sum(depth[u] ** 2)) >= min_influence]
        selected = set(counted)
        r_of = {u: float(np.sum(depth[u])) for u in counted}
        for group in velocity_components(systs, fitter._constraints.v2_constraints_by_uuid):
            members = [u for u in group if u in selected]
            if members:
                terms.append(VelocityComponent(members, sum(k_of[u] for u in members),
                                               sum(r_of[u] for u in members)))
    k = sum(k_of[u] for u in counted)

    residuals = np.full(len(spec.x), np.nan)
    residuals[pixels] = resid
    dof = n - k
    aic, aicc, bic = information_criteria(chi2, n, k)
    return ModelEvaluation(
        chi2=chi2, n_pixels=n, k=k, dof=dof,
        chi2_reduced=chi2 / dof if dof > 0 else None,
        aic=aic, aicc=aicc, bic=bic,
        **{name: spic(chi2, terms, f) for name, f in SPIC_WEIGHTS.items()},
        resid_mean=float(np.mean(finite)) if finite.size else None,
        frac_beyond_3sigma=(float(np.mean(np.abs(finite) > OUTLIER_SIGMA))
                            if finite.size else None),
        components=counted, n_not_counted=len(candidates) - len(counted),
        velocity_components=terms, residuals=residuals)
