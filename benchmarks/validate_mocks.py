"""Validation of the recipes on mocks (phase 4 of ROADMAP.md).

Runs the chain of the MCP pilot (Continuum > Auto-estimate, Absorbers >
Identify lines, Absorbers > Auto-detect doublets, with the defaults of the
GUI dialogs) on a grid of mock quasar spectra and measures the results
against the truth:

- continuum: ratio of the estimated to the true continuum, in the
  Lyman-alpha forest and on the red side;
- identification of the metal doublets (CIV, SiIV, MgII): purity and
  completeness, with duplicates (F49) and candidates beyond z_em (F55)
  counted apart;
- fit (``doublets_auto``): bias and pulls of z, logN and b of the matched
  components, spurious and missed components, components at a bound of
  the fit or with absurd errors (#517, F52).

Not a test (not collected by pytest): run it by hand, before and after a
change, and compare the summaries. Usage::

    python benchmarks/validate_mocks.py [--grid full|robustness|reduced|tiny]
        [--workers N] [--json out.json] [--compare] [--compare-baseline]
        [--save-reference]

The full grid (135 mocks: S/N 10, 30, 100; R 8900, 45000, 140000; z_em
2.5, 3.5, 4.5; 5 seeds) was agreed with Guido at the start of phase 4,
with a matching tolerance of 20 km/s.

Each mock covers 104-160 nm in the rest frame of the quasar (from the
Lyman-beta forest to beyond the CIV emission), with 2.5 pixels per
resolution element, a power-law continuum with emission lines, a
Lyman-alpha forest (the whole Lyman series for logN >= 15), metal systems
on the red side with their H I, the Lyman limit, and a noise that follows
the flux. Seeds 0 and 1 also have a DLA and a sub-DLA in the forest
(``strong`` in the records).

The grid ``robustness`` (phase 7, PR 7a.1; 252 mocks: the cases of
``ROBUSTNESS_CASES``, S/N 10, 30, 100; R 8900, 45000; z_em 3.5, 4.5, with
the tellurics in the window; 3 seeds) adds one case at a time to the same
forest and metals (without the DLA and the sub-DLA of ``full``): a DLA, a
sub-DLA, a BAL trough, the richer emission lines, the tellurics, the
proximity zone with associated absorbers, or all of them (``combined``).
The forest and metals of a grid point are the same for every case (the
draws of each case come from a generator of their own). The recipes get
a z_em perturbed by a draw of sigma 500 km/s (reproducible from the
seed; ``z_em_used`` in the records, ``z_em`` the true one).

The metrics of phase 7 (PR 7a.2, ``phase7`` in the records and the
summaries, on both grids) use the truth map of each mock
(:func:`astrocook.mocks.truth_map`): the continuum on the damped wings of
the DLAs and sub-DLAs, across the BAL troughs, on the peaks of the
emission lines and in the proximity zone (ratio, integrated deviation, and
the pull against the column ``dcont`` when it exists); what the truth
leaves to constrain the continuum in the forest; z_em of an estimator
against the truth (when the chain has one, F153); the false candidates on
tellurics; the absorbed pixels against the regions and identifications of
identify_lines (the Lyman-alpha by default of F51 apart); a confusion
table of the flags of the fitted components against their errors.

The reference, the current algorithms on both grids, is
``results/validate_mocks_robustness_baseline.json`` (``--save-reference``;
``--compare`` prints the numbers of phase 4 and of phase 7 of a run next to
it), regenerated in PR 7c.1 after the fix of the scoring of the doublets
(F68). The baseline of phase 4, ``results/validate_mocks_baseline.json``
(2.2.0), is kept as history (``--compare-baseline``). The truth map costs
15-60 s per mock.
"""
import argparse
import dataclasses
import json
import logging
import sys
import time
from itertools import pairwise
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS  # noqa: E402
from astrocook.fitting.component_flags import MEASUREMENT_FLAGS  # noqa: E402
from astrocook.fitting.voigt_fitter import COMPONENT_FLAGS  # noqa: E402
from astrocook.mocks import (  # noqa: E402
    DLA_LOGN, LYMAN_SERIES, QSO_EMISSION_LINES, QSO_EMISSION_LINES_FULL,
    SUBDLA_LOGN, BALTrough, MockTruth, TelluricModel, associated_absorbers,
    continuum, generate, lya_forest, metal_systems, optical_depth,
    strong_hi_systems,
)
from astrocook.mocks.physics import C_KMS, transitions  # noqa: E402
from astrocook.mocks.truthmap import (  # noqa: E402
    CATEGORIES, component_category, truth_map,
)

GRIDS = {
    "full": dict(snr=(10.0, 30.0, 100.0), resol=(8900.0, 45000.0, 140000.0),
                 z_em=(2.5, 3.5, 4.5), seeds=tuple(range(5))),
    "reduced": dict(snr=(10.0, 30.0), resol=(8900.0, 45000.0), z_em=(2.5, 3.5),
                    seeds=(0, 2)),
    "tiny": dict(snr=(30.0,), resol=(8900.0,), z_em=(2.5,), seeds=(2,)),
    "robustness": dict(snr=(10.0, 30.0, 100.0), resol=(8900.0, 45000.0),
                       z_em=(3.5, 4.5), seeds=(0, 1, 2),
                       cases=("dla", "subdla", "bal", "emission", "telluric",
                              "proximity", "combined")),
}
#: The cases of the grid robustness (phase 7, PR 7a.1)
ROBUSTNESS_CASES = GRIDS["robustness"]["cases"]
#: Sigma of the perturbation of the z_em given to the recipes in the grid
#: robustness (km/s; Guido, 13 November 2026)
Z_EM_SIGMA_KMS = 500.0
#: Velocity range of the associated absorbers from z_em (km/s, positive
#: beyond it), within the margin of the identification (F55, 3000 km/s),
#: and that of the one beyond the margin, in the mocks of seed 2
ASSOCIATED_KMS = (-2000.0, 3000.0)
ASSOCIATED_BEYOND_KMS = (3500.0, 6000.0)

#: Matching tolerance between candidates or fitted components and the truth
MATCH_KMS = 20.0
#: Rest-frame window of the mocks (nm)
REST_WINDOW_NM = (104.0, 160.0)
#: Pixels per resolution element (FWHM)
PIXELS_PER_FWHM = 2.5
#: Largest sub-pixel of the mock model (km/s); the narrowest lines have
#: b = 4 km/s
MAX_SUBPIXEL_KMS = 1.0
#: Variance at zero flux relative to the continuum (MockTruth.noise_floor)
NOISE_FLOOR = 0.1
#: Detection threshold of the weaker line of a doublet, in sigma of its
#: equivalent width (proposed in the plan of phase 4, to be confirmed)
DETECTION_SIGMA = 3.0
#: Metal systems per mock on the red side, and their species
METAL_SYSTEMS = {"CIV": 4, "SiIV": 2, "MgII": 2}
METALS = tuple(METAL_SYSTEMS)
#: Bounds of the fits, for the flags (VoigtModelConstraint.DEFAULT_BOUNDS
#: and the max_b and max_logN of doublets_auto; the upper bound on logN was
#: 23 before F62)
LOGN_BOUNDS = (10.0, 19.0)
B_BOUNDS = (1.0, 50.0)
LYA_NM, LYB_NM = 121.567, 102.572


# --- Mocks --------------------------------------------------------------

def _doublet_z_range(species, z_em):
    """Redshifts at which both lines of a doublet fall on the red side of
    the window (between the Lyman-alpha emission and the window end)."""
    waves = [ATOM_DATA[t]["wave"] / 10.0 for t in STANDARD_MULTIPLETS[species]]
    lo = LYA_NM * (1 + z_em) / min(waves) - 1
    hi = REST_WINDOW_NM[1] * (1 + z_em) / max(waves) - 1
    return lo + 0.01, min(hi, z_em) - 0.01


def _rng(snr, resol, z_em, seed, *extra):
    return np.random.default_rng([int(seed), int(snr), int(resol), int(z_em * 10),
                                  *extra])


def _window(z_em):
    return tuple(w * (1 + z_em) for w in REST_WINDOW_NM)


def _forest_and_metals(rng, z_em, proximity_kms=None):
    """The forest and the metal systems on the red side of a mock."""
    lam_min, _ = _window(z_em)
    forest = lya_forest(rng, lam_min / LYA_NM - 1, z_em, series="Ly_ab",
                        strong_series=LYMAN_SERIES, z_em=z_em,
                        proximity_kms=proximity_kms)
    comps, systems = list(forest), []
    for species, n in METAL_SYSTEMS.items():
        z_lo, z_hi = _doublet_z_range(species, z_em)
        for _ in range(n):
            system = metal_systems(rng, 1, z_lo, z_hi, species=(species,),
                                   hi_logN_range=(14.0, 16.0))
            metals = [c for c in system if c.series == species]
            weights = np.array([10**c.logN for c in metals])
            z_sys = float(np.sum(weights * [c.z for c in metals]) / weights.sum())
            systems.append(dict(species=species, z=z_sys, components=metals))
            comps += system
    return comps, systems


def _strong(rng, z_em, logN_range):
    """One strong H I absorber, its Lyman-alpha in the forest of the
    window, away from the edges."""
    lam_min, _ = _window(z_em)
    return strong_hi_systems(rng, 1, lam_min * 1.05 / LYA_NM - 1, z_em - 0.05,
                             logN_range=logN_range)


def _truth(comps, snr, resol, z_em, seed, **kwargs):
    lam_min, lam_max = _window(z_em)
    dv = C_KMS / resol / PIXELS_PER_FWHM
    kwargs.setdefault("emission_lines", QSO_EMISSION_LINES)
    return MockTruth(
        components=sorted(comps, key=lambda c: c.z), lam_min_nm=lam_min,
        lam_max_nm=lam_max, dv_kms=dv, resol=resol, snr=snr, z_em=z_em,
        seed=int(seed), oversample=max(1, int(np.ceil(dv / MAX_SUBPIXEL_KMS))),
        cont_slope=-1.0, lyman_limit=True, noise_floor=NOISE_FLOOR, **kwargs)


def make_truth(snr, resol, z_em, seed):
    """Mock of one grid point, and its metal systems.

    Returns
    -------
    truth : MockTruth
    systems : list of dict
        ``{"species", "z", "components": [MockComponent, ...]}``, one per
        metal system (without its H I).
    """
    rng = _rng(snr, resol, z_em, seed)
    comps, systems = _forest_and_metals(rng, z_em)
    strong = None
    if seed in (0, 1):
        strong = "DLA" if seed == 0 else "subDLA"
        comps += _strong(rng, z_em, DLA_LOGN if seed == 0 else SUBDLA_LOGN)
    return _truth(comps, snr, resol, z_em, seed), systems, strong


# --- Robustness grid (phase 7, PR 7a.1) ---------------------------------

def _bal_trough(rng):
    """A BAL trough: its start 2000-10000 km/s, its width 2000-10000 km/s,
    its CIV optical depth 0.3-3 (SiIV has 2.4 times that, at the same
    column: the troughs range from shallow to black)."""
    v_min = rng.uniform(2000.0, 10000.0)
    width = rng.uniform(2000.0, 10000.0)
    return BALTrough(float(v_min), float(v_min + width),
                     float(rng.uniform(0.3, 3.0)), edge_kms=min(1000.0, width / 4))


def _emission_lines(rng):
    """QSO_EMISSION_LINES_FULL with the amplitudes scaled by 0.7-1.3, the
    shifts drawn around theirs (sigma 300 km/s; CIV 0-2000 km/s to the
    blue, its asymmetry 0-0.3)."""
    lines = []
    for line in QSO_EMISSION_LINES_FULL:
        amp = line.amplitude * rng.uniform(0.7, 1.3)
        if abs(line.rest_nm - 154.9) < 0.1:
            shift, asym = -rng.uniform(0.0, 2000.0), rng.uniform(0.0, 0.3)
        else:
            shift, asym = line.shift_kms + 300.0 * rng.standard_normal(), line.asymmetry
        lines.append(dataclasses.replace(line, amplitude=float(amp),
                                         shift_kms=float(shift), asymmetry=float(asym)))
    return lines


def _telluric(rng):
    """Tellurics at airmass 1-2, PWV 0.5-10 mm, shifted by up to 2 km/s."""
    return TelluricModel(airmass=float(rng.uniform(1.0, 2.0)),
                         pwv_mm=float(rng.uniform(0.5, 10.0)),
                         shift_kms=float(rng.uniform(-2.0, 2.0)))


def z_em_perturbed(snr, resol, z_em, seed):
    """The z_em given to the recipes: the true one shifted by a draw of
    sigma Z_EM_SIGMA_KMS, reproducible from the grid point."""
    dv = Z_EM_SIGMA_KMS * _rng(snr, resol, z_em, seed, 99).standard_normal()
    return float(z_em + dv / C_KMS * (1 + z_em)), float(dv)


def make_robust_truth(snr, resol, z_em, seed, case):
    """Mock of one point of the grid robustness.

    Returns
    -------
    truth : MockTruth
    systems : list of dict
        The metal systems on the red side, as :func:`make_truth`.
    extra : dict
        What the case added: ``strong`` (``"DLA"``, ``"subDLA"``, both
        joined by ``+``, or None) and ``associated``, the associated
        systems (species, z, velocity from z_em, beyond the margin of
        F55 or not).
    """
    if case not in ROBUSTNESS_CASES:
        raise ValueError(f"Unknown case {case!r}")
    every = case == "combined"
    case_rng = _rng(snr, resol, z_em, seed, 100 + ROBUSTNESS_CASES.index(case))
    proximity_kms = (float(case_rng.uniform(1000.0, 3000.0))
                     if every or case == "proximity" else None)
    comps, systems = _forest_and_metals(_rng(snr, resol, z_em, seed), z_em,
                                        proximity_kms)
    kwargs, strong, associated = {}, [], []
    if every or case == "dla":
        comps += _strong(case_rng, z_em, DLA_LOGN)
        strong.append("DLA")
    if every or case == "subdla":
        comps += _strong(case_rng, z_em, SUBDLA_LOGN)
        strong.append("subDLA")
    if every or case == "bal":
        kwargs["bal_troughs"] = [_bal_trough(case_rng)]
    if every or case == "emission":
        kwargs["emission_lines"] = _emission_lines(case_rng)
    if every or case == "telluric":
        kwargs["telluric"] = _telluric(case_rng)
    if proximity_kms is not None:
        kwargs["proximity_kms"] = proximity_kms
        ranges = [ASSOCIATED_KMS] * 2 + ([ASSOCIATED_BEYOND_KMS] if seed == 2 else [])
        for v_range in ranges:
            system = associated_absorbers(case_rng, 1, z_em, v_range_kms=v_range)
            metals = [c for c in system if c.series != LYMAN_SERIES]
            z_sys = float(np.mean([c.z for c in metals]))
            v = (z_sys - z_em) / (1 + z_em) * C_KMS
            associated.append(dict(species=metals[0].series, z=z_sys, v_kms=v,
                                   beyond_margin=v_range is ASSOCIATED_BEYOND_KMS))
            comps += system
    truth = _truth(comps, snr, resol, z_em, seed, **kwargs)
    return truth, systems, dict(strong="+".join(strong) or None,
                                associated=associated)


# --- Detection limit ----------------------------------------------------

def weak_line(species):
    """The weaker line of a doublet (lowest f * lambda)."""
    return min(STANDARD_MULTIPLETS[species],
               key=lambda t: ATOM_DATA[t]["f"] * ATOM_DATA[t]["wave"])


def equivalent_width_ang(transition, z, logN, b):
    """Observed equivalent width (Angstrom) of one line of one component,
    from its optical depth on a fine local grid (the instrumental profile
    conserves it)."""
    lam0 = ATOM_DATA[transition]["wave"] * (1 + z)
    half = max(300.0, 10 * b)
    dv = np.linspace(-half, half, 6001)
    lam = lam0 * (1 + dv / C_KMS)
    tau = optical_depth(lam, transition, z, logN, b)
    return float(trapezoid(1 - np.exp(-tau), lam))


def is_detectable(system, truth):
    """Whether the weaker line of a system is detected at DETECTION_SIGMA.

    The equivalent width is summed over the components; its error is
    sigma_EW = dlambda_pix sqrt(N_pix) / snr over N_pix = 2 FWHM / dv
    pixels, with FWHM the instrumental one combined with the thermal one
    of the broadest component and the velocity spread of the system. An
    approximation, to be confirmed against the benchmark.
    """
    line = weak_line(system["species"])
    comps = system["components"]
    ew = sum(equivalent_width_ang(line, c.z, c.logN, c.b) for c in comps)
    zs = [c.z for c in comps]
    spread = (max(zs) - min(zs)) / (1 + system["z"]) * C_KMS
    fwhm = np.hypot(C_KMS / truth.resol, 1.665 * max(c.b for c in comps)) + spread
    lam = ATOM_DATA[line]["wave"] * (1 + system["z"])
    dlam_pix = lam * truth.dv_kms / C_KMS
    n_pix = max(1.0, 2 * fwhm / truth.dv_kms)
    sigma_ew = dlam_pix * np.sqrt(n_pix) / truth.snr
    return ew >= DETECTION_SIGMA * sigma_ew, ew, sigma_ew


# --- The chain ----------------------------------------------------------

def schema_params(menu, recipe):
    """The non-hidden parameters of a recipe with their GUI defaults, as
    strings (as the dialogs and the MCP server pass them)."""
    from astrocook.recipes.registry import get_recipe
    return {p.name: str(p.default) for p in get_recipe(menu, recipe).params
            if not p.gui_hidden and p.name != "debug_rating"}


#: The options of the chain (``--method`` and the parameters of the
#: identification by hypotheses, PR 7c.4c): set by :func:`set_chain`, also in
#: the workers of a Pool
CHAIN: dict = {}


def set_chain(options):
    """The options of the chain: ``method``, and for ``identify_lines`` the
    other keys (``multiplets``, ``prefilter``, ``screen``, ``criterion``,
    ``temper``); ``doublets_auto`` gets ``method`` (and ``multiplets`` with
    ``doublets_multiplets``)."""
    CHAIN.clear()
    CHAIN.update({k: v for k, v in options.items() if v is not None})
    from astrocook.fitting import hypotheses
    if "max_added" in CHAIN:
        # The most components the refinement of the hypotheses adds (PR
        # 7c.4c; -1: no fixed number; 0: the
        # refinement off, as ``refine=False``)
        n = int(CHAIN["max_added"])
        hypotheses.MAX_ADDED = None if n < 0 else n
    if "split_gaps" in CHAIN:
        # The split of the regions at their gaps (PR 7c.4c-bis)
        hypotheses.SPLIT_GAPS = bool(CHAIN["split_gaps"])


def chain_params(recipe):
    """The parameters of a recipe of the chain: the GUI defaults and the
    options of :data:`CHAIN`."""
    params = schema_params("absorbers", recipe)
    if "method" in CHAIN:
        params["method"] = CHAIN["method"]
    if recipe == "identify_lines":
        params.update({k: str(v) for k, v in CHAIN.items()
                       if k in ("multiplets", "prefilter", "screen", "criterion", "temper",
                                "sweep")})
    elif "doublets_multiplets" in CHAIN:
        params["multiplets"] = CHAIN["doublets_multiplets"]
    return params


def run_chain(spec):
    """Continuum, identification and doublet fit on a mock spectrum, as a
    user of the GUI would run them."""
    from astrocook.core.session import Session
    session = Session("mock", None, spec=spec)
    t0 = time.time()
    cont = session.continuum.estimate_auto(**schema_params("continuum", "estimate_auto"))
    t1 = time.time()
    ident = cont.absorbers.identify_lines(**chain_params("identify_lines"))
    t2 = time.time()
    fit = cont.absorbers.doublets_auto(**chain_params("doublets_auto"))
    t3 = time.time()
    times = dict(continuum=t1 - t0, identify=t2 - t1, doublets=t3 - t2)
    return cont, ident, fit, times


# --- Metrics ------------------------------------------------------------

def _dv(z1, z2):
    return abs(z1 - z2) / (1 + z2) * C_KMS


def continuum_metrics(truth, spec):
    """Ratio of the estimated to the true continuum, by region."""
    x = spec.x.value
    ratio = spec.cont.value / continuum(truth, x)
    rest = x / (1 + truth.z_em)
    out = {}
    for name, (lo, hi) in {"forest": (LYB_NM, LYA_NM),
                           "red": (LYA_NM + 1.0, REST_WINDOW_NM[1])}.items():
        r = ratio[(rest > lo) & (rest < hi) & np.isfinite(ratio)]
        out[name] = (dict(median=float(np.median(r)),
                          p05=float(np.percentile(r, 5)),
                          p95=float(np.percentile(r, 95)))
                     if r.size else None)
    return out


def _candidates(spec):
    meta = spec.meta
    def decode(key):
        value = meta.get(key) or {}
        return json.loads(value) if isinstance(value, str) else dict(value)
    series, z = decode("series_map_json"), decode("z_map_json")
    return [(series[k], float(z[k])) for k in series]


def identification_metrics(truth, systems, spec):
    """Purity and completeness of the metal candidates of identify_lines.

    A candidate is correct when a true system of the same species has a
    component within MATCH_KMS; a second candidate on the same system is
    a duplicate (F49); candidates beyond z_em are counted apart (F55).
    Of the false ones, those on a component of their species in the truth
    but not in ``systems`` (the associated absorbers of the grid
    robustness) are counted apart too (``n_false_associated``, PR 7a.2).
    """
    candidates = _candidates(spec)
    out = {"n_lya": sum(1 for s, _ in candidates if s.startswith("Ly")),
           "systems": []}
    for species in METALS:
        true = [s for s in systems if s["species"] == species]
        detection = [is_detectable(s, truth) for s in true]
        detectable = [d[0] for d in detection]
        cands = [z for s, z in candidates if s == species]
        found, duplicates, false, beyond = set(), 0, 0, 0
        false_associated = 0
        for z in cands:
            beyond += z > truth.z_em
            hits = [i for i, s in enumerate(true)
                    if min(_dv(z, c.z) for c in s["components"]) <= MATCH_KMS]
            if not hits:
                false += 1
                false_associated += any(c.series == species and _dv(z, c.z) <= MATCH_KMS
                                        for c in truth.components)
            elif hits[0] in found:
                duplicates += 1
            else:
                found.add(hits[0])
        out[species] = dict(
            n_true=len(true), n_detectable=int(sum(detectable)),
            n_candidates=len(cands), n_found=len(found),
            n_found_detectable=sum(1 for i in found if detectable[i]),
            n_duplicates=duplicates, n_false=false, n_beyond_z_em=beyond,
            n_false_associated=int(false_associated))
        # Significance of the weaker line and outcome, per system, to
        # check the detection limit
        out["systems"] += [dict(species=species, significance=d[1] / d[2],
                                found=i in found)
                           for i, d in enumerate(detection)]
    return out


# --- The credibility of the hypotheses (F164, PR 7c.4c) ---------------------

#: The bins of the weight of the reliability diagram
WEIGHT_BINS = tuple(np.linspace(0.0, 1.0, 11))
#: True components of one kind closer than this (km/s) are one system, for
#: the truth of the hypotheses
SYSTEM_KMS = 300.0


def _true_systems(truth):
    """The true components grouped into systems: {(kind, key): [components]}
    with the kind the species of a metal or ``"HI"`` (the Lyman series), the
    components of a kind closer than SYSTEM_KMS one system (its key the
    index of its first)."""
    out = {}
    by_kind = {}
    for c in truth.components:
        kind = "HI" if _is_lyman(c.series.split(",")[0]) else c.series
        by_kind.setdefault(kind, []).append(c)
    for kind, comps in by_kind.items():
        comps = sorted(comps, key=lambda c: c.z)
        key = 0
        for i, c in enumerate(comps):
            if i and _dv(c.z, comps[i - 1].z) > SYSTEM_KMS:
                key = i
            out.setdefault((kind, key), []).append(c)
    return out


def _matched_systems(series, redshifts, systems):
    """The true systems that components of ``series`` at ``redshifts``
    match (a true component of the kind within MATCH_KMS)."""
    kind = "HI" if _is_lyman(series) else series
    hits = set()
    for z in redshifts:
        for key, comps in systems.items():
            if key[0] == kind and any(_dv(z, c.z) <= MATCH_KMS for c in comps):
                hits.add(key)
    return hits


def hypothesis_calibration(truth, spec):
    """The hypotheses of the groups compared by the identification by
    hypotheses (those with a weight), each with its weight and whether it
    is right against the truth (F164, PR 7c.4c).

    An edge (a candidate pair) is true when the components of its species
    in the hypothesis it anchors alone match a true system of the species
    (H I for Ly_ab) within MATCH_KMS; a hypothesis is right when none of
    its anchors is false and its anchors match all the true systems that
    the true edges of the group match (the hypothesis of the default labels
    is right when the group has no true edge). A dropped hypothesis has
    weight 0. Empty for the current method."""
    from astrocook.core.identification import META_KEY
    structure = spec.meta.get(META_KEY)
    if structure is None or structure.method != "hypotheses":
        return []
    systems = _true_systems(truth)
    out = []
    for g in structure.groups:
        if len(g.hypotheses) < 2 or not g.edges:
            continue
        edge_hits = {}
        for i, e in enumerate(g.edges):
            alone = [h for h in g.hypotheses if tuple(h.anchors) == (i,)]
            zs = [c.z for h in alone for c in h.components if c.series == e.series]
            edge_hits[i] = _matched_systems(e.series, zs, systems)
        target = set().union(*[hits for hits in edge_hits.values() if hits])
        kind = "metal" if any(not _is_lyman(e.series) for e in g.edges) else "HI"
        for h in g.hypotheses:
            if h.weight is None:
                continue
            false = any(not edge_hits[a] for a in h.anchors)
            hits = set().union(*[edge_hits[a] for a in h.anchors]) if h.anchors else set()
            out.append(dict(weight=float(h.weight), right=bool(not false and hits == target),
                            kind=kind, anchored=bool(h.anchors), status=str(h.status)))
    return out


# --- The sweep of the anchored systems (PR 7c.4d) ------------------------

#: How far from the redshift of a system a true component of a species is
#: there, for the tests of the sweep that were not accepted (km/s): the
#: spread of the components of a system (+-50 km/s) and their b
SWEEP_PRESENT_KMS = 100.0


def _shares_lines(a, b):
    """Two series with a transition in common (``Ly_ab`` and the Lyman
    series of the truth, ``FeII`` and a list of FeII lines)."""
    from astrocook.core.atomic_data import series_lines
    return bool(set(series_lines(a)) & set(series_lines(b)))


def _sweep_test(h, truth):
    """One test of the sweep against the truth."""
    species, rest = h.label.split(" at z = ", 1)
    z_sys = float(rest.split()[0])
    own = [c for c in h.components
           if (species == "HI") == c.series.startswith("Ly_ab") and
           not (species != "HI" and c.series.startswith("Ly_"))]
    series = own[0].series if own else None
    true = [c for c in truth.components
            if series is not None and _shares_lines(series, c.series)]
    accepted = str(h.status) == "accepted"
    return dict(species=species, accepted=accepted, status=str(h.status),
                weight=float(h.weight or 0.0), z=z_sys,
                present=bool(any(_dv(c.z, z_sys) <= SWEEP_PRESENT_KMS for c in true)),
                right=bool(accepted and any(_dv(o.z, c.z) <= MATCH_KMS
                                            for o in own for c in true)),
                replaced_lya=int(h.criteria.get("replaced_lya", 0)),
                replaced_unidentified=int(h.criteria.get("replaced_unidentified", 0)),
                hi_z=[o.z for o in own] if species == "HI" and accepted else [])


def sweep_metrics(truth, systems, spec):
    """The tests of the sweep of the anchored systems (the groups of the
    structure with ``anchored_on``), each against the truth: ``present``
    when a true component that shares a line with its series is within
    SWEEP_PRESENT_KMS of the system, ``right`` when it is accepted and one of
    its components is within MATCH_KMS of such a component; the default
    lines it replaced. And the H I of the metal systems of ``systems`` whose
    Ly-a is in the window and which the identification found (a candidate of
    their species within MATCH_KMS): their H I (the true Lyman component
    nearest to the system, within SWEEP_PRESENT_KMS) found by an accepted
    H I of the sweep within MATCH_KMS. None without the structure."""
    from astrocook.core.identification import META_KEY
    structure = spec.meta.get(META_KEY)
    if structure is None or structure.method != "hypotheses":
        return None
    tests = [_sweep_test(g.hypotheses[0], truth) for g in structure.groups
             if g.anchored_on is not None]
    hi_z = [z for t in tests for z in t.pop("hi_z")]
    cands = _candidates(spec)
    lyman = [c for c in truth.components if _is_lyman(c.series.split(",")[0])]
    hi = dict(n_covered=0, n_found=0)
    for s in systems:
        if LYA_NM * (1 + s["z"]) <= truth.lam_min_nm:
            continue
        if not any(sp == s["species"] and min(_dv(z, c.z) for c in s["components"]) <= MATCH_KMS
                   for sp, z in cands):
            continue
        near = [c for c in lyman if _dv(c.z, s["z"]) <= SWEEP_PRESENT_KMS]
        if not near:
            continue
        c_hi = min(near, key=lambda c: _dv(c.z, s["z"]))
        hi["n_covered"] += 1
        hi["n_found"] += any(_dv(z, c_hi.z) <= MATCH_KMS for z in hi_z)
    return dict(tests=tests, hi=hi)


class _Meta:
    """A stand-in with the ``meta`` of a spectrum, for the metrics."""
    def __init__(self, meta):
        self.meta = meta


def _without_sweep(spec):
    """``spec`` (its ``meta``) without the candidate systems of the sweep
    (the components of the accepted tests, groups with ``anchored_on``)."""
    from astrocook.core.identification import META_KEY
    structure = spec.meta.get(META_KEY)
    out = {str(c.candidate) for g in structure.groups if g.anchored_on is not None
           for h in g.hypotheses for c in h.components if c.candidate is not None}
    meta = dict(spec.meta)
    for key in ("series_map_json", "z_map_json"):
        value = meta.get(key) or {}
        value = json.loads(value) if isinstance(value, str) else dict(value)
        meta[key] = json.dumps({k: v for k, v in value.items() if str(k) not in out})
    return _Meta(meta)


def sweep_summary(rows_by_record):
    """Pool the tests of the sweep: per species the tests, those accepted,
    right, with the species present; the purity of the rarer species (not
    the H I) and of the H I; the H I of the metal systems; the default lines
    replaced (by right and wrong tests); the calibration of the weights
    (right: the species present). None without tests."""
    recs = [r for r in rows_by_record if r]
    if not recs:
        return None
    tests = [t for r in recs for t in r["tests"]]
    by_species = {}
    for t in tests:
        d = by_species.setdefault(t["species"], dict(n=0, accepted=0, right=0, present=0,
                                                     present_accepted=0))
        d["n"] += 1
        d["accepted"] += t["accepted"]
        d["right"] += t["right"]
        d["present"] += t["present"]
        d["present_accepted"] += t["present"] and t["accepted"]

    def purity(sel):
        acc = [t for t in tests if t["accepted"] and sel(t)]
        right = sum(t["right"] for t in acc)
        return dict(n=len(acc), right=right, purity=right / len(acc) if acc else None)
    hi = {k: sum(r["hi"][k] for r in recs) for k in ("n_covered", "n_found")}
    hi["completeness"] = hi["n_found"] / hi["n_covered"] if hi["n_covered"] else None
    replaced = {f"{k}_{ok}": sum(t[f"replaced_{k}"] for t in tests
                                 if t["accepted"] and t["right"] == (ok == "right"))
                for k in ("lya", "unidentified") for ok in ("right", "wrong")}
    return dict(n_tests=len(tests), by_species=by_species,
                rare=purity(lambda t: t["species"] != "HI"),
                hi_tests=purity(lambda t: t["species"] == "HI"),
                hi_of_systems=hi, replaced=replaced,
                calibration=calibration_summary([dict(weight=t["weight"], right=t["present"])
                                                 for t in tests]))


def print_sweep(sw):
    """Print :func:`sweep_summary`."""
    if not sw:
        return
    fmt = (lambda v: "-" if v is None else f"{v:.3f}")
    print(f"sweep: {sw['n_tests']} tests; rarer species accepted {sw['rare']['n']}, "
          f"purity {fmt(sw['rare']['purity'])}; H I accepted {sw['hi_tests']['n']}, purity "
          f"{fmt(sw['hi_tests']['purity'])}; H I of the metal systems "
          f"{sw['hi_of_systems']['n_found']}/{sw['hi_of_systems']['n_covered']} "
          f"({fmt(sw['hi_of_systems']['completeness'])}); default lines replaced "
          f"{sw['replaced']}")
    print("  by species (tests/accepted/right/present/present and accepted): " + ", ".join(
        f"{k} {v['n']}/{v['accepted']}/{v['right']}/{v['present']}/{v['present_accepted']}"
        for k, v in sorted(sw["by_species"].items())))
    c = sw["calibration"]
    if c:
        print(f"  weights: ECE {c['ece']:.3f}, Brier {c['brier']:.3f}")


def calibration_summary(rows):
    """The reliability diagram (per bin of weight: number, mean weight,
    fraction right), its expected calibration error and the Brier score."""
    if not rows:
        return None
    w = np.array([r["weight"] for r in rows])
    y = np.array([r["right"] for r in rows], dtype=float)
    bins = []
    ece = 0.0
    for lo, hi in pairwise(WEIGHT_BINS):
        sel = (w >= lo) & ((w < hi) if hi < 1 else (w <= hi))
        n = int(sel.sum())
        if n:
            conf, acc = float(w[sel].mean()), float(y[sel].mean())
            ece += n / len(w) * abs(acc - conf)
            bins.append(dict(lo=float(lo), hi=float(hi), n=n, weight=conf, right=acc))
        else:
            bins.append(dict(lo=float(lo), hi=float(hi), n=0, weight=None, right=None))
    return dict(n=len(rows), ece=float(ece), brier=float(np.mean((w - y) ** 2)),
                frac_right=float(y.mean()), bins=bins)


def calibration_tables(rows):
    """:func:`calibration_summary` of all the rows, by kind (metal, HI), by
    S/N and by R (rows with ``snr`` and ``resol``)."""
    if not rows:
        return None
    def by(key):
        return {str(k): calibration_summary([c for c in rows if c[key] == k])
                for k in sorted({c[key] for c in rows})}
    return dict(all=calibration_summary(rows), by_kind=by("kind"), by_snr=by("snr"),
                by_resol=by("resol"))


def _flags(c):
    """The flags of the benchmark (the old ones: ``n_flagged``). The absurd
    errors of F52 are the flag ``large_errors`` of the fitter since phase 7,
    PR 7d.2 (one definition: before, dlogN > 1 or db > b computed here; the
    fitter adds the error of z above the width of the line), kept under the
    key ``absurd_errors`` for the comparison."""
    flags = []
    errors = [c.dz, c.dlogN, c.db]
    if any(e is None or not np.isfinite(e) for e in errors):
        flags.append("nan_errors")
    elif "large_errors" in (c.flags or ()):
        flags.append("absurd_errors")
    if abs(c.logN - LOGN_BOUNDS[0]) < 0.02 or abs(c.logN - LOGN_BOUNDS[1]) < 0.02:
        flags.append("logN_at_bound")
    if abs(c.b - B_BOUNDS[0]) < 0.02 or abs(c.b - B_BOUNDS[1]) < 0.02 * B_BOUNDS[1]:
        flags.append("b_at_bound")
    return flags


#: Categories of the matched components, for the pulls: "single" if the
#: true system has one component, "multi" otherwise; "clean" if only one
#: fitted component of the series is within CROWD_KMS of the true one,
#: "crowded" otherwise (more components fitted than resolved, e.g. F49)
PULL_CATEGORIES = ("single_clean", "single_crowded", "multi_clean", "multi_crowded")
CROWD_KMS = 50.0


def pull_category(true_comp, n_in_system, fitted):
    near = sum(1 for f in fitted if f.series == true_comp.series
               and _dv(f.z, true_comp.z) <= CROWD_KMS)
    return (f"{'single' if n_in_system == 1 else 'multi'}_"
            f"{'clean' if near <= 1 else 'crowded'}")


def system_column(system, fitted):
    """Total column density of a true system against the sum of the fitted
    components of its species within MATCH_KMS of its velocity range: a
    measure that does not depend on how the system is decomposed.

    The error of the fitted total propagates the errors of the components
    (as independent); ``None`` if any of them is not finite.
    """
    comps = system["components"]
    zs = [c.z for c in comps]
    lo = min(zs) - MATCH_KMS / C_KMS * (1 + min(zs))
    hi = max(zs) + MATCH_KMS / C_KMS * (1 + max(zs))
    inside = [f for f in fitted if f.series == system["species"] and lo <= f.z <= hi]
    true_logN = float(np.log10(np.sum([10**c.logN for c in comps])))
    out = dict(species=system["species"], n_true=len(comps), n_fitted=len(inside),
               true_logN=true_logN, fit_logN=None, dlogN=None,
               # The flags of 7d.2 of the fitted components (F59: the share
               # of the residual in the saturated or unresolved ones)
               measurement_flags=sorted({f for c in inside for f in (c.flags or ())
                                         if f in MEASUREMENT_FLAGS}))
    if inside:
        n = np.array([10**f.logN for f in inside])
        out["fit_logN"] = float(np.log10(n.sum()))
        errs = [f.dlogN for f in inside]
        if all(e is not None and np.isfinite(e) for e in errs):
            sigma_n = np.sqrt(np.sum((n * np.log(10) * np.array(errs)) ** 2))
            out["dlogN"] = float(sigma_n / (n.sum() * np.log(10)))
    return out


def fit_metrics(truth, systems, systs):
    """Components of doublets_auto against the true metal components.

    True and fitted components of the same series are matched one to one,
    closest first, within MATCH_KMS.
    """
    fitted = [c for c in (systs.components if systs else []) if c.series in METALS]
    true = [c for s in systems for c in s["components"]]
    n_in_system = [len(s["components"]) for s in systems for _ in s["components"]]
    pairs = sorted(((_dv(f.z, t.z), i, j) for i, f in enumerate(fitted)
                    for j, t in enumerate(true) if f.series == t.series),
                   key=lambda p: p[0])
    used_f, used_t, matches = set(), set(), []
    classes = ["spurious"] * len(fitted)
    for dv, i, j in pairs:
        if dv > MATCH_KMS:
            break
        if i not in used_f and j not in used_t:
            used_f.add(i)
            used_t.add(j)
            matches.append((fitted[i], true[j], j))
            classes[i] = _truth_class(fitted[i], true[j])
    pulls = {"z": [], "logN": [], "b": []}
    diffs = {"z": [], "logN": [], "b": []}
    by_category = {k: {"z": [], "logN": [], "b": []} for k in PULL_CATEGORIES}
    by_flag = {k: {"z": [], "logN": [], "b": []} for k in (*MEASUREMENT_FLAGS, "none")}
    for f, t, j in matches:
        category = pull_category(t, n_in_system[j], fitted)
        flag_keys = [k for k in MEASUREMENT_FLAGS if k in (f.flags or ())] or ["none"]
        for par, err in (("z", f.dz), ("logN", f.dlogN), ("b", f.db)):
            d = getattr(f, par) - getattr(t, par)
            if par == "z":
                d = d / (1 + t.z) * C_KMS  # km/s
                err = err / (1 + t.z) * C_KMS if err is not None else None
            diffs[par].append(d)
            if err is not None and np.isfinite(err) and err > 0:
                pulls[par].append(d / err)
                by_category[category][par].append(d / err)
                for k in flag_keys:
                    by_flag[k][par].append(d / err)
    flags = [_flags(c) for c in fitted]
    return dict(
        systems=[system_column(s, fitted) for s in systems],
        n_true=len(true), n_fitted=len(fitted), n_matched=len(matches),
        n_spurious=len(fitted) - len(matches), n_missed=len(true) - len(matches),
        n_flagged=sum(bool(f) for f in flags),
        flags={k: sum(k in f for f in flags)
               for k in ("nan_errors", "absurd_errors", "logN_at_bound", "b_at_bound")},
        # The flags of 7d.2, apart from n_flagged (the old ones), and the
        # pulls of the matched components by flag ("none": without any)
        measurement_flags={k: sum(k in (c.flags or ()) for c in fitted)
                           for k in MEASUREMENT_FLAGS},
        pulls_by_flag=by_flag,
        diffs=diffs, pulls=pulls, pulls_by_category=by_category,
        flag_confusion=flag_confusion(fitted, classes),
        chi2r=[float(c.chi2) for c in fitted if c.chi2 is not None])


# --- Metrics of phase 7 (PR 7a.2) --------------------------------------
#
# On the truth map of the mock (astrocook.mocks.truth_map): what absorbs
# each pixel, as a depth 1 - T per category (agreed with Guido, 16 November
# 2026).

#: Depth (1 - T) above which a pixel of the truth map is absorbed: the full
#: account of the absorption, detectable or not
ABSORBED_DEPTH = 0.01
#: True transmission of the H I Lyman-alpha line of a DLA or sub-DLA that
#: bounds its damped wings: between the two, the saturated core (below the
#: first) excluded
WING_T = (0.1, 0.95)
#: Telluric transmission below which a pixel is on a telluric line
TELLURIC_T = 0.9
#: True transmission above which a forest pixel is nearly unabsorbed (what
#: is left to constrain the continuum there)
FOREST_FREE_T = 0.97
#: Name of the column of the uncertainty of the continuum (fixed here for
#: PR 7b.2): the pull of the continuum is measured when it exists
DCONT_COLUMN = "dcont"
#: Regions of the continuum metrics
CONTINUUM_REGIONS = ("forest", "red", "dla_wings", "subdla_wings", "bal",
                     "emission_peaks", "proximity")
#: Outcomes of an absorbed pixel (pixel_metrics)
PIXEL_OUTCOMES = ("undetected", "unidentified", "lya_default_right",
                  "lya_default_wrong", "right", "wrong")
#: Rows of the pixel table: the category that absorbs a pixel most, and
#: ``none`` for the unabsorbed pixels inside a region
PIXEL_ROWS = (*CATEGORIES, "none")
#: Flags of the confusion table: those of the fit (#517) and those of the
#: benchmark (errors not finite or absurd, F52)
CONFUSION_FLAGS = (*COMPONENT_FLAGS, "absurd_errors")
#: Classes of the fitted components against the truth (confusion table)
TRUTH_CLASSES = ("good", "pull_gt3", "spurious")


def detectable_depth(snr):
    """Depth of a pixel detected at DETECTION_SIGMA over a resolution
    element (PIXELS_PER_FWHM pixels), not below ABSORBED_DEPTH."""
    return max(ABSORBED_DEPTH, DETECTION_SIGMA / (snr * np.sqrt(PIXELS_PER_FWHM)))


def total_depth(tmap):
    """Depth of all the categories together (the product of their
    transmissions)."""
    return 1.0 - np.prod([1.0 - d for d in tmap.depth.values()], axis=0)


def dominant_row(tmap):
    """Index in PIXEL_ROWS of each pixel: the category that absorbs it
    most, or ``none`` below ABSORBED_DEPTH."""
    depth = np.vstack([tmap.depth[c] for c in CATEGORIES])
    return np.where(total_depth(tmap) > ABSORBED_DEPTH, np.argmax(depth, axis=0),
                    len(CATEGORIES))


def _hi_rows(row):
    return np.isin(row, [CATEGORIES.index("forest"), CATEGORIES.index("strong_hi")])


def _strong_hi(truth, logN_range):
    return [c for c in truth.components if component_category(c) == "strong_hi"
            and logN_range[0] <= c.logN < logN_range[1]]


def continuum_regions(truth, tmap):
    """Masks of the regions of the continuum metrics, from the truth:

    - ``forest``, ``red``: as in phase 4 (rest frame of the true z_em);
    - ``dla_wings``, ``subdla_wings``: the damped wings of the strong H I
      systems, where the true transmission of their Lyman-alpha line alone
      is within WING_T (the saturated core excluded);
    - ``bal``: the pixels absorbed by the BAL troughs (ABSORBED_DEPTH);
    - ``emission_peaks``: within half a FWHM of the peak of each emission
      line (one FWHM in all), at its shift;
    - ``proximity``: from ``proximity_kms`` to the blue of the Lyman-alpha
      emission to the emission itself.
    """
    x = tmap.x_nm
    rest = x / (1 + truth.z_em)
    out = {"forest": (rest > LYB_NM) & (rest < LYA_NM),
           "red": (rest > LYA_NM + 1.0) & (rest < REST_WINDOW_NM[1])}
    for name, logN_range in (("dla_wings", (DLA_LOGN[0], np.inf)),
                             ("subdla_wings", SUBDLA_LOGN)):
        mask = np.zeros(len(x), dtype=bool)
        for c in _strong_hi(truth, logN_range):
            t = np.exp(-optical_depth(x * 10.0, "Ly_a", c.z, c.logN, c.b, c.btur))
            mask |= (t >= WING_T[0]) & (t <= WING_T[1])
        out[name] = mask
    out["bal"] = tmap.absorbed("bal", ABSORBED_DEPTH)
    mask = np.zeros(len(x), dtype=bool)
    for line in truth.emission_lines:
        dv = np.log(rest / line.rest_nm) * C_KMS - line.shift_kms
        mask |= np.abs(dv) <= line.fwhm_kms / 2
    out["emission_peaks"] = mask
    mask = np.zeros(len(x), dtype=bool)
    if truth.proximity_kms:
        dv = np.log(rest / LYA_NM) * C_KMS
        mask = (dv >= -truth.proximity_kms) & (dv <= 0)
    out["proximity"] = mask
    return out


def region_continuum_metrics(truth, spec, regions):
    """The continuum estimated against the true one in each region: the
    ratio (median, p05, p95), the integrated deviation (signed, in units of
    the true continuum integrated over the region: a bias of the
    equivalent widths; and absolute), and the pull against DCONT_COLUMN
    when the spectrum has it (None otherwise)."""
    x = spec.x.value
    true = continuum(truth, x)
    est = spec.cont.value
    width = x  # pixels uniform in velocity: their width goes as x
    dcont = (np.asarray(spec.get_column(DCONT_COLUMN).value, dtype=float)
             if spec.has_aux_column(DCONT_COLUMN) else None)
    out = {}
    for name, mask in regions.items():
        sel = mask & np.isfinite(est) & np.isfinite(true)
        if not sel.any():
            out[name] = None
            continue
        r = est[sel] / true[sel]
        norm = np.sum(width[sel] * true[sel])
        diff = width[sel] * (est[sel] - true[sel])
        pull = None
        if dcont is not None:
            ok = sel & np.isfinite(dcont) & (dcont > 0)
            pull = _stats((est[ok] - true[ok]) / dcont[ok])
        out[name] = dict(n=int(sel.sum()), median=float(np.median(r)),
                         p05=float(np.percentile(r, 5)),
                         p95=float(np.percentile(r, 95)),
                         dev=float(np.sum(diff) / norm),
                         abs_dev=float(np.sum(np.abs(diff)) / norm), pull=pull)
    return out


def forest_truth(truth, tmap):
    """What the truth leaves to constrain the continuum in the forest: the
    mean true transmission and the fraction of pixels with T above
    FOREST_FREE_T."""
    rest = tmap.x_nm / (1 + truth.z_em)
    t = 1.0 - total_depth(tmap)[(rest > LYB_NM) & (rest < LYA_NM)]
    return dict(mean_t=float(np.mean(t)), frac_free=float(np.mean(t > FOREST_FREE_T)))


def z_em_estimate(session):
    """The z_em estimated by the chain, if it has an estimator (F153, PR
    7b.3); None until then."""
    return None


def _decode_meta(meta, key):
    value = meta.get(key) or {}
    return json.loads(value) if isinstance(value, str) else dict(value)


def candidate_list(spec):
    """The candidates of identify_lines: series, z, region and score."""
    meta = spec.meta
    maps = {k: _decode_meta(meta, f"{k}_map_json")
            for k in ("series", "z", "region", "score")}
    return [dict(series=s, z=float(maps["z"][k]),
                 region=int(float(maps["region"].get(k, 0) or 0)),
                 score=float(maps["score"].get(k, np.nan)))
            for k, s in maps["series"].items()]


def _is_lyman(series):
    return series.startswith("Ly")


def _is_default_lya(cand):
    """The Lyman-alpha of F51: assigned to a region without a test (score
    0)."""
    return cand["series"] == "Ly_a" and cand["score"] == 0.0


def _matches_truth(cand, truth):
    """A metal candidate with a true component of its series within
    MATCH_KMS (associated absorbers included)."""
    return any(c.series == cand["series"] and _dv(cand["z"], c.z) <= MATCH_KMS
               for c in truth.components)


def _line_pixels(cand, x):
    """Pixels of the lines of a candidate (nearest pixel of each line in
    the spectrum)."""
    out = []
    for t in transitions(cand["series"]):
        lam = ATOM_DATA[t]["wave"] / 10.0 * (1 + cand["z"])
        if x[0] <= lam <= x[-1]:
            out.append(int(np.argmin(np.abs(x - lam))))
    return out


def telluric_false(truth, tmap, candidates):
    """False candidates on tellurics, by series: no true absorber of the
    series within MATCH_KMS (for H I: none of its pixels in an H I row of
    :func:`dominant_row`), and a line on a pixel of telluric transmission
    below TELLURIC_T. Counted apart from (and also in) the false ones of
    identification_metrics."""
    out = {}
    if not tmap.depth["telluric"].any():
        return out
    tell = tmap.depth["telluric"] > 1.0 - TELLURIC_T
    hi = _hi_rows(dominant_row(tmap))
    for cand in candidates:
        pixels = _line_pixels(cand, tmap.x_nm)
        if _is_lyman(cand["series"]):
            false = not any(hi[i] for i in pixels)
        else:
            false = not _matches_truth(cand, truth)
        if false and any(tell[i] for i in pixels):
            out[cand["series"]] = out.get(cand["series"], 0) + 1
    return out


def pixel_metrics(truth, tmap, spec, candidates):
    """The absorbed pixels of the truth map against the regions and the
    identifications of identify_lines.

    Each pixel is in the row of the category that absorbs it most (or
    ``none``, unabsorbed, counted only inside a region), and has one
    outcome: ``undetected`` (in no region), ``unidentified`` (in a region
    without candidates), ``lya_default_right`` or ``lya_default_wrong``
    (the region has only the Lyman-alpha of F51: right if the pixel is in
    the row ``forest`` or ``strong_hi``), ``right`` (a metal candidate of
    the region matches a true component of its series and the pixel is in
    the row ``metal``, or an H I candidate tested by the identification
    and the pixel is in an H I row) or ``wrong``. The row decides, not any
    absorption above ABSORBED_DEPTH: in the forest nearly every pixel has
    some H I, also those of a BAL trough or of a metal line. Two tables: the pixels with total depth above
    ABSORBED_DEPTH (``absorbed``) and above :func:`detectable_depth`
    (``detectable``).
    """
    n = len(tmap.x_nm)
    ids = (np.asarray(spec.get_column("abs_ids").value).astype(int)
           if spec.has_aux_column("abs_ids") else np.zeros(n, dtype=int))
    total = total_depth(tmap)
    row = dominant_row(tmap)
    metal_px = row == CATEGORIES.index("metal")
    hi_px = _hi_rows(row)
    by_region = {}
    for cand in candidates:
        by_region.setdefault(cand["region"], []).append(cand)
    outcome = np.zeros(n, dtype=int)  # undetected
    for rid in np.unique(ids[ids > 0]):
        inside = ids == rid
        cands = by_region.get(int(rid), [])
        if not cands:
            outcome[inside] = PIXEL_OUTCOMES.index("unidentified")
            continue
        if all(_is_default_lya(c) for c in cands):
            outcome[inside] = np.where(hi_px[inside],
                                       PIXEL_OUTCOMES.index("lya_default_right"),
                                       PIXEL_OUTCOMES.index("lya_default_wrong"))
            continue
        metal_ok = any(not _is_lyman(c["series"]) and _matches_truth(c, truth)
                       for c in cands)
        hi_ok = any(_is_lyman(c["series"]) and not _is_default_lya(c) for c in cands)
        right = (metal_ok & metal_px[inside]) | (hi_ok & hi_px[inside])
        outcome[inside] = np.where(right, PIXEL_OUTCOMES.index("right"),
                                   PIXEL_OUTCOMES.index("wrong"))
    unabsorbed_in_region = (row == len(CATEGORIES)) & (ids > 0)
    out = {}
    for name, depth_min in (("absorbed", ABSORBED_DEPTH),
                            ("detectable", detectable_depth(truth.snr))):
        sel = (total > depth_min) | unabsorbed_in_region
        table = {}
        for r, row_name in enumerate(PIXEL_ROWS):
            in_row = sel & (row == r)
            if row_name == "none":
                in_row &= ids > 0
            table[row_name] = {o: int(np.sum(in_row & (outcome == k)))
                               for k, o in enumerate(PIXEL_OUTCOMES)}
        out[name] = table
    return out


def _truth_class(f, t):
    """``good`` or ``pull_gt3`` for a matched component (any pull above 3,
    or an error not finite or not positive)."""
    for par, err in (("z", f.dz), ("logN", f.dlogN), ("b", f.db)):
        if err is None or not np.isfinite(err) or err <= 0:
            return "pull_gt3"
        if abs(getattr(f, par) - getattr(t, par)) / err > 3:
            return "pull_gt3"
    return "good"


def flag_confusion(fitted, classes):
    """Counts of the fitted components by flag (CONFUSION_FLAGS, plus
    ``any`` and ``none``) and class against the truth (TRUTH_CLASSES)."""
    table = {k: dict.fromkeys(TRUTH_CLASSES, 0) for k in (*CONFUSION_FLAGS, "any", "none")}
    for f, cls in zip(fitted, classes, strict=True):
        # The flags of 7d.2 are counted apart (measurement_flags)
        flags = (set(f.flags or ()) - set(MEASUREMENT_FLAGS)) | ({"absurd_errors"} & set(_flags(f)))
        flags |= {"nan_errors"} & set(_flags(f))
        for k in flags & set(CONFUSION_FLAGS):
            table[k][cls] += 1
        table["any" if flags else "none"][cls] += 1
    return table


def run_point(point):
    """Mock, chain and metrics of one grid point (a record).

    A point of the grid robustness has a fifth element, the case: its
    record adds ``case``, ``associated``, ``z_em_used`` (the z_em given to
    the recipes) and ``dz_em_kms`` (its difference from the true one).
    """
    snr, resol, z_em, seed = point[:4]
    logging.disable(logging.WARNING)
    t0 = time.time()
    extra = {}
    if len(point) > 4:
        truth, systems, added = make_robust_truth(*point)
        strong = added["strong"]
        z_used, dv = z_em_perturbed(snr, resol, z_em, seed)
        extra = dict(case=point[4], associated=added["associated"],
                     z_em_used=z_used, dz_em_kms=dv)
    else:
        truth, systems, strong = make_truth(snr, resol, z_em, seed)
    # The true continuum is not in the spectrum, as for an observed one
    spec = generate(truth, truth_columns=False)
    if extra:
        spec = spec.with_properties(z_em=extra["z_em_used"])
    t_mock = time.time() - t0
    cont, ident, fit, times = run_chain(spec)
    t0 = time.time()
    tmap = truth_map(truth)
    regions = continuum_regions(truth, tmap)
    candidates = candidate_list(ident.spec)
    z_est = z_em_estimate(fit)
    phase7 = dict(
        continuum=region_continuum_metrics(truth, cont.spec, regions),
        forest_truth=forest_truth(truth, tmap),
        dz_em_est_kms=(None if z_est is None else
                       float((z_est - z_em) / (1 + z_em) * C_KMS)),
        telluric_false=telluric_false(truth, tmap, candidates),
        pixels=pixel_metrics(truth, tmap, ident.spec, candidates))
    phase7["calibration"] = hypothesis_calibration(truth, ident.spec)
    phase7["sweep"] = sweep_metrics(truth, systems, ident.spec)
    if phase7["sweep"] is not None:
        # The identification without the candidates the sweep added: the
        # first pass is the same with and without it (PR 7c.4d)
        phase7["identification_first"] = identification_metrics(
            truth, systems, _without_sweep(ident.spec))
    times["metrics"] = time.time() - t0
    return dict(
        snr=snr, resol=resol, z_em=z_em, seed=seed, strong=strong, **extra,
        n_pixels=len(cont.spec.x), times=dict(mock=t_mock, **times),
        continuum=continuum_metrics(truth, cont.spec),
        identification=identification_metrics(truth, systems, ident.spec),
        fit=fit_metrics(truth, systems, fit.systs), phase7=phase7)


# --- Summary ------------------------------------------------------------

def _stats(values):
    """Median, robust sigma (half the 16-84% interval: a few components
    with errors of ~0 give pulls of 1e6 or more, which would dominate a
    standard deviation) and the fraction beyond 3 (for the pulls)."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if not v.size:
        return None
    p16, p84 = np.percentile(v, [16, 84])
    return dict(n=int(v.size), median=float(np.median(v)),
                sigma=float((p84 - p16) / 2), p16=float(p16), p84=float(p84),
                frac_gt3=float(np.mean(np.abs(v) > 3)))


def summarize(records, key=None):
    """Pool the records (all, or grouped by one grid parameter)."""
    if key is not None:
        groups = sorted({r[key] for r in records})
        return {str(g): summarize([r for r in records if r[key] == g])
                for g in groups}
    out = {"n_mocks": len(records)}
    for region in ("forest", "red"):
        meds = [r["continuum"][region]["median"] for r in records
                if r["continuum"][region]]
        p05 = [r["continuum"][region]["p05"] for r in records if r["continuum"][region]]
        p95 = [r["continuum"][region]["p95"] for r in records if r["continuum"][region]]
        out[f"continuum_{region}"] = dict(
            median_of_medians=float(np.median(meds)) if meds else None,
            median_p05=float(np.median(p05)) if p05 else None,
            median_p95=float(np.median(p95)) if p95 else None)
    ident = {}
    for species in METALS:
        tot = {k: sum(r["identification"][species][k] for r in records)
               for k in records[0]["identification"][species]}
        unique = tot["n_found"] + tot["n_false"]
        tot["purity"] = tot["n_found"] / unique if unique else None
        tot["completeness"] = (tot["n_found_detectable"] / tot["n_detectable"]
                               if tot["n_detectable"] else None)
        ident[species] = tot
    out["identification"] = ident
    fit = {k: sum(r["fit"][k] for r in records)
           for k in ("n_true", "n_fitted", "n_matched", "n_spurious", "n_missed",
                     "n_flagged")}
    fit["flags"] = {k: sum(r["fit"]["flags"][k] for r in records)
                    for k in records[0]["fit"]["flags"]}
    if "measurement_flags" in records[0]["fit"]:
        fit["measurement_flags"] = {k: sum(r["fit"]["measurement_flags"][k] for r in records)
                                    for k in records[0]["fit"]["measurement_flags"]}
        fit["pulls_by_flag"] = {
            k: {par: _stats([p for r in records for p in r["fit"]["pulls_by_flag"][k][par]])
                for par in ("z", "logN", "b")}
            for k in records[0]["fit"]["pulls_by_flag"]}
    for par in ("z", "logN", "b"):
        fit[f"diff_{par}"] = _stats([d for r in records for d in r["fit"]["diffs"][par]])
        fit[f"pull_{par}"] = _stats([p for r in records for p in r["fit"]["pulls"][par]])
    fit["chi2r"] = _stats([c for r in records for c in r["fit"]["chi2r"]])
    columns = [c for r in records for c in r["fit"]["systems"] if c["fit_logN"] is not None]
    fit["system_column"] = dict(
        n_systems=sum(len(r["fit"]["systems"]) for r in records),
        n_fitted=len(columns),
        diff=_stats([c["fit_logN"] - c["true_logN"] for c in columns]),
        pull=_stats([(c["fit_logN"] - c["true_logN"]) / c["dlogN"] for c in columns
                     if c["dlogN"]]),
        diff_multi=_stats([c["fit_logN"] - c["true_logN"] for c in columns
                           if c["n_true"] > 1]),
        pull_multi=_stats([(c["fit_logN"] - c["true_logN"]) / c["dlogN"]
                           for c in columns if c["dlogN"] and c["n_true"] > 1]))
    fit["pulls_by_category"] = {
        cat: {par: _stats([p for r in records
                           for p in r["fit"]["pulls_by_category"][cat][par]])
              for par in ("z", "logN", "b")}
        for cat in PULL_CATEGORIES}
    out["fit"] = fit
    out["completeness_by_significance"] = completeness_by_significance(records)
    out["completeness_by_significance_species"] = {
        species: completeness_by_significance(records, species) for species in METALS}
    if all("phase7" in r for r in records):
        out["phase7"] = summarize_phase7(records)
    out["seconds"] = {k: float(np.sum([r["times"][k] for r in records]))
                      for k in records[0]["times"]}
    return out


def _median(values):
    v = [x for x in values if x is not None]
    return float(np.median(v)) if v else None


def summarize_phase7(records):
    """Pool the metrics of phase 7 (PR 7a.2): medians over the mocks for
    the continuum and the forest, sums for the counts."""
    out = {}
    cont = {}
    for region in CONTINUUM_REGIONS:
        rs = [r["phase7"]["continuum"].get(region) for r in records]
        rs = [m for m in rs if m]
        pulls = [m["pull"] for m in rs if m["pull"]]
        cont[region] = dict(
            n_mocks=len(rs),
            **{f"median_{k}": _median([m[k] for m in rs])
               for k in ("median", "p05", "p95", "dev", "abs_dev")},
            pull=(dict(median=_median([p["median"] for p in pulls]),
                       sigma=_median([p["sigma"] for p in pulls]),
                       frac_gt3=_median([p["frac_gt3"] for p in pulls]))
                  if pulls else None))
    out["continuum"] = cont
    out["forest_truth"] = {k: _median([r["phase7"]["forest_truth"][k] for r in records])
                           for k in ("mean_t", "frac_free")}
    dz = [r["phase7"]["dz_em_est_kms"] for r in records
          if r["phase7"]["dz_em_est_kms"] is not None]
    out["z_em_estimate"] = (dict(**_stats(dz), frac_gt500=float(np.mean(np.abs(dz) > 500)))
                            if dz else None)
    tell = {}
    for r in records:
        for k, v in r["phase7"]["telluric_false"].items():
            tell[k] = tell.get(k, 0) + v
    out["telluric_false"] = tell
    pixels = {}
    for name in ("absorbed", "detectable"):
        table = {row: {o: sum(r["phase7"]["pixels"][name][row][o] for r in records)
                       for o in PIXEL_OUTCOMES} for row in PIXEL_ROWS}
        fractions = {}
        for row, counts in table.items():
            n = sum(counts.values())
            fractions[row] = ({o: c / n for o, c in counts.items()} if n else None)
        pixels[name] = dict(counts=table, fractions=fractions)
    out["pixels"] = pixels
    conf = {k: {c: sum(r["fit"]["flag_confusion"][k][c] for r in records)
                for c in TRUTH_CLASSES}
            for k in records[0]["fit"]["flag_confusion"]}
    n_bad = sum(conf[k][c] for k in ("any", "none") for c in ("pull_gt3", "spurious"))
    rates = {}
    for k, counts in conf.items():
        n = sum(counts.values())
        bad = counts["pull_gt3"] + counts["spurious"]
        rates[k] = dict(n=n, frac_bad=bad / n if n else None,
                        frac_of_bad=bad / n_bad if n_bad else None)
    out["flag_confusion"] = dict(counts=conf, rates=rates)
    cal = calibration_tables([dict(c, snr=r["snr"], resol=r["resol"]) for r in records
                              for c in r["phase7"].get("calibration", [])])
    if cal:
        out["calibration"] = cal
    sw = sweep_summary([r["phase7"].get("sweep") for r in records])
    if sw:
        out["sweep"] = sw
    return out


#: Bins of the significance of the weaker line (EW / sigma_EW)
SIGNIFICANCE_BINS = (0.0, 1.0, 2.0, 3.0, 5.0, 10.0, np.inf)


def completeness_by_significance(records, species=None):
    """Fraction of the true systems found, in bins of the significance of
    their weaker line: it should rise around DETECTION_SIGMA. All the
    species, or one (F57, by species since PR 7c.1)."""
    systems = [s for r in records for s in r["identification"]["systems"]
               if species is None or s["species"] == species]
    out = {}
    for lo, hi in pairwise(SIGNIFICANCE_BINS):
        sel = [s["found"] for s in systems if lo <= s["significance"] < hi]
        out[f"{lo:g}-{hi:g}"] = dict(n=len(sel),
                                     found=float(np.mean(sel)) if sel else None)
    return out


def grid_points(grid):
    g = GRIDS[grid]
    points = [(snr, resol, z_em, seed) for snr in g["snr"] for resol in g["resol"]
              for z_em in g["z_em"] for seed in g["seeds"]]
    if "cases" in g:
        points = [(*p, case) for case in g["cases"] for p in points]
    return points


def run(grid="full", workers=1):
    points = grid_points(grid)
    if workers > 1:
        with Pool(workers, initializer=set_chain, initargs=(dict(CHAIN),)) as pool:
            records = pool.map(run_point, points, chunksize=1)
    else:
        records = [run_point(p) for p in points]
    summary = dict(all=summarize(records), by_snr=summarize(records, "snr"),
                   by_resol=summarize(records, "resol"),
                   by_z_em=summarize(records, "z_em"))
    if "cases" in GRIDS[grid]:
        summary["by_case"] = summarize(records, "case")
    return dict(grid=grid, match_kms=MATCH_KMS, records=records, summary=summary)


def print_summary(result):
    s = result["summary"]["all"]
    print(f"{s['n_mocks']} mocks")
    for region in ("forest", "red"):
        c = s[f"continuum_{region}"]
        print(f"continuum {region}: median ratio {c['median_of_medians']:.3f} "
              f"(p05 {c['median_p05']:.3f}, p95 {c['median_p95']:.3f})")
    for species, t in s["identification"].items():
        print(f"{species}: purity {t['purity']}, completeness {t['completeness']}, "
              f"duplicates {t['n_duplicates']}, beyond z_em {t['n_beyond_z_em']}")
    f = s["fit"]
    print(f"fit: {f['n_matched']}/{f['n_true']} matched, {f['n_spurious']} spurious, "
          f"{f['n_flagged']} flagged {f['flags']}")
    sc = f["system_column"]
    if sc["diff"] and sc["pull"]:
        print(f"system column ({sc['n_fitted']}/{sc['n_systems']}): diff median "
              f"{sc['diff']['median']:.3f} sigma {sc['diff']['sigma']:.3f}; pull median "
              f"{sc['pull']['median']:.2f} sigma {sc['pull']['sigma']:.2f}")
    for cat, pars in f["pulls_by_category"].items():
        text = ", ".join(f"{par} {p['sigma']:.2f}" for par, p in pars.items() if p)
        n = pars["logN"]["n"] if pars["logN"] else 0
        print(f"  pull sigma, {cat} ({n}): {text}")
    print("completeness by significance: " + ", ".join(
        f"{k}: {v['found']:.2f} ({v['n']})" for k, v in
        s["completeness_by_significance"].items() if v["found"] is not None))
    for par in ("z", "logN", "b"):
        p = f[f"pull_{par}"]
        if p:
            print(f"  pull {par}: median {p['median']:.2f}, sigma {p['sigma']:.2f}, "
                  f"|pull|>3 {p['frac_gt3']:.2f}")


#: The baseline (version 2.2.0, after phase 4): the pooled summary of the
#: full grid, kept as history since PR 7c.1 (``--compare-baseline``; the
#: reference below is that of ``--compare``). That of 2.1.0, before the
#: fixes of phase 4, is kept in validate_mocks_baseline_2.1.0.json
BASELINE = ROOT / "benchmarks" / "results" / "validate_mocks_baseline.json"


def key_numbers(summary):
    """The numbers of a pooled summary compared against the baseline."""
    out = {}
    for region in ("forest", "red"):
        out[f"continuum_{region}"] = summary[f"continuum_{region}"]["median_of_medians"]
    for species, t in summary["identification"].items():
        out[f"{species}_purity"] = t["purity"]
        out[f"{species}_completeness"] = t["completeness"]
        out[f"{species}_duplicates"] = t["n_duplicates"]
    f = summary["fit"]
    out["fit_matched"] = f["n_matched"]
    out["fit_spurious"] = f["n_spurious"]
    out["fit_flagged"] = f["n_flagged"]
    for par in ("z", "logN", "b"):
        out[f"pull_{par}_sigma"] = f[f"pull_{par}"]["sigma"] if f[f"pull_{par}"] else None
    sc = f["system_column"]
    out["column_diff_sigma"] = sc["diff"]["sigma"] if sc["diff"] else None
    out["column_pull_sigma"] = sc["pull"]["sigma"] if sc["pull"] else None
    return out


def compare(summary, baseline):
    """Print the key numbers of a run next to those of the baseline."""
    now, before = key_numbers(summary), key_numbers(baseline)
    print(f"{'':28s}{'baseline':>12s}{'now':>12s}")
    for key, value in now.items():
        old = before.get(key)
        fmt = (lambda v: "-" if v is None else
               (f"{v:.3f}" if isinstance(v, float) else str(v)))
        mark = "" if old == value else "  *"
        print(f"{key:28s}{fmt(old):>12s}{fmt(value):>12s}{mark}")


#: The reference of phase 7 (PR 7a.2, regenerated in PR 7c.1, after F68):
#: the summaries of the current algorithms on the grids robustness and
#: full, each under its name; the reference of ``--compare`` for both grids
#: (the baseline of 2.2.0 above is history)
REFERENCE = ROOT / "benchmarks" / "results" / "validate_mocks_robustness_baseline.json"


def pixel_fractions(pixels):
    """Fraction of each outcome over the absorbed pixels of every
    category (the row ``none`` left out)."""
    counts = {o: sum(pixels["counts"][row][o] for row in CATEGORIES)
              for o in PIXEL_OUTCOMES}
    n = sum(counts.values())
    return {o: c / n if n else None for o, c in counts.items()}


def key_numbers_phase7(summary):
    """The numbers of the metrics of phase 7 compared against the
    reference."""
    p7 = summary["phase7"]
    out = {}
    for region, c in p7["continuum"].items():
        out[f"cont_{region}_median"] = c["median_median"]
        out[f"cont_{region}_dev"] = c["median_dev"]
    out["forest_mean_t"] = p7["forest_truth"]["mean_t"]
    out["forest_frac_free"] = p7["forest_truth"]["frac_free"]
    z = p7["z_em_estimate"]
    out["z_em_est_sigma"] = z["sigma"] if z else None
    out["telluric_false"] = sum(p7["telluric_false"].values())
    for name in ("absorbed", "detectable"):
        for o, f in pixel_fractions(p7["pixels"][name]).items():
            out[f"pix_{name[:6]}_{o}"] = f
    for k in ("any", "none"):
        out[f"flags_{k}_frac_bad"] = p7["flag_confusion"]["rates"][k]["frac_bad"]
    return out


def compare_phase7(summary, reference):
    """Print the numbers of phase 7 of a run next to the reference."""
    now, before = key_numbers_phase7(summary), key_numbers_phase7(reference)
    fmt = (lambda v: "-" if v is None else
           (f"{v:.3f}" if isinstance(v, float) else str(v)))
    print(f"{'':32s}{'reference':>12s}{'now':>12s}")
    for key, value in now.items():
        old = before.get(key)
        mark = "" if old == value else "  *"
        print(f"{key:32s}{fmt(old):>12s}{fmt(value):>12s}{mark}")


def print_calibration(cal):
    """Print the calibration of the credibility (:func:`calibration_tables`)."""
    if not cal:
        return
    a = cal["all"]
    print(f"credibility ({a['n']} hypotheses, {a['frac_right']:.2f} right): ECE "
          f"{a['ece']:.3f}, Brier {a['brier']:.3f}; right by weight: " + ", ".join(
              f"{b['lo']:.1f}-{b['hi']:.1f} {b['right']:.2f} ({b['n']})"
              for b in a["bins"] if b["n"]))
    for name in ("by_kind", "by_snr", "by_resol"):
        print(f"  {name}: " + ", ".join(
            f"{k} ECE {v['ece']:.3f} Brier {v['brier']:.3f} ({v['n']})"
            for k, v in cal[name].items() if v))


def print_phase7(summary):
    """Print the key numbers of the metrics of phase 7."""
    p7 = summary.get("phase7")
    if not p7:
        return
    for region, c in p7["continuum"].items():
        if c["n_mocks"]:
            print(f"continuum {region} ({c['n_mocks']} mocks): median ratio "
                  f"{c['median_median']:.3f} (p05 {c['median_p05']:.3f}, p95 "
                  f"{c['median_p95']:.3f}), deviation {c['median_dev']:+.3f} "
                  f"(abs {c['median_abs_dev']:.3f})")
    f = p7["forest_truth"]
    print(f"forest: mean T {f['mean_t']:.3f}, T > {FOREST_FREE_T}: {f['frac_free']:.3f}")
    print(f"false candidates on tellurics: {p7['telluric_false']}")
    for name in ("absorbed", "detectable"):
        print(f"pixels {name}: " + ", ".join(
            f"{o} {v:.3f}" for o, v in pixel_fractions(p7["pixels"][name]).items()
            if v is not None))
    print_calibration(p7.get("calibration"))
    print_sweep(p7.get("sweep"))
    rates = p7["flag_confusion"]["rates"]
    print("flags (n, bad fraction, fraction of the bad): " + ", ".join(
        f"{k} ({r['n']}, {r['frac_bad']:.2f}, {r['frac_of_bad']:.2f})"
        for k, r in rates.items() if r["n"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--grid", choices=sorted(GRIDS), default="full")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--json", help="write the records and the summary here")
    parser.add_argument("--save-baseline", action="store_true",
                        help="write the summary of a full run as the baseline")
    parser.add_argument("--compare", action="store_true",
                        help="compare the summary of a full or robustness run, and "
                             "its metrics of phase 7, with the reference")
    parser.add_argument("--compare-baseline", action="store_true",
                        help="compare the summary of a full run with the baseline "
                             "of 2.2.0 (history, PR 7c.1)")
    parser.add_argument("--save-reference", action="store_true",
                        help="write the summary of a full or robustness run in the "
                             "reference of phase 7")
    parser.add_argument("--method", choices=("doublets", "hypotheses"), default=None,
                        help="the method of identify_lines and doublets_auto (PR 7c.4c)")
    parser.add_argument("--multiplets", default=None,
                        help="the species of identify_lines (e.g. 'all,Ly_ab')")
    parser.add_argument("--doublets-multiplets", default=None,
                        help="the species of doublets_auto (e.g. 'all')")
    parser.add_argument("--prefilter", type=float, default=None,
                        help="the pre-filter of identify_lines (method hypotheses)")
    parser.add_argument("--screen", type=float, default=None,
                        help="the screen of identify_lines (method hypotheses)")
    parser.add_argument("--criterion", default=None,
                        help="the criterion of identify_lines (method hypotheses)")
    parser.add_argument("--temper", action="store_true", default=None,
                        help="the weights tempered (method hypotheses)")
    parser.add_argument("--max-added", type=int, default=None,
                        help="the most components added by the refinement of each "
                             "hypothesis (method hypotheses; 0: no refinement; -1: no "
                             "fixed number; the default 3)")
    parser.add_argument("--no-split-gaps", dest="split_gaps", action="store_false",
                        default=None,
                        help="the regions not split at their gaps (method hypotheses, "
                             "as PR 7c.4c)")
    parser.add_argument("--sweep", choices=["HI", "all", "none"], default=None,
                        help="the species of the sweep of the anchored systems (method "
                             "hypotheses; default that of identify_lines, HI; 'none' as "
                             "PR 7c.4c-bis)")
    args = parser.parse_args()
    set_chain(dict(method=args.method, multiplets=args.multiplets, max_added=args.max_added,
                   sweep=args.sweep,
                   split_gaps=args.split_gaps,
                   doublets_multiplets=args.doublets_multiplets, prefilter=args.prefilter,
                   screen=args.screen, criterion=args.criterion, temper=args.temper))
    t0 = time.time()
    result = run(args.grid, args.workers)
    result["chain"] = dict(CHAIN)
    result["seconds"] = time.time() - t0
    print_summary(result)
    print_phase7(result["summary"]["all"])
    print(f"{result['seconds']:.0f} s")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(result, f, indent=1)
    if args.compare_baseline and args.grid == "full":
        with open(BASELINE) as f:
            compare(result["summary"]["all"], json.load(f)["summary"]["all"])
    reference = json.loads(REFERENCE.read_text()) if REFERENCE.exists() else {}
    if args.compare and args.grid in reference:
        compare(result["summary"]["all"], reference[args.grid]["summary"]["all"])
        compare_phase7(result["summary"]["all"], reference[args.grid]["summary"]["all"])
    if args.save_reference and args.grid in ("full", "robustness"):
        from astrocook import __version__
        reference[args.grid] = dict(astrocook_version=__version__, match_kms=MATCH_KMS,
                                    summary=result["summary"])
        REFERENCE.write_text(json.dumps(reference, indent=1))
    if args.save_baseline and args.grid == "full":
        from astrocook import __version__
        BASELINE.parent.mkdir(exist_ok=True)
        with open(BASELINE, "w") as f:
            json.dump(dict(astrocook_version=__version__, grid=args.grid,
                           match_kms=MATCH_KMS, summary=result["summary"]),
                      f, indent=1)


if __name__ == "__main__":
    main()
