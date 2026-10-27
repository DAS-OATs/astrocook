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

    python benchmarks/validate_mocks.py [--grid full|reduced|tiny]
        [--workers N] [--json out.json]

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
"""
import argparse
import json
import logging
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS  # noqa: E402
from astrocook.mocks import (  # noqa: E402
    DLA_LOGN, LYMAN_SERIES, QSO_EMISSION_LINES, SUBDLA_LOGN, MockTruth,
    continuum, generate, lya_forest, metal_systems, optical_depth,
    strong_hi_systems,
)
from astrocook.mocks.physics import C_KMS  # noqa: E402

GRIDS = {
    "full": dict(snr=(10.0, 30.0, 100.0), resol=(8900.0, 45000.0, 140000.0),
                 z_em=(2.5, 3.5, 4.5), seeds=tuple(range(5))),
    "reduced": dict(snr=(10.0, 30.0), resol=(8900.0, 45000.0), z_em=(2.5, 3.5),
                    seeds=(0, 2)),
    "tiny": dict(snr=(30.0,), resol=(8900.0,), z_em=(2.5,), seeds=(2,)),
}

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
#: Bounds of the fits, for the flags (VoigtModelConstraintV2.DEFAULT_BOUNDS
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


def make_truth(snr, resol, z_em, seed):
    """Mock of one grid point, and its metal systems.

    Returns
    -------
    truth : MockTruth
    systems : list of dict
        ``{"species", "z", "components": [MockComponent, ...]}``, one per
        metal system (without its H I).
    """
    rng = np.random.default_rng([int(seed), int(snr), int(resol), int(z_em * 10)])
    lam_min, lam_max = (w * (1 + z_em) for w in REST_WINDOW_NM)
    dv = C_KMS / resol / PIXELS_PER_FWHM
    forest = lya_forest(rng, lam_min / LYA_NM - 1, z_em, series="Ly_ab",
                        strong_series=LYMAN_SERIES)
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
    strong = None
    if seed in (0, 1):
        # Its Lyman-alpha in the forest of the window, away from the edges
        strong = "DLA" if seed == 0 else "subDLA"
        comps += strong_hi_systems(rng, 1, lam_min * 1.05 / LYA_NM - 1, z_em - 0.05,
                                   logN_range=DLA_LOGN if seed == 0 else SUBDLA_LOGN)
    truth = MockTruth(
        components=sorted(comps, key=lambda c: c.z), lam_min_nm=lam_min,
        lam_max_nm=lam_max, dv_kms=dv, resol=resol, snr=snr, z_em=z_em,
        seed=int(seed), oversample=max(1, int(np.ceil(dv / MAX_SUBPIXEL_KMS))),
        emission_lines=QSO_EMISSION_LINES, cont_slope=-1.0, lyman_limit=True,
        noise_floor=NOISE_FLOOR)
    return truth, systems, strong


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
    from astrocook.recipes.absorbers import ABSORBERS_RECIPES_SCHEMAS
    from astrocook.recipes.continuum import CONTINUUM_RECIPES_SCHEMAS
    schema = {"continuum": CONTINUUM_RECIPES_SCHEMAS,
              "absorbers": ABSORBERS_RECIPES_SCHEMAS}[menu][recipe]
    return {p["name"]: str(p["default"]) for p in schema["params"]
            if not p.get("gui_hidden") and p["name"] != "debug_rating"}


def run_chain(spec):
    """Continuum, identification and doublet fit on a mock spectrum, as a
    user of the GUI would run them."""
    from astrocook.core.session import SessionV2
    session = SessionV2("mock", None, spec=spec)
    t0 = time.time()
    cont = session.continuum.estimate_auto(**schema_params("continuum", "estimate_auto"))
    t1 = time.time()
    ident = cont.absorbers.identify_lines(**schema_params("absorbers", "identify_lines"))
    t2 = time.time()
    fit = cont.absorbers.doublets_auto(**schema_params("absorbers", "doublets_auto"))
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
        for z in cands:
            beyond += z > truth.z_em
            hits = [i for i, s in enumerate(true)
                    if min(_dv(z, c.z) for c in s["components"]) <= MATCH_KMS]
            if not hits:
                false += 1
            elif hits[0] in found:
                duplicates += 1
            else:
                found.add(hits[0])
        out[species] = dict(
            n_true=len(true), n_detectable=int(sum(detectable)),
            n_candidates=len(cands), n_found=len(found),
            n_found_detectable=sum(1 for i in found if detectable[i]),
            n_duplicates=duplicates, n_false=false, n_beyond_z_em=beyond)
        # Significance of the weaker line and outcome, per system, to
        # check the detection limit
        out["systems"] += [dict(species=species, significance=d[1] / d[2],
                                found=i in found)
                           for i, d in enumerate(detection)]
    return out


def _flags(c):
    flags = []
    errors = [c.dz, c.dlogN, c.db]
    if any(e is None or not np.isfinite(e) for e in errors):
        flags.append("nan_errors")
    elif (c.dlogN > 1.0) or (c.db > c.b):
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
               true_logN=true_logN, fit_logN=None, dlogN=None)
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
    for dv, i, j in pairs:
        if dv > MATCH_KMS:
            break
        if i not in used_f and j not in used_t:
            used_f.add(i)
            used_t.add(j)
            matches.append((fitted[i], true[j], j))
    pulls = {"z": [], "logN": [], "b": []}
    diffs = {"z": [], "logN": [], "b": []}
    by_category = {k: {"z": [], "logN": [], "b": []} for k in PULL_CATEGORIES}
    for f, t, j in matches:
        category = pull_category(t, n_in_system[j], fitted)
        for par, err in (("z", f.dz), ("logN", f.dlogN), ("b", f.db)):
            d = getattr(f, par) - getattr(t, par)
            if par == "z":
                d = d / (1 + t.z) * C_KMS  # km/s
                err = err / (1 + t.z) * C_KMS if err is not None else None
            diffs[par].append(d)
            if err is not None and np.isfinite(err) and err > 0:
                pulls[par].append(d / err)
                by_category[category][par].append(d / err)
    flags = [_flags(c) for c in fitted]
    return dict(
        systems=[system_column(s, fitted) for s in systems],
        n_true=len(true), n_fitted=len(fitted), n_matched=len(matches),
        n_spurious=len(fitted) - len(matches), n_missed=len(true) - len(matches),
        n_flagged=sum(bool(f) for f in flags),
        flags={k: sum(k in f for f in flags)
               for k in ("nan_errors", "absurd_errors", "logN_at_bound", "b_at_bound")},
        diffs=diffs, pulls=pulls, pulls_by_category=by_category,
        chi2r=[float(c.chi2) for c in fitted if c.chi2 is not None])


def run_point(point):
    """Mock, chain and metrics of one grid point (a record)."""
    snr, resol, z_em, seed = point
    logging.disable(logging.WARNING)
    t0 = time.time()
    truth, systems, strong = make_truth(snr, resol, z_em, seed)
    # The true continuum is not in the spectrum, as for an observed one
    spec = generate(truth, truth_columns=False)
    t_mock = time.time() - t0
    cont, ident, fit, times = run_chain(spec)
    return dict(
        snr=snr, resol=resol, z_em=z_em, seed=seed, strong=strong,
        n_pixels=len(cont.spec.x), times=dict(mock=t_mock, **times),
        continuum=continuum_metrics(truth, cont.spec),
        identification=identification_metrics(truth, systems, ident.spec),
        fit=fit_metrics(truth, systems, fit.systs))


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
    out["seconds"] = {k: float(np.sum([r["times"][k] for r in records]))
                      for k in records[0]["times"]}
    return out


#: Bins of the significance of the weaker line (EW / sigma_EW)
SIGNIFICANCE_BINS = (0.0, 1.0, 2.0, 3.0, 5.0, 10.0, np.inf)


def completeness_by_significance(records):
    """Fraction of the true systems found, in bins of the significance of
    their weaker line: it should rise around DETECTION_SIGMA."""
    systems = [s for r in records for s in r["identification"]["systems"]]
    out = {}
    for lo, hi in zip(SIGNIFICANCE_BINS[:-1], SIGNIFICANCE_BINS[1:]):
        sel = [s["found"] for s in systems if lo <= s["significance"] < hi]
        out[f"{lo:g}-{hi:g}"] = dict(n=len(sel),
                                     found=float(np.mean(sel)) if sel else None)
    return out


def grid_points(grid):
    g = GRIDS[grid]
    return [(snr, resol, z_em, seed) for snr in g["snr"] for resol in g["resol"]
            for z_em in g["z_em"] for seed in g["seeds"]]


def run(grid="full", workers=1):
    points = grid_points(grid)
    if workers > 1:
        with Pool(workers) as pool:
            records = pool.map(run_point, points, chunksize=1)
    else:
        records = [run_point(p) for p in points]
    return dict(grid=grid, match_kms=MATCH_KMS, records=records,
                summary=dict(all=summarize(records),
                             by_snr=summarize(records, "snr"),
                             by_resol=summarize(records, "resol"),
                             by_z_em=summarize(records, "z_em")))


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
#: full grid. That of 2.1.0, before the fixes of phase 4, is kept in
#: validate_mocks_baseline_2.1.0.json
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


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--grid", choices=sorted(GRIDS), default="full")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--json", help="write the records and the summary here")
    parser.add_argument("--save-baseline", action="store_true",
                        help="write the summary of a full run as the baseline")
    parser.add_argument("--compare", action="store_true",
                        help="compare the summary of a full run with the baseline")
    args = parser.parse_args()
    t0 = time.time()
    result = run(args.grid, args.workers)
    result["seconds"] = time.time() - t0
    print_summary(result)
    print(f"{result['seconds']:.0f} s")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(result, f, indent=1)
    if args.compare and args.grid == "full":
        with open(BASELINE) as f:
            compare(result["summary"]["all"], json.load(f)["summary"]["all"])
    if args.save_baseline and args.grid == "full":
        from astrocook import __version__
        BASELINE.parent.mkdir(exist_ok=True)
        with open(BASELINE, "w") as f:
            json.dump(dict(astrocook_version=__version__, grid=args.grid,
                           match_kms=MATCH_KMS, summary=result["summary"]),
                      f, indent=1)


if __name__ == "__main__":
    main()
