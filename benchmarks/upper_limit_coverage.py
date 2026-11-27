"""The coverage of the upper limits on the mocks (phase 7, PR 7d.7; F156).

On the mocks of ``validate_mocks.py`` (the grid ``full`` by default), after
the chain of the benchmark (continuum, identification, ``doublets_auto``),
for each true metal system that the chain did not find (no fitted component
of its species within ``MATCH_KMS`` of one of its components), the upper
limit on its species at its redshift (the column-weighted one), with the
methods ``profile`` and ``ew``, at 3 sigma, with b fixed at the default
(the instrumental b) and at the effective b of the system (its velocity
dispersion, sigma = b / sqrt 2 per component plus the spread of the
components, as one Gaussian).

The coverage is the fraction of these systems whose true total column (the
sum of its components) is below the limit; the systems are split by
``is_detectable`` of the benchmark (the weaker line of the doublet at 3
sigma in the noiseless mock): those below the threshold are the ones a
limit is for (the plan of 7d.7: "the fraction of the true systems below
the threshold whose log N is below the limit"), those above it are misses
of the chain, which the limit should see as a detection (significance
above 3).

Not a test: run it by hand. Usage::

    python benchmarks/upper_limit_coverage.py [--grid full|reduced|tiny]
        [--workers N] [--json out.json]
"""
import argparse
import json
import logging
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "benchmarks"))

import validate_mocks as vm  # noqa: E402

from astrocook.fitting.upper_limits import upper_limit  # noqa: E402
from astrocook.mocks import generate  # noqa: E402

SIGMA = 3.0


def _effective_b(system):
    comps = system["components"]
    w = np.array([10 ** c.logN for c in comps])
    w = w / w.sum()
    v = np.array([(c.z - system["z"]) / (1 + system["z"]) * vm.C_KMS for c in comps])
    var = float(np.sum(w * (np.array([c.b for c in comps]) ** 2 / 2 + v ** 2)))
    return float(np.clip(np.sqrt(2 * var), 1.0, 300.0))


def one_mock(point):
    snr, resol, z_em, seed = point
    logging.disable(logging.WARNING)
    warnings.simplefilter("ignore")
    t0 = time.time()
    truth, systems, _ = vm.make_truth(snr, resol, z_em, seed)
    spec = generate(truth, truth_columns=False)
    _, _, fit, _ = vm.run_chain(spec)
    fitted = [c for c in fit.systs.components] if fit.systs else []
    rows = []
    for system in systems:
        species = system["species"]
        found = any(c.series == species and any(vm._dv(c.z, t.z) <= vm.MATCH_KMS
                                                for t in system["components"])
                    for c in fitted)
        if found:
            continue
        detectable, _, _ = vm.is_detectable(system, truth)
        logN_true = float(np.log10(sum(10 ** c.logN for c in system["components"])))
        row = dict(snr=snr, resol=resol, z_em=z_em, seed=seed, species=species,
                   z=system["z"], logN_true=logN_true, detectable=bool(detectable),
                   b_eff=_effective_b(system), n_comp=len(system["components"]))
        for method in ("profile", "ew"):
            for b_name, b in (("binst", None), ("beff", row["b_eff"])):
                key = f"{method}_{b_name}"
                try:
                    lim = upper_limit(fit.spec, fit.systs, species, system["z"], b=b,
                                      sigma=SIGMA, method=method)
                    row[key] = dict(lim=lim.logN_lim, significance=lim.significance,
                                    n_refitted=len(lim.refitted))
                except (ValueError, RuntimeError) as e:
                    row[key] = dict(error=str(e))
        rows.append(row)
    return dict(point=point, rows=rows, seconds=time.time() - t0)


def summary(rows):
    out = {}
    for detectable in (False, True):
        sel = [r for r in rows if r["detectable"] == detectable]
        part = {"n": len(sel)}
        for key in ("profile_binst", "profile_beff", "ew_binst", "ew_beff"):
            ok = [r for r in sel if "lim" in r[key]]
            if not ok:
                continue
            covered = [r["logN_true"] < r[key]["lim"] for r in ok]
            margin = [r[key]["lim"] - r["logN_true"] for r in ok]
            part[key] = dict(
                n=len(ok), failed=len(sel) - len(ok),
                coverage=float(np.mean(covered)),
                margin_median=float(np.median(margin)),
                detected=float(np.mean([r[key]["significance"] > SIGMA for r in ok])))
        out["detectable" if detectable else "below_threshold"] = part
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--grid", default="full", choices=("full", "reduced", "tiny"))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--json", default=None)
    args = parser.parse_args()
    g = vm.GRIDS[args.grid]
    points = [(s, r, z, seed) for s in g["snr"] for r in g["resol"] for z in g["z_em"]
              for seed in g["seeds"]]
    t0 = time.time()
    with Pool(args.workers) as pool:
        results = pool.map(one_mock, points, chunksize=1)
    rows = [r for res in results for r in res["rows"]]
    report = {"grid": args.grid, "sigma": SIGMA, "mocks": len(points),
              "minutes": (time.time() - t0) / 60, "all": summary(rows),
              "by_resol": {str(r): summary([x for x in rows if x["resol"] == r])
                           for r in g["resol"]},
              "by_snr": {str(s): summary([x for x in rows if x["snr"] == s])
                         for s in g["snr"]},
              "by_species": {sp: summary([x for x in rows if x["species"] == sp])
                             for sp in vm.METALS}}
    print(json.dumps(report, indent=1))
    if args.json:
        with open(args.json, "w") as f:
            json.dump({**report, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
