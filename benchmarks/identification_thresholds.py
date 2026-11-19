"""The thresholds of the identification by hypotheses on mocks (phase 7, PR 7c.4c).

The two new heuristics of ``identify_lines(method="hypotheses")``
(``astrocook.fitting.hypotheses``) measured on a subset of the grid
``full`` of ``validate_mocks.py``, before Guido decides them (18 November
2026): the threshold of the pre-filter (the R^2 of a pair above which it is
an edge of the graph) and that of the screen (the difference of the
criterion, by the apparent optical depth, beyond which a hypothesis is
dropped before the short fit). For each mock the continuum is estimated
once (``estimate_auto`` with the GUI defaults), then ``identify_lines`` runs
with each pair of thresholds, and the current method for reference.

Per run: the completeness of the metal systems detectable at 3 sigma and
by the significance of their weaker line, the purity (CIV, SiIV, MgII, the
species of the mocks; the candidates of the other species, with
``--multiplets all``, are all false), the calibration of the weights (F164:
the expected calibration error and the Brier score of
``validate_mocks.hypothesis_calibration``), the hypotheses compared and the
time of the identification.

Not a test: run it by hand. Usage::

    python benchmarks/identification_thresholds.py [--workers 4]
        [--prefilters 0 0.2 0.3] [--screens 10 20 40 inf]
        [--multiplets CIV,SiIV,MgII,Ly_ab] [--snr ...] [--resol ...]
        [--z-em 3.5] [--seeds 0 2] [--json out.json]
"""
import argparse
import json
import logging
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
for path in (ROOT, ROOT / "benchmarks"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import validate_mocks as vm  # noqa: E402

from astrocook.mocks import generate  # noqa: E402

OPTIONS: dict = {}


def _init(options):
    OPTIONS.clear()
    OPTIONS.update(options)
    if options.get("max_added") is not None:
        from astrocook.fitting import hypotheses
        hypotheses.MAX_ADDED = int(options["max_added"])


def _runs():
    """The runs of each mock: the current method, then each pair of
    thresholds for each list of species."""
    out = [dict(method="doublets")]
    for mult in OPTIONS["multiplets"]:
        for p in OPTIONS["prefilters"]:
            for s in OPTIONS["screens"]:
                out.append(dict(method="hypotheses", multiplets=mult, prefilter=p,
                                screen=s if np.isfinite(s) else None))
    return out


def run_point(point):
    from astrocook.session import Session
    logging.disable(logging.WARNING)
    snr, resol, z_em, seed = point
    truth, systems, _ = vm.make_truth(snr, resol, z_em, seed)
    spec = generate(truth, truth_columns=False)
    cont = Session("mock", None, spec=spec).continuum.estimate_auto(
        **vm.schema_params("continuum", "estimate_auto"))
    base = vm.schema_params("absorbers", "identify_lines")
    records = []
    for run in _runs():
        params = dict(base, **{k: str(v) for k, v in run.items()})
        t0 = time.time()
        ident = cont.absorbers.identify_lines(**params)
        seconds = time.time() - t0
        cands = vm._candidates(ident.spec)
        structure = ident.spec.meta.get("identification")
        compared = [g for g in structure.groups if len(g.hypotheses) > 1]
        records.append(dict(
            run=run, snr=snr, resol=resol, z_em=z_em, seed=seed, seconds=seconds,
            identification=vm.identification_metrics(truth, systems, ident.spec),
            n_other=sum(1 for s, _ in cands if s not in vm.METALS and not s.startswith("Ly")),
            n_groups_compared=len(compared),
            n_hypotheses=sum(len(g.hypotheses) for g in compared),
            calibration=[dict(c, snr=snr, resol=resol)
                         for c in vm.hypothesis_calibration(truth, ident.spec)]))
    return records


def _key(run):
    if run["method"] == "doublets":
        return "doublets"
    screen = "none" if run["screen"] is None else f"{run['screen']:g}"
    return f"{run['multiplets']} | prefilter {run['prefilter']:g} | screen {screen}"


def summarize(records):
    out = {}
    for key in dict.fromkeys(_key(r["run"]) for r in records):
        rs = [r for r in records if _key(r["run"]) == key]
        species = {}
        for sp in vm.METALS:
            tot = {k: sum(r["identification"][sp][k] for r in rs)
                   for k in ("n_found", "n_false", "n_detectable", "n_found_detectable")}
            uniq = tot["n_found"] + tot["n_false"]
            species[sp] = dict(purity=tot["n_found"] / uniq if uniq else None,
                               completeness=(tot["n_found_detectable"] / tot["n_detectable"]
                                             if tot["n_detectable"] else None),
                               n_false=tot["n_false"])
        secs = np.array([r["seconds"] for r in rs])
        rows = [c for r in rs for c in r["calibration"]]
        out[key] = dict(
            n_mocks=len(rs), species=species,
            n_other=sum(r["n_other"] for r in rs),
            by_significance=vm.completeness_by_significance(rs),
            n_groups_compared=sum(r["n_groups_compared"] for r in rs),
            n_hypotheses=sum(r["n_hypotheses"] for r in rs),
            seconds=dict(median=float(np.median(secs)), p95=float(np.percentile(secs, 95)),
                         total=float(secs.sum())),
            calibration=vm.calibration_tables(rows))
    return out


def print_table(summary):
    sig_bins = ("2-3", "3-5", "5-10", "10-inf")
    print("| run | CIV c/p | SiIV c/p | MgII c/p | other | "
          + " | ".join(f"found {b}" for b in sig_bins)
          + " | groups | hyps | ECE | Brier | s median | s p95 |")
    print("|" + "---|" * (11 + len(sig_bins)))
    for key, s in summary.items():
        cp = []
        for sp in vm.METALS:
            c, p = s["species"][sp]["completeness"], s["species"][sp]["purity"]
            cp.append(f"{c:.2f}/{p:.2f}" if c is not None and p is not None else "-")
        found = []
        for b in sig_bins:
            v = s["by_significance"].get(b, {})
            found.append(f"{v['found']:.2f}" if v.get("found") is not None else "-")
        cal = s["calibration"]["all"] if s["calibration"] else None
        ece = f"{cal['ece']:.3f}" if cal else "-"
        brier = f"{cal['brier']:.3f}" if cal else "-"
        print(f"| {key} | " + " | ".join(cp) + f" | {s['n_other']} | " + " | ".join(found)
              + f" | {s['n_groups_compared']} | {s['n_hypotheses']} | {ece} | {brier} | "
              f"{s['seconds']['median']:.2f} | {s['seconds']['p95']:.2f} |")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--prefilters", type=float, nargs="*", default=[0.0, 0.2, 0.3])
    parser.add_argument("--screens", type=float, nargs="*", default=[10.0, 20.0, 40.0, np.inf])
    parser.add_argument("--multiplets", nargs="*", default=["CIV,SiIV,MgII,Ly_ab"])
    parser.add_argument("--snr", type=float, nargs="*", default=[10.0, 30.0, 100.0])
    parser.add_argument("--resol", type=float, nargs="*", default=[8900.0, 45000.0, 140000.0])
    parser.add_argument("--z-em", type=float, nargs="*", default=[3.5])
    parser.add_argument("--seeds", type=int, nargs="*", default=[0, 2])
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument("--max-added", type=int, default=None,
                        help="components added by the refinement (0: none)")
    args = parser.parse_args(argv)
    options = dict(prefilters=args.prefilters, screens=args.screens,
                   multiplets=args.multiplets, max_added=args.max_added)
    _init(options)
    points = [(s, r, z, seed) for s in args.snr for r in args.resol for z in args.z_em
              for seed in args.seeds]
    t0 = time.time()
    with Pool(args.workers, initializer=_init, initargs=(options,)) as pool:
        records = [r for rs in pool.imap_unordered(run_point, points) for r in rs]
    summary = summarize(records)
    print(f"{len(points)} mocks in {time.time() - t0:.0f} s")
    print_table(summary)
    if args.json:
        args.json.write_text(json.dumps(dict(options=dict(options, screens=[
            str(s) for s in args.screens]), summary=summary, records=records), indent=1))


if __name__ == "__main__":
    main()
