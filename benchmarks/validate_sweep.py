"""The sweep of the anchored systems on mocks with rich systems (phase 7,
PR 7c.4d).

The metal systems of ``validate_mocks.py`` have one species each (and
their H I): on its grids every species that the sweep adds to a system
other than its H I is false, so they measure the false additions and the
H I of the systems, not whether the sweep finds the rarer lines. This
script makes the mocks of the grid ``full`` of ``validate_mocks.py`` (the
same forest, metals, DLA and sub-DLA), and gives half of the metal systems
(``RICH_FRACTION``, drawn per system from a generator of its own, so the
rest of the mock is unchanged) the lines of other ions at the redshift and
with the b of each of their components, at the column of the component
plus the offsets of ``RICH_IONS`` (indicative values of the low ions of
the absorbers of the forest and the red side, not a calibration): those of
the CIV and SiIV systems in the forest too (SiII 1190, 1193, 1260, 1304,
SiIII 1206, CII 1036 and 1334, OI 1039, 1302), on the red side SiII 1526
and AlII 1670; FeII and MgI with the MgII systems. It runs the continuum
and identify_lines by hypotheses (the species of the GUI, ``Ly_ab`` among
them, so that the forest has its assumed Ly-a), and measures with
``validate_mocks.sweep_metrics``: per species the tests of the sweep, those
accepted, right, with the species present, the purity of the rarer species
and of the H I, the H I of the metal systems, the default lines replaced,
the calibration of the weights; and the completeness of the ions of the
rich systems (a true ion of a system found by an accepted test).

Not a test: run it by hand. Usage::

    python benchmarks/validate_sweep.py [--workers N] [--json out.json]
        [--sweep HI|all|none]
"""
import argparse
import dataclasses
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

from astrocook.fitting.hypotheses_sweep import SWEEP_SPECIES  # noqa: E402
from astrocook.mocks import MockComponent, generate  # noqa: E402

GRID = dict(snr=(10.0, 30.0, 100.0), resol=(8900.0, 45000.0), z_em=(2.5, 3.5),
            seeds=(0, 1, 2))
#: The fraction of the metal systems with the other ions
RICH_FRACTION = 0.5
#: The other ions of a rich system of each species: (species of the sweep,
#: offset of logN from the component of the anchor's species), b the same
RICH_IONS = {
    "CIV": (("SiII", -0.6), ("CII", 0.0), ("SiIII", -0.4), ("OI", -0.8), ("AlII", -1.4)),
    "SiIV": (("SiII", -0.3), ("CII", 0.3), ("SiIII", 0.0), ("OI", -0.5), ("AlII", -1.1)),
    "MgII": (("FeII", -0.3), ("MgI", -1.2)),
}


def make_rich_truth(snr, resol, z_em, seed):
    """The mock of ``validate_mocks.make_truth`` with half of its metal
    systems rich. Returns the truth, the systems (with ``rich`` and the
    ions added) and the strong absorber."""
    truth, systems, strong = vm.make_truth(snr, resol, z_em, seed)
    rng = vm._rng(snr, resol, z_em, seed, 300)
    added = []
    for s in systems:
        s["rich"] = bool(rng.random() < RICH_FRACTION)
        s["ions"] = []
        if not s["rich"]:
            continue
        for species, offset in RICH_IONS[s["species"]]:
            s["ions"].append(species)
            added += [MockComponent(SWEEP_SPECIES[species], c.z, c.logN + offset, c.b)
                      for c in s["components"]]
    comps = sorted(list(truth.components) + added, key=lambda c: c.z)
    return dataclasses.replace(truth, components=comps), systems, strong


def ion_metrics(systems, sweep):
    """The ions of the rich systems found by the sweep: a test of the
    species at the system accepted and right."""
    out = {}
    tests = sweep["tests"] if sweep else []
    for s in systems:
        for species in s["ions"]:
            d = out.setdefault(species, dict(n=0, found=0))
            d["n"] += 1
            d["found"] += any(t["species"] == species and t["right"]
                              and vm._dv(t["z"], s["z"]) <= vm.SWEEP_PRESENT_KMS
                              for t in tests)
    return out


def run_point(point):
    snr, resol, z_em, seed = point
    logging.disable(logging.WARNING)
    truth, systems, strong = make_rich_truth(snr, resol, z_em, seed)
    spec = generate(truth, truth_columns=False)
    from astrocook.session import Session
    session = Session("mock", None, spec=spec)
    cont = session.continuum.estimate_auto(**vm.schema_params("continuum", "estimate_auto"))
    t0 = time.time()
    ident = cont.absorbers.identify_lines(**vm.chain_params("identify_lines"))
    t_ident = time.time() - t0
    sweep = vm.sweep_metrics(truth, systems, ident.spec)
    return dict(snr=snr, resol=resol, z_em=z_em, seed=seed, strong=strong,
                t_identify=t_ident, sweep=sweep, ions=ion_metrics(systems, sweep),
                identification=vm.identification_metrics(truth, systems, ident.spec),
                n_rich=sum(s["rich"] for s in systems), n_systems=len(systems))


def summarize(records):
    out = dict(n_mocks=len(records), t_identify=float(np.sum([r["t_identify"] for r in records])),
               t_identify_median=float(np.median([r["t_identify"] for r in records])),
               sweep=vm.sweep_summary([r["sweep"] for r in records]))
    ions = {}
    for r in records:
        for k, v in r["ions"].items():
            d = ions.setdefault(k, dict(n=0, found=0))
            d["n"] += v["n"]
            d["found"] += v["found"]
    for d in ions.values():
        d["completeness"] = d["found"] / d["n"] if d["n"] else None
    out["ions"] = ions
    metals = {}
    for sp in vm.METALS:
        n = {k: sum(r["identification"][sp][k] for r in records)
             for k in ("n_true", "n_found", "n_candidates", "n_false", "n_duplicates")}
        n["completeness"] = n["n_found"] / n["n_true"] if n["n_true"] else None
        n["purity"] = (1 - n["n_false"] / n["n_candidates"]) if n["n_candidates"] else None
        metals[sp] = n
    out["metals"] = metals
    return out


def run(workers=1):
    points = [(snr, resol, z_em, seed) for snr in GRID["snr"] for resol in GRID["resol"]
              for z_em in GRID["z_em"] for seed in GRID["seeds"]]
    if workers > 1:
        with Pool(workers, initializer=vm.set_chain, initargs=(dict(vm.CHAIN),)) as pool:
            records = pool.map(run_point, points, chunksize=1)
    else:
        records = [run_point(p) for p in points]
    summary = dict(all=summarize(records))
    for key in ("snr", "resol", "z_em"):
        summary[f"by_{key}"] = {str(v): summarize([r for r in records if r[key] == v])
                                for v in sorted({r[key] for r in records})}
    return dict(records=records, summary=summary)


def print_summary(s):
    fmt = (lambda v: "-" if v is None else f"{v:.3f}")
    print(f"{s['n_mocks']} mocks; identify_lines {s['t_identify']:.0f} s "
          f"(median {s['t_identify_median']:.1f} s)")
    print("  metals (completeness/purity): " + ", ".join(
        f"{k} {fmt(v['completeness'])}/{fmt(v['purity'])}" for k, v in s["metals"].items()))
    print("  ions of the rich systems found: " + ", ".join(
        f"{k} {v['found']}/{v['n']} ({fmt(v['completeness'])})"
        for k, v in sorted(s["ions"].items())))
    vm.print_sweep(s["sweep"])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--json", help="write the records and the summary here")
    parser.add_argument("--sweep", choices=["HI", "all", "none"], default="all",
                        help="the species of the sweep (default all: the rarer species "
                             "are what these mocks measure)")
    parser.add_argument("--multiplets", default="CIV,SiIV,MgII,Ly_ab",
                        help="the species of identify_lines")
    args = parser.parse_args()
    vm.set_chain(dict(method="hypotheses", multiplets=args.multiplets, sweep=args.sweep))
    t0 = time.time()
    result = run(args.workers)
    result["seconds"] = time.time() - t0
    print_summary(result["summary"]["all"])
    for key in ("snr", "resol", "z_em"):
        for value, s in result["summary"][f"by_{key}"].items():
            print(f"{key} = {value}:")
            print_summary(s)
    print(f"{result['seconds']:.0f} s")
    if args.json:
        Path(args.json).write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
