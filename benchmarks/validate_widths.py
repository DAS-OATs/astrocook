"""The b of the default lines of the forest against the truth (phase 7, PR
7c.4d-quater; the measure of F174 direction 4).

The first pass of the identification by hypotheses fits every region of
the forest that no pair links with assumed Ly-a lines, b free. Direction 4
of F174 would keep the lines too narrow for an H I as candidates of
metals: does the b of those fits separate the true metal lines of the
forest from the H I, and at which S/N and R?

The mocks are those of ``validate_sweep.py`` (the grid ``full`` of
``validate_mocks.py`` with half of the metal systems rich in low ions, whose
lines fall in the forest), at S/N 10, 30, 100, R 8900, 45000, 140000, z_em
3.5, seeds 0-2 (27 mocks); the chain of ``validate_mocks.py`` with the
method ``hypotheses`` and no sweep. Each assumed Ly-a of an accepted
hypothesis is matched to the true line nearest to its Ly-a within the
larger of half the instrumental FWHM and ``MATCH_KMS`` (the strongest, by N
f lambda, when several): ``HI`` (a Lyman line), ``metal`` or ``none``.
Reported per S/N and R: the medians of b and of its error, the fraction
below 12 km/s, the fraction fitted (with an error: in a group compared; the
others keep the AOD parameters of a region that no pair links), and the classification by the ratio of the priors of b of
:mod:`astrocook.fitting.evidence` (``B_METAL`` against ``B_FOREST``): the
fraction of the metals and of the H I called metal. The b of the mocks
(forest lognormal of median 25 km/s with a floor at 10, metals 4-15 km/s,
``astrocook.mocks.populations``) separate more than real spectra do: an
optimistic test.

Not a test: run it by hand. Usage::

    python benchmarks/validate_widths.py [--workers N] [--json out.json]
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
import validate_sweep as vs  # noqa: E402

from astrocook.core.atomic_data import ATOM_DATA, series_lines  # noqa: E402
from astrocook.fitting import evidence as E  # noqa: E402
from astrocook.mocks import generate  # noqa: E402
from astrocook.mocks.physics import C_KMS  # noqa: E402

GRID = dict(snr=(10.0, 30.0, 100.0), resol=(8900.0, 45000.0, 140000.0), z_em=(3.5,),
            seeds=(0, 1, 2))
#: The least distance of a match (km/s)
MATCH_KMS = 10.0
#: The threshold of "narrow" (km/s)
NARROW_KMS = 12.0


def _true_lines(truth):
    """(observed wavelength in A, N f lambda, class) of every line of the
    truth."""
    out = []
    for c in truth.components:
        cls = "HI" if c.series.startswith("Ly") else "metal"
        for s in c.series.split(","):
            for t in series_lines(s):
                if t in ATOM_DATA and ATOM_DATA[t]["f"] > 0:
                    a = ATOM_DATA[t]
                    out.append((a["wave"] * (1 + c.z), 10 ** c.logN * a["f"] * a["wave"], cls))
    return np.array([x[0] for x in out]), np.array([x[1] for x in out]), [x[2] for x in out]


def run_point(point):
    snr, resol, z_em, seed = point
    logging.disable(logging.WARNING)
    truth, _, _ = vs.make_rich_truth(snr, resol, z_em, seed)
    spec = generate(truth, truth_columns=False)
    from astrocook.core.identification import META_KEY
    from astrocook.session import Session
    session = Session("mock", None, spec=spec)
    cont = session.continuum.estimate_auto(**vm.schema_params("continuum", "estimate_auto"))
    ident = cont.absorbers.identify_lines(**vm.chain_params("identify_lines"))
    structure = ident.spec.meta[META_KEY]
    lam, strength, cls = _true_lines(truth)
    fwhm = C_KMS / resol
    tol = max(0.5 * fwhm, MATCH_KMS)
    lines = []
    for g in structure.groups:
        for h in g.hypotheses:
            if str(h.status) != "accepted":
                continue
            for c in h.components:
                if c.series != "Ly_a":
                    continue
                lo = ATOM_DATA["Ly_a"]["wave"] * (1 + c.z)
                dv = np.abs(np.log(lam / lo)) * C_KMS
                near = np.flatnonzero(dv <= tol)
                kind = cls[near[np.argmax(strength[near])]] if near.size else "none"
                lines.append(dict(b=float(c.b), db=None if c.db is None else float(c.db),
                                  kind=kind))
    return dict(snr=snr, resol=resol, z_em=z_em, seed=seed, lines=lines)


def _prior_says_metal(b):
    return E.B_METAL.terms(b)[0] > E.B_FOREST.terms(b)[0]


def summarize(lines):
    out = {}
    for kind in ("HI", "metal", "none"):
        sel = [x for x in lines if x["kind"] == kind]
        b = np.array([x["b"] for x in sel])
        db = np.array([x["db"] for x in sel if x["db"] is not None and np.isfinite(x["db"])])
        out[kind] = dict(
            n=len(sel), b_median=float(np.median(b)) if b.size else None,
            db_median=float(np.median(db)) if db.size else None,
            # With an error: fitted (a group compared); without: the AOD
            # parameters of a region that no pair links
            fitted=float(db.size / len(sel)) if sel else None,
            narrow=float(np.mean(b < NARROW_KMS)) if b.size else None,
            called_metal=float(np.mean([_prior_says_metal(x) for x in b])) if b.size else None)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--json")
    args = parser.parse_args()
    vm.set_chain(dict(method="hypotheses", multiplets="CIV,SiIV,MgII,Ly_ab", sweep="none"))
    points = [(s, r, z, k) for s in GRID["snr"] for r in GRID["resol"] for z in GRID["z_em"]
              for k in GRID["seeds"]]
    t0 = time.time()
    with Pool(args.workers, initializer=vm.set_chain, initargs=(dict(vm.CHAIN),)) as pool:
        records = pool.map(run_point, points, chunksize=1)
    print(f"{len(records)} mocks in {time.time() - t0:.0f} s")
    print("| S/N | R | H I n | b | db | <12 | called metal | fitted | metal n | b | db | <12 "
          "| called metal | fitted | none n |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    fmt = (lambda v: "-" if v is None else f"{v:.2f}")
    rows = {}
    for s in GRID["snr"]:
        for r in GRID["resol"]:
            lines = [x for rec in records if rec["snr"] == s and rec["resol"] == r
                     for x in rec["lines"]]
            d = summarize(lines)
            rows[f"{s:g}_{r:g}"] = d
            h, m = d["HI"], d["metal"]
            print(f"| {s:g} | {r:g} | {h['n']} | {fmt(h['b_median'])} | {fmt(h['db_median'])} | "
                  f"{fmt(h['narrow'])} | {fmt(h['called_metal'])} | {fmt(h['fitted'])} | {m['n']} | "
                  f"{fmt(m['b_median'])} | {fmt(m['db_median'])} | {fmt(m['narrow'])} | "
                  f"{fmt(m['called_metal'])} | {fmt(m['fitted'])} | {d['none']['n']} |")
    if args.json:
        Path(args.json).write_text(json.dumps(dict(records=records, summary=rows), indent=1))


if __name__ == "__main__":
    main()
