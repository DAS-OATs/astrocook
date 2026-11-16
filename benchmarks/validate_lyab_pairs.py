"""The Ly-a/Ly-b pairs of identify_lines on mocks (phase 7, PR 7c.1, F68).

The mocks of ``validate_mocks.py`` start at 104 nm in the rest frame of the
quasar, above the Ly-b emission: no Ly-b line of the forest falls in them,
so the grids ``full`` and ``robustness`` cannot measure the pairs. This
script makes the same mocks (forest with its Ly-b, metals, a DLA in seed 0
and a sub-DLA in seed 1) on a window that starts at ``LYAB_REST_NM[0]``,
inside the Ly-b forest, runs the continuum and identify_lines with the
defaults of the GUI (``Ly_ab`` among the multiplets) and measures:

- the pairs (the ``Ly_a`` candidates with a score: the Lyman-alpha by
  default of F51 has score 0) against the H I of the truth whose Ly-a and
  Ly-b both fall in the forest of the window: completeness by logN, over
  all such components and over those whose Ly-b line is detectable at
  ``DETECTION_SIGMA`` (as ``is_detectable`` of ``validate_mocks.py``);
  purity (a pair within ``MATCH_KMS`` of a true component with a
  detectable Ly-b), and the same match for the candidates moved by
  ``CHANCE_SHIFT_KMS``, the rate of a chance coincidence in the forest;
- the pixels where the Ly-b of the forest absorbs most (above the
  detectable depth of ``validate_mocks.py``): in no region, in a region
  labelled Ly-b by a pair (``region_identifications``), in a region with
  only the Lyman-alpha by default (wrong for these pixels), or other;
- the Ly-b lines with a Ly-a system on them (a ``Ly_a`` candidate within
  ``MATCH_KMS`` of z_b = (1 + z) lambda_b / lambda_a - 1 of a true
  component with a detectable Ly-b; also where a true Ly-a is there).

The candidates of a pair: identify_lines gives a ``Ly_a`` system, with the
score of the pair, at every peak of the region of the Ly-a line (the
de-blending), and none in the region of the Ly-b line, which is labelled
``Ly_b`` and no longer gets the Lyman-alpha by default.

Not a test: run it by hand. Usage::

    python benchmarks/validate_lyab_pairs.py [--workers N] [--json out.json]
        [--seeds 0 1 2]
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

from astrocook.mocks import generate  # noqa: E402
from astrocook.mocks.generate import (  # noqa: E402
    components_optical_depth, pixel_transmission, subpixel_grid,
)
from astrocook.mocks.physics import C_KMS, transitions  # noqa: E402

#: Rest-frame window of these mocks (nm): from the Ly-b forest (above the
#: Ly-g emission, 97.25 nm) to beyond the CIV emission, as in validate_mocks
LYAB_REST_NM = (97.5, vm.REST_WINDOW_NM[1])
# The window of the mocks of validate_mocks, changed here at import (also
# in the workers of a Pool, which import this module)
vm.REST_WINDOW_NM = LYAB_REST_NM

GRID = dict(snr=(10.0, 30.0, 100.0), resol=(8900.0, 45000.0), z_em=(2.5, 3.5, 4.5),
            seeds=(0, 1, 2))
#: Bins of logN of the completeness
LOGN_BINS = (12.0, 13.0, 13.5, 14.0, 15.0, 22.0)
#: Shift of the candidates for the chance coincidences (km/s, both signs)
CHANCE_SHIFT_KMS = 300.0
LAM_A, LAM_B = vm.LYA_NM, vm.LYB_NM


def _is_hi(comp):
    return all(t.startswith("Ly_") for t in transitions(comp.series))


def _pair_range(truth):
    """The H I components with both Ly-a and Ly-b in the forest of the
    window (between its start and the Ly-a emission)."""
    lam_min = truth.lam_min_nm
    return [c for c in truth.components
            if _is_hi(c) and "Ly_b" in transitions(c.series) and c.z < truth.z_em
            and LAM_B * (1 + c.z) > lam_min]


def lyb_detectable(comp, truth):
    """Whether the Ly-b line of one component is detected at
    DETECTION_SIGMA (the rule of ``vm.is_detectable`` for one component)."""
    ew = vm.equivalent_width_ang("Ly_b", comp.z, comp.logN, comp.b)
    fwhm = np.hypot(C_KMS / truth.resol, 1.665 * comp.b)
    lam = LAM_B * 10.0 * (1 + comp.z)
    sigma_ew = lam * truth.dv_kms / C_KMS * np.sqrt(max(1.0, 2 * fwhm / truth.dv_kms)) / truth.snr
    return ew >= vm.DETECTION_SIGMA * sigma_ew


def _without_lyb(comp):
    rest = [t for t in transitions(comp.series) if t != "Ly_b"]
    return dataclasses.replace(comp, series=",".join(rest)) if rest else None


def lyb_depths(truth):
    """Depth (1 - T) per pixel of the Ly-b lines of the forest alone and of
    everything else, as the truth map makes them."""
    x, sub = subpixel_grid(truth)
    lyb = [dataclasses.replace(c, series="Ly_b") for c in truth.components
           if _is_hi(c) and "Ly_b" in transitions(c.series)]
    rest = [c for c in (_without_lyb(c) if _is_hi(c) else c for c in truth.components)
            if c is not None]
    out = []
    for comps in (lyb, rest):
        tau = components_optical_depth(truth, sub, comps) if comps else np.zeros_like(sub)
        out.append(1.0 - pixel_transmission(truth, np.exp(-tau), len(x)))
    return x, out[0], out[1]


def _near(z, comps, kms=vm.MATCH_KMS):
    return [c for c in comps if vm._dv(z, c.z) <= kms]


def pair_metrics(truth, candidates):
    true_pairs = _pair_range(truth)
    detectable = [lyb_detectable(c, truth) for c in true_pairs]
    pairs = [c for c in candidates if c["series"] == "Ly_a" and c["score"] > 0]
    # Right: within MATCH_KMS of a true component with a detectable Ly-b
    # (any H I would do too often by chance in the dense forest: 3/4 of the
    # candidates moved by CHANCE_SHIFT_KMS on a mock at z_em 3.5)
    hi = [c for c, d in zip(true_pairs, detectable, strict=True) if d]
    found = [any(vm._dv(p["z"], c.z) <= vm.MATCH_KMS for p in pairs) for c in true_pairs]
    right = [bool(_near(p["z"], hi)) for p in pairs]
    chance = [bool(_near((1 + p["z"]) * (1 + s * CHANCE_SHIFT_KMS / C_KMS) - 1, hi))
              for p in pairs for s in (-1, 1)]
    # The Ly-b lines with a Ly-a system on them (the spurious systems of
    # F68, where nothing else is there)
    lya = [c for c in candidates if c["series"] == "Ly_a"]
    lyb_as_lya = 0
    for c in hi:
        z_b = (1 + c.z) * LAM_B / LAM_A - 1
        lyb_as_lya += any(vm._dv(cand["z"], z_b) <= vm.MATCH_KMS for cand in lya)
    by_logn = {}
    for lo, hi_ in zip(LOGN_BINS[:-1], LOGN_BINS[1:], strict=True):
        sel = [i for i, c in enumerate(true_pairs) if lo <= c.logN < hi_]
        by_logn[f"{lo:g}-{hi_:g}"] = dict(
            n=len(sel), found=sum(found[i] for i in sel),
            n_detectable=sum(detectable[i] for i in sel),
            found_detectable=sum(found[i] and detectable[i] for i in sel))
    return dict(n_true=len(true_pairs), n_detectable=int(sum(detectable)),
                n_found=int(sum(found)),
                n_found_detectable=int(sum(f and d for f, d in zip(found, detectable,
                                                                    strict=True))),
                n_pairs=len(pairs), n_right=int(sum(right)),
                chance_rate=float(np.mean(chance)) if chance else None,
                n_lyb_as_lya=lyb_as_lya, by_logn=by_logn,
                scores=[p["score"] for p in pairs])


PIXEL_CLASSES = ("undetected", "lyb_pair", "lya_default", "other")


def pixel_metrics(truth, spec, candidates):
    """The pixels where the Ly-b of the forest absorbs most, by outcome."""
    x, d_lyb, d_rest = lyb_depths(truth)
    sel = (d_lyb > vm.detectable_depth(truth.snr)) & (d_lyb > d_rest)
    ids = (np.asarray(spec.get_column("abs_ids").value).astype(int)
           if spec.has_aux_column("abs_ids") else np.zeros(len(x), dtype=int))
    labels = vm._decode_meta(spec.meta, "region_identifications")
    by_region = {}
    for c in candidates:
        by_region.setdefault(c["region"], []).append(c)
    out = dict.fromkeys(PIXEL_CLASSES, 0)
    for i in np.flatnonzero(sel):
        rid = int(ids[i])
        if rid == 0:
            out["undetected"] += 1
        elif any(line == "Ly_b" for line, _ in labels.get(str(rid), [])):
            out["lyb_pair"] += 1
        elif by_region.get(rid) and all(vm._is_default_lya(c) for c in by_region[rid]):
            out["lya_default"] += 1
        else:
            out["other"] += 1
    return out


def run_point(point):
    snr, resol, z_em, seed = point
    logging.disable(logging.WARNING)
    truth, _, strong = vm.make_truth(snr, resol, z_em, seed)
    spec = generate(truth, truth_columns=False)
    from astrocook.session import Session
    session = Session("mock", None, spec=spec)
    cont = session.continuum.estimate_auto(**vm.schema_params("continuum", "estimate_auto"))
    t0 = time.time()
    ident = cont.absorbers.identify_lines(**vm.schema_params("absorbers", "identify_lines"))
    t_ident = time.time() - t0
    candidates = vm.candidate_list(ident.spec)
    return dict(snr=snr, resol=resol, z_em=z_em, seed=seed, strong=strong,
                t_identify=t_ident, pairs=pair_metrics(truth, candidates),
                pixels=pixel_metrics(truth, ident.spec, candidates))


def summarize(records):
    keys = ("n_true", "n_detectable", "n_found", "n_found_detectable", "n_pairs",
            "n_right", "n_lyb_as_lya")
    tot = {k: sum(r["pairs"][k] for r in records) for k in keys}
    tot["completeness"] = tot["n_found"] / tot["n_true"] if tot["n_true"] else None
    tot["completeness_detectable"] = (tot["n_found_detectable"] / tot["n_detectable"]
                                      if tot["n_detectable"] else None)
    tot["purity"] = tot["n_right"] / tot["n_pairs"] if tot["n_pairs"] else None
    chance = [r["pairs"]["chance_rate"] for r in records
              if r["pairs"]["chance_rate"] is not None]
    weights = [r["pairs"]["n_pairs"] for r in records if r["pairs"]["chance_rate"] is not None]
    tot["chance_rate"] = (float(np.average(chance, weights=weights)) if chance else None)
    by_logn = {}
    for b in records[0]["pairs"]["by_logn"]:
        n = {k: sum(r["pairs"]["by_logn"][b][k] for r in records)
             for k in records[0]["pairs"]["by_logn"][b]}
        n["completeness"] = n["found"] / n["n"] if n["n"] else None
        n["completeness_detectable"] = (n["found_detectable"] / n["n_detectable"]
                                        if n["n_detectable"] else None)
        by_logn[b] = n
    tot["by_logn"] = by_logn
    scores = [s for r in records for s in r["pairs"]["scores"]]
    tot["score_median"] = float(np.median(scores)) if scores else None
    pixels = {k: sum(r["pixels"][k] for r in records) for k in PIXEL_CLASSES}
    n = sum(pixels.values())
    tot["pixels"] = dict(counts=pixels,
                         fractions={k: v / n if n else None for k, v in pixels.items()})
    tot["t_identify"] = float(np.sum([r["t_identify"] for r in records]))
    tot["n_mocks"] = len(records)
    return tot


def run(workers=1, seeds=GRID["seeds"]):
    points = [(snr, resol, z_em, seed) for snr in GRID["snr"] for resol in GRID["resol"]
              for z_em in GRID["z_em"] for seed in seeds]
    if workers > 1:
        with Pool(workers) as pool:
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
    print(f"{s['n_mocks']} mocks; pairs: {s['n_pairs']} candidates, purity "
          f"{fmt(s['purity'])} (chance {fmt(s['chance_rate'])}), completeness "
          f"{fmt(s['completeness'])} ({s['n_true']} true), Ly-b detectable "
          f"{fmt(s['completeness_detectable'])} ({s['n_detectable']}); Ly-b taken for "
          f"Ly-a: {s['n_lyb_as_lya']}; identify_lines {s['t_identify']:.0f} s")
    print("  by logN: " + ", ".join(
        f"{b}: {fmt(v['completeness'])} ({v['n']}), det. {fmt(v['completeness_detectable'])} "
        f"({v['n_detectable']})" for b, v in s["by_logn"].items()))
    print("  Ly-b pixels: " + ", ".join(f"{k} {fmt(v)}"
                                        for k, v in s["pixels"]["fractions"].items()))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--json", help="write the records and the summary here")
    parser.add_argument("--seeds", type=int, nargs="+", default=list(GRID["seeds"]))
    args = parser.parse_args()
    t0 = time.time()
    result = run(args.workers, tuple(args.seeds))
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
