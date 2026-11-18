"""The detection of the lines by matched filters on the mocks (phase 7,
PR 7c.4b, F57).

On the mocks of ``validate_mocks.py`` (grid ``full`` or ``robustness``),
after the continuum of the GUI (``estimate_auto`` with its defaults, the
common decision of 7b and 7c), three detections of the absorbed regions
are measured against the truth map of the mock
(:func:`astrocook.mocks.truth_map`):

- ``abs_mask``: the regions of ``identify_lines`` today (the mask of the
  continuum, ``detect_regions`` with 3 pixels at least, merged within
  10 km/s);
- ``mf_single``: the matched filters of
  :func:`astrocook.core.line_detection.detect_lines`, the lines alone;
- ``mf``: the same with the joint significance of the doublets CIV, SiIV,
  MgII (no doublet beyond z_em by more than 3000 km/s, the margin of F55);
- ``mf_j4``: the same with the doublets at 4 sigma instead of 3.

The measures:

- the absorbed pixels of the truth (total depth above ``ABSORBED_DEPTH``)
  in a region, by the class of their true significance (the matched
  filters applied to the noiseless depth of the truth with the error of
  the mock, the best width; ``SIGMA_BINS``), by zone (forest, red side)
  and by the category that absorbs them most (the rows of 7a.2);
- the true metal systems (``systems`` of the mocks) by the significance of
  their weaker line (``is_detectable``, the bins of F57): the weaker line
  in a region (a pixel within ``MATCH_KMS`` of it at the redshift of a
  component), and both lines;
- the false regions: those where the truth is below 1 sigma for every
  width (no pixel with true significance above 1: noise, or a continuum
  that is off), and those with no absorbed pixel at all; by zone;
- the pixels in the regions, by zone, and the time of the detection.

Not a test: run it by hand. Usage::

    python benchmarks/validate_detection.py [--grid full|robustness|reduced|tiny]
        [--workers N] [--json out.json]
"""
import argparse
import json
import logging
import sys
import time
from itertools import pairwise
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
for p in (ROOT, ROOT / "benchmarks"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS  # noqa: E402
from astrocook.core.line_detection import (  # noqa: E402
    DEFAULT_MULTIPLETS, detect_lines, matched_filter_detection,
)
from astrocook.core.spectrum_operations import (  # noqa: E402
    detect_regions, merge_regions_by_velocity,
)
from astrocook.mocks import continuum, generate  # noqa: E402
from astrocook.mocks.truthmap import CATEGORIES, truth_map  # noqa: E402
from validate_mocks import (  # noqa: E402
    ABSORBED_DEPTH, C_KMS, GRIDS, LYA_NM, LYB_NM, MATCH_KMS, SIGNIFICANCE_BINS,
    dominant_row, grid_points, is_detectable, make_robust_truth, make_truth,
    schema_params, total_depth, weak_line, z_em_perturbed,
)

#: The methods compared
METHODS = ("abs_mask", "mf_single", "mf", "mf_j4")
#: The threshold of the joint significance of ``mf_j4``, an alternative
#: measured next to the 3 sigma of the design
J4_SIGMA = 4.0
#: Classes of the true significance of a pixel
SIGMA_BINS = (0.0, 1.0, 3.0, 5.0, 10.0, np.inf)
#: Zones of the spectrum (rest frame of the true z_em)
ZONES = ("forest", "red")
#: Margin of the doublets beyond z_em (km/s), that of identify_lines (F55)
Z_EM_MARGIN_KMS = 3000.0
#: Rows of the pixels: the category that absorbs them most
ROWS = CATEGORIES


def _bin_labels(bins):
    return [f"{lo:g}-{hi:g}" for lo, hi in pairwise(bins)]


def _zones(truth, x):
    rest = x / (1 + truth.z_em)
    return {"forest": (rest > LYB_NM) & (rest < LYA_NM), "red": rest >= LYA_NM}


def abs_mask_regions(spec, z_em):
    """The regions of identify_lines from the mask of the continuum."""
    mask = np.asarray(spec.get_column("abs_mask").value).astype(bool)
    raw, _ = detect_regions(mask, 3)
    merged, _ = merge_regions_by_velocity(raw, spec.x, z_em, 10.0)
    return merged


def true_significance(truth, spec, depth):
    """The best S of the filters on the noiseless depth of the truth, with
    the error of the mock over the true continuum."""
    x = spec.x.to_value("nm")
    err = spec.dy.value / continuum(truth, x)
    det = matched_filter_detection(x, 1.0 - depth, err, spec.resolution(),
                                   multiplets=())
    return np.nan_to_num(det.max_significance, nan=0.0)


def pixel_table(regions, s_true, depth, row, zones):
    """Absorbed pixels and those in a region, by zone, row and class of
    the true significance."""
    absorbed = depth > ABSORBED_DEPTH
    cls = np.digitize(s_true, SIGMA_BINS[1:-1])
    inside = regions > 0
    out = {}
    for zone, zmask in zones.items():
        out[zone] = {}
        for r, name in enumerate(ROWS):
            sel = absorbed & zmask & (row == r)
            out[zone][name] = [[int(np.sum(sel & (cls == k))),
                                int(np.sum(sel & (cls == k) & inside))]
                               for k in range(len(SIGMA_BINS) - 1)]
    return out


def _line_pixels(x, transition, zs):
    lam0 = ATOM_DATA[transition]["wave"] / 10.0
    out = []
    for z in zs:
        lam = lam0 * (1 + z)
        dv = np.abs(np.log(x / lam)) * C_KMS
        out.extend(np.flatnonzero(dv <= MATCH_KMS))
    return np.array(out, dtype=int)


def system_table(truth, systems, x, regions_by_method):
    """Per true metal system: the significance of its weaker line and, per
    method, whether its weaker line and both lines are in a region."""
    out = []
    for s in systems:
        zs = [c.z for c in s["components"]]
        _, ew, sigma_ew = is_detectable(s, truth)
        weak = weak_line(s["species"])
        pix = {t: _line_pixels(x, t, zs) for t in STANDARD_MULTIPLETS[s["species"]]}
        rec = dict(species=s["species"], significance=float(ew / sigma_ew))
        for m, reg in regions_by_method.items():
            found = {t: bool(len(p) and (reg[p] > 0).any()) for t, p in pix.items()}
            rec[m] = [found[weak], all(found.values())]
        out.append(rec)
    return out


def region_table(regions, s_true, depth, zones):
    """Regions by zone (of their first pixel): all, without a pixel of true
    significance above 1, without an absorbed pixel."""
    out = {z: dict(n=0, noise=0, unabsorbed=0) for z in zones}
    n = int(regions.max(initial=0))
    if n == 0:
        return out
    idx = np.arange(1, n + 1)
    from scipy.ndimage import maximum, minimum_position
    s_max = np.asarray(maximum(s_true, regions, idx))
    d_max = np.asarray(maximum(depth, regions, idx))
    first = np.array([p[0] for p in minimum_position(np.arange(len(regions)),
                                                       regions, idx)])
    for z, zmask in zones.items():
        sel = zmask[first]
        out[z] = dict(n=int(sel.sum()), noise=int(np.sum(sel & (s_max <= 1.0))),
                      unabsorbed=int(np.sum(sel & (d_max <= ABSORBED_DEPTH))))
    return out


def run_point(point):
    from astrocook.session import Session
    logging.disable(logging.WARNING)
    snr, resol, z_em, seed = point[:4]
    extra = {}
    if len(point) > 4:
        truth, systems, _ = make_robust_truth(*point)
        z_used, _ = z_em_perturbed(snr, resol, z_em, seed)
        extra = dict(case=point[4])
    else:
        truth, systems, _ = make_truth(snr, resol, z_em, seed)
        z_used = z_em
    spec = generate(truth, truth_columns=False).with_properties(z_em=z_used)
    cont = Session("mock", None, spec=spec).continuum.estimate_auto(
        **schema_params("continuum", "estimate_auto")).spec
    z_max = z_used + Z_EM_MARGIN_KMS / C_KMS * (1 + z_used)
    t0 = time.time()
    single = detect_lines(cont, multiplets=())
    t1 = time.time()
    full = detect_lines(cont, multiplets=DEFAULT_MULTIPLETS, z_max=z_max)
    t2 = time.time()
    j4 = detect_lines(cont, multiplets=DEFAULT_MULTIPLETS, z_max=z_max,
                      multiplet_sigma=J4_SIGMA)
    regions = dict(abs_mask=abs_mask_regions(cont, z_used),
                   mf_single=single.regions, mf=full.regions, mf_j4=j4.regions)
    tmap = truth_map(truth)
    depth = total_depth(tmap)
    row = dominant_row(tmap)
    x = cont.x.to_value("nm")
    s_true = true_significance(truth, cont, depth)
    zones = _zones(truth, x)
    return dict(
        snr=snr, resol=resol, z_em=z_em, seed=seed, **extra,
        pixels={m: pixel_table(r, s_true, depth, row, zones) for m, r in regions.items()},
        systems=system_table(truth, systems, x, regions),
        regions={m: region_table(r, s_true, depth, zones) for m, r in regions.items()},
        masked={m: {z: int(np.sum((r > 0) & zm)) for z, zm in zones.items()}
                for m, r in regions.items()},
        n_zone={z: int(zm.sum()) for z, zm in zones.items()},
        km_s_zone={z: float(np.sum(zm) * truth.dv_kms) for z, zm in zones.items()},
        seconds=dict(single=t1 - t0, joint=t2 - t1))


# --- Summary ------------------------------------------------------------

def _frac(a, b):
    return float(a / b) if b else None


def summarize(records):
    out = {}
    labels = _bin_labels(SIGMA_BINS)
    for m in METHODS:
        pix = {}
        for zone in ZONES:
            tot = np.zeros((len(labels), 2), dtype=int)
            for r in records:
                for name in ROWS:
                    tot += np.array(r["pixels"][m][zone][name])
            pix[zone] = {lab: dict(n=int(t[0]), found=_frac(t[1], t[0]))
                         for lab, t in zip(labels, tot, strict=True)}
            pix[f"{zone}_metal"] = {
                lab: dict(n=int(t[0]), found=_frac(t[1], t[0]))
                for lab, t in zip(labels, np.sum(
                    [r["pixels"][m][zone]["metal"] for r in records], axis=0),
                    strict=True)}
        systems = [s for r in records for s in r["systems"]]
        sys_tab = {}
        for lo, hi in pairwise(SIGNIFICANCE_BINS):
            sel = [s[m] for s in systems if lo <= s["significance"] < hi]
            sys_tab[f"{lo:g}-{hi:g}"] = dict(
                n=len(sel), weak=_frac(sum(w for w, _ in sel), len(sel)),
                both=_frac(sum(b for _, b in sel), len(sel)))
        regs = {}
        for zone in ZONES:
            n = sum(r["regions"][m][zone]["n"] for r in records)
            noise = sum(r["regions"][m][zone]["noise"] for r in records)
            unabs = sum(r["regions"][m][zone]["unabsorbed"] for r in records)
            kms = sum(r["km_s_zone"][zone] for r in records)
            regs[zone] = dict(n=n, noise=noise, unabsorbed=unabs,
                              noise_per_1e4_kms=_frac(1e4 * noise, kms),
                              masked=_frac(sum(r["masked"][m][zone] for r in records),
                                           sum(r["n_zone"][zone] for r in records)))
        out[m] = dict(pixels=pix, systems=sys_tab, regions=regs)
    out["seconds"] = {k: float(np.median([r["seconds"][k] for r in records]))
                      for k in ("single", "joint")}
    return out


def _fmt(v):
    return "  -  " if v is None else f"{v:5.2f}"


def print_summary(summary, title):
    print(f"\n== {title}")
    labels = _bin_labels(SIGMA_BINS)
    for key in ("forest", "red", "red_metal"):
        print(f"absorbed pixels in a region, {key}, by true S: " + " ".join(f"{lab:>8}" for lab in labels))
        for m in METHODS:
            row = summary[m]["pixels"][key]
            print(f"  {m:10} " + " ".join(f"{_fmt(row[lab]['found']):>8}" for lab in labels))
        print("  n          " + " ".join(f"{summary[METHODS[0]]['pixels'][key][lab]['n']:>8}" for lab in labels))
    sys_labels = list(summary[METHODS[0]]["systems"])
    print("metal systems, weaker line / both lines in a region, by S of the weaker line:")
    print("             " + " ".join(f"{lab:>11}" for lab in sys_labels))
    for m in METHODS:
        tab = summary[m]["systems"]
        print(f"  {m:10} " + " ".join(f"{_fmt(tab[lab]['weak'])}/{_fmt(tab[lab]['both'])}" for lab in sys_labels))
    print("  n          " + " ".join(f"{summary[METHODS[0]]['systems'][lab]['n']:>11}" for lab in sys_labels))
    print("regions: n, noise (true S <= 1), unabsorbed, noise per 1e4 km/s, fraction of pixels in regions")
    for zone in ZONES:
        for m in METHODS:
            r = summary[m]["regions"][zone]
            print(f"  {zone:6} {m:10} {r['n']:7} {r['noise']:6} {r['unabsorbed']:6} "
                  f"{_fmt(r['noise_per_1e4_kms'])} {_fmt(r['masked'])}")
    print(f"median time of the detection: lines {summary['seconds']['single']:.2f} s, "
          f"with the doublets {summary['seconds']['joint']:.2f} s")


def run(grid, workers):
    points = grid_points(grid)
    if workers > 1:
        with Pool(workers) as pool:
            records = pool.map(run_point, points, chunksize=1)
    else:
        records = [run_point(p) for p in points]
    summary = dict(all=summarize(records))
    for key in ("snr", "resol"):
        for val in GRIDS[grid][key]:
            summary[f"{key}={val:g}"] = summarize([r for r in records if r[key] == val])
    for case in GRIDS[grid].get("cases", ()):
        summary[f"case={case}"] = summarize([r for r in records if r["case"] == case])
    return dict(grid=grid, records=records, summary=summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--grid", choices=sorted(GRIDS), default="full")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--json", help="write the records and the summary here")
    args = parser.parse_args()
    t0 = time.time()
    result = run(args.grid, args.workers)
    result["seconds"] = time.time() - t0
    for key, summary in result["summary"].items():
        print_summary(summary, f"{args.grid}, {key}")
    print(f"{result['seconds']:.0f} s")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(result, f, indent=1)


if __name__ == "__main__":
    main()
