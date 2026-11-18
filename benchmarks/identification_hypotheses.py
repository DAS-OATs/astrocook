"""Controlled pairs of alternative identifications on mocks (phase 7, PR 7c.2).

The experiment of the design document ``design/identification.md``
(section 4): which criterion of ``evaluate_model`` (chi2, AIC, AICc, BIC,
SpIC H, a, b) picks the true one of two alternative identifications of the
same absorption, by S/N, resolving power and significance of the lines.

Each trial is a small mock (a flat continuum, a window of a few nm around
the candidate, no forest) with a known truth, and two hypotheses for it, one
true and one false, both measured on the same pixels (the union of the
windows around the lines of both). The scenarios:

- ``doublet_vs_lya``: a CIV doublet; the alternatives a CIV system (true)
  and two independent Ly-a lines on the two lines (the interloper);
- ``doublet_vs_lyab``: as ``doublet_vs_lya``, the interlopers modelled
  with their Ly-b (series ``Ly_ab``), whose pixels are added (in these
  mocks the Ly-b falls on a clean continuum: in a spectrum it falls in the
  Ly-b forest, so this is the best case of the forward model of a series);
- ``coincidence``: two unrelated Ly-a lines at the separation of the CIV
  doublet (within ``COINCIDENCE_KMS``), with independent logN and b; the
  alternatives two Ly-a (true) and a CIV system (the chance coincidence);
- ``alias``: a CIV doublet; the alternatives a CIV system (true) and the
  CIV at the redshift that puts its 1548 on the true 1550, with a Ly-a on
  the true 1548 (the wrong pairing of a kinematic search);
- ``detection``: a CIV doublet, from below the detection limit; the
  alternatives a CIV system (true) and no absorber (unidentified noise).

Each hypothesis is evaluated twice: without a fit, with N and b from the
apparent optical depth (AOD) of its lines (``aod``), and after a short fit
started there (``fit``). Every component of a hypothesis is counted
(``min_influence=0``): a hypothesis is a model, and its parameters are its
cost whatever their influence.

The prior on the incidence of the ions: for the scenarios with a Ly-a
alternative, the decision of each criterion is also given with the prior
odds of the Ly-a hypothesis against the CIV one (``PRIOR_ODDS``), added as
2 ln(odds) to the difference of the criterion.

Not a test: run it by hand (about 2 minutes on 4 cores with the defaults,
40 trials per scenario and cell; 4.5 minutes with 100, the run of the
document).
Usage::

    python benchmarks/identification_hypotheses.py [--trials 40]
        [--workers 4] [--json out.json]
"""
import argparse
import dataclasses
import json
import logging
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
for path in (ROOT, ROOT / "benchmarks"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import validate_mocks as vm  # noqa: E402

from astrocook.core.atomic_data import ATOM_DATA  # noqa: E402
from astrocook.core.structures import SystemListData  # noqa: E402
from astrocook.core.system_list import SystemList  # noqa: E402
from astrocook.fitting.model_evaluation import evaluate_model  # noqa: E402
from astrocook.fitting.voigt_fitter import VoigtFitter  # noqa: E402
from astrocook.mocks import MockComponent, MockTruth, generate  # noqa: E402
from astrocook.mocks.physics import C_KMS  # noqa: E402

SCENARIOS = ("doublet_vs_lya", "doublet_vs_lyab", "coincidence", "alias", "detection")
CRITERIA = ("chi2", "aic", "aicc", "bic", "spic", "spic_a", "spic_b")
MODES = ("aod", "fit")
GRID = dict(snr=(10.0, 30.0, 100.0), resol=(8900.0, 45000.0))
#: Bins of the significance of the weaker line (sigma, as ``vm.is_detectable``)
SIG_BINS = (0.0, 3.0, 5.0, 10.0, np.inf)
SIG_LABELS = ("<3", "3-5", "5-10", ">10")
#: Prior odds Ly-a : CIV for the scenarios with a Ly-a alternative
PRIOR_ODDS = (1.0, 10.0, 100.0)

#: The candidate: CIV_1548 at about 500 nm, in the forest of a quasar at
#: z_em 3.5 (the forest itself is not simulated: no blends)
Z_CAND = 500.0 / (ATOM_DATA["CIV_1548"]["wave"] / 10.0) - 1.0
Z_EM = 3.5
LAM_1548 = ATOM_DATA["CIV_1548"]["wave"]
LAM_1550 = ATOM_DATA["CIV_1550"]["wave"]
LAM_LYA = ATOM_DATA["Ly_a"]["wave"]
#: Half width of the window of the pixels around every line of either
#: hypothesis (km/s)
PIXEL_HALF_KMS = 200.0
#: Half width of the AOD integration around a line, beyond the instrumental
#: FWHM (km/s)
AOD_HALF_KMS = 60.0
#: The largest offset of the second line of a chance coincidence from the
#: position of the CIV doublet (km/s): within the tolerance of a search
COINCIDENCE_KMS = 5.0
#: Ranges of the draws
LOGN_CIV = (11.8, 14.8)
B_CIV = (4.0, 20.0)
LOGN_LYA = (12.3, 14.5)
B_LYA = (15.0, 40.0)
#: Bounds of the short fit (those of the automatic pipelines, F62)
FIT_NFEV = 300
FIT_B_MAX = 100.0
FIT_LOGN_MAX = 19.0
#: AOD: N = AOD_CONST / (f lambda[A]) * int tau dv[km/s]
AOD_CONST = 3.768e14


@dataclasses.dataclass(frozen=True)
class Line:
    """One component of a hypothesis: its series and its redshift (the
    candidate position); logN and b come from the AOD or the fit."""
    series: str
    z: float


@dataclasses.dataclass(frozen=True)
class Trial:
    scenario: str
    snr: float
    resol: float
    seed: int


# --- The truth ------------------------------------------------------------

def z_of_lya_at(lam_ang: float) -> float:
    return lam_ang / LAM_LYA - 1.0


def draw_truth(trial: Trial) -> Tuple[List[MockComponent], Dict[str, List[Line]], str]:
    """The true components, the two hypotheses (true first) and the
    transition whose significance bins the trial."""
    rng = np.random.default_rng(1000 * SCENARIOS.index(trial.scenario) + trial.seed)
    z = Z_CAND + rng.uniform(-2e-4, 2e-4)
    lam1, lam2 = LAM_1548 * (1 + z), LAM_1550 * (1 + z)
    civ = Line("CIV", z)
    if trial.scenario == "coincidence":
        lam2 *= 1 + rng.uniform(-COINCIDENCE_KMS, COINCIDENCE_KMS) / C_KMS
        lya = [Line("Ly_a", z_of_lya_at(lam1)), Line("Ly_a", z_of_lya_at(lam2))]
        comps = [MockComponent("Ly_a", ln.z, rng.uniform(*LOGN_LYA), rng.uniform(*B_LYA))
                 for ln in lya]
        return comps, {"true": lya, "false": [civ]}, "Ly_a"
    comps = [MockComponent("CIV", z, rng.uniform(*LOGN_CIV), rng.uniform(*B_CIV))]
    if trial.scenario in ("doublet_vs_lya", "doublet_vs_lyab"):
        series = "Ly_a" if trial.scenario == "doublet_vs_lya" else "Ly_ab"
        false = [Line(series, z_of_lya_at(lam1)), Line(series, z_of_lya_at(lam2))]
    elif trial.scenario == "alias":
        false = [Line("CIV", lam2 / LAM_1548 - 1), Line("Ly_a", z_of_lya_at(lam1))]
    else:
        false = []
    return comps, {"true": [civ], "false": false}, "CIV_1550"


def make_truth(trial: Trial, comps: List[MockComponent], hyps) -> MockTruth:
    lams = [lam / 10.0 for ln in hyps["true"] + hyps["false"]
            for lam in line_waves(ln)]
    pad = 2 * PIXEL_HALF_KMS / C_KMS
    dv = C_KMS / trial.resol / vm.PIXELS_PER_FWHM
    return MockTruth(
        components=comps, lam_min_nm=min(lams) * (1 - pad), lam_max_nm=max(lams) * (1 + pad),
        dv_kms=dv, resol=trial.resol, snr=trial.snr, z_em=Z_EM, seed=trial.seed,
        oversample=max(1, int(np.ceil(dv / vm.MAX_SUBPIXEL_KMS))))


SERIES_LINES = {"CIV": ["CIV_1548", "CIV_1550"], "Ly_a": ["Ly_a"], "Ly_ab": ["Ly_a", "Ly_b"]}
#: The lines of a series whose AOD gives N and b (the Ly-b of an interloper
#: is too weak to help)
AOD_LINES = {"CIV": ["CIV_1548", "CIV_1550"], "Ly_a": ["Ly_a"], "Ly_ab": ["Ly_a"]}


def line_waves(line: Line) -> List[float]:
    """Observed wavelengths (Angstrom) of the lines of one component."""
    return [ATOM_DATA[n]["wave"] * (1 + line.z) for n in SERIES_LINES[line.series]]


def significance(comps: List[MockComponent], truth: MockTruth, transition: str) -> float:
    """EW / sigma_EW of the weaker true line (the rule of ``vm.is_detectable``)."""
    best = np.inf
    for c in comps:
        line = transition if c.series == "CIV" else "Ly_a"
        ew = vm.equivalent_width_ang(line, c.z, c.logN, c.b)
        fwhm = np.hypot(C_KMS / truth.resol, 1.665 * c.b)
        lam = ATOM_DATA[line]["wave"] * (1 + c.z)
        sigma_ew = (lam * truth.dv_kms / C_KMS * np.sqrt(max(1.0, 2 * fwhm / truth.dv_kms))
                    / truth.snr)
        best = min(best, ew / sigma_ew)
    return float(best)


# --- The hypotheses -------------------------------------------------------

def aod_profile(spec, lam_ang: float, half_kms: float, snr: float):
    """Velocity (km/s) and apparent optical depth around one line."""
    x = spec.x.to_value("nm") * 10.0
    v = C_KMS * np.log(x / lam_ang)
    sel = np.abs(v) <= half_kms
    flux = spec.y.value[sel] / spec.cont.value[sel]
    tau = -np.log(np.clip(flux, max(0.02, 1.0 / snr), None))
    return v[sel], tau


def aod_estimate(spec, line: Line, resol: float, snr: float) -> Tuple[float, float]:
    """logN and b of one component from the AOD of its lines: N the largest
    of the lines (the weaker one corrects the saturation of the stronger),
    b from the second moment of the stacked N(v), the instrumental width
    subtracted in quadrature."""
    fwhm = C_KMS / resol
    half = fwhm + AOD_HALF_KMS
    names = AOD_LINES[line.series]
    n_lines, stack, v_ref = [], None, None
    for name in names:
        lam = ATOM_DATA[name]["wave"] * (1 + line.z)
        v, tau = aod_profile(spec, lam, half, snr)
        nv = AOD_CONST / (ATOM_DATA[name]["f"] * ATOM_DATA[name]["wave"]) * tau
        n_lines.append(np.trapezoid(nv, v))
        if stack is None:
            v_ref, stack = v, nv
        else:
            stack = stack + np.interp(v_ref, v, nv, left=0.0, right=0.0)
    logn = float(np.log10(max(max(n_lines), 1e11)))
    weights = np.clip(stack, 0.0, None)
    if weights.sum() <= 0:
        return logn, 10.0
    mean = np.sum(weights * v_ref) / weights.sum()
    var = np.sum(weights * (v_ref - mean) ** 2) / weights.sum() - (fwhm / 2.3548) ** 2
    b_max = 40.0 if line.series == "CIV" else 60.0
    b = float(np.clip(np.sqrt(2 * max(var, 0.0)), 3.0, b_max))
    return logn, b


def hypothesis_systs(spec, lines: List[Line], resol: float, snr: float) -> SystemList:
    systs = SystemList(SystemListData(components=[]))
    for ln in lines:
        logn, b = aod_estimate(spec, ln, resol, snr)
        systs = systs.add_component(ln.series, ln.z, logN=logn, b=b)
    return systs


def pixel_mask(spec, hyps) -> np.ndarray:
    x = spec.x.to_value("nm") * 10.0
    mask = np.zeros(len(x), dtype=bool)
    for ln in hyps["true"] + hyps["false"]:
        for lam in line_waves(ln):
            mask |= np.abs(C_KMS * np.log(x / lam)) <= PIXEL_HALF_KMS
    return mask


def measure(spec, systs: SystemList, pixels: np.ndarray) -> Dict[str, Optional[float]]:
    ev = evaluate_model(spec, systs, pixels, min_influence=0.0)
    return {c: getattr(ev, c) for c in CRITERIA} | {"k": ev.k, "n": ev.n_pixels}


def short_fit(spec, systs: SystemList) -> SystemList:
    if not systs.components:
        return systs
    fitted, _, _ = VoigtFitter(spec, systs).fit(
        max_nfev=FIT_NFEV, b_max=FIT_B_MAX, logN_max=FIT_LOGN_MAX)
    return fitted


def run_trial(trial: Trial) -> Dict:
    comps, hyps, weak = draw_truth(trial)
    truth = make_truth(trial, comps, hyps)
    spec = generate(truth)
    pixels = pixel_mask(spec, hyps)
    out = {"scenario": trial.scenario, "snr": trial.snr, "resol": trial.resol,
           "seed": trial.seed, "sigma": significance(comps, truth, weak),
           "truth": [dataclasses.asdict(c) for c in comps]}
    t0 = time.perf_counter()
    for key in ("true", "false"):
        systs = hypothesis_systs(spec, hyps[key], trial.resol, trial.snr)
        out[f"{key}_aod"] = measure(spec, systs, pixels)
        out[f"{key}_fit"] = measure(spec, short_fit(spec, systs), pixels)
    out["seconds"] = time.perf_counter() - t0
    return out


# --- The summary ----------------------------------------------------------

def delta(rec: Dict, mode: str, crit: str) -> Optional[float]:
    """IC(false) - IC(true): positive when the criterion picks the truth."""
    t, f = rec[f"true_{mode}"][crit], rec[f"false_{mode}"][crit]
    return None if t is None or f is None else f - t


def lya_sign(scenario: str) -> int:
    """+1 when the true hypothesis is the Ly-a one, -1 when the false is,
    0 without a Ly-a alternative (the prior does not apply)."""
    return {"coincidence": 1, "doublet_vs_lya": -1, "doublet_vs_lyab": -1,
            "alias": -1}.get(scenario, 0)


def picks(records: List[Dict], mode: str, crit: str, odds: float = 1.0) -> Tuple[float, float, int]:
    """Fraction of trials where the criterion picks the truth, the median
    weight of the truth (1 / (1 + exp(-d/2))) and the number of trials."""
    ds = []
    for r in records:
        d = delta(r, mode, crit)
        if d is None:
            continue
        ds.append(d + 2 * np.log(odds) * lya_sign(r["scenario"]))
    if not ds:
        return np.nan, np.nan, 0
    arr = np.array(ds)
    weight = 1.0 / (1.0 + np.exp(-np.clip(arr, -700, 700) / 2))
    return float(np.mean(arr > 0)), float(np.median(weight)), len(arr)


def sig_bin(sigma: float) -> str:
    return SIG_LABELS[int(np.searchsorted(SIG_BINS, sigma, side="right")) - 1]


def table(records: List[Dict], mode: str, rows, row_of, title: str) -> str:
    head = f"| {title} | n | " + " | ".join(CRITERIA) + " |"
    lines = [head, "|" + "---|" * (len(CRITERIA) + 2)]
    for row in rows:
        sel = [r for r in records if row_of(r) == row]
        if not sel:
            continue
        cells = [f"{picks(sel, mode, c)[0]:.2f}" for c in CRITERIA]
        lines.append(f"| {row} | {len(sel)} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def summary(records: List[Dict]) -> str:
    out = []
    for scen in SCENARIOS:
        recs = [r for r in records if r["scenario"] == scen]
        for mode in MODES:
            out.append(f"\n### {scen}, {mode}: fraction of trials that pick the truth\n")
            out.append(table(recs, mode, SIG_LABELS, lambda r: sig_bin(r["sigma"]), "sigma"))
            out.append("")
            cells = [(s, R) for s in GRID["snr"] for R in GRID["resol"]]
            out.append(table(recs, mode, [f"S/N {s:g}, R {R:g}" for s, R in cells],
                             lambda r: f"S/N {r['snr']:g}, R {r['resol']:g}", "cell"))
    out.append("\n### Median weight of the truth, exp(-dIC/2) normalised, by scenario (fit)\n")
    out.append("| scenario | " + " | ".join(CRITERIA) + " |")
    out.append("|" + "---|" * (len(CRITERIA) + 1))
    for scen in SCENARIOS:
        recs = [r for r in records if r["scenario"] == scen and r["sigma"] >= 3]
        out.append(f"| {scen} | " + " | ".join(
            f"{picks(recs, 'fit', c)[1]:.2f}" for c in CRITERIA) + " |")
    out.append("\n### Prior odds Ly-a : CIV (fit, sigma >= 3), fraction that picks the truth\n")
    out.append("| scenario | odds | aic | bic | spic |")
    out.append("|---|---|---|---|---|")
    for scen in ("doublet_vs_lya", "doublet_vs_lyab", "coincidence", "alias"):
        recs = [r for r in records if r["scenario"] == scen and r["sigma"] >= 3]
        for odds in PRIOR_ODDS:
            vals = [f"{picks(recs, 'fit', c, odds)[0]:.2f}" for c in ("aic", "bic", "spic")]
            out.append(f"| {scen} | {odds:g} | " + " | ".join(vals) + " |")
    secs = np.array([r["seconds"] for r in records])
    out.append(f"\nTime per trial (two hypotheses, AOD and fit): median {np.median(secs):.2f} s,"
               f" p95 {np.percentile(secs, 95):.2f} s")
    return "\n".join(out)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--trials", type=int, default=40, help="trials per scenario and cell")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument("--scenarios", nargs="*", default=list(SCENARIOS))
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.WARNING)
    logging.getLogger("astrocook").setLevel(logging.ERROR)
    trials = [Trial(s, snr, R, seed) for s in args.scenarios for snr in GRID["snr"]
              for R in GRID["resol"] for seed in range(args.trials)]
    t0 = time.perf_counter()
    with Pool(args.workers) as pool:
        records = pool.map(run_trial, trials, chunksize=4)
    print(f"{len(records)} trials in {time.perf_counter() - t0:.0f} s")
    print(summary(records))
    if args.json:
        args.json.write_text(json.dumps(records, indent=1))


if __name__ == "__main__":
    main()
