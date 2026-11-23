"""The sweep of the anchored systems (phase 7, PR 7c.4d): the rarer and
single lines, and the H I of the metal systems.

The second pass of the identification by hypotheses
(:mod:`astrocook.fitting.hypotheses`, ``design/identification.md`` section
3.5). The anchors (the doublets and the pairs) find the systems; the other
lines of a system are too weak, too rare or single (CII 1334, SiIII 1206,
OI 1302, AlII 1670, the weaker FeII and SiII lines, the H I in the forest)
to be found by a pair. For each metal system accepted by the first pass,
and for each species of :data:`SWEEP_SPECIES` with a line in the coverage
other than that of the anchor, the hypothesis "Y at z":

- its components at the redshifts of the anchor's components (one for the
  H I, at the column-weighted redshift of the anchor: the H I of a system
  is broader than its metal components and one line of it is what the
  forest resolves), N from the apparent optical depth at its predicted
  lines (the least over its lines: that of a line blended with the forest
  is an upper limit), z and b free; the components without influence on the pixels
  (:func:`~astrocook.fitting.model_evaluation.evaluate_model`,
  ``min_influence`` 1) dropped after the first fit, the components that
  leave the profile of the anchor (:func:`_span`) make the hypothesis fail;
- compared with the labels its pixels already have: the free lines of the
  default labels (the assumed Ly-a, the unidentified lines) that fall on its
  predicted lines are replaced by it, the other default components of the
  regions within its span (its neighbours: a line of the forest blended
  with one of its lines) are free on both sides, the rest of the model
  (every accepted component of the first pass and of the sweep so far, the
  components of the session) is fixed. Both sides get the short fit on the
  same pixels, without the refinement by the residuals of the first pass
  (the neighbours free do its work here; with it a test of SiII in the
  forest took minutes), and Y is accepted when its criterion is not above
  that of the labels by more than :data:`SWEEP_TIE`;
- an accepted species becomes an edge of the graph (from the region of the
  anchor's strong line to the region of its strongest line), the regions it
  explains labelled with its lines, its components candidate systems with
  the weight of the comparison (exp(-dIC/2) normalised over the two).

What it does not do (``design/identification.md``): the upper limits where
nothing is there (F156, 7d.7; the sweep skips a species with no absorption
at its lines, :data:`SWEEP_SIGMA`); the systems of the Ly-a/Ly-b pairs of
the forest (their metals, the metals of the forest off by default:
decision 2 of the design); the tellurics (7c.5).
"""
from __future__ import annotations

import dataclasses
import logging
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

from astrocook.core import identification as ident
from astrocook.core.atomic_data import ATOM_DATA, series_lines
from astrocook.core.structures import ComponentData, SystemListData
from astrocook.core.system_list import SystemList
from astrocook.fitting import hypotheses as H
from astrocook.fitting.evaluation_window import EvaluationWindow

logger = logging.getLogger(__name__)

#: The species of the sweep, with the series of their model: the metals of
#: the anchors (each tested at the systems of the others: a SiIV too weak
#: for its pair at the redshift of a CIV), the single lines of the design
#: (section 3.5), the rarer ions of the damped systems, and the H I (Ly-a
#: and Ly-b). The series of CII and SiII are written as their ground-state
#: lines (a custom series, ``series_lines``): the standard multiplets also
#: hold the excited fine-structure lines CII*_1335 and SiII*_1264, which a
#: model would put at the column of the ground state (F173). The species of
#: the optical (NaI, CaII) and the lines below 912 A of the table left out:
#: never in the coverage of the forest and the red side of a UV spectrum, or
#: below the Lyman limit
SWEEP_SPECIES: Dict[str, str] = {
    "HI": "Ly_ab",
    "CIV": "CIV",
    "SiIV": "SiIV",
    "MgII": "MgII",
    "NV": "NV",
    "OVI": "OVI",
    "AlIII": "AlIII",
    "CII": "CII_1036,CII_1334",
    "CIII": "CIII_977",
    "SiII": "SiII_1190,SiII_1193,SiII_1260,SiII_1304,SiII_1526,SiII_1808",
    "SiIII": "SiIII_1206",
    "OI": "OI_1039,OI_1302",
    "AlII": "AlII_1670",
    "MgI": "MgI_2852",
    "FeII": "FeII",
    "NI": "NI",
    "SII": "SII",
    "NiII": "NiII",
    "CrII": "CrII",
    "ZnII": "ZnII",
    "MnII": "MnII",
}
#: The species of an anchor (``ANCHOR_PAIRS``, ``Ly_ab``) as a species of
#: the sweep: the anchor's own is not swept
_ANCHOR_SPECIES = {"Ly_ab": "HI"}
#: The least significance of the absorption at one of the predicted lines
#: of a species (the sum of 1 - F over the profile of the anchor, in sigma
#: of that sum) for it to be compared. Its purpose: the time; a species with
#: nothing at its lines cannot lower the criterion (its components only add
#: parameters), and its comparison is the upper limit of F156 (7d.7)
SWEEP_SIGMA = 3.0
#: The tie of the comparison, in the units of the criterion: "Y at z" is
#: accepted when its criterion is not above that of the labels by more than
#: this. Its purpose: where the H I of a system is the assumed Ly-a of a
#: line whose Ly-b is not covered, the two models are the same (one line,
#: the same three parameters) and their criteria differ by the noise of the
#: two fits (about 1e-6); the tie goes to the label with a known redshift,
#: the prior that makes a single line credible (design, section 3.5). 0.1:
#: below any difference that the weights tell apart (0.49 against 0.51)
SWEEP_TIE = 0.1
#: The least influence of a component of "Y at z" (the chi2 change of
#: removing it, as ``evaluate_model``) to be kept after its first fit
SWEEP_INFLUENCE = 1.0
#: The largest comparison, in the free parameters of the labels of its
#: pixels (the replaced components and the neighbours). Its purpose: the
#: time; the fits of both sides and their refinement cost about k times
#: FIT_NFEV models each (F171), and a test of SiII on a forest mock with 15
#: free lines took 5 minutes
SWEEP_MAX_K = 30
#: The start of the b of the H I of a system (km/s)
HI_B = 20.0


@dataclasses.dataclass
class SweepTest:
    """One hypothesis "Y at z" of the sweep, and its outcome."""
    group: int
    """The index of the group of the anchored system (in the results)."""
    anchor: int
    """The candidate (edge) of the system."""
    species: str
    series: str
    z: float
    """The column-weighted redshift of the anchor."""
    hyp: H._Hyp
    """Its components (Y, then its neighbours)."""
    labels: H._Hyp
    """The labels of its pixels (the replaced default components, refitted
    on the same pixels)."""
    accepted: bool = False
    weight: Optional[float] = None
    reason: str = ""
    regions: Tuple[int, ...] = ()
    """The regions its lines explain."""
    from_kinds: Tuple[str, ...] = ()
    """The roles of the default components it replaced (``"lya"``,
    ``"unidentified"``), one per component."""
    lines: Tuple[Tuple[int, str], ...] = ()
    """(region, transition) of the lines it explains."""


@dataclasses.dataclass
class SweepResult:
    """The sweep of a spectrum."""
    tests: List[SweepTest]
    replaced: Set[str]
    """The uuids of the default components replaced by an accepted species."""
    n_skipped: int
    """The species with a line in the coverage and nothing at their lines."""
    seconds: float = 0.0


@dataclasses.dataclass
class _System:
    group: int
    anchor: int
    species: str
    z: np.ndarray
    logn: np.ndarray
    b: np.ndarray
    region: int
    """The region of the strong line of the anchor."""
    line: str
    """The strong line of the anchor."""


class _Model:
    """The accepted components of the whole spectrum (the first pass and the
    sweep so far): their data by uuid, and the default ones with their
    role and region (the components the sweep can replace)."""

    def __init__(self, results: Sequence[H.GroupResult], fixed: Optional[SystemList]):
        self.comps: Dict[str, ComponentData] = {}
        self.defaults: Dict[str, H._Comp] = {}
        self.series: Dict[str, str] = {}
        self.session = list(fixed.components) if fixed is not None else []
        for g in results:
            h = g.hyps[g.accepted]
            last = h.fitted if h.fitted is not None else h.systs
            if last is None:
                continue
            data = {c.uuid: c for c in last.components}
            for comp, u in zip(h.comps, h.uuids, strict=True):
                if u not in data:
                    continue
                self.comps[u] = data[u]
                self.series[u] = comp.series
                if comp.role in ("lya", "unidentified"):
                    self.defaults[u] = comp

    def without(self, uuids: Set[str]) -> SystemList:
        """The fixed model: the session and every accepted component but
        ``uuids``."""
        comps = self.session + [c for u, c in self.comps.items() if u not in uuids]
        return SystemList(SystemListData(components=comps))

    def subset(self, uuids: Sequence[str]) -> SystemList:
        return SystemList(SystemListData(components=[self.comps[u] for u in uuids]))

    def replace(self, gone: Set[str], new: SystemList, comps: Sequence[H._Comp]) -> None:
        for u in gone:
            self.comps.pop(u, None)
            self.defaults.pop(u, None)
        for c, comp in zip(new.components, comps, strict=True):
            self.comps[c.uuid] = c
            self.series[c.uuid] = comp.series
            if comp.role in ("lya", "unidentified"):
                self.defaults[c.uuid] = comp


def _systems(results: Sequence[H.GroupResult], cands: Mapping[int, H._Candidate]
             ) -> List[_System]:
    """The metal systems of the accepted hypotheses: per anchor, its
    components (those of the refinement on its lines too)."""
    out = []
    for gi, g in enumerate(results):
        h = g.hyps[g.accepted]
        last = h.fitted if h.fitted is not None else h.systs
        if last is None or not h.anchors:
            continue
        data = {c.uuid: c for c in last.components}
        for a in h.anchors:
            c = cands[a]
            if c.species == H.LYAB:
                continue
            own = [data[u] for comp, u in zip(h.comps, h.uuids, strict=True)
                   if comp.role == "anchor" and comp.cand == a and u in data]
            if not own:
                continue
            out.append(_System(gi, a, c.species, np.array([x.z for x in own]),
                               np.array([x.logN for x in own]), np.array([x.b for x in own]),
                               c.regions[0], c.lines[0]))
    return out


def _span(s: _System, fwhm_kms: float) -> Tuple[float, float]:
    """The profile of the anchor, in km/s from its column-weighted
    redshift: from its first to its last component, widened by the larger
    of half the instrumental FWHM and the largest b of its components (the
    range of the z of "Y at z", and of the pixels of :data:`SWEEP_SIGMA`)."""
    zc = _zc(s)
    dv = H.C_KMS * np.log((1.0 + s.z) / (1.0 + zc))
    w = max(0.5 * fwhm_kms, float(np.max(s.b)))
    return float(dv.min() - w), float(dv.max() + w)


def _zc(s: _System) -> float:
    w = 10.0 ** (s.logn - s.logn.max())
    return float(np.sum(w * s.z) / np.sum(w))


def _covered(series: str, zc: float, v: np.ndarray, forest_lo: float,
             x_range: Optional[Tuple[float, float]] = None) -> List[str]:
    """The lines of ``series`` at ``zc`` within the spectrum, above the
    Lyman limit of the quasar and within ``x_range`` (nm), if given."""
    lo, hi = max(v[0], forest_lo), v[-1]
    if x_range is not None:
        lo, hi = max(lo, H.C_KMS * np.log(x_range[0])), min(hi, H.C_KMS * np.log(x_range[1]))
    out = []
    for t in series_lines(series):
        if t not in ATOM_DATA or ATOM_DATA[t]['f'] <= 0:
            continue
        vc = H.C_KMS * np.log(H._wave_nm(t) * (1.0 + zc))
        if lo < vc < hi:
            out.append(t)
    return out


def _significance(lines: Sequence[str], zc: float, span: Tuple[float, float],
                  v: np.ndarray, flux: np.ndarray, err: np.ndarray) -> float:
    """The largest significance of the absorption at the lines (the sum of
    1 - F over the span in sigma of that sum)."""
    best = 0.0
    for t in lines:
        vc = H.C_KMS * np.log(H._wave_nm(t) * (1.0 + zc))
        lo, hi = np.searchsorted(v, [vc + span[0], vc + span[1]])
        f, e = flux[lo:hi], err[lo:hi]
        ok = np.isfinite(f) & np.isfinite(e) & (e > 0)
        if ok.sum() < 2:
            continue
        best = max(best, float(np.sum(1.0 - f[ok]) / np.sqrt(np.sum(e[ok] ** 2))))
    return best


def _touched(lines: Sequence[str], zc: float, span: Tuple[float, float], v: np.ndarray,
             labels: np.ndarray) -> Set[int]:
    """The regions within the span of the predicted lines."""
    out: Set[int] = set()
    for t in lines:
        vc = H.C_KMS * np.log(H._wave_nm(t) * (1.0 + zc))
        lo, hi = np.searchsorted(v, [vc + span[0], vc + span[1]])
        out |= {int(r) for r in np.unique(labels[lo:hi + 1]) if r > 0}
    return out


def _near(model: _Model, s: _System, species: str, lines: Sequence[str], zc: float,
          span: Tuple[float, float], v: np.ndarray, fwhm: np.ndarray, labels: np.ndarray
          ) -> Tuple[List[str], List[str]]:
    """The default components at the predicted lines of the species: those
    on the core of one of its components (within the larger of half the
    instrumental FWHM and the b of the anchor's component; for the H I the
    span, :data:`HI_B` wider), which it replaces, and the others of the regions within its
    span, its neighbours, free in both sides of the comparison (a line of
    the forest blended with a line of the species; the components of a
    group without pairs keep their AOD parameters, and fixed they left a
    misfit of thousands in the chi2 of a region)."""
    if species == "HI":
        # One H I line for the whole system: its core the profile of the
        # anchor, HI_B wider
        cores = [(0.5 * (span[0] + span[1]), 0.5 * (span[1] - span[0]) + HI_B)]
    else:
        cores = list(zip(H.C_KMS * np.log((1.0 + s.z) / (1.0 + zc)), s.b, strict=True))
    touched = _touched(lines, zc, span, v, labels)
    replaced, neighbours = [], []
    for u, comp in model.defaults.items():
        vu = H.C_KMS * np.log(H._wave_nm(H.PLACEHOLDER) * (1.0 + model.comps[u].z))
        on_core, in_span = False, comp.region in touched
        for t in lines:
            vc = H.C_KMS * np.log(H._wave_nm(t) * (1.0 + zc))
            i = min(max(int(np.searchsorted(v, vc)), 0), len(v) - 1)
            half = 0.5 * fwhm[i]
            on_core |= any(abs(vu - vc - dv) <= max(half, b) for dv, b in cores)
            in_span |= vc + span[0] - half <= vu <= vc + span[1] + half
        if on_core:
            replaced.append(u)
        elif in_span:
            neighbours.append(u)
    return replaced, neighbours


def _start(s: _System, species: str, series: str, lines: Sequence[str], zc: float,
           span: Tuple[float, float], v, flux, err, fwhm_kms: float
           ) -> Tuple[SystemList, List[H._Comp]]:
    """The components of "Y at z" with their AOD guesses."""
    systs = H._empty()
    comps: List[H._Comp] = []
    if species == "HI":
        zs, ranges = [zc], [span]
    else:
        dv = H.C_KMS * np.log((1.0 + s.z) / (1.0 + zc))
        order = np.argsort(dv)
        zs, ranges = [], []
        for j, k in enumerate(order):
            lo = span[0] if j == 0 else 0.5 * (dv[order[j - 1]] + dv[k])
            hi = span[1] if j == len(order) - 1 else 0.5 * (dv[k] + dv[order[j + 1]])
            zs.append(float(s.z[k]))
            ranges.append((lo - dv[k], hi - dv[k]))
    for z, rng in zip(zs, ranges, strict=True):
        # The least column over the lines: the AOD of a line blended with the
        # forest is an upper limit (the largest, as for the anchors, started
        # the SiII of a system at the column of a forest line on its 1260)
        guesses = [H.aod_guess(v, flux, err, [t], z, rng, fwhm_kms) for t in lines]
        logn, b = min(guesses)
        if species == "HI":
            b = max(b, HI_B)
        systs = systs.add_component(series, z, logN=logn, b=b)
        comps.append(H._Comp(series, z, -1, -1, "sweep", None, rng))
    return systs, comps


def _window(ctx: H._Context, model: _Model, systs: SystemList, uuids: Sequence[str],
            lines: Sequence[str], zc: float, span: Tuple[float, float], v: np.ndarray
            ) -> EvaluationWindow:
    """The window of a comparison: that of its lines and of the labels it
    replaces (:meth:`EvaluationWindow.around`), without the pixels of the
    regions that neither touches (as the window of a group of the first
    pass, ``hypotheses._window``: the components of a group without pairs
    keep their AOD parameters, and the misfit of their pixels, the same for
    both sides, took the fits of the comparison).

    Raises
    ------
    ValueError
        If no pixel is left.
    """
    win = EvaluationWindow.around(ctx.spec, [systs, model.subset(uuids)],
                                  fixed=model.without(set(uuids)))
    ids = {model.defaults[u].region for u in uuids} | _touched(lines, zc, span, v, ctx.labels)
    keep = win.mask & ((ctx.labels == 0) | np.isin(ctx.labels, sorted(ids)))
    if keep.sum() == win.index.size:
        return win
    return EvaluationWindow.from_mask(ctx.spec, keep, fixed=model.without(set(uuids)))


def _fit_labels(ctx: H._Context, win: EvaluationWindow, model: _Model, uuids: Sequence[str],
                cands: Mapping[int, H._Candidate], v: np.ndarray, extra: int) -> H._Hyp:
    """The labels of the pixels (the replaced components and the
    neighbours), refitted on the window and refined by its residuals, with
    ``extra`` more components than the refinement of "Y at z" (its
    components beyond the default lines it replaces): the null hypothesis of
    the test is that its pixels need as many free lines as Y puts there. Its purpose: the regions of the
    groups without pairs keep the components of the de-blending, often fewer
    than their lines, and any species with a line on them won by the misfit
    alone (an H I where the mock has none, by 66 in the SpIC)."""
    h = H._Hyp((), [model.defaults[u] for u in uuids])
    h.systs = model.subset(uuids)
    h.uuids = list(uuids)
    if uuids:
        H._fit(win, h)
    else:
        h.fitted, h.fit = h.systs, H._measure(win, h.systs)
    H._refine(ctx, win, h, cands, v, extra)
    return h


def _drop(win: EvaluationWindow, h: H._Hyp, gone: Set[str]) -> None:
    """The components ``gone`` left out of ``h``, the rest fitted again."""
    assert h.fitted is not None
    data = {c.uuid: c for c in h.fitted.components}
    idx = [i for i, u in enumerate(h.uuids) if u not in gone]
    h.comps = [h.comps[i] for i in idx]
    h.uuids = [h.uuids[i] for i in idx]
    h.systs = SystemList(SystemListData(components=[data[u] for u in h.uuids]))
    H._fit(win, h)


def _own(h: H._Hyp) -> List[str]:
    """The uuids of the components of the species."""
    return [u for c, u in zip(h.comps, h.uuids, strict=True) if c.role == "sweep"]


def _outside(h: H._Hyp, zc: float, span: Tuple[float, float]) -> Set[str]:
    """The components of the species out of the profile of the anchor."""
    assert h.fitted is not None
    data = {c.uuid: c for c in h.fitted.components}
    return {u for u in _own(h)
            if not span[0] <= H.C_KMS * np.log((1.0 + data[u].z) / (1.0 + zc)) <= span[1]}


def _fit_species(win: EvaluationWindow, h: H._Hyp, zc: float, span: Tuple[float, float]
                 ) -> str:
    """The short fit of "Y at z": its components that leave the profile of
    the anchor, then those without influence (:data:`SWEEP_INFLUENCE`),
    dropped and the rest fitted again. An empty string, or why it failed."""
    H._fit(win, h)
    for step in ("out", "influence"):
        if h.fit is None:
            return h.reason
        if step == "out":
            gone = _outside(h, zc, span)
        else:
            assert h.fitted is not None
            counted = set(win.evaluate(h.fitted, min_influence=SWEEP_INFLUENCE).components)
            gone = {u for u in _own(h) if u not in counted}
        if len(gone) == len(_own(h)):
            return ("its components out of the profile of the anchor" if step == "out"
                    else "no component with an influence on the pixels")
        if gone:
            _drop(win, h, gone)
    if h.fit is None:
        return h.reason
    return "" if not _outside(h, zc, span) else "its components out of the profile of the anchor"


def _explained(h: H._Hyp, lines: Sequence[str], labels: np.ndarray, v: np.ndarray
               ) -> List[Tuple[int, str]]:
    """(region, line) of the lines of the species at its fitted components
    that fall in a region."""
    assert h.fitted is not None
    data = {c.uuid: c for c in h.fitted.components}
    out = []
    for comp, u in zip(h.comps, h.uuids, strict=True):
        if comp.role != "sweep":
            continue
        for t in lines:
            vc = H.C_KMS * np.log(H._wave_nm(t) * (1.0 + data[u].z))
            i = int(np.clip(np.searchsorted(v, vc), 0, len(v) - 1))
            if labels[i] > 0 and (int(labels[i]), t) not in out:
                out.append((int(labels[i]), t))
    return out


def _test(ctx: H._Context, model: _Model, s: _System, species: str, cands, v, flux, err,
          fwhm) -> Optional[SweepTest]:
    """The comparison of "Y at z" with the labels of its pixels; None when
    the species is not covered or has nothing at its lines."""
    series = SWEEP_SPECIES[species]
    zc = _zc(s)
    i = int(np.clip(np.searchsorted(v, H.C_KMS * np.log(H._wave_nm(s.line) * (1 + zc))),
                    0, len(v) - 1))
    fwhm_kms = float(fwhm[i])
    span = _span(s, fwhm_kms)
    lines = _covered(series, zc, v, ctx.forest[0], ctx.x_range)
    if not lines:
        return None
    if _significance(lines, zc, span, v, flux, err) < SWEEP_SIGMA:
        return SweepTest(s.group, s.anchor, species, series, zc, H._Hyp((), []),
                         H._Hyp((), []), reason="nothing at its lines")
    gone, near = _near(model, s, species, lines, zc, span, v, fwhm, ctx.labels)
    k = 3 * (len(gone) + len(near))
    if k > SWEEP_MAX_K:
        # As a group of the first pass (F171): the fits of the comparison
        # would take minutes
        return SweepTest(s.group, s.anchor, species, series, zc, H._Hyp((), []),
                         H._Hyp((), []), reason=f"not compared: the labels of its pixels "
                         f"have k = {k} free parameters, above {SWEEP_MAX_K}")
    systs, comps = _start(s, species, series, lines, zc, span, v, flux, err, fwhm_kms)
    both = SystemList(SystemListData(
        components=list(systs.components) + list(model.subset(near).components)))
    hyp = H._Hyp((), comps + [model.defaults[u] for u in near], systs=both,
                 uuids=[c.uuid for c in both.components])
    test = SweepTest(s.group, s.anchor, species, series, zc, hyp, H._Hyp((), []))
    try:
        win = _window(ctx, model, both, gone + near, lines, zc, span, v)
    except ValueError as e:
        test.reason = f"not evaluated: {e}"
        return test
    why = _fit_species(win, hyp, zc, span)
    if why:
        test.reason = why
        return test
    test.labels = _fit_labels(ctx, win, model, gone + near, cands, v,
                              max(0, len(_own(hyp)) - len(gone)))
    H._refine(ctx, win, hyp, cands, v)
    return _decide(ctx, model, test, gone, near, lines, v)


def _decide(ctx: H._Context, model: _Model, test: SweepTest, gone: Sequence[str],
            near: Sequence[str], lines: Sequence[str], v: np.ndarray) -> SweepTest:
    """The weight of "Y at z" against the labels, and the model updated
    when it is accepted."""
    crit = ctx.criterion
    h, lab = test.hyp, test.labels
    assert h.fit is not None and h.fitted is not None
    if lab.fit is None:
        test.reason = f"the labels could not be fitted: {lab.reason}"
        return test
    d = h.fit[crit] - lab.fit[crit]
    test.weight = float(1.0 / (1.0 + np.exp(np.clip(0.5 * d, -700, 700))))
    test.accepted = d <= SWEEP_TIE
    names = {"lya": "assumed Ly-a", "unidentified": "unidentified"}
    what = ", ".join(f"{sum(1 for u in gone if model.defaults[u].role == r)} {n}"
                     for r, n in names.items()
                     if any(model.defaults[u].role == r for u in gone)) or "no line"
    test.reason = f"{crit} {d:+.1f} against the labels of its pixels ({what})"
    if not test.accepted:
        return test
    test.from_kinds = tuple(model.defaults[u].role for u in gone)
    test.lines = tuple(_explained(h, lines, ctx.labels, v))
    test.regions = tuple(dict.fromkeys(r for r, _ in test.lines))
    data = {c.uuid: c for c in h.fitted.components}
    model.replace(set(gone) | set(near),
                  SystemList(SystemListData(components=[data[u] for u in h.uuids])), h.comps)
    return test


def sweep(ctx: H._Context, results: Sequence[H.GroupResult], cands: Mapping[int, H._Candidate],
          v: np.ndarray, flux: np.ndarray, err: np.ndarray, fwhm: np.ndarray) -> SweepResult:
    """The sweep of the metal systems of ``results`` (see the module)."""
    import time
    t0 = time.perf_counter()
    model = _Model(results, ctx.fixed)
    first = set(model.defaults)
    tests: List[SweepTest] = []
    skipped = 0
    for s in _systems(results, cands):
        own = _ANCHOR_SPECIES.get(s.species, s.species)
        for species in SWEEP_SPECIES:
            if species == own:
                continue
            t = _test(ctx, model, s, species, cands, v, flux, err, fwhm)
            if t is None:
                continue
            if t.reason == "nothing at its lines":
                skipped += 1
                continue
            tests.append(t)
    replaced = first - set(model.defaults)
    seconds = time.perf_counter() - t0
    n_acc = sum(t.accepted for t in tests)
    logger.info(f"Sweep of the anchored systems: {len(tests)} species compared "
                f"({skipped} with nothing at their lines), {n_acc} accepted, "
                f"{len(replaced)} default lines replaced, in {seconds:.1f} s.")
    return SweepResult(tests, replaced, skipped, seconds)


def structure(test: SweepTest, cands: Mapping[int, H._Candidate], gid: int,
              anchored_on: int, comps: Sequence[ident.Component]) -> ident.Group:
    """The group of the structure of a test: the regions of the anchor and
    those the species explains, an edge when it is accepted and explains a
    region, its hypothesis with the weight against the labels."""
    c = cands[test.anchor]
    edges: Tuple[ident.Edge, ...] = ()
    anchors: Tuple[int, ...] = ()
    if test.accepted and test.lines:
        rid, t = max(test.lines, key=lambda x: H._flam(x[1]))
        lag = float(H.C_KMS * np.log(ATOM_DATA[t]['wave'] / ATOM_DATA[c.lines[0]]['wave']))
        edges = (ident.Edge(regions=(c.regions[0], rid), transitions=(c.lines[0], t),
                            series=test.series, z=float(test.z), lag_kms=lag),)
        anchors = (0,)
    h = test.hyp
    stats = h.fit or {}
    criteria = {k: val for k, val in stats.items() if k not in ("k", "n", "chi2", "chi2_reduced")}
    if test.labels.fit is not None:
        criteria.update({f"labels_{k}": val for k, val in test.labels.fit.items()
                         if k in ("chi2", "spic", "bic", "aic")})
    # The default lines it replaced, by their role (an accepted test)
    for role in ("lya", "unidentified"):
        criteria[f"replaced_{role}"] = float(sum(r == role for r in test.from_kinds))
    regions = tuple(dict.fromkeys((c.regions[0],) + test.regions))
    status = (ident.Status.ACCEPTED if test.accepted else
              ident.Status.ALTERNATIVE if test.weight is not None else ident.Status.DROPPED)
    hyp = ident.Hypothesis(
        label=f"{test.species} at z = {test.z:.5f} (sweep of {c.species})",
        regions=regions, components=tuple(comps),
        stage=ident.Stage.FIT if h.fit is not None else ident.Stage.AOD,
        status=status, reason=test.reason,
        weight=test.weight if test.weight is not None else 0.0,
        chi2=stats.get("chi2"), k=H._int(stats.get("k")), n=H._int(stats.get("n")),
        criteria=criteria, anchors=anchors)
    return ident.Group(id=gid, regions=regions, edges=edges, hypotheses=(hyp,),
                       anchored_on=anchored_on)
