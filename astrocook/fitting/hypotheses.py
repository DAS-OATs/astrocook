"""Identification as a comparison between hypotheses (phase 7, PR 7c.4c).

The prototype of ``design/identification.md`` (sections 3, 5, 6, 10-12),
selected by ``identify_lines(method="hypotheses")`` next to the current
method (``"doublets"``, :meth:`~astrocook.core.spectrum.Spectrum.identify_lines`):

1. *Regions*: the lines found by the matched filters of PR 7c.4b
   (:func:`~astrocook.core.line_detection.detect_lines`), not the pixels that
   the clipping of the continuum rejects; each region is split into its
   peaks (the de-blending of the current method, on the normalised flux).
2. *Edges*: two regions that a transition pair of a requested species puts
   at the same redshift (the anchors of :data:`ANCHOR_PAIRS`: the doublets,
   two pairs of FeII and of SiII, and Ly-a/Ly-b in the forest), kept when
   the R^2 of F_1 = F_2^gamma (the score of the current method, now a cheap
   pre-filter) is above ``prefilter``. The metals only on the red side of
   the Ly-a emission (the metals of the forest off: decision 2 of the
   design), the Ly-a/Ly-b pairs only in the forest, no system beyond
   ``z_em`` by more than the margin of F55.
3. *Groups*: the connected components of the graph; a group of more than
   ``max_group`` regions is split at its weakest edges (lowest R^2).
4. *Hypotheses*: per group, each edge (a candidate) as the anchor of one
   hypothesis, the other peaks of the group with their default label, and
   the hypothesis "all default": n + 1 for n candidates. The default label
   is the Ly-a assumed in the forest when the Ly-a is asked (``Ly_a`` or
   ``Ly_ab`` among the species, F51), the unidentified line otherwise: one
   free line per peak, the Voigt profile of a placeholder transition (the
   Ly-a at its wavelength, k = 3). Two candidates that win alone against
   "all default" and share no region are combined in one hypothesis.
5. *Evaluation*: on an :class:`EvaluationWindow` of the group (PR 7c.4a),
   the same pixels for all its hypotheses, every component counted
   (``min_influence = 0``). First with N and b from the apparent optical
   depth (AOD); the hypotheses behind the best by more than ``screen`` are
   dropped; the others (the finalists) get a short fit (:data:`FIT_NFEV`
   evaluations) and are evaluated again.
6. *Credibility*: w_i = exp(-dIC_i / 2) normalised over the finalists of
   the group, IC the SpIC_H by default (decision 1 of the design), BIC and
   the others reported; optionally tempered by the chi2 per degree of
   freedom of the accepted model (``temper``, section 3.3: off by default).

The result has the same form as that of the current method, so that the
rest of the chain works unchanged: the ``abs_ids`` column (here every
region, the unidentified ones too), the maps of the candidate systems (one
per component of the accepted hypotheses that is not an unidentified
line; their score the weight of the hypothesis, 0 for the assumed Ly-a as
in the current method) with the parameters of the short fit as their
guesses, the labels of the regions, the dropped candidates, and the
structure of the identification (``core/identification.py``, PR 7c.3) with
every group, its edges and all its hypotheses.

The higher orders of the Lyman series are not used (design, section 12);
the assumed Ly-a has no Ly-b in its model (a Ly-b is modelled only by a
Ly-a/Ly-b hypothesis, anchored by its edge); the tellurics (7c.5) and the
prior of the incidence are not in the prototype.

7. *Sweep* (PR 7c.4d, ``sweep``, the H I alone by default, the other
   species with ``sweep="all"``): the H I and the other species tested at
   the redshift of each accepted metal system, against the labels of their
   pixels (:mod:`astrocook.fitting.hypotheses_sweep`); each
   test a group of its own in the structure (``anchored_on``).
"""
from __future__ import annotations

import dataclasses
import logging
import time
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Optional, Sequence, Tuple

import astropy.units as au
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths

from astrocook.core import identification as ident
from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS, xem_d
from astrocook.core.line_detection import detect_lines
from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import SystemListData
from astrocook.core.system_list import SystemList
from astrocook.fitting.evaluation_window import EvaluationWindow
from astrocook.fitting.voigt_fitter import NoPixelsToFit

if TYPE_CHECKING:
    from astrocook.fitting.hypotheses_sweep import SweepResult

logger = logging.getLogger(__name__)

C_KMS = 299792.458

#: The anchor pairs of each species: the two lines most likely to be
#: detected together (design, section 8): the doublets of the table, and
#: two pairs each for FeII and SiII (their strongest lines that are not
#: blended with each other). An edge of the graph is a pair of regions at
#: the lag of a pair. The pair of CII is CII_1334/CII_1036 since F173 (PR
#: 7c.4d-ter; before, CII_1334 and the excited CII*_1335)
ANCHOR_PAIRS: Dict[str, Tuple[Tuple[str, str], ...]] = {
    **{name: ((lines[0], lines[1]),) for name, lines in STANDARD_MULTIPLETS.items()
       if len(lines) == 2 and not name.startswith("Ly")},
    "FeII": (("FeII_2382", "FeII_2600"), ("FeII_2344", "FeII_2586")),
    "SiII": (("SiII_1260", "SiII_1526"), ("SiII_1304", "SiII_1526")),
}
#: The metal anchors of ``"all"`` (the default of the design, section 12)
ALL_METALS: Tuple[str, ...] = tuple(ANCHOR_PAIRS)
#: The Ly-a/Ly-b pair, an anchor in the forest
LYAB = "Ly_ab"
LYAB_PAIR = ("Ly_a", "Ly_b")

#: The pre-filter of the edges: the R^2 of F_1 = F_2^gamma of the pair (the
#: score of the current method) at or above which a candidate becomes an
#: edge of the graph. Its purpose: in the forest, and on the red side with
#: the false regions of the bias of the continuum (F168), every region is at
#: the lag of some other; without a pre-filter the groups percolate, every
#: coincidence costs a fit and many win by chance. Measured in PR 7c.4c (18
#: mocks of the grid full, no screen): at 0 the purity of CIV and MgII fell
#: to 0.56 and 0.52 (ECE 0.20, 36 s per spectrum, 10 minutes at R 140000);
#: at 0.2 and 0.3 0.96 and 0.91-0.94, with the same completeness (ECE 0.02-
#: 0.03, 4 s). 0.3 (chosen with Guido's leave, 19 November 2026). Measured
#: again in PR 7c.4c-bis with the score on the peaks (F169), the regions
#: split at their gaps and the refinement without a fixed number (the same
#: 18 mocks): at 0 the purity of CIV, SiIV, MgII 0.50, 0.79, 0.39 (ECE 0.23,
#: 262 s per spectrum, median); at 0.2 and 0.3 0.98, 1.00, 0.95, the same
#: hypotheses (ECE 0.044 and 0.050, 9.4 and 8.8 s). 0.3 (Guido, 20 November
#: 2026)
PREFILTER_R2 = 0.3
#: The screen: the hypotheses behind the best by more than this (in the
#: criterion, evaluated with the AOD parameters) are dropped before the
#: short fit; None: no screen, every hypothesis fitted. Its purpose: the
#: time (a fit per hypothesis). Measured in PR 7c.4c at 10, 20 and 40: the
#: AOD model of a strong line, close in logN and b to the fit (MgII logN
#: 13.42, b 6.4 against 13.50, 8.0), has a chi2 of thousands at S/N 100
#: (7645 against 88), so every threshold of a few tens dropped about half
#: of the true CIV and SiIV systems (completeness 0.44 against 0.77, ECE
#: 0.31 against 0.03) for 2.5 s per spectrum saved. None (chosen with
#: Guido's leave, 19 November 2026)
SCREEN_DELTA: Optional[float] = None
#: The largest group, in regions (decision 5 of the design): above it the
#: group is split at its weakest edges, so that the n + 1 hypotheses stay
#: few and the windows short
MAX_GROUP = 6
#: The criteria of the weights: SpIC_H (default), BIC, AIC
CRITERIA = ("spic", "bic", "aic")
#: The short fit of the finalists (the experiment of the design, section 4:
#: 300 evaluations; the bounds those of the experiment, wider in b than
#: the 50 km/s of the pipelines, as the unidentified lines can be blends)
FIT_NFEV = 300
FIT_B_MAX = 100.0
FIT_LOGN_MAX = 19.0
#: The largest hypothesis that gets the short fit, in free parameters (40
#: components of three). Its purpose: the time. A group with a hypothesis
#: above it is not compared: its default labels are accepted with their AOD
#: parameters, as a group without edges. 60 at first (Guido, 19 November
#: 2026): the Jacobian was by finite differences, one model per parameter,
#: and in the trough of a DLA (mock "combined" of the grid robustness at
#: S/N 10) hypotheses of 94 components, k = 281, kept one mock running for
#: more than 100 minutes. 120 since the analytic Jacobian (Guido, 25
#: November 2026, PR 7c.4d-bis): a short fit of k = 120 in a DLA trough
#: takes 6 s (107 s before), k = 282 18 s; on the grid robustness 120
#: compares one hypothesis more than 60, with the same numbers, at 13% more
#: time to identify; without a cap the grid did not end in 2 hours (a few
#: groups far above 120). The room is for the joint hypotheses of F174
MAX_FIT_K = 120
#: The refinement of the finalists by their residuals: after its short fit
#: a hypothesis gets a component where the residuals of its model are most
#: negative (a 3-pixel sum below -REFINE_SIGMA): one of its species when the
#: pixel is in a region of one of its anchors, an unidentified line (or the
#: assumed Ly-a) elsewhere, kept when it lowers the criterion. Its purpose:
#: the components of a hypothesis are the peaks of the de-blending, fewer
#: than the true ones in a saturated or blended system; with too few, a
#: metal hypothesis cannot fit both lines with one ratio, and the free lines
#: of "all default" win on the misfit (on the mocks of the grid the CIV,
#: SiIV, MgII systems of 3-4 components lost by 20-60 in the SpIC). Every
#: finalist gets it, so that each identification is compared with its best
#: model. An option of the recipe (``refine``), on by default (F170, Guido,
#: 20 November 2026)
#: The most components the refinement adds to a hypothesis (it also stops
#: when the next component does not lower the criterion of the group, or at
#: :data:`MAX_FIT_K`; None: no fixed number). 3, as in 7c.4c (Guido, 21
#: November 2026): on the grid robustness no fixed number gave the same
#: completeness, a lower purity (SiIV 0.71 against 0.75, MgII 0.70 against
#: 0.76) and 2.6 times the time of the identification (PR 7c.4c-bis)
MAX_ADDED: Optional[int] = 3
#: Where the refinement looks: the most negative 3-pixel sum of the
#: residuals, in sigma
REFINE_SIGMA = 3.0
#: The split of the regions of the matched filters at their gaps (F169, PR
#: 7c.4c-bis): two adjacent peaks of a region are in two regions of the
#: identification when the smoothed flux between them comes back to within
#: GAP_SIGMA times the median error from the continuum. Its purpose: the
#: regions of the filters are wide (F168), and in the forest they join
#: lines that the flux separates (at z_em 3.5, R 45000, S/N 30 a region of
#: 8500 km/s and 49 peaks, back to the continuum between 190 of 318 pairs
#: of adjacent peaks): a pair of regions is then a pair of blends, its
#: group has hypotheses above MAX_FIT_K, and one candidate stands for many
#: systems
SPLIT_GAPS = True
GAP_SIGMA = 2.0
#: The start of an added component
ADDED_LOGN = 12.5
ADDED_B = 10.0
#: The placeholder transition of an unidentified line (and of the assumed
#: Ly-a): one line, k = 3
PLACEHOLDER = "Ly_a"
#: AOD: N = AOD_CONST / (f lambda[A]) int tau dv[km/s] (Savage & Sembach 1991)
AOD_CONST = 3.768e14
#: The lowest normalised flux of the AOD (a saturated pixel is not infinite)
AOD_FLOOR = 0.02
#: The bounds of the AOD guesses
AOD_LOGN = (11.0, 19.0)
AOD_B = (3.0, 60.0)


# --- The pieces ------------------------------------------------------------

@dataclasses.dataclass(frozen=True)
class _Peak:
    index: int
    """The pixel of the peak, in the parent spectrum."""
    width_kms: float
    """Its full width at half its depth."""


@dataclasses.dataclass(frozen=True)
class _Region:
    id: int
    lo: int
    hi: int
    """The pixels lo..hi (inclusive) of the parent spectrum."""
    zone: str
    """``"forest"`` (Lyman limit to Ly-a), ``"red"`` (redward of Ly-a) or
    ``"blue"`` (below the Lyman limit)."""
    significance: float
    peaks: Tuple[_Peak, ...]


@dataclasses.dataclass(frozen=True)
class _Candidate:
    """An edge: a pair of regions at the lag of an anchor pair."""
    index: int
    species: str
    lines: Tuple[str, str]
    """(strong, partner): the line with the larger f lambda first."""
    regions: Tuple[int, int]
    """The regions of the strong line and of the partner."""
    z: float
    score: float


@dataclasses.dataclass(frozen=True)
class _Comp:
    """A component of a hypothesis, before its parameters."""
    series: str
    z: float
    region: int
    peak: int
    role: str
    """``"anchor"`` (of the candidate ``cand``), ``"lya"`` (the assumed
    Ly-a) or ``"unidentified"``."""
    cand: Optional[int] = None
    aod_range: Tuple[float, float] = (-30.0, 30.0)
    """The range of its AOD, in km/s from the peak."""
    width_kms: float = 20.0
    """The FWHM of its peak."""
    fwhm_kms: float = 0.0
    """The instrumental FWHM at its peak."""


@dataclasses.dataclass
class _Hyp:
    anchors: Tuple[int, ...]
    comps: List[_Comp]
    explained: frozenset = frozenset()
    """The peaks (region, index) its anchors explain."""
    uuids: List[str] = dataclasses.field(default_factory=list)
    """The uuid of the component of each of ``comps``."""
    systs: Optional[SystemList] = None
    fitted: Optional[SystemList] = None
    aod: Optional[Dict[str, float]] = None
    fit: Optional[Dict[str, float]] = None
    status: str = ident.Status.ALTERNATIVE
    reason: str = ""
    weight: Optional[float] = None


def _empty() -> SystemList:
    return SystemList(SystemListData(components=[]))


def _flam(t: str) -> float:
    return float(ATOM_DATA[t]['f'] * ATOM_DATA[t]['wave'])


def _wave_nm(t: str) -> float:
    return float(ATOM_DATA[t]['wave']) / 10.0


def anchor_pairs(species: Sequence[str]) -> List[Tuple[str, Tuple[str, str]]]:
    """(species, (strong line, partner)) of every anchor pair of ``species``."""
    out = []
    for s in species:
        pairs = (LYAB_PAIR,) if s == LYAB else ANCHOR_PAIRS.get(s, ())
        for a, b in pairs:
            out.append((s, (a, b) if _flam(a) >= _flam(b) else (b, a)))
    return out


def parse_species(multiplets: Sequence[str]) -> Tuple[List[str], bool]:
    """The anchors and whether the Ly-a is asked, from the species of the
    recipe: ``"all"`` is every metal anchor; ``Ly_ab`` asks the Ly-a and
    its pairs with the Ly-b, ``Ly_a`` the Ly-a alone (the assumed Ly-a of
    the forest).

    Raises
    ------
    ValueError
        For a species without an anchor.
    """
    anchors: List[str] = []
    lya = False
    for m in multiplets:
        if m == "all":
            anchors += [s for s in ALL_METALS if s not in anchors]
        elif m in ("Ly_a", LYAB):
            lya = True
            if m == LYAB and LYAB not in anchors:
                anchors.append(LYAB)
        elif m in ANCHOR_PAIRS:
            if m not in anchors:
                anchors.append(m)
        else:
            raise ValueError(f"No anchor for '{m}': the species searched are "
                             f"{', '.join(ALL_METALS)}, Ly_a and Ly_ab (or 'all').")
    return anchors, lya


# --- Regions and peaks -----------------------------------------------------

def _peaks(v: np.ndarray, flux: np.ndarray, err: np.ndarray, lo: int, hi: int
           ) -> Tuple[_Peak, ...]:
    """The peaks of the absorption of pixels lo..hi: the de-blending of the
    current method (``Spectrum.identify_lines``), on the normalised flux."""
    seg = 1.0 - flux[lo:hi + 1]
    seg = np.where(np.isfinite(seg), seg, 0.0)
    dv = np.abs(np.gradient(v[lo:hi + 1])) if hi > lo else np.array([1.0])
    if len(seg) < 3:
        return (_Peak(lo + int(np.argmax(seg)), float(np.sum(dv))),)
    smooth = gaussian_filter1d(seg, sigma=1.0)
    e = err[lo:hi + 1]
    e = e[np.isfinite(e) & (e > 0)]
    prominence = max(0.02, 3.0 * float(np.median(e))) if e.size else 0.02
    idx, _ = find_peaks(smooth, prominence=prominence, width=1.0)
    if len(idx) == 0:
        return (_Peak(lo + int(np.argmax(smooth)), float(np.sum(dv))),)
    widths = peak_widths(smooth, idx, rel_height=0.5)[0]
    return tuple(_Peak(lo + int(i), float(w * dv[i])) for i, w in zip(idx, widths, strict=True))


def _gaps(v: np.ndarray, flux: np.ndarray, err: np.ndarray, lo: int, hi: int,
          peaks: Sequence[_Peak]) -> List[int]:
    """The pixels of lo..hi where the region is split (F169, PR 7c.4c-bis):
    between two adjacent peaks, the highest pixel of the smoothed flux when
    it comes back to within :data:`GAP_SIGMA` times the median error of the
    region from the continuum."""
    if len(peaks) < 2:
        return []
    smooth = gaussian_filter1d(np.where(np.isfinite(flux[lo:hi + 1]), flux[lo:hi + 1], 1.0),
                               sigma=1.0)
    e = err[lo:hi + 1]
    e = e[np.isfinite(e) & (e > 0)]
    level = 1.0 - GAP_SIGMA * (float(np.median(e)) if e.size else 0.0)
    out = []
    for p, q in zip(peaks[:-1], peaks[1:], strict=True):
        k = p.index - lo + int(np.argmax(smooth[p.index - lo:q.index - lo + 1]))
        if smooth[k] >= level and p.index < lo + k < q.index:
            out.append(lo + k)
    return out


def _regions(spec: Spectrum, labels: np.ndarray, significance: np.ndarray,
             v: np.ndarray, flux: np.ndarray, err: np.ndarray, z_em: float,
             x_range: Optional[Tuple[float, float]]
             ) -> Tuple[List[_Region], np.ndarray]:
    """The regions of the identification and their labels (ids from 1): the
    runs of the matched filters, each split at its gaps (:func:`_gaps`)
    when :data:`SPLIT_GAPS`."""
    x_nm = np.exp(v / C_KMS)
    lya = (1.0 + z_em) * xem_d['Ly_a'].to_value(au.nm)
    lylim = (1.0 + z_em) * xem_d['Ly_lim'].to_value(au.nm)
    out: List[_Region] = []
    new_labels = np.zeros(len(labels), dtype=int)
    n = int(labels.max(initial=0))
    if n == 0:
        return out, new_labels
    starts = np.full(n + 1, -1)
    ends = np.full(n + 1, -1)
    idx = np.flatnonzero(labels)
    lab = labels[idx]
    # The regions of scipy.ndimage.label are contiguous runs
    first = np.r_[True, lab[1:] != lab[:-1]]
    last = np.r_[lab[1:] != lab[:-1], True]
    starts[lab[first]] = idx[first]
    ends[lab[last]] = idx[last]
    for rid in range(1, n + 1):
        lo, hi = int(starts[rid]), int(ends[rid])
        if lo < 0:
            continue
        centre = float(np.sqrt(x_nm[lo] * x_nm[hi]))
        if x_range is not None and not x_range[0] <= centre <= x_range[1]:
            continue
        peaks = _peaks(v, flux, err, lo, hi)
        cuts = _gaps(v, flux, err, lo, hi, peaks) if SPLIT_GAPS else []
        for s_lo, s_hi in zip([lo] + [c + 1 for c in cuts], cuts + [hi], strict=True):
            centre = float(np.sqrt(x_nm[s_lo] * x_nm[s_hi]))
            zone = "red" if centre >= lya else ("forest" if centre > lylim else "blue")
            sig = significance[s_lo:s_hi + 1]
            sig = sig[np.isfinite(sig)]
            seg = tuple(p for p in peaks if s_lo <= p.index <= s_hi)
            new_id = len(out) + 1
            new_labels[s_lo:s_hi + 1] = new_id
            out.append(_Region(id=new_id, lo=s_lo, hi=s_hi, zone=zone,
                               significance=float(sig.max()) if sig.size else float('nan'),
                               peaks=seg or _peaks(v, flux, err, s_lo, s_hi)))
    return out, new_labels


# --- Edges ---------------------------------------------------------------

def _r2(v: np.ndarray, norm: np.ndarray, sel: np.ndarray, lag: float, gamma: float) -> float:
    """The R^2 of 1 - F_s against 1 - F_p^gamma on the pixels ``sel`` of the
    strong line, the partner read at ``v + lag``; 0 when it cannot be rated."""
    if np.count_nonzero(sel) < 3:
        return 0.0
    y_data = 1.0 - norm[sel]
    y_model = 1.0 - np.power(np.interp(v[sel] + lag, v, norm, left=1.0, right=1.0), gamma)
    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
    if ss_tot == 0:
        return 0.0
    return float(max(0.0, 1.0 - np.sum((y_data - y_model) ** 2) / ss_tot))


def _lag_gamma(lines: Tuple[str, str]) -> Tuple[float, float]:
    strong, partner = lines
    lag = C_KMS * np.log(ATOM_DATA[partner]['wave'] / ATOM_DATA[strong]['wave'])
    return float(lag), _flam(strong) / _flam(partner)


def _norm(flux: np.ndarray) -> np.ndarray:
    return np.clip(np.where(np.isfinite(flux), flux, 1.0), 0.0, 1.5)


def pair_score(v: np.ndarray, flux: np.ndarray, region_strong: Tuple[int, int],
               region_partner: Tuple[int, int], lines: Tuple[str, str],
               z: float) -> float:
    """The R^2 of F_s = F_p^gamma of a pair on its two regions (the score of
    the current method,
    :func:`~astrocook.core.spectrum_operations.rate_doublet_candidate`, for
    any two transitions): the absorption of the strong line against that of
    the partner raised to gamma = (f lambda)_s / (f lambda)_p, on the
    velocities where the two regions overlap at the lag, 30 km/s more on
    each side. 0 when it cannot be rated. The score of PR 7c.4c, replaced
    in the identification by :func:`peak_score` (F169); kept as the
    reference of the current method."""
    lag, gamma = _lag_gamma(lines)
    if gamma <= 0 or not np.isfinite(gamma):
        return 0.0
    s_lo, s_hi = region_strong
    p_lo, p_hi = region_partner
    if s_hi - s_lo < 1 or p_hi - p_lo < 1:
        return 0.0
    start = max(v[s_lo], v[p_lo] - lag) - 30.0
    end = min(v[s_hi], v[p_hi] - lag) + 30.0
    if start >= end:
        return 0.0
    return _r2(v, _norm(flux), (v >= start) & (v <= end), lag, gamma)


def peak_score(c: _Candidate, by_id: Mapping[int, _Region], v: np.ndarray,
               flux: np.ndarray, fwhm: np.ndarray,
               matched: Optional[Sequence[Tuple[int, int]]] = None) -> float:
    """The R^2 of F_s = F_p^gamma of a pair on the peaks of the de-blending
    that it matches (F169, PR 7c.4c-bis): for each peak of the strong line
    whose partner falls on a peak of the partner's region
    (:func:`anchor_peaks`), the velocities of the share of both peaks
    (:func:`_range`, at least +-FWHM/2 of the instrument), the union of
    them. Not the whole regions (:func:`pair_score`): the regions of the
    matched filters are wide (F168), and a true doublet whose two lines
    fall in one region, or whose regions also hold other lines, scored
    <= 0 there. 0 when no peak matches or it cannot be rated."""
    if matched is None:
        matched = anchor_peaks(c, by_id, v, fwhm)
    if not matched:
        return 0.0
    lag, gamma = _lag_gamma(c.lines)
    if gamma <= 0 or not np.isfinite(gamma):
        return 0.0
    a, b = by_id[c.regions[0]], by_id[c.regions[1]]
    sel = np.zeros(len(v), dtype=bool)
    for k, kk in matched:
        vs = v[a.peaks[k].index]
        vp = v[b.peaks[kk].index] - lag
        rs, rp = _range(a, k, v, fwhm), _range(b, kk, v, fwhm)
        start, end = max(vs + rs[0], vp + rp[0]), min(vs + rs[1], vp + rp[1])
        if start < end:
            lo, hi = np.searchsorted(v, [start, end], side="left")
            sel[lo:hi] = True
    return _r2(v, _norm(flux), sel, lag, gamma)


def _candidates(regions: List[_Region], labels: np.ndarray, v: np.ndarray,
                flux: np.ndarray, fwhm: np.ndarray, anchors: Sequence[str], z_max: float
                ) -> List[_Candidate]:
    """Every pair of regions at the lag of an anchor pair (the overlap of
    the regions of the strong line with those of the partner moved back by
    the lag, as the kinematic search of the current method), with its
    R^2 on the peaks it matches (:func:`peak_score`)."""
    by_id = {r.id: r for r in regions}
    out: List[_Candidate] = []
    for species, (strong, partner) in anchor_pairs(anchors):
        zone = "forest" if species == LYAB else "red"
        lag = C_KMS * np.log(ATOM_DATA[partner]['wave'] / ATOM_DATA[strong]['wave'])
        at = v + lag
        inside = (at >= v[0]) & (at <= v[-1])
        j = np.clip(np.searchsorted(v, at), 1, len(v) - 1)
        j = np.where(np.abs(at - v[j - 1]) <= np.abs(v[j] - at), j - 1, j)
        shifted = np.where(inside, labels[j], 0)
        both = (labels > 0) & (shifted > 0)
        if not np.any(both):
            continue
        pairs = np.unique(np.stack((labels[both], shifted[both])), axis=1)
        for rs, rp in pairs.T:
            a, b = by_id.get(int(rs)), by_id.get(int(rp))
            # (a and b can be one region: the two lines of a strong or a
            # broad system in one run of the filters)
            if a is None or b is None or a.zone != zone or b.zone != zone:
                continue
            over = np.flatnonzero((labels == rs) & (shifted == rp))
            centre = float(np.exp(0.5 * (v[over[0]] + v[over[-1]]) / C_KMS))
            z = centre / _wave_nm(strong) - 1.0
            if z < 0 or z > z_max:
                continue
            cand = _Candidate(len(out), species, (strong, partner), (a.id, b.id), z, 0.0)
            score = peak_score(cand, by_id, v, flux, fwhm)
            out.append(dataclasses.replace(cand, score=score))
    return out


def _root(parent: Dict[int, int], r: int) -> int:
    while parent[r] != r:
        parent[r] = parent[parent[r]]
        r = parent[r]
    return r


def _components_of(regions: List[_Region], edges: List[_Candidate]
                   ) -> Tuple[Dict[int, int], Dict[int, List[int]]]:
    """The root of each region in the graph of ``edges``, and the regions of
    each root (union-find)."""
    parent: Dict[int, int] = {r.id: r.id for r in regions}
    for e in edges:
        a, b = _root(parent, e.regions[0]), _root(parent, e.regions[1])
        if a != b:
            parent[max(a, b)] = min(a, b)
    members: Dict[int, List[int]] = {}
    for r in regions:
        members.setdefault(_root(parent, r.id), []).append(r.id)
    return {r.id: _root(parent, r.id) for r in regions}, members


def _groups(regions: List[_Region], edges: List[_Candidate], max_group: int
            ) -> Tuple[List[Tuple[List[int], List[_Candidate]]], List[_Candidate]]:
    """The connected components of the graph of the regions and the edges,
    those above ``max_group`` regions split at their weakest edges (one at a
    time, the lowest R^2 first); the edges cut."""
    kept = sorted(edges, key=lambda e: -e.score)
    cut: List[_Candidate] = []
    while True:
        root, members = _components_of(regions, kept)
        big = {k for k, m in members.items() if len(m) > max_group}
        if not big:
            break
        for k in big:
            inside = [e for e in kept if root[e.regions[0]] == k]
            weakest = min(inside, key=lambda e: (e.score, -e.index))
            kept.remove(weakest)
            cut.append(weakest)
    by_root: Dict[int, List[_Candidate]] = {}
    for e in kept:
        by_root.setdefault(root[e.regions[0]], []).append(e)
    groups = [(sorted(m), sorted(by_root.get(k, []), key=lambda e: e.index))
              for k, m in members.items()]
    groups.sort(key=lambda g: g[0][0])
    return groups, cut


# --- Hypotheses ----------------------------------------------------------

def _range(region: _Region, k: int, v: np.ndarray, fwhm: np.ndarray) -> Tuple[float, float]:
    """The AOD range of peak ``k`` of a region, in km/s from the peak: its
    share of the region, from the midpoint with the previous peak (or the
    start of the region) to that with the next (or its end), at least
    +-FWHM/2 of the instrument."""
    p = region.peaks[k]
    vp = v[p.index]
    lo, hi = v[region.lo] - vp, v[region.hi] - vp
    for q in region.peaks:
        d = v[q.index] - vp
        if d < 0:
            lo = max(lo, 0.5 * d)
        elif d > 0:
            hi = min(hi, 0.5 * d)
    half = 0.5 * fwhm[p.index]
    return (min(lo, -half), max(hi, half))


def _match(region: _Region, v_at: float, v: np.ndarray, fwhm: np.ndarray,
           width_kms: float = 0.0) -> Optional[int]:
    """The peak of ``region`` nearest to the velocity ``v_at``, within half
    the larger of its width, ``width_kms`` (that of the peak projected) and
    the instrumental FWHM; None if none."""
    best, best_d = None, np.inf
    for k, q in enumerate(region.peaks):
        d = abs(v[q.index] - v_at)
        if d <= 0.5 * max(q.width_kms, width_kms, fwhm[q.index]) and d < best_d:
            best, best_d = k, d
    return best


def anchor_peaks(c: _Candidate, by_id: Mapping[int, _Region], v: np.ndarray,
                 fwhm: np.ndarray) -> List[Tuple[int, int]]:
    """The peaks of the region of the strong line of a candidate whose
    partner line falls on a peak of the region of the partner (pairs of
    their indices): where the pair is a doublet. None for a candidate whose
    regions overlap at the lag but whose peaks do not: its pair says
    nothing of the species."""
    a, b = by_id[c.regions[0]], by_id[c.regions[1]]
    lag = C_KMS * np.log(ATOM_DATA[c.lines[1]]['wave'] / ATOM_DATA[c.lines[0]]['wave'])
    out = []
    for k, p in enumerate(a.peaks):
        kk = _match(b, v[p.index] + lag, v, fwhm, p.width_kms)
        if kk is not None and not (a.id == b.id and kk == k):
            out.append((k, kk))
    return out


def _components(group_regions: Sequence[_Region], anchors: Sequence[_Candidate],
                v: np.ndarray, fwhm: np.ndarray, lya: bool) -> Tuple[List[_Comp], frozenset]:
    """The components of the hypothesis with ``anchors``, and the peaks the
    anchors explain: each anchor a component of its species at every peak
    of the region of its strong line whose partner falls on a peak of the
    region of the partner line (:func:`anchor_peaks`); every other peak of
    the group its default label (a peak of a metal line that is not that
    metal is left to its own line, not to a component that the partner
    line contradicts)."""
    by_id = {r.id: r for r in group_regions}
    explained = set()
    comps: List[_Comp] = []
    for c in anchors:
        a, b = by_id[c.regions[0]], by_id[c.regions[1]]
        for k, kk in anchor_peaks(c, by_id, v, fwhm):
            p = a.peaks[k]
            z = float(np.exp(v[p.index] / C_KMS)) / _wave_nm(c.lines[0]) - 1.0
            comps.append(_Comp(c.species, z, a.id, k, "anchor", c.index, _range(a, k, v, fwhm),
                               p.width_kms, float(fwhm[p.index])))
            explained.add((a.id, k))
            explained.add((b.id, kk))
    anchored = frozenset(explained)
    for r in group_regions:
        default = "lya" if (lya and r.zone == "forest") else "unidentified"
        for k, p in enumerate(r.peaks):
            if (r.id, k) in explained:
                continue
            z = float(np.exp(v[p.index] / C_KMS)) / _wave_nm(PLACEHOLDER) - 1.0
            comps.append(_Comp(PLACEHOLDER, z, r.id, k, default, None, _range(r, k, v, fwhm),
                               p.width_kms, float(fwhm[p.index])))
    return comps, anchored


def aod_guess(v: np.ndarray, flux: np.ndarray, err: np.ndarray, lines: Sequence[str],
              z: float, aod_range: Tuple[float, float], fwhm_kms: float
              ) -> Tuple[float, float]:
    """logN and b of a component from the apparent optical depth of its
    lines within ``aod_range`` (km/s from each line): N the largest over the
    lines (the weaker one corrects the saturation of the stronger), b from
    the FWHM of the apparent optical depth of that line, the instrumental
    one (``fwhm_kms``) subtracted in quadrature (a second moment of N(v) over
    the range of a blend took the noise of its wings; the width of the flux
    of a saturated line was too wide). Within :data:`AOD_LOGN` and
    :data:`AOD_B`."""
    n_best, b = 0.0, 10.0
    for t in lines:
        vc = C_KMS * np.log(_wave_nm(t) * (1.0 + z))
        if vc < v[0] or vc > v[-1]:
            continue
        sel = (v >= vc + aod_range[0]) & (v <= vc + aod_range[1]) & np.isfinite(flux)
        if np.count_nonzero(sel) < 2:
            continue
        floor = np.maximum(AOD_FLOOR, np.where(np.isfinite(err[sel]), err[sel], AOD_FLOOR))
        tau = -np.log(np.clip(flux[sel], floor, None))
        vel = v[sel] - vc
        n_t = float(np.trapezoid(AOD_CONST / _flam(t) * tau, vel))
        if n_t > n_best:
            n_best = n_t
            peak = float(np.max(tau))
            if peak > 0:
                dv = np.abs(np.gradient(vel)) if len(vel) > 1 else np.ones(1)
                width = float(np.sum(dv[tau >= 0.5 * peak]))
                b = np.sqrt(max(width ** 2 - fwhm_kms ** 2, 0.0)) / (2.0 * np.sqrt(np.log(2.0)))
    logn = float(np.clip(np.log10(n_best) if n_best > 0 else AOD_LOGN[0], *AOD_LOGN))
    return logn, float(np.clip(b, *AOD_B))


def _aod_lines(comp: _Comp, cands: Mapping[int, _Candidate]) -> Tuple[str, ...]:
    if comp.role == "anchor" and comp.cand is not None:
        c = cands[comp.cand]
        # The Ly-b falls in the Ly-b forest: its AOD would be that of the
        # blends (the experiment of the design, section 4)
        return (c.lines[0],) if c.species == LYAB else c.lines
    return (PLACEHOLDER,)


def _build(h: _Hyp, cands: Mapping[int, _Candidate], v, flux, err) -> None:
    """The system list of a hypothesis, with the AOD guesses."""
    systs = _empty()
    for c in h.comps:
        logn, b = aod_guess(v, flux, err, _aod_lines(c, cands), c.z, c.aod_range,
                            c.fwhm_kms)
        systs = systs.add_component(c.series, c.z, logN=logn, b=b)
    h.systs = systs
    h.uuids = [c.uuid for c in systs.components]


def _measure(win: EvaluationWindow, systs: SystemList) -> Dict[str, float]:
    ev = win.evaluate(systs, min_influence=0.0)
    out = {"chi2": ev.chi2, "aic": ev.aic, "aicc": ev.aicc, "bic": ev.bic,
           "spic": ev.spic, "spic_a": ev.spic_a, "spic_b": ev.spic_b,
           "k": ev.k, "n": ev.n_pixels, "chi2_reduced": ev.chi2_reduced}
    return {k: float(val) for k, val in out.items() if val is not None}


def _label(h: _Hyp, cands: Mapping[int, _Candidate]) -> str:
    parts = [f"{cands[a].species} at z = {cands[a].z:.5f}" for a in h.anchors]
    rest = [c for c in h.comps if c.role != "anchor"]
    if rest:
        n_lya = sum(1 for c in rest if c.role == "lya")
        n_un = len(rest) - n_lya
        if n_lya:
            parts.append(f"{n_lya} Ly_a (assumed)")
        if n_un:
            parts.append(f"{n_un} unidentified")
    return " + ".join(parts) if parts else "no line"


@dataclasses.dataclass
class GroupResult:
    """The outcome of a group (for the structure and the output)."""
    regions: List[int]
    edges: List[_Candidate]
    hyps: List[_Hyp]
    accepted: int


class _Context:
    """What the evaluation of the groups shares."""

    def __init__(self, spec: Spectrum, fixed: Optional[SystemList], criterion: str,
                 screen: Optional[float], temper: bool, lya: bool, labels: np.ndarray,
                 forest: Tuple[float, float], refine: bool = True):
        self.refine = refine
        self.labels = labels
        self.forest = forest
        """The forest, Lyman limit to Ly-a (velocity, c ln lambda[nm])."""
        self.spec = spec
        self.fixed = fixed
        self.criterion = criterion
        self.screen = screen
        self.temper = temper
        self.lya = lya
        self.x_range: Optional[Tuple[float, float]] = None
        """The range of the regions (nm), when limited: the sweep tests only
        the lines in it (the pixels outside have no regions, and their
        absorption no model)."""


def _window(ctx: _Context, systs: Sequence[SystemList], ids: Sequence[int]
            ) -> EvaluationWindow:
    """The window of a group (:meth:`EvaluationWindow.around`) without the
    pixels of the regions of other groups (PR 7c.4c-bis): their absorption
    is in no hypothesis of the group, and the hypotheses with more free
    lines took it (with the regions split at their gaps, F169, the window
    of a group of the forest reached its neighbours: chi2 of hundreds per
    pixel).

    Raises
    ------
    ValueError
        If no pixel is left.
    """
    win = EvaluationWindow.around(ctx.spec, systs, fixed=ctx.fixed)
    keep = win.mask & ((ctx.labels == 0) | np.isin(ctx.labels, list(ids)))
    if keep.sum() == win.index.size:
        return win
    return EvaluationWindow.from_mask(ctx.spec, keep, fixed=ctx.fixed)


def _evaluate_group(ctx: _Context, regions: List[_Region], edges: List[_Candidate],
                    v, flux, err, fwhm) -> GroupResult:
    """The hypotheses of a group, evaluated and weighted (module, steps 4-6)."""
    cands = {e.index: e for e in edges}
    ids = [r.id for r in regions]
    hyps = [_Hyp((e.index,), *_components(regions, [e], v, fwhm, ctx.lya)) for e in edges]
    if edges:
        # Nothing: no line where the noise made the regions; refined by its
        # residuals, the fewest free lines the data need
        hyps.append(_Hyp((), []))
    hyps.append(_Hyp((), *_components(regions, [], v, fwhm, ctx.lya)))
    if not edges:
        # One hypothesis: nothing to compare, its AOD parameters for the output
        h = hyps[0]
        _build(h, cands, v, flux, err)
        h.status = ident.Status.ACCEPTED
        h.reason = ("assumed: a region of the forest that no pair links" if any(
            c.role == "lya" for c in h.comps) else "no pair links the region")
        return GroupResult(ids, edges, hyps, 0)
    for h in hyps:
        _build(h, cands, v, flux, err)
    crit = ctx.criterion
    try:
        win = _window(ctx, [h.systs for h in hyps if h.systs is not None], ids)
        for h in hyps:
            assert h.systs is not None
            h.aod = _measure(win, h.systs)
    except ValueError as e:
        # No valid pixel around the group: the default labels, unevaluated
        logger.debug(f"A group of regions {ids} cannot be evaluated ({e}).")
        default = hyps[-1]
        default.status = ident.Status.ACCEPTED
        default.reason = f"not evaluated: {e}"
        for h in hyps[:-1]:
            h.status = ident.Status.DROPPED
            h.reason = f"not evaluated: {e}"
        return GroupResult(ids, edges, hyps, len(hyps) - 1)

    k_max = max(h.aod.get("k", 0.0) for h in hyps if h.aod is not None)
    if k_max > MAX_FIT_K:
        default = hyps[-1]
        default.status = ident.Status.ACCEPTED
        default.reason = (f"not compared: a hypothesis of k = {k_max:.0f} free parameters, "
                          f"above {MAX_FIT_K} for the short fit; the default labels")
        for h in hyps[:-1]:
            h.status = ident.Status.DROPPED
            h.reason = f"not compared: the group is too large for the short fit (k = {k_max:.0f})"
            h.weight = 0.0
        return GroupResult(ids, edges, hyps, len(hyps) - 1)

    best_aod = min(h.aod[crit] for h in hyps if h.aod is not None)
    for h in hyps:
        assert h.aod is not None
        d = h.aod[crit] - best_aod
        if ctx.screen is not None and d > ctx.screen:
            h.status = ident.Status.DROPPED
            h.reason = f"screen: {crit} {d:.1f} above the best by the AOD (threshold {ctx.screen:g})"
            h.weight = 0.0
            continue
        _fit(win, h)
        _refine(ctx, win, h, cands, v)

    finalists = [h for h in hyps if h.fit is not None]
    # The anchored finalists that explain different peaks: combined, from
    # the best one, each added when it lowers the criterion (design, section
    # 3.2; by the peaks, not the regions: a line of one system can be
    # blended with a line of another in one region, and then neither wins
    # alone against "all default", which fits the blend with free lines)
    default = hyps[-1]
    anchored = sorted((h for h in finalists if h.anchors),
                      key=lambda h: h.fit[crit])  # type: ignore[index]
    if len(anchored) > 1:
        current = anchored[0]
        for h in anchored[1:]:
            if current.explained & h.explained:
                continue
            anchors = list(current.anchors) + list(h.anchors)
            trial = _Hyp(tuple(anchors),
                         *_components(regions, [cands[a] for a in anchors], v, fwhm, ctx.lya))
            _build(trial, cands, v, flux, err)
            assert trial.systs is not None
            trial.aod = _measure(win, trial.systs)
            if trial.aod.get("k", 0.0) > MAX_FIT_K:
                continue
            _fit(win, trial)
            _refine(ctx, win, trial, cands, v)
            hyps.insert(len(hyps) - 1, trial)
            if trial.fit is not None and trial.fit[crit] < current.fit[crit]:  # type: ignore[index]
                current = trial
        finalists = [h for h in hyps if h.fit is not None]

    if not finalists:
        default.status = ident.Status.ACCEPTED
        default.reason = "no hypothesis could be fitted: the default labels"
        return GroupResult(ids, edges, hyps, hyps.index(default))
    values = np.array([h.fit[crit] for h in finalists])  # type: ignore[index]
    i_best = int(np.argmin(values))
    scale = 1.0
    if ctx.temper:
        chi2_nu = finalists[i_best].fit.get("chi2_reduced")  # type: ignore[union-attr]
        scale = max(1.0, chi2_nu) if chi2_nu is not None and np.isfinite(chi2_nu) else 1.0
    delta = (values - values[i_best]) / scale
    w = np.exp(-0.5 * delta)
    w /= w.sum()
    for h, wi, d in zip(finalists, w, delta, strict=True):
        h.weight = float(wi)
        h.status = ident.Status.ALTERNATIVE
        h.reason = f"{crit} {d:+.1f} from the accepted"
    accepted = finalists[i_best]
    accepted.status = ident.Status.ACCEPTED
    others = sorted((h for h in finalists if h is not accepted), key=lambda h: -(h.weight or 0))
    accepted.reason = (f"lowest {crit} (w = {accepted.weight:.2f}" +
                       (f"; next {_label(others[0], cands)}, w = {others[0].weight:.2f}"
                        if others else "") + ")")
    logger.debug(f"{_label(accepted, cands)} (w = {accepted.weight:.2f})"
                 + "".join(f" against {_label(h, cands)} ({h.weight:.2f})" for h in others)
                 + f", {crit}")
    return GroupResult(ids, edges, hyps, hyps.index(accepted))


def _added(ctx: _Context, h: _Hyp, cands: Mapping[int, _Candidate], i: int,
           v: np.ndarray) -> List[_Comp]:
    """The components that can be added to ``h`` at pixel ``i`` by the
    refinement (:func:`_refine`): its default label, and, when the pixel is in a
    region of one of its anchors, the species of the anchor on the line of
    that region (a region of a metal line can also hold a blend, or a
    depression of the continuum)."""
    rid = int(ctx.labels[i])
    lam = float(np.exp(v[i] / C_KMS))
    in_forest = ctx.forest[0] < v[i] < ctx.forest[1]
    role = "lya" if (ctx.lya and in_forest) else "unidentified"
    out = [_Comp(PLACEHOLDER, lam / _wave_nm(PLACEHOLDER) - 1.0, rid, -1, role)]
    for a in h.anchors:
        c = cands[a]
        for r, t in zip(c.regions, c.lines, strict=True):
            if rid == r:
                out.append(_Comp(c.species, lam / _wave_nm(t) - 1.0, rid, -1, "anchor", a))
                return out
    return out


def _refine(ctx: _Context, win: EvaluationWindow, h: _Hyp,
            cands: Mapping[int, _Candidate], v: np.ndarray, extra: int = 0) -> None:
    """Components added to a fitted hypothesis where its residuals are most
    negative, each the best of :func:`_added`, while it lowers the criterion
    (:data:`MAX_ADDED`, ``extra`` more; :data:`MAX_FIT_K`); nothing when
    ``ctx.refine`` is off."""
    crit = ctx.criterion
    if not ctx.refine:
        return
    added = 0
    while MAX_ADDED is None or added < MAX_ADDED + extra:
        added += 1
        if h.fitted is None or h.fit is None:
            return
        res = win.evaluate(h.fitted, min_influence=0.0).residuals[win.index]
        res = np.where(np.isfinite(res), res, 0.0)
        s3 = (res + np.r_[res[1:], 0.0] + np.r_[0.0, res[:-1]]) / np.sqrt(3.0)
        k = int(np.argmin(s3))
        if s3[k] > -REFINE_SIGMA or h.fit.get("k", 0.0) + 3 > MAX_FIT_K:
            return
        best = None
        for comp in _added(ctx, h, cands, int(win.index[k]), v):
            trial = h.fitted.add_component(comp.series, comp.z, logN=ADDED_LOGN, b=ADDED_B)
            try:
                fitted = win.fit(trial, max_nfev=FIT_NFEV, b_max=FIT_B_MAX,
                                 logN_max=FIT_LOGN_MAX)
                measured = _measure(win, fitted)
            except (NoPixelsToFit, ValueError, np.linalg.LinAlgError):
                continue
            if best is None or measured[crit] < best[2][crit]:
                best = (comp, fitted, measured, trial.components[-1].uuid)
        if best is None or best[2][crit] >= h.fit[crit]:
            return
        comp, h.fitted, h.fit, uuid = best
        h.comps.append(comp)
        h.uuids.append(uuid)


def _fit(win: EvaluationWindow, h: _Hyp) -> None:
    """The short fit of a finalist and its evaluation after it."""
    assert h.systs is not None
    try:
        h.fitted = win.fit(h.systs, max_nfev=FIT_NFEV, b_max=FIT_B_MAX,
                           logN_max=FIT_LOGN_MAX)
        h.fit = _measure(win, h.fitted)
    except (NoPixelsToFit, ValueError, np.linalg.LinAlgError) as e:
        # A fit that fails leaves the hypothesis out of the comparison
        h.status = ident.Status.DROPPED
        h.reason = f"the short fit failed: {e}"
        h.weight = 0.0


# --- The whole spectrum ----------------------------------------------------

@dataclasses.dataclass(frozen=True)
class HypothesesResult:
    """The spectrum with the identification, and its summary."""
    spectrum: Spectrum
    seconds: float
    n_regions: int
    n_groups: int
    n_compared: int
    n_hypotheses: int


def identify_hypotheses(spec: Spectrum, multiplets: Sequence[str],
                        fixed: Optional[SystemList] = None,
                        criterion: str = "spic", prefilter: float = PREFILTER_R2,
                        screen: Optional[float] = SCREEN_DELTA, max_group: int = MAX_GROUP,
                        temper: bool = False, refine: bool = True, sweep: str = "HI",
                        z_em_margin_kms: float = 3000.0,
                        x_range: Optional[Tuple[float, float]] = None,
                        params: Optional[Mapping[str, Any]] = None) -> HypothesesResult:
    """The identification as a comparison between hypotheses (see the
    module) of the regions of ``spec``.

    Parameters
    ----------
    spec : Spectrum
        With a continuum, the resolution and ``z_em``.
    multiplets : sequence of str
        The species searched (:func:`parse_species`).
    fixed : SystemList, optional
        Components that stay fixed in every evaluation (those already in the
        session).
    criterion : str
        ``"spic"`` (SpIC_H), ``"bic"`` or ``"aic"``: the weights.
    prefilter : float
        The R^2 of a pair above which it is an edge (:data:`PREFILTER_R2`).
    screen : float or None
        The screen of the AOD stage (:data:`SCREEN_DELTA`); None for none
        (every hypothesis fitted).
    max_group : int
        The largest group, in regions (:data:`MAX_GROUP`).
    temper : bool
        The differences of the criterion divided by max(1, chi2_nu) of the
        accepted model of the group (design, section 3.3).
    refine : bool
        The refinement of the finalists by their residuals (:func:`_refine`,
        F170).
    sweep : str
        The sweep of the anchored systems (PR 7c.4d,
        :mod:`astrocook.fitting.hypotheses_sweep`) at the redshift of each
        accepted metal system: ``"HI"`` (the H I, the default), ``"all"``
        (the H I and the other species of ``SWEEP_ALL``) or ``"none"``
        (:data:`~astrocook.fitting.hypotheses_sweep.SWEEP_CHOICES`).
    z_em_margin_kms : float
        No system beyond ``z_em`` by more than this (F55).
    x_range : (float, float), optional
        Only the regions whose centre is in this range of wavelength (nm).
    params : mapping, optional
        The parameters of the recipe, for the structure.

    Raises
    ------
    ValueError
        Without ``z_em``, a continuum or the resolution, for an unknown
        criterion, species or choice of the sweep.
    """
    t0 = time.perf_counter()
    if criterion not in CRITERIA:
        raise ValueError(f"Unknown criterion '{criterion}': one of {', '.join(CRITERIA)}.")
    from astrocook.fitting.hypotheses_sweep import SWEEP_CHOICES
    if sweep not in SWEEP_CHOICES:
        raise ValueError(f"Unknown sweep '{sweep}': one of {', '.join(SWEEP_CHOICES)}.")
    sweep_species = SWEEP_CHOICES[sweep]
    z_em = spec._data.z_em
    if not z_em:
        raise ValueError("z_em must be set before running identify_lines.")
    anchors, lya = parse_species(multiplets)
    z_max = z_em + z_em_margin_kms / C_KMS * (1.0 + z_em)
    metals = [s for s in anchors if s != LYAB and len(STANDARD_MULTIPLETS.get(s, ())) == 2]
    det = detect_lines(spec, multiplets=metals, z_max=z_max)

    from astrocook.core.spectrum_lines import _normalized_flux, _wavelength_nm
    x_nm = _wavelength_nm(spec)
    v = C_KMS * np.log(x_nm)
    flux = _normalized_flux(spec, "the spectrum")
    norm = np.asarray(au.Quantity(spec.norm).to_value(spec.y.unit), dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        err = np.asarray(spec.dy.to_value(spec.y.unit), dtype=float) / norm
    fwhm = C_KMS / spec.resolution()

    regions, labels = _regions(spec, det.regions, det.max_significance, v, flux, err, z_em,
                               x_range)
    all_cands = _candidates(regions, labels, v, flux, fwhm, anchors, z_max)
    by_id = {r.id: r for r in regions}
    unmatched = [c for c in all_cands if not anchor_peaks(c, by_id, v, fwhm)]
    edges = [c for c in all_cands if c.score >= prefilter and c not in unmatched]
    low = [c for c in all_cands if c.score < prefilter and c not in unmatched]
    groups, cut = _groups(regions, edges, max_group)

    forest = (C_KMS * np.log((1.0 + z_em) * xem_d['Ly_lim'].to_value(au.nm)),
              C_KMS * np.log((1.0 + z_em) * xem_d['Ly_a'].to_value(au.nm)))
    ctx = _Context(spec, fixed if fixed is not None and fixed.components else None,
                   criterion, screen, temper, lya, labels, forest, refine)
    ctx.x_range = x_range
    results = [_evaluate_group(ctx, [by_id[i] for i in g_ids], g_edges, v, flux, err, fwhm)
               for g_ids, g_edges in groups]
    swept = None
    if sweep_species:
        from astrocook.fitting.hypotheses_sweep import sweep as run_sweep
        swept = run_sweep(ctx, results, {c.index: c for c in all_cands}, v, flux, err, fwhm,
                          species=sweep_species)
    seconds = time.perf_counter() - t0
    run_params = dict(params or {})
    run_params.update(multiplets=list(multiplets), criterion=criterion, prefilter=prefilter,
                      screen=screen, max_group=max_group, temper=temper, refine=refine,
                      sweep=sweep, max_added=MAX_ADDED, split_gaps=SPLIT_GAPS,
                      z_em_margin_kms=z_em_margin_kms, seconds=round(seconds, 3))
    if swept is not None:
        run_params["sweep_seconds"] = round(swept.seconds, 3)
    new_spec = _output(spec, x_nm, regions, labels, all_cands, results, low, cut, unmatched,
                       prefilter, run_params, swept)
    n_compared = sum(1 for g in results if g.edges)
    n_hyps = sum(len(g.hyps) for g in results if g.edges)
    n_acc = sum(1 for g in results if g.hyps[g.accepted].anchors)
    logger.info(f"{len(regions)} regions, {len(edges)} pairs above the pre-filter (R^2 > "
                f"{prefilter:g}) in {n_compared} groups compared ({n_hyps} hypotheses): "
                f"{n_acc} identified by their pairs"
                + (f", {sum(t.accepted for t in swept.tests)} species added by the sweep"
                   if swept is not None else "") + f", in {seconds:.1f} s.")
    return HypothesesResult(new_spec, seconds, len(regions), len(results), n_compared, n_hyps)


def _output(spec: Spectrum, x_nm: np.ndarray, regions: List[_Region], labels: np.ndarray,
            cands: List[_Candidate], results: List[GroupResult], low: List[_Candidate],
            cut: List[_Candidate], unmatched: List[_Candidate], prefilter: float,
            params: Mapping[str, Any], swept: Optional["SweepResult"] = None) -> Spectrum:
    """The spectrum with the ``abs_ids`` column, the maps of the candidate
    systems, the labels of the regions, the dropped candidates and the
    structure (the form of the current method)."""
    import json
    from copy import deepcopy

    from astrocook.core.structures import DataColumn

    by_cand = {c.index: c for c in cands}
    by_id = {r.id: r for r in regions}
    replaced = swept.replaced if swept is not None else set()
    # The values of the maps: those after the sweep (its accepted tests fit
    # the default components around their lines again)
    values = swept.values if swept is not None else {}
    series_map: Dict[int, str] = {}
    z_map: Dict[int, float] = {}
    logn_map: Dict[int, float] = {}
    b_map: Dict[int, float] = {}
    region_map: Dict[int, int] = {}
    score_map: Dict[int, float] = {}
    region_info: Dict[int, List[List[Any]]] = {}
    kinds: Dict[int, ident.RegionKind] = {}
    transitions: Dict[int, List[str]] = {}
    s_groups = []
    for gi, g in enumerate(results):
        s_hyps = []
        for hi, h in enumerate(g.hyps):
            # The components in the order of h.comps (by their uuids)
            last = h.fitted if h.fitted is not None else h.systs
            final = {c.uuid: c for c in last.components} if last is not None else {}
            comps = [final[u] for u in h.uuids]
            is_acc = hi == g.accepted
            s_comps = []
            for c, sc in zip(h.comps, comps, strict=True):
                sid = None
                if is_acc and c.role != "unidentified" and sc.uuid not in replaced:
                    sid = len(series_map) + 1
                    now = values.get(sc.uuid, sc)
                    series_map[sid] = c.series
                    z_map[sid] = float(now.z)
                    logn_map[sid] = float(now.logN)
                    b_map[sid] = float(now.b)
                    region_map[sid] = c.region
                    score_map[sid] = float(h.weight) if c.role == "anchor" and h.weight else 0.0
                s_comps.append(ident.Component(
                    series=c.series, z=float(sc.z), logN=float(sc.logN), b=float(sc.b),
                    dz=_finite(sc.dz), dlogN=_finite(sc.dlogN), db=_finite(sc.db),
                    candidate=sid))
            stage = (ident.Stage.FIT if h.fit is not None else ident.Stage.AOD)
            stats = h.fit if h.fit is not None else (h.aod or {})
            criteria = {k: val for k, val in stats.items()
                        if k not in ("k", "n", "chi2", "chi2_reduced")}
            if h.anchors:
                criteria["score"] = float(max(by_cand[a].score for a in h.anchors))
            if h.aod is not None and h.fit is not None:
                criteria.update({f"aod_{k}": val for k, val in h.aod.items()
                                 if k in ("chi2", "spic", "bic", "aic")})
            regs = tuple(sorted({c.region for c in h.comps}
                                | {r for a in h.anchors for r in by_cand[a].regions}))
            edge_pos = {e.index: i for i, e in enumerate(g.edges)}
            s_hyps.append(ident.Hypothesis(
                anchors=tuple(edge_pos[a] for a in h.anchors),
                label=_label(h, by_cand), regions=regs, components=tuple(s_comps),
                stage=stage, status=ident.Status(h.status), reason=h.reason, weight=h.weight,
                chi2=stats.get("chi2"), k=_int(stats.get("k")), n=_int(stats.get("n")),
                criteria=criteria))
            if is_acc:
                for a in h.anchors:
                    cand = by_cand[a]
                    w = float(h.weight) if h.weight is not None else 0.0
                    kind = (ident.RegionKind.FOREST if cand.species == LYAB
                            else ident.RegionKind.METAL)
                    for rid, t in zip(cand.regions, cand.lines, strict=True):
                        kinds[rid] = kind
                        transitions.setdefault(rid, []).append(t)
                        region_info.setdefault(rid, []).append([t, w])
                for c, u in zip(h.comps, h.uuids, strict=True):
                    if c.role == "lya" and c.region not in kinds and u not in replaced:
                        kinds[c.region] = ident.RegionKind.FOREST_ASSUMED
                        transitions[c.region] = ["Ly_a"]
                        region_info[c.region] = [["Ly_a", 0.0]]
        edges = tuple(ident.Edge(regions=e.regions, transitions=e.lines, series=e.species,
                                 z=float(e.z), lag_kms=_lag(e.lines), prefilter=float(e.score))
                      for e in g.edges)
        s_groups.append(ident.Group(id=gi + 1, regions=tuple(g.regions), edges=edges,
                                    hypotheses=tuple(s_hyps)))
    if swept is not None:
        s_groups += _sweep_output(swept, by_cand, len(s_groups), series_map, z_map, logn_map,
                                  b_map, region_map, score_map, region_info, kinds, transitions)
    s_regions = tuple(ident.Region(
        id=r.id, x_min=float(x_nm[r.lo]), x_max=float(x_nm[r.hi]),
        kind=kinds.get(r.id, ident.RegionKind.UNIDENTIFIED),
        transitions=tuple(dict.fromkeys(transitions.get(r.id, []))),
        significance=_finite(r.significance)) for r in regions)
    s_dropped = tuple(ident.DroppedCandidate(
        series=c.species, z=float(c.z), score=float(c.score), threshold=float(prefilter),
        x_min=float(x_nm[by_id[c.regions[0]].lo]), x_max=float(x_nm[by_id[c.regions[0]].hi]),
        reason=why) for cs, why in ((low, "below the pre-filter"),
                                     (unmatched, "no peak of its lines at the lag"),
                                     (cut, "cut to split a group too large")) for c in cs)
    # The dropped candidates also as the current method writes them (F119)
    dropped = [dict(series=d.series, z=d.z, score=d.score, threshold=d.threshold,
                    x_nm=[d.x_min, d.x_max], reason=d.reason) for d in s_dropped]
    structure = ident.Identification(
        method="hypotheses", params=dict(params), regions=s_regions,
        groups=tuple(s_groups), dropped=s_dropped, x_unit="nm")

    meta = dict(spec.meta)
    meta['series_map_json'] = json.dumps(series_map)
    meta['z_map_json'] = json.dumps(z_map)
    meta['logN_map_json'] = json.dumps(logn_map)
    meta['b_map_json'] = json.dumps(b_map)
    meta['region_map_json'] = json.dumps(region_map)
    meta['score_map_json'] = json.dumps(score_map)
    meta['rejected_candidates_json'] = json.dumps(dropped)
    meta['region_identifications'] = json.dumps(region_info)
    meta[ident.META_KEY] = structure
    aux = deepcopy(spec._data.aux_cols)
    aux['abs_ids'] = DataColumn(values=labels.astype(int), unit=au.dimensionless_unscaled,
                                description="Absorption regions with identifications")
    data = dataclasses.replace(spec._data, meta=meta, aux_cols=aux)
    n_sys = len(series_map)
    return spec.__class__(data=data, history=spec.history + [
        f"Identified {n_sys} systems (hypotheses, {len(regions)} regions)"])


def _sweep_output(swept: "SweepResult", by_cand: Mapping[int, _Candidate], n_groups: int,
                  series_map: Dict[int, str], z_map: Dict[int, float],
                  logn_map: Dict[int, float], b_map: Dict[int, float],
                  region_map: Dict[int, int], score_map: Dict[int, float],
                  region_info: Dict[int, List[List[Any]]],
                  kinds: Dict[int, ident.RegionKind], transitions: Dict[int, List[str]]
                  ) -> List[ident.Group]:
    """The tests of the sweep in the maps of the candidate systems (the
    components of an accepted species, with the weight of its test as their
    score, and those its refinement added), the regions it explains
    relabelled, and one group of the structure per test."""
    from astrocook.fitting.hypotheses_sweep import structure
    out: List[ident.Group] = []
    for t in swept.tests:
        h = t.hyp
        last = h.fitted if h.fitted is not None else h.systs
        final = {c.uuid: c for c in last.components} if last is not None else {}
        s_comps = []
        for c, u in zip(h.comps, h.uuids, strict=True):
            if u not in final:
                continue
            sc = final[u]
            sid = None
            if t.accepted and c.role != "unidentified":
                sid = len(series_map) + 1
                series_map[sid] = c.series
                z_map[sid] = float(sc.z)
                logn_map[sid] = float(sc.logN)
                b_map[sid] = float(sc.b)
                region_map[sid] = c.region if c.region > 0 else by_cand[t.anchor].regions[0]
                score_map[sid] = float(t.weight or 0.0) if c.role == "sweep" else 0.0
            s_comps.append(ident.Component(
                series=c.series, z=float(sc.z), logN=float(sc.logN), b=float(sc.b),
                dz=_finite(sc.dz), dlogN=_finite(sc.dlogN), db=_finite(sc.db), candidate=sid))
        if t.accepted:
            kind = ident.RegionKind.FOREST if t.species == "HI" else ident.RegionKind.METAL
            for rid, line in t.lines:
                if kinds.get(rid) in (None, ident.RegionKind.FOREST_ASSUMED,
                                      ident.RegionKind.UNIDENTIFIED):
                    kinds[rid] = kind
                    transitions[rid] = []
                    region_info[rid] = []
                transitions.setdefault(rid, []).append(line)
                region_info.setdefault(rid, []).append([line, float(t.weight or 0.0)])
        out.append(structure(t, by_cand, n_groups + len(out) + 1, t.group + 1, s_comps))
    return out


def _lag(lines: Tuple[str, str]) -> float:
    return float(C_KMS * np.log(ATOM_DATA[lines[1]]['wave'] / ATOM_DATA[lines[0]]['wave']))


def _finite(x) -> Optional[float]:
    return float(x) if x is not None and np.isfinite(x) else None


def _int(x) -> Optional[int]:
    return int(x) if x is not None else None
