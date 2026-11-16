"""Truth map of a mock: what absorbs each pixel (phase 7, PR 7a.1).

From the truth alone, for the metrics of the benchmark (PR 7a.2). Each
category of absorber gets its own absorption depth per pixel, 1 - T, with
T its transmission alone, convolved with the instrumental profile and
averaged over the pixel as the flux of :func:`~astrocook.mocks.generate`.
The components are classified by a rule on their series and column:

- ``strong_hi``: H I only (Lyman lines), logN >= 17.2 (the lower bound of
  :data:`~astrocook.mocks.LLS_LOGN`), with its Lyman limit;
- ``forest``: H I only, below that, the H I of the metal systems
  included;
- ``metal``: any other series, the associated absorbers included;
- ``bal``: the BAL troughs; ``telluric``: the tellurics.
"""
import dataclasses
from typing import Dict

import numpy as np

from astrocook.mocks.generate import (
    bal_tau, components_optical_depth, pixel_transmission, subpixel_grid,
)
from astrocook.mocks.physics import transitions
from astrocook.mocks.populations import LLS_LOGN
from astrocook.mocks.telluric import telluric_optical_depth
from astrocook.mocks.truth import MockComponent, MockTruth

#: The categories, in the order of their bits in :meth:`TruthMap.mask`
CATEGORIES = ("forest", "metal", "strong_hi", "bal", "telluric")
#: Bit of each category in :meth:`TruthMap.mask`
BITS = {name: 1 << i for i, name in enumerate(CATEGORIES)}


def component_category(comp: MockComponent) -> str:
    """``'forest'``, ``'strong_hi'`` or ``'metal'`` (see the module)."""
    if all(t.startswith("Ly_") for t in transitions(comp.series)):
        return "strong_hi" if comp.logN >= LLS_LOGN[0] else "forest"
    return "metal"


@dataclasses.dataclass(frozen=True)
class TruthMap:
    """Absorption depth per pixel of each category.

    Attributes
    ----------
    x_nm : np.ndarray
        Pixel centres, as in the mock spectrum.
    depth : dict
        Category -> 1 - T per pixel (zeros where the category is absent).
    """
    x_nm: np.ndarray
    depth: Dict[str, np.ndarray]

    def mask(self, threshold: float = 0.01) -> np.ndarray:
        """Bits (:data:`BITS`) of the categories that absorb each pixel by
        more than ``threshold`` (the threshold of the metrics is a choice
        of PR 7a.2)."""
        out = np.zeros(len(self.x_nm), dtype=np.int16)
        for name, depth in self.depth.items():
            out[depth > threshold] |= BITS[name]
        return out

    def absorbed(self, category: str, threshold: float = 0.01) -> np.ndarray:
        """True where ``category`` absorbs by more than ``threshold``."""
        return self.depth[category] > threshold


def truth_map(truth: MockTruth) -> TruthMap:
    """The :class:`TruthMap` of ``truth``."""
    x, sub = subpixel_grid(truth)
    taus = {name: np.zeros_like(sub) for name in CATEGORIES}
    groups: Dict[str, list] = {"forest": [], "metal": [], "strong_hi": []}
    for comp in truth.components:
        groups[component_category(comp)].append(comp)
    for name, comps in groups.items():
        if comps:
            taus[name] = components_optical_depth(truth, sub, comps)
    if truth.bal_troughs:
        taus["bal"] = bal_tau(truth, sub)
    if truth.telluric is not None:
        taus["telluric"] = telluric_optical_depth(sub, truth.telluric)
    depth = {name: (1.0 - pixel_transmission(truth, np.exp(-tau), len(x))
                    if tau.any() else np.zeros(len(x)))
             for name, tau in taus.items()}
    return TruthMap(x, depth)
