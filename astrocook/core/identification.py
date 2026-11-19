"""The identification of the absorbed regions as a structure (phase 7, PR 7c.3).

The requirements are in ``design/identification.md``, section 9: per group
of regions, the regions and the edges that made it, the hypotheses that
explain it (their components, stage, statistics, weight, status and
reason), the candidates dropped by the score, and the parameters of the
run. The structure is immutable (frozen dataclasses, tuples) and lives in
the metadata of the spectrum under :data:`META_KEY`, next to the ``abs_ids``
column that it describes, so that it follows the history of the session
and the undo with no extra code. In the archive it is a file of its own,
``<name>_identification.json`` (``astrocook.io.archive_data``): a FITS
header cannot hold it.

:class:`RegionKind` is the vocabulary of the regions shared with the
emission model of 7b (``design/identification.md``, "Common decisions with
7b", item 1): 7b writes the regions it owns (BAL troughs, DLA cores and
wings, tellurics), the identification labels the others.

The current method of :meth:`~astrocook.core.spectrum.Spectrum.identify_lines`
(``method="doublets"``) fills the structure in a minimal form: every region
it finds, the accepted pairs as edges, one hypothesis per accepted pair
(stage :attr:`Stage.SCORE`, the R^2 of the pair as its only statistic) and
one per region of the forest assumed to be Ly-a, the dropped candidates.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any, Dict, Iterable, Iterator, Mapping, Optional, Tuple

#: The key of the structure in the metadata of the spectrum
META_KEY = "identification"
#: The version of the schema of :meth:`Identification.to_dict`: raised when
#: a field changes meaning or is removed (a new optional field does not
#: need it, the reader fills its default)
SCHEMA_VERSION = 1
#: The name of the schema, written in the file
SCHEMA_NAME = "astrocook.identification"


class RegionKind(StrEnum):
    """What absorbs in a region: the vocabulary shared with 7b."""
    BAL = "bal"
    """A trough of a broad absorption line (written by 7b)."""
    DLA = "dla"
    """The core or the wings of a damped Ly-a system (written by 7b)."""
    TELLURIC = "telluric"
    """Telluric absorption (from the ``telluric`` column of 7b.6)."""
    FOREST = "forest"
    """Ly-a of the forest, identified (by its Ly-b, or by a comparison)."""
    FOREST_ASSUMED = "forest_assumed"
    """Ly-a of the forest by assumption: no identification supports it."""
    METAL = "metal"
    """A metal line."""
    UNIDENTIFIED = "unidentified"
    """No identification."""


#: The category of the truth map of the mocks (``astrocook.mocks.truthmap``,
#: PR 7a.1) of each kind of region, for the metrics; ``None`` where the truth
#: map has no category (what is unidentified is compared with the truth)
TRUTH_CATEGORY: Dict[RegionKind, Optional[str]] = {
    RegionKind.BAL: "bal",
    RegionKind.DLA: "strong_hi",
    RegionKind.TELLURIC: "telluric",
    RegionKind.FOREST: "forest",
    RegionKind.FOREST_ASSUMED: "forest",
    RegionKind.METAL: "metal",
    RegionKind.UNIDENTIFIED: None,
}


class Stage(StrEnum):
    """How far a hypothesis was evaluated."""
    SCORE = "score"
    """Only the score of the pre-filter (the R^2 of a pair): the current method."""
    AOD = "aod"
    """The model evaluated with N and b from the apparent optical depth."""
    FIT = "fit"
    """The model evaluated after a short fit."""


class Status(StrEnum):
    """What became of a hypothesis."""
    ACCEPTED = "accepted"
    """Its components are those of the system list."""
    ALTERNATIVE = "alternative"
    """Evaluated and kept, with its weight, not accepted."""
    DROPPED = "dropped"
    """Dropped by the screen before the end of the evaluation."""


@dataclass(frozen=True)
class Region:
    """An absorbed region. ``id`` is its value in the ``abs_ids`` column
    when the region is there (a region that no identification labels is
    not: its ``id`` is unique in the structure only)."""
    id: int
    x_min: float
    x_max: float
    kind: RegionKind
    transitions: Tuple[str, ...] = ()
    """The transitions that label it (several for a region shared by
    candidates)."""
    significance: Optional[float] = None
    """The significance of the detection (sigma); ``None`` for the current
    method, which has none."""


@dataclass(frozen=True)
class Edge:
    """Two regions that a transition pair puts at the same redshift."""
    regions: Tuple[int, int]
    transitions: Tuple[str, str]
    """The transition of each region, in the order of ``regions``."""
    series: str
    """The multiplet of the pair (``"CIV"``, ``"Ly_ab"``, ...)."""
    z: float
    lag_kms: Optional[float] = None
    """The separation of the two transitions in velocity."""
    prefilter: Optional[float] = None
    """The value of the pre-filter (the R^2 of F1 = F2^gamma)."""


@dataclass(frozen=True)
class Component:
    """A component of a hypothesis."""
    series: str
    z: float
    logN: Optional[float] = None
    b: Optional[float] = None
    dz: Optional[float] = None
    dlogN: Optional[float] = None
    db: Optional[float] = None
    candidate: Optional[int] = None
    """The key of the candidate system in the maps of the spectrum
    (``series_map_json``, ...) that the current method writes."""
    uuid: Optional[str] = None
    """The component of the system list made from it (accepted hypotheses,
    after :meth:`~astrocook.recipes.absorbers.RecipeAbsorbers.populate_from_identification`)."""


@dataclass(frozen=True)
class Hypothesis:
    """One explanation of the regions of a group."""
    label: str
    """In words, e.g. ``"CIV at z = 2.22950"``."""
    regions: Tuple[int, ...]
    components: Tuple[Component, ...]
    stage: Stage
    status: Status
    reason: str = ""
    """Why it has its status, in words (for the log)."""
    weight: Optional[float] = None
    """Its credibility, exp(-dIC/2) normalised over the group."""
    chi2: Optional[float] = None
    k: Optional[int] = None
    """The number of free parameters."""
    n: Optional[int] = None
    """The number of pixels."""
    criteria: Mapping[str, float] = field(default_factory=dict)
    """Every information criterion computed (``"aic"``, ``"aicc"``,
    ``"bic"``, ``"spic"`` (SpIC_H), ...), and the ``"score"`` of the pre-filter."""
    anchors: Tuple[int, ...] = ()
    """The edges of the group (their index in :attr:`Group.edges`) that
    anchor it: none for the hypothesis of the default labels (PR 7c.4c)."""


@dataclass(frozen=True)
class Group:
    """Regions linked by the transitions that can explain them together."""
    id: int
    regions: Tuple[int, ...]
    edges: Tuple[Edge, ...]
    hypotheses: Tuple[Hypothesis, ...]


@dataclass(frozen=True)
class DroppedCandidate:
    """A candidate pair scored below its threshold (F119)."""
    series: str
    z: float
    score: float
    threshold: float
    x_min: float
    x_max: float
    reason: str = ""
    """Why it was dropped, when not by its score alone (PR 7c.4c: below
    the pre-filter, or cut to split a group too large)."""


@dataclass(frozen=True)
class Identification:
    """The identification of the absorbed regions of a spectrum.

    The positions (``x_min``, ``x_max`` of the regions and of the dropped
    candidates) are in ``x_unit``.
    """
    method: str
    params: Mapping[str, Any]
    """The parameters of the run, to replay it."""
    regions: Tuple[Region, ...] = ()
    groups: Tuple[Group, ...] = ()
    dropped: Tuple[DroppedCandidate, ...] = ()
    x_unit: str = "nm"

    def region(self, region_id: int) -> Optional[Region]:
        """The region ``region_id``, or ``None``."""
        return next((r for r in self.regions if r.id == region_id), None)

    def hypotheses(self, status: Optional[Status] = None) -> Iterator[Hypothesis]:
        """The hypotheses of every group (with ``status``, if given)."""
        for g in self.groups:
            for h in g.hypotheses:
                if status is None or h.status == status:
                    yield h

    def with_uuids(self, uuids: Mapping[int, str]) -> "Identification":
        """The structure with the uuids of the components made from the
        candidate systems (``{candidate: uuid}``)."""
        def comp(c: Component) -> Component:
            u = uuids.get(c.candidate) if c.candidate is not None else None
            return dataclasses.replace(c, uuid=u) if u is not None else c
        groups = tuple(dataclasses.replace(g, hypotheses=tuple(
            dataclasses.replace(h, components=tuple(comp(c) for c in h.components))
            for h in g.hypotheses)) for g in self.groups)
        return dataclasses.replace(self, groups=groups)

    def merged(self, new: "Identification", offset: int,
               old_column: Iterable[int], new_column: Iterable[int]) -> "Identification":
        """This identification and ``new`` of a later step of a pipeline, as
        :func:`~astrocook.recipes.absorbers.merge_region_ids` merges the
        ``abs_ids`` columns: the regions of ``new`` shifted by ``offset``.

        ``old_column`` and ``new_column`` are the region IDs of the two
        columns before the merge. The regions in neither (those that no
        identification labels, or that an earlier merge left out of the
        column) have IDs of the structure only: they are numbered again past
        the others, so that no two regions share one. The parameters are
        those of ``new`` (each run is in the log).
        """
        old_fixed = {int(i) for i in old_column if i > 0}
        new_fixed = {int(i) + offset for i in new_column if i > 0}
        new = new._remapped({r.id: r.id + offset for r in new.regions})
        old_fixed &= {r.id for r in self.regions}
        new_fixed &= {r.id for r in new.regions}
        next_id = max(old_fixed | new_fixed, default=0) + 1
        maps = []
        for ident, fixed in ((self, old_fixed), (new, new_fixed)):
            m = {}
            for r in ident.regions:
                if r.id in fixed:
                    m[r.id] = r.id
                else:
                    m[r.id] = next_id
                    next_id += 1
            maps.append(m)
        old, new = self._remapped(maps[0]), new._remapped(maps[1])
        n_groups = max((g.id for g in old.groups), default=0)
        groups = old.groups + tuple(dataclasses.replace(g, id=g.id + n_groups)
                                    for g in new.groups)
        return dataclasses.replace(
            new, regions=tuple(sorted(old.regions + new.regions, key=lambda r: r.id)),
            groups=groups, dropped=old.dropped + new.dropped)

    def _remapped(self, ids: Mapping[int, int]) -> "Identification":
        """The structure with the IDs of its regions changed by ``ids``."""
        def m(seq):
            return tuple(ids.get(i, i) for i in seq)
        regions = tuple(dataclasses.replace(r, id=ids.get(r.id, r.id)) for r in self.regions)
        groups = tuple(dataclasses.replace(
            g, regions=m(g.regions),
            edges=tuple(dataclasses.replace(e, regions=m(e.regions)) for e in g.edges),
            hypotheses=tuple(dataclasses.replace(h, regions=m(h.regions))
                             for h in g.hypotheses)) for g in self.groups)
        return dataclasses.replace(self, regions=regions, groups=groups)

    # --- JSON ---

    def to_dict(self) -> Dict[str, Any]:
        """The structure as a dictionary that ``json`` can write, with the
        name and the version of its schema."""
        d = _plain(dataclasses.asdict(self))
        return {"schema": SCHEMA_NAME, "schema_version": SCHEMA_VERSION, **d}

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "Identification":
        """The structure from :meth:`to_dict`. A field unknown to this
        version (written by a later one, same schema version) is ignored.

        Raises
        ------
        ValueError
            If the schema is not this one, or its version is later than
            :data:`SCHEMA_VERSION`.
        """
        if d.get("schema") != SCHEMA_NAME:
            raise ValueError(f"Not an identification: schema {d.get('schema')!r}.")
        version = d.get("schema_version")
        if not isinstance(version, int) or version > SCHEMA_VERSION:
            raise ValueError(f"The version {version} of the schema of the identification "
                             f"is not known to this Astrocook (it reads up to "
                             f"{SCHEMA_VERSION}).")

        def build(kls, raw, **convert):
            names = {f.name for f in dataclasses.fields(kls)}
            kw = {k: v for k, v in raw.items() if k in names}
            for k, fn in convert.items():
                if k in kw:
                    kw[k] = fn(kw[k])
            return kls(**kw)

        def edge(e):
            return build(Edge, e, regions=tuple, transitions=tuple)

        def comp(c):
            return build(Component, c)

        def hyp(h):
            return build(Hypothesis, h, regions=tuple, stage=Stage, status=Status, anchors=tuple,
                         components=lambda cs: tuple(comp(c) for c in cs),
                         criteria=dict)

        def group(g):
            return build(Group, g, regions=tuple,
                         edges=lambda es: tuple(edge(e) for e in es),
                         hypotheses=lambda hs: tuple(hyp(h) for h in hs))

        return build(
            cls, d, params=dict,
            regions=lambda rs: tuple(build(Region, r, kind=RegionKind,
                                           transitions=tuple) for r in rs),
            groups=lambda gs: tuple(group(g) for g in gs),
            dropped=lambda ds: tuple(build(DroppedCandidate, x) for x in ds))


def _plain(obj: Any) -> Any:
    """``obj`` with tuples as lists, enums as their values and NumPy
    scalars as Python numbers, for ``json``."""
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, StrEnum):
        return obj.value
    if hasattr(obj, "item") and callable(obj.item):  # NumPy scalar
        return obj.item()
    return obj
