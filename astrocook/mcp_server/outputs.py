"""Structured outputs of the tools.

The SDK derives the output schema of each tool from these models, and the
client shows it to the model next to the tool description: the field
descriptions are part of the interface.
"""
import json
from typing import Any, Optional

import astropy.units as au
import numpy as np
from pydantic import BaseModel, Field

from astrocook.core.session import SessionV2


class SpectrumSummary(BaseModel):
    """What the agent needs to know about a spectrum, without the data."""

    n_pixels: int = Field(description="Number of pixels.")
    x_min: float = Field(description="First pixel centre, in x_unit.")
    x_max: float = Field(description="Last pixel centre, in x_unit.")
    x_unit: str = Field(description="Unit of the x axis (e.g. nm).")
    y_unit: str = Field(description="Unit of the flux.")
    z_em: float = Field(description="Emission redshift; 0 if not set.")
    resol: float = Field(description=(
        "Resolving power R of the spectrum (median over the pixels, from "
        "its resol column); 0 if not set."))
    has_continuum: bool = Field(description="Whether a continuum is present.")
    n_components: int = Field(description="Number of absorption components.")


class SessionInfo(BaseModel):
    """A session after opening it or running a recipe on it."""

    session_id: str = Field(description=(
        "Identifier of the session, to pass to the other tools."))
    spectrum: SpectrumSummary
    log_entries: int = Field(description=(
        "Number of operations recorded in the session log."))
    saved_with_version: Optional[str] = Field(None, description=(
        "Astrocook version that saved the archive, if opened from one "
        "that records it."))
    warnings: list[str] = Field(default_factory=list, description=(
        "Warnings logged by Astrocook during the call."))
    summary: str = Field("", description=(
        "What the recipe did, in one line (e.g. the components added and "
        "removed, the columns added)."))
    notes: list[str] = Field(default_factory=list, description=(
        "What else is worth knowing: the components removed and why, those "
        "the recipe made and discarded, those whose fit failed, the steps "
        "run first to provide what the session lacked."))
    provided: list[str] = Field(default_factory=list, description=(
        "Recipes run first, with their defaults, to provide what the "
        "session lacked (provide_missing=true); each is in the session "
        "log before the recipe."))


class SaveResult(BaseModel):
    """Where a session was saved."""

    session_id: str
    path: str = Field(description="Absolute path of the .acs2 archive.")
    log_entries: int = Field(description=(
        "Number of operations recorded in the saved log."))
    warnings: list[str] = Field(default_factory=list, description=(
        "Warnings logged by Astrocook during the call."))



def _median_resolution(spec) -> float:
    """Median resolving power of the pixels (#526), 0 if unknown."""
    r = spec.resolution()
    r = r[np.isfinite(r)]
    return float(np.median(r)) if r.size else 0.0


def spectrum_summary(session: SessionV2) -> SpectrumSummary:
    """Summarize the spectrum and the components of a session."""
    spec = session.spec
    x = spec.x
    return SpectrumSummary(
        n_pixels=len(x),
        x_min=float(np.nanmin(x.value)),
        x_max=float(np.nanmax(x.value)),
        x_unit=x.unit.to_string(),
        y_unit=spec.y.unit.to_string(),
        z_em=float(spec._data.z_em),
        resol=_median_resolution(spec),
        has_continuum=spec.cont is not None,
        n_components=len(session.systs.components) if session.systs else 0,
    )


# --- Continuum ---------------------------------------------------------

LYA_REST_NM = 121.567


class ContinuumSummary(BaseModel):
    """The continuum and the mask of absorbed pixels just estimated."""

    absorbed_fraction: float = Field(description=(
        "Fraction of the pixels flagged as absorbed (column abs_mask)."))
    absorbed_fraction_forest: Optional[float] = Field(None, description=(
        "The same, bluewards of the Lyman-alpha emission (z_em); null if "
        "the spectrum has no such pixels or z_em is not set."))
    absorbed_fraction_red: Optional[float] = Field(None, description=(
        "The same, redwards of the Lyman-alpha emission."))
    snr_median: Optional[float] = Field(None, description=(
        "Median signal-to-noise ratio per pixel at the continuum level "
        "(continuum / error)."))


class ContinuumInfo(SessionInfo):
    continuum: ContinuumSummary


def _lya_split(session: SessionV2) -> Optional[np.ndarray]:
    """True for the pixels bluewards of the Lyman-alpha emission, or None
    if z_em is not set or the x axis is not a wavelength."""
    spec = session.spec
    z_em = float(spec._data.z_em)
    if z_em <= 0:
        return None
    try:
        lya = (LYA_REST_NM * (1 + z_em) * au.nm).to_value(spec.x.unit)
    except au.UnitConversionError:
        return None
    return spec.x.value < lya


def _fraction(mask: np.ndarray) -> Optional[float]:
    return float(np.mean(mask)) if mask.size else None


def continuum_summary(session: SessionV2) -> ContinuumSummary:
    """Summarize the continuum and the absorbed mask of a session."""
    spec = session.spec
    absorbed = spec.get_column("abs_mask").value.astype(bool)
    out = {"absorbed_fraction": float(np.mean(absorbed))}
    blue = _lya_split(session)
    if blue is not None:
        out["absorbed_fraction_forest"] = _fraction(absorbed[blue])
        out["absorbed_fraction_red"] = _fraction(absorbed[~blue])
    if spec.cont is not None:
        dy = spec.dy.value
        ok = np.isfinite(dy) & (dy > 0)
        if ok.any():
            out["snr_median"] = float(np.nanmedian(spec.cont.value[ok] / dy[ok]))
    return ContinuumSummary(**out)


# --- Identification ----------------------------------------------------

class Candidate(BaseModel):
    transition: str = Field(description="Transition, e.g. CIV_1548.")
    score: float = Field(description=(
        "R^2 score of the multiplet test (0-1); 0 for a Lyman-alpha "
        "assigned to a forest region by default, without a test."))


class Region(BaseModel):
    region_id: int
    x_min: float = Field(description="Start of the region, in x_unit.")
    x_max: float = Field(description="End of the region, in x_unit.")
    candidates: list[Candidate]


class SystemCandidate(BaseModel):
    series: str = Field(description="Series, e.g. CIV or Ly_a.")
    z: float
    logN_guess: float = Field(description=(
        "Initial guess of log10 column density (cm^-2), not a fit."))
    b_guess: float = Field(description=(
        "Initial guess of the Doppler parameter (km/s), not a fit."))


class IdentificationSummary(BaseModel):
    """The identifications just made; no component is created."""

    n_regions: int = Field(description="Number of identified regions.")
    n_systems: int = Field(description=(
        "Number of candidate systems (a region can hold several after "
        "de-blending)."))
    systems_by_series: dict[str, int]
    absorbed_fraction_identified: Optional[float] = Field(None, description=(
        "Fraction of the absorbed pixels (abs_mask) that lie in an "
        "identified region; the rest is absorption not accounted for."))
    regions: list[Region] = Field(description=(
        "Identified regions, best score first (at most max_listed)."))
    systems: list[SystemCandidate] = Field(description=(
        "Candidate systems, metals first, then by redshift (at most "
        "max_listed)."))
    truncated: bool = Field(description=(
        "True if regions or systems were cut to max_listed."))


class IdentificationInfo(SessionInfo):
    identification: IdentificationSummary


def _json_meta(meta: dict, key: str) -> dict:
    """A map stored by identify_lines in the metadata, as JSON or already
    decoded (after an .acs2 round trip, F27)."""
    value = meta.get(key) or {}
    return json.loads(value) if isinstance(value, str) else dict(value)


def identification_summary(session: SessionV2, max_listed: int
                           ) -> IdentificationSummary:
    """Summarize the output of identify_lines, which lives in the
    ``abs_ids`` column and in JSON strings in the metadata."""
    spec = session.spec
    meta = spec.meta
    x = spec.x.value
    ids = spec.get_column("abs_ids").value.astype(int)
    region_info = _json_meta(meta, "region_identifications")
    series = _json_meta(meta, "series_map_json")
    z = _json_meta(meta, "z_map_json")
    logN = _json_meta(meta, "logN_map_json")
    b = _json_meta(meta, "b_map_json")

    regions = []
    for rid, cands in region_info.items():
        sel = ids == int(rid)
        if not sel.any():
            continue
        regions.append(Region(
            region_id=int(rid), x_min=float(x[sel].min()),
            x_max=float(x[sel].max()),
            candidates=[Candidate(transition=t, score=s) for t, s in cands]))
    regions.sort(key=lambda r: -max((c.score for c in r.candidates), default=0))

    systems = [SystemCandidate(series=series[k], z=z[k], logN_guess=logN[k],
                               b_guess=b[k]) for k in series]
    systems.sort(key=lambda s: (s.series.startswith("Ly"), s.z))
    by_series: dict[str, int] = {}
    for s in systems:
        by_series[s.series] = by_series.get(s.series, 0) + 1

    identified = None
    if spec.has_aux_column("abs_mask"):
        absorbed = spec.get_column("abs_mask").value.astype(bool)
        if absorbed.any():
            identified = float(np.mean(ids[absorbed] > 0))
    return IdentificationSummary(
        n_regions=len(regions), n_systems=len(systems),
        systems_by_series=by_series,
        absorbed_fraction_identified=identified,
        regions=regions[:max_listed], systems=systems[:max_listed],
        truncated=len(regions) > max_listed or len(systems) > max_listed)


# --- Fit ---------------------------------------------------------------

class Component(BaseModel):
    uuid: str = Field(description="Identifier of the component.")
    series: str
    z: float
    dz: Optional[float]
    logN: float = Field(description="log10 column density (cm^-2).")
    dlogN: Optional[float]
    b: float = Field(description="Doppler parameter (km/s).")
    db: Optional[float]
    chi2r: Optional[float] = Field(description=(
        "Reduced chi2 of the fit group of the component."))
    flags: list[str] = Field(description=(
        "Warnings left by the fit: 'z_at_window', 'logN_at_lower_bound', "
        "'logN_at_upper_bound', 'b_at_lower_bound', 'b_at_upper_bound' (the "
        "parameter ended at a bound of the fit, so its value and error are "
        "not reliable; the upper bounds of the automatic fits are max_logN "
        "and max_b), 'nan_errors' (the uncertainties could not be computed)."))


class FitSummary(BaseModel):
    n_components: int = Field(description=(
        "Number of components in the session after the call."))
    n_added: int = Field(description="Components added by the call.")
    n_removed: int = Field(description=(
        "Components that were in the session before the call and are no "
        "longer (the recipe clears the region it fits)."))
    n_flagged: int = Field(description=(
        "Components added by the call with at least one flag."))
    components: list[Component] = Field(description=(
        "Components added by the call, by redshift (at most max_listed)."))
    truncated: bool = Field(description=(
        "True if the components were cut to max_listed."))


class FitInfo(SessionInfo):
    fit: FitSummary


def component_flags(c) -> list[str]:
    """Flags left on a component by its last fit (#517, see
    :class:`Component`)."""
    return list(c.flags)


def _optional(value) -> Optional[float]:
    return None if value is None else float(value)


def fit_summary(before: SessionV2, after: SessionV2, max_listed: int) -> FitSummary:
    """Compare the components before and after a fitting recipe."""
    old = {c.uuid for c in before.systs.components} if before.systs else set()
    comps = after.systs.components if after.systs else []
    new_uuids = {c.uuid for c in comps}
    added = sorted((c for c in comps if c.uuid not in old), key=lambda c: c.z)
    listed = []
    for c in added:
        listed.append(Component(
            uuid=c.uuid, series=c.series, z=c.z, dz=_optional(c.dz),
            logN=c.logN, dlogN=_optional(c.dlogN), b=c.b, db=_optional(c.db),
            chi2r=_optional(c.chi2), flags=component_flags(c)))
    return FitSummary(
        n_components=len(comps), n_added=len(added),
        n_removed=len(old - new_uuids),
        n_flagged=sum(bool(c.flags) for c in listed),
        components=listed[:max_listed], truncated=len(listed) > max_listed)


# --- The result of any recipe (phase 6-ter, PR 1) -----------------------
#
# Built from the RecipeResult in the same way for every recipe: nothing
# here names a recipe. The summaries of the pilot (continuum,
# identification, fit) are chosen by what the recipe changed or provides.

class ComponentRemoved(BaseModel):
    uuid: str
    series: str
    z: float
    logN: float
    b: float
    reason: Optional[str] = Field(None, description=(
        "Why the recipe removed it, when it says so."))


class ComponentChanged(BaseModel):
    uuid: str
    series: str
    z: float = Field(description="Redshift after the call.")
    changed: dict[str, list] = Field(description=(
        "The parameters that changed, each as [before, after]."))


class ChangesSummary(BaseModel):
    """What the recipe changed in the session, listed (at most max_listed
    items per list; the counts are complete)."""

    n_components_added: int = 0
    n_components_removed: int = 0
    n_components_changed: int = 0
    components_removed: list[ComponentRemoved] = Field(default_factory=list, description=(
        "Components that were in the session and are no longer, with the "
        "reason when the recipe gives one."))
    components_changed: list[ComponentChanged] = Field(default_factory=list)
    columns_added: list[str] = Field(default_factory=list)
    columns_removed: list[str] = Field(default_factory=list)
    columns_changed: list[str] = Field(default_factory=list)
    properties_changed: dict[str, list] = Field(default_factory=dict, description=(
        "Properties of the session that changed (name, z_em, resolution "
        "as the median R, units, number of pixels), each as [before, after]."))
    truncated: bool = Field(False, description=(
        "True if a list was cut to max_listed."))


class NewSession(BaseModel):
    """A session opened by a recipe that branches (e.g. a region extracted
    from the spectrum)."""

    session_id: str = Field(description="Identifier of the new session.")
    spectrum: SpectrumSummary
    summary: str = Field(description="What changed, compared with the original.")


class RecipeInfo(SessionInfo):
    """The session after a recipe, and what the recipe did."""

    recipe: str = Field(description="The recipe that ran, as menu.recipe.")
    changes: Optional[ChangesSummary] = Field(None, description=(
        "What the recipe changed in this session; null if it opened new "
        "sessions instead (new_sessions) or returned a value (value)."))
    new_sessions: list[NewSession] = Field(default_factory=list, description=(
        "Sessions opened by the recipe; this session is unchanged."))
    value: Optional[Any] = Field(None, description=(
        "The result of a recipe that returns a value instead of changing "
        "the session (e.g. a redshift)."))
    not_fitted: list[str] = Field(default_factory=list, description=(
        "Components whose fit group failed: their values are those from "
        "before the call."))
    continuum: Optional[ContinuumSummary] = Field(None, description=(
        "The continuum and the absorbed mask, when the recipe made or "
        "changed them."))
    identification: Optional[IdentificationSummary] = Field(None, description=(
        "The identified regions and candidate systems, when the recipe "
        "identified lines."))
    fit: Optional[FitSummary] = Field(None, description=(
        "The components added by the recipe, with errors, reduced chi2 and "
        "flags, when it added any."))
    details: dict[str, Any] = Field(default_factory=dict, description=(
        "What is specific to the recipe, as it reports it: e.g. what "
        "optimize_system tried and why it stopped (optimization), the "
        "doublets that identify_lines scored and dropped (dropped_candidates)."))
    undo_hint: Optional[str] = Field(None, description=(
        "How to have back what the recipe removed, when it removed "
        "components or columns."))


def _jsonable(value):
    """A value of a change, as JSON can carry it."""
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return str(value)


def changes_summary(changes, removed_reasons: dict, max_listed: int) -> ChangesSummary:
    """The :class:`~astrocook.recipes.result.Changes` of a recipe, listed."""
    removed = [ComponentRemoved(uuid=c.uuid, series=c.series, z=c.z, logN=c.logN,
                                b=c.b, reason=removed_reasons.get(c.uuid))
               for c in changes.components_removed]
    changed = [ComponentChanged(uuid=c.component.uuid, series=c.component.series,
                                z=c.component.z,
                                changed={k: [_jsonable(a), _jsonable(b)]
                                         for k, (a, b) in c.changed.items()})
               for c in changes.components_changed]
    return ChangesSummary(
        n_components_added=len(changes.components_added),
        n_components_removed=len(removed), n_components_changed=len(changed),
        components_removed=removed[:max_listed],
        components_changed=changed[:max_listed],
        columns_added=list(changes.columns_added),
        columns_removed=list(changes.columns_removed),
        columns_changed=list(changes.columns_changed),
        properties_changed={k: [_jsonable(a), _jsonable(b)]
                            for k, (a, b) in changes.properties_changed.items()},
        truncated=len(removed) > max_listed or len(changed) > max_listed)


#: Keys of ``RecipeResult.extra`` that the output says in its own fields
_EXTRA_SAID_ELSEWHERE = frozenset({"removed", "discarded", "not_fitted", "candidates",
                                   "regions"})


def details(extra: dict) -> dict:
    """The rest of ``RecipeResult.extra``: what a recipe reports of its own,
    as JSON, for every recipe in the same way."""
    return {k: _deep_jsonable(v) for k, v in extra.items()
            if k not in _EXTRA_SAID_ELSEWHERE}


def _deep_jsonable(value):
    if isinstance(value, dict):
        return {str(k): _deep_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_deep_jsonable(v) for v in value]
    return _jsonable(value)


def _value(value):
    """The value of a recipe that does not return a session."""
    if isinstance(value, (np.generic,)):
        value = value.item()
    return _jsonable(value)


def summaries(result, before: SessionV2, after: SessionV2, provides,
              max_listed: int) -> dict:
    """The summaries of the pilot that fit what the recipe did: the
    continuum when it provides one or changed it, the identification when
    it identified lines, the fit when it added components."""
    out = {}
    changes = result.changes
    touched = set(changes.columns_added) | set(changes.columns_changed)
    if after.spec.cont is not None and after.spec.has_aux_column("abs_mask") and (
            "continuum" in provides or touched & {"cont", "abs_mask"}):
        out["continuum"] = continuum_summary(after)
    if ("candidates" in result.extra or "regions" in result.extra) \
            and after.spec.has_aux_column("abs_ids") \
            and "region_identifications" in after.spec.meta:
        out["identification"] = identification_summary(after, max_listed)
    if changes.components_added:
        out["fit"] = fit_summary(before, after, max_listed)
    return out


# --- The sessions (phase 6-ter, PR 2) -------------------------------------

class SessionListItem(BaseModel):
    session_id: str
    origin: str = Field(description=(
        "'opened' (from a file), 'branched' (by a recipe that opens new "
        "sessions) or 'duplicated' (session_duplicate)."))
    parent: Optional[str] = Field(None, description=(
        "The session it was made from, for a branch or a copy (it may be "
        "closed since)."))
    path: Optional[str] = Field(None, description="The file it was opened from.")
    spectrum: SpectrumSummary
    n_states: int = Field(description=(
        "States in its history: the session as opened, then one per recipe."))
    current_state: int = Field(description=(
        "Index of the current state (0: as opened); lower than n_states - 1 "
        "after session_undo."))
    unsaved: bool = Field(description=(
        "True if the current state, or the log (a note), is not the one last "
        "saved or opened."))
    running: Optional[str] = Field(None, description=(
        "The recipe running on it and for how long, if any: until it ends, "
        "the tools read and save the state before it."))


class SessionList(BaseModel):
    sessions: list[SessionListItem]


class ColumnInfo(BaseModel):
    name: str
    unit: str = Field(description="Unit; empty if dimensionless.")


class NoteInfo(BaseModel):
    kind: str = Field(description=(
        "'note' or 'gap': written with session_note; 'gap_seen': a gap seen "
        "by the server (a prerequisite that no recipe estimates, a tool that "
        "failed on a valid request, a recipe undone)."))
    text: str


def note_info(note) -> NoteInfo:
    return NoteInfo(kind=note.kind, text=note.text)


class LogEntryInfo(BaseModel):
    notes_before: list[NoteInfo] = Field(default_factory=list, description=(
        "The notes written at the state this recipe ran on, before it."))
    recipe: str = Field(description="The recipe, as in the log (its name).")
    params: dict[str, Any] = Field(description="The values it used.")
    reason: Optional[str] = Field(None, description=(
        "Why it was run, if given (the 'reason' of a recipe tool)."))


class StateInfo(BaseModel):
    index: int = Field(description="0: the session as opened.")
    recipe: Optional[str] = Field(None, description=(
        "The recipe that made this state from the previous one; null for "
        "the first."))
    params: dict[str, Any] = Field(default_factory=dict)
    reason: Optional[str] = Field(None, description=(
        "Why the recipe was run, if given (the 'reason' of a recipe tool)."))
    summary: str = Field("", description=(
        "What it changed, compared with the previous state."))
    notes: list[NoteInfo] = Field(default_factory=list, description=(
        "The notes written at this state (session_note)."))
    current: bool = Field(False, description="True for the current state.")


class NoteResult(BaseModel):
    """A note added to the log of a session."""

    session_id: str
    state: int = Field(description="The state the note was written at (the current one).")
    note: NoteInfo
    unsaved: bool = Field(description=(
        "True: the note is in the session, saved with it by session_save."))


class MethodsResult(BaseModel):
    """The log of a session as the text of a methods section."""

    session_id: str
    text: str = Field(description=(
        "The text: an opening paragraph, one per recipe of the log (with its "
        "reason, what it did and changed, and the notes written there), the "
        "fit and the atomic data, and the references."))
    paragraphs: list[str] = Field(description="The paragraphs of the text, in order.")
    references: list[str] = Field(description="The works the text cites.")


class SessionDescription(SessionListItem):
    """A session in full, without the data."""

    columns: list[ColumnInfo] = Field(description="The columns of the spectrum.")
    components_by_series: dict[str, int] = Field(description=(
        "Number of absorption components per series (e.g. CIV)."))
    log_before: list[LogEntryInfo] = Field(default_factory=list, description=(
        "The log of the archive the session was opened from: what was done "
        "before its first state, for reading only (it is not replayed, and "
        "session_save does not save it again: an archive holds the log from "
        "its original data to its state)."))
    saved_path: Optional[str] = Field(None, description=(
        "Where it was last saved or opened from."))


class HistoryInfo(BaseModel):
    session_id: str
    states: list[StateInfo] = Field(description=(
        "The states of the session, oldest first; session_undo and "
        "session_redo move between them."))
    can_undo: bool
    can_redo: bool
    truncated: bool = Field(False, description=(
        "True if only the last max_listed states are listed."))


class UndoResult(SessionInfo):
    """The session after session_undo or session_redo."""

    moved: list[str] = Field(description=(
        "The recipes undone (or redone), in the order they were moved over."))
    current_state: int
    n_states: int


class CloseResult(BaseModel):
    session_id: str
    closed: bool
    open_sessions: list[str]


class GuiResult(BaseModel):
    """The GUI launched on a session (``gui_open``)."""

    session_id: str
    path: str = Field(description=(
        "The archive the GUI follows, saved again after every change of the "
        "session; not the user's own save (use session_save for that)."))
    command: list[str] = Field(description="The command that launched the GUI.")


class GapGroup(BaseModel):
    kind: str = Field(description=(
        "'declared' (by an agent), 'prerequisite' (needed and no recipe "
        "estimates it), 'failed' (a tool failed on a valid request), "
        "'undone' (a recipe whose result was undone)."))
    text: str = Field(description="The first text of the group.")
    count: int
    sessions: int = Field(description="How many different sessions.")
    versions: list[str]
    first: str = Field(description="Date of the first.")
    last: str = Field(description="Date of the last.")
    examples: list[str] = Field(description="Up to three different texts of the group.")


class GapsSummary(BaseModel):
    """The gaps of the collection of recipes, grouped."""

    enabled: bool = Field(description="False if the record is turned off.")
    groups: list[GapGroup]
    n_groups: int = Field(description="Groups before max_groups cut the list.")
    skipped_lines: int = Field(description="Lines of the file that could not be read.")
