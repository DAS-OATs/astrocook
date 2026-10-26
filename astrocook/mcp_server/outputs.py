"""Structured outputs of the tools.

The SDK derives the output schema of each tool from these models, and the
client shows it to the model next to the tool description: the field
descriptions are part of the interface.
"""
import json
from typing import Optional

import astropy.units as au
import numpy as np
from pydantic import BaseModel, Field

from astrocook.core.session import SessionV2
from astrocook.fitting.voigt_model import VoigtModelConstraintV2


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

# Bounds of the fits: logN and the lower bound on b from the constraint
# model (VoigtModelConstraintV2.DEFAULT_BOUNDS), the upper bounds on b and
# logN from the max_b and max_logN of the recipe
LOGN_BOUNDS = VoigtModelConstraintV2.DEFAULT_BOUNDS["logN"]
B_MIN = VoigtModelConstraintV2.DEFAULT_BOUNDS["b"][0]
#: Relative tolerance to call a parameter "at a bound"
AT_BOUND_RTOL = 1e-3


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
        "Warnings on the fit: 'nan_errors'; 'logN_at_lower_bound', "
        "'logN_at_upper_bound', 'logN_at_max_logN', 'b_at_lower_bound', "
        "'b_at_max_b' (the "
        "parameter ended at a bound of the fit, so its value and error "
        "are not reliable)."))


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


def _close(value: float, bound: float) -> bool:
    return abs(value - bound) <= AT_BOUND_RTOL * max(1.0, abs(bound))


def component_flags(c, max_b: float, max_logN: Optional[float] = None) -> list[str]:
    """Flags of a fitted component (see :class:`Component`)."""
    flags = []
    errors = (c.dz, c.dlogN, c.db)
    if any(e is not None and not np.isfinite(e) for e in errors):
        flags.append("nan_errors")
    if _close(c.logN, LOGN_BOUNDS[0]):
        flags.append("logN_at_lower_bound")
    if _close(c.logN, LOGN_BOUNDS[1]):
        flags.append("logN_at_upper_bound")
    elif max_logN is not None and _close(c.logN, max_logN):
        flags.append("logN_at_max_logN")
    if _close(c.b, B_MIN):
        flags.append("b_at_lower_bound")
    if _close(c.b, max_b):
        flags.append("b_at_max_b")
    return flags


def _optional(value) -> Optional[float]:
    return None if value is None else float(value)


def fit_summary(before: SessionV2, after: SessionV2, max_b: float,
                max_listed: int, max_logN: Optional[float] = None) -> FitSummary:
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
            chi2r=_optional(c.chi2), flags=component_flags(c, max_b, max_logN)))
    return FitSummary(
        n_components=len(comps), n_added=len(added),
        n_removed=len(old - new_uuids),
        n_flagged=sum(bool(c.flags) for c in listed),
        components=listed[:max_listed], truncated=len(listed) > max_listed)
