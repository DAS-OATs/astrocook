"""The tools to look into a session (phase 6-ter, PR 3): the components,
the spectrum in a range, a plot.

The pilot showed that the agents need to see the data to judge a result,
and worked around the server to do it (Python on the FITS file: section
4.2 of ``ROADMAP.md``). These tools read the current state of a session;
they change nothing and log nothing.
"""
import io
from typing import Annotated, Optional

import astropy.units as au
import numpy as np
from mcp.server.mcpserver import Image, MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from pydantic import BaseModel, Field

from astrocook.core.system_list import get_all_w_rest
from astrocook.mcp_server.outputs import _jsonable
from astrocook.mcp_server.session_tools import SESSION_ID
from astrocook.mcp_server.sessions import SessionStore
from astrocook.recipes.absorbers import FIT_EXPANSION_DEPTH

C_KMS = 299792.458


# --- Components ----------------------------------------------------------------

class ComponentInfo(BaseModel):
    uuid: str
    series: str
    z: float
    dz: Optional[float]
    logN: float = Field(description="log10 column density (cm^-2).")
    dlogN: Optional[float]
    b: float = Field(description="Doppler parameter (km/s).")
    db: Optional[float]
    btur: float = Field(description="Turbulent broadening (km/s); 0 if not used.")
    chi2r: Optional[float] = Field(description=(
        "Reduced chi2 of the fit group of the component, at its last fit."))
    flags: list[str] = Field(description=(
        "Warnings left by the last fit (a parameter at a bound, errors "
        "that could not be computed, 'not_fitted')."))
    group: int = Field(description=(
        "The fit group: components with the same number are fitted "
        "together (their lines overlap or they are linked), as Refit All "
        "groups them."))
    constraints: dict[str, str] = Field(description=(
        "The parameters that are not free: 'fixed', or 'linked: <expression>' "
        "(e.g. the z of a component linked to another)."))
    dv_kms: Optional[float] = Field(None, description=(
        "Velocity offset (km/s) from the first component of its group with "
        "the same series."))


class ComponentList(BaseModel):
    session_id: str
    n_components: int = Field(description="Components in the session.")
    n_matched: int = Field(description="Components that match the filters.")
    n_groups: int = Field(description="Fit groups in the session.")
    components: list[ComponentInfo] = Field(description=(
        "The components that match, by redshift (at most max_listed, from "
        "offset)."))
    truncated: bool


def fit_groups(session) -> dict[str, int]:
    """The fit group of each component: the groups of Refit All before it
    splits the large ones (seeds by increasing z, the friends of friends of
    ``get_connected_group`` at the depth of the fits)."""
    systs = session.systs
    if not systs or not systs.components:
        return {}
    comps = sorted(systs.components, key=lambda c: (c.z, c.series))
    group_of: dict[str, int] = {}
    n = 0
    for c in comps:
        if c.uuid in group_of:
            continue
        members = systs.get_connected_group([c.uuid], max_depth=FIT_EXPANSION_DEPTH,
                                            spec=session.spec)
        for u in members | {c.uuid}:
            group_of.setdefault(u, n)
        n += 1
    return group_of


def _constraints(systs, uuid: str) -> dict[str, str]:
    out = {}
    for param, constr in systs._data.v2_constraints_map.get(uuid, {}).items():
        if constr.target_uuid or constr.expression:
            out[param] = f"linked: {constr.expression or constr.target_uuid}"
        elif not constr.is_free:
            out[param] = "fixed"
    return out


def _optional(value) -> Optional[float]:
    return None if value is None or not np.isfinite(value) else float(value)


# --- Spectrum ------------------------------------------------------------------

class SpectrumRange(BaseModel):
    session_id: str
    x_unit: str
    x_min: float
    x_max: float
    n_pixels: int = Field(description="Pixels of the spectrum in the range.")
    binned_by: int = Field(description=(
        "Pixels averaged into each point listed (1: every pixel); for "
        "reading only, the session is not rebinned."))
    columns: dict[str, list[Optional[float]]] = Field(description=(
        "The values of the columns asked, one per point; null where not "
        "finite. 'x' is always given."))
    units: dict[str, str]


DEFAULT_COLUMNS = ("y", "dy", "cont", "model")


def _range_mask(spec, x_min, x_max) -> np.ndarray:
    x = spec.x.value
    lo = np.nanmin(x) if x_min is None else x_min
    hi = np.nanmax(x) if x_max is None else x_max
    if lo >= hi:
        raise ToolError(f"Empty range: x_min ({lo}) must be less than x_max ({hi}).")
    sel = (x >= lo) & (x <= hi)
    if not sel.any():
        raise ToolError(f"No pixel between {lo} and {hi} {spec.x.unit}: the spectrum "
                        f"covers {np.nanmin(x):.4g}-{np.nanmax(x):.4g} {spec.x.unit}.")
    return sel


def _binned(values: np.ndarray, n: int, error: bool) -> np.ndarray:
    """Means over blocks of ``n`` (an error column: the error of the mean)."""
    if n == 1:
        return values
    m = len(values) // n * n
    blocks = values[:m].reshape(-1, n)
    rest = values[m:]
    if error:
        out = np.sqrt(np.nanmean(blocks ** 2, axis=1) / n)
        tail = [np.sqrt(np.nanmean(rest ** 2) / len(rest))] if len(rest) else []
    else:
        out = np.nanmean(blocks, axis=1)
        tail = [np.nanmean(rest)] if len(rest) else []
    return np.concatenate([out, tail])


def register_inspect_tools(server: MCPServer, store: SessionStore) -> None:
    """Add the tools that read a session to ``server``."""

    @server.tool()
    def components_list(
        session_id: SESSION_ID,
        series: Annotated[Optional[str], Field(description=(
            "Only these series, comma-separated (e.g. 'CIV, SiIV')."))] = None,
        z_min: Annotated[Optional[float], Field(description="Only z >= z_min.")] = None,
        z_max: Annotated[Optional[float], Field(description="Only z <= z_max.")] = None,
        flagged_only: Annotated[bool, Field(description=(
            "Only the components with a flag."))] = False,
        offset: Annotated[int, Field(ge=0, description=(
            "Skip the first matching components (to read a long list in "
            "pages)."))] = 0,
        max_listed: Annotated[int, Field(ge=0, description=(
            "Maximum number of components listed."))] = 50,
    ) -> ComponentList:
        """List the absorption components of a session, with their errors,
        the reduced chi2 and flags of their last fit, their fit group and
        their constraints, by redshift."""
        with store.lock:
            session = store.get(session_id).current_state
            systs = session.systs
            comps = list(systs.components) if systs else []
            groups = fit_groups(session)
            wanted = ({s.strip() for s in series.split(",") if s.strip()}
                      if series else None)
            matched = [c for c in sorted(comps, key=lambda c: (c.z, c.series))
                       if (wanted is None or c.series in wanted)
                       and (z_min is None or c.z >= z_min)
                       and (z_max is None or c.z <= z_max)
                       and (not flagged_only or c.flags)]
            first: dict[tuple, float] = {}
            for c in sorted(comps, key=lambda c: c.z):
                first.setdefault((groups.get(c.uuid), c.series), c.z)
            listed = []
            for c in matched[offset:offset + max_listed]:
                z0 = first[(groups.get(c.uuid), c.series)]
                listed.append(ComponentInfo(
                    uuid=c.uuid, series=c.series, z=c.z, dz=_optional(c.dz),
                    logN=c.logN, dlogN=_optional(c.dlogN), b=c.b, db=_optional(c.db),
                    btur=c.btur, chi2r=_optional(c.chi2), flags=list(c.flags),
                    group=groups.get(c.uuid, -1), constraints=_constraints(systs, c.uuid),
                    dv_kms=(c.z - z0) / (1 + z0) * C_KMS))
            return ComponentList(
                session_id=session_id, n_components=len(comps), n_matched=len(matched),
                n_groups=len(set(groups.values())), components=listed,
                truncated=offset + max_listed < len(matched))

    @server.tool()
    def spectrum_read(
        session_id: SESSION_ID,
        x_min: Annotated[Optional[float], Field(description=(
            "Start of the range, in the unit of the x axis (spectrum.x_unit "
            "of the session, e.g. nm); null: the start of the spectrum."))] = None,
        x_max: Annotated[Optional[float], Field(description=(
            "End of the range; null: the end of the spectrum."))] = None,
        columns: Annotated[Optional[list[str]], Field(description=(
            "Columns to read (session_describe lists them); null: y, dy, "
            "cont and model, those present."))] = None,
        max_points: Annotated[int, Field(ge=10, le=5000, description=(
            "Maximum number of points: a longer range is averaged in blocks "
            "of pixels."))] = 200,
    ) -> SpectrumRange:
        """Read the columns of the spectrum of a session in a range of x, as
        numbers. A long range is averaged in blocks of pixels (said in
        binned_by), so that it fits in the context of a model; for a look
        at the shape, spectrum_plot is cheaper."""
        with store.lock:
            spec = store.get(session_id).current_state.spec
            names = spec.t.colnames
            if columns is None:
                columns = [c for c in DEFAULT_COLUMNS if c in names]
            missing = [c for c in columns if c not in names]
            if missing:
                raise ToolError(f"No column {', '.join(missing)} in '{session_id}'. "
                                f"Its columns: {', '.join(names)}.")
            sel = _range_mask(spec, x_min, x_max)
            n = int(sel.sum())
            binned_by = int(np.ceil(n / max_points))
            out, units = {}, {}
            for name in ["x", *[c for c in columns if c != "x"]]:
                values = np.asarray(spec.t[name], dtype=float)[sel]
                values = _binned(values, binned_by, error=name.startswith("d"))
                out[name] = [_jsonable(float(v)) for v in values]
                units[name] = str(spec.t[name].unit or "")
            x = spec.x.value[sel]
            return SpectrumRange(session_id=session_id, x_unit=str(spec.x.unit),
                                 x_min=float(x.min()), x_max=float(x.max()), n_pixels=n,
                                 binned_by=binned_by, columns=out, units=units)

    @server.tool()
    def spectrum_plot(
        session_id: SESSION_ID,
        x_min: Annotated[Optional[float], Field(description=(
            "Start of the range, in the unit of the x axis; null: the start "
            "of the spectrum."))] = None,
        x_max: Annotated[Optional[float], Field(description=(
            "End of the range; null: the end of the spectrum."))] = None,
        normalized: Annotated[bool, Field(description=(
            "Divide by the continuum (if any)."))] = False,
    ) -> list:
        """Plot the spectrum of a session in a range of x, as a PNG image:
        the flux and its error, the continuum and the model (if any), and a
        tick at each line of each component, labelled with its series."""
        with store.lock:
            session = store.get(session_id).current_state
            return [Image(data=_plot(session, x_min, x_max, normalized), format="png"),
                    f"{session_id}: spectrum plotted"
                    + (" (normalized)" if normalized and session.spec.cont is not None else "")]


def _plot(session, x_min, x_max, normalized) -> bytes:
    """The PNG of :func:`spectrum_plot` (matplotlib, without a display)."""
    from matplotlib.figure import Figure

    spec = session.spec
    sel = _range_mask(spec, x_min, x_max)
    x = spec.x.value[sel]
    cont = spec.cont.value[sel] if spec.cont is not None else None
    norm = cont if (normalized and cont is not None) else 1.0
    fig = Figure(figsize=(10, 3.6), dpi=100)
    ax = fig.add_subplot()
    ax.step(x, spec.y.value[sel] / norm, where="mid", lw=0.7, color="black", label="flux")
    ax.step(x, spec.dy.value[sel] / norm, where="mid", lw=0.6, color="tab:gray",
            label="error")
    if cont is not None:
        ax.plot(x, cont / norm, lw=1, color="tab:blue", label="continuum")
    if spec.has_aux_column("model"):
        ax.plot(x, spec.get_column("model").value[sel] / norm, lw=1, color="tab:red",
                label="model")
    top = np.nanpercentile(spec.y.value[sel] / norm, 99) if sel.any() else 1
    comps = session.systs.components if session.systs else []
    for c in comps:
        for w in get_all_w_rest(c.series):
            xo = (w * (1 + c.z) * au.AA).to_value(spec.x.unit, equivalencies=au.spectral())
            if x.min() <= xo <= x.max():
                ax.axvline(xo, ymin=0.9, ymax=0.97, color="tab:green", lw=0.8)
                ax.text(xo, top * 1.12, c.series, fontsize=6, rotation=90,
                        ha="center", va="bottom", color="tab:green")
    ax.set_xlim(x.min(), x.max())
    ax.set_xlabel(f"x ({spec.x.unit})")
    ax.set_ylabel("normalized flux" if np.ndim(norm) else f"flux ({spec.y.unit})")
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    return buf.getvalue()
