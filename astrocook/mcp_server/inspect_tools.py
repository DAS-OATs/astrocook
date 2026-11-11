"""The tools to look into a session (phase 6-ter, PR 3): the components,
the spectrum in a range, a plot.

The pilot showed that the agents need to see the data to judge a result,
and worked around the server to do it (Python on the FITS file: section
4.2 of ``ROADMAP_HISTORY.md``). These tools read the current state of a session;
they change nothing and log nothing. What they compute is in the package,
for scripts too (phase 6-ter, PR 3, Guido): the fit groups in
:meth:`~astrocook.core.system_list.SystemList.fit_groups` (those of
Refit All), the reading of a range in
:meth:`~astrocook.core.spectrum.Spectrum.columns_in_range`, the plot in
:mod:`astrocook.plotting`.
"""
from typing import Annotated, Optional

import numpy as np
from mcp.server.mcpserver import Image, MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from pydantic import BaseModel, Field

from astrocook.mcp_server.outputs import _jsonable
from astrocook.recipes._spectrum import spec_of
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
    splits the large ones (:meth:`SystemList.fit_groups` at the depth of
    the fits). A component pulled in by two groups gets the first."""
    systs = session.systs
    if not systs or not systs.components:
        return {}
    groups = systs.fit_groups([c.uuid for c in systs.components],
                              max_depth=FIT_EXPANSION_DEPTH, spec=session.spec)
    group_of: dict[str, int] = {}
    for n, members in enumerate(groups):
        for u in members:
            group_of.setdefault(u, n)
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
        with store.state_lock:
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
            "Columns to read (session_describe lists them); the names of "
            "the GUI work too (F = y, dF = dy, λ = x); null: y, dy, cont "
            "and model, those present."))] = None,
        max_points: Annotated[int, Field(ge=10, le=5000, description=(
            "Maximum number of points: a longer range is averaged in blocks "
            "of pixels."))] = 200,
    ) -> SpectrumRange:
        """Read the columns of the spectrum of a session in a range of x, as
        numbers. A long range is averaged in blocks of pixels (said in
        binned_by), so that it fits in the context of a model; for a look
        at the shape, spectrum_plot is cheaper."""
        with store.state_lock:
            spec = spec_of(store.get(session_id).current_state)
            if columns is None:
                columns = [c for c in DEFAULT_COLUMNS if c in spec.t.colnames]
            try:
                data = spec.columns_in_range(columns, x_min, x_max, max_points)
            except KeyError as e:
                raise ToolError(f"{e.args[0]} (session '{session_id}')") from None
            except ValueError as e:
                raise ToolError(str(e)) from None
            x = data["columns"]["x"]
            return SpectrumRange(
                session_id=session_id, x_unit=str(spec.x.unit),
                x_min=float(np.nanmin(x)), x_max=float(np.nanmax(x)),
                n_pixels=data["n_pixels"], binned_by=data["binned_by"],
                columns={k: [_jsonable(float(v)) for v in vals]
                         for k, vals in data["columns"].items()},
                units=data["units"])

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
        from astrocook.plotting import figure_png, plot_spectrum
        with store.state_lock:
            session = store.get(session_id).current_state
            try:
                fig = plot_spectrum(session, x_min, x_max, normalized)
            except ValueError as e:
                raise ToolError(str(e)) from None
            return [Image(data=figure_png(fig), format="png"),
                    f"{session_id}: spectrum plotted"
                    + (" (normalized)" if normalized and spec_of(session).cont is not None else "")]
