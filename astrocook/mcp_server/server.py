"""The MCP server and its tools (pilot, ``ROADMAP.md`` section 4.1).

Each tool is named after what it calls: ``<menu>_<recipe>`` for a recipe,
``session_<method>`` for a method of :class:`~astrocook.core.session.SessionV2`.
The docstrings of the tools are their descriptions for the model.

The defaults and the parameter descriptions of the recipe tools are read
from the recipe schemas (``*_RECIPES_SCHEMAS``, the ones the GUI dialogs
use), so that the tools and the GUI share them (F40: the method
signatures may differ). Parameters that need the GUI (debug plots) or
the internals (column names) are left out.

Launch with ``astrocook-mcp`` (``stdio`` transport); a client starts it as
a subprocess, e.g. ``claude mcp add astrocook -- astrocook-mcp``.
"""
import os
from pathlib import Path
from typing import Annotated, Optional

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from pydantic import Field

from astrocook import __version__
from astrocook.core.session import SessionV2
from astrocook.mcp_server.errors import run_recipe
from astrocook.mcp_server.outputs import (
    ContinuumInfo, FitInfo, IdentificationInfo, SaveResult, SessionInfo,
    continuum_summary, fit_summary, identification_summary, spectrum_summary)
from astrocook.mcp_server.sessions import SessionStore
from astrocook.recipes.absorbers import ABSORBERS_RECIPES_SCHEMAS
from astrocook.recipes.continuum import CONTINUUM_RECIPES_SCHEMAS

INSTRUCTIONS = """\
Astrocook analyses quasar spectra: continuum, identification of absorption
systems, Voigt fitting. Open a spectrum with session_open and use the
session_id it returns in the other tools; each tool returns a summary of
the session, not the data. Save the result with session_save, so that a
person can review it in the Astrocook GUI. Paths are on the machine where
the server runs; give them absolute.
"""


SCHEMAS = {"continuum": CONTINUUM_RECIPES_SCHEMAS,
           "absorbers": ABSORBERS_RECIPES_SCHEMAS}

SESSION_ID = Annotated[str, Field(description=(
    "Identifier returned by session_open."))]


def _resolve(path: str) -> Path:
    return Path(os.path.expanduser(path)).resolve()


def _param(menu: str, recipe: str, name: str) -> dict:
    """The schema entry of a recipe parameter."""
    for p in SCHEMAS[menu][recipe]["params"]:
        if p["name"] == name:
            return p
    raise KeyError(f"{menu}.{recipe} has no parameter {name}")


def default(menu: str, recipe: str, name: str):
    """Default of a recipe parameter, as in the GUI dialog."""
    return _param(menu, recipe, name)["default"]


def doc(menu: str, recipe: str, name: str) -> str:
    """Description of a recipe parameter, as in the GUI dialog."""
    return _param(menu, recipe, name)["doc"]


def description(menu: str, recipe: str, more: str) -> str:
    """Tool description: the recipe brief and details, then ``more``."""
    schema = SCHEMAS[menu][recipe]
    return f"{schema['brief']} {schema['details']}\n\n{more}"


def build_server() -> MCPServer:
    """Return a new server, with its own (empty) set of sessions."""
    server = MCPServer("astrocook", instructions=INSTRUCTIONS,
                       version=__version__, log_level="WARNING")
    store = SessionStore()

    def info(session_id: str, warnings: list[str]) -> SessionInfo:
        history = store.get(session_id)
        state = history.current_state
        return SessionInfo(
            session_id=session_id,
            spectrum=spectrum_summary(state),
            log_entries=len(history.log_manager.to_dict()["entries"]),
            saved_with_version=history.states[0].saved_with_version,
            warnings=warnings,
        )

    def run(session_id: str, menu: str, recipe: str, params: dict):
        """Run a recipe on the current state of a session, record it in
        the log and return the new state with the warnings."""
        current = store.get(session_id).current_state
        new_state, warnings = run_recipe(
            f"{menu}_{recipe}", getattr(getattr(current, menu), recipe),
            **params)
        store.apply(session_id, new_state, recipe, params)
        return new_state, warnings

    @server.tool()
    def session_open(
        path: Annotated[str, Field(description=(
            "Absolute path of the file: a spectrum (FITS, ESO/MIDAS, ...) "
            "or an Astrocook archive (.acs2, or V1 .acs)."))],
        z_em: Annotated[Optional[float], Field(description=(
            "Emission redshift of the quasar, if known; needed by the "
            "continuum. Leave out to keep the value in the file."))] = None,
        resol: Annotated[Optional[float], Field(description=(
            "Resolving power R (e.g. 45000), if known; needed by the "
            "fits. Leave out to keep the value in the file."))] = None,
    ) -> SessionInfo:
        """Open a spectrum or an Astrocook archive as a new session.

        The format is detected from the file. If z_em or resol are given,
        they are set with the recipe edit.set_properties, which is
        recorded in the session log.
        """
        file_path = _resolve(path)
        if not file_path.is_file():
            raise ToolError(f"No such file: {file_path}")
        with store.lock:
            session, warnings = run_recipe(
                "session_open", SessionV2.open_new, str(file_path),
                file_path.stem, None, "auto")
            session_id = store.add(session, str(file_path))
            if z_em is not None or resol is not None:
                params = {}
                if z_em is not None:
                    params["z_em"] = str(z_em)
                if resol is not None:
                    params["resol"] = str(resol)
                warnings += run(session_id, "edit", "set_properties",
                                params)[1]
            return info(session_id, warnings)

    @server.tool()
    def session_save(
        session_id: SESSION_ID,
        path: Annotated[str, Field(description=(
            "Absolute path of the archive; .acs2 is appended if missing."))],
        overwrite: Annotated[bool, Field(description=(
            "Replace the file if it exists."))] = False,
    ) -> SaveResult:
        """Save a session as an Astrocook archive (.acs2).

        The archive holds the current state, the original data and the
        session log, and opens in the Astrocook GUI.
        """
        file_path = _resolve(path)
        if file_path.suffix.lower() != ".acs2":
            file_path = file_path.with_name(file_path.name + ".acs2")
        if file_path.exists() and not overwrite:
            raise ToolError(f"{file_path} exists; pass overwrite=true to "
                            "replace it.")
        if not file_path.parent.is_dir():
            raise ToolError(f"No such folder: {file_path.parent}")
        with store.lock:
            history = store.get(session_id)
            # save() returns 0 on success, unlike the recipes
            _, warnings = run_recipe(
                "session_save", lambda: history.current_state.save(
                    str(file_path), initial_session=history.states[0]) or True)
            return SaveResult(
                session_id=session_id, path=str(file_path),
                log_entries=len(history.log_manager.to_dict()["entries"]),
                warnings=warnings)

    @server.tool(description=description("continuum", "estimate_auto", (
        "Needs z_em (set it with session_open), to tell the Lyman-alpha "
        "forest from the red side. Replaces any existing continuum and "
        "mask. Returns the fraction of absorbed pixels and the S/N.")))
    def continuum_estimate_auto(
        session_id: SESSION_ID,
        smooth_len_lya: Annotated[float, Field(description=doc(
            "continuum", "estimate_auto", "smooth_len_lya"))] = default("continuum", "estimate_auto", "smooth_len_lya"),
        smooth_len_out: Annotated[float, Field(description=doc(
            "continuum", "estimate_auto", "smooth_len_out"))] = default("continuum", "estimate_auto", "smooth_len_out"),
        kappa: Annotated[float, Field(description=doc(
            "continuum", "estimate_auto", "kappa"))] = default("continuum", "estimate_auto", "kappa"),
        fudge: Annotated[str, Field(description=doc(
            "continuum", "estimate_auto", "fudge"))] = default("continuum", "estimate_auto", "fudge"),
        smooth_fwhm: Annotated[float, Field(description=doc(
            "continuum", "estimate_auto", "smooth_fwhm"))] = default("continuum", "estimate_auto", "smooth_fwhm"),
        template: Annotated[bool, Field(description=doc(
            "continuum", "estimate_auto", "template"))] = default("continuum", "estimate_auto", "template"),
    ) -> ContinuumInfo:
        params = {"smooth_len_lya": str(smooth_len_lya),
                  "smooth_len_out": str(smooth_len_out), "kappa": str(kappa),
                  "fudge": str(fudge), "smooth_fwhm": str(smooth_fwhm),
                  "template": str(template)}
        with store.lock:
            state, warnings = run(session_id, "continuum", "estimate_auto",
                                  params)
            return ContinuumInfo(**info(session_id, warnings).model_dump(),
                                 continuum=continuum_summary(state))

    @server.tool(description=description("absorbers", "identify_lines", (
        "Needs z_em; uses the absorbed mask of continuum_estimate_auto (or "
        "makes one with the default parameters). Creates no component: "
        "returns the identified regions with their candidate transitions "
        "and scores, and the candidate systems with initial guesses of "
        "logN and b.")))
    def absorbers_identify_lines(
        session_id: SESSION_ID,
        multiplets: Annotated[str, Field(description=doc(
            "absorbers", "identify_lines", "multiplets"))] = default("absorbers", "identify_lines", "multiplets"),
        min_pix_region: Annotated[int, Field(description=doc(
            "absorbers", "identify_lines", "min_pix_region"))] = default("absorbers", "identify_lines", "min_pix_region"),
        merge_dv: Annotated[float, Field(description=doc(
            "absorbers", "identify_lines", "merge_dv"))] = default("absorbers", "identify_lines", "merge_dv"),
        score_threshold: Annotated[float, Field(description=doc(
            "absorbers", "identify_lines", "score_threshold"))] = default("absorbers", "identify_lines", "score_threshold"),
        bypass_scoring: Annotated[bool, Field(description=doc(
            "absorbers", "identify_lines", "bypass_scoring"))] = default("absorbers", "identify_lines", "bypass_scoring"),
        max_listed: Annotated[int, Field(ge=0, description=(
            "Maximum number of regions and of systems listed in the "
            "output; the counts are always complete."))] = 50,
    ) -> IdentificationInfo:
        # mask_col (a column name) and debug_rating (plots in the GUI)
        # are left out
        params = {"multiplets": multiplets,
                  "min_pix_region": str(min_pix_region),
                  "merge_dv": str(merge_dv),
                  "score_threshold": str(score_threshold),
                  "bypass_scoring": str(bypass_scoring)}
        with store.lock:
            state, warnings = run(session_id, "absorbers", "identify_lines",
                                  params)
            return IdentificationInfo(
                **info(session_id, warnings).model_dump(),
                identification=identification_summary(state, max_listed))

    @server.tool(description=description("absorbers", "doublets_auto", (
        "Identifies the doublets redwards of the Lyman-alpha emission, adds "
        "a component for each candidate and fits them (Voigt profiles, "
        "least squares). Needs z_em and the resolution (session_open) and "
        "a continuum (made with the default parameters if missing). It "
        "makes a new absorbed mask with the default parameters, replacing "
        "the one of continuum_estimate_auto, and removes the components "
        "already on the red side before fitting. Can take minutes on a "
        "large spectrum. Returns the components added, with errors, "
        "reduced chi2 and flags for the ones at a bound of the fit.")))
    def absorbers_doublets_auto(
        session_id: SESSION_ID,
        multiplets: Annotated[str, Field(description=doc(
            "absorbers", "doublets_auto", "multiplets"))] = default("absorbers", "doublets_auto", "multiplets"),
        score_threshold: Annotated[float, Field(description=doc(
            "absorbers", "doublets_auto", "score_threshold"))] = default("absorbers", "doublets_auto", "score_threshold"),
        z_window_kms: Annotated[float, Field(description=doc(
            "absorbers", "doublets_auto", "z_window_kms"))] = default("absorbers", "doublets_auto", "z_window_kms"),
        max_b: Annotated[float, Field(description=doc(
            "absorbers", "doublets_auto", "max_b"))] = default("absorbers", "doublets_auto", "max_b"),
        max_listed: Annotated[int, Field(ge=0, description=(
            "Maximum number of components listed in the output; the "
            "counts are always complete."))] = 50,
    ) -> FitInfo:
        params = {"multiplets": multiplets,
                  "score_threshold": str(score_threshold),
                  "z_window_kms": str(z_window_kms), "max_b": str(max_b)}
        with store.lock:
            before = store.get(session_id).current_state
            state, warnings = run(session_id, "absorbers", "doublets_auto",
                                  params)
            return FitInfo(**info(session_id, warnings).model_dump(),
                           fit=fit_summary(before, state, max_b, max_listed))

    return server


def main() -> None:
    """Entry point of ``astrocook-mcp``: serve over stdio."""
    build_server().run()


if __name__ == "__main__":
    main()
