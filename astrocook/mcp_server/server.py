"""The MCP server and its tools (``ROADMAP.md``, section 4; phase 6-ter).

Each tool is named after what it calls: ``<menu>_<recipe>`` for a recipe,
``session_<method>`` for the sessions. The tools of the recipes are
generated from the recipe registry (:mod:`astrocook.mcp_server.tools`,
phase 6-ter PR 1): their parameters, descriptions and outputs come from the
registry and from :class:`~astrocook.recipes.result.RecipeResult`, so that
a recipe added to the registry reaches the agents with no code here. The
tools of the sessions are written here; the docstrings of the tools are
their descriptions for the model.

What a session lacks (z_em, the resolution, a continuum) is said, not
provided, unless the agent asks with ``provide_missing`` (phase 6-bis,
PR 8: the agent decides).

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
from astrocook.mcp_server.outputs import SaveResult, SessionInfo
from astrocook.mcp_server.sessions import SessionStore
from astrocook.mcp_server.tools import (
    menus_text, register_recipe_tools, run_tool, session_info)
from astrocook.recipes.registry import get_recipe

INSTRUCTIONS_HEAD = """\
Astrocook analyses quasar spectra: continuum, identification of absorption
systems, Voigt fitting. Open a spectrum with session_open and use the
session_id it returns in the other tools; each tool returns a summary of
the session and of what changed, not the data. Save the result with
session_save, so that a person can review it in the Astrocook GUI. Paths
are on the machine where the server runs; give them absolute.

Each recipe of Astrocook is a tool, named <menu>_<recipe>, as in the menus
of the GUI:
"""

INSTRUCTIONS_TAIL = """
A recipe runs on the current state of a session and makes the next one,
recorded in the session log with its parameters (the log replays). A
recipe that opens new sessions (e.g. a region extracted) leaves the
session unchanged and lists the new ones, each with its session_id.

A recipe that needs something the session lacks (the emission redshift,
the resolution, a continuum) does not run: the tool says what is missing
and which tool provides it. z_em and the resolution are set with
session_open or edit_set_properties; a continuum is made with
continuum_estimate_auto, or, with provide_missing=true, with its default
parameters before the recipe (the step is reported in 'provided' and
recorded in the session log).
"""


def instructions() -> str:
    """The instructions of the server: the menus of the recipes come from
    the registry."""
    return INSTRUCTIONS_HEAD + menus_text() + "\n" + INSTRUCTIONS_TAIL


SESSION_ID = Annotated[str, Field(description=(
    "Identifier returned by session_open."))]


def _resolve(path: str) -> Path:
    return Path(os.path.expanduser(path)).resolve()


def build_server() -> MCPServer:
    """Return a new server, with its own (empty) set of sessions, and a
    tool for each recipe of the registry at that time."""
    server = MCPServer("astrocook", instructions=instructions(),
                       version=__version__, log_level="WARNING")
    store = SessionStore()

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
                file_path.stem, store.context, "auto")
            session_id = store.add(session, str(file_path))
        params = {}
        if z_em is not None:
            params["z_em"] = z_em
        if resol is not None:
            params["resol"] = str(resol)
        if params:
            # Through the tool of the recipe (it takes the lock): logged, as
            # in the GUI
            result = run_tool(store, get_recipe("edit", "set_properties"), session_id,
                              params, False, 0)
            warnings += result.warnings
        return session_info(store, session_id).model_copy(update={"warnings": warnings})

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

    register_recipe_tools(server, store)
    return server


def main() -> None:
    """Entry point of ``astrocook-mcp``: serve over stdio."""
    build_server().run()


if __name__ == "__main__":
    main()
