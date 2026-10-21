"""The MCP server and its tools (pilot, ``ROADMAP.md`` section 4.1).

Each tool is named after what it calls: ``<menu>_<recipe>`` for a recipe,
``session_<method>`` for a method of :class:`~astrocook.core.session.SessionV2`.
The docstrings of the tools are their descriptions for the model.

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
from astrocook.mcp_server.outputs import SaveResult, SessionInfo, spectrum_summary
from astrocook.mcp_server.sessions import SessionStore

INSTRUCTIONS = """\
Astrocook analyses quasar spectra: continuum, identification of absorption
systems, Voigt fitting. Open a spectrum with session_open and use the
session_id it returns in the other tools; each tool returns a summary of
the session, not the data. Save the result with session_save, so that a
person can review it in the Astrocook GUI. Paths are on the machine where
the server runs; give them absolute.
"""


def _resolve(path: str) -> Path:
    return Path(os.path.expanduser(path)).resolve()


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
                current = store.get(session_id).current_state
                new_state, more = run_recipe(
                    "edit_set_properties", current.edit.set_properties,
                    **params)
                store.apply(session_id, new_state, "set_properties", params)
                warnings += more
            return info(session_id, warnings)

    @server.tool()
    def session_save(
        session_id: Annotated[str, Field(description=(
            "Identifier returned by session_open."))],
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

    return server


def main() -> None:
    """Entry point of ``astrocook-mcp``: serve over stdio."""
    build_server().run()


if __name__ == "__main__":
    main()
