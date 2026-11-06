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
from mcp.server.mcpserver import MCPServer

from astrocook import __version__
from astrocook.mcp_server.inspect_tools import register_inspect_tools
from astrocook.mcp_server.session_tools import register_session_tools
from astrocook.mcp_server.sessions import SessionStore
from astrocook.mcp_server.tools import menus_text, register_recipe_tools

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

Every state of a session is kept: session_history lists them, with what
each recipe changed; session_undo goes back (e.g. to have back the
components a recipe removed) and session_redo forward, without running
anything again. To try something on a copy, use session_duplicate.
session_list lists the open sessions; session_close closes one, and
refuses if it has changes not saved.

To look into a session: session_describe (its columns), components_list
(the components with their errors, flags, fit groups and constraints),
spectrum_read (the columns in a range, as numbers) and spectrum_plot (a
plot of a range, with the continuum, the model and the components). They
change nothing.
"""


def instructions() -> str:
    """The instructions of the server: the menus of the recipes come from
    the registry."""
    return INSTRUCTIONS_HEAD + menus_text() + "\n" + INSTRUCTIONS_TAIL


def build_server() -> MCPServer:
    """Return a new server, with its own (empty) set of sessions, and a
    tool for each recipe of the registry at that time."""
    server = MCPServer("astrocook", instructions=instructions(),
                       version=__version__, log_level="WARNING")
    store = SessionStore()

    register_session_tools(server, store)
    register_inspect_tools(server, store)
    register_recipe_tools(server, store)
    return server


def main() -> None:
    """Entry point of ``astrocook-mcp``: serve over stdio."""
    build_server().run()


if __name__ == "__main__":
    main()
