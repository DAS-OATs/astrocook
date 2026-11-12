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
from astrocook.mcp_server.gui_tools import register_gui_tools
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

The log is what the person who reviews your work reads: give each recipe
tool a 'reason' (why this recipe, why these values), and write with
session_note what you did or decided outside the tools (e.g. a value
measured with your own code, a range you chose by eye). Notes make no new
state and the replay skips them. session_methods writes the log as the
text of a methods section, with your reasons and notes in it.

Whenever you work outside the tools because no recipe does what the analysis
needs (e.g. z_em measured with your own code), declare a gap: session_note
with kind='gap', saying what was needed and what you did instead. These
gaps, with the tools that failed on a valid request and the recipes undone
(give session_undo a 'reason'), are collected in a local file of the user,
never sent anywhere, to decide which recipes to add; gaps_summary lists them.

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

A long recipe sends its progress, and a line every 20 s when it is silent.
To stop one, cancel its call; if your client cannot, session_cancel does
it. A recipe cancelled stops at its next step and leaves the session as it
was before it.

To look into a session: session_describe (its columns), components_list
(the components with their errors, flags, fit groups and constraints),
spectrum_read (the columns in a range, as numbers) and spectrum_plot (a
plot of a range, with the continuum, the model and the components). They
change nothing.

The columns have the names of the GUI everywhere (reads, expressions,
parameters): λ (wavelength), F (flux), dF (its error), cont, model, ...;
λ can be written lam or x, and x, y, dy (the names in the data and in
the CSV files) are accepted too.

To let the user watch the analysis, gui_open shows a session in the Astrocook
GUI and keeps it up to date after every change (read only for the user;
nothing they do there reaches you). It is not available where there is no
display.
"""

#: Typical workflows on the components. The tools of Absorbers overlap, and
#: what an agent most needs to know is which of them fit and which do not
#: (read from the recipes, phase 6-quater PR 12: update_component and
#: update_constraint change values and links and run no fit; the model drawn
#: stays as it was until the next fit, F152). Each example is a call written
#: as ``tool(parameter=value, ...)``: a test reads them and checks the names
#: of the tools and of the parameters against the registry.
INSTRUCTIONS_WORKFLOWS = """
Components, by hand. These tools fit, with the components blended with the
one they name: absorbers_add_component, absorbers_add_linked_system,
absorbers_fit_component; absorbers_refit_all fits every group;
absorbers_optimize_system adds components where the residuals call for them,
keeping each only if the AIC improves. These do not fit: absorbers_update_component,
absorbers_update_components, absorbers_update_constraint,
absorbers_update_constraints_column change values and links only, and the
model stays as it was until a fit: follow them with absorbers_fit_component.
- Adapt a system: components_list(series='CIV') gives the uuids and the fit
  of each component; absorbers_update_component(uuid='<uuid>', logN=14.2);
  absorbers_fit_component(uuid='<uuid>'); judge it with components_list and
  spectrum_plot; absorbers_optimize_system(uuid='<uuid>') if the residuals
  need more components.
- Try a component: absorbers_add_component(series='CIV', z=1.9, logN=13.5,
  b=10.0) adds and fits it; keep it, or session_undo to drop it (a try on a
  copy: session_duplicate).
- Link ions: absorbers_add_linked_system(series_list='CIV,SiIV', z=1.9) adds
  new components with the redshifts tied to the first; to tie a parameter of
  components that exist, absorbers_update_constraint(uuid='<uuid of SiIV>',
  param='b', expression="p['<uuid of CIV>'].b"), then
  absorbers_fit_component(uuid='<uuid of SiIV>'); is_free=True frees it again.
"""


def instructions() -> str:
    """The instructions of the server: the menus of the recipes come from
    the registry."""
    return (INSTRUCTIONS_HEAD + menus_text() + "\n" + INSTRUCTIONS_TAIL
            + INSTRUCTIONS_WORKFLOWS)


def build_server() -> MCPServer:
    """Return a new server, with its own (empty) set of sessions, and a
    tool for each recipe of the registry at that time."""
    server = MCPServer("astrocook", instructions=instructions(),
                       version=__version__, log_level="WARNING")
    store = SessionStore()

    register_session_tools(server, store)
    register_inspect_tools(server, store)
    register_gui_tools(server, store)
    register_recipe_tools(server, store)
    return server


def main() -> None:
    """Entry point of ``astrocook-mcp``: serve over stdio."""
    build_server().run()


if __name__ == "__main__":
    main()
