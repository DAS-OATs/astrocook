"""The tools of the recipes, generated from the registry (phase 6-ter, PR 1).

One tool per recipe of :mod:`astrocook.recipes.registry`, named
``<menu>_<recipe>`` (agreed in the pilot), with:

- the schema of its parameters from :func:`~astrocook.mcp_server.schemas.parameters_schema`
  (the parameters offered to agents), plus ``session_id``,
  ``provide_missing`` and ``max_listed``;
- its description from the ``brief`` and ``details`` of the registry, and
  what it requires and provides;
- its output, a :class:`~astrocook.mcp_server.outputs.RecipeInfo`, built
  from the :class:`~astrocook.recipes.result.RecipeResult` in the same way
  for every recipe.

Nothing here names a recipe, except :data:`EXCLUDED`: a recipe added to the
registry reaches the agents with no code of the server (the principle of
the phase, ``ROADMAP.md``, phase 6-ter). The values of the parameters are
checked and converted by the registry, as for the GUI and the scripts
(``typed_recipe``): a wrong value is refused with the message of the recipe.
"""
from typing import Annotated, Any, Callable

from mcp.server.mcpserver import Context, MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.server.mcpserver.utilities.func_metadata import ArgModelBase
from pydantic import ConfigDict, Field

from astrocook.mcp_server.errors import run_recipe
from astrocook.mcp_server.outputs import (
    NewSession, RecipeInfo, SessionInfo, _value, changes_summary,
    spectrum_summary, summaries)
from astrocook.mcp_server.progress import report_progress
from astrocook.mcp_server.schemas import parameters_schema
from astrocook.mcp_server.sessions import SessionStore
from astrocook.recipes.errors import MissingPrerequisite, RecipeError
from astrocook.recipes.prerequisites import missing as missing_prerequisites
from astrocook.recipes.registry import RecipeSpec, all_recipes
from astrocook.recipes.result import run as run_registry

#: The recipes that are not tools, each with its reason. A recipe of the
#: registry that is neither a tool nor here fails a test.
EXCLUDED: dict[str, str] = {}

#: The parameters of every recipe tool, besides those of the recipe
SERVER_PARAMS = {
    "session_id": {"type": "string", "description": (
        "Identifier of the session to run the recipe on, as returned by "
        "session_open or by a recipe that opened a new session.")},
    "provide_missing": {"type": "boolean", "default": False, "description": (
        "Provide what the session lacks and a recipe can make with its "
        "default parameters (e.g. a continuum), before running this one; "
        "the step is reported in 'provided' and logged. False: the tool "
        "says what is missing instead.")},
    "max_listed": {"type": "integer", "minimum": 0, "default": 50, "description": (
        "Maximum number of items in each list of the output (components, "
        "regions, systems); the counts are always complete.")},
    "reason": {"type": ["string", "null"], "default": None, "description": (
        "Why you run this recipe with these parameters (e.g. 'z_em from the "
        "peak of Ly-a'), in a sentence or two: stored with the step in the "
        "session log, for the person who reviews it. A recipe that only "
        "returns a value logs nothing.")},
}

#: In the result of a recipe that removed components or columns
UNDO_HINT = ("session_undo restores the state before this recipe, with what it "
             "removed.")

#: Added to the message of a prerequisite that a recipe can provide
PROVIDE_HINT = (" Run it first, or call this tool again with provide_missing=true to "
                "provide it with its defaults.")

_NEEDS = {"z_em": "the emission redshift (z_em)", "resolution": "the resolution",
          "continuum": "a continuum"}


def as_tools(message: str) -> str:
    """A message of Astrocook with the recipes it names, ``(menu.recipe)``,
    as the tools that run them, ``(tool menu_recipe)``."""
    for spec in recipe_specs():
        message = message.replace(f"({spec.category}.{spec.name})",
                                  f"(tool {tool_name(spec)})")
    return message


def tool_name(spec: RecipeSpec) -> str:
    """``<menu>_<recipe>``, e.g. ``continuum_estimate_auto``."""
    return f"{spec.category}_{spec.name}"


def is_tool(spec: RecipeSpec) -> bool:
    return f"{spec.category}.{spec.name}" not in EXCLUDED


def tool_description(spec: RecipeSpec) -> str:
    """The brief and details of the registry, then what the recipe
    requires and provides, and whether it opens new sessions."""
    text = f"{spec.brief} {spec.details}".strip()
    more = []
    if spec.requires:
        needs = ", ".join(_NEEDS.get(r, r) for r in sorted(spec.requires))
        more.append(f"Needs {needs}: if the session lacks it, the tool says what "
                    f"is missing and which tool provides it.")
    if spec.provides:
        gives = ", ".join(_NEEDS.get(p, p) for p in sorted(spec.provides))
        more.append(f"Provides {gives}.")
    if spec.branching:
        more.append("Opens new sessions (listed in new_sessions, each with its "
                    "session_id); the session it runs on does not change.")
    return text + ("\n\n" + " ".join(more) if more else "")


def input_schema(spec: RecipeSpec) -> dict:
    """The parameters of the recipe offered to agents, with those of the
    server."""
    schema = parameters_schema(spec)
    clash = set(schema["properties"]) & set(SERVER_PARAMS)
    if clash:
        raise ValueError(f"The recipe {spec.category}.{spec.name} has parameters "
                         f"named as those of the server: {sorted(clash)}")
    schema["properties"] = {"session_id": SERVER_PARAMS["session_id"],
                            **schema["properties"],
                            "provide_missing": SERVER_PARAMS["provide_missing"],
                            "max_listed": SERVER_PARAMS["max_listed"],
                            "reason": SERVER_PARAMS["reason"]}
    schema["required"] = ["session_id", *schema["required"]]
    return schema


class _RecipeArgs(ArgModelBase):
    """The arguments of a recipe tool: those of the server are fields, the
    parameters of the recipe are kept as given (the registry checks them)."""

    model_config = ConfigDict(extra="allow", arbitrary_types_allowed=True)

    def model_dump_one_level(self) -> dict[str, Any]:
        out = super().model_dump_one_level()
        out["recipe_params"] = dict(self.model_extra or {})
        return out


def _register(server: MCPServer, fn: Callable, name: str, description: str,
              schema: dict) -> None:
    """Add ``fn`` as the tool ``name``, with ``schema`` as the schema of its
    parameters, and the given parameters not declared by ``fn`` passed to it
    as ``recipe_params``.

    The SDK derives the schema from the signature of the function; the
    tools of the recipes take theirs from the registry, so the registered
    tool is adjusted after (``_tool_manager``, private in ``mcp`` 2.x: the
    tests of the tools fail if it changes)."""
    server.add_tool(fn, name=name, description=description, structured_output=True)
    tool = server._tool_manager.get_tool(name)
    tool.parameters = schema
    base = tool.fn_metadata.arg_model
    tool.fn_metadata.arg_model = type(base.__name__, (base, _RecipeArgs),
                                      {"model_config": _RecipeArgs.model_config})


def _record_gap(store: SessionStore, spec: RecipeSpec, session_id: str,
                session, error: ToolError) -> None:
    """Record the failure of a recipe tool as a gap, if it is one (phase
    6-ter, PR 5b, decision 6). What is a mistake of the agent is not:
    unknown parameters and sessions (refused before), a value the registry
    refuses (``InvalidParameter``), what a recipe of the collection provides
    (``MissingPrerequisite``, the agent runs it) and what the data do not
    allow (``RecipeFailed``, with a reason the agent can act on). A gap is
    a prerequisite that no recipe estimates from the data (z_em, the
    resolution: they must be given), and a failure that is none of the
    above: the request was valid and the recipe stopped with an unexpected
    error (usually a bug). The message of the error is not recorded (it may
    hold a path)."""
    cause = error.__cause__
    tool = tool_name(spec)
    if isinstance(cause, MissingPrerequisite):
        for gap in missing_prerequisites(session, spec):
            if gap.automatic is None:
                store.record_gap(
                    session_id, "prerequisite",
                    f"{tool} needs {_NEEDS.get(gap.name, gap.name)}, which no recipe "
                    f"estimates from the data.", once=True)
    elif cause is not None and not isinstance(cause, (RecipeError, ToolError)):
        store.record_gap(
            session_id, "failed",
            f"{tool} failed on a valid request: {type(cause).__name__}.", once=True)


def run_tool(store: SessionStore, spec: RecipeSpec, session_id: str,
             params: dict, provide_missing: bool, max_listed: int,
             ctx: Context | None = None, reason: str | None = None) -> RecipeInfo:
    """Run the recipe of a tool and build its output; ``ctx``, the context
    of the call, receives the progress of the recipe (``progress.py``);
    ``reason`` goes to the log entry of the recipe (phase 6-ter, PR 5a)."""
    offered = {p.name for p in spec.params if p.agents}
    unknown = sorted(set(params) - offered)
    if unknown:
        raise ToolError(f"Unknown parameters for {tool_name(spec)}: "
                        f"{', '.join(unknown)}. Its parameters: "
                        f"{', '.join(sorted(offered)) or 'none'}.")
    recipe = f"{spec.category}.{spec.name}"
    with store.lock:
        history = store.get(session_id)
        before = history.current_state

        def call():
            try:
                return run_registry(before, recipe,
                                    on_missing="provide" if provide_missing else "raise",
                                    **params)
            except MissingPrerequisite as e:
                hint = PROVIDE_HINT if "provides it with its defaults" in str(e) else ""
                raise ToolError(as_tools(str(e)) + hint) from e

        try:
            with report_progress(ctx, tool_name(spec)):
                result, _ = run_recipe(tool_name(spec), call)
        except ToolError as e:
            _record_gap(store, spec, session_id, before, e)
            raise
        removed = dict(result.extra.get("removed", {}))
        out: dict[str, Any] = {"recipe": recipe,
                               "not_fitted": [str(u) for u in result.extra.get("not_fitted", ())]}
        if result.branching and result.session is not None:
            all_changes = (result.changes if isinstance(result.changes, tuple)
                           else (result.changes,))
            new_ids = store.branch(session_id, result.sessions, spec.name, result.params,
                                   provided=result.provided, reason=reason)
            out["new_sessions"] = [
                NewSession(session_id=i, spectrum=spectrum_summary(s),
                           summary="; ".join(c.describe()) or "no change")
                for i, s, c in zip(new_ids, result.sessions, all_changes, strict=True)]
            if len(new_ids) == 1:
                out.update(summaries(result, before, result.session, spec.provides,
                                     max_listed))
        elif result.session is not None:
            store.apply(session_id, result.session, spec.name, result.params,
                        provided=result.provided, reason=reason)
            out["changes"] = changes_summary(result.changes, removed, max_listed)
            if result.changes.components_removed or result.changes.columns_removed:
                # Decision 1 of the phase: what a recipe removed can be had
                # back from the history of the session
                out["undo_hint"] = UNDO_HINT
            out.update(summaries(result, before, result.session, spec.provides,
                                 max_listed))
        else:
            # A value (e.g. a redshift): the session does not change, and
            # nothing is logged, as in the GUI
            out["value"] = _value(result.value)
        return RecipeInfo(**session_info(store, session_id, result).model_dump(), **out)


def session_info(store: SessionStore, session_id: str, result=None) -> SessionInfo:
    """The session after a call, with what the recipe said."""
    history = store.get(session_id)
    return SessionInfo(
        session_id=session_id,
        spectrum=spectrum_summary(history.current_state),
        log_entries=len(history.log_manager.to_dict()["entries"]),
        saved_with_version=history.states[0].saved_with_version,
        warnings=list(result.warnings) if result is not None else [],
        summary=result.summary if result is not None else "",
        notes=list(result.notes) if result is not None else [],
        provided=[p.recipe for p in result.provided] if result is not None else [],
    )


def _make_tool(store: SessionStore, spec: RecipeSpec) -> Callable:
    def tool(session_id: Annotated[str, Field()],
             provide_missing: bool = False,
             max_listed: int = 50,
             reason: str | None = None,
             recipe_params: dict | None = None,
             ctx: Context | None = None) -> RecipeInfo:
        return run_tool(store, spec, session_id, recipe_params or {},
                        provide_missing, max_listed, ctx, reason)
    tool.__name__ = tool_name(spec)
    return tool


def recipe_specs():
    """The recipes offered as tools, in the order of the registry."""
    return [s for s in all_recipes() if is_tool(s)]


def register_recipe_tools(server: MCPServer, store: SessionStore) -> list[str]:
    """Add a tool for each recipe of the registry that is not excluded, and
    return their names."""
    names = []
    for spec in recipe_specs():
        name = tool_name(spec)
        _register(server, _make_tool(store, spec), name, tool_description(spec),
                  input_schema(spec))
        names.append(name)
    return names


def menus_text() -> str:
    """The menus of the recipes and their tools, for the instructions of
    the server."""
    by_menu: dict[str, list[str]] = {}
    for spec in recipe_specs():
        by_menu.setdefault(spec.category, []).append(spec.name)
    return "\n".join(f"- {menu}: {', '.join(f'{menu}_{n}' for n in names)}"
                     for menu, names in by_menu.items())
