"""The tools of the sessions (``session_<method>``).

Open and save (the pilot), and since phase 6-ter (PR 2) list, describe,
duplicate and close, and move in the history of a session (decision 1 of
the phase: the server keeps every state, as the GUI does, so that what a
recipe removed can be had back with ``session_undo``, without running
anything again).

These tools are written here, not generated: they are about the sessions,
not the recipes, which reach the server from the registry
(:mod:`astrocook.mcp_server.tools`).
"""
import os
from collections import Counter
from pathlib import Path
from typing import Annotated, Literal, Optional

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from pydantic import Field

from astrocook import gaps
from astrocook.core.session import SessionV2
from astrocook.mcp_server.errors import run_recipe
from astrocook.mcp_server.outputs import (
    CloseResult, ColumnInfo, GapGroup, GapsSummary, HistoryInfo, LogEntryInfo, NoteResult, SaveResult,
    SessionDescription, SessionInfo, SessionList, SessionListItem, StateInfo,
    UndoResult, _jsonable, note_info, spectrum_summary)
from astrocook.mcp_server.sessions import SessionStore
from astrocook.mcp_server.tools import run_tool, session_info, tool_name
from astrocook.recipes.registry import find_recipe, get_recipe
from astrocook.recipes.result import compare_sessions

SESSION_ID = Annotated[str, Field(description=(
    "Identifier of the session, as returned by session_open, by a recipe "
    "that opened a new session, or by session_duplicate."))]


def _resolve(path: str) -> Path:
    return Path(os.path.expanduser(path)).resolve()


def _params(params: dict) -> dict:
    return {k: _jsonable(v) for k, v in params.items()}


def list_item(store: SessionStore, session_id: str) -> dict:
    """The fields of :class:`SessionListItem` for a session."""
    history = store.get(session_id)
    meta = store.meta[session_id]
    return dict(session_id=session_id, origin=meta.origin, parent=meta.parent,
                path=meta.path, spectrum=spectrum_summary(history.current_state),
                n_states=len(history.states), current_state=history.current_index,
                unsaved=store.unsaved(session_id))


def _moved_recipe(store: SessionStore, session_id: str, index: int) -> str:
    entry = store.entry_of(session_id, index)
    return entry["recipe_name"] if entry else "?"


def register_session_tools(server: MCPServer, store: SessionStore) -> None:
    """Add the tools of the sessions to ``server``."""

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
            store.mark_saved(session_id, str(file_path))
            return SaveResult(
                session_id=session_id, path=str(file_path),
                log_entries=len(history.log_manager.to_dict()["entries"]),
                warnings=warnings)

    @server.tool()
    def session_list() -> SessionList:
        """List the open sessions: where each comes from, its spectrum in a
        few numbers, its history, and whether it has changes not saved."""
        with store.lock:
            return SessionList(sessions=[SessionListItem(**list_item(store, i))
                                         for i in store.ids])

    @server.tool()
    def session_describe(session_id: SESSION_ID) -> SessionDescription:
        """Describe a session: the columns of the spectrum, the components
        per series, where it comes from and the log that came with it. The
        states made in this server are listed by session_history."""
        with store.lock:
            history = store.get(session_id)
            meta = store.meta[session_id]
            state = history.current_state
            spec = state.spec
            columns = [ColumnInfo(name=c, unit=str(spec.t[c].unit or ""))
                       for c in spec.t.colnames]
            comps = state.systs.components if state.systs else []
            log = history.log_manager
            before = [LogEntryInfo(notes_before=[note_info(n) for n in log.notes_at(i)],
                                   recipe=e.recipe_name, params=_params(e.params),
                                   reason=e.reason)
                      for i, e in enumerate(log.entries[:meta.first_entry])]
            return SessionDescription(
                **list_item(store, session_id), columns=columns,
                components_by_series=dict(Counter(c.series for c in comps)),
                log_before=before, saved_path=meta.saved_path)

    @server.tool()
    def session_history(
        session_id: SESSION_ID,
        max_listed: Annotated[int, Field(ge=1, description=(
            "Maximum number of states listed, the most recent."))] = 50,
    ) -> HistoryInfo:
        """List the states of a session: the session as opened, then one per
        recipe run on it, with the recipe, its parameters and what it
        changed. session_undo and session_redo move between them."""
        with store.lock:
            history = store.get(session_id)
            states = []
            for i, state in enumerate(history.states):
                entry = store.entry_of(session_id, i)
                summary = ""
                if i > 0:
                    summary = "; ".join(
                        compare_sessions(history.states[i - 1], state).describe()
                    ) or "nothing changed"
                states.append(StateInfo(
                    index=i, recipe=entry["recipe_name"] if entry else None,
                    params=_params(entry["params"]) if entry else {},
                    reason=entry["reason"] if entry else None, summary=summary,
                    notes=[note_info(n) for n in store.notes_of(session_id, i)],
                    current=i == history.current_index))
            return HistoryInfo(session_id=session_id, states=states[-max_listed:],
                               can_undo=history.can_undo(), can_redo=history.can_redo(),
                               truncated=len(states) > max_listed)

    def undone_gap(session_id: str, index: int, why: Optional[str]) -> None:
        """Record the result of a recipe that the agent undid (a step provided
        before a recipe is not a recipe of the agent)."""
        entry = store.entry_of(session_id, index)
        if entry is None or store.was_provided(session_id, index):
            return
        spec = find_recipe(entry["recipe_name"])
        text = f"{tool_name(spec) if spec else entry['recipe_name']} undone."
        if entry["reason"]:
            text += f" It was run because: {entry['reason']}"
        if why:
            text += f" It was undone because: {why}"
        store.record_gap(session_id, "undone", text)

    def move(session_id: str, steps: int, back: bool,
             why: Optional[str] = None) -> UndoResult:
        with store.lock:
            history = store.get(session_id)
            can = history.can_undo if back else history.can_redo
            if not can():
                raise ToolError(f"Nothing to {'undo' if back else 'redo'} in "
                                f"'{session_id}': it is at state {history.current_index} "
                                f"of {len(history.states) - 1}.")
            moved, undone = [], []
            for _ in range(steps):
                if not can():
                    break
                index = history.current_index if back else history.current_index + 1
                moved.append(_moved_recipe(store, session_id, index))
                if back:
                    undone.append(index)
                history.undo() if back else history.redo()
            for index in undone:
                undone_gap(session_id, index, why)
            info = session_info(store, session_id)
            return UndoResult(**info.model_dump(), moved=moved,
                              current_state=history.current_index,
                              n_states=len(history.states))

    STEPS = Annotated[int, Field(ge=1, description="How many states to move.")]

    @server.tool()
    def session_note(
        session_id: SESSION_ID,
        text: Annotated[str, Field(min_length=1, description=(
            "The note: what you did or decided outside the tools (e.g. 'z_em "
            "= 2.345 from the peak of Ly-a, measured on the flux with my own "
            "code'), or what you saw in the session. For a gap, what you "
            "needed and what you did instead."))],
        kind: Annotated[Literal["note", "gap"], Field(description=(
            "'note' (default), or 'gap' when you work outside the tools "
            "because no recipe does what the analysis needs."))] = "note",
    ) -> NoteResult:
        """Write a note in the log of a session, at its current state: for
        what you did or decided outside the tools, so that the person who
        reviews the session can read it. It makes no new state and is
        skipped when the log is replayed; it is undone with its state, and
        saved with the session by session_save. To say why you run a
        recipe, give its tool a 'reason' instead.

        Declare a gap (kind='gap') every time you do outside the tools
        something that a recipe should do: what the analysis needed, and
        what you did instead (e.g. 'Needed z_em: no recipe estimates it.
        Measured the peak of Ly-a on the flux with my own code: 2.345').
        Gaps are collected, with the failures and the recipes undone, to
        decide which recipes to add; they are kept in a local file of the
        user and never sent anywhere. gaps_summary lists them."""
        with store.lock:
            history = store.get(session_id)
            try:
                note = store.note(session_id, text, kind)
            except ValueError as e:
                raise ToolError(str(e)) from e
            return NoteResult(session_id=session_id, state=history.current_index,
                              note=note_info(note), unsaved=store.unsaved(session_id))

    @server.tool()
    def session_undo(
        session_id: SESSION_ID, steps: STEPS = 1,
        reason: Annotated[Optional[str], Field(description=(
            "Why you undo it (e.g. 'it removed the CIV system I need'), in a "
            "sentence: recorded with the gap that an undone recipe is."))] = None,
    ) -> UndoResult:
        """Go back to the state before the last recipe (or the last
        ``steps``): its components, columns and log are those of that state.
        The states undone are kept until another recipe runs, and
        session_redo goes forward again."""
        return move(session_id, steps, back=True, why=reason)

    @server.tool()
    def session_redo(session_id: SESSION_ID, steps: STEPS = 1) -> UndoResult:
        """Go forward again over the states undone by session_undo."""
        return move(session_id, steps, back=False)

    @server.tool()
    def session_duplicate(session_id: SESSION_ID) -> SessionListItem:
        """Make a copy of a session, with its history up to the current
        state and its log, under a new identifier: to try something without
        touching the original."""
        with store.lock:
            new_id = store.duplicate(session_id)
            return SessionListItem(**list_item(store, new_id))

    @server.tool()
    def session_close(
        session_id: SESSION_ID,
        discard_changes: Annotated[bool, Field(description=(
            "Close it even if it has changes not saved, which are lost."))] = False,
    ) -> CloseResult:
        """Close a session. A session with changes not saved is not closed
        unless discard_changes is true: save it first with session_save."""
        with store.lock:
            if store.unsaved(session_id) and not discard_changes:
                raise ToolError(
                    f"'{session_id}' has changes not saved: save it with "
                    f"session_save, or pass discard_changes=true to lose them.")
            store.close(session_id)
            return CloseResult(session_id=session_id, closed=True,
                               open_sessions=store.ids)

    @server.tool()
    def gaps_summary(
        since: Annotated[Optional[str], Field(description=(
            "Only the gaps from this date on, YYYY-MM-DD."))] = None,
        kind: Annotated[Optional[Literal["declared", "prerequisite", "failed", "undone"]],
                        Field(description="Only this kind.")] = None,
        max_groups: Annotated[int, Field(ge=1, description=(
            "Maximum number of groups listed, the most frequent."))] = 50,
    ) -> GapsSummary:
        """Summarize the gaps of the collection of recipes recorded in the
        local file of the user, grouped: what analyses needed and no recipe
        provided (declared by agents, prerequisites that no recipe
        estimates, tools that failed on a valid request, recipes undone).
        The texts of the declared ones are free: join the similar groups
        yourself."""
        path = gaps.gaps_path()
        if path is None:
            return GapsSummary(enabled=False, groups=[], n_groups=0, skipped_lines=0)
        records, skipped = gaps.read_gaps(path)
        try:
            groups = gaps.summarize(records, since=since, kind=kind)
        except ValueError as e:
            raise ToolError(str(e)) from e
        return GapsSummary(enabled=True, groups=[GapGroup(**g) for g in groups[:max_groups]],
                           n_groups=len(groups), skipped_lines=skipped)
