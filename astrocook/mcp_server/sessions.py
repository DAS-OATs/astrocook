"""The sessions open in the server.

Each session is a :class:`~astrocook.core.session_manager.SessionHistory`,
as in the GUI, under an identifier that the tools return and take. The log
entries are added here, after each recipe, as the GUI does: the recipes
do not write them (``ROADMAP.md``, section 4.1).

The recipes that read another open session by its name (a parameter of
kind ``session`` in the registry, e.g. ``reference_session`` of
``equalize``) look for it in ``session._gui.session_histories``, by
``display_name``: the GUI is that context in the application. Here the
store is (:attr:`SessionStore.context`), and the name of a session is its
identifier (phase 6-ter, PR 1).

Each state of a history is the result of one log entry (the first state
is the session as opened, after the entries of an archive's own log): the
store keeps where the entries of this server begin (``first_entry``), so
that the history and the log can be read together (phase 6-ter, PR 2). It
also keeps the state last saved (or opened), so that a session with
changes not saved is not closed without being asked.
"""
import threading
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from mcp.server.mcpserver.exceptions import ToolError

from astrocook.core.session import SessionV2
from astrocook.core.session_manager import SessionHistory
from astrocook.core.structures import HistoryLogV2
from astrocook.recipes.params import resolve_params


class _Context:
    """What the recipes look for in ``session._gui``: the open sessions."""

    def __init__(self, histories: dict):
        self._histories = histories

    @property
    def session_histories(self) -> list:
        return list(self._histories.values())


@dataclass
class SessionMeta:
    """What the store knows of a session besides its history."""
    #: The file it was opened from, or None for a session made in the
    #: server (a branch, a copy)
    path: Optional[str] = None
    #: The session it was made from (a branch, a copy), or None
    parent: Optional[str] = None
    #: How it was made: 'opened', 'branched' (by the recipe ``made_by``),
    #: 'duplicated'
    origin: str = "opened"
    made_by: Optional[str] = None
    #: The index of the log entry of the second state (the entries before it
    #: came with an archive, or from the parent)
    first_entry: int = 0
    #: The state last saved or opened; None if never saved
    saved_state: Optional[SessionV2] = None
    #: Where it was last saved
    saved_path: Optional[str] = None


class SessionStore:
    """Sessions by identifier.

    ``lock`` serializes the tool calls that run recipes: the SDK runs the
    calls in worker threads, and the recipes were written for one call at
    a time (e.g. module-level state, F12).
    """

    def __init__(self):
        self._histories: dict[str, SessionHistory] = {}
        self.meta: dict[str, SessionMeta] = {}
        self.lock = threading.Lock()
        #: Passed to the sessions as their ``gui``, so that the recipes
        #: find the other open sessions
        self.context = _Context(self._histories)

    @property
    def ids(self) -> list[str]:
        """The identifiers of the open sessions, in the order they were
        opened."""
        return list(self._histories)

    def _new_id(self, stem: str) -> str:
        """``stem``, with a suffix (``-2``, ``-3``, ...) if it is taken."""
        session_id, n = stem, 1
        while session_id in self._histories:
            n += 1
            session_id = f"{stem}-{n}"
        return session_id

    def _child_id(self, parent: str) -> str:
        """``<parent>-1``, ``<parent>-2``, ...: the first one free."""
        n = 1
        while f"{parent}-{n}" in self._histories:
            n += 1
        return f"{parent}-{n}"

    def add(self, session: SessionV2, path: str) -> str:
        """Store a newly opened session and return its identifier.

        The identifier is the file name without extension, with a suffix
        (``-2``, ``-3``, ...) if that is already taken. A V2 log restored
        from an ``.acs2`` archive is kept, as in the GUI.
        """
        session_id = self._new_id(Path(path).stem or "session")
        log = session.log_manager
        if not isinstance(log, HistoryLogV2):
            log = HistoryLogV2()
        session._gui = self.context
        history = SessionHistory(session, log)
        history.display_name = session_id
        self._histories[session_id] = history
        self.meta[session_id] = SessionMeta(
            path=str(path), first_entry=log.current_index + 1,
            saved_state=session, saved_path=str(path))
        return session_id

    def get(self, session_id: str) -> SessionHistory:
        """Return the history of a session.

        Raises
        ------
        ToolError
            If no session has this identifier.
        """
        try:
            return self._histories[session_id]
        except KeyError:
            known = ", ".join(self._histories) or "none"
            raise ToolError(f"No session '{session_id}'. Open sessions: "
                            f"{known}; open one with session_open.") from None

    def unsaved(self, session_id: str) -> bool:
        """Whether the current state differs from the one last saved (or
        opened)."""
        return self.get(session_id).current_state is not self.meta[session_id].saved_state

    def mark_saved(self, session_id: str, path: str) -> None:
        meta = self.meta[session_id]
        meta.saved_state = self.get(session_id).current_state
        meta.saved_path = str(path)

    def entry_of(self, session_id: str, state_index: int):
        """The log entry that made the state ``state_index`` (a dict with
        ``recipe_name`` and ``params``), or None for the first state."""
        if state_index == 0:
            return None
        log = self.get(session_id).log_manager
        i = self.meta[session_id].first_entry + state_index - 1
        if 0 <= i < len(log.entries):
            entry = log.entries[i]
            return {"recipe_name": entry.recipe_name, "params": dict(entry.params)}
        return None

    def _add_provided(self, history: SessionHistory, provided) -> None:
        """The recipes run first to provide what the session lacked
        (:class:`~astrocook.recipes.prerequisites.ProvidedStep`), each with
        its state and log entry, as in the GUI (phase 6-bis, PR 8): the
        replay of the log runs them."""
        for step in provided:
            history.log_manager.add_entry(step.recipe.split('.')[1], dict(step.params))
            step.session.log_manager = history.log_manager
            history.add_state(step.session)

    def apply(self, session_id: str, new_state: SessionV2, recipe_name: str,
              params: dict[str, Any], provided=()) -> SessionV2:
        """Record a recipe in the log and make its result the current state."""
        history = self.get(session_id)
        self._add_provided(history, provided)
        # The values the recipe used, defaults included (#573, item 3)
        history.log_manager.add_entry(recipe_name, resolve_params(recipe_name, params))
        new_state.log_manager = history.log_manager
        history.add_state(new_state)
        return new_state

    def _copy(self, session_id: str, new_id: str) -> SessionHistory:
        """A new history with a copy of the states of ``session_id`` up to the
        current one, and of its log, as the GUI does for a branch
        (``update_gui_session_state``)."""
        parent = self.get(session_id)
        log = deepcopy(parent.log_manager)
        states = []
        for state in parent.states[:parent.current_index + 1]:
            copy = state.with_new_spectrum(state.spec).with_new_system_list(state.systs)
            copy.name = new_id
            copy._gui = self.context
            states.append(copy)
        history = SessionHistory(states[0], log)
        for state in states[1:]:
            history.add_state(state)
        history.display_name = new_id
        return history

    def branch(self, session_id: str, new_states, recipe_name: str,
               params: dict[str, Any], provided=()) -> list[str]:
        """Open the sessions made by a branching recipe and return their
        identifiers (``<parent>-1``, ``<parent>-2``, ...).

        As in the GUI (``update_gui_session_state``): the steps provided
        first are states of the parent; each new session starts as a copy
        of the parent's history and log, then the recipe is logged in it,
        not in the parent, which does not change (F41, #575).
        """
        parent = self.get(session_id)
        self._add_provided(parent, provided)
        entry = resolve_params(recipe_name, params)
        ids = []
        for new_state in new_states:
            new_id = self._child_id(session_id)
            history = self._copy(session_id, new_id)
            history.log_manager.add_entry(recipe_name, entry)
            new_state.name = new_id
            new_state._gui = self.context
            history.add_state(new_state)
            self._histories[new_id] = history
            self.meta[new_id] = SessionMeta(
                parent=session_id, origin="branched", made_by=recipe_name,
                first_entry=self.meta[session_id].first_entry)
            ids.append(new_id)
        return ids

    def duplicate(self, session_id: str) -> str:
        """A copy of a session, with its history up to the current state and
        its log, under a new identifier (``<id>-1``, ...)."""
        new_id = self._child_id(session_id)
        self._histories[new_id] = self._copy(session_id, new_id)
        self.meta[new_id] = SessionMeta(parent=session_id, origin="duplicated",
                                        first_entry=self.meta[session_id].first_entry)
        return new_id

    def close(self, session_id: str) -> None:
        self.get(session_id)
        del self._histories[session_id]
        del self.meta[session_id]
