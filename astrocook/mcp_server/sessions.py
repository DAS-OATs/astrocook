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
"""
import threading
from copy import deepcopy
from pathlib import Path
from typing import Any

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


class SessionStore:
    """Sessions by identifier.

    ``lock`` serializes the tool calls that run recipes: the SDK runs the
    calls in worker threads, and the recipes were written for one call at
    a time (e.g. module-level state, F12).
    """

    def __init__(self):
        self._histories: dict[str, SessionHistory] = {}
        self.lock = threading.Lock()
        #: Passed to the sessions as their ``gui``, so that the recipes
        #: find the other open sessions
        self.context = _Context(self._histories)

    def _new_id(self, stem: str) -> str:
        """``stem``, with a suffix (``-2``, ``-3``, ...) if it is taken."""
        session_id, n = stem, 1
        while session_id in self._histories:
            n += 1
            session_id = f"{stem}-{n}"
        return session_id

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
            known = ", ".join(sorted(self._histories)) or "none"
            raise ToolError(f"No session '{session_id}'. Open sessions: "
                            f"{known}; open one with session_open.") from None

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
            n = 1
            while f"{session_id}-{n}" in self._histories:
                n += 1
            new_id = f"{session_id}-{n}"
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
            log.add_entry(recipe_name, entry)
            new_state.name = new_id
            new_state._gui = self.context
            history.add_state(new_state)
            self._histories[new_id] = history
            ids.append(new_id)
        return ids
