"""The sessions open in the server.

Each session is a :class:`~astrocook.core.session_manager.SessionHistory`,
as in the GUI, under an identifier that the tools return and take. The log
entries are added here, after each recipe, as the GUI does: the recipes
do not write them (``ROADMAP.md``, section 4.1).
"""
import threading
from pathlib import Path
from typing import Any

from mcp.server.mcpserver.exceptions import ToolError

from astrocook.core.session import SessionV2
from astrocook.core.session_manager import SessionHistory
from astrocook.core.structures import HistoryLogV2


class SessionStore:
    """Sessions by identifier.

    ``lock`` serializes the tool calls that run recipes: the SDK runs the
    calls in worker threads, and the recipes were written for one call at
    a time (e.g. module-level state, F12).
    """

    def __init__(self):
        self._histories: dict[str, SessionHistory] = {}
        self.lock = threading.Lock()

    def add(self, session: SessionV2, path: str) -> str:
        """Store a newly opened session and return its identifier.

        The identifier is the file name without extension, with a suffix
        (``-2``, ``-3``, ...) if that is already taken. A V2 log restored
        from an ``.acs2`` archive is kept, as in the GUI.
        """
        stem = Path(path).stem or "session"
        session_id, n = stem, 1
        while session_id in self._histories:
            n += 1
            session_id = f"{stem}-{n}"
        log = session.log_manager
        if not isinstance(log, HistoryLogV2):
            log = HistoryLogV2()
        self._histories[session_id] = SessionHistory(session, log)
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

    def apply(self, session_id: str, new_state: SessionV2, recipe_name: str,
              params: dict[str, Any]) -> SessionV2:
        """Record a recipe in the log and make its result the current state."""
        history = self.get(session_id)
        history.log_manager.add_entry(recipe_name, params)
        new_state.log_manager = history.log_manager
        history.add_state(new_state)
        return new_state
