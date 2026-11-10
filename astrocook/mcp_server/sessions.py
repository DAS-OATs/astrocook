"""The sessions open in the server.

Each session is a :class:`~astrocook.core.session_manager.SessionHistory`,
as in the GUI, under an identifier that the tools return and take. The log
entries are added here, after each recipe, as the GUI does: the recipes
do not write them (``ROADMAP_HISTORY.md``, section 4.1).

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

A session can be followed by the GUI (``gui_open``, phase 6-ter PR 5-bis):
the store then saves it to an archive of its own folder after every change of
its state or log (a recipe, the steps provided before it, undo and redo, a
note), with the marker of :mod:`astrocook.follow`, and the GUI reloads it.
The folder is removed when the store is collected or the server exits.
"""
import logging
import shutil
import tarfile
import tempfile
import threading
import time
import weakref
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Optional

from mcp.server.mcpserver.exceptions import ToolError

from astrocook import follow, gaps
from astrocook.core.session import SessionV2
from astrocook.core.session_manager import SessionHistory
from astrocook.core.structures import HistoryLogV2
from astrocook.recipes.params import resolve_params

logger = logging.getLogger(__name__)


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
    #: The index of the log entry of the second state; 0 since the log of an
    #: archive opened is no longer continued (F109)
    first_entry: int = 0
    #: The state last saved or opened; None if never saved
    saved_state: Optional[SessionV2] = None
    #: Where it was last saved
    saved_path: Optional[str] = None
    #: The log as last saved or opened (``HistoryLogV2.to_dict``): a note
    #: changes the log, not the state (phase 6-ter, PR 5a)
    saved_log: Optional[dict] = None
    #: The indices of the log entries of the steps provided before a
    #: recipe (a continuum made for it), not run by the agent (PR 5b)
    provided: set = field(default_factory=set)
    #: The gaps seen by the server in this session, recorded once each
    seen_gaps: set = field(default_factory=set)
    #: The archive followed by the GUI (``gui_open``), saved after every
    #: change; None if the session is not followed
    gui_path: Optional[Path] = None
    #: Why the archive followed by the GUI could not be saved, for the next
    #: output of a tool on the session
    gui_warnings: list = field(default_factory=list)


class SessionStore:
    """Sessions by identifier.

    Two locks (phase 6-ter, F113). ``recipe_lock`` serializes the tool calls
    that run recipes or change the states: the SDK runs the calls in worker
    threads, and the recipes were written for one call at a time (e.g.
    module-level state, F12). ``state_lock`` guards the histories and the
    records of the store, and is held only briefly: the tools that read a
    session (and save it, or write a note) take only this one, so that they
    answer while a recipe runs. They see the last state made: the states
    are immutable, and a recipe adds its own only when it ends. A tool that
    changes the states takes both (:meth:`writing`), always in this order.
    """

    def __init__(self):
        self._histories: dict[str, SessionHistory] = {}
        self.meta: dict[str, SessionMeta] = {}
        self.recipe_lock = threading.Lock()
        self.state_lock = threading.RLock()
        #: The recipe running on a session, and when it started (monotonic)
        self.running: dict[str, tuple[str, float]] = {}
        #: Passed to the sessions as their ``gui``, so that the recipes
        #: find the other open sessions
        self.context = _Context(self._histories)
        #: The folder of the archives followed by the GUI, made when the
        #: first one is
        self._gui_folder: Optional[Path] = None

    @contextmanager
    def writing(self) -> Iterator[None]:
        """Both locks: for a tool that changes the states of a session."""
        with self.recipe_lock, self.state_lock:
            yield

    def running_text(self, session_id: str) -> Optional[str]:
        """The recipe running on a session and for how long, or None."""
        if session_id not in self.running:
            return None
        recipe, start = self.running[session_id]
        return f"{recipe}, for {time.monotonic() - start:.0f} s"

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
        (``-2``, ``-3``, ...) if that is already taken. The log starts
        empty, as in the GUI: the log of an ``.acs2`` archive is its
        ``previous``, read only (``log_before`` of ``session_describe``),
        never replayed nor saved again (F109).
        """
        session_id = self._new_id(Path(path).stem or "session")
        log = HistoryLogV2.opened_from(session.log_manager)
        session._gui = self.context
        history = SessionHistory(session, log)
        history.display_name = session_id
        self._histories[session_id] = history
        self.meta[session_id] = SessionMeta(
            path=str(path), first_entry=log.current_index + 1,
            saved_state=session, saved_path=str(path), saved_log=log.to_dict())
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
        """Whether the current state, or the log (a note), differs from the
        one last saved (or opened)."""
        history = self.get(session_id)
        meta = self.meta[session_id]
        return (history.current_state is not meta.saved_state
                or history.log_manager.to_dict() != meta.saved_log)

    def mark_saved(self, session_id: str, path: str) -> None:
        meta = self.meta[session_id]
        history = self.get(session_id)
        meta.saved_state = history.current_state
        meta.saved_log = history.log_manager.to_dict()
        meta.saved_path = str(path)

    def entry_of(self, session_id: str, state_index: int):
        """The log entry that made the state ``state_index`` (a dict with
        ``recipe_name``, ``params`` and ``reason``), or None for the first
        state."""
        if state_index == 0:
            return None
        log = self.get(session_id).log_manager
        i = self.meta[session_id].first_entry + state_index - 1
        if 0 <= i < len(log.entries):
            entry = log.entries[i]
            return {"recipe_name": entry.recipe_name, "params": dict(entry.params),
                    "reason": entry.reason}
        return None

    def notes_of(self, session_id: str, state_index: int) -> list:
        """The notes written at the state ``state_index`` (``LogNoteV2``);
        those of a state undone are kept, as the state is, for a redo."""
        log = self.get(session_id).log_manager
        notes = log.notes_at(self.meta[session_id].first_entry + state_index)
        if state_index == 0 and log.previous is not None:
            # Those written at the end of the log of the archive (F109)
            notes = log.previous.notes_at(len(log.previous.entries)) + notes
        return notes

    def note(self, session_id: str, text: str, kind: str = "note"):
        """Add a note at the current state of a session (no new state). A
        gap (``kind="gap"``, declared by the agent) is also recorded in the
        file of the gaps of the user."""
        note = self.get(session_id).log_manager.add_note(text, kind)
        if kind == "gap":
            gaps.record_gap("declared", note.text, session_id)
        self.refresh_gui(session_id)
        return note

    def record_gap(self, session_id: str, kind: str, text: str,
                   once: bool = False) -> None:
        """Record a gap seen by the server (``kind`` of :data:`astrocook.gaps.KINDS`):
        a note of the session log, kept with the session, and a line of the
        file of the gaps. With ``once``, a gap of the same kind and text is
        recorded once per session. The note does not make the session
        unsaved: nothing the agent did is lost if it is closed."""
        meta = self.meta[session_id]
        if once:
            if (kind, text) in meta.seen_gaps:
                return
            meta.seen_gaps.add((kind, text))
        was_unsaved = self.unsaved(session_id)
        self.get(session_id).log_manager.add_note(text, "gap_seen")
        if not was_unsaved:
            meta.saved_log = self.get(session_id).log_manager.to_dict()
        gaps.record_gap(kind, text, session_id)
        self.refresh_gui(session_id)

    def was_provided(self, session_id: str, state_index: int) -> bool:
        """Whether the state ``state_index`` is a step provided before a
        recipe, not run by the agent."""
        return (self.meta[session_id].first_entry + state_index - 1
                in self.meta[session_id].provided)

    def _prune_provided(self, session_id: str) -> None:
        """Forget the provided steps whose entries a new entry drops (those
        after the current one)."""
        meta, log = self.meta[session_id], self.get(session_id).log_manager
        meta.provided = {i for i in meta.provided if i <= log.current_index}

    def _add_provided(self, history: SessionHistory, provided) -> None:
        """The recipes run first to provide what the session lacked
        (:class:`~astrocook.recipes.prerequisites.ProvidedStep`), each with
        its state and log entry, as in the GUI (phase 6-bis, PR 8): the
        replay of the log runs them."""
        meta = self.meta.get(history.display_name)
        for step in provided:
            log = history.log_manager
            if meta is not None:
                # Entries after the current one are dropped by the new one
                meta.provided = {i for i in meta.provided if i <= log.current_index}
            log.add_entry(step.recipe.split('.')[1], dict(step.params))
            if meta is not None:
                meta.provided.add(log.current_index)
            step.session.log_manager = history.log_manager
            history.add_state(step.session)

    def apply(self, session_id: str, new_state: SessionV2, recipe_name: str,
              params: dict[str, Any], provided=(), reason: Optional[str] = None) -> SessionV2:
        """Record a recipe in the log, with the reason given for it, and make
        its result the current state."""
        history = self.get(session_id)
        self._add_provided(history, provided)
        self._prune_provided(session_id)
        # The values the recipe used, defaults included (#573, item 3)
        history.log_manager.add_entry(recipe_name, resolve_params(recipe_name, params),
                                      reason=reason)
        new_state.log_manager = history.log_manager
        history.add_state(new_state)
        self.refresh_gui(session_id)
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
               params: dict[str, Any], provided=(),
               reason: Optional[str] = None) -> list[str]:
        """Open the sessions made by a branching recipe and return their
        identifiers (``<parent>-1``, ``<parent>-2``, ...).

        As in the GUI (``update_gui_session_state``): the steps provided
        first are states of the parent; each new session starts as a copy
        of the parent's history and log, then the recipe is logged in it,
        not in the parent, which does not change (F41, #575).
        """
        parent = self.get(session_id)
        self._add_provided(parent, provided)
        self._prune_provided(session_id)
        entry = resolve_params(recipe_name, params)
        ids = []
        for new_state in new_states:
            new_id = self._child_id(session_id)
            history = self._copy(session_id, new_id)
            history.log_manager.add_entry(recipe_name, entry, reason=reason)
            new_state.name = new_id
            new_state._gui = self.context
            history.add_state(new_state)
            self._histories[new_id] = history
            self.meta[new_id] = SessionMeta(
                parent=session_id, origin="branched", made_by=recipe_name,
                first_entry=self.meta[session_id].first_entry,
                provided=set(self.meta[session_id].provided))
            ids.append(new_id)
        if provided:
            # The steps provided are states of the parent; the new sessions
            # are not followed
            self.refresh_gui(session_id)
        return ids

    def duplicate(self, session_id: str) -> str:
        """A copy of a session, with its history up to the current state and
        its log, under a new identifier (``<id>-1``, ...)."""
        new_id = self._child_id(session_id)
        self._histories[new_id] = self._copy(session_id, new_id)
        self.meta[new_id] = SessionMeta(parent=session_id, origin="duplicated",
                                        first_entry=self.meta[session_id].first_entry,
                                        provided=set(self.meta[session_id].provided))
        return new_id

    def close(self, session_id: str) -> None:
        self.get(session_id)
        self.unfollow(session_id)
        del self._histories[session_id]
        del self.meta[session_id]

    # --- The GUI that follows a session (phase 6-ter, PR 5-bis) ---

    def gui_folder(self) -> Path:
        """The folder of the archives followed by the GUI (a temporary one,
        removed when the store is collected or the server exits)."""
        if self._gui_folder is None:
            self._gui_folder = Path(tempfile.mkdtemp(prefix="astrocook-gui-"))
            weakref.finalize(self, shutil.rmtree, self._gui_folder, True)
        return self._gui_folder

    def follow(self, session_id: str) -> Path:
        """Save a session to its archive for the GUI, with the marker that
        says it is followed, and keep it saved after every change; return
        the path of the archive.

        Raises
        ------
        OSError, ValueError
            If the archive cannot be written (it is then not followed).
        """
        history = self.get(session_id)
        meta = self.meta[session_id]
        path = meta.gui_path or self.gui_folder() / f"{session_id}.acs2"
        follow.write_archive(history, path)
        follow.write_marker(path, session_id)
        meta.gui_path = path
        return path

    def unfollow(self, session_id: str) -> None:
        """Stop following a session: its marker and archive are removed (the
        GUI says that the session is no longer followed)."""
        meta = self.meta[session_id]
        if meta.gui_path is not None:
            follow.remove(meta.gui_path)
            meta.gui_path = None

    def refresh_gui(self, session_id: str) -> None:
        """Save the archive followed by the GUI again, if the session is
        followed. A failure does not stop the tool that changed the session:
        it becomes a warning of its output, and the next change saves the
        archive again."""
        meta = self.meta.get(session_id)
        if meta is None or meta.gui_path is None:
            return
        try:
            follow.write_archive(self.get(session_id), meta.gui_path)
        except (OSError, ValueError, TypeError, KeyError, tarfile.TarError) as e:
            # The GUI is a view of the session: what stops its archive (a disk
            # full, a file held open on Windows, a table that cannot be
            # written) must not lose the result of the recipe
            message = (f"The archive followed by the GUI was not updated ({e}); "
                       "it is saved again after the next change.")
            logger.warning(message)
            meta.gui_warnings.append(message)

    def pop_gui_warnings(self, session_id: str) -> list:
        """The warnings of :meth:`refresh_gui` not yet returned by a tool."""
        meta = self.meta.get(session_id)
        if meta is None:
            return []
        warnings, meta.gui_warnings = meta.gui_warnings, []
        return warnings
