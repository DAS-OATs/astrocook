"""The sessions followed by the main window (phase 6-ter, PR 5-bis; decision 7
of the phase, option b).

An archive opened with a marker next to it (:mod:`astrocook.follow`) is
*followed*: another program (the MCP server, through ``gui_open``) saves it
again after every change of the session, and the window reloads it, keeping
the view (the range of the plot, the System Inspector and its selection).

How a half-written archive is never read: the writer replaces the archive
only when it is complete (``os.replace`` of a temporary file in the same
folder, atomic on the same file system), so the path always holds a whole
archive, old or new; the window reacts to the path of the archive only
(not to the temporary files in its folder), waits for the events to stop
(``FOLLOW_DEBOUNCE_MS``), reloads only if the time and size of the file
changed since the last load, in a worker thread, and if a load fails anyway
it keeps the state it shows and tries once more (``FOLLOW_RETRY_MS``).

The folder is watched as well as the file: on Linux a file replaced is a new
file, and the watch on the old one is lost; it is added again at each event.

A followed session is read only: the recipes that change it, the continuum
editor, undo and redo, the scripts and renaming are refused, with the offer
of an editable copy (with the log of the session). The recipes that open new
sessions (branching) are allowed: their sessions are new, editable and not
followed (Guido, plan of the PR). Saving and exporting are allowed. When the
marker disappears (the agent closed the session, the server stopped) the
session is no longer followed, and stays read only.
"""
import logging
import os
import time
from copy import deepcopy
from dataclasses import dataclass
from typing import List, Optional, Tuple, cast

from PySide6.QtCore import QFileSystemWatcher, Qt, QTimer
from PySide6.QtWidgets import (
    QHBoxLayout, QLabel, QMessageBox, QPushButton, QToolBar, QWidget,
)

from astrocook import follow
from astrocook.core.session_manager import SessionHistory
from astrocook.core.structures import HistoryLogV2
from astrocook.recipes.registry import get_recipe

from .qt_workers import SessionLoadWorker
from .theme import ACCENT, FOLLOW_BANNER, themed
from astrocook.gui._window_host import HostBase

logger = logging.getLogger(__name__)

#: How long the events of the archive must stop before it is read: a writer
#: may replace it several times in a row
FOLLOW_DEBOUNCE_MS = 300
#: After a load that failed, how long before the one more try
FOLLOW_RETRY_MS = 1000


@dataclass
class _Followed:
    """A session history followed, and what is known of its archive."""
    history: SessionHistory
    #: The archive, with the links resolved
    path: str
    session_id: str
    #: (mtime_ns, size) of the archive at the last load
    stamp: Optional[Tuple[int, int]]
    #: True when the marker is gone (or its writer is dead)
    ended: bool = False
    loading: bool = False
    #: Another change arrived while loading
    again: bool = False
    retried: bool = False
    #: The worker that loads it, kept alive while it runs
    worker: object = None


def _stamp(path: str) -> Optional[Tuple[int, int]]:
    try:
        st = os.stat(path)
    except OSError:
        return None
    return (st.st_mtime_ns, st.st_size)


class MainFollowMixin(HostBase):
    """The sessions followed by ``MainWindowV2``.

    Uses the window's ``session_histories``, ``active_history``, the thread
    pool, ``_load_session``, ``_update_view_for_session``,
    ``_on_duplicate_session``, ``_show_custom_message``, ``show_status`` and
    ``_reposition_floating_widgets``.
    """

    def _setup_follow(self):
        """The watcher of the archives, its timer and the banner of a
        followed session (hidden)."""
        self._followed: List[_Followed] = []
        self._follow_watcher = QFileSystemWatcher(self)
        self._follow_watcher.fileChanged.connect(self._on_followed_path_changed)
        self._follow_watcher.directoryChanged.connect(self._on_followed_path_changed)
        self._follow_timer = QTimer(self)
        self._follow_timer.setSingleShot(True)
        self._follow_timer.setInterval(FOLLOW_DEBOUNCE_MS)
        self._follow_timer.timeout.connect(self._check_followed)

        banner = QToolBar("Followed Session", self)
        banner.setObjectName("followBanner")
        banner.setMovable(False)
        banner.setFloatable(False)
        banner.toggleViewAction().setVisible(False)  # not in the context menu
        content = QWidget()
        layout = QHBoxLayout(content)
        layout.setContentsMargins(8, 2, 8, 2)
        self.follow_banner_label = QLabel()
        # One line: the rest is in the tooltip (a label that wraps made the
        # toolbar several lines tall)
        self.follow_banner_label.setWordWrap(False)
        self.follow_banner_label.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(self.follow_banner_label, 1)
        self.follow_copy_button = QPushButton("Make Editable Copy")
        self.follow_copy_button.setToolTip(
            "A new session with the current state and log of this one, which "
            "you can change; it is not followed.")
        self.follow_copy_button.clicked.connect(lambda: self._on_make_editable_copy())
        layout.addWidget(self.follow_copy_button)
        banner.addWidget(content)
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, banner)
        banner.hide()
        self.follow_banner = banner
        self._apply_follow_style()

    def _apply_follow_style(self):
        """The colours of the banner, for the light or dark palette."""
        self.follow_banner.setStyleSheet(
            f"QToolBar#followBanner {{ background: {themed(FOLLOW_BANNER)}; border: none; "
            f"border-bottom: 2px solid {ACCENT}; padding: 2px; }}")

    def _banner_room(self) -> int:
        """The height of the banner, when it is shown (the cards float below
        it)."""
        banner = getattr(self, 'follow_banner', None)
        if banner is None or banner.isHidden():
            return 0
        return max(banner.height(), banner.sizeHint().height())

    # --- Which sessions are followed ---

    def _followed_of(self, history: Optional[SessionHistory]) -> Optional[_Followed]:
        for followed in getattr(self, '_followed', []):
            if followed.history is history:
                return followed
        return None

    def is_followed(self, history: Optional[SessionHistory] = None) -> bool:
        """Whether a session history (the active one by default) is
        followed, and so read only (also when it is no longer updated)."""
        return self._followed_of(history or self.active_history) is not None

    def _follow_if_marked(self, history: SessionHistory, path: Optional[str]):
        """Follow ``history`` if the archive it was opened from has a
        marker."""
        if not path or self._followed_of(history) is not None:
            return
        marker = follow.read_marker(path)
        if marker is None:
            return
        real = os.path.realpath(path)
        followed = _Followed(history, real, marker.session_id, _stamp(real),
                             ended=not marker.writer_alive())
        self._followed.append(followed)
        if followed.ended:
            logger.info(f"{real} was followed, but the program that wrote it has "
                        "stopped: it is read only and no longer updated.")
        else:
            self._watch(followed)
            logger.info(f"Following {real} (session '{marker.session_id}').")
        self._update_follow_banner()
        self._update_undo_redo_actions()

    def _show_if_followed(self, path: str) -> bool:
        """If ``path`` is the archive of a session already followed (and still
        updated), show that session and return True: ``gui_open`` called
        again sends the same archive to the application already open, on
        macOS."""
        real = os.path.realpath(path)
        for followed in self._followed:
            if followed.path == real and not followed.ended:
                index = self.session_histories.index(followed.history)
                self.active_history = followed.history
                self._update_view_for_session(followed.history.current_state,
                                              set_current_list_item=True,
                                              target_list_index=index)
                self._update_undo_redo_actions()
                self.activateWindow()
                self.raise_()
                return True
        return False

    def _watch(self, followed: _Followed):
        paths = [p for p in (followed.path, os.path.dirname(followed.path))
                 if os.path.exists(p) and p not in self._follow_watcher.files()
                 and p not in self._follow_watcher.directories()]
        if paths:
            self._follow_watcher.addPaths(paths)

    def _unwatch(self, followed: _Followed):
        others = [f for f in self._followed if f is not followed and not f.ended]
        for path in (followed.path, os.path.dirname(followed.path)):
            if any(path in (f.path, os.path.dirname(f.path)) for f in others):
                continue
            if path in self._follow_watcher.files() or path in self._follow_watcher.directories():
                self._follow_watcher.removePath(path)

    def _forget_followed(self, history: SessionHistory):
        """A followed session closed in the window."""
        followed = self._followed_of(history)
        if followed is not None:
            self._unwatch(followed)
            self._followed.remove(followed)

    # --- Reloading ---

    def _on_followed_path_changed(self, path: str):
        """An event of a watched archive or folder: watch the archives again
        (a replaced file is a new file), and check them when the events
        stop."""
        for followed in self._followed:
            if not followed.ended:
                self._watch(followed)
        self._follow_timer.start()  # restarted by each event

    def _check_followed(self):
        """Reload the followed archives that changed; end those whose marker
        is gone."""
        for followed in list(self._followed):
            if followed.ended:
                continue
            if followed.history not in self.session_histories:
                self._forget_followed(followed.history)
                continue
            marker = follow.read_marker(followed.path)
            if marker is None or marker.session_id != followed.session_id:
                self._end_follow(followed)
                continue
            stamp = _stamp(followed.path)
            if stamp is not None and stamp != followed.stamp:
                self._reload_followed(followed, stamp)

    def _reload_followed(self, followed: _Followed, stamp):
        if followed.loading:
            followed.again = True
            return
        followed.loading = True
        worker = SessionLoadWorker([followed.path], self._load_session)
        worker.signals.loaded.connect(
            lambda i, path, session: self._on_followed_loaded(followed, stamp, session))
        worker.signals.failed.connect(
            lambda i, path, message: self._on_followed_failed(followed, message))
        followed.worker = worker
        self.thread_pool.start(worker)

    def _on_followed_loaded(self, followed: _Followed, stamp, session):
        followed.loading = False
        followed.retried = False
        if followed not in self._followed or followed.history not in self.session_histories:
            return  # closed while it loaded
        followed.stamp = stamp
        self._replace_followed_state(followed.history, session)
        if followed.again:
            followed.again = False
            self._check_followed()

    def _on_followed_failed(self, followed: _Followed, message: str):
        """A load that failed: the state shown is kept; one more try."""
        followed.loading = False
        if followed not in self._followed:
            return
        if followed.again or not followed.retried:
            followed.again = False
            followed.retried = True
            QTimer.singleShot(FOLLOW_RETRY_MS, self._check_followed)
        else:
            logger.warning(f"The followed session could not be reloaded from "
                           f"{followed.path} ({message}): it shows the last state read.")

    def _replace_followed_state(self, history: SessionHistory, session):
        """The state read from the archive replaces the one shown (a followed
        session has one state: its history is in its log), with the log of
        the archive; the view is kept."""
        # The log of the agent, read only (F109: the state read is the start)
        log = HistoryLogV2.opened_from(session.log_manager)
        session.__dict__.pop('_source_path', None)
        session._gui = self
        session.name = history.current_state.name
        history.states = [session]
        history.current_index = 0
        history.log_manager = log
        session.log_manager = log
        if history is self.active_history:
            # Not autoscaled: a zoomed view is kept, and the System Inspector
            # keeps its selection (by the identifiers of the components)
            self._update_view_for_session(session, set_current_list_item=False,
                                          force_autoscale=False)
            self._update_undo_redo_actions()
        self.show_status(f"'{history.display_name}' updated by the agent at "
                         f"{time.strftime('%H:%M:%S')}.")

    def _end_follow(self, followed: _Followed):
        followed.ended = True
        self._unwatch(followed)
        logger.info(f"{followed.path} is no longer followed.")
        self._update_follow_banner()
        self.show_status(f"'{followed.history.display_name}' is no longer followed.")

    # --- The banner and the read-only guard ---

    def _update_follow_banner(self):
        """Show the banner when the active session is followed, with what it
        says; hide it otherwise."""
        banner = getattr(self, 'follow_banner', None)
        if banner is None:
            return
        followed = self._followed_of(self.active_history)
        was_visible = not banner.isHidden()
        if followed is None:
            banner.hide()
        else:
            if followed.ended:
                text = ("<b>No longer followed.</b> The agent closed this session, "
                        "or its server stopped. Read only.")
                tip = ("This session was changed by an AI agent through the MCP server; "
                       "it is not updated any more.")
            else:
                text = ("<b>Followed session.</b> Updated by an AI agent after every "
                        "change. Read only.")
                tip = ("An AI agent changes this session through the MCP server, and the "
                       "window shows each change. Nothing done here reaches the agent: "
                       "make an editable copy to work on it.")
            self.follow_banner_label.setText(text)
            self.follow_banner_label.setToolTip(tip)
            banner.show()
        if was_visible != (not banner.isHidden()):
            QTimer.singleShot(0, self._reposition_floating_widgets)

    def _refuse_if_followed(self, action: str) -> bool:
        """If the active session is followed, say that it is read only and
        offer an editable copy, and return True: the caller does not go on."""
        history = cast(SessionHistory, self.active_history)
        if self._followed_of(history) is None:
            return False
        ret = self._show_custom_message(
            title="Read-Only Session",
            header=f"'{history.display_name}' is read only.",
            text=(f"{action} would change this session, which is followed: an AI agent "
                  "changes it through the MCP server, and the window shows its "
                  "changes. Make an editable copy to work on it."),
            buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            default_btn=QMessageBox.StandardButton.Yes,
            button_texts={QMessageBox.StandardButton.Yes: "Make Editable Copy"})
        if ret == QMessageBox.StandardButton.Yes:
            self._on_make_editable_copy(history)
        return True

    def _refuse_recipe_if_followed(self, category: str, name: str) -> bool:
        """:meth:`_refuse_if_followed` for a recipe that changes the session;
        the branching ones (a new session) are allowed, and the export (it
        writes files, and changes nothing)."""
        if name == "export_ascii":
            return False
        try:
            spec = get_recipe(category, name)
        except KeyError:
            return self._refuse_if_followed("This operation")
        if spec.branching:
            return False
        return self._refuse_if_followed(spec.brief.rstrip('.'))

    def _on_make_editable_copy(self, history: Optional[SessionHistory] = None):
        """A new session with the current state of a followed one and a copy
        of its log (what the agent did, with its reasons and notes), not
        followed."""
        history = history or self.active_history
        if history is None:
            return
        log = deepcopy(history.log_manager)
        self._on_duplicate_session(history)
        copy = cast(SessionHistory, self.active_history)
        if copy is history:
            return
        copy.log_manager = log
        copy.current_state.log_manager = log
        if self.log_scripter_dialog and self.log_scripter_dialog.isVisible():
            self.log_scripter_dialog.set_log_object(log)
        self.show_status(f"'{copy.display_name}': an editable copy of "
                         f"'{history.display_name}'.")
