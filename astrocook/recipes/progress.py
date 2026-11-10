"""The progress of a running recipe, as the recipes report it in the log.

A long recipe writes ``>> PROGRESS: n`` (n from 0 to 100) with
``logger.info``, next to the lines that say what it is doing. The reader of
those lines is whoever shows the progress: the GUI (the dialog with its
bar, ``gui/recipe_progress.py``), the MCP server (progress notifications to
the client, phase 6-ter PR 4) or a script (:func:`watch_progress`). The
format is read here only, so that the readers cannot differ.

The same places are where a recipe can be cancelled (phase 6-ter, PR
7-quater): whoever runs a recipe in a thread runs it inside
:func:`cancellable`, with a ``threading.Event`` that it sets from another
thread to cancel it (the Stop of the GUI, the cancellation of a client of
the MCP server, its ``session_cancel`` tool); the recipe calls
:func:`check_cancelled`, or :func:`progress`, which also logs the progress
line, and stops with :class:`~astrocook.recipes.errors.RecipeCancelled`.
The event belongs to the thread, not to the session, so that the sessions
stay immutable and the steps of a pipeline (each on a new session) see it.
"""
import logging
import threading
from contextlib import contextmanager
from typing import Callable, Iterator, Optional

from astrocook.recipes.errors import RecipeCancelled

#: What a recipe writes before the percentage
PROGRESS_MARKER = ">> PROGRESS:"


def is_progress_line(message: str) -> bool:
    """Whether a log message is a progress line (its value may be unreadable)."""
    return PROGRESS_MARKER in message


def parse_progress(message: str) -> int | None:
    """The percentage of a progress line, or None if the message is not
    one or its value is not an integer."""
    if not is_progress_line(message):
        return None
    try:
        return int(message.split(PROGRESS_MARKER)[1].strip())
    except ValueError:
        return None


class _ProgressHandler(logging.Handler):
    """Calls ``callback(percent, status)`` for each progress line logged by
    one thread (a synchronous call of the server runs in a worker thread,
    possibly next to others); ``status`` is the last other line logged."""

    def __init__(self, thread_id: int, callback: Callable[[int, str], None]):
        super().__init__(level=logging.INFO)
        self.thread_id = thread_id
        self.callback = callback
        self.status = ""

    def emit(self, record: logging.LogRecord) -> None:
        if record.thread != self.thread_id:
            return
        message = record.getMessage()
        if not is_progress_line(message):
            self.status = message.strip()
            return
        percent = parse_progress(message)
        if percent is not None:
            self.callback(percent, self.status)


@contextmanager
def watch_progress(callback: Callable[[int, str], None],
                   logger_name: str = "astrocook") -> Iterator[None]:
    """Call ``callback(percent, status)`` for each progress line that the
    current thread logs, while the block runs.

    The records are INFO: the level of the logger ``logger_name`` (the
    package, by default) is lowered to INFO if it is higher, and restored
    after, so that the lines reach the handler whatever the application
    configured for the root logger.
    """
    handler = _ProgressHandler(threading.get_ident(), callback)
    logger = logging.getLogger(logger_name)
    level = logger.level
    if logger.getEffectiveLevel() > logging.INFO:
        logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    try:
        yield
    finally:
        logger.removeHandler(handler)
        logger.setLevel(level)


# --- Cancellation (phase 6-ter, PR 7-quater) ------------------------------

_current = threading.local()


@contextmanager
def cancellable(event: Optional[threading.Event] = None) -> Iterator[threading.Event]:
    """Run the block (a recipe) so that setting ``event`` cancels it.

    The event (a new one if None) is the one that :func:`check_cancelled`
    reads in the current thread while the block runs; it is yielded, to be
    set from another thread.
    """
    event = event if event is not None else threading.Event()
    previous = getattr(_current, "event", None)
    _current.event = event
    try:
        yield event
    finally:
        _current.event = previous


def cancel_requested() -> bool:
    """Whether the recipe running in the current thread was cancelled."""
    event = getattr(_current, "event", None)
    return event is not None and event.is_set()


def check_cancelled() -> None:
    """Stop the recipe running in the current thread if it was cancelled.

    Raises
    ------
    RecipeCancelled
        If the event of :func:`cancellable` is set.
    """
    if cancel_requested():
        raise RecipeCancelled("cancelled while it ran: the session is as before it")


def progress(logger: logging.Logger, percent: int) -> None:
    """Log the progress line of a recipe, then stop it if it was cancelled."""
    logger.info(f"{PROGRESS_MARKER} {percent}")
    check_cancelled()
