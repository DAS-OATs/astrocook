"""The progress of a running recipe, as the recipes report it in the log.

A long recipe writes ``>> PROGRESS: n`` (n from 0 to 100) with
``logger.info``, next to the lines that say what it is doing. The reader of
those lines is whoever shows the progress: the GUI (the dialog with its
bar, ``gui/recipe_progress.py``), the MCP server (progress notifications to
the client, phase 6-ter PR 4) or a script (:func:`watch_progress`). The
format is read here only, so that the readers cannot differ.
"""
import logging
import threading
from contextlib import contextmanager
from typing import Callable, Iterator

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
