"""Run a recipe and turn its log records into MCP results.

The recipes report a failure by logging an error and returning 0, and a
problem that does not stop them by logging a warning. The model sees
neither: it sees only what a tool returns. :func:`run_recipe` collects the
records logged while a recipe runs; a failure becomes a
:class:`~mcp.server.mcpserver.exceptions.ToolError` carrying the logged
messages, the warnings are returned to be put in the tool output.

An unexpected exception becomes a ``ToolError`` too: the SDK would
otherwise send the model only "Error executing tool <name>", without the
message (``mcp`` 2.3).
"""
import logging
import threading
from contextlib import contextmanager
from typing import Any, Callable, Iterator

from mcp.server.mcpserver.exceptions import ToolError


class _Collector(logging.Handler):
    """Keep the records of one thread (the SDK runs each call of a
    synchronous tool in a worker thread, possibly next to other calls)."""

    def __init__(self, thread_id: int):
        super().__init__(level=logging.WARNING)
        self.thread_id = thread_id
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        if record.thread == self.thread_id:
            self.records.append(record)


@contextmanager
def collect_log() -> Iterator[list[logging.LogRecord]]:
    """Collect the warnings and errors logged by the current thread."""
    collector = _Collector(threading.get_ident())
    root = logging.getLogger()
    root.addHandler(collector)
    try:
        yield collector.records
    finally:
        root.removeHandler(collector)


def _messages(records: list[logging.LogRecord], min_level: int,
              max_level: int = logging.CRITICAL) -> list[str]:
    """Messages of the records within the levels, without repetitions."""
    out: list[str] = []
    for r in records:
        msg = r.getMessage()
        if min_level <= r.levelno <= max_level and msg not in out:
            out.append(msg)
    return out


def run_recipe(name: str, fn: Callable[..., Any], *args, **kwargs
               ) -> tuple[Any, list[str]]:
    """Run a recipe and return its result with the warnings it logged.

    Parameters
    ----------
    name : str
        The recipe name, used in the error messages.
    fn : callable
        The recipe (or any function following the recipe convention:
        0 on failure, with the reason logged).

    Returns
    -------
    tuple
        The result of ``fn`` and the list of the warning messages.

    Raises
    ------
    ToolError
        If ``fn`` returns 0 or ``None`` or raises; the message lists the
        errors it logged.
    """
    with collect_log() as records:
        try:
            result = fn(*args, **kwargs)
        except Exception as e:
            errors = _messages(records, logging.ERROR)
            raise ToolError("; ".join(
                errors + [f"{name} failed: {type(e).__name__}: {e}"])) from e
    warnings = _messages(records, logging.WARNING, logging.WARNING)
    if result is None or (isinstance(result, int) and result == 0):
        errors = _messages(records, logging.ERROR) or warnings
        detail = "; ".join(errors) if errors else "no reason was logged"
        raise ToolError(f"{name} failed: {detail}")
    return result, warnings
