"""Progress notifications of a running tool (phase 6-ter, PR 4).

The recipes write ``>> PROGRESS: n`` in the log while they run
(:mod:`astrocook.recipes.progress`, read by the GUI too); a tool sends them
to the client as MCP progress notifications, with the last line of the log
as the message, the status the GUI shows under its bar. A recipe that
writes none sends two: one when it starts, one when it ends. Nothing here
names a recipe. A client that did not ask for progress (no progress token)
receives nothing: the SDK ignores the call.

Nothing is written to stdout: the notifications go through the transport of
the SDK, and the log to stderr.
"""
import logging
import time
from contextlib import contextmanager
from typing import Iterator

import anyio
import anyio.from_thread
from mcp.server.mcpserver import Context

from astrocook.recipes.progress import watch_progress

logger = logging.getLogger(__name__)

#: Progress is a percentage
TOTAL = 100


@contextmanager
def report_progress(ctx: Context | None, name: str) -> Iterator[None]:
    """Send ``ctx`` the progress of the recipe that runs in the block, in
    the thread that calls this (a synchronous tool runs in a worker thread
    of the SDK, from which the notifications are sent to its event loop).

    Parameters
    ----------
    ctx : Context or None
        The context of the call of the tool; None sends nothing.
    name : str
        What is running, for the messages of the start and the end.
    """
    if ctx is None:
        yield
        return
    last = -1
    open_ = True

    def send(percent: int, message: str) -> None:
        # The value of a notification must increase: repeated lines (e.g.
        # several recipes of a pipeline starting at 0) are not sent
        nonlocal last, open_
        percent = min(max(percent, 0), TOTAL)
        if not open_ or percent <= last:
            return
        try:
            anyio.from_thread.run(ctx.report_progress, percent, TOTAL, message)
        except (anyio.BrokenResourceError, anyio.ClosedResourceError):
            # The client went away: the recipe goes on, and says it once
            open_ = False
            logger.info("%s: the client is gone, no more progress is sent.", name)
            return
        last = percent

    start = time.monotonic()
    send(0, f"Running {name}")
    # The lines of the recipe stop at 99: 100 is the end of the tool, which
    # always has its message (a recipe writes 100 before it returns its result)
    with watch_progress(lambda percent, message: send(min(percent, TOTAL - 1), message)):
        yield
    send(TOTAL, f"{name} done in {time.monotonic() - start:.1f} s")
