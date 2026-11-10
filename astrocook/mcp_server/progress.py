"""Progress notifications of a running tool (phase 6-ter, PR 4).

The recipes write ``>> PROGRESS: n`` in the log while they run
(:mod:`astrocook.recipes.progress`, read by the GUI too); a tool sends them
to the client as MCP progress notifications, with the last line of the log
as the message, the status the GUI shows under its bar. A recipe that
writes none sends two: one when it starts, one when it ends. Nothing here
names a recipe. A client that did not ask for progress (no progress token)
receives nothing: the SDK ignores the call.

A recipe that logs nothing for :data:`HEARTBEAT_S` seconds gets a heartbeat
(phase 6-ter, section 4.3 of ``ROADMAP.md``, F113): a notification that it
is still running, with the time elapsed, so that neither the agent nor the
user mistakes a long silent recipe (e.g. ``metals_auto``) for a hang, and a
client that resets its timeout on progress does not give up. The value of a
notification must increase, so a heartbeat adds a part of what is left
before the next percent (1/2, 3/4, ...): a line of the recipe is never
below it.

Nothing is written to stdout: the notifications go through the transport of
the SDK, and the log to stderr.
"""
import functools
import logging
import threading
import time
from contextlib import contextmanager
from typing import Any, Callable, Iterator

import anyio
import anyio.from_thread
import anyio.lowlevel
import anyio.to_thread
from mcp.server.mcpserver import Context

from astrocook.recipes.progress import watch_progress

logger = logging.getLogger(__name__)

#: Progress is a percentage
TOTAL = 100
#: Seconds of silence of a recipe before a heartbeat (F113): below the
#: timeouts of the clients, usually 60 s
HEARTBEAT_S = 20.0

#: What a notification raises when the client or the loop is gone
_GONE = tuple(e for e in (getattr(anyio, n, None) for n in (
    "BrokenResourceError", "ClosedResourceError", "RunFinishedError")) if e)


async def run_cancellable(fn: Callable[[threading.Event], Any]) -> Any:
    """Run ``fn(event)`` in a worker thread, and set ``event`` if the call
    of the tool is cancelled meanwhile (phase 6-ter, PR 7-quater).

    The SDK cancels the scope of a call when its client sends
    ``notifications/cancelled`` (or gives up), but a synchronous tool runs
    in a thread that cannot be interrupted: the cancellation arrived when
    the recipe ended, hours later on a long spectrum (F113). Here the call
    stops waiting for the thread when it is cancelled, and sets the event
    that the recipe checks (:func:`astrocook.recipes.progress.cancellable`):
    the recipe stops at its next check, and its result, if any, is dropped.
    """
    event = threading.Event()
    try:
        return await anyio.to_thread.run_sync(functools.partial(fn, event),
                                              abandon_on_cancel=True)
    except anyio.get_cancelled_exc_class():
        event.set()
        raise


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
    # The token of the event loop, so that the thread of the heartbeat can
    # send too (only the worker thread of the tool knows it)
    token = anyio.from_thread.run_sync(anyio.lowlevel.current_token)
    last = -1        # the last percent of the recipe
    sent = -1.0      # the last value sent (a heartbeat is above ``last``)
    beats = 0        # heartbeats since the last percent
    status = ""      # the last line of the recipe
    last_time = time.monotonic()
    open_ = True
    guard = threading.Lock()

    def notify(value: float, message: str) -> bool:
        nonlocal open_, sent, last_time
        try:
            anyio.from_thread.run(ctx.report_progress, value, TOTAL, message, token=token)
        except _GONE:
            # The client went away: the recipe goes on, and says it once
            open_ = False
            logger.info("%s: the client is gone, no more progress is sent.", name)
            return False
        sent, last_time = value, time.monotonic()
        return True

    def send(percent: int, message: str) -> None:
        # The value of a notification must increase: repeated lines (e.g.
        # several recipes of a pipeline starting at 0) are not sent
        nonlocal last, beats, status
        percent = min(max(percent, 0), TOTAL)
        with guard:
            status = message
            if not open_ or percent <= last:
                return
            if notify(percent, message):
                last, beats = percent, 0

    start = time.monotonic()
    done = threading.Event()

    def heartbeat() -> None:
        nonlocal beats
        while not done.wait(min(HEARTBEAT_S, 1.0)):
            with guard:
                if not open_ or done.is_set():
                    return
                if time.monotonic() - last_time < HEARTBEAT_S:
                    continue
                beats += 1
                value = last + 1 - 0.5 ** beats
                if value <= sent:    # below what floats tell apart
                    continue
                elapsed = time.monotonic() - start
                message = f"Still running {name}: {elapsed:.0f} s"
                if status:
                    message += f"; last: {status}"
                notify(value, message)

    send(0, f"Running {name}")
    beat = threading.Thread(target=heartbeat, name=f"heartbeat-{name}", daemon=True)
    beat.start()
    try:
        # The lines of the recipe stop at 99: 100 is the end of the tool, which
        # always has its message (a recipe writes 100 before it returns its result)
        with watch_progress(lambda percent, message: send(min(percent, TOTAL - 1), message)):
            yield
    finally:
        with guard:
            done.set()
        beat.join()
    send(TOTAL, f"{name} done in {time.monotonic() - start:.1f} s")
