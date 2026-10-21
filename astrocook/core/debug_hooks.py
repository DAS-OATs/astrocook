"""Hook for the debug plots requested by the core.

The core does not depend on the GUI (findings register F6): it hands the
plot data to a handler, which the GUI sets to the signal that draws them on
the main thread. Without a handler (scripts, tests) the request is dropped,
as the GUI signal did when nothing was connected to it.
"""
import logging
from typing import Callable, Optional

_debug_plot_handler: Optional[Callable[[dict], None]] = None


def set_debug_plot_handler(handler: Optional[Callable[[dict], None]]) -> None:
    """Set the function that receives the debug plot requests.

    Parameters
    ----------
    handler : callable or None
        Called with a dictionary of plot data (``'title'`` and the arrays to
        plot); ``None`` removes the current handler.
    """
    global _debug_plot_handler
    _debug_plot_handler = handler


def request_debug_plot(plot_data: dict) -> None:
    """Pass ``plot_data`` to the debug plot handler, if one is set."""
    if _debug_plot_handler is None:
        logging.debug(f"No handler for the debug plot: {plot_data.get('title', '')}")
        return
    _debug_plot_handler(plot_data)
