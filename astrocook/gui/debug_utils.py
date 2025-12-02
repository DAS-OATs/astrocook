from PySide6.QtCore import QObject, Signal

class DebugPlotter(QObject):
    """
    A global QObject to emit plot requests from worker threads
    to the main GUI thread.
    """
    # Signal emits a dictionary of data to be plotted
    plot_requested = Signal(dict)

# Create a single, global instance that can be imported anywhere
GLOBAL_PLOTTER = DebugPlotter()