from PySide6.QtCore import QObject, Signal

from astrocook.recipes.result import run

class AppController(QObject):
    """
    The Public API exposed to In-App Scripts.
    This acts as the 'app' object in the script context.
    """
    # Signals to update GUI from script thread
    log_message = Signal(str)
    refresh_request = Signal()

    def __init__(self, main_window):
        super().__init__()
        self._mw = main_window

    def get_session(self, name: str):
        """Retrieves a loaded session by its tab name."""
        for i in range(self._mw.tab_widget.count()):
            # Assuming main_window stores sessions in a retrievable way
            # This depends on your specific MainWindow implementation
            # Ideally, main_window has a mapping or list of active wrappers
            wrapper = self._mw.get_tab_wrapper(i) 
            if wrapper and wrapper.name == name:
                return wrapper.session
        return None

    def get_active_session(self):
        """Returns the currently visible session."""
        return self._mw.get_current_session()

    def run_recipe(self, session, recipe: str, method: str, **kwargs):
        """
        Executes a recipe on a session and returns its result.
        Example: app.run_recipe(s, 'absorbers', 'optimize_system', uuid=..., max_components=5)

        The recipe runs here, on ``session`` and in the thread of the
        script, through the registry (``astrocook.recipes.run``), as in the
        MCP server: it does not go through the dispatcher of the window (a
        worker thread, a progress dialog, the active session), which
        returns nothing to the caller (F101). The result is the new
        session (the first, for a recipe that returns several), or the
        value of a recipe that returns none; the window and its history
        are not changed until the script calls :meth:`create_tab`.
        """
        self.log_message.emit(f"Running {recipe}.{method}...")
        result = run(session, f"{recipe}.{method}", **kwargs)
        return result.session if result.session is not None else result.value

    def create_tab(self, session, name: str):
        """Adds a result session back to the GUI as a new tab."""
        # Must use QMetaObject.invokeMethod if running from a thread,
        # but for simple scripts, direct call might work if main thread blocked
        # (Better to emit a signal for thread safety)
        self._mw.add_session_tab(session, name)
        self.refresh_request.emit()