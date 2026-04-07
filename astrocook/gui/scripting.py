import logging
import threading
from PySide6.QtCore import QObject, Signal, Slot

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
        Executes a recipe on a session.
        Example: app.run_recipe(s, 'absorbers', 'optimize_system', uuid=..., max_components=5)
        """
        self.log_message.emit(f"Running {recipe}.{method}...")
        
        # We reuse the MainWindow's central dispatcher
        # This handles parameter parsing and recipe instantiation
        result_session = self._mw._on_recipe_requested(
            recipe, method, kwargs, session_override=session
        )
        
        if result_session:
            return result_session
        else:
            raise RuntimeError(f"Recipe {recipe}.{method} failed.")

    def create_tab(self, session, name: str):
        """Adds a result session back to the GUI as a new tab."""
        # Must use QMetaObject.invokeMethod if running from a thread,
        # but for simple scripts, direct call might work if main thread blocked
        # (Better to emit a signal for thread safety)
        self._mw.add_session_tab(session, name)
        self.refresh_request.emit()