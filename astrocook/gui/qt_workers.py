import logging
import re
import threading
from PySide6.QtCore import QObject, QRunnable, Signal
from typing import TYPE_CHECKING, Dict

from astrocook.core.session import SessionV2
from astrocook.core.session_manager import SessionHistory
from astrocook.core.structures import HistoryLogV2
from astrocook.gui.log_script import ScriptComments
from astrocook.recipes.errors import RecipeCancelled, RecipeError
from astrocook.recipes.params import resolve_params
from astrocook.recipes.prerequisites import collect_provided
from astrocook.recipes.progress import cancellable
from astrocook.recipes.registry import find_recipe
from astrocook.recipes.result import run

if TYPE_CHECKING:
    from .main_window import MainWindowV2

logger = logging.getLogger(__name__)

# Define the regex for parsing script lines
SCRIPT_LINE_REGEX = re.compile(r'(\w+)\((.*)\)')


class WorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.
    - finished: Emits the new SessionHistory object on success.
    - error: Emits a tuple of (line_number, line_text, error_message).
    - progress: Emits a string to update the user (e.g., "Running line 5...")
    - stopped: Emits the name of the recipe stopped by the user (it left
      no new state)
    """
    finished = Signal(object) 
    error = Signal(tuple)
    progress = Signal(str)
    stopped = Signal(str)


class RecipeWorker(QRunnable):
    """
    Worker thread for running a single recipe.
    """
    def __init__(self, 
                 session: SessionV2,
                 category: str,
                 recipe_name: str,
                 params: dict,
                 alias_map: dict):
        super().__init__()
        
        self.signals = WorkerSignals()
        #: The RecipeResult of the recipe, once it has run (phase 6-bis, PR 4)
        self.result = None
        self.session = session
        self.category = category
        self.recipe_name = recipe_name
        self.params = params
        self.alias_map = alias_map
        #: Set by the Stop: the recipe, run in :func:`cancellable`, stops at
        #: its next check (phase 6-ter, PR 7-quater; it was a ``_stop_flag``
        #: written in the session)
        self._cancel = threading.Event()

    def stop(self):
        """Thread-safe method to request a stop."""
        logger.info(f"Stopping {self.recipe_name}.")
        self._cancel.set()
        
    def run(self):
        """The main work. This runs on the background thread."""
        try:
            logger.debug(f"Running {self.category}.{self.recipe_name}.")

            # 1. Find the recipe method
            recipe_instance = getattr(self.session, self.category, None)
            if not recipe_instance:
                raise AttributeError(f"Session object has no attribute '{self.category}'")

            recipe_method = getattr(recipe_instance, self.recipe_name, None)
            if not callable(recipe_method):
                raise AttributeError(f"Recipe instance has no callable method '{self.recipe_name}'")

            # 2. Add alias_map if needed
            if self.recipe_name in ("apply_expression", "mask_expression"):
                self.params['alias_map'] = self.alias_map

            # 3. Run the recipe, with its structured result (phase 6-bis, PR 4)
            with cancellable(self._cancel):
                if find_recipe(self.recipe_name) is not None:
                    self.result = run(self.session, f"{self.category}.{self.recipe_name}",
                                      **self.params)
                    new_session_state = (self.result.session if self.result.session is not None
                                         else self.result.value)
                    if isinstance(new_session_state, tuple):
                        new_session_state = list(new_session_state)
                else:
                    new_session_state = recipe_method(**self.params)

            if self._cancel.is_set():
                # Stopped after its last check: its result is not kept
                raise RecipeCancelled("stopped after its last step: the session is as before it")

            if new_session_state is None or isinstance(new_session_state, int):
                # A recipe raises when it fails (phase 6-bis, PR 3): no result
                # is a bug, shown with its traceback
                raise RuntimeError(f"{self.recipe_name} gave no result ({new_session_state!r}).")
            
            # 4. Success!
            self.signals.finished.emit((new_session_state, self.recipe_name, self.params))
            
        except Exception as e:
            # 5. Failure! Emit the error
            
            import traceback
            
            if isinstance(e, RecipeCancelled):
                # The Stop of the user: no new state, no error dialog
                logger.warning(f"{self.recipe_name} was stopped.")
                self.signals.stopped.emit(self.recipe_name)
                return

            if isinstance(e, RecipeError):
                # The recipe could not do what was asked: the message names
                # the recipe (and the step of a pipeline) and the reason.
                # Logged here, once (phase 6-bis, PR 3)
                logger.error(str(e))
                self.signals.error.emit(("Recipe Aborted", str(e), ""))
                return

            # Check if this is a standard "user error" (like empty split, bad input)
            is_user_error = isinstance(e, (ValueError, TypeError))
            
            if is_user_error:
                # Case A: User Error. Log mildly, send NO traceback.
                logger.warning(f"{self.recipe_name} stopped: {e}")
                
                # Signal: (Title, Message, Traceback)
                self.signals.error.emit((
                    "Recipe Aborted",  # Milder title
                    str(e),           # Simple message
                    ""                # Empty traceback
                ))
            else:
                # Case B: Real Bug. Log severely, send FULL traceback.
                logger.error(f"Unexpected error in {self.recipe_name}: {e}", exc_info=True)
                
                self.signals.error.emit((
                    f"Recipe Error: {self.recipe_name}", # Scary title
                    f"An unexpected error occurred:\n{e}",
                    traceback.format_exc() # Full traceback
                ))
            
class ScriptWorker(QRunnable):
    """
    Worker thread for running a full analysis script.
    Inherits from QRunnable to run on the QThreadPool.
    """
    def __init__(self, 
                 script_text: str, 
                 initial_session_state: SessionV2, 
                 main_window: 'MainWindowV2'):
        super().__init__()
        
        self.signals = WorkerSignals()
        self.script_text = script_text
        self.initial_state = initial_session_state
        self.main_window = main_window # For context (aliases, stubs)

    def _build_alias_map(self) -> Dict[str, str]:
        """Generates a map of short aliases (s1, s2) to full session names."""
        alias_map = {}
        s_index = 1
        if self.main_window and self.main_window.active_history:
            current_session_name = self.main_window.active_history.display_name
            for hist in self.main_window.session_histories:
                if hist.display_name != current_session_name:
                    alias = f"s{s_index}"
                    alias_map[alias] = hist.display_name
                    s_index += 1
        return alias_map

    def run(self):
        """The main work. This runs on the background thread."""
        try:
            logger.debug("Script worker started.")

            # 1. Create a new history to build into; the log of the archive
            # the session was opened from stays read only (F109)
            new_log = HistoryLogV2()
            new_log.previous = getattr(getattr(self.initial_state, 'log_manager', None),
                                       'previous', None)
            # We must create a true copy of the initial state
            initial_state_copy = self.initial_state.with_new_spectrum(
                self.initial_state.spec
            ).with_new_system_list(
                self.initial_state.systs
            )
            # Re-link the new state to the GUI
            initial_state_copy._gui = self.main_window
            
            new_history = SessionHistory(initial_state_copy, new_log)
            current_processing_state = initial_state_copy
            
            alias_map = self._build_alias_map() # Build alias map once

            # 2. Process the script line by line
            lines = self.script_text.splitlines()
            # The reasons and the notes are not steps, but are kept in the
            # new log (phase 6-ter, PR 5a)
            comments = ScriptComments()
            for i, line in enumerate(lines):
                line = line.strip()
                if not line or line.startswith('#'):
                    comments.read(line)
                    continue # Skip comments and empty lines
                for kind, text in comments.notes():
                    new_log.add_note(text, kind=kind)
                reason = comments.reason()

                self.signals.progress.emit(f"Running line {i+1}: {line}")

                # 3a. Parse the line
                match = SCRIPT_LINE_REGEX.match(line)
                if not match:
                    raise ValueError("Invalid syntax. Expected 'recipe_name(args)'.")
                
                recipe_name = match.group(1)
                args_str = match.group(2)
                
                params_dict = {}
                if args_str:
                    arg_pairs = re.findall(r"(\w+)\s*=\s*([^,]+)", args_str)
                    for key, val in arg_pairs:
                        key = key.strip()
                        val = val.strip()
                        if (val.startswith("'") and val.endswith("'")) or \
                           (val.startswith('"') and val.endswith('"')):
                            params_dict[key] = val[1:-1]
                        else:
                            params_dict[key] = val 

                # 3b. Find and run the recipe
                category = self.main_window.recipe_category_map.get(recipe_name)
                if not category:
                    raise ValueError(f"Recipe '{recipe_name}' not found.")
                
                recipe_instance = getattr(current_processing_state, category)
                recipe_method = getattr(recipe_instance, recipe_name)

                # Add alias map if needed
                if recipe_name in ("apply_expression", "mask_expression"):
                    params_dict['alias_map'] = alias_map
                
                with collect_provided() as provided:
                    new_session_state = recipe_method(**params_dict)
                
                if new_session_state is None or isinstance(new_session_state, int):
                    # A recipe raises when it fails (phase 6-bis, PR 3)
                    raise RuntimeError(f"{recipe_name} gave no result ({new_session_state!r}).")
                    
                # 3c. Success! Add to the new history
                params_to_log = params_dict.copy()
                params_to_log.pop('alias_map', None)
                # What was provided first: its own state and entry (PR 5)
                for step in provided:
                    new_log.add_entry(step.recipe.split('.')[1], dict(step.params))
                    new_history.add_state(step.session)
                new_log.add_entry(recipe_name, resolve_params(recipe_name, params_to_log),
                                  reason=reason)

                new_history.add_state(new_session_state)
                current_processing_state = new_session_state

            for kind, text in comments.notes():
                new_log.add_note(text, kind=kind)

            # 4. Success! Emit the final history
            self.signals.finished.emit(new_history)

        except Exception as e:
            # 5. Failure! Emit the error
            line_info = f"line {i+1}: {line}" if 'i' in locals() else "N/A"
            # A RecipeError names the recipe and the reason: no traceback
            # (phase 6-bis, PR 3); another error is a bug, with it
            logger.error(f"The script failed at {line_info}: {e}",
                         exc_info=not isinstance(e, RecipeError))
            self.signals.error.emit((i+1, line, str(e)))

class SessionLoadSignals(QObject):
    """The signals of SessionLoadWorker: each file starting (index, path),
    loaded (index, path, session) or failed (index, path, message), and the
    end (loaded, failed, cancelled)."""
    started = Signal(int, str)
    loaded = Signal(int, str, object)
    failed = Signal(int, str, str)
    done = Signal(int, int, bool)


class SessionLoadWorker(QRunnable):
    """Loads several files one after the other in a background thread, so
    that the window stays responsive and can show the progress (phase
    5-bis: the files were loaded in the thread of the interface, with only
    the wait cursor). ``cancel`` skips the files not yet started."""

    def __init__(self, paths, loader):
        super().__init__()
        self.paths = list(paths)
        self.loader = loader  # path -> session, or raises / returns 0
        self.signals = SessionLoadSignals()
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def run(self):
        loaded = failed = 0
        for i, path in enumerate(self.paths):
            if self._cancelled:
                break
            self.signals.started.emit(i, path)
            try:
                session = self.loader(path)
                if session is None or session == 0:
                    raise ValueError("the loader could not read it")
                loaded += 1
                self.signals.loaded.emit(i, path, session)
            except Exception as e:
                failed += 1
                logger.error(f"Could not open {path}: {e}", exc_info=not isinstance(e, ValueError))
                self.signals.failed.emit(i, path, str(e))
        self.signals.done.emit(loaded, failed, self._cancelled)
