import logging
import re
import numpy as np
from PySide6.QtCore import QObject, QRunnable, Signal
from typing import TYPE_CHECKING, Dict

# --- Import all schemas and the map ---
from astrocook.recipes.absorbers import ABSORBERS_RECIPES_SCHEMAS
from astrocook.recipes.continuum import CONTINUUM_RECIPES_SCHEMAS
from astrocook.recipes.edit import EDIT_RECIPES_SCHEMAS
from astrocook.recipes.file import FILE_RECIPES_SCHEMAS
from astrocook.recipes.flux import FLUX_RECIPES_SCHEMAS
from astrocook.core.session import SessionV2
from astrocook.core.session_manager import SessionHistory
from astrocook.core.structures import HistoryLogV2
from astrocook.legacy.gui_log import GUILog
from astrocook.legacy.defaults import Defaults

if TYPE_CHECKING:
    from .main_window import MainWindowV2

# Define a map of all known recipe modules
ALL_SCHEMAS = {
    'file': FILE_RECIPES_SCHEMAS,
    'edit': EDIT_RECIPES_SCHEMAS,
    'flux': FLUX_RECIPES_SCHEMAS,
    'continuum': CONTINUUM_RECIPES_SCHEMAS,
    'absorbers': ABSORBERS_RECIPES_SCHEMAS
}
# Define the regex for parsing script lines
SCRIPT_LINE_REGEX = re.compile(r'(\w+)\((.*)\)')


class WorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.
    - finished: Emits the new SessionHistory object on success.
    - error: Emits a tuple of (line_number, line_text, error_message).
    - progress: Emits a string to update the user (e.g., "Running line 5...")
    """
    finished = Signal(object) 
    error = Signal(tuple)
    progress = Signal(str)


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
        self.session = session
        self.category = category
        self.recipe_name = recipe_name
        self.params = params
        self.alias_map = alias_map
        self._is_stopped = False

    def stop(self):
        """Thread-safe method to request a stop."""
        logging.info(f"Stopping recipe '{self.recipe_name}'...")
        self._is_stopped = True
        
        # Inject the stop flag into the session object so recipes can see it
        # This acts as a 'shared variable' between GUI thread and Physics thread
        self.session._stop_flag = True
        
    def run(self):
        """The main work. This runs on the background thread."""
        try:
            logging.info(f"RecipeWorker: Running {self.category}.{self.recipe_name}")

            # Reset stop flag on session start
            self.session._stop_flag = False
            
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

            # 3. Run the recipe
            new_session_state = recipe_method(**self.params)

            if self._is_stopped or getattr(self.session, '_stop_flag', False):
                logging.warning(f"Recipe '{self.recipe_name}' was stopped by user.")
                
                # Use the returned session (likely the original one) or fallback
                result_session = new_session_state if new_session_state else self.session
                
                # [CRITICAL] Tag the RESULT object so MainWindow knows to abort history update
                result_session._stop_flag = True
                
                self.signals.finished.emit((result_session, self.recipe_name, self.params))
                return

            if not new_session_state or new_session_state == 0:
                raise ValueError(f"Recipe failed to execute (returned 0). Check logs.")
            
            # 4. Success!
            self.signals.finished.emit((new_session_state, self.recipe_name, self.params))
            
        except Exception as e:
            # 5. Failure! Emit the error
            
            import traceback
            
            # Check if this is a standard "user error" (like empty split, bad input)
            is_user_error = isinstance(e, (ValueError, TypeError))
            
            if is_user_error:
                # Case A: User Error. Log mildly, send NO traceback.
                logging.warning(f"Recipe '{self.recipe_name}' aborted: {e}")
                
                # Signal: (Title, Message, Traceback)
                self.signals.error.emit((
                    "Recipe Aborted",  # Milder title
                    str(e),           # Simple message
                    ""                # Empty traceback
                ))
            else:
                # Case B: Real Bug. Log severely, send FULL traceback.
                logging.error(f"RecipeWorker crashed on {self.recipe_name}: {e}", exc_info=True)
                
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
            logging.info("ScriptWorker: Starting script run...")

            # 1. Create a new history to build into
            new_log = HistoryLogV2()
            # We must create a true copy of the initial state
            initial_state_copy = self.initial_state.with_new_spectrum(
                self.initial_state.spec
            ).with_new_system_list(
                self.initial_state.systs
            )
            # Re-link the new state to the GUI and V1 stubs
            initial_state_copy._gui = self.main_window
            initial_state_copy.log = GUILog(self.main_window.mock_gui_context) 
            initial_state_copy.defs = Defaults(self.main_window.mock_gui_context)
            
            new_history = SessionHistory(initial_state_copy, new_log)
            current_processing_state = initial_state_copy
            
            alias_map = self._build_alias_map() # Build alias map once

            # 2. Process the script line by line
            lines = self.script_text.splitlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue # Skip comments and empty lines

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
                
                new_session_state = recipe_method(**params_dict)
                
                if not new_session_state or new_session_state == 0:
                    raise ValueError(f"Recipe failed to execute (returned 0).")
                    
                # 3c. Success! Add to the new history
                params_to_log = params_dict.copy()
                params_to_log.pop('alias_map', None)
                new_log.add_entry(recipe_name, params_to_log) 

                new_history.add_state(new_session_state)
                current_processing_state = new_session_state

            # 4. Success! Emit the final history
            self.signals.finished.emit(new_history)

        except Exception as e:
            # 5. Failure! Emit the error
            line_info = f"line {i+1}: {line}" if 'i' in locals() else "N/A"
            logging.error(f"ScriptWorker: Failed on {line_info}\nError: {e}", exc_info=True)
            self.signals.error.emit((i+1, line, str(e)))