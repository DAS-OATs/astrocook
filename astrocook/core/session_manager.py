import logging
from typing import List, Optional, Union

# --- V2 Imports ---
from astrocook.core.session import SessionV2
from astrocook.core.structures import HistoryLogV2, V1LogArtifact
# --- V1 Imports (for type hinting only) ---
from astrocook.legacy.gui_log import GUILog 

# Define the type for the log manager
LogManager = Union[HistoryLogV2, V1LogArtifact, GUILog]

class SessionHistory:
    """Manages the list of states and current index for one session lineage."""
    
    def __init__(self, initial_state: SessionV2, log_object: LogManager):
        if not isinstance(initial_state, SessionV2):
            raise TypeError("SessionHistory must be initialized with a SessionV2 object.")
        if not isinstance(log_object, (HistoryLogV2, V1LogArtifact, GUILog)):
             raise TypeError(f"SessionHistory initialized with invalid log object: {type(log_object)}")
        
        self.states: List[SessionV2] = [initial_state]
        self.current_index: int = 0 # Points to the active state
        self.log_manager: LogManager = log_object
        self._display_name: str = self.states[0].name if self.states else "Unnamed Session"

        # Link the log manager to the session state for recipes to access
        # We use a new attribute 'log_manager' to avoid V1 conflicts
        initial_state.log_manager = self.log_manager
        # Also set the old .log for V1 recipes if they're still used
        if isinstance(log_object, GUILog):
            initial_state.log = log_object

    @property
    def current_state(self) -> SessionV2:
        """Returns the currently active SessionV2 state."""
        if 0 <= self.current_index < len(self.states):
            return self.states[self.current_index]
        logging.error("History index out of bounds!")
        return self.states[-1] if self.states else None # Fallback

    @property
    def display_name(self) -> str:
        """Returns the name of the session (e.g., from the initial state)."""
        return self._display_name
    
    @display_name.setter
    def display_name(self, new_name: str):
        """ Allows the GUI to set the display name. """
        # --- THIS IS THE FIX ---
        self._display_name = new_name
        # --- END FIX ---

    def add_state(self, new_state: SessionV2):
        """
        Adds a new state, truncating redo history.
        NOTE: This method assumes the log entry has *already* been added
        by the recipe that created new_state.
        """
        
        # 1. Truncate stale future states
        if self.current_index < len(self.states) - 1:
            logging.debug(f"Truncating SessionHistory states from index {self.current_index + 1}")
            del self.states[self.current_index + 1:]
        
        # 2. Ensure the new state has the correct log manager reference
        if not hasattr(new_state, 'log_manager') or new_state.log_manager is not self.log_manager:
            logging.warning("SessionHistory.add_state: New state had incorrect log reference. Forcibly linking.")
            new_state.log_manager = self.log_manager
            # Also re-link the old .log attribute
            if isinstance(self.log_manager, GUILog):
                new_state.log = self.log_manager

        # 3. Append the new state
        self.states.append(new_state)
        
        # 4. Update the index
        self.current_index = len(self.states) - 1
        
        # 5. Log manager index is assumed to be handled by the recipe's
        #    add_entry call (which truncates the log)

    def can_undo(self) -> bool:
        return self.current_index > 0

    def can_redo(self) -> bool:
        return self.current_index < len(self.states) - 1

    def undo(self) -> Optional[SessionV2]:
        """Moves index back and returns the previous state."""
        if self.can_undo():
            self.current_index -= 1
            
            # --- Sync with V2 Log ---
            if isinstance(self.log_manager, HistoryLogV2):
                self.log_manager.undo()
            # ------------------------
                
            return self.current_state
        return None

    def redo(self) -> Optional[SessionV2]:
        """Moves index forward and returns the next state."""
        if self.can_redo():
            self.current_index += 1

            # --- Sync with V2 Log ---
            if isinstance(self.log_manager, HistoryLogV2):
                self.log_manager.redo()
            # ------------------------

            return self.current_state
        return None