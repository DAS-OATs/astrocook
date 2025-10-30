import logging
from typing import List, Optional

from ..session import SessionV2
from ...v1.gui_log import GUILog

class SessionHistory:
    """Manages the list of states and current index for one session lineage."""
    def __init__(self, initial_state: SessionV2, gui_log: GUILog):
        if not isinstance(initial_state, SessionV2):
            raise TypeError("SessionHistory must be initialized with a SessionV2 object.")
        if not isinstance(gui_log, GUILog): # <<< Add check
            raise TypeError("SessionHistory must be initialized with a GUILog object.")
        
        self.states: List[SessionV2] = [initial_state]
        self.current_index: int = 0 # Points to the active state
        self.log: GUILog = gui_log

        # Ensure the initial state uses this history's log instance
        initial_state.log = self.log

    @property
    def current_state(self) -> SessionV2:
        """Returns the currently active SessionV2 state."""
        if 0 <= self.current_index < len(self.states):
            return self.states[self.current_index]
        # Should ideally not happen if index is managed correctly
        logging.error("History index out of bounds!")
        return self.states[-1] if self.states else None # Fallback

    @property
    def display_name(self) -> str:
        """Returns the name of the session (e.g., from the initial state)."""
        # Could potentially be made editable later
        return self.states[0].name if self.states else "Unnamed Session"

    def add_state(self, new_state: SessionV2):
        """Adds a new state, truncating redo history."""
        # Truncate states after the current index
        if self.current_index < len(self.states) - 1:
            del self.states[self.current_index + 1:]
        
        # Ensure the new state (which should have a ref) points to this log
        if new_state.log is not self.log:
            logging.warning("SessionHistory.add_state: New state had incorrect log reference. Forcibly linking.")
            new_state.log = self.log

        # Append the new state
        self.states.append(new_state)
        # Update the index to point to the new state
        self.current_index = len(self.states) - 1

    def can_undo(self) -> bool:
        return self.current_index > 0

    def can_redo(self) -> bool:
        return self.current_index < len(self.states) - 1

    def undo(self) -> Optional[SessionV2]:
        """Moves index back and returns the previous state."""
        if self.can_undo():
            self.current_index -= 1
            return self.current_state
        return None

    def redo(self) -> Optional[SessionV2]:
        """Moves index forward and returns the next state."""
        if self.can_redo():
            self.current_index += 1
            return self.current_state
        return None