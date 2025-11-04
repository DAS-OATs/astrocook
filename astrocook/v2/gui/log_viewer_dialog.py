import logging
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QListWidget, QPushButton, QDialogButtonBox, QListWidgetItem
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor # <<< *** ADD THIS IMPORT ***
from typing import TYPE_CHECKING, Optional, Union

# --- V2 Imports ---
from ..structures import HistoryLogV2, V1LogArtifact
# --- V1 Imports ---
from ...v1.gui_log import GUILog

# Define the type for the log manager
LogManager = Union[HistoryLogV2, V1LogArtifact, GUILog]


class LogViewerDialog(QDialog):
    """
    A non-modal dialog to display the session's log history
    in a human-readable format.
    Supports both V2 (HistoryLogV2) and V1 (V1LogArtifact) logs.
    """
    def __init__(self, log_object: LogManager, session_name: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Log History: {session_name}")
        # self.setWindowFlags(self.windowFlags() | Qt.Window) # Removed to fix bugs
        
        self.log_object = log_object # Save the log manager

        self.layout = QVBoxLayout(self)
        self.list_widget = QListWidget()
        self.layout.addWidget(self.list_widget)

        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        button_box.rejected.connect(self.reject)
        self.layout.addWidget(button_box)

        self._populate_list() # Populate on init
        self.resize(600, 600)

    def set_log_object(self, new_log_object: Optional[LogManager]):
        """
        Updates the log viewer to point to a new log manager
        (e.g., when the active session changes).
        """
        if new_log_object is None:
            # Clear the view if no session is active
            self.log_object = None
            self.list_widget.clear()
            self.list_widget.addItem("No active session.")
            self.setWindowTitle("Log History")
        elif new_log_object is not self.log_object:
            # Find the new session's name from the main window
            new_title = "Log History"
            if hasattr(self, 'parent') and hasattr(self.parent(), 'active_history'):
                 if self.parent().active_history.log_manager is new_log_object:
                    new_title = f"Log History: {self.parent().active_history.display_name}"
            
            self.log_object = new_log_object
            self.setWindowTitle(new_title) # Update window title
            self.refresh() # Refresh the view with the new object

    def _populate_list(self):
        """
        Parses the log_object and populates the list widget.
        """
        self.list_widget.clear()
        
        # --- *** START V2/V1 LOGIC *** ---
        
        # Check if this is a V2 Log
        if hasattr(self.log_object, 'is_v2_log') and self.log_object.is_v2_log:
            self._populate_v2_list()
        
        # Check if this is a V1 Artifact
        elif hasattr(self.log_object, 'is_v2_log') and not self.log_object.is_v2_log:
            self._populate_v1_list(self.log_object.v1_json)
            
        # Fallback for raw V1 GUILog (from branching)
        elif isinstance(self.log_object, GUILog):
            logging.warning("LogViewerDialog received raw GUILog, parsing .json")
            self._populate_v1_list(self.log_object.json)
            
        else:
            self.list_widget.addItem("Error: Unknown log object type.")
            logging.error(f"LogViewerDialog: Unknown log object type: {type(self.log_object)}")

        # --- *** END V2/V1 LOGIC *** ---

        # Scroll to the bottom to show the most recent command
        self.list_widget.scrollToBottom()

    def _populate_v2_list(self):
        """Populates the list from a HistoryLogV2 object."""
        
        entries = self.log_object.entries
        current_idx = self.log_object.current_index
        
        if not entries:
            self.list_widget.addItem("Log is empty.")
            return

        gray_color = QColor(Qt.GlobalColor.gray)

        for i, entry in enumerate(entries):
            # Format the string
            rec_display = entry.recipe_name.capitalize()
            params_str = ", ".join(f"{k}={v}" for k, v in entry.params.items())
            display_string = f"{rec_display} ({params_str})"
            
            item = QListWidgetItem(display_string)
            
            # --- *** THIS IS THE SHADING LOGIC *** ---
            if i > current_idx:
                item.setForeground(gray_color)
            # ----------------------------------
                
            self.list_widget.addItem(item)

    def _populate_v1_list(self, v1_json: dict):
        """Populates the list from a V1 log JSON dictionary."""
        try:
            commands = v1_json.get('set_menu', [])
        except Exception as e:
            logging.error(f"Failed to parse V1 log JSON: {e}")
            self.list_widget.addItem("Error: Could not parse V1 log.")
            return

        if not commands:
            self.list_widget.addItem("V1 Log is empty.")
            return

        for cmd in commands:
            if cmd.get('recipe') == '_refresh':
                continue
            rec = cmd.get('recipe', 'N/A')
            rec_display = rec.capitalize() 
            params_dict = cmd.get('params', {})
            params_str = ", ".join(f"{k}={v}" for k, v in params_dict.items())
            display_string = f"{rec_display} ({params_str})"
            
            # V1 logs are always "active", no shading
            self.list_widget.addItem(display_string)

    def refresh(self):
        """
        Public method to clear and re-populate the log list.
        """
        logging.debug(f"Refreshing LogViewerDialog for {self.windowTitle()}")
        self.list_widget.clear()
        self._populate_list() # Uses self.log_object