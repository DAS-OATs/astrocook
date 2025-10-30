import logging
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QListWidget, QPushButton, QDialogButtonBox
)
from PySide6.QtCore import Qt
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...v1.gui_log import GUILog

class LogViewerDialog(QDialog):
    """
    A non-modal dialog to display the session's V1 log history
    in a human-readable format.
    """
    def __init__(self, log_object: 'GUILog', session_name: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Log History: {session_name}")
        self.setWindowFlags(self.windowFlags() | Qt.Window) # Ensure it's a modeless window

        self.log_object = log_object

        self.layout = QVBoxLayout(self)
        
        # 1. The List Widget to show parsed commands
        self.list_widget = QListWidget()
        self.layout.addWidget(self.list_widget)

        # 2. Close Button
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        button_box.rejected.connect(self.reject) # Close connects to reject
        self.layout.addWidget(button_box)

        # 3. Populate the list
        self._populate_list(log_object)
        
        # Set a reasonable default size
        self.resize(600, 600)

    def _populate_list(self, log_object: 'GUILog'):
        """
        Parses the GUILog's JSON and populates the list widget.
        """
        self.list_widget.clear()
        
        try:
            # Get the list of commands
            commands = log_object.json.get('set_menu', [])
        except Exception as e:
            logging.error(f"Failed to parse log JSON: {e}")
            self.list_widget.addItem("Error: Could not parse log history.")
            return

        if not commands:
            self.list_widget.addItem("Log is empty.")
            return

        # Parse each command dictionary
        for cmd in commands:
            # Skip V1 internal refresh commands
            if cmd.get('recipe') == '_refresh':
                continue

            # Get recipe name
            rec = cmd.get('recipe', 'N/A')
            # Capitalize first letter for readability
            rec_display = rec.capitalize() 

            # Format params: (key1=val1, key2=val2)
            params_dict = cmd.get('params', {})
            params_str = ", ".join(f"{k}={v}" for k, v in params_dict.items())
            
            # The human-readable string (prefix omitted)
            display_string = f"{rec_display} ({params_str})"
            
            self.list_widget.addItem(display_string)

        # Scroll to the bottom to show the most recent command
        self.list_widget.scrollToBottom()

    def refresh(self):
        """
        Public method to clear and re-populate the log list.
        """
        logging.debug(f"Refreshing LogViewerDialog for {self.windowTitle()}")
        self.list_widget.clear()
        self._populate_list(self.log_object)