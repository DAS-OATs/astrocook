import logging
import json
from PySide6.QtWidgets import (
    QDialog, QMessageBox, QVBoxLayout, QHBoxLayout, QTextEdit, QPushButton, QDialogButtonBox, 
    QListWidgetItem, QFileDialog, QMenu
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QKeySequence, QTextCursor
from typing import TYPE_CHECKING, Union, Optional

# --- V2 Imports ---
from ..structures import HistoryLogV2, V1LogArtifact
from ..utils import get_recipe_schema # For inserting recipes
# --- V1 Imports ---
from ...v1.gui_log import GUILog

if TYPE_CHECKING:
    from .main_window import MainWindowV2 # Import main window for typing

# Define the type for the log manager
LogManager = Union[HistoryLogV2, V1LogArtifact, GUILog]


class LogScripterDialog(QDialog):
    """
    A non-modal dialog to display, edit, and save the session's log history
    as a runnable script.
    """
    def __init__(self, log_object: LogManager, session_name: str, parent: 'MainWindowV2'):
        super().__init__(parent)
        self.main_window = parent # Store reference to MainWindowV2
        self.setWindowTitle(f"Log Scripter: {session_name}")
        
        self.log_object = log_object # Save the log manager

        self.layout = QVBoxLayout(self)
        
        # --- 1. New Button Bar ---
        button_layout = QHBoxLayout()
        
        self.undo_button = QPushButton("Undo")
        self.undo_button.setToolTip("Undo last action in the main session")
        self.undo_button.setShortcut(QKeySequence.Undo)
        self.undo_button.clicked.connect(self._on_undo)
        button_layout.addWidget(self.undo_button)

        self.redo_button = QPushButton("Redo")
        self.redo_button.setToolTip("Redo last action in the main session")
        self.redo_button.setShortcut(QKeySequence.Redo)
        self.redo_button.clicked.connect(self._on_redo)
        button_layout.addWidget(self.redo_button)
        
        button_layout.addStretch()

        self.run_all_button = QPushButton("Run All")
        self.run_all_button.setToolTip("Re-run the entire script from the original data")
        self.run_all_button.setStyleSheet(
            #"QPushButton { background-color: #4CAF50; color: white; border-radius: 5px; }"
            "QPushButton { font-weight: bold; }"
        )
        self.run_all_button.clicked.connect(self._on_run_all)
        button_layout.addWidget(self.run_all_button)

        self.insert_recipe_button = QPushButton("Insert Recipe...")
        self.insert_recipe_button.setToolTip("Insert a recipe template at the cursor")
        self.insert_recipe_button.clicked.connect(self._on_insert_recipe)
        button_layout.addWidget(self.insert_recipe_button)

        self.save_log_button = QPushButton("Save Log As...")
        self.save_log_button.setToolTip("Save the content of this scripter to a .json file")
        self.save_log_button.clicked.connect(self._on_save_log)
        button_layout.addWidget(self.save_log_button)

        self.layout.addLayout(button_layout)

        # --- 2. Scripter Text Edit ---
        self.scripter_edit = QTextEdit()
        self.scripter_edit.setReadOnly(False) # Make it editable
        #self.scripter_edit.setFontFamily("Monospace")
        # Connect textChanged to update the button state
        self.scripter_edit.textChanged.connect(self._update_button_state)
        self.layout.addWidget(self.scripter_edit)

        # --- 3. Close Button ---
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        button_box.rejected.connect(self.reject)
        self.layout.addWidget(button_box)

        self._populate_scripter() # Populate on init
        self.resize(700, 500)
        self._update_button_state() # Set initial state for undo/redo
    
    def _build_alias_map(self) -> dict:
        """Generates a map of short aliases (s1, s2) to full session names."""
        alias_map = {}
        s_index = 1
        if self.main_window and self.main_window.active_history:
            # Use the active_history's name as the one to exclude
            current_session_name = self.main_window.active_history.display_name
            for hist in self.main_window.session_histories:
                if hist.display_name != current_session_name:
                    alias = f"s{s_index}"
                    alias_map[alias] = hist.display_name
                    s_index += 1
        return alias_map

    def set_log_object(self, new_log_object: Optional[LogManager]):
        """
        Updates the log scripter to point to a new log manager
        (e.g., when the active session changes).
        """
        new_title = "Log Scripter"
        if new_log_object is None:
            self.log_object = None
            self.scripter_edit.setPlainText("# No active session.")
            self.setWindowTitle(new_title)
        elif new_log_object is not self.log_object:
            if self.main_window and self.main_window.active_history:
                 if self.main_window.active_history.log_manager is new_log_object:
                    new_title = f"Log Scripter: {self.main_window.active_history.display_name}"
            
            self.log_object = new_log_object
            self.setWindowTitle(new_title) 
            self.refresh() 
        
        self._update_button_state() # Update buttons on change

    def _populate_scripter(self):
        """
        Parses the log_object and populates the text edit.
        """
        self.scripter_edit.clear()
        
        if self.log_object is None:
            self.scripter_edit.setPlainText("# No log available.")
            return

        # --- *** START V2/V1 LOGIC *** ---
        if hasattr(self.log_object, 'is_v2_log') and self.log_object.is_v2_log:
            self._populate_v2_script()
        
        elif hasattr(self.log_object, 'is_v2_log') and not self.log_object.is_v2_log:
            self._populate_v1_script(self.log_object.v1_json)
            
        elif isinstance(self.log_object, GUILog):
            logging.warning("LogScripterDialog received raw GUILog, parsing .json")
            self._populate_v1_script(self.log_object.json)
            
        else:
            self.scripter_edit.setPlainText(f"# Error: Unknown log object type: {type(self.log_object)}")
            logging.error(f"LogScripterDialog: Unknown log object type: {type(self.log_object)}")
        # --- *** END V2/V1 LOGIC *** ---
        
        self.scripter_edit.moveCursor(QTextCursor.MoveOperation.End)

    def _populate_v2_script(self):
        """Populates the scripter from a HistoryLogV2 object."""
        
        entries = self.log_object.entries
        current_idx = self.log_object.current_index
        
        if not entries:
            #self.scripter_edit.setPlainText("# Log is empty.")
            return

        script_lines = []
        for i, entry in enumerate(entries):
            # Format as "recipe_name(param1=value1, param2='value2')"
            params_str = ", ".join(f"{k}='{v}'" if isinstance(v, str) else f"{k}={v}" 
                                   for k, v in entry.params.items())
            line = f"{entry.recipe_name}({params_str})"
            
            # Add comment for "undone" actions
            if i > current_idx:
                line = f"# {line}"
                
            script_lines.append(line)

        self.scripter_edit.setPlainText("\n".join(script_lines))

    def _populate_v1_script(self, v1_json: dict):
        """Populates the scripter from a V1 log JSON dictionary."""
        try:
            commands = v1_json.get('set_menu', [])
        except Exception as e:
            logging.error(f"Failed to parse V1 log JSON: {e}")
            self.scripter_edit.setPlainText("# Error: Could not parse V1 log.")
            return

        if not commands:
            self.scripter_edit.setPlainText("# V1 Log is empty.")
            return
            
        script_lines = []
        for cmd in commands:
            if cmd.get('recipe') == '_refresh':
                continue
            rec = cmd.get('recipe', 'N/A')
            params_dict = cmd.get('params', {})
            params_str = ", ".join(f"{k}='{v}'" if isinstance(v, str) else f"{k}={v}" 
                                   for k, v in params_dict.items())
            line = f"{rec}({params_str})"
            script_lines.append(line)
        
        self.scripter_edit.setPlainText("\n".join(script_lines))

    def refresh(self):
        """
        Public method to clear and re-populate the log list.
        """
        logging.debug(f"Refreshing LogScripterDialog for {self.windowTitle()}")
        self._populate_scripter()
        self._update_button_state()
        
    def _update_button_state(self):
        """Syncs the Undo/Redo/Run buttons with the main window's state."""
        can_undo = False
        can_redo = False
        can_run = False # <-- New check
        
        if self.main_window and self.main_window.active_history:
            can_undo = self.main_window.active_history.can_undo()
            can_redo = self.main_window.active_history.can_redo()

        # Check the *text content* of the scripter, not the log object
        if self.scripter_edit:
            script_text = self.scripter_edit.toPlainText().strip()
            # Check if there is any non-comment, non-empty text
            if script_text and not script_text.startswith("#"):
                can_run = True
        
        self.undo_button.setEnabled(can_undo)
        self.redo_button.setEnabled(can_redo)
        self.run_all_button.setEnabled(can_run)
        
    def _on_undo(self):
        if self.main_window:
            self.main_window._undo_last_action()
            # The main window's undo method will refresh us
            
    def _on_redo(self):
        if self.main_window:
            self.main_window._redo_last_action()
            # The main window's redo method will refresh us

    def _on_run_all(self):
        """Tells the main window to re-run the *entire* script."""
        if not self.main_window:
            return
            
        script_text = self.scripter_edit.toPlainText()
        
        reply = QMessageBox.question(
            self,
            "Run Script",
            "This will replace your current session history by re-running this script from the original data.\n\nAre you sure you want to proceed?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel
        )
        
        if reply == QMessageBox.StandardButton.Yes:
            self.main_window.run_script(script_text)

    def _on_save_log(self):
        """Saves the scripter's text content to a file."""
        if not self.log_object:
            QMessageBox.warning(self, "No Log", "Nothing to save.")
            return
            
        # Suggest a filename
        default_name = "astrocook_log.json"
        if self.main_window and self.main_window.active_history:
            default_name = f"{self.main_window.active_history.display_name}_log.json"
            
        file_name, _ = QFileDialog.getSaveFileName(
            self, "Save Log Scripter", default_name, "JSON Log (*.json);;Text File (*.txt)"
        )
        
        if file_name:
            try:
                text_content = self.scripter_edit.toPlainText()
                # Here we could save as JSON, or just plain text.
                # For now, let's just save the text.
                # To save as JSON, we'd have to parse the text.
                with open(file_name, 'w') as f:
                    f.write(text_content)
                logging.info(f"Log scripter content saved to {file_name}")
            except Exception as e:
                logging.error(f"Failed to save log file: {e}")
                QMessageBox.critical(self, "Save Error", f"Could not save file:\n{e}")

    def _on_insert_recipe(self):
        """Shows a menu of available recipes to insert as a template."""
        
        # This is a complex task. For now, I'll provide a placeholder
        # In the future, this would scan all recipe schemas.
        
        menu = QMenu(self)
        
        # --- Create a simple menu ---
        try:
            # Get Edit schemas
            from ..recipes.edit import EDIT_RECIPES_SCHEMAS
            edit_menu = menu.addMenu("Edit")
            for name in EDIT_RECIPES_SCHEMAS.keys():
                action = edit_menu.addAction(name)
                action.triggered.connect(lambda checked=False, n=name: self._insert_template(n))
                
            # Get Flux schemas
            from ..recipes.flux import FLUX_RECIPES_SCHEMAS
            flux_menu = menu.addMenu("Flux")
            for name in FLUX_RECIPES_SCHEMAS.keys():
                action = flux_menu.addAction(name)
                action.triggered.connect(lambda checked=False, n=name: self._insert_template(n))

        except Exception as e:
            logging.error(f"Failed to build recipe insert menu: {e}")
            action = menu.addAction("Error building menu")
            action.setEnabled(False)

        # Show the menu at the button's position
        cursor_pos = self.insert_recipe_button.mapToGlobal(self.insert_recipe_button.rect().bottomLeft())
        menu.exec(cursor_pos)
        
    def _insert_template(self, recipe_name: str):
        """Inserts a text template for the given recipe."""
        
        # Try to find the schema
        schema = None
        params_str = "..." # Default
        try:
            from ..recipes.edit import EDIT_RECIPES_SCHEMAS
            if recipe_name in EDIT_RECIPES_SCHEMAS:
                schema = EDIT_RECIPES_SCHEMAS[recipe_name]
            else:
                from ..recipes.flux import FLUX_RECIPES_SCHEMAS
                if recipe_name in FLUX_RECIPES_SCHEMAS:
                    schema = FLUX_RECIPES_SCHEMAS[recipe_name]
            
            if schema and 'params' in schema:
                params = []
                for p in schema['params']:
                    params.append(f"{p['name']}={p['default']}")
                params_str = ", ".join(params)
                
        except Exception:
            pass # Use default "..."
            
        template = f"{recipe_name}({params_str})\n"
        self.scripter_edit.textCursor().insertText(template)
        self.scripter_edit.setFocus()

        # Manually update the button state after inserting text
        self._update_button_state()