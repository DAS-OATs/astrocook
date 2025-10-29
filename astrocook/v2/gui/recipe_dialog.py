import logging
import webbrowser
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QFormLayout, QLabel, QLineEdit, QCheckBox,
    QPushButton, QDialogButtonBox, QWidget, QComboBox, QMessageBox # Added QMessageBox
)
from PySide6.QtGui import QIntValidator, QDoubleValidator, QFont
from PySide6.QtCore import Qt, QLocale

from ..session import SessionV2
from ..utils import get_recipe_schema, is_branching_recipe
from astrocook import settings

class RecipeDialog(QDialog):
    """
    A dynamic dialog box for configuring and running Astrocook V2 recipes.
    """
    def __init__(self, recipe_category: str, recipe_name: str,
                 session: SessionV2, parent=None):
        super().__init__(parent)
        self.recipe_category = recipe_category
        self.recipe_name = recipe_name
        self.session = session
        self.parent_window = parent # Reference to MainWindowV2

        self.input_widgets = {} # Dictionary to store input widgets {param_name: widget}

        try:
            self.schema = get_recipe_schema(recipe_category, recipe_name)
        except Exception as e:
            logging.error(f"Failed to get schema for {recipe_category}.{recipe_name}: {e}")
            QMessageBox.critical(parent, "Schema Error", f"Could not load recipe schema:\n{e}")
            self.schema = None # Mark as invalid

        self._setup_ui()

    def _setup_ui(self):
        """Creates the dialog's layout and widgets."""
        if not self.schema:
            # If schema loading failed in init, don't proceed
            # Close the dialog after the error message shown in init
            # Use QTimer to reject after the constructor finishes
            from PySide6.QtCore import QTimer
            QTimer.singleShot(0, self.reject)
            return

        self.setWindowTitle(f"Recipe: {self.schema.get('brief', self.recipe_name)}")
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(15) # Add spacing between elements

        # 1. Description (Optional)
        if self.schema.get('details'):
            desc_label = QLabel(self.schema['details'])
            desc_label.setWordWrap(True)
            main_layout.addWidget(desc_label)

        # 2. Parameters Form
        form_widget = QWidget()
        form_layout = QFormLayout(form_widget)
        form_layout.setRowWrapPolicy(QFormLayout.DontWrapRows) # Labels on left
        form_layout.setLabelAlignment(Qt.AlignLeft)
        form_layout.setSpacing(10) # Space between rows

        params = self.schema.get('params', [])
        if not params:
            form_layout.addRow(QLabel("This recipe takes no parameters."))
        else:
            for param in params:
                param_name = param['name']
                param_type = param['type']
                param_doc = param.get('doc', param_name) # Use name if no doc
                param_default = str(param.get('default', '')) # Ensure string default

                label = QLabel(f"{param_doc}:")
                label.setToolTip(param_doc) # Add tooltip to label too

                widget = None

                # --- Widget Creation based on Type ---
                if param_type == str:
                    widget = QLineEdit(param_default)
                    widget.setToolTip(param_doc)
                elif param_type == int:
                    widget = QLineEdit(param_default)
                    widget.setValidator(QIntValidator())
                    widget.setToolTip(f"{param_doc} (integer)")
                elif param_type == float:
                    widget = QLineEdit(param_default)
                    validator = QDoubleValidator()
                    validator.setLocale(QLocale.C) # Use period for decimal
                    validator.setNotation(QDoubleValidator.StandardNotation)
                    widget.setValidator(validator)
                    widget.setToolTip(f"{param_doc} (decimal)")
                elif param_type == bool:
                    # Handle boolean defaults (True/False or "True"/"False")
                    default_checked = param.get('default', False)
                    if isinstance(default_checked, str):
                         default_checked = default_checked.lower() == 'true'
                    widget = QCheckBox()
                    widget.setChecked(default_checked)
                    widget.setToolTip(param_doc)
                    # Use label text directly for CheckBox
                    label.setText("") # Clear label text
                    widget.setText(f"{param_doc}")
                # TODO: Add cases for 'filepath' (QLineEdit + QPushButton)
                # TODO: Add cases for 'choice' (QComboBox using param['choices'])
                else:
                    logging.warning(f"Unsupported parameter type '{param_type}' for '{param_name}'. Using QLineEdit.")
                    widget = QLineEdit(param_default) # Fallback
                    widget.setToolTip(param_doc)
                # ------------------------------------

                if widget:
                    self.input_widgets[param_name] = widget
                    form_layout.addRow(label, widget)

        main_layout.addWidget(form_widget)

        # 3. Standard Buttons (OK/Run, Cancel, Help)
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel | QDialogButtonBox.Help
        )
        button_box.button(QDialogButtonBox.Ok).setText("Run") # Rename OK to Run

        button_box.accepted.connect(self.accept) # Connect OK/Run to accept slot
        button_box.rejected.connect(self.reject) # Connect Cancel to reject slot
        button_box.helpRequested.connect(self._show_help) # Connect Help

        main_layout.addWidget(button_box)

    def accept(self):
        """Gathers parameters, calls the recipe, and handles the result."""
        if not self.schema:
            super().reject() # Just close if schema failed
            return

        params_to_pass = {}
        logging.debug(f"Gathering parameters for {self.recipe_category}.{self.recipe_name}")
        try:
            for param_name, widget in self.input_widgets.items():
                if isinstance(widget, QLineEdit):
                    params_to_pass[param_name] = widget.text()
                elif isinstance(widget, QCheckBox):
                    # Pass boolean state as string for V1 recipe compatibility for now
                    params_to_pass[param_name] = str(widget.isChecked())
                elif isinstance(widget, QComboBox):
                    params_to_pass[param_name] = widget.currentText() # Or currentData()
                # TODO: Handle other widget types (e.g., file path)
                logging.debug(f" - {param_name}: {params_to_pass[param_name]}")

            # --- Find and Call the V2 Recipe Method ---
            # Get the recipe class instance from the session (e.g., session.flux)
            recipe_instance = getattr(self.session, self.recipe_category, None)
            if not recipe_instance:
                 raise AttributeError(f"Session object has no attribute '{self.recipe_category}'")

            # Get the actual method from the instance
            recipe_method = getattr(recipe_instance, self.recipe_name, None)
            if not callable(recipe_method):
                 raise AttributeError(f"Recipe instance has no callable method '{self.recipe_name}'")

            logging.info(f"Running recipe: {self.recipe_category}.{self.recipe_name} with params: {params_to_pass}")
            new_session_state = recipe_method(**params_to_pass)
            # -----------------------------------------

            # --- Handle Recipe Result ---
            if new_session_state == 0:
                # V1 compatibility: 0 often means validation failure within the recipe
                logging.warning(f"Recipe {self.recipe_name} indicated failure (returned 0). Check logs.")
                # Optionally show a message box here
                # QMessageBox.warning(self, "Recipe Warning", "Recipe failed validation. Check console logs.")
                # Don't close the dialog on failure, let user correct parameters
                return # Keep dialog open

            elif isinstance(new_session_state, SessionV2):
                # Success: Update the main window's state
                logging.info(f"Recipe {self.recipe_name} completed successfully.")
                if self.parent_window and hasattr(self.parent_window, 'update_gui_session_state'):
                     # Find index of old session to replace/insert after
                     try:
                         current_index = self.parent_window.active_sessions.index(self.session)
                     except ValueError:
                         current_index = -1 # Or handle as error/append case

                     # Determine if it's a branching recipe
                     branching = is_branching_recipe(self.recipe_name)

                     self.parent_window.update_gui_session_state(
                         new_session_state,
                         original_session_index=current_index,
                         is_branching=branching
                     )
                else:
                     logging.warning("Parent window reference missing or lacks update method.")

                super().accept() # Close the dialog on success

            else:
                 # Unexpected return type
                 logging.error(f"Recipe {self.recipe_name} returned unexpected type: {type(new_session_state)}")
                 QMessageBox.critical(self, "Recipe Error", "Recipe returned an unexpected result.")
                 return # Keep dialog open

        except Exception as e:
            logging.error(f"Error running recipe {self.recipe_name}: {e}", exc_info=True) # Log traceback
            QMessageBox.critical(self, "Recipe Execution Error", f"An error occurred:\n{e}")
            # Keep dialog open on error
            return


    def reject(self):
        """Closes the dialog without running the recipe."""
        logging.debug(f"Recipe dialog {self.recipe_name} cancelled.")
        super().reject() # Call the base class reject (closes dialog)

    def _show_help(self):
        """Opens the documentation URL in a web browser."""
        if self.schema and self.schema.get('url'):
            url = self.schema['url']
            logging.info(f"Opening help URL: {url}")
            try:
                # Construct full URL if needed (assuming schema url is relative)
                # from astrocook.vars import docs_url # If needed
                # full_url = docs_url + url
                webbrowser.open(url) # Or full_url
            except Exception as e:
                logging.error(f"Could not open help URL {url}: {e}")
                QMessageBox.warning(self, "Help Error", f"Could not open help URL:\n{e}")
        else:
            logging.warning("No help URL found in schema.")
            QMessageBox.information(self, "No Help", "No documentation URL is available for this recipe.")