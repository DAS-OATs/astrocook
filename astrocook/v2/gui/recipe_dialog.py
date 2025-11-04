import logging
import webbrowser
from PySide6.QtWidgets import (
    QDialog, QFormLayout, QFrame, QGridLayout, QLabel, QLineEdit, QCheckBox,
    QPushButton, QDialogButtonBox, QWidget, QComboBox, QMessageBox, QVBoxLayout
)
from PySide6.QtGui import QIntValidator, QDoubleValidator, QFont
from PySide6.QtCore import Qt, QLocale

from ..session import SessionV2
from ..structures import HistoryLogV2  # <<< *** ADD THIS IMPORT ***
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
        
        # 1. Operators (from your code)
        self.operator_list = ["+", "-", "*", "/", "**"]
        
        # 2. Functions (now with stats)
        self.function_list = ["log", "log10", "exp", "sqrt", "abs",]
                              #"sin", "cos", "tan", "asin", "acos",
                              #"min", "max", "median", "mean", "std"] # Added mean and std
    
        # 3. Booleans (for masking)
        self.boolean_list = ['&', '|', '~', '(', ')']

        # 4. Session alias map
        self.session_alias_map = self._build_alias_map()

        try:
            self.schema = get_recipe_schema(recipe_category, recipe_name)
        except Exception as e:
            logging.error(f"Failed to get schema for {recipe_category}.{recipe_name}: {e}")
            QMessageBox.critical(parent, "Schema Error", f"Could not load recipe schema:\n{e}")
            self.schema = None # Mark as invalid

        self._setup_ui()

    def _get_available_columns(self) -> list[str]:
        """Helper to get a sorted list of available column names."""
        if not self.session or not self.session.spec:
            return []
        try:
            # Get columns from the TableAdapter's raw dictionary
            cols = list(self.session.spec.t._data_dict.keys())
            cols.sort()
            return cols
        except Exception as e:
            logging.warning(f"Could not get column list for dialog: {e}")
            return []

    def _build_alias_map(self) -> dict:
        """Generates a map of short aliases (s1, s2) to full session names."""
        alias_map = {}
        s_index = 1
        if self.parent_window and hasattr(self.parent_window, 'session_histories'):
            current_session_name = self.session.name
            for hist in self.parent_window.session_histories:
                if hist.display_name != current_session_name:
                    alias = f"s{s_index}"
                    alias_map[alias] = hist.display_name
                    s_index += 1
        return alias_map

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

                if param_default == "_current_":
                    try:
                        # Get the value from the session's data core
                        current_val = getattr(self.session.spec._data, param_name)
                        param_default = str(current_val)
                    except AttributeError:
                        logging.warning(f"Could not find current value for '{param_name}'.")
                        param_default = "0.0" # Fallback

                label = QLabel(f"{param_doc}:")
                label.setToolTip(param_doc) # Add tooltip to label too

                widget = None

                if self.recipe_name == 'resample' and param_name == 'target_session':
                    widget = QComboBox()

                    other_names = list(self.session_alias_map.values())
                    other_names.sort()

                    if not other_names:
                        # No other sessions, disable
                        widget.addItem("No other sessions open")
                        widget.setEnabled(False)
                    else:
                        widget.addItem("None")
                        widget.addItems(other_names)
                    
                    # Connect signal to check for oversampling
                    widget.currentTextChanged.connect(self._on_resample_target_changed)

                # --- Widget Creation based on Type ---
                elif param_type == str:
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

        if self.recipe_name == 'resample':
            self.oversample_warning_label = QLabel("⚠️ **Warning:** Target grid is finer. This will oversample the data.")
            self.oversample_warning_label.setStyleSheet("color: #D35400;") # Orange
            self.oversample_warning_label.setVisible(False)
            form_layout.addRow(self.oversample_warning_label)

        main_layout.addWidget(form_widget)

        # --- *** NEW USABILITY PANE *** ---
        if self.recipe_name in ("apply_expression", "mask_expression"):
            
            available_cols = self._get_available_columns()
            
            if available_cols:
                ref_container = QFrame(self)
                ref_container.setObjectName("ReferencePane")
                ref_container.setFrameShape(QFrame.StyledPanel) 
                ref_layout = QVBoxLayout(ref_container)
                ref_layout.setSpacing(5)
                ref_layout.setContentsMargins(10, 10, 10, 10)

                # --- Available Columns (Buttons - Always Show) ---
                cols_title = QLabel("<b>Available Columns:</b>")
                ref_layout.addWidget(cols_title)
                # ... (column button grid logic) ...
                col_button_grid_widget = QWidget()
                col_button_grid = QGridLayout(col_button_grid_widget)
                col_button_grid.setSpacing(5)
                col_button_grid.setContentsMargins(0, 0, 0, 0)
                
                max_cols_per_row = 5 
                for i, col_name in enumerate(available_cols):
                    row = i // max_cols_per_row
                    col = i % max_cols_per_row
                    
                    button = QPushButton(col_name)
                    button.setStyleSheet("QPushButton { text-align: center; padding: 3px 5px; }")
                    button.clicked.connect(
                        lambda checked=False, text=col_name: self._insert_text_into_expression(text)
                    )
                    col_button_grid.addWidget(button, row, col)
                
                ref_layout.addWidget(col_button_grid_widget)

                # --- *** NEW: Available Sessions Buttons *** ---
                if self.session_alias_map:
                    sess_title = QLabel("<b>Available Sessions (e.g., s1.y):</b>")
                    sess_title.setContentsMargins(0, 8, 0, 0)
                    ref_layout.addWidget(sess_title)
                    
                    sess_button_widget = QWidget()
                    sess_button_layout = QVBoxLayout(sess_button_widget)
                    sess_button_layout.setSpacing(5)
                    sess_button_layout.setContentsMargins(0, 0, 0, 0)

                    # Loop through the map: {'s1': 'full_name', 's2': 'full_name_2'}
                    for alias, full_name in self.session_alias_map.items():
                        
                        # --- *** THIS IS THE FIX *** ---
                        # Button text is now "full_name (s1)"
                        button_text = f"{full_name} ({alias})"
                        button = QPushButton(button_text) 
                        # --- *** END FIX *** ---

                        button.setStyleSheet("QPushButton { text-align: left; padding: 4px; }")
                        button.setToolTip(f"Inserts '{alias}.' into expression") 
                        
                        # Action still sends the ALIAS ('s1')
                        button.clicked.connect(
                            lambda checked=False, text=alias: self._insert_text_into_expression(text)
                        )
                        sess_button_layout.addWidget(button)
                    
                    ref_layout.addWidget(sess_button_widget)

                # --- Operators (Buttons - Always Show for both) ---
                op_title = QLabel("<b>Operators:</b>")
                op_title.setContentsMargins(0, 8, 0, 0)
                ref_layout.addWidget(op_title)

                op_button_grid_widget = QWidget()
                op_button_grid = QGridLayout(op_button_grid_widget)
                op_button_grid.setSpacing(5)
                op_button_grid.setContentsMargins(0, 0, 0, 0)
                
                for i, op_name in enumerate(self.operator_list):
                    button = QPushButton(op_name)
                    button.setStyleSheet("QPushButton { text-align: center; padding: 3px 5px; }")
                    button.clicked.connect(
                        lambda checked=False, text=op_name: self._insert_text_into_expression(text)
                    )
                    op_button_grid.addWidget(button, 0, i) # Single row

                ref_layout.addWidget(op_button_grid_widget)

                # --- *** START CONDITIONAL LOGIC *** ---
                if self.recipe_name == "apply_expression":
                    # --- Common Functions (Buttons) ---
                    funcs_title = QLabel("<b>Functions:</b>")
                    funcs_title.setContentsMargins(0, 8, 0, 0) 
                    ref_layout.addWidget(funcs_title)
    
                    func_button_grid_widget = QWidget()
                    func_button_grid = QGridLayout(func_button_grid_widget)
                    func_button_grid.setSpacing(5)
                    func_button_grid.setContentsMargins(0, 0, 0, 0)
    
                    for i, func_name in enumerate(self.function_list):
                        row = i // max_cols_per_row
                        col = i % max_cols_per_row
                        
                        button = QPushButton(func_name)
                        button.setStyleSheet("QPushButton { text-align: center; padding: 3px 5px; }")
                        button.clicked.connect(
                            lambda checked=False, text=func_name: self._insert_text_into_expression(text)
                        )
                        func_button_grid.addWidget(button, row, col)
    
                    ref_layout.addWidget(func_button_grid_widget)

                elif self.recipe_name == "mask_expression":
                    # --- Boolean Operators (Buttons) ---
                    bool_title = QLabel("<b>Boolean Operators:</b>")
                    bool_title.setContentsMargins(0, 8, 0, 0)
                    ref_layout.addWidget(bool_title)
                    
                    bool_button_grid_widget = QWidget()
                    bool_button_grid = QGridLayout(bool_button_grid_widget)
                    bool_button_grid.setSpacing(5)
                    bool_button_grid.setContentsMargins(0, 0, 0, 0)

                    for i, bool_name in enumerate(self.boolean_list):
                        button = QPushButton(bool_name)
                        button.setStyleSheet("QPushButton { text-align: center; padding: 3px 5px; }")
                        button.clicked.connect(
                            lambda checked=False, text=bool_name: self._insert_text_into_expression(text)
                        )
                        bool_button_grid.addWidget(button, 0, i) # Single row

                    ref_layout.addWidget(bool_button_grid_widget)

                # --- *** END CONDITIONAL LOGIC *** ---

                main_layout.addWidget(ref_container)

        # 3. Standard Buttons (OK/Run, Cancel, Help)
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel | QDialogButtonBox.Help
        )
        button_box.button(QDialogButtonBox.Ok).setText("Run") # Rename OK to Run

        button_box.accepted.connect(self.accept) # Connect OK/Run to accept slot
        button_box.rejected.connect(self.reject) # Connect Cancel to reject slot
        button_box.helpRequested.connect(self._show_help) # Connect Help

        main_layout.addWidget(button_box)

    def _on_resample_target_changed(self, target_name: str):
        """Slot to check for oversampling when QComboBox changes."""
        if not self.oversample_warning_label:
            return

        if target_name == "None" or target_name == "":
            self.oversample_warning_label.setVisible(False)
            return
            
        try:
            # 1. Find the target session
            target_hist = None
            for hist in self.parent_window.session_histories:
                if hist.display_name == target_name:
                    target_hist = hist
                    break
            
            if not target_hist:
                self.oversample_warning_label.setVisible(False)
                return

            # 2. Compare data points
            current_points = len(self.session.spec.x)
            target_points = len(target_hist.current_state.spec.x)

            if target_points > current_points:
                self.oversample_warning_label.setVisible(True)
            else:
                self.oversample_warning_label.setVisible(False)

        except Exception as e:
            logging.warning(f"Could not check for oversampling: {e}")
            self.oversample_warning_label.setVisible(False)

    def _insert_text_into_expression(self, text: str):
        """
        Slot to insert text into the 'expression' QLineEdit.
        Intelligently adds parentheses for functions.
        """
        
        # 1. Find the QLineEdit widget for 'expression'
        target_widget = None
        if self.recipe_name == "apply_expression":
            target_widget = self.input_widgets.get('expression')
        elif self.recipe_name == "mask_expression":
            target_widget = self.input_widgets.get('expression')
        
        if not (target_widget and isinstance(target_widget, QLineEdit)):
            logging.warning(f"Could not find 'expression' QLineEdit to insert text.")
            return
            
        if text.startswith("_") or text in self.function_list:
            # It's a function: add () and move cursor
            target_widget.insert(f" {text}() ")
            target_widget.cursorBackward(False, 2) 
        elif text in self.session_alias_map:
            # It's a session alias: add . and wait
            target_widget.insert(f" {text}. ")
        elif text in self.operator_list or text in self.boolean_list:
            # It's an operator
            target_widget.insert(f" {text} ") 
        else:
            # It's a column name
            target_widget.insert(f" {text} ") 
            
        target_widget.setFocus()

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
        except Exception as e:
            logging.error(f"Error gathering parameters: {e}")
            QMessageBox.critical(self, "Parameter Error", f"Could not read parameters:\n{e}")
            return # Keep dialog open

        if self.recipe_name in ("apply_expression", "mask_expression"):
            params_to_pass['alias_map'] = self.session_alias_map

        try:
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
                    # Find the index of the SessionHistory that was active when dialog opened
                    original_history_index = -1
                    try:
                        # Find the history manager whose *current state* matches the session
                        # the dialog was launched with. This is safer than relying on active_history directly.
                        for i, history in enumerate(self.parent_window.session_histories):
                             # Compare object identity
                             if history.current_state is self.session:
                                 original_history_index = i
                                 break
                        if original_history_index == -1:
                             # Fallback or error if the original session wasn't found
                             logging.error("Could not find original session history manager!")
                             # Maybe default to current active_history index?
                             original_history_index = self.parent_window.session_histories.index(self.parent_window.active_history)

                    except (ValueError, AttributeError, IndexError) as find_err:
                        logging.error(f"Error finding original history index: {find_err}. Using active history.")
                        try: # Fallback
                             original_history_index = self.parent_window.session_histories.index(self.parent_window.active_history)
                        except (ValueError, AttributeError):
                             original_history_index = -1 # Should not happen

                    if original_history_index != -1:
                        try:
                            active_log_manager = self.parent_window.session_histories[original_history_index].log_manager
                            if isinstance(active_log_manager, HistoryLogV2):
                                active_log_manager.add_entry(
                                    recipe_name=self.recipe_name, 
                                    params=params_to_pass
                                )
                                logging.debug(f"Logged successful recipe: {self.recipe_name}")
                        except Exception as e:
                            logging.error(f"Failed to log successful recipe: {e}")

                    if original_history_index != -1:
                        branching = is_branching_recipe(self.recipe_name)
                        # Call MainWindowV2's update method
                        self.parent_window.update_gui_session_state(
                            new_session_state,
                            original_session_index=original_history_index, # Pass index in session_histories
                            is_branching=branching
                        )
                    else:
                         logging.error("Could not trigger GUI update: Original history index not found.")

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