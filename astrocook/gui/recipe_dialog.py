import logging
import webbrowser
from PySide6.QtWidgets import (
    QDialog, QFormLayout, QFrame, QGridLayout, QLabel, QLineEdit, QCheckBox,
    QPushButton, QDialogButtonBox, QWidget, QComboBox, QMessageBox, QVBoxLayout, QApplication
)
from PySide6.QtGui import QIntValidator, QDoubleValidator, QFont
from PySide6.QtCore import Qt, QLocale, Signal

from astrocook.core.session import SessionV2
from astrocook.core.structures import HistoryLogV2 
from astrocook.core.utils import get_recipe_schema, is_branching_recipe
from astrocook import settings

class RecipeDialog(QDialog):
    """
    A dynamic dialog box for configuring and running Astrocook V2 recipes.
    """

    # Emits the parameters needed for the worker thread
    recipe_requested = Signal(str, str, dict, dict) # category, recipe_name, params, alias_map

    def __init__(self, recipe_category: str, recipe_name: str,
                 session: SessionV2, parent=None):
        super().__init__(parent)
        self.recipe_category = recipe_category
        self.recipe_name = recipe_name
        self.session = session
        self.parent_window = parent 

        self.input_widgets = {} 
                
        # 1. Operators 
        self.operator_list = ["+", "-", "*", "/", "**"]
        
        # 2. Functions 
        self.function_list = ["log", "log10", "exp", "sqrt", "abs",]
    
        # 3. Booleans 
        self.boolean_list = ['&', '|', '~', '(', ')']

        # 4. Session alias map
        self.session_alias_map = self._build_alias_map()

        try:
            self.schema = get_recipe_schema(recipe_category, recipe_name)
        except Exception as e:
            logging.error(f"Failed to get schema for {recipe_category}.{recipe_name}: {e}")
            QMessageBox.critical(parent, "Schema Error", f"Could not load recipe schema:\n{e}")
            self.schema = None 

        # [CHANGE] Apply Rounded Style Globally to this Dialog
        self._apply_styles()
        
        self._setup_ui()

    def _apply_styles(self):
        """Applies consistent rounded styling to inputs."""
        pal = QApplication.palette()
        text_col = pal.color(pal.ColorRole.Text).name()
        base_col = pal.color(pal.ColorRole.Base).name()
        border_col = "#aaa" # Neutral gray

        self.setStyleSheet(f"""
            QLineEdit, QComboBox {{
                padding: 4px;
                border-radius: 4px;
                background-color: {base_col};
                color: {text_col};
            }}
            QLineEdit:focus, QComboBox:focus {{
                border: 1px solid #296bff; /* Highlight color */
            }}
            /* Slight padding for the form layout labels */
            QLabel {{ margin-right: 5px; }}
        """)

    def _get_available_columns(self) -> list[str]:
        # ... (Method unchanged) ...
        if not self.session or not self.session.spec:
            return []
        try:
            cols = list(self.session.spec.t._data_dict.keys())
            cols.sort()
            return cols
        except Exception as e:
            return []

    def _build_alias_map(self) -> dict:
        # ... (Method unchanged) ...
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
            from PySide6.QtCore import QTimer
            QTimer.singleShot(0, self.reject)
            return

        self.setWindowTitle(f"Recipe: {self.schema.get('brief', self.recipe_name).rstrip('.')}")
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(15) 

        # 1. Description (Optional)
        if self.schema.get('details'):
            desc_label = QLabel(self.schema['details'])
            desc_label.setWordWrap(True)
            main_layout.addWidget(desc_label)

        # 2. Parameters Form
        form_widget = QWidget()
        form_layout = QFormLayout(form_widget)
        form_layout.setRowWrapPolicy(QFormLayout.DontWrapRows) 
        form_layout.setLabelAlignment(Qt.AlignLeft)
        form_layout.setSpacing(5) 

        params = self.schema.get('params', [])
        if not params:
            form_layout.addRow(QLabel("This recipe takes no parameters."))
        else:
            for param in params:
                if param.get('gui_hidden', False):
                    continue
                
                param_name = param['name']
                param_type = param['type']
                param_doc = param.get('doc', param_name) 
                param_default = str(param.get('default', '')) 

                if param_default == "_current_":
                    try:
                        current_val = getattr(self.session.spec._data, param_name)
                        param_default = str(current_val)
                    except AttributeError:
                        param_default = "0.0" 

                label = QLabel(f"{param_doc}.strip('.'):")
                label.setToolTip(param_doc) 

                widget = None

                if self.recipe_name == 'resample' and param_name == 'target_session':
                    widget = QComboBox()
                    other_names = list(self.session_alias_map.values())
                    other_names.sort()

                    if not other_names:
                        widget.addItem("No other sessions open")
                        widget.setEnabled(False)
                    else:
                        widget.addItem("None")
                        widget.addItems(other_names)
                    
                    widget.currentTextChanged.connect(self._on_resample_target_changed)

                elif param_type == str:
                    widget = QLineEdit(param_default)
                    widget.setToolTip(param_doc)
                    widget.setMinimumWidth(200)
                elif param_type == int:
                    widget = QLineEdit(param_default)
                    widget.setValidator(QIntValidator())
                    widget.setToolTip(f"{param_doc} (integer)")
                    widget.setMinimumWidth(200)
                elif param_type == float:
                    widget = QLineEdit(param_default)
                    validator = QDoubleValidator()
                    validator.setLocale(QLocale.C) 
                    validator.setNotation(QDoubleValidator.StandardNotation)
                    widget.setValidator(validator)
                    widget.setToolTip(f"{param_doc} (decimal)")
                    widget.setMinimumWidth(200)
                elif param_type == bool:
                    default_checked = param.get('default', False)
                    if isinstance(default_checked, str):
                         default_checked = default_checked.lower() == 'true'
                    widget = QCheckBox()
                    widget.setChecked(default_checked)
                    widget.setToolTip(param_doc)
                    label.setText("") 
                    widget.setText(f"{param_doc}")
                else:
                    widget = QLineEdit(param_default) 
                    widget.setToolTip(param_doc)

                if widget:
                    self.input_widgets[param_name] = widget
                    form_layout.addRow(label, widget)

        if self.recipe_name == 'resample':
            self.oversample_warning_label = QLabel("⚠️ **Warning:** Target grid is finer. This will oversample the data.")
            self.oversample_warning_label.setStyleSheet("color: #D35400;") 
            self.oversample_warning_label.setVisible(False)
            form_layout.addRow(self.oversample_warning_label)

        main_layout.addWidget(form_widget)

        # --- *** NEW USABILITY PANE *** ---
        # (This section remains unchanged from your provided file)
        if self.recipe_name in ("apply_expression", "mask_expression", "split"):
            # ... (Existing logic for helper buttons) ...
            # I am keeping this section brief in the snippet to save space, 
            # but assume the exact logic from your uploaded file is preserved here.
            self._setup_usability_pane(main_layout) 

        # 3. Standard Buttons
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel | QDialogButtonBox.Help
        )
        button_box.button(QDialogButtonBox.Ok).setText("Run") 

        button_box.accepted.connect(self.accept) 
        button_box.rejected.connect(self.reject) 
        button_box.helpRequested.connect(self._show_help) 

        main_layout.addWidget(button_box)
        
    def _setup_usability_pane(self, main_layout):
        """Helper to set up the buttons for expressions."""
        # (This is the code block from your file, extracted for cleanliness)
        available_cols = self._get_available_columns()
        if not available_cols: return

        ref_container = QFrame(self)
        ref_container.setObjectName("ReferencePane")
        ref_container.setFrameShape(QFrame.StyledPanel) 
        ref_layout = QVBoxLayout(ref_container)
        ref_layout.setSpacing(5)
        ref_layout.setContentsMargins(10, 10, 10, 10)

        cols_title = QLabel("<b>Available Columns:</b>")
        ref_layout.addWidget(cols_title)
        
        col_button_grid_widget = QWidget()
        col_button_grid = QGridLayout(col_button_grid_widget)
        col_button_grid.setSpacing(5)
        col_button_grid.setContentsMargins(0, 0, 0, 0)
        
        max_cols_per_row = 5 
        for i, col_name in enumerate(available_cols):
            row = i // max_cols_per_row
            col = i % max_cols_per_row
            button = QPushButton(col_name)
            # Remove styling that conflicts with general stylesheet or refine it
            button.setStyleSheet("QPushButton { text-align: center; padding: 3px 5px; }") 
            button.clicked.connect(lambda checked=False, text=col_name: self._insert_text_into_expression(text))
            col_button_grid.addWidget(button, row, col)
        
        ref_layout.addWidget(col_button_grid_widget)

        if self.session_alias_map:
            sess_title = QLabel("<b>Available Sessions (e.g., s1.y):</b>")
            sess_title.setContentsMargins(0, 8, 0, 0)
            ref_layout.addWidget(sess_title)
            
            sess_button_widget = QWidget()
            sess_button_layout = QVBoxLayout(sess_button_widget)
            sess_button_layout.setSpacing(5)
            sess_button_layout.setContentsMargins(0, 0, 0, 0)

            for alias, full_name in self.session_alias_map.items():
                button_text = f"{full_name} ({alias})"
                button = QPushButton(button_text) 
                button.setStyleSheet("QPushButton { text-align: left; padding: 4px; }")
                button.setToolTip(f"Inserts '{alias}.' into expression") 
                button.clicked.connect(lambda checked=False, text=alias: self._insert_text_into_expression(text))
                sess_button_layout.addWidget(button)
            
            ref_layout.addWidget(sess_button_widget)

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
            button.clicked.connect(lambda checked=False, text=op_name: self._insert_text_into_expression(text))
            op_button_grid.addWidget(button, 0, i)

        ref_layout.addWidget(op_button_grid_widget)

        if self.recipe_name == "apply_expression":
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
                button.clicked.connect(lambda checked=False, text=func_name: self._insert_text_into_expression(text))
                func_button_grid.addWidget(button, row, col)

            ref_layout.addWidget(func_button_grid_widget)

        elif self.recipe_name in ("mask_expression", "split"):
            bool_title = QLabel("<b>Boolean Operators:</b>")
            bool_title.setContentsMargins(0, 8, 0, 0)
            ref_layout.addWidget(bool_title)
            
            bool_button_grid_widget = QWidget()
            bool_button_grid = QGridLayout(bool_button_grid_widget)
            bool_button_grid.setSpacing(5)
            bool_button_grid.setContentsMargins(0, 0, 0, 0)

            for i, bool_name in enumerate(self.boolean_list):
                display_text = bool_name.replace('&', '&&')
                button = QPushButton(display_text)
                button.setStyleSheet("QPushButton { text-align: center; padding: 3px 5px; }")
                button.clicked.connect(lambda checked=False, text=bool_name: self._insert_text_into_expression(text))
                bool_button_grid.addWidget(button, 0, i)

            ref_layout.addWidget(bool_button_grid_widget)

        main_layout.addWidget(ref_container)

    def _on_resample_target_changed(self, target_name: str):
        # ... (Method unchanged) ...
        if not self.oversample_warning_label: return
        if target_name == "None" or target_name == "":
            self.oversample_warning_label.setVisible(False)
            return
        try:
            target_hist = None
            for hist in self.parent_window.session_histories:
                if hist.display_name == target_name:
                    target_hist = hist
                    break
            if not target_hist:
                self.oversample_warning_label.setVisible(False)
                return
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
        # ... (Method unchanged) ...
        target_widget = None
        if self.recipe_name == "apply_expression":
            target_widget = self.input_widgets.get('expression')
        elif self.recipe_name == "mask_expression":
            target_widget = self.input_widgets.get('expression')
        elif self.recipe_name in ("mask_expression", "split"):
            target_widget = self.input_widgets.get('expression')
        
        if not (target_widget and isinstance(target_widget, QLineEdit)):
            return
            
        if text.startswith("_") or text in self.function_list:
            target_widget.insert(f" {text}() ")
            target_widget.cursorBackward(False, 2) 
        elif text in self.session_alias_map:
            target_widget.insert(f" {text}. ")
        elif text in self.operator_list or text in self.boolean_list:
            target_widget.insert(f" {text} ") 
        else:
            target_widget.insert(f" {text} ") 
        target_widget.setFocus()

    def accept(self):
        # ... (Method unchanged) ...
        if not self.schema:
            super().reject(); return

        params_to_pass = {}
        try:
            for param_name, widget in self.input_widgets.items():
                if isinstance(widget, QLineEdit):
                    params_to_pass[param_name] = widget.text()
                elif isinstance(widget, QCheckBox):
                    params_to_pass[param_name] = str(widget.isChecked())
                elif isinstance(widget, QComboBox):
                    params_to_pass[param_name] = widget.currentText() 
        except Exception as e:
            logging.error(f"Error gathering parameters: {e}")
            QMessageBox.critical(self, "Parameter Error", f"Could not read parameters:\n{e}")
            return 

        self.recipe_requested.emit(
            self.recipe_category,
            self.recipe_name,
            params_to_pass,
            self.session_alias_map
        )
        super().accept()

    def reject(self):
        super().reject() 

    def _show_help(self):
        # ... (Method unchanged) ...
        if self.schema and self.schema.get('url'):
            url = self.schema['url']
            try:
                webbrowser.open(url) 
            except Exception as e:
                QMessageBox.warning(self, "Help Error", f"Could not open help URL:\n{e}")
        else:
            QMessageBox.information(self, "No Help", "No documentation URL is available for this recipe.")