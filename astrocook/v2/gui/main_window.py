# astrocook/v2/gui/main_window.py

import logging
import os
from PySide6.QtCore import ( # <<< Modify this import
    Qt,
    QItemSelectionModel, # <<< Add QItemSelectionModel
    QLocale,
    QPropertyAnimation, QEasingCurve, QRect, QPoint,
    QParallelAnimationGroup,
    QSize, QStringListModel, QTimer
)
from PySide6.QtGui import QAction, QDoubleValidator, QKeySequence 
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QDialog, QFileDialog, 
    QMainWindow, QWidget, QVBoxLayout, QFormLayout, QLabel, QLineEdit, QListView, QMessageBox,
    QPushButton, QSizePolicy, QSpacerItem, QStackedWidget, QStyle
)

from .pyside_plot import SpectrumPlotWidget
from .recipe_dialog import RecipeDialog
from ..session import SessionV2, load_session_from_file
try:
    from ...v1.functions import trans_parse
    from ...v1.vars import xem_d
    V1_FUNCTIONS_AVAILABLE = True
except ImportError:
    logging.error("Could not import V1 functions (trans_parse, xem_d) needed for redshift cursor.")
    V1_FUNCTIONS_AVAILABLE = False

# --- Constants for Sidebar Widths ---
LEFT_SIDEBAR_WIDTH = 300
RIGHT_SIDEBAR_WIDTH = 200
ANIMATION_DURATION = 150 # ** Speed up animation **
BUTTON_WIDTH = 20
BUTTON_HEIGHT = 30

class MainWindowV2(QMainWindow):
    def __init__(self, session):
        super().__init__()
        self.active_sessions = []
        self.session_manager = session
        self.session_model = QStringListModel()

        # ** Initialize animation attributes to None **
        self.left_animation_group = None
        self.right_animation_group = None

        # ** Add History Index **
        self.history_index = -1 # Index of the currently active session in active_sessions

        self.setGeometry(100, 100, 450, 150) # Initial small size
        screen_geometry = QApplication.primaryScreen().geometry()
        x = (screen_geometry.width() - self.width()) // 2
        y = (screen_geometry.height() - self.height()) // 2
        self.move(x, y)

        # --- ** Central Widget is NOW the Stack ** ---
        self.central_stack = QStackedWidget()
        self._setup_plot_view()
        self._setup_empty_view()
        self.setCentralWidget(self.central_stack) # Stack fills the window initially

        # --- ** Sidebars are Children of the MainWindow, floating above ** ---
        self._setup_left_sidebar()  # Creates self.left_sidebar_widget
        self._setup_right_sidebar() # Creates self.right_sidebar_widget

        # --- ** Buttons are also Children of MainWindow ** ---
        self._setup_collapse_buttons()

        self._create_menubar()
        self._apply_styles()

        # --- Initial State ---
        is_initial_session_valid = bool(self.session_manager.spec and len(self.session_manager.spec.x) > 0)
        # Add the initial session (if valid) *without* triggering history truncation
        if is_initial_session_valid:
            self._add_session_internal(session, is_initial=True) # Use internal add
        else:
            self._update_ui_state(False, is_startup=True) # Show empty state

        self._update_undo_redo_actions() # Set initial state

    def resizeEvent(self, event):
        """Handle window resizing to reposition floating elements."""
        super().resizeEvent(event)
        self._reposition_floating_widgets()

    def _reposition_floating_widgets(self):
        """Positions sidebars and buttons based on current window size and visibility."""
        window_height = self.height()
        window_width = self.width()
        menubar_height = self.menuBar().height() if self.menuBar() else 0
        content_y_start = menubar_height
        content_height = window_height - content_y_start

        # ** Calculate centered Y position for buttons **
        button_y = content_y_start + (content_height - BUTTON_HEIGHT) // 2

        # --- Left Sidebar & Button ---
        left_width = self.left_sidebar_widget.width()
        button_x = left_width - BUTTON_WIDTH if left_width > 0 else left_width
        # ** Use calculated button_y **
        self.session_collapse_button.setGeometry(button_x, button_y, BUTTON_WIDTH, BUTTON_HEIGHT)
        # Sidebar still fills the full content height
        self.left_sidebar_widget.setGeometry(0, content_y_start, left_width, content_height)

        # --- Right Sidebar & Button ---
        right_width = self.right_sidebar_widget.width()
        sidebar_x = window_width - right_width
        button_x = sidebar_x if right_width > 0 else sidebar_x - BUTTON_WIDTH
         # ** Use calculated button_y **
        self.plot_controls_collapse_button.setGeometry(button_x, button_y, BUTTON_WIDTH, BUTTON_HEIGHT)
        # Sidebar still fills the full content height
        self.right_sidebar_widget.setGeometry(sidebar_x, content_y_start, right_width, content_height)
        
    def _setup_plot_view(self):
        self.plot_viewer = SpectrumPlotWidget(self.session_manager, self)
        self.central_stack.addWidget(self.plot_viewer)

    def _setup_empty_view(self):
        empty_widget = QWidget()
        layout = QVBoxLayout(empty_widget)
        label = QLabel("Welcome to Astrocook v2.\n\nUse 'File > Open Spectrum...'")
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        font = label.font() # Get the current font
        font.setPointSize(18) # Set a specific size: font.setPointSize(14)
        label.setFont(font) # Apply the modified font

        layout.addWidget(label)
        self.central_stack.addWidget(empty_widget)

    def _setup_left_sidebar(self):
        """Creates the left sidebar widget as a child of the main window."""
        self.left_sidebar_widget = QWidget(self) # ** Parent is main window **
        sidebar_layout = QVBoxLayout(self.left_sidebar_widget)
        sidebar_layout.setContentsMargins(0, 0, 0, 0); sidebar_layout.setSpacing(0)

        self.session_list_view = QListView()
        # ... (setup model, font, connection) ...
        self.session_list_view.setModel(self.session_model)
        font = self.session_list_view.font(); font.setPointSize(14); self.session_list_view.setFont(font)
        self.session_list_view.clicked.connect(self._on_session_switched)

        sidebar_layout.addWidget(self.session_list_view)
        self.left_sidebar_widget.setObjectName("SessionContainer")
        self.session_list_view.setObjectName("SessionListView")

        # Initial geometry set later, start hidden/collapsed
        self.left_sidebar_widget.resize(0, self.height())
        self.left_sidebar_widget.setVisible(False)

    def _setup_right_sidebar(self):
        """Creates the right sidebar widget with plot toggles and cursor controls."""
        self.right_sidebar_widget = QWidget(self)
        sidebar_layout = QVBoxLayout(self.right_sidebar_widget)
        sidebar_layout.setContentsMargins(10, 10, 10, 10)
        sidebar_layout.setSpacing(15) # Increase spacing between sections

        # --- Plot Element Toggles ---
        plot_toggles_layout = QVBoxLayout() # Group toggles
        plot_toggles_layout.setSpacing(8)
        plot_toggles_layout.addWidget(QLabel("<b>Plot Elements:</b>")) # Section title

        self.error_checkbox = QCheckBox("1-sigma error"); self.error_checkbox.setChecked(True)
        self.error_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.error_checkbox)

        self.continuum_checkbox = QCheckBox("Continuum"); self.continuum_checkbox.setChecked(True)
        self.continuum_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.continuum_checkbox)

        self.model_checkbox = QCheckBox("Absorption model"); self.model_checkbox.setChecked(True)
        self.model_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.model_checkbox)

        self.systems_checkbox = QCheckBox("Systems"); self.systems_checkbox.setChecked(True) # Often useful default
        self.systems_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.systems_checkbox)

        # self.lines_checkbox = QCheckBox("Lines"); self.lines_checkbox.setChecked(False) # Add if needed
        # self.lines_checkbox.stateChanged.connect(self._trigger_replot)
        # plot_toggles_layout.addWidget(self.lines_checkbox)

        sidebar_layout.addLayout(plot_toggles_layout) # Add group to main layout
    
        # --- Redshift Cursor Controls ---
        cursor_layout = QVBoxLayout() # Group cursor controls
        cursor_layout.setSpacing(8)
        cursor_layout.addWidget(QLabel("<b>Redshift Cursor:</b>")) # Section title

        # Use QFormLayout for label-input pairs
        form_layout = QFormLayout()
        form_layout.setRowWrapPolicy(QFormLayout.DontWrapRows)
        form_layout.setLabelAlignment(Qt.AlignLeft) # Align labels left
        form_layout.setSpacing(5)

        self.cursor_series_input = QLineEdit("Ly_a") # Default series
        self.cursor_z_input = QLineEdit("0.0") # Default redshift
        # Add validator for redshift input
        z_validator = QDoubleValidator()
        # ** Set locale to one using period (e.g., C locale or en_US) **
        z_validator.setLocale(QLocale.C) # Force C locale (period decimal separator)
        z_validator.setNotation(QDoubleValidator.StandardNotation)
        self.cursor_z_input.setValidator(z_validator)

        form_layout.addRow("Series:", self.cursor_series_input)
        form_layout.addRow("Redshift:", self.cursor_z_input)
        cursor_layout.addLayout(form_layout) # Add form to cursor group

        self.cursor_show_checkbox = QCheckBox("Show Cursor Lines"); self.cursor_show_checkbox.setChecked(False)
        cursor_layout.addWidget(self.cursor_show_checkbox)

        sidebar_layout.addLayout(cursor_layout) # Add group to main layout

        # Connect signals for cursor updates
        self.cursor_series_input.editingFinished.connect(self._update_cursor_and_replot)
        self.cursor_z_input.editingFinished.connect(self._update_cursor_and_replot)
        self.cursor_show_checkbox.stateChanged.connect(self._trigger_replot) # Just replot on show/hide

        # --- Spacer ---
        spacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        sidebar_layout.addItem(spacer)

        self.right_sidebar_widget.setObjectName("PlotControlsContainer")
        # Assign object names for styling if needed
        self.error_checkbox.setObjectName("PlotControlCheckbox")
        self.continuum_checkbox.setObjectName("PlotControlCheckbox")
        self.model_checkbox.setObjectName("PlotControlCheckbox")
        self.systems_checkbox.setObjectName("PlotControlCheckbox")
        self.cursor_show_checkbox.setObjectName("PlotControlCheckbox")

        self.right_sidebar_widget.setVisible(False)

    def _setup_collapse_buttons(self):
        """Creates and positions the collapse buttons as children of main window."""
        # Left Button
        self.session_collapse_button = QPushButton(self) # ** Parent is main window **
        self.session_collapse_button.setToolTip("Collapse/Expand Session Panel")
        self.session_collapse_button.setFixedSize(QSize(BUTTON_WIDTH, BUTTON_HEIGHT))
        self.session_collapse_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed) # Fill height
        self.session_collapse_button.clicked.connect(self._toggle_left_sidebar)
        self.session_collapse_button.setObjectName("CollapseButton")
        self.session_collapse_button.setVisible(False) # Start hidden

        # Right Button
        self.plot_controls_collapse_button = QPushButton(self) # ** Parent is main window **
        self.plot_controls_collapse_button.setToolTip("Collapse/Expand Plot Controls")
        self.plot_controls_collapse_button.setFixedSize(QSize(BUTTON_WIDTH, BUTTON_HEIGHT))
        self.plot_controls_collapse_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed) # Fill height
        self.plot_controls_collapse_button.clicked.connect(self._toggle_right_sidebar)
        self.plot_controls_collapse_button.setObjectName("CollapseButton")
        self.plot_controls_collapse_button.setVisible(False) # Start hidden

    def _trigger_replot(self):
        """Slot to be called when any plot toggle checkbox changes state."""
        if self.plot_viewer: # Check if plot viewer exists
            # Call the plot function, which will now read checkbox states
            self.plot_viewer.plot_spectrum()

    def _update_cursor_and_replot(self):
        """Slot called when cursor series or redshift changes."""
        # Optional: Add validation for series input here if needed
        # (e.g., check if series exists in xem_d)
        if self.cursor_show_checkbox.isChecked(): # Only replot if cursor is visible
            self._trigger_replot()

    def _create_menubar(self):
        menu_bar = self.menuBar()
        
        # --- Define all 6 Scientific Menus ---
        
        # 1. Standard Application Menus
        file_menu = menu_bar.addMenu("&File")
        view_menu = menu_bar.addMenu("&View")

        # 2. Scientific Interaction Menus
        edit_menu = menu_bar.addMenu("&Edit")
        flux_menu = menu_bar.addMenu("&Flux")
        continuum_menu = menu_bar.addMenu("&Continuum")
        absorbers_menu = menu_bar.addMenu("&Absorbers")
        
        # --- File and View operations ---

        # Open file
        open_action = QAction("&Open Spectrum...", self); open_action.setShortcut("Ctrl+O"); open_action.triggered.connect(self._on_open_spectrum)
        file_menu.addAction(open_action)
        save_action = QAction("&Save Session...", self); save_action.setShortcut("Ctrl+S"); save_action.triggered.connect(self._on_save_session); 
        save_action.setEnabled(False)
        file_menu.addAction(save_action)

        # --- View Menu Actions ---
        toggle_left_action = QAction("Toggle Session Panel", self); toggle_left_action.triggered.connect(self._toggle_left_sidebar)
        view_menu.addAction(toggle_left_action)
        toggle_right_action = QAction("Toggle Plot Controls", self); toggle_right_action.triggered.connect(self._toggle_right_sidebar)
        view_menu.addAction(toggle_right_action)

        self.toggle_left_action = toggle_left_action
        self.toggle_right_action = toggle_right_action

        # --- Add QActions for V2 Recipes ---
        
        # ** Undo Action **
        self.undo_action = QAction("&Undo", self)
        self.undo_action.setShortcut(QKeySequence.Undo) # Standard shortcut (Ctrl+Z / Cmd+Z)
        self.undo_action.triggered.connect(self._undo_last_action)
        edit_menu.addAction(self.undo_action)

        # ** Redo Action **
        self.redo_action = QAction("&Redo", self)
        self.redo_action.setShortcut(QKeySequence.Redo) # Standard shortcut (Ctrl+Y / Cmd+Shift+Z)
        self.redo_action.triggered.connect(self._redo_last_action)
        edit_menu.addAction(self.redo_action)

        edit_menu.addSeparator()

        # RECIPES FOR 'EDIT' MENU (x_convert, y_convert)
        
        # x_convert Action
        x_convert_action = QAction("Convert &X Axis...", self)
        x_convert_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "x_convert"))
        edit_menu.addAction(x_convert_action)
        self.x_convert_action = x_convert_action; x_convert_action.setEnabled(False)

        # y_convert Action
        y_convert_action = QAction("Convert &Y Axis...", self)
        y_convert_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "y_convert"))
        edit_menu.addAction(y_convert_action)
        self.y_convert_action = y_convert_action; y_convert_action.setEnabled(False)

        # RECIPES FOR 'FLUX' MENU (rebin)
        
        # rebin Action
        rebin_action = QAction("&Rebin Spectrum...", self)
        rebin_action.triggered.connect(lambda: self._launch_recipe_dialog("flux", "rebin"))
        flux_menu.addAction(rebin_action)
        self.rebin_action = rebin_action; rebin_action.setEnabled(False)

        # ... (Other menu item actions will be added here later) ...

        self._update_undo_redo_actions()

    def _launch_recipe_dialog(self, category, name):
        """Creates and executes the RecipeDialog for the specified recipe."""
        if not self.session_manager or not self.session_manager.spec:
            QMessageBox.warning(self, "No Session", "Please load a spectrum before running a recipe.")
            return

        logging.info(f"Launching dialog for recipe: {category}.{name}")
        dialog = RecipeDialog(category, name, self.session_manager, self)
        # .exec() runs the dialog modally and returns QDialog.Accepted or QDialog.Rejected
        result = dialog.exec()

        if result == QDialog.Accepted:
            logging.debug(f"Recipe dialog {name} accepted (state update handled by dialog's accept method).")
            # The dialog's accept() method now calls self.update_gui_session_state
        else:
            logging.debug(f"Recipe dialog {name} cancelled.")

    def _apply_styles(self):
        # ... (Get palette colors) ...
        # ... (Define sidebar_bg, item_selected_bg, etc.) ...
        try: # Added try-except for palette
            palette = QApplication.palette()
            # Get RGB components from the desired background color
            sidebar_base_color = palette.color(palette.ColorRole.Button)
            sidebar_r = sidebar_base_color.red()
            sidebar_g = sidebar_base_color.green()
            sidebar_b = sidebar_base_color.blue()
            sidebar_a = 0.7 # Alpha value (0.0 to 1.0), e.g., 90% opaque

            sidebar_bg = palette.color(palette.ColorRole.Button).name()
            item_selected_bg = palette.color(palette.ColorRole.Highlight).name()
            item_selected_text = palette.color(palette.ColorRole.HighlightedText).name()
            #item_hover_bg = palette.color(palette.ColorRole.Highlight).name()
            border_color = palette.color(palette.ColorRole.Mid).name()
            button_fg = palette.color(palette.ColorRole.ButtonText).name()
            text_color = palette.color(palette.ColorRole.Text).name()
        except Exception as e:
            logging.warning(f"Could not query palette: {e}. Using fallback colors.")
            sidebar_r, sidebar_g, sidebar_b = 224, 224, 224 # RGB for #E0E0E0
            sidebar_a = 0.9
            sidebar_bg="#E0E0E0"; item_selected_bg="#808080"; item_selected_text="#FFFFFF";
            item_hover_bg="#D0D0D0"; border_color="#B0B0B0"; button_fg="#000000";

        qss = f"""
            /* Sidebar Containers: Transparent background, NO side borders */
            QWidget#SessionContainer, QWidget#PlotControlsContainer {{
                background-color: rgba({sidebar_r}, {sidebar_g}, {sidebar_b}, {sidebar_a});
                border: none; /* ** Remove all borders ** */
            }}
            
            /* Make ListView background fully transparent */
            QListView#SessionListView {{ background-color: transparent; border: none; }}
            /* Ensure list items are NOT transparent */
            QListView#SessionListView::item {{
                padding: 6px 10px; border: none;
                background-color: transparent; /* Inherit sidebar bg initially */
                color: {button_fg}; /* Use button text color for readability */
             }}
            /* Selected/Hover items should be opaque */
            QListView#SessionListView::item:selected {{ background-color: {item_selected_bg}; color: {item_selected_text}; }}
            
            /* Right Sidebar Checkboxes & Labels */
            QWidget#PlotControlsContainer QLabel {{ /* Style section labels */
                margin-bottom: 4px; /* Space below label */
            }}
            QCheckBox#PlotControlCheckbox {{
                color: {text_color}; spacing: 5px; padding: 2px 0px; /* Reduced padding */
                background-color: transparent; padding-left: 15px; /* Indentation */
            }}
            /* Style LineEdits for Cursor */
            QWidget#PlotControlsContainer QLineEdit {{
                padding: 3px;
                border: 1px solid {border_color};
                border-radius: 3px;
                background-color: {palette.color(palette.ColorRole.Base).name() if 'palette' in locals() else '#FFFFFF'};
                color: {text_color};
            }}
            /* Style Form Layout labels */
            QWidget#PlotControlsContainer QFormLayout QLabel {{ /* More specific selector */
                color: {text_color};
                margin-bottom: 0px; /* Override default label margin */
                padding-top: 4px; /* Align better with LineEdit */
            }}

            /* Collapse Buttons: Blend, NO borders */
            QPushButton#CollapseButton {{
                background-color: none;
                color: {button_fg}; /* Icon color */
                border: none; /* ** Remove all borders ** */
                margin: 0px; padding: 0px;
                icon-size: 16px;
                /* Maybe add slight rounding */
                border-radius: 3px;
            }}
        """
        self.setStyleSheet(qss)
        # Set object names after creation if not done elsewhere

    # --- ** NEW Toggle & Animation Methods ** ---

    def _toggle_left_sidebar(self):
        """Animates the left sidebar overlay and its button together."""
        is_currently_visible = self.left_sidebar_widget.width() > 0 # Use width check
        start_sidebar_geom = QRect(self.left_sidebar_widget.geometry()) # Get current geometry
        end_sidebar_geom = QRect(start_sidebar_geom) # Copy for end state
        start_button_x = self.session_collapse_button.x()
        end_button_x = 0 # Default button X when closed

        if is_currently_visible: # Closing
            end_sidebar_geom.setWidth(0) # Animate width to 0
            # Button X already calculated (end_button_x = 0)
            target_visible = False
        else: # Opening
            # Start geom should reflect the hidden state (width 0)
            start_sidebar_geom.setWidth(0)
            end_sidebar_geom.setWidth(LEFT_SIDEBAR_WIDTH) # Animate width open
            end_button_x = LEFT_SIDEBAR_WIDTH - BUTTON_WIDTH # Button X moves out
            target_visible = True
            self.left_sidebar_widget.setVisible(True) # Show before animating

        # --- Animation Group ---
        # Stop previous animation if running
        if self.left_animation_group and self.left_animation_group.state() == QParallelAnimationGroup.Running:
            self.left_animation_group.stop()
        self.left_animation_group = QParallelAnimationGroup(self)

        # ** Animate sidebar GEOMETRY **
        sidebar_anim = QPropertyAnimation(self.left_sidebar_widget, b"geometry")
        sidebar_anim.setDuration(ANIMATION_DURATION)
        sidebar_anim.setStartValue(start_sidebar_geom) # Start geometry
        sidebar_anim.setEndValue(end_sidebar_geom)     # End geometry
        sidebar_anim.setEasingCurve(QEasingCurve.InOutQuad)
        self.left_animation_group.addAnimation(sidebar_anim)

        # Animate button X position
        button_anim = QPropertyAnimation(self.session_collapse_button, b"pos")
        button_anim.setDuration(ANIMATION_DURATION)
        # Use current Y, animate X
        button_anim.setStartValue(QPoint(start_button_x, self.session_collapse_button.y()))
        button_anim.setEndValue(QPoint(end_button_x, self.session_collapse_button.y()))
        button_anim.setEasingCurve(QEasingCurve.InOutQuad)
        self.left_animation_group.addAnimation(button_anim)

        # Hide sidebar widget AFTER closing animation finishes
        if not target_visible:
            # Ensure lambda captures the correct widget state if needed
            self.left_animation_group.finished.connect(lambda: self.left_sidebar_widget.setVisible(False))

        self.left_animation_group.start()
        self._update_sidebar_button_icon(self.session_collapse_button, target_visible, is_left=True)


    def _toggle_right_sidebar(self):
        # ... (get start/end width, window_width, start/end sidebar X, start/end button X) ...
        start_width = self.right_sidebar_widget.width(); end_width = 0 if start_width > 0 else RIGHT_SIDEBAR_WIDTH
        window_width = self.width()
        start_sidebar_x = self.right_sidebar_widget.x(); end_sidebar_x = window_width if end_width == 0 else window_width - end_width
        start_button_x = self.plot_controls_collapse_button.x(); end_button_x = window_width - BUTTON_WIDTH if end_width == 0 else window_width - end_width
        target_visible = end_width > 0

        # ** Calculate the constant button Y **
        button_y = self.plot_controls_collapse_button.y() # Keep current Y

        if target_visible: self.right_sidebar_widget.setVisible(True)
        self.right_animation_group = QParallelAnimationGroup(self)
        # Sidebar Animation (geometry)
        sidebar_geom_anim = QPropertyAnimation(self.right_sidebar_widget, b"geometry"); sidebar_geom_anim.setDuration(ANIMATION_DURATION)
        sidebar_geom_anim.setStartValue(QRect(start_sidebar_x, self.right_sidebar_widget.y(), start_width, self.right_sidebar_widget.height()))
        sidebar_geom_anim.setEndValue(QRect(end_sidebar_x, self.right_sidebar_widget.y(), end_width, self.right_sidebar_widget.height()))
        sidebar_geom_anim.setEasingCurve(QEasingCurve.InOutQuad)
        self.right_animation_group.addAnimation(sidebar_geom_anim)
        # Button Animation (pos)
        button_anim = QPropertyAnimation(self.plot_controls_collapse_button, b"pos"); button_anim.setDuration(ANIMATION_DURATION)
        button_anim.setStartValue(QPoint(start_button_x, button_y)) # Start at current Y
        button_anim.setEndValue(QPoint(end_button_x, button_y))     # End at same Y
        button_anim.setEasingCurve(QEasingCurve.InOutQuad)
        self.right_animation_group.addAnimation(button_anim)

        if not target_visible: self.right_animation_group.finished.connect(lambda: self.right_sidebar_widget.setVisible(False))
        self.right_animation_group.start()
        self._update_sidebar_button_icon(self.plot_controls_collapse_button, target_visible, is_left=False)

    def _update_sidebar_button_icon(self, button, sidebar_is_visible, is_left):
        """Updates collapse button icon based on sidebar visibility."""
        # Icon logic adjusted: Left button shows '>' when sidebar is hidden (points right to expand)
        # Right button shows '<' when sidebar is hidden (points left to expand)
        if is_left:
            icon_name = QStyle.SP_ArrowLeft if sidebar_is_visible else QStyle.SP_ArrowRight
            tooltip = "Collapse Session Panel" if sidebar_is_visible else "Expand Session Panel"
        else: # Right sidebar button
            icon_name = QStyle.SP_ArrowRight if sidebar_is_visible else QStyle.SP_ArrowLeft
            tooltip = "Collapse Plot Controls" if sidebar_is_visible else "Expand Plot Controls"

        icon = self.style().standardIcon(icon_name)
        button.setIcon(icon); button.setToolTip(tooltip)

    def _update_ui_state(self, is_valid_session, is_startup=False):
        is_valid_session = bool(is_valid_session)

        # ... (Enable/Disable View menu actions) ...
        if hasattr(self, 'toggle_left_action'): self.toggle_left_action.setEnabled(is_valid_session)
        if hasattr(self, 'toggle_right_action'): self.toggle_right_action.setEnabled(is_valid_session)
        # --- Update Session Logic ---

        # ** Enable/Disable Recipe Actions **
        enable_recipes = is_valid_session
        # Check if actions exist before enabling/disabling
        if hasattr(self, 'x_convert_action'): self.x_convert_action.setEnabled(enable_recipes)
        if hasattr(self, 'y_convert_action'): self.y_convert_action.setEnabled(enable_recipes)
        if hasattr(self, 'rebin_action'): self.rebin_action.setEnabled(enable_recipes)
        # ... enable/disable other recipe actions ...

        # Enable Save action only if valid session
        if hasattr(self, 'save_action'): self.save_action.setEnabled(is_valid_session)

        if is_valid_session:
            # --- State when a valid session IS loaded ---
            if self.central_stack.currentIndex() != 0: # If switching from empty
                self.central_stack.setCurrentIndex(0)
                # Resize only on first valid load
                # Check if *any* previous session existed and was valid
                was_previously_empty = not any(s and s.spec and len(s.spec.x) > 0 for s in self.active_sessions[:-1])
                if is_startup or was_previously_empty:
                    logging.debug("Resizing and centering window for first valid session.")
                    self.resize(1400, 900)
                    # Recenter after resize
                    screen_geometry = QApplication.primaryScreen().geometry()
                    x = (screen_geometry.width() - self.width()) // 2
                    y = (screen_geometry.height() - self.height()) // 2
                    self.move(x, y)

            # Show buttons
            self.session_collapse_button.setVisible(True)
            self.plot_controls_collapse_button.setVisible(True)

            # Update icons based on *current width* (more reliable during animation)
            self._update_sidebar_button_icon(self.session_collapse_button, self.left_sidebar_widget.width() > 0, is_left=True)
            self._update_sidebar_button_icon(self.plot_controls_collapse_button, self.right_sidebar_widget.width() > 0, is_left=False)

            # Ensure sidebars are visible if they have width > 0
            # (Animation handles the actual hiding/showing based on width)
            if self.left_sidebar_widget.width() > 0: self.left_sidebar_widget.setVisible(True)
            if self.right_sidebar_widget.width() > 0: self.right_sidebar_widget.setVisible(True)         

        else:
            if self.central_stack.currentIndex() != 1: self.central_stack.setCurrentIndex(1)

            # Hide buttons
            self.session_collapse_button.setVisible(False)
            self.plot_controls_collapse_button.setVisible(False)

            # Hide sidebars immediately (stop animations if running)
            if self.left_animation_group and self.left_animation_group.state() == QParallelAnimationGroup.Running: self.left_animation_group.stop()
            if self.right_animation_group and self.right_animation_group.state() == QParallelAnimationGroup.Running: self.right_animation_group.stop()
            self.left_sidebar_widget.setVisible(False); self.left_sidebar_widget.setGeometry(0, self.left_sidebar_widget.y(), 0, self.left_sidebar_widget.height()) # Reset geom
            self.right_sidebar_widget.setVisible(False); self.right_sidebar_widget.setGeometry(self.width(), self.right_sidebar_widget.y(), 0, self.right_sidebar_widget.height()) # Reset geom

        # Reposition elements after visibility/state changes might affect layout needs
        # QTimer.singleShot(0, self._reposition_floating_widgets) # Schedule reposition slightly later
        self._reposition_floating_widgets() # Try immediate reposition first

    def update_gui_session_state(self, new_session: SessionV2, original_session_index: int, is_branching: bool):
        """
        Updates the GUI state after a recipe successfully returns a new session.
        Called by RecipeDialog.accept().
        """
        if not isinstance(new_session, SessionV2):
            logging.error("update_gui_session_state received invalid session object.")
            return

        # ** 1. Truncate history beyond the current point **
        # If user did Undo, then performed an action, remove the 'redoable' states
        current_len = len(self.active_sessions)
        if self.history_index < current_len - 1:
            logging.debug(f"Truncating history: removing {current_len - 1 - self.history_index} states.")
            del self.active_sessions[self.history_index + 1:]
            # Update model (remove rows from the end)
            self.session_model.removeRows(self.history_index + 1, current_len - 1 - self.history_index)

        # ** 2. Append the new state (Works for both linear and branching after truncation) **
        insert_pos = self.history_index + 1
        logging.debug(f"Adding new session state at index {insert_pos} (Branching={is_branching})")
        self.active_sessions.insert(insert_pos, new_session)
        self.session_model.insertRow(insert_pos)
        self.session_model.setData(self.session_model.index(insert_pos), new_session.name)
        self.history_index = insert_pos # New state becomes current

        # ** 3. Update the view to show the new state **
        self._update_view_for_session(new_session, set_current_list_item=True)
        self._update_undo_redo_actions()

    def _add_session_internal(self, new_session, is_initial=False):
        """Internal helper to add a session to the list and model, updating history index."""
        if new_session is None or isinstance(new_session, int): return

        # On initial load or explicit add, truncate any potential redo history
        if not is_initial and self.history_index < len(self.active_sessions) - 1:
            del self.active_sessions[self.history_index + 1:]
            self.session_model.removeRows(self.history_index + 1, len(self.active_sessions) - 1 - self.history_index)

        # Append
        self.active_sessions.append(new_session)
        new_idx = len(self.active_sessions) - 1
        self.session_model.insertRow(new_idx)
        self.session_model.setData(self.session_model.index(new_idx), new_session.name)
        self.history_index = new_idx # Make the newly added session current

    def _update_view_for_session(self, session_to_show: SessionV2, set_current_list_item=False):
        """Updates the central plot widget and UI state for the given session."""
        self.session_manager = session_to_show # Update the reference
        is_valid = bool(session_to_show and session_to_show.spec and len(session_to_show.spec.x) > 0)

        # Update UI visibility etc.
        self._update_ui_state(is_valid)

        if is_valid:
            self.plot_viewer.update_plot(session_to_show) # Update plot content
            if set_current_list_item:
                # Select the corresponding item in the list view
                try:
                    # Index in history IS the index in the list/model now
                    q_model_index = self.session_model.index(self.history_index)
                    selection_model = self.session_list_view.selectionModel()
                    selection_flag = QItemSelectionModel.SelectionFlag.ClearAndSelect
                    selection_model.setCurrentIndex(q_model_index, selection_flag)
                except Exception as e:
                    logging.warning(f"Could not select session index {self.history_index} in list view: {e}")
        else:
             logging.info("Updating to an empty session view.")

    def _undo_last_action(self):
        """Switches the view to the previous state in the history."""
        if self.history_index > 0:
            self.history_index -= 1
            logging.debug(f"Undo: Switching to history index {self.history_index}")
            session_to_show = self.active_sessions[self.history_index]
            self._update_view_for_session(session_to_show, set_current_list_item=True)
            self._update_undo_redo_actions()
        else:
            logging.debug("Undo: Already at oldest state.")

    def _redo_last_action(self):
        """Switches the view to the next state in the history."""
        if self.history_index < len(self.active_sessions) - 1:
            self.history_index += 1
            logging.debug(f"Redo: Switching to history index {self.history_index}")
            session_to_show = self.active_sessions[self.history_index]
            self._update_view_for_session(session_to_show, set_current_list_item=True)
            self._update_undo_redo_actions()
        else:
            logging.debug("Redo: Already at newest state.")

    def _update_undo_redo_actions(self):
        """Enables/disables Undo/Redo actions based on history index."""
        can_undo = self.history_index > 0
        can_redo = self.history_index < len(self.active_sessions) - 1

        # ** Add Debug Print **
        print(f"DEBUG _update_undo_redo: history_index={self.history_index}, "
              f"len(active_sessions)={len(self.active_sessions)}, "
              f"Can Undo={can_undo}, Can Redo={can_redo}")
        # ********************

        if hasattr(self, 'undo_action'):
            self.undo_action.setEnabled(can_undo)
        if hasattr(self, 'redo_action'):
            self.redo_action.setEnabled(can_redo)

    def _on_open_spectrum(self):
        """Launches the file dialog and initiates V2 loading."""
        file_name, _ = QFileDialog.getOpenFileName(
            self, 
            "Open Spectrum File", 
            os.getcwd(), # Start in current working directory
            "Astrocook Sessions (*.acs *.acs2);;FITS Files (*.fits);;All Files (*)" # Updated filter
        )
        
        if file_name:
            # We assume a fixed format name for now, as V1 auto-detection is complex
            format_name = 'generic_spectrum'
            # Extract just the filename stem for the session name
            session_name = os.path.splitext(os.path.basename(file_name))[0]


            # Make sure a valid GUI context placeholder exists
            # If session_manager might be None initially, handle it
            gui_context = getattr(self.session_manager, '_gui', None)
            if gui_context is None:
                 # Create a temporary mock if no session exists yet
                 from tests.launch_pyside_app import MockGUI # Assuming mock is accessible
                 gui_context = MockGUI()
                 logging.warning("No active session, using temporary GUI context for loading.")


            try:
                # CRITICAL FIX: Call the utility function with 'archive_path'
                new_session = load_session_from_file(
                    archive_path=file_name, # <<< CORRECT ARGUMENT NAME
                    name=session_name, # Use extracted name
                    format_name=format_name,
                    gui_context=gui_context
                )

                if new_session == 0: # Check for V1-style failure
                     raise RuntimeError("load_session_from_file returned failure code 0.")

                # Add the new session to the application state
                self.add_session(new_session)

            except Exception as e:
                logging.error(f"Failed to load file via V2 adapter: {e}")
                # Consider showing an error message box to the user here
                # from PySide6.QtWidgets import QMessageBox
                # QMessageBox.critical(self, "Error Loading File", f"Could not load {file_name}:\n{e}")

    def _on_save_session(self):
        """Handles the File > Save Session action."""
        if not self.session_manager or not self.session_manager.spec:
            QMessageBox.warning(self, "No Session", "No active session to save.")
            return

        # Suggest filename based on session name, default to .acs2
        default_name = self.session_manager.name + ".acs2"
        # Use QFileDialog to get save path
        file_name, selected_filter = QFileDialog.getSaveFileName(
            self,
            "Save Astrocook Session",
            default_name, # Default filename
            "Astrocook V2 Session (*.acs2);;Astrocook V1 Session (*.acs)" # Filters
        )

        if file_name:
            logging.info(f"Saving session to: {file_name}")
            try:
                # Call the SessionV2 save method
                # The 'models' argument might be irrelevant for V2, check save() signature
                result = self.session_manager.save(file_path=file_name) # Pass file_path
                if result != 0: # Check for V1-style error code
                     raise RuntimeError("Session save method returned non-zero error code.")
                logging.info("Session saved successfully.")
                # Optionally update window title or status bar
            except Exception as e:
                logging.error(f"Failed to save session to {file_name}: {e}", exc_info=True)
                QMessageBox.critical(self, "Save Error", f"Could not save session:\n{e}")

    def _launch_x_convert_dialog(self):
        # Placeholder function called by the QAction
        print("Launching X Convert Dialog...")
        
    def _launch_y_convert_dialog(self):
        # Placeholder function called by the QAction
        print("Launching Y Convert Dialog...")
        
    def _launch_rebin_dialog(self):
        # Placeholder function called by the QAction
        print("Launching Rebin Dialog...")
        

# --- Modify add_session (used by _on_open_spectrum) ---
    def add_session(self, new_session, initial_load=False):
        """Adds a new session, potentially clearing redo history, and updates view."""
        self._add_session_internal(new_session, is_initial=initial_load)
        # Update view to show the newly added session
        self._update_view_for_session(new_session, set_current_list_item=True)
        self._update_undo_redo_actions() # Update enabled state

    # --- Modify _on_session_switched ---
    def _on_session_switched(self, index):
        """Sets the active history index when the list view is clicked."""
        try:
            new_history_index = index.row()
            if 0 <= new_history_index < len(self.active_sessions):
                 if new_history_index != self.history_index: # Only update if changed
                     logging.debug(f"Session list clicked: Setting history index to {new_history_index}")
                     self.history_index = new_history_index
                     session_to_show = self.active_sessions[self.history_index]
                     self._update_view_for_session(session_to_show, set_current_list_item=False) # View update, no need to re-select list item
                     self._update_undo_redo_actions()
            else:
                 logging.warning(f"Invalid index {new_history_index} clicked in session list.")
        except Exception as e:
            logging.error(f"Error switching session via list click: {e}")