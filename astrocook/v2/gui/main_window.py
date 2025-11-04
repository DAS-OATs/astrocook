import json
import logging
import os
from PySide6.QtCore import (
    Qt,
    QItemSelectionModel,
    QLocale,
    QPropertyAnimation, QEasingCurve, QRect, QPoint,
    QParallelAnimationGroup,
    QSize, QStringListModel, QTimer
)
from PySide6.QtGui import QAction, QDoubleValidator, QKeySequence 
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialog, QFileDialog, 
    QMainWindow, QWidget, QVBoxLayout, QFormLayout, QLabel, QLineEdit, QListView, 
    QMenu, QMessageBox,
    QPushButton, QSizePolicy, QSpacerItem, QStackedWidget, QStyle
)
import re
from typing import List, Optional

from .log_scripter_dialog import LogScripterDialog
from .pyside_plot import SpectrumPlotWidget
from .recipe_dialog import RecipeDialog
from ..session_manager import SessionHistory
from ..session import SessionV2, load_session_from_file, LogManager
from ..structures import HistoryLogV2, V1LogArtifact
from ..utils import guarded_deepcopy_v1_state, get_recipe_schema, is_branching_recipe # Import recipe helpers
from ...v1.gui_log import GUILog
from ...v1.defaults import Defaults
try:
    from ...v1.functions import trans_parse
    from ...v1.vars import xem_d
    V1_FUNCTIONS_AVAILABLE = True
except ImportError:
    logging.error("Could not import V1 functions (trans_parse, xem_d) needed for redshift cursor.")
    V1_FUNCTIONS_AVAILABLE = False

# --- *** RECIPE LOOKUP MAP *** ---
from ..recipes.continuum import CONTINUUM_RECIPES_SCHEMAS
from ..recipes.edit import EDIT_RECIPES_SCHEMAS
from ..recipes.flux import FLUX_RECIPES_SCHEMAS
RECIPE_CATEGORY_MAP = {}
for name in CONTINUUM_RECIPES_SCHEMAS: 
    RECIPE_CATEGORY_MAP[name] = 'continuum'
for name in EDIT_RECIPES_SCHEMAS:
    RECIPE_CATEGORY_MAP[name] = 'edit'
for name in FLUX_RECIPES_SCHEMAS:
    RECIPE_CATEGORY_MAP[name] = 'flux'
# (Add other recipe categories here as they are created)

# --- Constants for Sidebar Widths ---
LEFT_SIDEBAR_WIDTH = 250
RIGHT_SIDEBAR_WIDTH = 250
ANIMATION_DURATION = 150 # ** Speed up animation **
BUTTON_WIDTH = 20
BUTTON_HEIGHT = 30

class MockV1GUIContext:
    """Provides the minimal methods expected by V1 loaders."""
    def __init__(self):
        # V1 loaders might expect _flags attribute, even if empty
        self._flags = []
    def _flags_cond(self, flag):
        # V2 GUI doesn't use these flags, so always return False
        return False
    def _flags_extr(self, flag):
        # V2 GUI doesn't use these flags, so return None
        return None
class MainWindowV2(QMainWindow):
    def __init__(self, initial_session: SessionV2, initial_log_object: Optional[LogManager]):
        super().__init__()
        self.session_histories: List[SessionHistory] = [] # List of history managers
        self.active_history: Optional[SessionHistory] = None # Reference to the selected manager
        self.session_model = QStringListModel()
        
        self.log_scripter_dialog: Optional[LogScripterDialog] = None
        self.active_recipe_dialog: Optional[QDialog] = None

        self.recipe_category_map = RECIPE_CATEGORY_MAP

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
        
        # --- ** Sidebars are Children of the MainWindow, floating above ** ---
        self._setup_left_sidebar()  # Creates self.left_sidebar_widget
        self._setup_right_sidebar() # Creates self.right_sidebar_widget
        self._setup_collapse_buttons()

        # --- ** Central Stack Views ** ---
        self._setup_plot_view(initial_session)
        self._setup_empty_view()
        
        # --- 3. Set Central Widget THIRD ---
        self.setCentralWidget(self.central_stack) # Stack fills the window initially

        # --- 4. CRITICAL FIX: Raise floating widgets LAST ---
        # (This brings them visually on top of the central widget)
        self.left_sidebar_widget.raise_()
        self.right_sidebar_widget.raise_()
        self.session_collapse_button.raise_()
        self.plot_controls_collapse_button.raise_()

        # --- Menubar, Styles, Initial State ---
        self._create_menubar()
        self._apply_styles()

        # --- Initial State ---
        is_initial_session_valid = bool(initial_session and initial_session.spec and len(initial_session.spec.x) > 0)
        
        # ** Create the Mock GUI context once for all loggers **
        self.mock_gui_context = MockV1GUIContext()

        if is_initial_session_valid:
            log_object = HistoryLogV2() # <<< ALWAYS create a new, empty log
            
            initial_history = SessionHistory(initial_session, log_object)
            
            # Link the *real* GUI (self) to the session for recipes
            initial_session._gui = self
            
            initial_session.log = GUILog(self.mock_gui_context)
            initial_session.defs = Defaults(self.mock_gui_context)
            # We no longer check for V1LogArtifact here, it's always new

            self.session_histories.append(initial_history)
            self.active_history = initial_history
            self.session_model.setStringList([h.display_name for h in self.session_histories])
            self._update_view_for_session(initial_history.current_state, set_current_list_item=True, is_startup=True)
        else:
            self.active_history = None
            self._update_view_for_session(None, set_current_list_item=False, is_startup=True)

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
        

        self.central_stack.lower()

        # 1. Raise the sidebars (middle layer)
        self.left_sidebar_widget.raise_()
        self.right_sidebar_widget.raise_()
        
        # 2. Raise the buttons (top layer)
        self.session_collapse_button.raise_()
        self.plot_controls_collapse_button.raise_()

    def _setup_plot_view(self, session_for_plot: SessionV2):
        self.plot_viewer = SpectrumPlotWidget(session_for_plot, self)
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

        # 1. Set the context menu policy
        self.session_list_view.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        
        # 2. Connect the signal
        self.session_list_view.customContextMenuRequested.connect(self._on_session_list_context_menu)

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
    
        # --- ** Axis & View Controls ** ---
        view_layout = QVBoxLayout()
        view_layout.setSpacing(8)
        view_layout.addWidget(QLabel("<b>View Controls:</b>"))

        # X-Axis Unit
        form_layout_x = QFormLayout()
        self.x_unit_combo = QComboBox()
        self.x_unit_options = ["nm", "Angstrom", "micron"] # Simplified units
        self.x_unit_combo.addItems(self.x_unit_options)
        self.x_unit_combo.setToolTip("Change X-axis display units (data remains nm).")
        form_layout_x.addRow("X-Unit:", self.x_unit_combo)
        view_layout.addLayout(form_layout_x)
        
        # X/Y Toggles
        self.norm_y_checkbox = QCheckBox("Normalize flux")
        self.norm_y_checkbox.setToolTip("Plot Y / Continuum and set Y-limits.")
        view_layout.addWidget(self.norm_y_checkbox)

        self.log_x_checkbox = QCheckBox("Logarithmic X-Axis")
        view_layout.addWidget(self.log_x_checkbox)
        
        self.log_y_checkbox = QCheckBox("Logarithmic Y-Axis")
        #view_layout.addWidget(self.log_y_checkbox)
        
        sidebar_layout.addLayout(view_layout)
        # -----------------------------

        # Connect new signals
        self.x_unit_combo.currentTextChanged.connect(self._trigger_replot)
        self.norm_y_checkbox.stateChanged.connect(self._trigger_replot)
        self.log_x_checkbox.stateChanged.connect(self._trigger_replot)
        self.log_y_checkbox.stateChanged.connect(self._trigger_replot)

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
        
        self.x_unit_combo.setObjectName("XUnitCombo")
        self.norm_y_checkbox.setObjectName("PlotControlCheckbox")
        self.log_x_checkbox.setObjectName("PlotControlCheckbox")
        self.log_y_checkbox.setObjectName("PlotControlCheckbox")
        
        self.cursor_show_checkbox.setObjectName("PlotControlCheckbox")
        
        self.right_sidebar_widget.resize(0, self.height())
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
            # ** 1. Get the currently active session state **
            session_state = None
            if self.active_history:
                session_state = self.active_history.current_state
            
            # ** 2. Pass the state to plot_spectrum **
            # plot_spectrum itself handles the None case
            self.plot_viewer.plot_spectrum(session_state=session_state)

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

        # --- File Menu Actions
        open_action = QAction("&Open Spectrum...", self); open_action.setShortcut("Ctrl+O"); open_action.triggered.connect(self._on_open_spectrum)
        file_menu.addAction(open_action)
        save_action = QAction("&Save Session...", self); save_action.setShortcut("Ctrl+S"); save_action.triggered.connect(self._on_save_session); 
        save_action.setEnabled(False)
        file_menu.addAction(save_action)
        self.save_action = save_action # Assign to self
        file_menu.addSeparator()
        close_action = QAction("&Close Session", self)
        close_action.setShortcut("Ctrl+W")
        close_action.triggered.connect(self._on_close_session_requested)
        close_action.setEnabled(False)
        file_menu.addAction(close_action)
        self.close_session_action = close_action

        # 'VIEW' MENU ACTIONS
        toggle_left_action = QAction("Toggle Session Panel", self); toggle_left_action.triggered.connect(self._toggle_left_sidebar)
        view_menu.addAction(toggle_left_action)
        toggle_right_action = QAction("Toggle Plot Controls", self); toggle_right_action.triggered.connect(self._toggle_right_sidebar)
        view_menu.addAction(toggle_right_action)
        self.toggle_left_action = toggle_left_action
        self.toggle_right_action = toggle_right_action

        view_menu.addSeparator()
        self.view_log_action = QAction("View Session &Log", self)
        self.view_log_action.triggered.connect(self._on_view_log)
        self.view_log_action.setEnabled(False) # Disabled by default
        view_menu.addAction(self.view_log_action)

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

        # RECIPES FOR 'EDIT' MENU

        # set_properties Action
        set_props_action = QAction("Set &Properties...", self)
        set_props_action.setToolTip("Set core session properties (e.g., z_em)")
        set_props_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "set_properties"))
        edit_menu.addAction(set_props_action)
        self.set_properties_action = set_props_action; self.set_properties_action.setEnabled(False)
        edit_menu.addSeparator()
        
        # apply_expression Action
        apply_expr_action = QAction("Apply &Expression...", self)
        apply_expr_action.setToolTip("Apply a numerical expression to columns (e.g., 'y / 2.0')")
        apply_expr_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "apply_expression"))
        edit_menu.addAction(apply_expr_action)
        self.apply_expression_action = apply_expr_action; self.apply_expression_action.setEnabled(False)
        
        # mask_expression Action
        mask_expr_action = QAction("&Mask by Expression...", self)
        mask_expr_action.setToolTip("Mask a column using a boolean expression (e.g., 'x < 300')")
        mask_expr_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "mask_expression"))
        edit_menu.addAction(mask_expr_action)
        self.mask_expression_action = mask_expr_action; self.mask_expression_action.setEnabled(False)

        # split Action
        split_action = QAction("S&plit Spectrum...", self)
        split_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "split"))
        edit_menu.addAction(split_action)
        self.split_action = split_action; split_action.setEnabled(False)

        # RECIPES FOR 'FLUX' MENU
        
        # rebin Action
        rebin_action = QAction("&Rebin Spectrum...", self)
        rebin_action.triggered.connect(lambda: self._launch_recipe_dialog("flux", "rebin"))
        flux_menu.addAction(rebin_action)
        self.rebin_action = rebin_action; rebin_action.setEnabled(False)

        # resample Action
        resample_action = QAction("Re&sample on Grid...", self)
        resample_action.setToolTip("Resample this spectrum onto another session's grid")
        resample_action.triggered.connect(lambda: self._launch_recipe_dialog("flux", "resample"))
        flux_menu.addAction(resample_action)
        self.resample_action = resample_action; self.resample_action.setEnabled(False)

        # RECIPES FOR 'CONTINUUM' MENU

        auto_cont_action = QAction("&Auto-estimate Continuum...", self)
        auto_cont_action.setToolTip("Find unabsorbed regions and fit a continuum")
        auto_cont_action.triggered.connect(lambda: self._launch_recipe_dialog("continuum", "estimate_auto"))
        continuum_menu.addAction(auto_cont_action)
        self.auto_cont_action = auto_cont_action; self.auto_cont_action.setEnabled(False)

        continuum_menu.addSeparator()

        find_unabs_action = QAction("&Find Unabsorbed Regions...", self)
        find_unabs_action.setToolTip("Create a 'mask_unabs' column using V1 'clip_flux' logic")
        find_unabs_action.triggered.connect(lambda: self._launch_recipe_dialog("continuum", "find_unabsorbed"))
        continuum_menu.addAction(find_unabs_action)
        self.find_unabs_action = find_unabs_action; self.find_unabs_action.setEnabled(False)

        fit_cont_action = QAction("Fit &Continuum to Mask...", self)
        fit_cont_action.setToolTip("Fit a continuum to the 'mask_unabs' column (V1 logic)")
        fit_cont_action.triggered.connect(lambda: self._launch_recipe_dialog("continuum", "fit_continuum"))
        continuum_menu.addAction(fit_cont_action)
        self.fit_cont_action = fit_cont_action; self.fit_cont_action.setEnabled(False)

        self._update_undo_redo_actions()

    def _launch_recipe_dialog(self, category, name):
        """Creates and executes the RecipeDialog for the specified recipe."""
        if self.active_recipe_dialog:
            self.active_recipe_dialog.activateWindow()
            return
        
        # ** Get current state from the *active history manager* **
        if not self.active_history or not self.active_history.current_state.spec:
            QMessageBox.warning(self, "No Session", "Please load a spectrum before running a recipe.")
            return

        logging.info(f"Launching dialog for recipe: {category}.{name}")
        
        # ** Pass the active history's *current state* to the dialog **
        current_state = self.active_history.current_state

        # 1. Use parent=None. This fixes the sidebar event bug.
        # We assume RecipeDialog has its __init__ updated to no longer
        # require the 'main_window_ref' argument.
        dialog = RecipeDialog(category, name, current_state, self) 

        # 2. Run the modal dialog.
        result = dialog.exec()

        QTimer.singleShot(0, self._force_restack_floating_widgets)

        if result == QDialog.Accepted:
            logging.debug(f"Recipe dialog {name} accepted.")
            # Note: The dialog's accept() method now handles
            # calling self.update_gui_session_state()
        else:
            logging.debug(f"Recipe dialog {name} cancelled.")

    def _force_restack_floating_widgets(self):
        """
        Slot to be called by QTimer.
        This is the programmatic equivalent of the user manually
        resizing the window, which we know fixes the bug.
        """
        logging.debug("QTimer: Forcing native resize event to fix widget stacking.")
        try:
            # 1. Get the current size
            current_size = self.size()
            
            # 2. "Jiggle" the size: resize by +1 pixel
            # This is invisible but forces the native resizeEvent to run
            self.resize(current_size.width() + 1, current_size.height())
            
            # 3. Immediately resize back to the original size
            self.resize(current_size)
            
            # The resizeEvent itself will call _reposition_floating_widgets,
            # which will handle the final Z-stacking.
            
        except Exception as e:
            logging.warning(f"Failed to force resize: {e}")

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
            sidebar_a = 0.9 # Alpha value (0.0 to 1.0), e.g., 90% opaque

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

    def update_gui_session_state(self, new_session: SessionV2, original_session_index: int, # original_index no longer needed
                                 is_branching: bool):
        """
        Updates the GUI state AND history after a recipe returns a new session.
        Operates on the currently active SessionHistory.
        """
        if not isinstance(new_session, SessionV2):
            logging.error("update_gui_session_state received invalid session.")
            return
            
        # ** Use the active_history as the target for the update **
        target_history = self.active_history
        
        if target_history is None:
            logging.error("Cannot update state: No active session history.")
            return # Should not happen if a session is loaded

        # Ensure the new state also points to the real GUI
        new_session._gui = self

        if is_branching:
            logging.debug(f"Branching: Creating new SessionHistory from current state.")
            # 1. Get the source log and create a V1-safe deep copy
            source_log_manager = target_history.log_manager
            new_log_copy = guarded_deepcopy_v1_state(source_log_manager)

            if new_log_copy is None:
                logging.error("Failed to deepcopy GUILog for branching, aborting.")
                QMessageBox.critical(self, "Branching Error", "Could not copy session log state.")
                return

            # 2. Update the new state from the recipe to point to the *new* log
            new_session.log_manager = new_log_copy
            
            # ** FIX: The new state's .log attribute also needs the *new* log object **
            if isinstance(new_log_copy, GUILog):
                new_session.log = new_log_copy
            elif isinstance(new_log_copy, V1LogArtifact):
                 # Create a *new* GUILog stub populated from the V1 artifact
                 new_session.log = GUILog(self.mock_gui_context)
                 try: new_session.log.str = json.dumps(new_log_copy.v1_json)
                 except Exception: pass
            else: # It's a HistoryLogV2
                 # We still need a V1 .log stub for V1 recipes
                 new_session.log = GUILog(self.mock_gui_context)


            # 3. Create the new history manager
            # ** FIX: The new branch *name* should be distinct **
            # We can get the base name from the original state
            base_name = target_history.states[0].name
            # Find a unique name (e.g., "Session-1", "Session-2")
            i = 1
            new_branch_name = f"{base_name}-{i}"
            existing_names = [h.display_name for h in self.session_histories]
            while new_branch_name in existing_names:
                i += 1
                new_branch_name = f"{base_name}-{i}"
            
            # ** FIX: The *new* session state needs the *new* name **
            new_session.name = new_branch_name
            
            # ** Get a copy of the *initial* state to start the branch **
            # We must copy it so we can change its name and log links
            initial_state_for_branch = target_history.states[0]
            
            # We need a *full copy* of the session state, not just a reference
            # Let's try using the .with_new_spectrum() trick
            initial_state_copy = initial_state_for_branch.with_new_spectrum(
                initial_state_for_branch.spec
            )
            # And copy the systs list too
            initial_state_copy = initial_state_copy.with_new_system_list(
                initial_state_for_branch.systs
            )

            initial_state_copy._gui = self # Ensure copy also has real GUI

            initial_state_copy.name = new_branch_name
            initial_state_copy.log_manager = new_log_copy
            initial_state_copy.log = new_session.log # Use the new log stub

            # 4. Create the new history manager with the *copied* initial state
            new_history = SessionHistory(initial_state_copy, new_log_copy)

            # 5. Copy the *rest* of the states, re-linking their log reference
            # We must copy *all* intermediate states
            other_states_to_copy = []
            for state in target_history.states[1:target_history.current_index + 1]:
                # Create a copy of each state
                state_copy = state.with_new_spectrum(state.spec)
                state_copy = state_copy.with_new_system_list(state.systs)
                
                state_copy._gui = self # Ensure copy also has real GUI

                state_copy.name = new_branch_name # Ensure name consistency
                state_copy.log_manager = new_log_copy # Re-link to the new log
                state_copy.log = new_session.log # Re-link V1 stub
                other_states_to_copy.append(state_copy)
            
            # 6. Manually set the full state list for the new history
            new_history.states.extend(other_states_to_copy) # Add the re-linked states
            new_history.current_index = target_history.current_index # Set index to match

            # 7. Add the *new* state from the recipe to this *new* history
            new_history.add_state(new_session) # This correctly handles truncation and index update

            # Add the new history manager to the main list
            self.session_histories.append(new_history)
            self.active_history = new_history # New branch becomes active
            # Update the sidebar model
            self.session_model.setStringList([h.display_name for h in self.session_histories])
            new_list_index = len(self.session_histories) - 1
            # Update view for the state in the new history
            self._update_view_for_session(new_history.current_state, set_current_list_item=True, target_list_index=new_list_index)

        else: # Linear update
            logging.debug(f"Linear update: Adding state to active SessionHistory.")
            # Add state to the *current* history manager
            target_history.add_state(new_session) # Handles truncation and index update
            # Update the name in the sidebar model *if* the active history corresponds to the last item
            # (This keeps the name updated for the current linear path)

            # Refresh the *single* log viewer if it's open
            if self.log_scripter_dialog and self.log_scripter_dialog.isVisible():
                try:
                    self.log_scripter_dialog.refresh()
                except Exception as e:
                    logging.error(f"Failed to refresh log viewer: {e}")

            try:
                active_list_index = self.session_histories.index(target_history)
                self.session_model.setData(self.session_model.index(active_list_index), target_history.display_name) # Update name if needed
            except (ValueError, IndexError): pass
            # Update view for the new state in the current history
            self._update_view_for_session(target_history.current_state, set_current_list_item=True) # List item selection doesn't change

        self._update_undo_redo_actions()
    
    def _add_session_internal(self, new_session, is_initial=False):
        """Creates and adds a new SessionHistory manager."""
        if new_session is None or isinstance(new_session, int): return

        # Per Point 3, we *always* create a new, empty log
        log_object = HistoryLogV2()
        
        # 1. Always create a new history for a newly loaded session
        new_history = SessionHistory(new_session, log_object)
        
        # Overwrite the mock GUI context with the real one
        new_session._gui = self
        
        # 2. Manually set the V1 stubs (log and defs) on the new session
        new_session.log = GUILog(self.mock_gui_context)
        new_session.defs = Defaults(self.mock_gui_context)
        
        # We no longer need to check for V1LogArtifact
        
        self.session_histories.append(new_history)
        self.active_history = new_history # Newly loaded becomes active
        # Update model
        self.session_model.setStringList([h.display_name for h in self.session_histories])

    def add_session(self, new_session: SessionV2, initial_load=False):
        """Adds a new session history and updates view."""
        # The log_object parameter has been removed
        self._add_session_internal(new_session, is_initial=initial_load)
        if self.active_history:
            new_list_index = len(self.session_histories) - 1
            self._update_view_for_session(self.active_history.current_state, set_current_list_item=True, target_list_index=new_list_index)
        self._update_undo_redo_actions()

    def _update_view_for_session(self, session_state_to_show: Optional[SessionV2],
                                 set_current_list_item=False, target_list_index=None, is_startup=False):
        """Updates the central plot widget and UI state for the given session state."""
        self.session_manager = session_state_to_show # Keep for plot widget compatibility
        is_valid = bool(session_state_to_show and session_state_to_show.spec and len(session_state_to_show.spec.x) > 0)
        
        # Add refresh call for the single log viewer
        if self.log_scripter_dialog and self.log_scripter_dialog.isVisible():
            if session_state_to_show and self.active_history: # Check for active_history
                self.log_scripter_dialog.set_log_object(self.active_history.log_manager)
            else:
                self.log_scripter_dialog.set_log_object(None) # Clear viewer

        # Update general UI visibility etc. based on validity
        #self._update_ui_state(is_valid, is_startup=is_startup)

        # ** Update X Unit Combo Box **
        # This part is different. We don't read from the session (which is always nm).
        # We *could* store the user's *preference* somewhere, but for now, let's just
        # reset it to "nm" when a new session is shown, or leave it as is.
        # Let's leave it as is, the plot function will just obey it.
        # We DO need to reset the Normalize/Log toggles to their defaults.
        
        # ** FIX: Only reset toggles if the *history* object changed **
        # Check if the session_state_to_show belongs to the same history as the plot
        # (This is tricky. Let's just reset if it's a *new* session state)
        # We need a better way to track "view preferences" per session.
        # For now, let's reset them *unless* it's an undo/redo.
        
        # Let's simplify: _update_view_for_session is called for *any*
        # view change. We *should* reset the toggles unless we store
        # the view state somewhere.
        
        if is_valid:
            # Check if the plot viewer's *current* session matches the *new* one
            # If they match, it's just a replot, don't reset toggles.
            # (This check is imperfect, plot_viewer doesn't store session)
            
            # Let's stick to the previous logic: always reset view toggles
            # when showing a state, as we don't persist view preferences.
            self.norm_y_checkbox.blockSignals(True)
            self.log_x_checkbox.blockSignals(True)
            self.log_y_checkbox.blockSignals(True)
            self.norm_y_checkbox.setChecked(False)
            self.log_x_checkbox.setChecked(False)
            self.log_y_checkbox.setChecked(False)
            self.norm_y_checkbox.blockSignals(False)
            self.log_x_checkbox.blockSignals(False)
            self.log_y_checkbox.blockSignals(False)
            # Let x_unit_combo keep its value.

        # Update Plot Widget
        if is_valid:
            self.plot_viewer.update_plot(session_state_to_show)
        else:
            self.plot_viewer.update_plot(None) # Clear plot if session is None

        #    Call this *after* the plot draw has completed,
        #    so that _reposition_floating_widgets is the *last*
        #    operation, resetting any corrupted state.
        self._update_ui_state(is_valid, is_startup=is_startup)

        # Update List View Selection
        if set_current_list_item:
            list_index_to_select = target_list_index # Use specific index if provided (for branching/add)
            if list_index_to_select is None: # Otherwise, find index of active history
                 try: list_index_to_select = self.session_histories.index(self.active_history)
                 except (ValueError, AttributeError): list_index_to_select = -1

            if list_index_to_select >= 0:
                try:
                    q_model_index = self.session_model.index(list_index_to_select)
                    selection_model = self.session_list_view.selectionModel()
                    selection_flag = QItemSelectionModel.SelectionFlag.ClearAndSelect
                    selection_model.setCurrentIndex(q_model_index, selection_flag)
                except Exception as e:
                    logging.warning(f"Could not select session index {list_index_to_select} in list view: {e}")

        self.left_sidebar_widget.raise_()
        self.right_sidebar_widget.raise_()
        self.session_collapse_button.raise_()
        self.plot_controls_collapse_button.raise_()

    def _undo_last_action(self):
        """Switches the view to the previous state in the active history."""
        if self.active_history:
            previous_state = self.active_history.undo()
            if previous_state:
                logging.debug(f"Undo: Switched to state index {self.active_history.current_index}")
                self._update_view_for_session(previous_state, set_current_list_item=False) 
                self._update_undo_redo_actions()
                
                # --- MODIFIED CHECK ---
                if self.log_scripter_dialog and self.log_scripter_dialog.isVisible():
                    self.log_scripter_dialog.refresh()
                # --- END MODIFIED CHECK ---
            else:
                logging.debug("Undo: Already at oldest state for this session.")
        else:
             logging.debug("Undo: No active session.")

    def _redo_last_action(self):
        """Switches the view to the next state in the active history."""
        if self.active_history:
            next_state = self.active_history.redo()
            if next_state:
                logging.debug(f"Redo: Switched to state index {self.active_history.current_index}")
                self._update_view_for_session(next_state, set_current_list_item=False)
                self._update_undo_redo_actions()

                # --- MODIFIED CHECK ---
                if self.log_scripter_dialog and self.log_scripter_dialog.isVisible():
                    self.log_scripter_dialog.refresh()
                # --- END MODIFIED CHECK ---
            else:
                logging.debug("Redo: Already at newest state for this session.")
        else:
             logging.debug("Redo: No active session.")

    def _update_undo_redo_actions(self):
        """Enables/disables Undo/Redo based on the active history manager."""
        can_undo = self.active_history.can_undo() if self.active_history else False
        can_redo = self.active_history.can_redo() if self.active_history else False

        if hasattr(self, 'undo_action'): self.undo_action.setEnabled(can_undo)
        if hasattr(self, 'redo_action'): self.redo_action.setEnabled(can_redo)


    def _on_session_switched(self, index):
        """Sets the active SessionHistory manager when the list view is clicked."""
        try:
            new_history_list_index = index.row()
            if 0 <= new_history_list_index < len(self.session_histories):
                new_active_history = self.session_histories[new_history_list_index]
                if new_active_history != self.active_history:
                    logging.debug(f"Session list clicked: Setting active history to index {new_history_list_index}")
                    self.active_history = new_active_history
                    # Show the state pointed to by the NEW active history's index
                    self._update_view_for_session(self.active_history.current_state, set_current_list_item=False)
                    self._update_undo_redo_actions()

                    # Update the log viewer if it's open
                    if self.log_scripter_dialog and self.log_scripter_dialog.isVisible():
                        self.log_scripter_dialog.set_log_object(self.active_history.log_manager)
                else:
                    # Clicked on the already active session, maybe ensure view reflects current index?
                    # This can happen if user undid then clicked the list item again.
                    self._update_view_for_session(self.active_history.current_state, set_current_list_item=False)
            else:
                 logging.warning(f"Invalid index {new_history_list_index} clicked.")
        except Exception as e:
            logging.error(f"Error switching session via list click: {e}")

    def _on_open_spectrum(self):
        """Launches the file dialog and initiates V2 loading."""
        file_name, _ = QFileDialog.getOpenFileName(
            self, 
            "Open Spectrum File", 
            os.getcwd(),
            "Astrocook Sessions (*.acs *.acs2);;FITS Files (*.fits);;All Files (*)"
        )
        
        if file_name:
            format_name = 'generic_spectrum'
            session_name = os.path.splitext(os.path.basename(file_name))[0]
            gui_context = self.mock_gui_context

            try:
                # load_session_from_file now only returns the session (Point 3)
                new_session = load_session_from_file(
                    archive_path=file_name, 
                    name=session_name, 
                    format_name=format_name,
                    gui_context=gui_context
                )

                if new_session == 0:
                    raise RuntimeError("load_session_from_file returned failure code 0.")

                # add_session no longer takes a log_object (Point 3)
                self.add_session(new_session)

            except Exception as e:
                logging.error(f"Failed to load file via V2 adapter: {e}")
                QMessageBox.critical(self, "Error Loading File", f"Could not load {file_name}:\n{e}")

    def _on_session_list_context_menu(self, pos: QPoint):
        """
        Handles the right-click on the session list.
        """
        index = self.session_list_view.indexAt(pos)
        if not index.isValid():
            return # Clicked on empty space

        row = index.row()
        if 0 <= row < len(self.session_histories):
            # Get the *clicked* history item
            """
        Handles the right-click on the session list.
        """
        index = self.session_list_view.indexAt(pos)
        if not index.isValid():
            return 

        row = index.row()
        if not (0 <= row < len(self.session_histories)):
            return
            
        history_item = self.session_histories[row]
        session_name = history_item.display_name
        menu = QMenu(self)

        # --- View Log action ---
        view_action = QAction(f"View Log for '{session_name}'", self)
        view_action.triggered.connect(
            lambda: self._launch_log_scripter(history_item)
        )
        menu.addAction(view_action)

        # --- Save Session Action ---
        save_action = QAction(f"Save '{session_name}' as...", self)
        save_action.triggered.connect(
            lambda: self._on_save_session_context(history_item)
        )
        menu.addAction(save_action)

        # --- Close Session Action ---
        menu.addSeparator()
        close_action = QAction(f"Close '{session_name}'", self)
        close_action.triggered.connect(
            lambda: self._on_close_session_requested(history_item)
        )
        menu.addAction(close_action)
        
        menu.exec(self.session_list_view.mapToGlobal(pos))

    def _on_view_log(self):
        """
        Handles the "View > View Session Log" menu action.
        Shows/raises the *single, persistent* log scripter.
        """
        if self.active_history:
            self._launch_log_scripter(self.active_history)
        else:
            logging.warning("View Log called with no active history.")

    def _on_save_session_context(self, history_item: SessionHistory):
        """
        Slot for the context menu's "Save" action.
        """
        # We don't need to switch the active session, just save the one clicked
        self._save_session(history_item)

    def _on_save_session(self):
        """Handles the File > Save Session action for the *active* session."""
        self._save_session(self.active_history)

    def _save_session(self, history_to_save: Optional[SessionHistory]):
        """Saves the specified SessionHistory."""
        if not history_to_save or not history_to_save.current_state.spec:
            QMessageBox.warning(self, "No Session", "No active session to save.")
            return False # Indicate failure

        session_to_save_final = history_to_save.current_state
        session_to_save_initial = history_to_save.states[0]
        
        default_name = session_to_save_final.name + ".acs2"
        
        file_name, selected_filter = QFileDialog.getSaveFileName(
            self,
            "Save Astrocook Session",
            default_name,
            "Astrocook V2 Session (*.acs2);;Astrocook V1 Session (*.acs)"
        )

        if file_name:
            logging.info(f"Saving session to: {file_name}")
            try:
                # Pass *both* sessions to the save method
                result = session_to_save_final.save(
                    file_path=file_name,
                    initial_session=session_to_save_initial
                )
                logging.info("Session saved successfully.")
                return True # Indicate success
            except Exception as e:
                logging.error(f"Failed to save session to {file_name}: {e}", exc_info=True)
                QMessageBox.critical(self, "Save Error", f"Could not save session:\n{e}")
                return False
        return False

    # --- *** START NEW "CLOSE" METHODS *** ---
    def _on_close_session_requested(self, history_to_close: Optional[SessionHistory] = None):
        """
        Handles the request to close a session, either from the menu (active)
        or context menu (specific).
        """
        if history_to_close is None:
            history_to_close = self.active_history

        if not history_to_close:
            logging.warning("Close request with no session specified or active.")
            return

        session_name = history_to_close.display_name
        
        # Check for unsaved changes (basic check: has an undo history)
        # A more robust check would be a "is_dirty" flag
        is_dirty = history_to_close.can_undo() or history_to_close.can_redo()
        
        if is_dirty:
            msg_box = QMessageBox(self)
            msg_box.setWindowTitle("Close Session")
            msg_box.setText(f"Session '{session_name}' has unsaved changes.")
            msg_box.setInformativeText("Do you want to save your changes before closing?")
            msg_box.setStandardButtons(QMessageBox.StandardButton.Save |
                                       QMessageBox.StandardButton.Discard |
                                       QMessageBox.StandardButton.Cancel)
            msg_box.setDefaultButton(QMessageBox.StandardButton.Save)
            
            ret = msg_box.exec()

            if ret == QMessageBox.StandardButton.Save:
                saved_ok = self._save_session(history_to_close)
                if saved_ok:
                    self._close_session(history_to_close) # Close only if save succeeded
            elif ret == QMessageBox.StandardButton.Discard:
                self._close_session(history_to_close) # Close without saving
            elif ret == QMessageBox.StandardButton.Cancel:
                return # Do nothing
        else:
            # Not dirty, just close it
            self._close_session(history_to_close)

    def _close_session(self, history_to_close: SessionHistory):
        """
        Performs the actual removal of the session from the GUI.
        """
        try:
            row_to_remove = self.session_histories.index(history_to_close)
        except ValueError:
            logging.error(f"Could not find session {history_to_close.display_name} to close.")
            return
            
        logging.info(f"Closing session: {history_to_close.display_name}")

        # Remove from list
        self.session_histories.pop(row_to_remove)
        
        # Update model
        self.session_model.removeRow(row_to_remove)

        # Clean up the single log viewer IF it was showing the closing session
        if self.log_scripter_dialog and self.log_scripter_dialog.log_object is history_to_close.log_manager:
            self.log_scripter_dialog.close()
        
        # --- Handle switching the active view ---
        new_active_history = None
        if history_to_close is not self.active_history:
            # We closed a background session. Active history is unchanged.
            new_active_history = self.active_history
        elif self.session_histories:
            # We closed the active session, but others remain.
            # Select the one at the same index, or the last one.
            new_index = min(row_to_remove, len(self.session_histories) - 1)
            new_active_history = self.session_histories[new_index]
        else:
            # We closed the last session.
            new_active_history = None

        self.active_history = new_active_history
        
        if self.active_history:
            # Update view to show the new active session
            self._update_view_for_session(self.active_history.current_state, set_current_list_item=True)
        else:
            # Show empty view
            self._update_view_for_session(None, set_current_list_item=False)

        self._update_undo_redo_actions()
                
    def _launch_log_scripter(self, history_object: 'SessionHistory'):
        """
        Creates and shows the *single, persistent* LogScripterDialog,
        or activates it if it already exists.
        """
        log_object_to_view = history_object.log_manager
        session_name = history_object.display_name

        if self.log_scripter_dialog:
            self.log_scripter_dialog.set_log_object(log_object_to_view) 
            self.log_scripter_dialog.setWindowTitle(f"Log Scripter: {session_name}")
            self.log_scripter_dialog.raise_()
            self.log_scripter_dialog.activateWindow()
            logging.debug(f"Activating existing log scripter for {session_name}")
            return

        try:
            dialog = LogScripterDialog(log_object_to_view, session_name, self)
            self.log_scripter_dialog = dialog
            dialog.finished.connect(
                lambda: self._on_log_scripter_closed() 
            )
            dialog.show()
            logging.debug(f"Creating new single log scripter for {session_name}")
        except Exception as e:
            logging.error(f"Failed to launch log scripter: {e}", exc_info=True)

    def _on_log_scripter_closed(self):
        """
        Slot called when the single LogScripterDialog is closed.
        Cleans up the reference.
        """
        self.log_scripter_dialog = None
        logging.debug("Single log scripter closed, reference removed.")
        
        try:
            self.session_collapse_button.raise_()
            self.plot_controls_collapse_button.raise_()
        except Exception as e:
            logging.warning(f"Could not re-raise buttons: {e}")

    # Simple regex to find recipe(args)
    _SCRIPT_LINE_REGEX = re.compile(r'(\w+)\((.*)\)')
    
    def run_script(self, script_text: str):
        """Re-runs an analysis script, replacing the active history."""
        
        if not self.active_history:
            QMessageBox.warning(self, "No Session", "Cannot run script: No active session.")
            return

        logging.info("Starting script run...")

        # 1. Get the original, pristine state
        initial_state = self.active_history.states[0]
        
        # 2. Create a new history to build into
        new_log = HistoryLogV2()
        # We must create a true copy of the initial state
        initial_state_copy = initial_state.with_new_spectrum(
            initial_state.spec
        ).with_new_system_list(
            initial_state.systs
        )
        # Re-link the new state to the GUI
        initial_state_copy._gui = self
        initial_state_copy.log = GUILog(self.mock_gui_context) 
        initial_state_copy.defs = Defaults(self.mock_gui_context)
        
        new_history = SessionHistory(initial_state_copy, new_log)
        current_processing_state = initial_state_copy

        # 3. Process the script line by line
        lines = script_text.splitlines()
        for i, line in enumerate(lines):
            line = line.strip()
            if not line or line.startswith('#'):
                continue # Skip comments and empty lines

            try:
                # 3a. Parse the line
                match = self._SCRIPT_LINE_REGEX.match(line)
                if not match:
                    raise ValueError("Invalid syntax. Expected 'recipe_name(args)'.")
                
                recipe_name = match.group(1)
                args_str = match.group(2)
                
                params_dict = {}
                if args_str:
                    # Use a more robust way to parse args than simple splitting
                    # This finds all 'key=value' pairs
                    arg_pairs = re.findall(r"(\w+)\s*=\s*([^,]+)", args_str)
                    for key, val in arg_pairs:
                        key = key.strip()
                        val = val.strip()
                        # Un-quote strings
                        if (val.startswith("'") and val.endswith("'")) or \
                           (val.startswith('"') and val.endswith('"')):
                            params_dict[key] = val[1:-1]
                        else:
                            params_dict[key] = val # Recipes expect strings

                # 3b. Find and run the recipe
                category = self.recipe_category_map.get(recipe_name)
                if not category:
                    raise ValueError(f"Recipe '{recipe_name}' not found.")
                
                recipe_instance = getattr(current_processing_state, category)
                recipe_method = getattr(recipe_instance, recipe_name)

                # Check for multi-session and add alias map if needed
                if recipe_name in ("apply_expression", "mask_expression"):
                    if not 'alias_map' in params_dict:
                        # Build the alias map relative to the *active* session
                        params_dict['alias_map'] = self.log_scripter_dialog._build_alias_map()
                
                new_session_state = recipe_method(**params_dict)
                
                if not new_session_state or new_session_state == 0:
                    raise ValueError(f"Recipe failed to execute.")
                    
                # 3c. Success! Add to the new history
                new_log.add_entry(recipe_name, params_dict)
                new_history.add_state(new_session_state)
                current_processing_state = new_session_state

            except Exception as e:
                logging.error(f"Failed to run script line {i+1}: {line}\nError: {e}", exc_info=True)
                QMessageBox.critical(
                    self, 
                    "Script Error",
                    f"Failed on line {i+1}:\n{line}\n\nError: {e}"
                )
                return # Stop processing

        # 4. Swap the old history with the new one
        try:
            old_history_index = self.session_histories.index(self.active_history)
            self.session_histories[old_history_index] = new_history
            self.active_history = new_history
            self.session_model.setData(self.session_model.index(old_history_index), new_history.display_name)
        except ValueError:
            logging.error("Could not find old history to replace. Appending new history.")
            self.session_histories.append(new_history)
            self.active_history = new_history
            self.session_model.setStringList([h.display_name for h in self.session_histories])


        # 5. Refresh everything
        logging.info("Script run completed successfully.")
        self._update_view_for_session(self.active_history.current_state, set_current_list_item=True)
        self._update_undo_redo_actions()
        if self.log_scripter_dialog:
            self.log_scripter_dialog.set_log_object(new_log)

    def _launch_x_convert_dialog(self):
        # Placeholder function called by the QAction
        print("Launching X Convert Dialog...")
        
    def _launch_y_convert_dialog(self):
        # Placeholder function called by the QAction
        print("Launching Y Convert Dialog...")
        
    def _launch_rebin_dialog(self):
        # Placeholder function called by the QAction
        print("Launching Rebin Dialog...")

    def _update_ui_state(self, is_valid_session, is_startup=False):
        is_valid_session = bool(is_valid_session)

        # ... (Enable/Disable File menu actions) ...
        if hasattr(self, 'save_action'): self.save_action.setEnabled(is_valid_session)
        if hasattr(self, 'close_session_action'): self.close_session_action.setEnabled(is_valid_session)

        # ... (Enable/Disable View menu actions) ...
        if hasattr(self, 'toggle_left_action'): self.toggle_left_action.setEnabled(is_valid_session)
        if hasattr(self, 'toggle_right_action'): self.toggle_right_action.setEnabled(is_valid_session)
        if hasattr(self, 'view_log_action'): self.view_log_action.setEnabled(is_valid_session)

        # ** Enable/Disable Recipe Actions **
        enable_recipes = is_valid_session
        # Check if actions exist before enabling/disabling
        if hasattr(self, 'set_properties_action'): self.set_properties_action.setEnabled(enable_recipes)
        if hasattr(self, 'apply_expression_action'): self.apply_expression_action.setEnabled(enable_recipes)
        if hasattr(self, 'mask_expression_action'): self.mask_expression_action.setEnabled(enable_recipes)
        if hasattr(self, 'split_action'): self.split_action.setEnabled(enable_recipes)
        
        if hasattr(self, 'rebin_action'): self.rebin_action.setEnabled(enable_recipes)
        if hasattr(self, 'resample_action'): self.resample_action.setEnabled(enable_recipes)
        
        if hasattr(self, 'auto_cont_action'): self.auto_cont_action.setEnabled(enable_recipes)
        if hasattr(self, 'find_unabs_action'): self.find_unabs_action.setEnabled(enable_recipes)
        if hasattr(self, 'fit_cont_action'): self.fit_cont_action.setEnabled(enable_recipes)
        # --- *** END NEW RECIPES *** ---
        
        # ... enable/disable other recipe actions ...

        # Enable Save action only if valid session
        if hasattr(self, 'save_action'): self.save_action.setEnabled(is_valid_session)

        if is_valid_session:
            # --- State when a valid session IS loaded ---
            if self.central_stack.currentIndex() != 0: # If switching from empty
                self.central_stack.setCurrentIndex(0)
            
            # 2. Resize logic (DE-NESTED from the stack index check)
            # Check if this is the very first session being loaded
            was_previously_empty = len(self.session_histories) <= 1
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

        self._update_undo_redo_actions()