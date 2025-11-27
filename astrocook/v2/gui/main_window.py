from copy import deepcopy
import json
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
from PySide6.QtCore import (
    Qt,
    QItemSelectionModel,
    QLocale,
    QPropertyAnimation, QEasingCurve, QRect, QPoint,
    QParallelAnimationGroup,
    QSize, QStringListModel, QThreadPool, QTimer
)
from PySide6.QtGui import QAction, QDoubleValidator, QKeySequence, QIcon
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog, QInputDialog,
    QHBoxLayout, QMainWindow, QWidget, QVBoxLayout, QFormLayout, QLabel, QLineEdit, QListView, 
    QMenu, QMessageBox,
    QPushButton, QProgressDialog, QSizePolicy, QSpacerItem, QStackedWidget, QStyle, QTextEdit,
)
import re
from typing import Any, Dict, List, Optional

from .debug_utils import GLOBAL_PLOTTER
from .identification_viewer_dialog import IdentificationViewerDialog
from .log_scripter_dialog import LogScripterDialog
from ..photometry import STANDARD_FILTERS
from .pyside_plot import SpectrumPlotWidget
from .qt_workers import RecipeWorker, ScriptWorker
from .recipe_dialog import RecipeDialog
from ..session_manager import SessionHistory
from ..session import SessionV2, load_session_from_file, LogManager
from ..structures import HistoryLogV2, V1LogArtifact
from .system_inspector import SystemInspector
from ..utils import guarded_deepcopy_v1_state, get_recipe_schema, is_branching_recipe, resource_path # Import recipe helpers
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
from ..recipes.absorbers import ABSORBERS_RECIPES_SCHEMAS
from ..recipes.continuum import CONTINUUM_RECIPES_SCHEMAS
from ..recipes.edit import EDIT_RECIPES_SCHEMAS
from ..recipes.flux import FLUX_RECIPES_SCHEMAS
RECIPE_CATEGORY_MAP = {}
for name in ABSORBERS_RECIPES_SCHEMAS:
    RECIPE_CATEGORY_MAP[name] = 'absorbers'
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
        self.identification_viewer_dialog: Optional[IdentificationViewerDialog] = None # <<< *** ADD ***
        self.active_recipe_dialog: Optional[QDialog] = None

        self.recipe_category_map = RECIPE_CATEGORY_MAP

        # Initialize threadpool and progress dialog
        self.thread_pool = QThreadPool()
        logging.info(f"Main thread pool started with {self.thread_pool.maxThreadCount()} threads.")
        self.progress_dialog: Optional[QProgressDialog] = None

        # ** Initialize animation attributes to None **
        self.left_animation_group = None
        self.right_animation_group = None

        # ** Add History Index **
        self.history_index = -1 # Index of the currently active session in active_sessions

        self._pending_recipe_on_properties_set: Optional[tuple] = None

        self._last_attempted_recipe: Optional[dict] = None

        #self.setGeometry(100, 100, 450, 150) # Initial small size
        self.resize(1400,900)
        # Centratura robusta (Fix per macOS)
        screen = QApplication.primaryScreen()
        if screen:
            # geometry() include la barra dei menu, availableGeometry() no (meglio!)
            screen_rect = screen.availableGeometry() 
            
            # Prendiamo il rettangolo della nostra finestra
            window_rect = self.frameGeometry()
            
            # Spostiamo il centro del rettangolo finestra al centro dello schermo
            window_rect.moveCenter(screen_rect.center())
            
            # Muoviamo la finestra vera e propria nella nuova posizione calcolata
            self.move(window_rect.topLeft())

        try:
            # Usa pure logo_3d.png per ora. 
            # Se in futuro crei un 'logo_small.png' o 'logo_bw.png', cambia qui il nome.
            icon_path = resource_path(os.path.join("assets", "logo_3d.png"))
            self.setWindowIcon(QIcon(icon_path))
        except Exception as e:
            logging.warning(f"Could not set window icon: {e}")

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

        GLOBAL_PLOTTER.plot_requested.connect(self._on_debug_plot_requested)

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
        label = QLabel("Welcome to Astrocook v2!\n\nUse 'File > Open Spectrum...'")
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

        # 3a. Enable in-place editing
        self.session_list_view.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked | 
            QAbstractItemView.EditTrigger.EditKeyPressed
        )
        
        # 3b. Connect the signal that fires *after* editing is done
        self.session_model.dataChanged.connect(self._on_session_name_changed)

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
        sidebar_layout.setContentsMargins(10, 10, 15, 10)
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

        # --- *** NEW: Auxiliary Column Plotter *** ---
        # (Using a form layout for clean alignment)
        form_layout_aux = QFormLayout()
        form_layout_aux.setSpacing(5)
        self.aux_col_combo = QComboBox()
        self.aux_col_combo.setToolTip("Select an auxiliary column to plot (e.g., a mask)")
        self.aux_col_combo.currentTextChanged.connect(self._trigger_replot)
        form_layout_aux.addRow("Aux Column:", self.aux_col_combo)
        plot_toggles_layout.addLayout(form_layout_aux)
        # --- *** END NEW *** ---

        self.strong_lines_checkbox = QCheckBox("Strong emission lines only")
        self.strong_lines_checkbox.setToolTip("Show only major emission lines when z_em is active.")
        self.strong_lines_checkbox.setChecked(True) # Default to clean view
        self.strong_lines_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.strong_lines_checkbox)

        sidebar_layout.addLayout(plot_toggles_layout) # Add group to main layout

        self.isomag_checkbox = QCheckBox("Show Iso-Mag Grid")
        self.isomag_checkbox.setToolTip("Show lines of constant AB magnitude")
        self.isomag_checkbox.toggled.connect(lambda b: self.plot_viewer.toggle_isomag_grid(b))
        self.isomag_checkbox.setObjectName("PlotControlCheckbox") # Ensure styling
        plot_toggles_layout.addWidget(self.isomag_checkbox)
    
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
        self.norm_y_checkbox.toggled.connect(self._on_view_toggle) # <-- CHANGE THIS
        view_layout.addWidget(self.norm_y_checkbox)

        self.snr_checkbox = QCheckBox("Show SNR (y / error)")
        self.snr_checkbox.setToolTip("Plot Y / Error Column.")
        self.snr_checkbox.toggled.connect(self._on_view_toggle) # <-- NEW
        view_layout.addWidget(self.snr_checkbox)
        form_layout_snr = QFormLayout()
        form_layout_snr.setSpacing(5)
        self.snr_col_combo = QComboBox()
        self.snr_col_combo.setToolTip("Select column to use as error for SNR")
        self.snr_col_combo.currentTextChanged.connect(self._trigger_replot)
        form_layout_snr.addRow("  Error Col:", self.snr_col_combo)
        view_layout.addLayout(form_layout_snr)

        self.log_x_checkbox = QCheckBox("Logarithmic X-Axis")
        view_layout.addWidget(self.log_x_checkbox)
        
        self.log_y_checkbox = QCheckBox("Logarithmic Y-Axis")
        #view_layout.addWidget(self.log_y_checkbox)
        
        sidebar_layout.addLayout(view_layout)
        # -----------------------------

        # Connect new signals
        self.x_unit_combo.currentTextChanged.connect(self._trigger_replot)
        #self.norm_y_checkbox.stateChanged.connect(self._trigger_replot)
        self.log_x_checkbox.stateChanged.connect(self._trigger_replot)
        self.log_y_checkbox.stateChanged.connect(self._trigger_replot)

        # --- *** NEW: Axis Limit Controls *** ---
        limits_layout = QVBoxLayout()
        limits_layout.setSpacing(8)
        limits_layout.addWidget(QLabel("<b>Axis Limits:</b>"))

        # Use QFormLayout for label-input pairs
        form_layout_limits = QFormLayout()
        form_layout_limits.setRowWrapPolicy(QFormLayout.DontWrapRows)
        form_layout_limits.setLabelAlignment(Qt.AlignLeft)
        form_layout_limits.setSpacing(5)

        # Create validators (use C locale for decimal points)
        validator = QDoubleValidator()
        validator.setLocale(QLocale.C)
        validator.setNotation(QDoubleValidator.StandardNotation)

        self.xmin_input = QLineEdit("0.0")
        self.xmin_input.setValidator(validator)
        self.xmax_input = QLineEdit("1.0")
        self.xmax_input.setValidator(validator)
        self.ymin_input = QLineEdit("0.0")
        self.ymin_input.setValidator(validator)
        self.ymax_input = QLineEdit("1.0")
        self.ymax_input.setValidator(validator)

        self.xmin_input.editingFinished.connect(self._on_set_custom_limits)
        self.xmax_input.editingFinished.connect(self._on_set_custom_limits)
        self.ymin_input.editingFinished.connect(self._on_set_custom_limits)
        self.ymax_input.editingFinished.connect(self._on_set_custom_limits)
        
        form_layout_limits.addRow("X Min:", self.xmin_input)
        form_layout_limits.addRow("X Max:", self.xmax_input)
        form_layout_limits.addRow("Y Min:", self.ymin_input)
        form_layout_limits.addRow("Y Max:", self.ymax_input)
        
        limits_layout.addLayout(form_layout_limits)

        #self.set_limits_button = QPushButton("Set Custom Limits")
        #self.set_limits_button.clicked.connect(self._on_set_custom_limits)
        #limits_layout.addWidget(self.set_limits_button)
        
        sidebar_layout.addLayout(limits_layout)

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
        self.aux_col_combo.setObjectName("AuxColumnCombo") # <-- NEW
        self.strong_lines_checkbox.setObjectName("PlotControlCheckbox")
        self.isomag_checkbox.setObjectName("PlotControlCheckbox")

        self.x_unit_combo.setObjectName("XUnitCombo")
        self.norm_y_checkbox.setObjectName("PlotControlCheckbox")
        self.snr_checkbox.setObjectName("PlotControlCheckbox")
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
        self.view_system_inspector_action = QAction("View System &Inspector", self)
        self.view_system_inspector_action.triggered.connect(self._on_view_system_inspector)
        self.view_system_inspector_action.setEnabled(False) # Disabled by default
        view_menu.addAction(self.view_system_inspector_action)

        view_menu.addSeparator()
        self.view_log_action = QAction("View Session &Log", self)
        self.view_log_action.triggered.connect(self._on_view_log)
        self.view_log_action.setEnabled(False) # Disabled by default
        view_menu.addAction(self.view_log_action)

        self.view_identifications_action = QAction("View &Identifications", self)
        self.view_identifications_action.triggered.connect(self._on_view_identifications)
        self.view_identifications_action.setEnabled(False) # Disabled by default
        view_menu.addAction(self.view_identifications_action)

        # 'EDIT' MENU ACTIONS

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

        # smooth_column Action
        smooth_col_action = QAction("Smooth &Column...", self)
        smooth_col_action.setToolTip("Apply Gaussian smoothing to a single column")
        smooth_col_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "smooth_column"))
        edit_menu.addAction(smooth_col_action)
        self.smooth_column_action = smooth_col_action; self.smooth_column_action.setEnabled(False)

        edit_menu.addSeparator()

        # split Action
        split_action = QAction("S&plit Spectrum...", self)
        split_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "split"))
        edit_menu.addAction(split_action)
        self.split_action = split_action; split_action.setEnabled(False)

        # extract_preset Action
        extract_preset_action = QAction("Extract &Preset Region...", self)
        extract_preset_action.setToolTip("Extract standard regions (Forest, etc.) based on z_em")
        extract_preset_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "extract_preset"))
        edit_menu.addAction(extract_preset_action)
        self.extract_preset_action = extract_preset_action; self.extract_preset_action.setEnabled(False)

        # RECIPES FOR 'FLUX' MENU
        
        std_action = QAction("Calculate Running &StdDev...", self)
        std_action.setToolTip("Calculate running standard deviation on a column (e.g., y)")
        std_action.triggered.connect(
            lambda: self._launch_recipe_dialog("flux", "calculate_running_std")
        )
        flux_menu.addAction(std_action)
        self.calculate_running_std_action = std_action; self.calculate_running_std_action.setEnabled(False)

        flux_menu.addSeparator()

        # smooth Action
        smooth_action = QAction("&Smooth Spectrum...", self)
        smooth_action.setToolTip("Apply Gaussian smoothing to the flux")
        smooth_action.triggered.connect(lambda: self._launch_recipe_dialog("flux", "smooth"))
        flux_menu.addAction(smooth_action)
        self.smooth_action = smooth_action; self.smooth_action.setEnabled(False)

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

        flux_menu.addSeparator()

        # flux calibrate Action
        calib_action = QAction("Flux &Calibrate...", self)
        calib_action.setToolTip("Rescale spectrum to match a photometric magnitude")
        calib_action.triggered.connect(lambda: self._launch_recipe_dialog("flux", "calibrate_from_magnitudes"))
        flux_menu.addAction(calib_action)
        self.calibrate_action = calib_action; self.calibrate_action.setEnabled(False)

        # RECIPES FOR 'CONTINUUM' MENU

        auto_cont_action = QAction("&Auto-estimate Continuum...", self)
        auto_cont_action.setToolTip("Find Absorbed regions and fit a continuum")
        auto_cont_action.triggered.connect(lambda: self._launch_recipe_dialog("continuum", "estimate_auto"))
        continuum_menu.addAction(auto_cont_action)
        self.auto_cont_action = auto_cont_action; self.auto_cont_action.setEnabled(False)

        continuum_menu.addSeparator()

        find_abs_action = QAction("&Find Aabsorbed Regions...", self)
        find_abs_action.setToolTip("Create a 'abs_mask' column using V1 'clip_flux' logic")
        find_abs_action.triggered.connect(lambda: self._launch_recipe_dialog("continuum", "find_absorbed"))
        continuum_menu.addAction(find_abs_action)
        self.find_abs_action = find_abs_action; self.find_abs_action.setEnabled(False)

        fit_cont_action = QAction("Fit &Continuum to Mask...", self)
        fit_cont_action.setToolTip("Fit a continuum to the 'abs_mask' column (V1 logic)")
        fit_cont_action.triggered.connect(lambda: self._launch_recipe_dialog("continuum", "fit_continuum"))
        continuum_menu.addAction(fit_cont_action)
        self.fit_cont_action = fit_cont_action; self.fit_cont_action.setEnabled(False)

        continuum_menu.addSeparator()

        fit_pl_action = QAction("Fit &Power-Law...", self)
        fit_pl_action.setToolTip("Fit a power-law to specified rest-frame regions")
        fit_pl_action.triggered.connect(lambda: self._launch_recipe_dialog("continuum", "fit_powerlaw"))
        continuum_menu.addAction(fit_pl_action)
        self.fit_powerlaw_action = fit_pl_action; self.fit_powerlaw_action.setEnabled(False)

        # RECIPES FOR 'ABSORBERS' MENU

        identify_action = QAction("&Identify Absorption Lines...", self)
        identify_action.setToolTip("Automatically identify absorption regions using correlation signals")
        identify_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "identify_lines"))
        absorbers_menu.addAction(identify_action)
        self.identify_lines_action = identify_action; self.identify_lines_action.setEnabled(False)
        
        self._update_undo_redo_actions()

    def _on_set_custom_limits(self):
        """ Reads values from limit boxes and applies them to the plot. """
        if not self.plot_viewer:
            return
            
        try:
            # Read values as-typed (they are in display units)
            display_xmin = float(self.xmin_input.text())
            display_xmax = float(self.xmax_input.text())
            ymin = float(self.ymin_input.text())
            ymax = float(self.ymax_input.text())
        except ValueError:
            logging.warning("Invalid custom limit values entered.")
            return # Ignore if values are not valid floats
            
        selected_unit = self.x_unit_combo.currentText()
        
        if selected_unit == 'Angstrom':
            xmin_nm = display_xmin / 10.0
            xmax_nm = display_xmax / 10.0
        elif selected_unit == 'micron':
            xmin_nm = display_xmin * 1000.0
            xmax_nm = display_xmax * 1000.0
        else: # 'nm'
            xmin_nm = display_xmin
            xmax_nm = display_xmax

        if xmax_nm <= xmin_nm or ymax <= ymin:
            logging.warning("Invalid custom limit range (max <= min).")
            return

        # Set limits on the Matplotlib axes
        # This will automatically trigger the 'on_lim_changed' callback
        # in pyside_plot.py, which schedules a new plot_spectrum draw.
        self.plot_viewer.canvas.axes.set_xlim(xmin_nm, xmax_nm)
        self.plot_viewer.canvas.axes.set_ylim(ymin, ymax)

    def _update_limit_boxes_from_plot(self):
        """ Updates the limit text boxes with the plot's current limits. """
        if not self.plot_viewer:
            return
            
        try:
            xlim_nm = self.plot_viewer.canvas.axes.get_xlim() # These are always in nm
            ylim_plot = self.plot_viewer.canvas.axes.get_ylim()
            
            # --- NEW: Convert nm limits to selected display unit ---
            selected_unit = self.x_unit_combo.currentText()

            if selected_unit == 'Angstrom':
                display_xlim = (xlim_nm[0] * 10.0, xlim_nm[1] * 10.0)
            elif selected_unit == 'micron':
                display_xlim = (xlim_nm[0] / 1000.0, xlim_nm[1] / 1000.0)
            else: # 'nm'
                display_xlim = xlim_nm
            # --- END NEW ---
            
            # Block signals to prevent a feedback loop
            self.xmin_input.blockSignals(True)
            self.xmax_input.blockSignals(True)
            self.ymin_input.blockSignals(True)
            self.ymax_input.blockSignals(True)
            
            # Set text using the converted display values
            self.xmin_input.setText(f"{display_xlim[0]:.4g}")
            self.xmax_input.setText(f"{display_xlim[1]:.4g}")
            self.ymin_input.setText(f"{ylim_plot[0]:.4g}")
            self.ymax_input.setText(f"{ylim_plot[1]:.4g}")
            
        except Exception as e:
            logging.warning(f"Failed to update limit boxes: {e}")
        finally:
            # Always unblock signals
            self.xmin_input.blockSignals(False)
            self.xmax_input.blockSignals(False)
            self.ymin_input.blockSignals(False)
            self.ymax_input.blockSignals(False)

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
                padding: 16px 10px; border: none;
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
                border-radius: 5px;
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

    def _update_window_title(self, session_name: Optional[str] = None):
        """Updates the main window title."""
        app_name = "Astrocook"
        if session_name:
            self.setWindowTitle(f"{session_name} - {app_name}")
        else:
            self.setWindowTitle(app_name)

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
                                 is_branching: bool, auto_show_aux: Optional[str] = None,
                                 force_autoscale: bool = False):
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
            new_log_copy = None
            # Check if it's a V1-style log that needs special handling
            if hasattr(source_log_manager, '_gui'):
                logging.debug("Branching: Found V1-style log, using guarded_deepcopy_v1_state.")
                new_log_copy = guarded_deepcopy_v1_state(source_log_manager)
            else:
                # It's a V2-style log (like HistoryLogV2), which is safe to deepcopy
                logging.debug("Branching: Found V2-style log, using standard deepcopy.")
                try:
                    new_log_copy = deepcopy(source_log_manager)
                except Exception as e:
                    logging.error(f"Standard deepcopy failed on log manager: {e}")
                    new_log_copy = None # Ensure it's None on failure

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
            self._update_view_for_session(
                new_history.current_state, 
                set_current_list_item=True, 
                target_list_index=new_list_index,
                auto_show_aux=auto_show_aux,
                force_autoscale=True
            )

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
            self._update_view_for_session(
                target_history.current_state, 
                set_current_list_item=True,
                auto_show_aux=auto_show_aux,
                force_autoscale=force_autoscale
            ) # List item selection doesn't change

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
            self._update_view_for_session(self.active_history.current_state, set_current_list_item=True, target_list_index=new_list_index, force_autoscale=True)
        self._update_undo_redo_actions()

    def _update_view_for_session(self, session_state_to_show: Optional[SessionV2],
                                 set_current_list_item=False, target_list_index=None, is_startup=False,
                                 force_autoscale: bool = False, auto_show_aux: Optional[str] = None):
        """Updates the central plot widget and UI state for the given session state."""
        self.session_manager = session_state_to_show # Keep for plot widget compatibility
        is_valid = bool(session_state_to_show and session_state_to_show.spec and len(session_state_to_show.spec.x) > 0)
        
        # Add refresh call for the single log viewer
        if self.log_scripter_dialog and self.log_scripter_dialog.isVisible():
            if session_state_to_show and self.active_history: # Check for active_history
                self.log_scripter_dialog.set_log_object(self.active_history.log_manager)
            else:
                self.log_scripter_dialog.set_log_object(None) # Clear viewer

        if session_state_to_show:
            self._update_window_title(session_state_to_show.name)
        else:
            self._update_window_title(None)

        # --- *** NEW: Update Aux Column Combo *** ---
        try:
            self.aux_col_combo.blockSignals(True)
            current_selection = self.aux_col_combo.currentText()
            self.aux_col_combo.clear()
            self.aux_col_combo.addItem("None")
            
            if is_valid:
                all_cols = list(session_state_to_show.spec.t._data_dict.keys())
                core_cols = {'x', 'xmin', 'xmax', 'y', 'dy'}
                aux_col_names = sorted([c for c in all_cols if c not in core_cols])
                if aux_col_names:
                    self.aux_col_combo.addItems(aux_col_names)
            

            target_index = 0
            if auto_show_aux:
                idx = self.aux_col_combo.findText(auto_show_aux)
                if idx != -1:
                    target_index = idx
            
            if target_index == 0 and current_selection != "None":
                 idx = self.aux_col_combo.findText(current_selection)
                 if idx != -1:
                     target_index = idx
            
            self.aux_col_combo.setCurrentIndex(target_index)
        
        except Exception as e:
            logging.warning(f"Failed to update Aux Column combobox: {e}")
        finally:
            self.aux_col_combo.blockSignals(False)
        # --- *** END NEW *** ---

        # --- *** NEW: Update SNR Column Combo *** ---
        try:
            self.snr_col_combo.blockSignals(True)
            current_snr_col = self.snr_col_combo.currentText()
            self.snr_col_combo.clear()
            self.snr_col_combo.addItem("dy") # 'dy' is always the default
            
            if is_valid:
                # Add other potential error columns (e.g., 'running_std', 'cont_err')
                all_cols = list(session_state_to_show.spec.t._data_dict.keys())
                # Find columns that are not 'dy' but might be errors
                potential_err_cols = sorted([
                    c for c in all_cols 
                    if ('err' in c or 'running_std' in c) and c != 'dy'
                ])
                if potential_err_cols:
                    self.snr_col_combo.addItems(potential_err_cols)
            
            index = self.snr_col_combo.findText(current_snr_col)
            if index != -1:
                self.snr_col_combo.setCurrentIndex(index)
            else:
                self.snr_col_combo.setCurrentIndex(0) # Default to 'dy'
        
        except Exception as e:
            logging.warning(f"Failed to update SNR Column combobox: {e}")
        finally:
            self.snr_col_combo.blockSignals(False)
            
        # Also set the visibility of the combo box
        #self.snr_col_combo.parentWidget().setVisible(self.snr_checkbox.isChecked())

        # Update Plot Widget
        if is_valid:
            self.plot_viewer.update_plot(session_state_to_show, force_autoscale=force_autoscale)
        else:
            self.plot_viewer.update_plot(None, force_autoscale=False) # Clear plot if session is None

        #    Call this *after* the plot draw has completed,
        #    so that _reposition_floating_widgets is the *last*
        #    operation, resetting any corrupted state.
        self._update_ui_state(is_valid, is_startup=is_startup)

        # Update System Inspector if open
        if hasattr(self, 'system_inspector') and self.system_inspector and self.system_inspector.isVisible():
            self.system_inspector.set_session(session_state_to_show)

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

    def _on_view_toggle(self, checked: bool):
        """ Handles toggling the main view (Norm vs SNR vs Raw). """
        sender = self.sender()
        
        # 1. Enforce mutual exclusion
        if checked:
            if sender is self.norm_y_checkbox:
                self.snr_checkbox.blockSignals(True)
                self.snr_checkbox.setChecked(False)
                self.snr_checkbox.blockSignals(False)
            elif sender is self.snr_checkbox:
                self.norm_y_checkbox.blockSignals(True)
                self.norm_y_checkbox.setChecked(False)
                self.norm_y_checkbox.blockSignals(False)
        
        # 2. Show/hide the SNR error column selector
        #self.snr_col_combo.parentWidget().setVisible(self.snr_checkbox.isChecked())
        
        # 3. Trigger the redraw
        self._trigger_replot()

    def _undo_last_action(self):
        """Switches the view to the previous state in the active history."""
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
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
                    QTimer.singleShot(0, self._force_restack_floating_widgets)
                else:
                    logging.debug("Undo: Already at oldest state for this session.")
            else:
                 logging.debug("Undo: No active session.")
        finally:
            QApplication.restoreOverrideCursor()

    def _redo_last_action(self):
        """Switches the view to the next state in the active history."""
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
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
                    QTimer.singleShot(0, self._force_restack_floating_widgets)
                else:
                    logging.debug("Redo: Already at newest state for this session.")
            else:
                 logging.debug("Redo: No active session.")
        finally:
            QApplication.restoreOverrideCursor()

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
                    self._update_view_for_session(self.active_history.current_state, set_current_list_item=False, force_autoscale=False)
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

        # --- *** NEW: Info action *** ---
        info_action = QAction(f"Info", self)
        info_action.triggered.connect(
            # Use lambda to pass the specific history item
            lambda checked=False, item=history_item: self._on_session_info(item)
        )
        menu.addAction(info_action)

        # --- View Log action ---
        view_action = QAction(f"View Log", self)
        view_action.triggered.connect(
            lambda: self._launch_log_scripter(history_item)
        )
        menu.addAction(view_action)

        view_ids_action = QAction(f"View Identifications", self)
        view_ids_action.triggered.connect(
            lambda: self._launch_identification_viewer(history_item)
        )
        # Enable only if identifications exist on this item's *current* state
        has_ids = False
        if history_item and history_item.current_state and history_item.current_state.spec:
            has_ids = history_item.current_state.spec.meta.get('region_identifications') is not None
        view_ids_action.setEnabled(has_ids)
        menu.addAction(view_ids_action)

        menu.addSeparator()
        
        # --- Save Session Action ---
        save_action = QAction(f"Save as...", self)
        save_action.triggered.connect(
            lambda: self._on_save_session_context(history_item)
        )
        menu.addAction(save_action)

        # --- Close Session Action ---
        close_action = QAction(f"Close", self)
        close_action.triggered.connect(
            lambda: self._on_close_session_requested(history_item)
        )
        menu.addAction(close_action)
        
        menu.exec(self.session_list_view.mapToGlobal(pos))
    
    def _on_session_info(self, history_item: SessionHistory):
        """ Displays an info box for the selected session. """
        if not history_item:
            return
        
        try:
            state = history_item.current_state
            spec = state.spec
            systs = state.systs
            meta = spec.meta if spec else {}

            # --- Build an HTML table for alignment ---
            # We can control the font and style here reliably.
            info_html = (
                "<style>"
                "table { border: none; font-size: 13px; }"
                "td { border: none; padding-right: 15px; }"
                "b { font-weight: bold; }"
                "</style>"
                "<table>"
            )

            # Helper for clean table rows
            def add_row(key: str, value, sub_key: bool = False):
                indent = "&nbsp;" * (4 if sub_key else 0)
                return (f"<tr>"
                        f"<td>{indent}<b>{key}:</b></td>"
                        f"<td>{value}</td>"
                        f"</tr>")
            
            def add_header(title: str):
                 return (f"<tr>"
                        f"<td colspan='2'><br><b>{title}:</b></td>"
                        f"</tr>")

            info_html += add_row("Session", history_item.display_name)
            
            info_html += add_header("FITS Header")
            info_html += add_row("Object", meta.get('OBJECT', 'N/A'), sub_key=True)
            info_html += add_row("Instrument", meta.get('INSTRUME', 'N/A'), sub_key=True)
            info_html += add_row("Date Obs", meta.get('DATE-OBS', 'N/A'), sub_key=True)
            
            info_html += add_header("Astrocook Properties")
            info_html += add_row("z_em", f"{spec._data.z_em:.5f}", sub_key=True)
            info_html += add_row("z_rf", f"{spec._data.z_rf:.5f}", sub_key=True)
            
            info_html += add_header("Data Summary")
            info_html += add_row("Data Points", len(spec.x), sub_key=True)
            info_html += add_row("Components", len(systs.components), sub_key=True)
            
            info_html += "</table>"
            # --- End HTML string ---

            # --- Create the custom dialog ---
            dialog = QDialog(self)
            dialog.setWindowTitle(f"Session Info: {history_item.display_name}")
            
            layout = QVBoxLayout(dialog)
            
            text_edit = QTextEdit()
            text_edit.setReadOnly(True)
            text_edit.setHtml(info_html) # Set the HTML content
            
            button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
            button_box.rejected.connect(dialog.reject)
            
            layout.addWidget(text_edit)
            layout.addWidget(button_box)
            
            # --- Set the size you want ---
            dialog.resize(450, 350)
            dialog.exec() # Show the modal dialog

        except Exception as e:
            logging.error(f"Could not generate session info: {e}")
            QMessageBox.warning(self, "Error", f"Could not retrieve session info:\n{e}")

    def _on_session_name_changed(self, index_top_left, index_bottom_right):
        """ Slot called when data in the session_model changes (i.e., rename). """
        
        row = index_top_left.row()
        if not (0 <= row < len(self.session_histories)):
            return
            
        try:
            # 1. Get the history item
            history_item = self.session_histories[row]
            
            # 2. Get the new name from the model
            new_name = self.session_model.data(index_top_left, Qt.ItemDataRole.DisplayRole)
            
            if new_name and new_name != history_item.display_name:
                old_name = history_item.display_name
                
                # 3. Update the history object's display name
                history_item.display_name = new_name
                
                # 4. If this is the active session, update the log scripter title
                if history_item is self.active_history and self.log_scripter_dialog:
                    self.log_scripter_dialog.setWindowTitle(f"Log Scripter: {new_name}")

                    self._update_window_title(new_name)
                    
                logging.info(f"Renamed session '{old_name}' to '{new_name}'")
        
        except AttributeError as e:
            # This is the error you were seeing
            logging.error(f"Error during in-place session rename: {e}")
            QMessageBox.critical(self, "Rename Error", 
                                 f"Failed to set display_name property.\n"
                                 f"Please check 'session_manager.py'.\nError: {e}")
        except Exception as e:
            logging.error(f"Error during in-place session rename: {e}")

    def _on_view_system_inspector(self):
        if not hasattr(self, 'system_inspector') or self.system_inspector is None:
            self.system_inspector = SystemInspector(self)
        
        # Initialize with current session
        if self.active_history:
            self.system_inspector.set_session(self.active_history.current_state)
            
        self.system_inspector.show()
        self.system_inspector.raise_()

    def _on_view_log(self):
        """
        Handles the "View > View Session Log" menu action.
        Shows/raises the *single, persistent* log scripter.
        """
        if self.active_history:
            self._launch_log_scripter(self.active_history)
        else:
            logging.warning("View Log called with no active history.")
    
    def _on_view_identifications(self):
        """
        Handles the "View > View Identifications" menu action.
        """
        if self.active_history:
            self._launch_identification_viewer(self.active_history)
        else:
            logging.warning("View Identifications called with no active history.")

    def _launch_identification_viewer(self, history_object: 'SessionHistory'):
        """
        Creates and shows the *single, persistent* IdentificationViewerDialog.
        """
        if not history_object or not history_object.current_state.spec:
            QMessageBox.warning(self, "No Spectrum", "No spectrum data to analyze.")
            return

        spec_to_view = history_object.current_state.spec
        session_name = history_object.display_name

        # Close the old dialog if it exists (it's read-only, safer to recreate)
        if self.identification_viewer_dialog:
            try:
                self.identification_viewer_dialog.close()
            except Exception:
                pass
            self.identification_viewer_dialog = None

        try:
            dialog = IdentificationViewerDialog(
                spec_to_view, 
                session_name, 
                self
            )
            self.identification_viewer_dialog = dialog
            dialog.finished.connect(
                lambda: self._on_identification_viewer_closed() 
            )
            dialog.show()
            logging.debug(f"Creating new identification viewer for {session_name}")
        except Exception as e:
            logging.error(f"Failed to launch identification viewer: {e}", exc_info=True)

    def _on_identification_viewer_closed(self):
        """
        Slot called when the IdentificationViewerDialog is closed.
        """
        self.identification_viewer_dialog = None
        logging.debug("Identification viewer closed, reference removed.")

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

    # --- DEFINE THE LIST OF RECIPES THAT NEED Z_EM ---
    _RECIPES_REQUIRING_Z_EM = {
        'fit_powerlaw', 
        'estimate_auto', 
        'find_absorbed',
        'extract_preset',
        'identify_lines'
    }

    def _launch_recipe_dialog(self, category, name, initial_params: dict = None):
        if self.active_recipe_dialog:
            self.active_recipe_dialog.activateWindow()
            return
        if not self.active_history or not self.active_history.current_state.spec:
            QMessageBox.warning(self, "No Session", "Please load a spectrum before running a recipe.")
            return

        # 1. Check if this recipe needs z_em
        if name in self._RECIPES_REQUIRING_Z_EM:
            current_z_em = self.active_history.current_state.spec._data.z_em
            
            # 2. Check if z_em is not set
            if current_z_em == 0.0:
                reply = QMessageBox.question(self,
                                             "Emission Redshift Required",
                                             "This recipe requires an emission redshift (z_em) to run.\n\n"
                                             "Do you want to open the 'Set Properties' dialog to set it now?",
                                             QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
                                             QMessageBox.StandardButton.Yes)
                
                if reply == QMessageBox.StandardButton.Yes:
                    # 3. SET THE PENDING ACTION FLAG
                    logging.info(f"Setting '{name}' as pending, launching 'set_properties' first.")
                    self._pending_recipe_on_properties_set = (category, name)
                    
                    # 4. Launch 'set_properties'
                    self._launch_recipe_dialog("edit", "set_properties")
                else:
                    logging.info(f"User cancelled '{name}' due to missing z_em.")
                return # Stop here

        # --- Fix for Zoom/Pan Bug ---
        if self.plot_viewer and self.plot_viewer.toolbar:
            try:
                self.plot_viewer.toolbar.pan(False)
                self.plot_viewer.toolbar.zoom(False)
            except Exception as e:
                logging.warning(f"Could not reset toolbar state: {e}")

        logging.info(f"Launching dialog for recipe: {category}.{name}")
        
        current_state = self.active_history.current_state
        dialog = RecipeDialog(category, name, current_state, self) 

        if initial_params:
            for param_name, value in initial_params.items():
                if param_name in dialog.input_widgets:
                    widget = dialog.input_widgets[param_name]
                    if isinstance(widget, QLineEdit):
                        widget.setText(str(value))
                    # (Add other widget types here if needed later, e.g., QCheckBox)
        
        # --- Connect the new signal to our worker-launching slot ---
        dialog.recipe_requested.connect(self._on_recipe_requested)
        dialog.finished.connect(self._on_recipe_dialog_finished)
        # ---
        
        # Run non-modally so the main window isn't blocked
        dialog.show() 
        self.active_recipe_dialog = dialog
        
        # We no longer care about the .exec() result
        #QTimer.singleShot(0, self._force_restack_floating_widgets)

    def _on_recipe_dialog_finished(self, result: int):
        """
        Slot called when a RecipeDialog is closed (Accepted, Rejected, or 'X').
        This is crucial for cleanup.
        """
        # If the user clicked "Cancel" or "X" (result=Rejected) AND
        # we were waiting for them to set properties, clear the pending recipe.
        if (result == QDialog.Rejected and 
            self._pending_recipe_on_properties_set is not None and
            self.active_recipe_dialog and 
            self.active_recipe_dialog.recipe_name == 'set_properties'):
            
            logging.info("User cancelled 'set_properties', clearing pending recipe.")
            self._pending_recipe_on_properties_set = None

        # We must set this to None so a new dialog can be opened,
        # regardless of how it was closed.
        logging.debug(f"Recipe dialog closed (result: {result}), clearing active dialog lock.")
        self.active_recipe_dialog = None


    def _ask_to_renormalize_model(self, recipe_name: str, params: Dict[str, Any]):
        """
        Checks if a recipe will modify 'cont' and asks the user
        if they want to re-normalize 'model'. Modifies params in-place.
        """
        target_col = params.get('target_col')
        
        # Check for:
        # 1. smooth_column targeting 'cont'
        is_cont_target = (recipe_name == 'smooth_column' and target_col == 'cont')
        # 2. Any of these recipes, which *always* create/modify 'cont'
        is_cont_recipe = recipe_name in ['fit_continuum', 'estimate_auto']

        if (is_cont_target or is_cont_recipe):
            # If it's a cont recipe, check if a model exists
            if self.active_history.current_state.spec.has_aux_column('model'):
                reply = QMessageBox.question(self, 
                                             "Re-normalize Model?",
                                             "This will modify the 'cont' column. Do you also want to re-normalize the 'model' column to this new continuum?",
                                             QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                                             QMessageBox.StandardButton.Yes)
                
                # Set the hidden parameter based on the user's answer
                if reply == QMessageBox.StandardButton.Yes:
                    params['renorm_model'] = 'True'
                else:
                    params['renorm_model'] = 'False'

    # --- *** 4. NEW: Worker-launching slot for single recipes *** ---
    def _on_recipe_requested(self, category: str, recipe_name: str, 
                             params: dict, alias_map: dict):
        """
        Slot called when a RecipeDialog emits its 'recipe_requested' signal.
        Launches the RecipeWorker on a background thread.
        """
        if not self.active_history:
            return # Should not happen

        self._last_attempted_recipe = {
            'category': category,
            'recipe_name': recipe_name,
            # We don't save 'params' here because we want the dialog 
            # to reload fresh defaults or current values, not the failed ones.
            # If you DO want to preserve their typed-in values, it's much harder.
            # Re-opening the dialog standard way is usually sufficient.
        }

        # 1. Close the dialog that sent the signal
        #if self.active_recipe_dialog:
        #    self.active_recipe_dialog.close()
        #    self.active_recipe_dialog = None

        self._ask_to_renormalize_model(recipe_name, params)

        # 2. Show the *same* progress dialog as "Run All"
        if self.progress_dialog:
            self.progress_dialog.cancel()
        self.progress_dialog = QProgressDialog(f"Running {recipe_name}...", "Cancel", 0, 0, self)
        self.progress_dialog.setWindowTitle("Recipe Running")
        self.progress_dialog.setCancelButton(None) 
        self.progress_dialog.setModal(True)
        self.progress_dialog.setFixedWidth(400) # Fix 2
        self.progress_dialog.show()

        # 3. Set up and run the worker
        worker = RecipeWorker(
            session=self.active_history.current_state,
            category=category,
            recipe_name=recipe_name,
            params=params,
            alias_map=alias_map
        )
        
        # Connect to the *same* slots as the ScriptWorker
        worker.signals.finished.connect(self._on_recipe_finished)
        worker.signals.error.connect(self._on_recipe_error)
        
        self.thread_pool.start(worker)

    # --- *** 5. NEW: Callback slots for single recipes *** ---
    def _on_recipe_finished(self, result_data: tuple):
        """Slot called when the RecipeWorker succeeds."""
        self._safely_close_progress_dialog()

        # Set the "Wait" cursor for the whole application
        QApplication.setOverrideCursor(Qt.WaitCursor)
        
        # Schedule the heavy GUI update to run in 10ms.
        # This gives the GUI time to breathe and show the cursor
        # before it freezes for the plot redraw.
        QTimer.singleShot(10, lambda: self._process_recipe_result(result_data))

    def _process_recipe_result(self, result_data: tuple):
        """The actual GUI update, now called by a QTimer."""
        try:
            new_session_state, recipe_name, params = result_data
            
            # 1. Log the successful action
            try:
                if isinstance(self.active_history.log_manager, HistoryLogV2):
                    params_to_log = params.copy()
                    params_to_log.pop('alias_map', None) 
                    self.active_history.log_manager.add_entry(
                        recipe_name=recipe_name, 
                        params=params_to_log
                    )
                    logging.debug(f"Logged successful recipe: {recipe_name}")
            except Exception as e:
                logging.error(f"Failed to log successful recipe: {e}")

            auto_show_col = None
            # Special handling for identify_lines
            if recipe_name == 'identify_lines':
                # Show an info dialog instead of auto-plotting the mask
                try:
                    num_identified = new_session_state.spec.meta.get('num_regions_identified', 0)
                    total_regions = new_session_state.spec.meta.get('num_regions_merged', 0)
                    
                    if total_regions == 0:
                        # Fallback just in case, though num_merged should always be >= 0
                        total_regions = new_session_state.spec.meta.get('num_regions_raw', 0)
                        
                    QMessageBox.information(
                        self, 
                        "Identification Complete", 
                        f"Found identifications for <b>{num_identified}</b> of <b>{total_regions}</b> total absorption regions.\n\n"
                        "Hover over regions for details or use 'View > View Identifications' for a full list."
                    )
                except Exception as e:
                    logging.warning(f"Could not show identification summary: {e}")
                
                auto_show_col = None # Explicitly prevent auto-plotting
            
            # Original logic for all other recipes
            else:
                if self.active_history:
                    old_spec = self.active_history.current_state.spec
                    if old_spec:
                        old_cols = set(old_spec._data.aux_cols.keys())
                        new_cols = set(new_session_state.spec._data.aux_cols.keys())
                        added_cols = new_cols - old_cols

                        if added_cols:
                            # Pick one to show. Prefer 'cont_pl', 'model', 'cont' in that order,
                            # otherwise just pick an arbitrary new one.
                            preferred_order = ['cont_pl', 'model', 'cont']
                            for pref in preferred_order:
                                if pref in added_cols:
                                    auto_show_col = pref
                                    break
                            else:
                                # If none of the preferred ones are new, just take the first one
                                auto_show_col = list(added_cols)[0]

                            logging.info(f"Recipe '{recipe_name}' added column '{auto_show_col}', auto-displaying.")

            # List recipes that drastically change vertical scale
            RECIPES_REQUIRING_AUTOSCALE = {'calibrate_from_magnitudes'}
            should_autoscale = recipe_name in RECIPES_REQUIRING_AUTOSCALE

            # 2. Update the GUI state
            original_history_index = self.session_histories.index(self.active_history)
            branching = is_branching_recipe(recipe_name) #or recipe_name == 'resample'
            
            self.update_gui_session_state(
                new_session_state,
                original_session_index=original_history_index, 
                is_branching=branching,
                auto_show_aux=auto_show_col,
                force_autoscale=should_autoscale
            )
            
            # 3. Check if a recipe was pending on this action
            if recipe_name == 'set_properties' and self._pending_recipe_on_properties_set:
                pending_category, pending_name = self._pending_recipe_on_properties_set
                self._pending_recipe_on_properties_set = None # Clear the flag
                
                # Check if user *actually* set z_em
                new_z_em = new_session_state.spec._data.z_em
                if new_z_em != 0.0:
                    logging.info(f"'set_properties' complete, now launching pending recipe: {pending_name}")
                    # Use a QTimer to launch the dialog in the next event loop
                    QTimer.singleShot(0, lambda: self._launch_recipe_dialog(pending_category, pending_name))
                else:
                    logging.warning(f"User ran 'set_properties' but left z_em=0.0. Aborting pending recipe.")

        except Exception as e:
            logging.error(f"Failed to process recipe result: {e}", exc_info=True)
            QMessageBox.critical(self, "GUI Error", f"Failed to update GUI after recipe:\n{e}")

        finally:
            # --- *** ALWAYS restore the cursor *** ---
            QApplication.restoreOverrideCursor()

            QTimer.singleShot(0, self._force_restack_floating_widgets)

    def _on_recipe_error(self, error_data: tuple):
        """Slot called when the RecipeWorker fails."""
        title, message, trace = error_data
        
        self._safely_close_progress_dialog()
            
        if trace:
            # Critical bug: just show standard error
            QMessageBox.critical(self, title, message)
        else:
            # User error: Offer to try again
            # --- *** START MODIFICATION *** ---
            msg_box = QMessageBox(self)
            msg_box.setIcon(QMessageBox.Warning)
            msg_box.setWindowTitle(title)
            msg_box.setText(message)
            
            # Add standard buttons + a custom "Try Again"
            try_again_btn = msg_box.addButton("Try Again", QMessageBox.AcceptRole)
            msg_box.addButton(QMessageBox.Cancel)
            msg_box.setDefaultButton(try_again_btn)
            
            msg_box.exec()
            
            if msg_box.clickedButton() == try_again_btn:
                if self._last_attempted_recipe:
                    logging.info("User clicked 'Try Again', re-launching dialog.")
                    # Re-launch the dialog
                    self._launch_recipe_dialog(
                        self._last_attempted_recipe['category'],
                        self._last_attempted_recipe['recipe_name']
                    )
                
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
            dialog = LogScripterDialog(
                log_object_to_view, 
                session_name, 
                self, 
                self.recipe_category_map # Pass the map
            )
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
        """
        Launches the ScriptWorker on a background thread to run the script.
        """
        
        if not self.active_history:
            QMessageBox.warning(self, "No Session", "Cannot run script: No active session.")
            return

        # 1. Ask for confirmation
        reply = QMessageBox.question(
            self,
            "Run Script",
            "This will replace your current session history by re-running this script from the original data.\n\nAre you sure you want to proceed?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel
        )
        if reply != QMessageBox.StandardButton.Yes:
            return

        # --- *** FIX 1: Reset toolbar *** ---
        if self.plot_viewer and self.plot_viewer.toolbar:
            try:
                self.plot_viewer.toolbar.pan(False)
                self.plot_viewer.toolbar.zoom(False)
            except Exception as e:
                logging.warning(f"Could not reset toolbar state: {e}")
        # --- *** END FIX 1 *** ---

        logging.info("Starting script run...")

        # 2. Get the original, pristine state
        initial_state = self.active_history.states[0]
        
        # 3. Set up the progress dialog
        if self.progress_dialog:
            self.progress_dialog.cancel()
        self.progress_dialog = QProgressDialog("Running script...", "Cancel", 0, 0, self)
        self.progress_dialog.setWindowTitle("Script Running")
        self.progress_dialog.setCancelButton(None) 
        self.progress_dialog.setModal(True)
        # --- *** FIX 2: Set fixed width *** ---
        self.progress_dialog.setFixedWidth(400)
        # --- *** END FIX 2 *** ---
        self.progress_dialog.show()

        # 4. Set up and run the worker
        worker = ScriptWorker(script_text, initial_state, self)
        
        # Connect to *different* slots than the recipe worker
        worker.signals.finished.connect(self._on_script_finished)
        worker.signals.error.connect(self._on_script_error)
        worker.signals.progress.connect(self._on_script_progress)
        
        self.thread_pool.start(worker)

    def _on_script_progress(self, message: str):
        """Slot to update the progress dialog."""
        if self.progress_dialog:
            self.progress_dialog.setLabelText(message)

    def _on_script_error(self, error_data: tuple):
        """Slot called when the ScriptWorker fails."""
        line_num, line, error_msg = error_data
        
        self._safely_close_progress_dialog()
            
        logging.error(f"Failed to run script line {line_num}: {line}\nError: {error_msg}")
        QMessageBox.critical(
            self, 
            "Script Error",
            f"Failed on line {line_num}:\n{line}\n\nError: {error_msg}"
        )
        if self.log_scripter_dialog:
            self.log_scripter_dialog._update_button_state()

    def _on_script_finished(self, new_history: SessionHistory):
        """Slot called when the ScriptWorker succeeds."""
        self._safely_close_progress_dialog()
        
        QApplication.setOverrideCursor(Qt.WaitCursor)
        # Schedule the heavy GUI update
        QTimer.singleShot(10, lambda: self._process_script_result(new_history))
            
    def _process_script_result(self, new_history: SessionHistory):
        """The actual GUI update for 'Run All', called by a QTimer."""
        # 1. Swap the old history with the new one
        try:
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

            # 2. Refresh everything
            logging.info("Script run completed successfully.")
            self._update_view_for_session(self.active_history.current_state, set_current_list_item=True)
            self._update_undo_redo_actions()
            if self.log_scripter_dialog:
                self.log_scripter_dialog.set_log_object(new_history.log_manager)
        finally:
            # --- *** ALWAYS restore the cursor *** ---
            QApplication.restoreOverrideCursor()

    def _launch_x_convert_dialog(self):
        # Placeholder function called by the QAction
        print("Launching X Convert Dialog...")
        
    def _launch_y_convert_dialog(self):
        # Placeholder function called by the QAction
        print("Launching Y Convert Dialog...")
        
    def _launch_rebin_dialog(self):
        # Placeholder function called by the QAction
        print("Launching Rebin Dialog...")

    def launch_split_from_region(self, xmin: float, xmax: float):
        """
        Launches the 'split' recipe dialog pre-filled with a region expression.
        """
        # 1. Turn off the selection mode in the plot widget
        if self.plot_viewer:
             self.plot_viewer.toggle_region_selector()
             
        # 2. Build the expression string
        # Using the current x-axis unit label from the plot for clarity,
        # although the data 'x' is always native (likely nm).
        # Assuming xmin/xmax passed here are already in native units from the plot.
        expression = f"(x > {xmin:.4f}) & (x < {xmax:.4f})"
        
        # 3. Launch the dialog with these params pre-filled
        self._launch_recipe_dialog(
            "edit", 
            "split", 
            initial_params={"expression": expression} # We need to support this!
        )

    def _update_ui_state(self, is_valid_session, is_startup=False):
        is_valid_session = bool(is_valid_session)

        # ... (Enable/Disable File menu actions) ...
        if hasattr(self, 'save_action'): self.save_action.setEnabled(is_valid_session)
        if hasattr(self, 'close_session_action'): self.close_session_action.setEnabled(is_valid_session)

        # ... (Enable/Disable View menu actions) ...
        if hasattr(self, 'toggle_left_action'): self.toggle_left_action.setEnabled(is_valid_session)
        if hasattr(self, 'toggle_right_action'): self.toggle_right_action.setEnabled(is_valid_session)
        if hasattr(self, 'view_system_inspector_action'): self.view_system_inspector_action.setEnabled(is_valid_session)
        if hasattr(self, 'view_log_action'): self.view_log_action.setEnabled(is_valid_session)
        has_ids = False
        if is_valid_session and self.active_history.current_state.spec:
            # Check if the key exists AND is not None
            has_ids = self.active_history.current_state.spec.meta.get('region_identifications') is not None
        if hasattr(self, 'view_identifications_action'): 
            self.view_identifications_action.setEnabled(has_ids)

        # ** Enable/Disable Recipe Actions **
        enable_recipes = is_valid_session
        # Check if actions exist before enabling/disabling
        if hasattr(self, 'set_properties_action'): self.set_properties_action.setEnabled(enable_recipes)
        if hasattr(self, 'apply_expression_action'): self.apply_expression_action.setEnabled(enable_recipes)
        if hasattr(self, 'mask_expression_action'): self.mask_expression_action.setEnabled(enable_recipes)
        if hasattr(self, 'smooth_column_action'): self.smooth_column_action.setEnabled(enable_recipes)
        if hasattr(self, 'split_action'): self.split_action.setEnabled(enable_recipes)
        if hasattr(self, 'extract_preset_action'): self.extract_preset_action.setEnabled(enable_recipes)
        
        if hasattr(self, 'calculate_running_std_action'): self.calculate_running_std_action.setEnabled(enable_recipes)
        if hasattr(self, 'smooth_action'): self.smooth_action.setEnabled(enable_recipes)
        if hasattr(self, 'rebin_action'): self.rebin_action.setEnabled(enable_recipes)
        if hasattr(self, 'resample_action'): self.resample_action.setEnabled(enable_recipes)
        if hasattr(self, 'calibrate_action'): self.calibrate_action.setEnabled(enable_recipes)

        if hasattr(self, 'auto_cont_action'): self.auto_cont_action.setEnabled(enable_recipes)
        if hasattr(self, 'find_abs_action'): self.find_abs_action.setEnabled(enable_recipes)
        if hasattr(self, 'fit_cont_action'): self.fit_cont_action.setEnabled(enable_recipes)
        if hasattr(self, 'fit_powerlaw_action'): self.fit_powerlaw_action.setEnabled(enable_recipes)

        if hasattr(self, 'identify_lines_action'): self.identify_lines_action.setEnabled(enable_recipes)
        # --- *** END NEW RECIPES *** ---
        
        # ... enable/disable other recipe actions ...

        # Enable Save action only if valid session
        if hasattr(self, 'save_action'): self.save_action.setEnabled(is_valid_session)

        if is_valid_session:
            # --- State when a valid session IS loaded ---
            if self.central_stack.currentIndex() != 0: # If switching from empty
                self.central_stack.setCurrentIndex(0)
            
            """
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
            """

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

    def _safely_close_progress_dialog(self):
        """
        Helper to close the progress dialog safely, deferring the action
        to avoid macOS 'modalSession exited prematurely' warnings.
        """
        if self.progress_dialog:
            # 1. Grab a local reference and clear the main attribute immediately
            dlg = self.progress_dialog
            self.progress_dialog = None
            
            # 2. Defer the actual .close() call to the next event loop iteration
            QTimer.singleShot(0, dlg.close)

    def _on_debug_plot_requested(self, plot_data: dict):
        """
        SLOT: Receives data from a worker thread and plots it in
        a new Matplotlib window. Runs on the main GUI thread.
        """
        try:
            logging.info(f"Received debug plot request: {plot_data['title']}")
            
            # Create a new, separate figure
            fig, ax = plt.subplots()
            
            # Plot the data
            ax.plot(plot_data['v_compare'], plot_data['Y_data'], 
                    label='Y_data (Blue Line AOD)', drawstyle='steps-mid')
            ax.plot(plot_data['v_compare'], plot_data['Y_model'], 
                    label='Y_model (Red Line AOD * f_ratio)', linestyle='--', drawstyle='steps-mid')
            
            # Format the plot
            ax.set_title(plot_data['title'])
            ax.set_xlabel("Velocity (km/s) [relative to primary line]")
            ax.set_ylabel("Mean-Subtracted AOD")
            ax.legend()
            ax.grid(True, linestyle=':')
            
            # Show the plot (non-blocking)
            plt.show()

        except Exception as e:
            logging.error(f"Failed to create debug plot: {e}", exc_info=True)

    def open_session_from_path(self, file_path: str):
        """
        Public method to load a session directly from a path string.
        Used by macOS FileOpen events and Drag & Drop.
        """
        if not os.path.exists(file_path):
            logging.error(f"Cannot open file: {file_path} does not exist.")
            return

        logging.info(f"Direct load requested for: {file_path}")
        
        format_name = 'generic_spectrum'
        session_name = os.path.splitext(os.path.basename(file_path))[0]
        gui_context = self.mock_gui_context

        try:
            # Mostra cursore di attesa
            QApplication.setOverrideCursor(Qt.WaitCursor)
            
            new_session = load_session_from_file(
                archive_path=file_path, 
                name=session_name, 
                format_name=format_name,
                gui_context=gui_context
            )

            if new_session and new_session != 0:
                self.add_session(new_session)
                # Se la finestra era in stato "vuoto", aggiorniamo UI
                self._update_ui_state(is_valid_session=True)
            else:
                logging.error("Session load returned failure.")
                QMessageBox.warning(self, "Load Error", "Failed to load session.")

        except Exception as e:
            logging.error(f"Failed to load file direct: {e}")
            QMessageBox.critical(self, "Error Loading File", f"Could not load {file_path}:\n{e}")
        finally:
            QApplication.restoreOverrideCursor()
            # Forza il focus sulla finestra (utile su macOS se l'app era in background)
            self.activateWindow()
            self.raise_()