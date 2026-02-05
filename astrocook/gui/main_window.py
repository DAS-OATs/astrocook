from copy import deepcopy
import dataclasses
import json
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
from PySide6.QtCore import (
    Qt, QObject, Signal,
    QItemSelectionModel,
    QLocale,
    QPropertyAnimation, QEasingCurve, QRect, QPoint,
    QParallelAnimationGroup,
    QSize, QStringListModel, QThreadPool, QTimer, QUrl
)
from PySide6.QtGui import QAction, QDoubleValidator, QKeySequence, QIcon, QPixmap, QDesktopServices
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog, QInputDialog,
    QHBoxLayout, QMainWindow, QWidget, QVBoxLayout, QGridLayout, QFormLayout, QLabel, QLineEdit, QListView, 
    QMenu, QMessageBox, QSlider, QGroupBox, QProgressBar,
    QPushButton, QProgressDialog, QSizePolicy, QSpacerItem, QStackedWidget, QStyle, QTextEdit,
)
import re
from typing import Any, Dict, List, Optional, Union

from .debug_utils import GLOBAL_PLOTTER
from .identification_viewer_dialog import IdentificationViewerDialog
from .log_scripter_dialog import LogScripterDialog
from astrocook.core.photometry import STANDARD_FILTERS
from .pyside_plot import SpectrumPlotWidget
from .qt_workers import RecipeWorker, ScriptWorker
from .recipe_dialog import RecipeDialog
from .script_dialog import ScriptDialog
from astrocook.core.session_manager import SessionHistory
from astrocook.core.session import SessionV2, load_session_from_file, LogManager
from astrocook.core.structures import HistoryLogV2, V1LogArtifact
from .system_inspector import SystemInspector
from astrocook.core.utils import guarded_deepcopy_v1_state, get_recipe_schema, is_branching_recipe, resource_path # Import recipe helpers
from astrocook.legacy.gui_log import GUILog
from astrocook.legacy.defaults import Defaults
try:
    from astrocook.legacy.functions import trans_parse
    from astrocook.legacy.vars import xem_d
    V1_FUNCTIONS_AVAILABLE = True
except ImportError:
    logging.error("Could not import V1 functions (trans_parse, xem_d) needed for redshift cursor.")
    V1_FUNCTIONS_AVAILABLE = False

# --- *** RECIPE LOOKUP MAP *** ---
from astrocook.recipes.absorbers import ABSORBERS_RECIPES_SCHEMAS
from astrocook.recipes.continuum import CONTINUUM_RECIPES_SCHEMAS
from astrocook.recipes.edit import EDIT_RECIPES_SCHEMAS
from astrocook.recipes.flux import FLUX_RECIPES_SCHEMAS
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
    def __init__(self, initial_sessions: Union[SessionV2, List[SessionV2], None], initial_log_object: Optional[LogManager]):
        super().__init__()
        # Forces the app name to be "Astrocook" instead of "Python"
        app = QApplication.instance()
        if app:
            app.setApplicationName("Astrocook")
            app.setApplicationDisplayName("Astrocook")

        self.session_histories: List[SessionHistory] = [] # List of history managers
        self.active_history: Optional[SessionHistory] = None # Reference to the selected manager
        self.session_model = QStringListModel()
        
        self.log_scripter_dialog: Optional[LogScripterDialog] = None
        self.identification_viewer_dialog: Optional[IdentificationViewerDialog] = None # <<< *** ADD ***
        self.active_recipe_dialog: Optional[QDialog] = None
        self.script_dialog: Optional[ScriptDialog] = None

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

        self._suppress_refit_warning = True  # The warnings are annoying; can be re-enabled for debugging

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
        initial_session = initial_sessions[0] if isinstance(initial_sessions, List) and len(initial_sessions) > 0 else initial_sessions
        valid_session = initial_session if isinstance(initial_session, SessionV2) else None
        
        self._setup_plot_view(valid_session)
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
        # Normalize input to a list
        if isinstance(initial_sessions, SessionV2):
            init_list = [initial_sessions]
        elif isinstance(initial_sessions, list):
            init_list = initial_sessions
        else:
            init_list = []
        
        # ** Create the Mock GUI context once for all loggers **
        self.mock_gui_context = MockV1GUIContext()

        # Iterate and add all passed sessions
        for i, session in enumerate(init_list):
            if isinstance(session, SessionV2) and session.spec and len(session.spec.x) > 0:
                
                # Per previous logic: Always create new history/log for new sessions
                log_object = HistoryLogV2()
                new_history = SessionHistory(session, log_object)
                
                # Link real GUI
                session._gui = self
                session.log = GUILog(self.mock_gui_context)
                session.defs = Defaults(self.mock_gui_context)
                
                self.session_histories.append(new_history)

        # Set active history (default to the last one loaded, or the first?)
        # Usually user wants to see the first one, or the last one. Let's pick the first.
        if self.session_histories:
            self.active_history = self.session_histories[0]
            self.session_model.setStringList([h.display_name for h in self.session_histories])
            self._update_view_for_session(self.active_history.current_state, set_current_list_item=True, is_startup=True)
        else:
            self.active_history = None
            self._update_view_for_session(None, set_current_list_item=False, is_startup=True)

        # Enable drag and drop of files onto the main window
        self.setAcceptDrops(True)

        self._update_undo_redo_actions()
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
        label = QLabel("Welcome to Astrocook V2!\n\nUse 'File > Open Session...'")
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
        
        # ENABLE MULTI-SELECTION
        self.session_list_view.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        
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
        sidebar_layout.setContentsMargins(15, 10, 15, 10)
        sidebar_layout.setSpacing(15)

        # --- Helper: Expand Field ---
        def _expand_field(widget):
            widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            return widget

        # --- Plot Element Toggles ---
        # (This section remains unchanged)
        plot_toggles_layout = QVBoxLayout()
        plot_toggles_layout.setSpacing(8)
        plot_toggles_layout.addWidget(QLabel("<b>Plot Elements:</b>"))
        
        self.error_checkbox = QCheckBox("1-sigma error"); self.error_checkbox.setChecked(True)
        self.error_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.error_checkbox)

        # --- Error Style Toggle ---
        self.error_line_checkbox = QCheckBox("Show error as line")
        self.error_line_checkbox.setToolTip("If checked, plots the error value (dy) as a line instead of a band.")
        self.error_line_checkbox.setChecked(True) # Default to Band
        # Indent to show hierarchy
        self.error_line_checkbox.setStyleSheet("margin-left: 20px;")
        
        # Connect to main replot
        self.error_line_checkbox.stateChanged.connect(self._trigger_replot)
        # Also connect to System Inspector update if open
        self.error_line_checkbox.stateChanged.connect(
            lambda: self.system_inspector.vel_plot._update_plot() 
            if hasattr(self, 'system_inspector') and self.system_inspector and self.system_inspector.isVisible() 
            else None
        )
        
        plot_toggles_layout.addWidget(self.error_line_checkbox)

        self.continuum_checkbox = QCheckBox("Normalization"); self.continuum_checkbox.setChecked(True)
        self.continuum_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.continuum_checkbox)

        self.model_checkbox = QCheckBox("Absorption model"); self.model_checkbox.setChecked(True)
        self.model_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.model_checkbox)

        self.systems_checkbox = QCheckBox("Systems"); self.systems_checkbox.setChecked(True)
        self.systems_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.systems_checkbox)

        # Aux Column
        form_layout_aux = QFormLayout()
        form_layout_aux.setSpacing(5)
        self.aux_col_combo = QComboBox()
        self.aux_col_combo.setToolTip("Select an auxiliary column to plot")
        self.aux_col_combo.currentTextChanged.connect(self._trigger_replot)
        form_layout_aux.addRow("Aux. Column:", self.aux_col_combo) 
        plot_toggles_layout.addLayout(form_layout_aux)

        self.strong_lines_checkbox = QCheckBox("Strong emission lines only")
        self.strong_lines_checkbox.setToolTip("Show only major emission lines when z_em is active.")
        self.strong_lines_checkbox.setChecked(True)
        self.strong_lines_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.strong_lines_checkbox)

        sidebar_layout.addLayout(plot_toggles_layout)

        self.isomag_checkbox = QCheckBox("Show Iso-Mag Grid")
        self.isomag_checkbox.setToolTip("Show lines of constant AB magnitude")
        self.isomag_checkbox.toggled.connect(lambda b: self.plot_viewer.toggle_isomag_grid(b))
        self.isomag_checkbox.setObjectName("PlotControlCheckbox")
        plot_toggles_layout.addWidget(self.isomag_checkbox)
    
        # --- ** Axis & View Controls ** ---
        view_layout = QVBoxLayout()
        view_layout.setSpacing(8)
        view_layout.addWidget(QLabel("<b>View Controls:</b>"))

        # X-Axis Unit
        form_layout_x = QFormLayout()
        self.x_unit_combo = QComboBox()
        self.x_unit_options = ["nm", "Angstrom", "micron"]
        self.x_unit_combo.addItems(self.x_unit_options)
        self.x_unit_combo.setToolTip("Change λ display units (data remains nm).")
        form_layout_x.addRow("λ Unit:", self.x_unit_combo)
        view_layout.addLayout(form_layout_x)
        
        # X/Y Toggles
        self.norm_y_checkbox = QCheckBox("Normalize F")
        self.norm_y_checkbox.setToolTip("Plot F / Continuum and set F limits.")
        self.norm_y_checkbox.toggled.connect(self._on_view_toggle)
        view_layout.addWidget(self.norm_y_checkbox)

        self.snr_checkbox = QCheckBox("Show SNR (F / error)")
        self.snr_checkbox.setToolTip("Plot F / Error Column.")
        self.snr_checkbox.toggled.connect(self._on_view_toggle)
        view_layout.addWidget(self.snr_checkbox)
        
        form_layout_snr = QFormLayout()
        form_layout_snr.setSpacing(5)
        self.snr_col_combo = QComboBox()
        self.snr_col_combo.setToolTip("Select column to use as error for SNR")
        self.snr_col_combo.currentTextChanged.connect(self._trigger_replot)
        form_layout_snr.addRow("  Error Col.:", self.snr_col_combo)
        view_layout.addLayout(form_layout_snr)

        self.log_x_checkbox = QCheckBox("Logarithmic λ")
        view_layout.addWidget(self.log_x_checkbox)
        self.log_y_checkbox = QCheckBox("Logarithmic F")
        
        sidebar_layout.addLayout(view_layout)

        self.x_unit_combo.currentTextChanged.connect(self._trigger_replot)
        self.log_x_checkbox.stateChanged.connect(self._trigger_replot)
        self.log_y_checkbox.stateChanged.connect(self._trigger_replot)

        # --- *** UNIFIED GRID FOR CONTROLS *** ---
        # [CHANGE] Use a SINGLE grid for both sections to guarantee alignment
        unified_grid = QGridLayout()
        unified_grid.setContentsMargins(0, 0, 0, 0)
        unified_grid.setVerticalSpacing(5)
        unified_grid.setHorizontalSpacing(10) # A bit more space between label and box
        
        # Column 0: Labels (Fixed min width)
        unified_grid.setColumnMinimumWidth(0, 80)
        # Column 1: Boxes (Stretch)
        unified_grid.setColumnStretch(1, 1)

        # -- SECTION 1: Axis Limits --
        # Header (Span 2 columns)
        header_limits = QLabel("<b>Axis Limits:</b>")
        # Add some top margin to separate from previous section
        header_limits.setStyleSheet("margin-top: 10px; margin-bottom: 4px;")
        unified_grid.addWidget(header_limits, 0, 0, 1, 2)

        validator = QDoubleValidator()
        validator.setLocale(QLocale.C)
        validator.setNotation(QDoubleValidator.StandardNotation)

        self.xmin_input = _expand_field(QLineEdit("0.0")); self.xmin_input.setValidator(validator)
        self.xmax_input = _expand_field(QLineEdit("1.0")); self.xmax_input.setValidator(validator)
        self.ymin_input = _expand_field(QLineEdit("0.0")); self.ymin_input.setValidator(validator)
        self.ymax_input = _expand_field(QLineEdit("1.0")); self.ymax_input.setValidator(validator)

        self.xmin_input.editingFinished.connect(self._on_set_custom_limits)
        self.xmax_input.editingFinished.connect(self._on_set_custom_limits)
        self.ymin_input.editingFinished.connect(self._on_set_custom_limits)
        self.ymax_input.editingFinished.connect(self._on_set_custom_limits)
        
        # Rows 1-4
        self.lbl_lmin = QLabel("λ Min:")
        self.lbl_lmax = QLabel("λ Max:")
        self.lbl_fmin = QLabel("F Min:")
        self.lbl_fmax = QLabel("F Max:")

        unified_grid.addWidget(self.lbl_lmin, 1, 0)
        unified_grid.addWidget(self.xmin_input, 1, 1)
        unified_grid.addWidget(self.lbl_lmax, 2, 0)
        unified_grid.addWidget(self.xmax_input, 2, 1)
        unified_grid.addWidget(self.lbl_fmin, 3, 0)
        unified_grid.addWidget(self.ymin_input, 3, 1)
        unified_grid.addWidget(self.lbl_fmax, 4, 0)
        unified_grid.addWidget(self.ymax_input, 4, 1)

        # -- SECTION 2: Redshift Cursor --
        # Header (Span 2 columns)
        header_cursor = QLabel("<b>Redshift Cursor:</b>")
        header_cursor.setStyleSheet("margin-top: 10px; margin-bottom: 4px;")
        unified_grid.addWidget(header_cursor, 5, 0, 1, 2)

        self.cursor_series_input = _expand_field(QLineEdit("Ly_a"))
        self.cursor_z_input = _expand_field(QLineEdit("0.0"))
        z_validator = QDoubleValidator()
        z_validator.setLocale(QLocale.C) 
        z_validator.setNotation(QDoubleValidator.StandardNotation)
        self.cursor_z_input.setValidator(z_validator)

        # Rows 6-7
        unified_grid.addWidget(QLabel("Trans.:"), 6, 0)
        unified_grid.addWidget(self.cursor_series_input, 6, 1)
        unified_grid.addWidget(QLabel("Redshift:"), 7, 0)
        unified_grid.addWidget(self.cursor_z_input, 7, 1)
        
        # Add the Unified Grid to the sidebar
        #sidebar_layout.addLayout(unified_grid)

        self.cursor_show_checkbox = QCheckBox("Show Cursor Lines"); self.cursor_show_checkbox.setChecked(False)
        # Add some margin to the checkbox
        self.cursor_show_checkbox.setStyleSheet("margin-top: 5px;")
        unified_grid.addWidget(self.cursor_show_checkbox, 8, 0, 1, 2)

        self.cursor_series_input.editingFinished.connect(self._update_cursor_and_replot)
        self.cursor_z_input.editingFinished.connect(self._update_cursor_and_replot)
        self.cursor_show_checkbox.stateChanged.connect(self._trigger_replot) 
        self.cursor_series_input.returnPressed.connect(lambda: self.cursor_show_checkbox.setChecked(True))

        # Continuum Tools
        header_continuum = QLabel("<b>Edit Continuum:</b>")
        # Add some top margin to separate from previous section
        header_continuum.setStyleSheet("margin-top: 10px; margin-bottom: 4px;")
        unified_grid.addWidget(header_continuum, 9, 0, 1, 2)

        # ROW 10: Slider (Spanning both columns for better control)
        self.stride_slider = QSlider(Qt.Horizontal)
        self.stride_slider.setMinimum(0)
        self.stride_slider.setMaximum(100)
        self.stride_slider.setEnabled(False)
        self.stride_slider.valueChanged.connect(self._on_stride_change)
        self.stride_slider.sliderPressed.connect(self._on_stride_slider_pressed)
        self.stride_slider.sliderReleased.connect(self._on_stride_slider_released)
        self.stride_slider.setToolTip("Adjust density. Drag to smooth/roughen.")
        self._set_slider_from_stride(500)
        
        unified_grid.addWidget(self.stride_slider, 10, 0, 1, 2) # Span 2 columns

        # ROW 11: Spacing Label + Entry
        self.lbl_spacing = QLabel("Spacing:")
        self.lbl_spacing.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        unified_grid.addWidget(self.lbl_spacing, 11, 0, 1, 1)

        self.stride_entry = QLineEdit()
        stride_validator = QDoubleValidator(0.0, 20000.0, 0)
        stride_validator.setLocale(QLocale.C)
        stride_validator.setNotation(QDoubleValidator.StandardNotation)
        self.stride_entry.setValidator(stride_validator)
        self.stride_entry.setAlignment(Qt.AlignLeft)
        self.stride_entry.setEnabled(False)
        self.stride_entry.setToolTip("Enter spacing in km/s")
        self.stride_entry.editingFinished.connect(self._on_stride_text_changed)
        
        unified_grid.addWidget(self.stride_entry, 11, 1, 1, 1)

        # ROW 12: Buttons (Start/Save and Reset) inside a sub-layout for equal sizing
        button_layout = QHBoxLayout()
        button_layout.setSpacing(10) # Gap between buttons
        button_layout.setContentsMargins(0, 0, 0, 0)

        self.edit_cont_button = QPushButton("Start")
        self.edit_cont_button.clicked.connect(self._toggle_continuum_editor)
        # SizePolicy: Expanding allows them to share available space equally
        self.edit_cont_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        button_layout.addWidget(self.edit_cont_button)

        self.reset_cont_button = QPushButton("Reset")
        self.reset_cont_button.setEnabled(False)
        self.reset_cont_button.setToolTip("Discard changes and reset to default spacing.")
        self.reset_cont_button.clicked.connect(self._on_reset_continuum_clicked)
        self.reset_cont_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        button_layout.addWidget(self.reset_cont_button)

        # Add the sub-layout to the main grid, spanning both columns
        unified_grid.addLayout(button_layout, 12, 0, 1, 2)

        # Add the group to your main sidebar layout
        sidebar_layout.addLayout(unified_grid)

        spacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        sidebar_layout.addItem(spacer)

        self.right_sidebar_widget.setObjectName("PlotControlsContainer")
        self.error_checkbox.setObjectName("PlotControlCheckbox")
        self.error_line_checkbox.setObjectName("PlotControlCheckbox")
        self.continuum_checkbox.setObjectName("PlotControlCheckbox")
        self.model_checkbox.setObjectName("PlotControlCheckbox")
        self.systems_checkbox.setObjectName("PlotControlCheckbox")
        self.aux_col_combo.setObjectName("AuxColumnCombo")
        self.strong_lines_checkbox.setObjectName("PlotControlCheckbox")
        self.isomag_checkbox.setObjectName("PlotControlCheckbox")
        self.x_unit_combo.setObjectName("XUnitCombo")
        self.norm_y_checkbox.setObjectName("PlotControlCheckbox")
        self.snr_checkbox.setObjectName("PlotControlCheckbox")
        self.log_x_checkbox.setObjectName("PlotControlCheckbox")
        self.log_y_checkbox.setObjectName("PlotControlCheckbox")
        self.cursor_show_checkbox.setObjectName("PlotControlCheckbox")
        
        self._populate_snr_combo()

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

        # 3. Help Menu
        help_menu = menu_bar.addMenu("&Help")
        
        # --- File and View operations ---

        # --- File Menu Actions
        open_action = QAction("&Open Session...", self); open_action.setShortcut("Ctrl+O"); open_action.triggered.connect(self._on_open_spectrum)
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

        self.view_script_console_action = QAction("View Script &Console", self)
        self.view_script_console_action.triggered.connect(self._on_view_script_console)
        self.view_script_console_action.setEnabled(False) # Disabled until session loaded
        view_menu.addAction(self.view_script_console_action)

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

        # import_systems Action
        import_sys_action = QAction("&Import Systems...", self)
        import_sys_action.setToolTip("Copy absorption systems from another open session.")
        # Note: We call the 'absorbers' category because that's where the recipe logic lives
        import_sys_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "import_systems"))
        edit_menu.addAction(import_sys_action)
        self.import_systems_action = import_sys_action; self.import_systems_action.setEnabled(False)

        # delete Action
        delete_action = QAction("&Delete Elements...", self)
        delete_action.setToolTip("Delete columns or clear the line list")
        delete_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "delete"))
        edit_menu.addAction(delete_action)
        self.delete_action = delete_action; self.delete_action.setEnabled(False)

        edit_menu.addSeparator()

        # split Action
        split_action = QAction("S&plit Out Region...", self)
        split_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "split"))
        edit_menu.addAction(split_action)
        self.split_action = split_action; split_action.setEnabled(False)

        # extract_preset Action
        extract_preset_action = QAction("Extract &Preset Region...", self)
        extract_preset_action.setToolTip("Extract standard regions (Forest, etc.) based on z_em")
        extract_preset_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "extract_preset"))
        edit_menu.addAction(extract_preset_action)
        self.extract_preset_action = extract_preset_action; self.extract_preset_action.setEnabled(False)

        edit_menu.addSeparator()

        # coadd Action
        coadd_action = QAction("Co-&add Spectra...", self)
        coadd_action.setToolTip("Combine multiple sessions and rebin to a single grid")
        coadd_action.triggered.connect(lambda: self._launch_recipe_dialog(
            "edit", "coadd",
            initial_params={'session_names': self.active_history.display_name if self.active_history else ""}
        ))
        edit_menu.addAction(coadd_action)
        self.coadd_action = coadd_action; self.coadd_action.setEnabled(False)

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

        find_abs_action = QAction("&Find Absorbed Regions...", self)
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

        # --- AUTOMATED PIPELINES (Prominent at top) ---
        doublets_auto_action = QAction("Auto-detect &Doublets...", self)
        doublets_auto_action.setToolTip("Automated pipeline: Identify, populate, and fit metal doublets.")
        doublets_auto_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "doublets_auto"))
        absorbers_menu.addAction(doublets_auto_action)
        self.doublets_auto_action = doublets_auto_action; self.doublets_auto_action.setEnabled(False)
        
        lya_auto_action = QAction("Auto-fit &Ly-alpha Forest...", self)
        lya_auto_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "lya_auto"))
        lya_auto_action.setToolTip("Automated pipeline: Identify, populate, and fit Ly-alpha lines.")
        absorbers_menu.addAction(lya_auto_action)
        self.lya_auto_action = lya_auto_action; self.lya_auto_action.setEnabled(False)

        forest_auto_action = QAction("Auto-fit &Full Lyman Series...", self)
        forest_auto_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "forest_auto"))
        forest_auto_action.setToolTip("Automated pipeline: Identify, populate, and fit the Lyman series.")
        #absorbers_menu.addAction(forest_auto_action)
        self.forest_auto_action = forest_auto_action; self.forest_auto_action.setEnabled(False)

        metals_auto_action = QAction("Check &Associated Metals...", self)
        metals_auto_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "metals_auto"))
        metals_auto_action.setToolTip("Scans all existing systems and attempts to add other ions at the same redshift.")
        #absorbers_menu.addAction(metals_auto_action)
        self.metals_auto_action = metals_auto_action; self.metals_auto_action.setEnabled(False)

        absorbers_menu.addSeparator()

        identify_action = QAction("&Identify Absorption Lines...", self)
        identify_action.setToolTip("Automatically identify absorption regions using correlation signals")
        identify_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "identify_lines"))
        absorbers_menu.addAction(identify_action)
        self.identify_lines_action = identify_action; self.identify_lines_action.setEnabled(False)
        
        refit_all_action = QAction("&Refit All Systems...", self)
        refit_all_action.setToolTip("Optimizes all components in the spectrum by fitting disjoint groups.")
        refit_all_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "refit_all"))
        absorbers_menu.addAction(refit_all_action)
        self.refit_all_action = refit_all_action; self.refit_all_action.setEnabled(False)

        clean_negligible_action = QAction("&Clean Negligible Components...", self)
        clean_negligible_action.setToolTip("Removes components with negligible column densities or broadening.")
        clean_negligible_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "clean_negligible"))
        absorbers_menu.addAction(clean_negligible_action)
        self.clean_negligible_action = clean_negligible_action; self.clean_negligible_action.setEnabled(False)
        self._update_undo_redo_actions()

        # --- HELP MENU ACTIONS ---
        
        # 1. Documentation Action (Keeps the menu alive on macOS)
        doc_action = QAction("&Online Documentation", self)
        doc_action.setShortcut(QKeySequence.HelpContents) # Standard F1 / Cmd+? shortcut
        doc_action.triggered.connect(self._open_docs)
        help_menu.addAction(doc_action)

        help_menu.addSeparator()

        # 2. About Action (Will move to App Menu on macOS)
        about_action = QAction("&About Astrocook...", self)
        about_action.triggered.connect(self._on_about_requested)
        help_menu.addAction(about_action)

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
            QLineEdit:focus {{
                border: 1px solid #296bff; /* Highlight color */
            }}
            
            
            /* Style Form Layout labels */
            QWidget#PlotControlsContainer QFormLayout QLabel {{ /* More specific selector */
                color: {text_color};
                margin-bottom: 0px; /* Override default label margin */
                padding-top: 4px; /* Align better with LineEdit */
                min-width: 80px;
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

    def _show_custom_message(self, title, header, text, icon_name="icon_3d_HR.png", 
                             buttons=QMessageBox.StandardButton.Ok, 
                             default_btn=None, 
                             parent=None,
                             checkbox_text=None): # <--- New argument
        """
        Creates a custom Message Box. 
        If checkbox_text is provided, returns (button, is_checked).
        Otherwise, returns button.
        """
        target_parent = parent if parent else self
        
        msg_box = QMessageBox(target_parent)
        msg_box.setWindowTitle(title)
        
        formatted_text = (
            f"<p style='font-size: 13pt; font-weight: bold; margin: 0px;'>{header}</p>"
            f"<p style='font-size: 12pt; font-weight: normal; margin-top: 8px;'>{text}</p>"
        )
        msg_box.setText(formatted_text)
        
        try:
            icon_path = resource_path(os.path.join("assets", icon_name))
            pixmap = QPixmap(icon_path)
            if not pixmap.isNull():
                pixmap = pixmap.scaled(64, 64, Qt.KeepAspectRatio, Qt.SmoothTransformation)
                msg_box.setIconPixmap(pixmap)
        except Exception:
            msg_box.setIcon(QMessageBox.Icon.Information)

        msg_box.setStandardButtons(buttons)
        if default_btn:
            msg_box.setDefaultButton(default_btn)

        msg_box.setWindowModality(Qt.WindowModal)
        
        # Style with optional checkbox margin
        msg_box.setStyleSheet("""
            QLabel { padding-top: 10px; }
            QCheckBox { margin-top: 10px; font-size: 12pt; }
        """)
        
        # [NEW] Add Checkbox if requested
        cb = None
        if checkbox_text:
            cb = QCheckBox(checkbox_text)
            msg_box.setCheckBox(cb)
        
        ret = msg_box.exec()
        
        if cb:
            return ret, cb.isChecked()
        return ret
     

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

            # If the recipe changed the session name (e.g. set_properties), 
            # update the history wrapper immediately so the Sidebar reflects it.
            if target_history.display_name != new_session.name:
                target_history.display_name = new_session.name

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
        is_valid = False
        if session_state_to_show and session_state_to_show.spec:
             # Check if x exists and has data points
             x_col = session_state_to_show.spec.x
             if hasattr(x_col, 'values'):
                 is_valid = len(x_col.values) > 0
             else:
                 # Fallback if it's a Quantity or Array directly (legacy safety)
                 is_valid = len(x_col) > 0

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
        self._populate_snr_combo()
            
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
        
        # Determine symbol
        is_norm = self.norm_y_checkbox.isChecked()
        sym = "ƒ" if is_norm else "F"
        
        # Update Labels
        self.lbl_fmin.setText(f"{sym} Min:")
        self.lbl_fmax.setText(f"{sym} Max:")

        # 3. Trigger the redraw
        self._trigger_replot()

    # Helper to populate the SNR combo with "dF" instead of "dy"
    # Call this whenever you load a session or refresh columns
    def _populate_snr_combo(self):
        self.snr_col_combo.blockSignals(True)
        current_text = self.snr_col_combo.currentText() 
        self.snr_col_combo.clear()
        
        # 1. Always add the standard error as "dF"
        self.snr_col_combo.addItem("dF") 
        
        # 2. Add filtered aux columns
        if self.active_history and self.active_history.current_state.spec:
            spec = self.active_history.current_state.spec
            
            if hasattr(spec, '_data') and hasattr(spec._data, 'aux_cols'):
                for col in spec._data.aux_cols:
                    # [CHANGE] Strict filtering: only 'running_std' or '*err*'
                    is_error_like = (col == 'running_std') or ('err' in col.lower())
                    
                    if is_error_like and col != 'dy': 
                        self.snr_col_combo.addItem(col)

        # Restore selection if it still exists
        index = self.snr_col_combo.findText(current_text)
        if index >= 0:
            self.snr_col_combo.setCurrentIndex(index)
        else:
            self.snr_col_combo.setCurrentIndex(0) # Default to dF

        self.snr_col_combo.blockSignals(False)

    def _toggle_continuum_editor(self):
        # Safety check
        if not hasattr(self, 'plot_viewer'): return

        is_editing = getattr(self.plot_viewer, '_edit_mode_active', False)

        if not is_editing:
            # --- STARTING ---
            if not self.active_history: return

            # Check if continuum exists
            spec = self.active_history.current_state.spec
            # Check if .cont is None or (if it's a Column) has no values
            has_cont = spec.cont is not None

            if not has_cont:
                ret = self._show_custom_message(
                    title="Continuum Missing",
                    header="No continuum found to edit.",
                    text="You need to estimate a continuum before you can edit it manually.\nDo you want to run 'Auto-estimate Continuum' now?",
                    buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
                    default_btn=QMessageBox.StandardButton.Yes
                )
                
                if ret == QMessageBox.StandardButton.Yes:
                    self._launch_recipe_dialog("continuum", "estimate_auto")
                
                # Stop here regardless of choice (don't start editor on empty data)
                return

            slider_pos = self.stride_slider.value()
            initial_stride = self._get_log_stride(slider_pos)
            self._default_continuum_stride = initial_stride
            self.plot_viewer.start_continuum_edit(initial_stride=initial_stride)

            self.edit_cont_button.setText("Save")
            self.stride_slider.setEnabled(True)
            self.stride_entry.setEnabled(True)
            self.reset_cont_button.setEnabled(True)

            # Update label with units
            current_stride = self._get_log_stride(self.stride_slider.value())
            self.stride_entry.setText(f"{current_stride:.0f}")
            
        else:
            # --- STOPPING & SAVING ---
            
            # 1. Retrieve the manual knots
            knots_x, knots_y = self.plot_viewer.stop_continuum_edit(save=True)

            # UI Reset
            self.edit_cont_button.setText("Start")
            self.edit_cont_button.setStyleSheet("") 
            self.stride_slider.setEnabled(False)
            self.stride_entry.setEnabled(False)
            self.reset_cont_button.setEnabled(False)
            self.stride_entry.clear()

            # Clear the saved default
            if hasattr(self, '_default_continuum_stride'):
                del self._default_continuum_stride

            # 2. PREPARE RECIPE PARAMETERS
            if knots_x is not None and knots_y is not None and self.active_history:
                
                recipe_name = "update_from_knots"
                
                # A. Fetch default params from Schema
                # This ensures we respect the defaults defined in continuum.py
                schema = CONTINUUM_RECIPES_SCHEMAS.get(recipe_name, {})
                params = {p['name']: p['default'] for p in schema.get('params', [])}
                
                # B. Inject our manual data
                # (We keep these hidden in the schema so the dialog won't show huge lists)
                params['knots_x'] = knots_x
                params['knots_y'] = knots_y

                # 3. ASK USER (Using the Unified Helper)
                # This mimics the standard recipe flow: it checks if 'renorm_model' 
                # is in the params and prompts the user if needed.
                if hasattr(self, '_ask_to_renormalize_model'):
                    # Pass the name and the fully populated params dict
                    updated_params = self._ask_to_renormalize_model(recipe_name, params)
                    if updated_params is None: # User Cancelled
                        return
                    params = updated_params
                
                # 4. COMMIT TO BACKEND
                try:
                    current_session = self.active_history.current_state
                    
                    # Call the recipe using the dictionary unpacking
                    # This maps keys 'knots_x', 'knots_y', 'renorm_model' to arguments
                    new_session = current_session.continuum.update_from_knots(**params)
                    
                    self.active_history.add_state(new_session)
                    logging.info("Continuum update saved to history.")
                    
                    # 5. FORCE REFRESH
                    if hasattr(self, 'session_list_view'):
                        selection_model = self.session_list_view.selectionModel()
                        current_index = selection_model.currentIndex()
                        if current_index.isValid():
                            self.session_list_view.clicked.emit(current_index)
                        else:
                            self.update_gui_from_state(new_session)
                            
                except Exception as e:
                    logging.error(f"Failed to save continuum: {e}")
                        
    def _on_stride_change(self):
        # 1. Safety Check: If plot_viewer doesn't exist yet (during startup), do nothing.
        if not hasattr(self, 'plot_viewer') or self.plot_viewer is None:
            return
        
        """Handles logarithmic slider changes."""
        slider_pos = self.stride_slider.value()
        stride = self._get_log_stride(slider_pos)

        # Update the LineEdit text
        # Check focus to avoid overwriting user while they are typing (though slider shouldn't move if they are typing)
        if not self.stride_entry.hasFocus():
            self.stride_entry.setText(f"{stride:.0f}")

        # Update the plot (using your previously fixed method)
        if hasattr(self.plot_viewer, '_edit_mode_active') and self.plot_viewer._edit_mode_active:
            self.plot_viewer.update_continuum_stride(stride)

    def _on_stride_text_changed(self):
        """Updates slider and plot when user edits the spacing text manually."""
        try:
            val = float(self.stride_entry.text())
        except ValueError:
            return

        # 1. Clamp value
        val = max(100.0, min(val, 20000.0))
        
        # 2. Update Plot DIRECTLY (High Precision)
        # We do this manually because the slider (0-100) is too coarse 
        # to represent specific km/s values accurately.
        if hasattr(self, 'plot_viewer') and self.plot_viewer._edit_mode_active:
             self.plot_viewer.update_continuum_stride(val)

        # 3. Sync Slider (Visual Only)
        # We block signals so the slider doesn't fire back and 
        # round your precise '1250' down to '1200' via the log conversion.
        self.stride_slider.blockSignals(True)
        self._set_slider_from_stride(val)
        self.stride_slider.blockSignals(False)

        # 4. Update text formatting and drop focus
        # (Dropping focus is important so the slider works again immediately after)
        self.stride_entry.setText(f"{val:.0f}")
        self.stride_entry.clearFocus()
    
    def _on_stride_slider_pressed(self):
        """Lock the current knot shape before resampling."""
        if hasattr(self, 'plot_viewer'):
            self.plot_viewer.begin_stride_interaction()

    def _on_stride_slider_released(self):
        """Release the lock."""
        if hasattr(self, 'plot_viewer'):
            self.plot_viewer.end_stride_interaction()

    def _on_reset_continuum_clicked(self):
        """Resets spacing to default and reloads knots from the original data."""
        if not hasattr(self, 'plot_viewer') or not self.plot_viewer._edit_mode_active:
            return

        # 1. Retrieve the default stride
        default_stride = getattr(self, '_default_continuum_stride', 500.0)
        
        # 2. Block signals to prevent intermediate updates (the "flicker")
        self.stride_slider.blockSignals(True)
        self.stride_entry.blockSignals(True)
        
        try:
            # 3. Reset UI controls without triggering their slots
            self.stride_entry.setText(f"{default_stride:.0f}")
            self._set_slider_from_stride(default_stride)
        finally:
            # 4. Always unblock signals
            self.stride_slider.blockSignals(False)
            self.stride_entry.blockSignals(False)
        
        # 5. Perform the single, authoritative reset update
        self.plot_viewer.reset_continuum_to_original(default_stride)

    def _get_log_stride(self, slider_value):
        """
        Maps linear slider (0-100) to logarithmic stride in km/s.
        Range: ~100 km/s to ~20,000 km/s.
        """
        min_val = 100.0   # Minimum spacing (High detail)
        max_val = 20000.0 # Maximum spacing (Smooth / Power law)
        
        # Logarithmic mapping
        # Formula: y = min * (max/min)^(x/100)
        stride = min_val * (max_val / min_val) ** (slider_value / 100.0)
        return stride # Returns float (km/s)
    
    def _set_slider_from_stride(self, stride):
        """Inverse mapping: Stride (km/s) -> Slider Position."""
        min_val = 100.0
        max_val = 20000.0
        
        # Clamp values
        if stride < min_val: stride = min_val
        if stride > max_val: stride = max_val

        # Formula: x = 100 * log(y/min) / log(max/min)
        position = 100.0 * np.log(stride / min_val) / np.log(max_val / min_val)
        self.stride_slider.setValue(int(position))

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
        file_names, _ = QFileDialog.getOpenFileNames(
            self, 
            "Open Spectrum File(s)", 
            os.getcwd(),
            "All Supported (*.acs *.acs2 *.fits *.txt *.dat);;Astrocook Sessions (*.acs *.acs2);;FITS Files (*.fits);;Text Files (*.txt *.dat);;All Files (*)"
        )
        
        if file_names:
            # Show a wait cursor if loading multiple files
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                for file_name in file_names:
                    logging.info(f"Opening: {file_name}")
                    
                    # 'auto' format detection is now standard
                    session_name = os.path.splitext(os.path.basename(file_name))[0]
                    gui_context = self.mock_gui_context

                    try:
                        new_session = load_session_from_file(
                            archive_path=file_name, 
                            name=session_name, 
                            format_name='auto', # Use auto detection
                            gui_context=gui_context
                        )

                        if new_session == 0:
                            logging.error(f"Failed to load {file_name}")
                            continue

                        self.add_session(new_session)

                    except Exception as e:
                        logging.error(f"Failed to load file {file_name}: {e}")
                        # Don't show a popup for every single failure in a batch, just log it.
            finally:
                QApplication.restoreOverrideCursor()

    def _on_session_list_context_menu(self, pos: QPoint):
        """
        Handles the right-click on the session list.
        Supports Multi-Selection for Stitching.
        """
        # 1. Get all selected indexes
        indexes = self.session_list_view.selectionModel().selectedIndexes()
        
        if not indexes:
            return

        menu = QMenu(self)

        # --- MULTI-SELECTION CASE ---
        if len(indexes) > 1:
            # Sort rows to determine primary (the topmost selected)
            selected_rows = sorted([i.row() for i in indexes])
            
            if not all(0 <= r < len(self.session_histories) for r in selected_rows):
                return

            primary_idx = selected_rows[0]
            primary_history = self.session_histories[primary_idx]
            
            others_names = []
            for r in selected_rows[1:]:
                others_names.append(self.session_histories[r].display_name)
            
            action_text = f"Stitch {len(indexes)} sessions..."
            stitch_action = QAction(action_text, self)
            stitch_action.setToolTip("Combine selected sessions into a new, separate session entry.")
            
            stitch_action.triggered.connect(
                lambda: self._on_stitch_requested(primary_history, others_names)
            )
            menu.addAction(stitch_action)

            coadd_action = QAction(f"Co-add {len(indexes)} sessions...", self)
            coadd_action.setToolTip("Combine and rebin selected sessions into a new, separate session entry.")
            coadd_action.triggered.connect(
                lambda: self._on_coadd_requested(primary_history, others_names)
            )
            menu.addAction(coadd_action)
            
            menu.addSeparator()
        
            # Close Multiple Sessions
            close_multi_action = QAction(f"Close {len(indexes)} sessions", self)
            close_multi_action.triggered.connect(lambda: self._on_close_multiple_sessions(selected_rows))
            menu.addAction(close_multi_action)

            menu.exec(self.session_list_view.mapToGlobal(pos))
            return

        # --- SINGLE SELECTION CASE ---
        row = indexes[0].row()
        if not (0 <= row < len(self.session_histories)):
            return
            
        history_item = self.session_histories[row]

        # --- *** NEW: Info action *** ---
        info_action = QAction(f"View Info", self)
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

        # Duplicate Session
        duplicate_action = QAction("Duplicate", self)
        duplicate_action.triggered.connect(lambda: self._on_duplicate_session(history_item))
        menu.addAction(duplicate_action)

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

    def _on_stitch_requested(self, primary_history, others_names: List[str]):
        """
        Callback for the Stitch menu action.
        Triggers the recipe on the primary session, passing names of others.
        """
        # We must switch active history to the primary one so the RecipeWorker works on the correct object
        if self.active_history != primary_history:
            self.active_history = primary_history
            self._update_view_for_session(primary_history.current_state, set_current_list_item=True)

        # Call the recipe via the standard pipeline
        # We pass a list of STRINGS (names). The hybrid recipe will resolve them.
        params = {'other_sessions': others_names}
        
        self._on_recipe_requested("edit", "stitch", params, {})

    def _on_coadd_requested(self, primary_history, others_names: List[str]):
        """
        Callback for the Co-add context menu action.
        Pre-fills the 'others_names' in the dialog so the user can adjust grid parameters.
        """
        if self.active_history != primary_history:
            self.active_history = primary_history
            self._update_view_for_session(primary_history.current_state, set_current_list_item=True)

        # 1. Combine primary + others into a full list
        all_names = [primary_history.display_name] + others_names
        
        # 2. Join into a single string for the text input
        all_names_str = ", ".join(all_names)
        
        # 3. Launch dialog with correct parameter key 'session_names'
        self._launch_recipe_dialog(
            "edit", "coadd", initial_params={'session_names': all_names_str}
        )

    def _on_duplicate_session(self, history_item: SessionHistory):
        """Creates a new session entry from the current state of another with unique naming."""
        if not history_item:
            return

        current_state = history_item.current_state

        # 1. Deep copy the session state structures
        new_state = current_state.with_new_spectrum(current_state.spec)
        new_state = new_state.with_new_system_list(current_state.systs)

        # 2. Logic to find a unique name
        base_name = current_state.name
        existing_names = [h.display_name for h in self.session_histories]

        # Check if we are already duplicating a copy
        if "_copy" in base_name:
            base_name = base_name.split("_copy")[0]

        counter = 1
        new_name = f"{base_name}_copy{counter}"
        while new_name in existing_names:
            counter += 1
            new_name = f"{base_name}_copy{counter}"

        new_state.name = new_name

        # 3. Add to the GUI session histories
        self.add_session(new_state)
        logging.info(f"Session '{current_state.name}' duplicated as '{new_name}'.")
    
    def _on_close_multiple_sessions(self, rows: List[int]):
        """Closes multiple sessions, handling dirty checks for each."""
        # We must iterate in reverse order so that removing index 5 
        # doesn't shift index 2 into a different position
        for row in sorted(rows, reverse=True):
            history_to_close = self.session_histories[row]
            # This calls your existing method which handles the Save/Discard/Cancel dialog
            self._on_close_session_requested(history_to_close)
    
    def _on_session_info(self, history_item: SessionHistory):
        """ 
        Displays an editable info/inspector box for the selected session. 
        """
        if not history_item: return
        
        state = history_item.current_state
        spec = state.spec
        systs = state.systs
        meta = spec.meta if spec else {}

        # --- Create Custom Dialog ---
        dialog = QDialog(self)
        dialog.setWindowTitle(f"Session Inspector: {history_item.display_name}")
        dialog.setMinimumWidth(400)
        
        dialog.setStyleSheet("""
            QLineEdit {
                padding: 3px;
                border-radius: 5px;
                background-color: {palette.color(palette.ColorRole.Base).name() if 'palette' in locals() else '#FFFFFF'};
                color: {text_color};
            }
            QLineEdit:focus {{
                border: 1px solid #296bff; /* Highlight color */
            }}
        """)
        
        layout = QVBoxLayout(dialog)
        
        # --- 1. Editable Properties Section ---
        form_layout = QFormLayout()
        form_layout.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        form_layout.setVerticalSpacing(5)

        def _expand_field(widget):
            widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            return widget

        # Name
        self.info_name_edit = _expand_field(QLineEdit(history_item.display_name))
        form_layout.addRow("Session Name:", self.info_name_edit)
        
        # Object
        current_obj = meta.get('OBJECT', '')
        obj_display = str(current_obj) if current_obj else "undefined"
        self.info_object_edit = _expand_field(QLineEdit(obj_display))
        form_layout.addRow("Object Name:", self.info_object_edit)
        
        # z_em
        self.info_zem_edit = _expand_field(QLineEdit(f"{spec._data.z_em:.5f}"))
        dbl_validator = QDoubleValidator(); dbl_validator.setLocale(QLocale.C)
        self.info_zem_edit.setValidator(dbl_validator)
        form_layout.addRow("Emission Redshift (z_em):", self.info_zem_edit)

        # Resolution
        if spec._data.resol > 0:
            current_resol_str = f"{spec._data.resol:.0f}"
        else:
            current_resol_str = "undefined"
            
        self.info_resol_edit = _expand_field(QLineEdit(current_resol_str))
        self.info_resol_edit.setPlaceholderText("Resolving Power R")
        form_layout.addRow("Resolution (R):", self.info_resol_edit)
        
        layout.addLayout(form_layout)
        
        # --- 2. Read-Only Stats Section (HTML) ---
        x_vals = spec.x.value
        x_range = "Empty"
        snr_str = "N/A"
        step_str = "N/A"
        
        if len(x_vals) > 0:
            x_range = f"{np.min(x_vals):.2f} - {np.max(x_vals):.2f} {spec.x.unit}"
            try:
                snr_arr = spec.y.value / spec.dy.value
                valid_snr = snr_arr[(spec.y.value > 0) & (spec.dy.value > 0)]
                if len(valid_snr) > 0:
                    snr_str = f"{np.nanmedian(valid_snr):.2f}"
            except: pass

            # --- Step Detection Logic ---
            if len(x_vals) > 1:
                dx = np.diff(x_vals)
                med_dx = np.median(dx)
                
                # 1. Check Linear Constancy (within 0.01%)
                is_linear = np.allclose(dx, med_dx, rtol=1e-4)
                
                if is_linear:
                    step_str = f"{med_dx:.4f} {spec.x.unit} (constant)"
                else:
                    # 2. Check Log Constancy
                    is_log = False
                    if np.all(x_vals > 0):
                        dlog = np.diff(np.log10(x_vals))
                        med_dlog = np.median(dlog)
                        is_log = np.allclose(dlog, med_dlog, rtol=1e-4)
                    
                    if is_log:
                        # Convert to Velocity Step: dv = c * ln(10) * dlog10(lambda)
                        c_kms = 299792.458
                        dv = c_kms * med_dlog * np.log(10)
                        step_str = f"{med_dx:.4f} {spec.x.unit} (log; dv: {dv:.2f} km/s)"
                    else:
                        step_str = f"{med_dx:.4f} {spec.x.unit} (variable)"

        stats_html = (
            "<p style='margin-bottom: 5px;'><b>Other info:</b></p>"
            "<style>td { padding-left: 15px; padding-top: 5px; }</style>"
            "<table width='100%'>"
            f"<tr><td>Wavelength Range:</td><td>{x_range}</td></tr>"
            f"<tr><td>Median Step:</td><td>{step_str}</td></tr>" # Added Row
            f"<tr><td>Median SNR:</td><td>{snr_str}</td></tr>"
            f"<tr><td>Number of data points:</td><td>{len(spec.x)}</td></tr>"
            f"<tr><td>Number of components:</td><td>{len(systs.components)}</td></tr>"
            "</table>"
        )
        
        stats_label = QLabel(stats_html)
        stats_label.setTextFormat(Qt.RichText)
        layout.addWidget(stats_label)
        
        # --- 3. Buttons ---
        btn_box = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        btn_box.accepted.connect(lambda: self._apply_session_info_changes(
            dialog, history_item
        ))
        btn_box.rejected.connect(dialog.reject)
        layout.addWidget(btn_box)
        
        dialog.exec()

    def _apply_session_info_changes(self, dialog, history_item):
        """
        Callback to apply changes from the info dialog via the set_properties recipe.
        """
        # Get values from the widgets
        new_name = self.info_name_edit.text()
        new_obj = self.info_object_edit.text()
        new_zem = self.info_zem_edit.text()
        
        # Handle Object "undefined"
        raw_obj = self.info_object_edit.text().strip()
        if raw_obj.lower() == "undefined" or not raw_obj:
            new_obj = "" # Set to empty string internally
        else:
            new_obj = raw_obj
            
        # Handle z_em safety (replace comma with dot just in case)
        new_zem = self.info_zem_edit.text().replace(',', '.')
        
        # Handle Resolution "undefined"
        raw_resol = self.info_resol_edit.text().strip()
        if raw_resol.lower() == "undefined" or not raw_resol:
            new_resol = "0.0"
        else:
            new_resol = raw_resol

        # 1. Switch active history to the one being edited 
        if self.active_history != history_item:
            self.active_history = history_item
            self._update_view_for_session(history_item.current_state, set_current_list_item=True)

        # 2. Build params dictionary
        params = {
            'name': new_name,
            'object': new_obj,
            'z_em': new_zem,
            'resol': new_resol # [CHANGE] Pass new resolution
        }

        # 3. Call the recipe through the standard pipeline
        self._on_recipe_requested("edit", "set_properties", params, {})
        
        dialog.accept()

    def _on_session_name_changed(self, index_top_left, index_bottom_right):
        """ 
        Slot called when the user renames a session directly in the sidebar list.
        """
        row = index_top_left.row()
        if not (0 <= row < len(self.session_histories)):
            return
            
        try:
            history_item = self.session_histories[row]
            
            # 1. Get the new name from the Qt Model (user input)
            new_name = self.session_model.data(index_top_left, Qt.ItemDataRole.DisplayRole)
            
            # 2. Check against the INTERNAL session state name
            current_internal_name = history_item.current_state.name
            
            if new_name and new_name != current_internal_name:
                logging.info(f"Sidebar rename: '{current_internal_name}' -> '{new_name}'.")
                
                # --- [FIX START] Hybrid Update ---
                
                # A. Update the wrapper IMMEDIATELY.
                # This ensures if you open the Info window 1ms later, it shows the new name.
                history_item.display_name = new_name
                
                # B. Switch active history if needed (recipe runs on active)
                if self.active_history != history_item:
                    self.active_history = history_item
                
                # C. Trigger recipe to sync the internal SessionV2 state and Log.
                params = {
                    'name': new_name,
                    'object': '_current_',
                    'z_em': '_current_',
                    'resol': '_current_'
                }
                self._on_recipe_requested("edit", "set_properties", params, {})
                
                # --- [FIX END] ---

            elif new_name != history_item.display_name:
                # Sync wrapper if model changed but internal state matched (rare edge case)
                history_item.display_name = new_name

        except Exception as e:
            logging.error(f"Error during in-place session rename: {e}")

    def _on_about_requested(self):
        """Launches the About dialog."""
        dlg = AboutDialog(self)
        dlg.exec()

    def _open_docs(self):
        """Opens the online documentation in the default browser."""
        # Replace with your actual documentation URL (e.g., ReadTheDocs or GitHub Wiki)
        url = "https://das-oats.github.io/astrocook/dev/" 
        QDesktopServices.openUrl(QUrl(url))

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
    
    def _on_view_script_console(self):
        """Launches the Script Console dialog."""
        if not self.script_dialog:
            self.script_dialog = ScriptDialog(self)
        
        self.script_dialog.show()
        self.script_dialog.raise_()
        self.script_dialog.activateWindow()

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
            ret = self._show_custom_message(
                title="Close Session",
                header=f"Do you want to save changes to '{session_name}'?",  # Grassetto
                text="Your changes will be lost if you don't save them.",     # Normale
                buttons=QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
                default_btn=QMessageBox.StandardButton.Save
            )

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
        'doublets_auto',
        'lya_auto',
        'forest_auto',
        'metals_auto',
        'identify_lines'
    }

    # Recipes that need Resolution info
    _RECIPES_REQUIRING_RESOL = {
        'fit_component',
        'add_component',
        'doublets_auto',
        'lya_auto',
        'forest_auto',
        'metals_auto',
        'identify_lines'  # Useful for matching linewidths
    }

    def _check_and_handle_requirements(self, category, name, params=None, mode='dialog') -> bool:
        """
        Centralized check for z_em and resolution.
        Returns True if requirements are met.
        If not, prompts user, sets pending state, launches set_properties, and returns False.
        """
        if not self.active_history or not self.active_history.current_state.spec:
            return True # Let validation fail downstream if no session

        spec = self.active_history.current_state.spec
        
        # 1. Check z_em
        if name in self._RECIPES_REQUIRING_Z_EM:
            if spec._data.z_em == 0.0:
                ret = self._show_custom_message(
                    title="Emission Redshift Required",
                    header="Emission Redshift (z_em) is required.",
                    text="This recipe cannot run without z_em.\nDo you want to open 'Set Properties' to set it now?",
                    buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
                    default_btn=QMessageBox.StandardButton.Yes
                )
                if ret == QMessageBox.StandardButton.Yes:
                    self._set_pending_action(category, name, params, mode)
                    self._launch_recipe_dialog("edit", "set_properties")
                return False

        # 2. Check Resolution
        if name in self._RECIPES_REQUIRING_RESOL:
            # Check Scalar R or Meta FWHM
            has_resol = (spec._data.resol > 0) or (spec.meta.get('resol_fwhm', 0.0) > 0)
            
            # Check Column (must be non-empty and positive)
            if not has_resol and spec.has_aux_column('resol'):
                try:
                    col = spec.get_column('resol')
                    if col is not None:
                        val = np.nanmedian(col.value)
                        if np.isfinite(val) and val > 0:
                            has_resol = True
                except Exception:
                    pass
            
            if not has_resol:
                ret = self._show_custom_message(
                    title="Resolution Required",
                    header="Instrument Resolution is missing.",
                    text="This recipe needs resolution (FWHM or R) for accurate fitting.\nDo you want to open 'Set Properties' to set it now?",
                    buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
                    default_btn=QMessageBox.StandardButton.Yes
                )
                if ret == QMessageBox.StandardButton.Yes:
                    self._set_pending_action(category, name, params, mode)
                    self._launch_recipe_dialog("edit", "set_properties")
                return False
                
        return True

    def _set_pending_action(self, category, name, params, mode):
        """Stores the intended action to resume after properties are set."""
        logging.info(f"Setting pending action: {name} ({mode})")
        self._pending_recipe_on_properties_set = {
            'category': category,
            'name': name,
            'params': params,
            'mode': mode
        }

    def _launch_recipe_dialog(self, category, name, initial_params: dict = None):
        if self.active_recipe_dialog:
            self.active_recipe_dialog.activateWindow()
            return
        if not self.active_history: # Simple check
            QMessageBox.warning(self, "No Session", "Please load a spectrum.")
            return
        
        if name == 'set_properties':
            current_params = {}
            # Name
            current_params['name'] = self.active_history.display_name
            
            # Object, z_em, Resol
            spec = self.active_history.current_state.spec
            if spec:
                meta = spec.meta
                current_params['object'] = meta.get('OBJECT', '')
                current_params['z_em'] = f"{spec._data.z_em:.5f}"
                
                if spec._data.resol > 0:
                    current_params['resol'] = f"{spec._data.resol:.0f}"
                else:
                    current_params['resol'] = "0" # Or leave blank
            
            # Merge with any passed initial_params (which take precedence)
            if initial_params:
                current_params.update(initial_params)
            initial_params = current_params

        # [FIX] Use centralized check (mode='dialog')
        if not self._check_and_handle_requirements(category, name, initial_params, mode='dialog'):
            return

        # ... (Existing Fix for Zoom/Pan Bug) ...
        if self.plot_viewer and self.plot_viewer.toolbar:
            try:
                self.plot_viewer.toolbar.pan(False)
                self.plot_viewer.toolbar.zoom(False)
            except Exception: pass

        logging.info(f"Launching dialog for recipe: {category}.{name}")
        
        current_state = self.active_history.current_state
        dialog = RecipeDialog(category, name, current_state, self) 

        if initial_params:
            for param_name, value in initial_params.items():
                if param_name in dialog.input_widgets:
                    widget = dialog.input_widgets[param_name]
                    if isinstance(widget, QLineEdit):
                        widget.setText(str(value))
        
        dialog.recipe_requested.connect(self._on_recipe_requested)
        dialog.finished.connect(self._on_recipe_dialog_finished)
        dialog.show() 
        self.active_recipe_dialog = dialog
    
    def _on_stop_recipe(self):
        """Slot called when Stop button is clicked."""
        if hasattr(self, '_current_worker') and self._current_worker:
            self._current_worker.stop()
            
            # Update UI immediately
            if self.progress_dialog and hasattr(self.progress_dialog, 'stop_btn'):
                self.progress_dialog.stop_btn.setEnabled(False)
                self.progress_dialog.stop_btn.setText("Stopping...")
                self.progress_dialog.status_label.setText("Aborting process...")
        else:
            # If state is weird (no worker), just close
            self._cleanup_progress_ui()

    def _on_recipe_requested(self, category: str, recipe_name: str, 
                             params: dict, alias_map: dict):
        """
        Slot called when a recipe is requested (Dialog OR Context Menu).
        """
        if not self.active_history: return

        # [Check Requirements logic remains the same...]
        if not self._check_and_handle_requirements(category, recipe_name, params, mode='direct'):
            return

        self._last_attempted_recipe = {
            'category': category,
            'recipe_name': recipe_name,
        }

        self._ask_to_renormalize_model(recipe_name, params)

        # Attach Warning Capturing (for ALL recipes)
        self.warning_capture_handler = WarningCaptureHandler()
        # We assume standard formatting is fine, or set a formatter if you wish
        logging.getLogger().addHandler(self.warning_capture_handler)

        # 1. Clean up old dialogs
        if self.progress_dialog: 
            self.progress_dialog.close()
            self.progress_dialog = None
        
        # 2. Decide: Simple Dialog vs. Detailed Dialog
        long_running_recipes = {
            'doublets_auto', 'lya_auto', 'forest_auto', 
            'refit_all', 'metals_auto', 'add_component', 'fit_component'
        }
        
        is_long_running = recipe_name in long_running_recipes

        if recipe_name == 'import_systems':
            # Check if refit is requested (handles both boolean True and string 'True')
            if str(params.get('refit', False)) == 'True':
                is_long_running = True

        if is_long_running:
            # A. Create Detailed Dialog
            self.progress_dialog = RecipeProgressDialog(f"{recipe_name}", self)

            if recipe_name in {'add_component', 'fit_component'}:
                self.progress_dialog.progress_bar.setRange(0, 0) # Indeterminate mode
                self.progress_dialog.progress_bar.setTextVisible(False) # Hide "%" text

            self.progress_dialog.stop_requested.connect(self._on_stop_recipe)
            self.progress_dialog.show()

            # B. Setup Log Tapping
            self.log_bridge = LogSignalBridge()
            self.log_bridge.new_log_message.connect(self.progress_dialog.update_log)
            
            self.qt_log_handler = QtLogHandler(self.log_bridge)
            self.qt_log_handler.setLevel(logging.INFO)
            
            # Attach to Root Logger to capture everything from core/recipes
            logging.getLogger().addHandler(self.qt_log_handler)
            
        else:
            # Standard "Short" Dialog
            self.progress_dialog = QProgressDialog(f"Running {recipe_name}...", "Cancel", 0, 0, self)
            self.progress_dialog.setWindowTitle("Recipe Running")
            self.progress_dialog.setCancelButton(None) 
            self.progress_dialog.setModal(True)
            self.progress_dialog.setFixedWidth(400)
            self.progress_dialog.show()
            
            # No log tapping for short recipes
            self.qt_log_handler = None


        # Start Worker (Existing Code)
        worker = RecipeWorker(
            session=self.active_history.current_state,
            category=category,
            recipe_name=recipe_name,
            params=params,
            alias_map=alias_map
        )
        self._current_worker = worker
        worker.signals.finished.connect(self._on_recipe_finished)
        worker.signals.error.connect(self._on_recipe_error)
        self.thread_pool.start(worker)

    def _process_recipe_result(self, result_data: tuple):
        """The actual GUI update."""
        try:
            new_session_payload, recipe_name, params = result_data
            
            # --- 1. HANDLE MULTI-RETURN (LISTS) ---
            # Normalize payload to a list so we can loop
            if isinstance(new_session_payload, list):
                new_sessions = new_session_payload
            else:
                new_sessions = [new_session_payload]

            # Store the ORIGINAL active history (the parent of these operations)
            # This is critical for "Sibling" branching (A->B, A->C) instead of "Chain" (A->B->C)
            original_parent_history = self.active_history

            # Loop through all returned sessions
            for i, new_session_state in enumerate(new_sessions):
                
                # [SAFETY CHECK] If operation was cancelled
                if getattr(new_session_state, '_stop_flag', False):
                    logging.info("Operation cancelled. Reverting.")
                    if i == len(new_sessions) - 1: # Only update view on last one
                        self._update_view_for_session(self.active_history.current_state)
                    continue

                # --- 2. LOGGING (Only once or per item?) ---
                # Let's log for every item to be safe, or just once if params are shared.
                # Since params are shared, maybe just log once? 
                # Existing logic logs every time this runs.
                try:
                    if i == 0 and isinstance(original_parent_history.log_manager, HistoryLogV2):
                        # Log the command once on the parent
                        params_to_log = params.copy()
                        params_to_log.pop('alias_map', None) 
                        original_parent_history.log_manager.add_entry(recipe_name, params_to_log)
                except Exception: pass

                # --- 3. IDENTIFY LINES LOGIC ---
                auto_show_col = None
                if recipe_name in ['identify_lines', 'doublets_auto', 'populate_from_identification']:
                    if new_session_state.spec.has_aux_column('abs_ids'):
                        auto_show_col = 'abs_ids'
                    elif new_session_state.spec.has_aux_column('abs_mask'):
                        auto_show_col = 'abs_mask'

                elif recipe_name == 'fit_powerlaw':
                    if new_session_state.spec.has_aux_column('cont_pl'):
                        auto_show_col = 'cont_pl'

                should_autoscale = recipe_name in {'calibrate_from_magnitudes'}

                if recipe_name == 'set_properties':

                    if 'resol' in params:
                        raw_val = str(params['resol'])

                        if raw_val != '_current_':
                            try:
                                new_res = float(raw_val)

                                if new_session_state.systs and new_session_state.systs.components:

                                    # 1. Force Update All Components
                                    updated_comps = []
                                    for i, c in enumerate(new_session_state.systs.components):
                                        # Create a fresh copy
                                        new_c = dataclasses.replace(c, resol=new_res)
                                        updated_comps.append(new_c)
                                        # Log the first one to verify

                                    # 2. Rebuild Session
                                    if updated_comps:
                                        # Replace data core
                                        new_syst_data = dataclasses.replace(new_session_state.systs._data, components=updated_comps)

                                        # Re-instantiate SystemListV2
                                        # IMPORTANT: Ensure SystemListV2 is imported or accessible via class
                                        SystemListClass = new_session_state.systs.__class__
                                        new_systs = SystemListClass(new_syst_data)

                                        # Restore constraints if present
                                        if new_session_state.systs.constraint_model:
                                            new_systs.constraint_model = new_session_state.systs.constraint_model

                                        # Update Session
                                        new_session_state = new_session_state.with_new_system_list(new_systs)
                            except Exception as e:
                                logging.error(f"Propagation crashed: {e}", exc_info=True)

                original_history_index = self.session_histories.index(original_parent_history)
                branching = is_branching_recipe(recipe_name)
                
                # --- 4. PREPARE PARENT FOR BRANCHING ---
                # If we are branching multiple times (e.g. extract_preset="lya, red"), 
                # we must ensure the 'active_history' is reset to the PARENT 
                # before creating the next branch.
                if branching and len(new_sessions) > 1:
                    self.active_history = original_parent_history

                # --- 5. UPDATE GUI ---
                self.update_gui_session_state(
                    new_session_state,
                    original_session_index=original_history_index, 
                    is_branching=branching,
                    auto_show_aux=auto_show_col,
                    force_autoscale=should_autoscale
                )

            # After loop, the *last* created session will be active (set by update_gui_session_state).
            # This is usually the desired behavior.

            # --- 6. WARNINGS (Run once at the end) ---
            # (Keep existing warning logic using self.active_history or original params)
            warn_resolution = False
            warn_refit = False

            if recipe_name == 'set_properties' \
                and 'resol' in params \
                and str(params['resol']) != '_current_' \
                and self.active_history.current_state.spec.has_aux_column('model') \
                and not self._suppress_refit_warning:
                warn_resolution = True
            elif recipe_name == 'update_component' and 'resol' in params and not self._suppress_refit_warning:
                 warn_resolution = True
            elif recipe_name in ['update_component', 'update_constraint'] and not self._suppress_refit_warning:
                 warn_refit = True

            if warn_resolution or warn_refit:
                msg_parent = self
                if hasattr(self, 'system_inspector') and self.system_inspector and self.system_inspector.isVisible():
                    msg_parent = self.system_inspector

                title = "Resolution Updated" if warn_resolution else "Parameters Updated"
                header = "Resolution updated internally." if warn_resolution else "Component parameters updated."
                text = "To see the effect on the absorption model (green line), please refit the components."
                
                ret, is_checked = self._show_custom_message(
                    title=title,
                    header=header,
                    text=text,
                    buttons=QMessageBox.StandardButton.Ok,
                    parent=msg_parent,
                    checkbox_text="Don't show this again"
                )
                
                if is_checked:
                    self._suppress_refit_warning = True
                    logging.info("User suppressed refit warnings for this session.")

            # [NEW] Show Captured Warnings (if any)
            if hasattr(self, '_captured_warnings') and self._captured_warnings:
                # Filter out duplicate messages to be clean
                unique_warnings = list(dict.fromkeys(self._captured_warnings))

                # Don't show if it's just the "Recipe Aborted" warning we handle elsewhere
                real_warnings = [w for w in unique_warnings if "Recipe" not in w and "Aborted" not in w]

                if real_warnings:
                    count = len(real_warnings)
                    msg_text = f"The recipe finished, but generated {count} warning(s):<ul>"
                    for w in real_warnings:
                        msg_text += f"<li>{w}</li>"
                    msg_text += "</ul>"

                    self._show_custom_message(
                        title="Warnings Issued",
                        header="Process completed with warnings",
                        text=msg_text,
                        buttons=QMessageBox.StandardButton.Ok,
                        icon_name="icon_3d_HR.png" # or a warning icon
                    )

            # --- Resume Pending Action ---
            if recipe_name == 'set_properties' and self._pending_recipe_on_properties_set:
                pending = self._pending_recipe_on_properties_set
                self._pending_recipe_on_properties_set = None 
                
                logging.info(f"Resuming pending action: {pending['name']}")
                
                if pending['mode'] == 'dialog':
                    QTimer.singleShot(0, lambda: self._launch_recipe_dialog(
                        pending['category'], pending['name'], pending.get('params')
                    ))
                elif pending['mode'] == 'direct':
                    QTimer.singleShot(0, lambda: self._on_recipe_requested(
                        pending['category'], pending['name'], pending['params'], {}
                    ))

        except Exception as e:
            logging.error(f"Failed to process recipe result: {e}", exc_info=True)
            QMessageBox.critical(self, "GUI Error", f"Failed to update GUI:\n{e}")

        finally:
            QApplication.restoreOverrideCursor()
            QTimer.singleShot(0, self._force_restack_floating_widgets)

    def _show_refit_warning_with_checkbox(self, parent_widget):
        """Shows a warning with a 'Don't show again' checkbox, styled like custom messages."""
        msg_box = QMessageBox(parent_widget)
        msg_box.setWindowTitle("Parameters Updated")
        
        # --- 1. Custom Icon (Same as _show_custom_message) ---
        try:
            icon_path = resource_path(os.path.join("assets", "icon_3d_HR.png"))
            pixmap = QPixmap(icon_path)
            if not pixmap.isNull():
                pixmap = pixmap.scaled(64, 64, Qt.KeepAspectRatio, Qt.SmoothTransformation)
                msg_box.setIconPixmap(pixmap)
        except Exception:
            msg_box.setIcon(QMessageBox.Icon.Information)

        # --- 2. HTML Formatting (Same as _show_custom_message) ---
        header = "Component parameters have been updated."
        text = "To see the effect on the absorption model (green line), please <b>refit</b> the components."
        
        formatted_text = (
            f"<p style='font-size: 13pt; font-weight: bold; margin: 0px;'>{header}</p>"
            f"<p style='font-size: 12pt; font-weight: normal; margin-top: 8px;'>{text}</p>"
        )
        msg_box.setText(formatted_text)
        
        # --- 3. CSS Styling (Same as _show_custom_message) ---
        msg_box.setStyleSheet("""
            QLabel { padding-top: 10px; }
            QCheckBox { margin-top: 10px; font-size: 12pt; }
        """)

        # Add Checkbox
        cb = QCheckBox("Don't show this again")
        msg_box.setCheckBox(cb)
        
        msg_box.setStandardButtons(QMessageBox.Ok)
        msg_box.exec()
        
        if cb.isChecked():
            self._suppress_refit_warning = True
            logging.info("User suppressed future refit warnings for this session.")

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
        is_cont_recipe = recipe_name in ['fit_continuum', 'estimate_auto', 'update_from_knots']

        if (is_cont_target or is_cont_recipe):
            # If it's a cont recipe, check if a model exists
            if self.active_history.current_state.spec.has_aux_column('model'):
                ret = self._show_custom_message(
                    title="Re-normalize Model?",
                    header="Update Model normalization?",
                    text="You are modifying the continuum. Do you want to automatically re-normalize the 'model' column to match?", # Normale
                    buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    default_btn=QMessageBox.StandardButton.Yes
                )

                # Set the hidden parameter based on the user's answer
                if ret == QMessageBox.StandardButton.Yes:
                    params['renorm_model'] = 'True'
                else:
                    params['renorm_model'] = 'False'
            
        return params


    # --- *** 5. NEW: Callback slots for single recipes *** ---
    def _cleanup_progress_ui(self):
        # [NEW] Retrieve and Detach Warnings
        self._captured_warnings = []
        if hasattr(self, 'warning_capture_handler') and self.warning_capture_handler:
            self._captured_warnings = self.warning_capture_handler.warnings[:] # Copy list
            logging.getLogger().removeHandler(self.warning_capture_handler)
            self.warning_capture_handler = None

        try:
            # 1. Detach Log Handler
            if hasattr(self, 'qt_log_handler') and self.qt_log_handler:
                # Remove from root logger to be sure
                logging.getLogger().removeHandler(self.qt_log_handler)
                self.qt_log_handler = None
                self.log_bridge = None

            # 2. Close Dialog safely
            self._safely_close_progress_dialog()
            
        except Exception as e:
            logging.error(f"Error cleaning up progress UI: {e}")
            # Fallback: force close if something failed
            if self.progress_dialog:
                self.progress_dialog.close()
                self.progress_dialog = None

    def _on_recipe_finished(self, result_data: tuple):
        """Slot called when the RecipeWorker succeeds."""
        # 1. Cleanup UI immediately
        self._cleanup_progress_ui()
        
        # 2. Restore Cursor (Important!)
        QApplication.restoreOverrideCursor()

        new_session_state, recipe_name, params = result_data
        
        # [NEW] Check if operation was cancelled (returned same object or flagged)
        if getattr(new_session_state, '_stop_flag', False):
            logging.info("Operation cancelled. Reverting to previous state.")
            # Do NOT update history.
            # Just ensure GUI reflects the ACTIVE history state (which hasn't changed)
            self._update_view_for_session(self.active_history.current_state)
            return

        # 3. Schedule the data processing
        # We pass the data to the processor
        QTimer.singleShot(10, lambda: self._process_recipe_result(result_data))

    def _on_recipe_error(self, error_data: tuple):
        """Slot called when the RecipeWorker fails."""
        title, message, trace = error_data
        
        # 1. Cleanup UI
        self._cleanup_progress_ui()
            
        # [NEW] Append warnings to the message
        if hasattr(self, '_captured_warnings') and self._captured_warnings:
            unique_warnings = list(dict.fromkeys(self._captured_warnings))
            if unique_warnings:
                message += "<br><br><b>Previous Warnings:</b><ul>"
                for w in unique_warnings:
                    message += f"<li>{w}</li>"
                message += "</ul>"

        if trace:
            # --- CASE A: System Crash (Show Traceback) ---
            # We format the traceback as a code block for readability
            # Replace newlines with <br> for HTML rendering if not using <pre>
            # But <pre> is safer for code.
            
            # Truncate trace if massive
            if len(trace) > 2000:
                display_trace = trace[:2000] + "\n... (traceback truncated)"
            else:
                display_trace = trace

            full_text = (
                f"{message}<br><br>"
                f"<b>Traceback:</b>"
                f"<pre style='font-size: 10pt; background-color: #f0f0f0; padding: 5px;'>{display_trace}</pre>"
            )

            self._show_custom_message(
                title="Internal Error",
                header=title, # e.g. "Recipe Error: fit_continuum"
                text=full_text,
                buttons=QMessageBox.StandardButton.Ok,
                icon_name="icon_3d_HR.png" # Or a warning icon if you have one
            )

        else:
            # --- CASE B: User Error (Offer Retry) ---
            # e.g. "Param X is invalid"
            
            ret = self._show_custom_message(
                title="Recipe Failed",
                header=title,   # e.g. "Recipe Aborted"
                text=message,   # e.g. "Value must be positive."
                buttons=QMessageBox.StandardButton.Retry | QMessageBox.StandardButton.Cancel,
                default_btn=QMessageBox.StandardButton.Retry
            )
            
            # If user clicked Retry (mapped from "Try Again"), re-launch
            if ret == QMessageBox.StandardButton.Retry:
                if self._last_attempted_recipe:
                    logging.info("User clicked Retry, re-launching dialog.")
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
        ret = self._show_custom_message(
            title="Run Script",
            header="Run script and replace history?", # Grassetto
            text="This action will clear your current session history and re-run the script from the original data.\n\nAre you sure?", # Normale
            buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            default_btn=QMessageBox.StandardButton.Cancel
        )
        if ret != QMessageBox.StandardButton.Yes:
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

    def launch_split_from_current_view(self):
        """
        Launches the 'split' recipe using the current plot view limits (zoom).
        """
        if not self.plot_viewer or not self.active_history:
            return

        # 1. Get current X-limits from the plot (always in nm in Astrocook V2?)
        # We check _update_limit_boxes_from_plot to confirm logic: 
        # "xlim_nm = self.plot_viewer.canvas.axes.get_xlim()"
        try:
            xlim_nm = self.plot_viewer.canvas.axes.get_xlim()
            xmin, xmax = xlim_nm[0], xlim_nm[1]
        except Exception as e:
            logging.error(f"Could not get plot limits: {e}")
            return

        # 2. Build the expression string
        # Note: 'x' in the data is typically nm. 
        expression = f"(x > {xmin:.4f}) & (x < {xmax:.4f})"
        
        logging.info(f"Splitting from view: {expression}")

        # 3. Launch the dialog with these params pre-filled
        self._launch_recipe_dialog(
            "edit", 
            "split", 
            initial_params={"expression": expression}
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
        if hasattr(self, 'view_script_console_action'): self.view_script_console_action.setEnabled(is_valid_session)
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
        if hasattr(self, 'import_systems_action'): self.import_systems_action.setEnabled(enable_recipes)
        if hasattr(self, 'delete_action'): self.delete_action.setEnabled(enable_recipes)
        if hasattr(self, 'split_action'): self.split_action.setEnabled(enable_recipes)
        if hasattr(self, 'extract_preset_action'): self.extract_preset_action.setEnabled(enable_recipes)
        if hasattr(self, 'coadd_action'): self.coadd_action.setEnabled(enable_recipes)
        
        if hasattr(self, 'calculate_running_std_action'): self.calculate_running_std_action.setEnabled(enable_recipes)
        if hasattr(self, 'smooth_action'): self.smooth_action.setEnabled(enable_recipes)
        if hasattr(self, 'rebin_action'): self.rebin_action.setEnabled(enable_recipes)
        if hasattr(self, 'resample_action'): self.resample_action.setEnabled(enable_recipes)
        if hasattr(self, 'calibrate_action'): self.calibrate_action.setEnabled(enable_recipes)

        if hasattr(self, 'auto_cont_action'): self.auto_cont_action.setEnabled(enable_recipes)
        if hasattr(self, 'find_abs_action'): self.find_abs_action.setEnabled(enable_recipes)
        if hasattr(self, 'fit_cont_action'): self.fit_cont_action.setEnabled(enable_recipes)
        if hasattr(self, 'fit_powerlaw_action'): self.fit_powerlaw_action.setEnabled(enable_recipes)

        if hasattr(self, 'doublets_auto_action'): self.doublets_auto_action.setEnabled(enable_recipes)
        if hasattr(self, 'lya_auto_action'): self.lya_auto_action.setEnabled(enable_recipes)
        if hasattr(self, 'forest_auto_action'): self.forest_auto_action.setEnabled(enable_recipes)
        if hasattr(self, 'metals_auto_action'): self.metals_auto_action.setEnabled(enable_recipes)
        if hasattr(self, 'identify_lines_action'): self.identify_lines_action.setEnabled(enable_recipes)
        if hasattr(self, 'refit_all_action'): self.refit_all_action.setEnabled(enable_recipes)
        if hasattr(self, 'clean_negligible_action'): self.clean_negligible_action.setEnabled(enable_recipes)
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
        Aggressively close and destroy the progress dialog.
        Uses a QTimer to ensure this happens on the next event loop cycle,
        clearing any pending update signals first.
        """
        if not self.progress_dialog:
            return

        dialog_ref = self.progress_dialog  # Capture local reference
        self.progress_dialog = None        # Clear instance ref immediately
        
        def _destroy():
            try:
                logging.debug("Force-closing progress dialog...")
                
                # 1. Break Log Connections (Stop incoming text)
                if hasattr(self, 'log_bridge') and self.log_bridge:
                    try:
                        self.log_bridge.new_log_message.disconnect(dialog_ref.update_log)
                    except:
                        pass 

                # 2. Hard Close
                dialog_ref.hide()               # Visual removal
                dialog_ref.done(QDialog.Rejected) # Logical rejection
                dialog_ref.close()              # Event close
                dialog_ref.deleteLater()        # Memory cleanup
                
            except Exception as e:
                logging.error(f"Error destroying dialog: {e}")

        # Execute on the next frame of the GUI loop (0ms delay)
        # This allows pending log signals to flush safely before we kill the receiver.
        QTimer.singleShot(0, _destroy)

    def _on_debug_plot_requested(self, plot_data: dict):
        """
        SLOT: Receives data from a worker thread and plots it in
        a new Matplotlib window. Runs on the main GUI thread.
        """
        try:
            import matplotlib.pyplot as plt
            
            logging.info(f"Received debug plot request: {plot_data.get('title', 'Debug Plot')}")
            
            # --- OPTIMIZATION PLOT TYPE ---
            if plot_data.get('type') == 'optimization_step':
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
                
                # Extract Data
                x = plot_data['x']
                mask = plot_data['mask']
                
                # Zoom view to the relevant window
                x_view = x[mask]
                if len(x_view) > 0:
                    pad = (x_view[-1] - x_view[0]) * 0.2
                    ax1.set_xlim(x_view[0] - pad, x_view[-1] + pad)
                
                # Plot Flux
                ax1.step(x, plot_data['y'], color='black', label='Data', where='mid')
                ax1.plot(x, plot_data['model'], color='red', label='Model')
                ax1.axvline(plot_data['candidate_x'], color='green', linestyle='--', label='New Candidate')
                
                # Plot Residuals
                ax2.step(x, plot_data['resid'], color='blue', where='mid')
                ax2.axhline(0, color='gray', linestyle=':')
                ax2.axhline(-3, color='orange', linestyle='--') # Typical threshold
                
                ax1.set_title(plot_data['title'])
                ax2.set_xlabel("Wavelength")
                ax2.set_ylabel("Resid (Sigma)")
                ax1.legend()
                
                plt.show()
                return
        
            else:
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
        
        session_name = os.path.splitext(os.path.basename(file_path))[0]
        gui_context = self.mock_gui_context

        try:
            new_session = load_session_from_file(
                archive_path=file_path, 
                name=session_name, 
                format_name='auto',  # <--- Let the loader decide
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
    
    def open_inspector_on_component(self, uuid: str):
        """
        Opens the System Inspector and focuses on the specific component UUID.
        """
        # 1. Ensure Inspector is created and visible
        self._on_view_system_inspector()
        
        # 2. Delegate focus logic to the inspector
        if self.system_inspector:
            self.system_inspector.focus_on_component(uuid, force_group_view=True)

    def update_session_highlights(self, components: list):
        """
        Updates the main plot to highlight the specified components.
        """
        if not self.plot_viewer:
            return
            
        # Extract UUIDs for efficient storage
        uuids = {c.uuid for c in components} if components else set()
        self.plot_viewer.set_highlights(uuids)

    def dragEnterEvent(self, event):
        """Check if the drag contains files."""
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
        else:
            event.ignore()

    def dropEvent(self, event):
        """Process the dropped files."""
        urls = event.mimeData().urls()
        if not urls:
            return

        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            for url in urls:
                file_path = url.toLocalFile()
                if file_path:
                    logging.info(f"Dropped file detected: {file_path}")
                    # Use your existing public method
                    self.open_session_from_path(file_path)
        finally:
            QApplication.restoreOverrideCursor()
        
        event.acceptProposedAction()

class LogSignalBridge(QObject):
    """Thread-safe bridge to emit log messages as Qt signals."""
    new_log_message = Signal(str)

class QtLogHandler(logging.Handler):
    """Python Logging Handler that pipes messages to a Qt Signal."""
    def __init__(self, bridge):
        super().__init__()
        self.bridge = bridge
    
    def emit(self, record):
        # Format the message
        msg = self.format(record)
        # Emit signal (thread-safe way to talk to GUI)
        self.bridge.new_log_message.emit(msg)

class RecipeProgressDialog(QDialog):
    """
    A minimal progress dialog for long-running recipes.
    Shows a determinate bar and the last status message.
    """
    stop_requested = Signal()

    def __init__(self, title, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Running {title}")
        self.setModal(True)
        # Resize to be compact but wide enough for text
        self.resize(400, 120)
        
        # Remove [X] to prevent unsafe closing
        self.setWindowFlags(self.windowFlags() & ~Qt.WindowCloseButtonHint)

        layout = QVBoxLayout(self)
        layout.setSpacing(15)
        layout.setContentsMargins(20, 20, 20, 20)
        
        # Header
        self.header_label = QLabel(f"Recipe Running...")
        self.header_label.setStyleSheet("font-size: 14px;")
        layout.addWidget(self.header_label)

        # Progress Bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        # Optional: Make it look a bit thicker/modern
        self.progress_bar.setStyleSheet("QProgressBar { min-height: 20px; }")
        layout.addWidget(self.progress_bar)

        # Status Label (Last log line)
        self.status_label = QLabel("Initializing...")
        self.status_label.setStyleSheet("font-size: 11px;")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        # Add a stretch to keep things tight at the top
        layout.addStretch()

        # Use a horizontal layout to manage alignment
        btn_layout = QHBoxLayout()
        btn_layout.setContentsMargins(0, 0, 0, 0)

        # 1. Add a "spring" to push the button to the right
        btn_layout.addStretch()

        # [NEW] Stop Button
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setMaximumWidth(100)
        self.stop_btn.setToolTip("Stop execution gracefully")
        self.stop_btn.clicked.connect(self._on_stop_clicked)
        
        # 3. Add button to the horizontal row
        btn_layout.addWidget(self.stop_btn)

        # 4. Add the row to the main vertical layout
        layout.addLayout(btn_layout)

    def update_log(self, message):
        """
        Parses logs. 
        - '>> PROGRESS: 50' -> Updates bar.
        - Other text -> Updates status label.
        """
        if ">> PROGRESS:" in message:
            try:
                parts = message.split(">> PROGRESS:")
                val = int(parts[1].strip())
                self.progress_bar.setValue(val)
            except ValueError:
                pass
            return

        clean_msg = message.strip()
        if clean_msg:
            # Update the text label
            self.status_label.setText(clean_msg)
            # Force UI update immediately so text doesn't lag behind heavy calculations
            QApplication.processEvents()
    
    def _on_stop_clicked(self):
        self.stop_btn.setEnabled(False)
        self.stop_btn.setText("Stopping...")
        self.status_label.setText("Attempting to stop gracefully...")
        # Emit signal to parent
        self.stop_requested.emit()
class WarningCaptureHandler(logging.Handler):
    """Captures warning messages to a list for GUI display."""
    def __init__(self):
        super().__init__()
        self.warnings = []
        
    def emit(self, record):
        # Only capture meaningful warnings (skip debug/info)
        if record.levelno >= logging.WARNING:
            self.warnings.append(self.format(record))

class AboutDialog(QDialog):
    """
    Dialog to display application information, credits, and citation.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("About Astrocook")
        self.setFixedWidth(400)
        
        layout = QVBoxLayout(self)
        layout.setSpacing(5)
        layout.setContentsMargins(20, 20, 20, 20)

        layout.addSpacing(20)

        # --- Logo ---
        try:
            # Reusing the resource path logic from MainWindow
            icon_path = resource_path(os.path.join("assets", "logo_3d_HR.png"))
            if os.path.exists(icon_path):
                logo_lbl = QLabel()
                pixmap = QPixmap(icon_path).scaled(160, 160, Qt.KeepAspectRatio, Qt.SmoothTransformation)
                logo_lbl.setPixmap(pixmap)
                logo_lbl.setAlignment(Qt.AlignCenter)
                layout.addWidget(logo_lbl)
        except Exception: 
            pass

        # --- Title & Version ---
        title_lbl = QLabel("Astrocook V2")
        title_lbl.setStyleSheet("font-size: 20pt; font-weight: bold; color: #296bff;")
        title_lbl.setAlignment(Qt.AlignCenter)
        #layout.addWidget(title_lbl)

        copy_lbl = QLabel("© 2017-2026 Guido Cupani")
        copy_lbl.setStyleSheet("font-size: 14pt; color: #296bff;")
        copy_lbl.setAlignment(Qt.AlignCenter)
        layout.addWidget(copy_lbl)

        # Try to import version, else fallback
        try:
            from astrocook import __version__ as ver
        except ImportError:
            ver = "v2.0.0-beta.2"

        ver_lbl = QLabel(f"Version: {ver}")
        ver_lbl.setStyleSheet("font-weight: bold; color: #f5a100;")
        ver_lbl.setAlignment(Qt.AlignCenter)
        layout.addWidget(ver_lbl)

        layout.addSpacing(20)

        tagline_lbl = QLabel("A thousand recipes to cook a spectrum")
        tagline_lbl.setStyleSheet("font-size: 16pt; font-weight: bold;")
        tagline_lbl.setAlignment(Qt.AlignCenter)
        layout.addWidget(tagline_lbl)

        layout.addSpacing(20)

        # --- Description / Credits ---
        credits_html = (
            "<p align='center'>Developed by <b>Guido Cupani</b><br>at INAF–Osservatorio Astronomico di Trieste.</p>"
            "<p align='center'><b>Many thanks</b> to (in alphabetical order):<br>Giorgio Calderone, Stefano Cristiani, Simona Di Stefano,<br>Valentina D'Odorico, Francesco Guarneri, Elena Marcuzzo,<br> Stefano Alberto Russo, Andrea Trost.</p>"
            "<p align='center'>Logo and icon by <a href='https://claudiacupani.com'><b>Claudia Cupani</b></a>.</p>"
            "<p align='center'>Partly funded by the PRIN PNRR project<br>“Next generation computing & data technologies<br>to probe the cosmic metal content”.</p>"
        )
        credits_lbl = QLabel(credits_html)
        credits_lbl.setWordWrap(True)
        credits_lbl.setAlignment(Qt.AlignCenter)
        credits_lbl.setOpenExternalLinks(True)
        layout.addWidget(credits_lbl)

        layout.addSpacing(30)

        # --- Buttons ---
        # Create the standard Dialog Button Box
        btn_box = QDialogButtonBox(QDialogButtonBox.Ok)
        btn_box.accepted.connect(self.accept)

        # 1. Report Issue Button (GitHub)
        issue_btn = QPushButton("Report Issue")
        issue_btn.setToolTip("Open the GitHub Issues page to report bugs or request features")
        issue_btn.clicked.connect(self._open_github_issues)
        btn_box.addButton(issue_btn, QDialogButtonBox.ActionRole)
        
        # 2. Cite Button (ADS)
        cite_btn = QPushButton("Cite on ADS")
        cite_btn.setToolTip("Open the NASA ADS citation export page for Astrocook")
        cite_btn.clicked.connect(self._open_ads_citation)
        
        # Add "Cite" to the button box with ActionRole 
        # (This usually places it to the left of the OK button)
        btn_box.addButton(cite_btn, QDialogButtonBox.ActionRole)

        layout.addWidget(btn_box)

    def _open_ads_citation(self):
        """Opens the ADS export citation page in the default browser."""
        # URL for the 2020 Astronomy and Computing paper
        url = "https://ui.adsabs.harvard.edu/abs/2020SPIE11452E..1UC/exportcitation"
        QDesktopServices.openUrl(QUrl(url))

    def _open_github_issues(self):
        """Opens the GitHub Issues page in the default browser."""
        # Update this URL if the repo is under a different organization
        url = "https://github.com/DAS-OATs/astrocook/issues"
        QDesktopServices.openUrl(QUrl(url))