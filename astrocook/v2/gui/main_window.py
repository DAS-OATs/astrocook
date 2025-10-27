# astrocook/v2/gui/main_window.py

import logging
import os
from PySide6.QtCore import ( # <<< Modify this import
    Qt, QStringListModel, QSize,
    QItemSelectionModel # <<< Add QItemSelectionModel
)
from PySide6.QtGui import QAction, QIcon, QPalette, QFont # Import QPalette
from PySide6.QtWidgets import (
    QCheckBox, 
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QStackedWidget, QLabel, QDockWidget, QListView, QFileDialog, QApplication, 
    QPushButton, QSizePolicy, QSpacerItem, QStyle
)

from .pyside_plot import SpectrumPlotWidget
from ..session import load_session_from_file

class MainWindowV2(QMainWindow):
    def __init__(self, session):
        super().__init__()
        self.active_sessions = []
        self.session_manager = session
        self.session_model = QStringListModel()

        self.setGeometry(100, 100, 450, 150) # Initial small size
        screen_geometry = QApplication.primaryScreen().geometry()
        x = (screen_geometry.width() - self.width()) // 2
        y = (screen_geometry.height() - self.height()) // 2
        self.move(x, y)

        # --- Central Widget Setup ---
        # Main widget to hold the layout
        self.main_widget = QWidget()
        self.main_layout = QHBoxLayout(self.main_widget)
        self.main_layout.setContentsMargins(0, 0, 0, 0)
        self.main_layout.setSpacing(0) # No space between dock area and content

        # Area for dock widgets (will contain the session panel)
        # QMainWindow handles dock areas implicitly, but we manage the central part
        self.setCentralWidget(self.main_widget)

        # --- Content Stack (Plot or Empty View) ---
        self.central_stack = QStackedWidget()
        self._setup_plot_view()
        self._setup_empty_view()
        # Add stack AFTER dock setup if button is between dock and stack

        # --- Setup Dockable Session Panel ---
        self._setup_session_panel() # This ADDS the dock widget to the main window

        # --- ** Collapse Button (Now part of main layout) ** ---
        self.session_collapse_button = QPushButton() # No text initially
        self.session_collapse_button.setToolTip("Collapse/Expand Session Panel")
        self.session_collapse_button.setFixedSize(QSize(20, 40)) # Narrower, taller
        self.session_collapse_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding) # Fill vertically
        self.session_collapse_button.clicked.connect(self._toggle_session_dock)
        self.main_layout.addWidget(self.session_collapse_button)
        
        # ** Add button to the main layout, between dock and central stack **
        # (QMainWindow manages dock areas, so we add stack to main layout)
        self.main_layout.addWidget(self.central_stack, 1) # Stack takes remaining space

        # --- ** Right Collapse Button ** ---
        self.plot_controls_collapse_button = QPushButton() # NEW button
        self.plot_controls_collapse_button.setToolTip("Collapse/Expand Plot Controls")
        self.plot_controls_collapse_button.setFixedSize(QSize(20, 40))
        self.plot_controls_collapse_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        self.plot_controls_collapse_button.clicked.connect(self._toggle_plot_controls_dock) # New slot
        self.plot_controls_collapse_button.setObjectName("CollapseButton") # Keep common ID for style
        self.main_layout.addWidget(self.plot_controls_collapse_button) # Add RIGHT button

        # --- Setup Dockable Plot Controls Panel (Right) ---
        self._setup_plot_controls_panel() # Adds the right dock widget

        self._create_menubar()

        # --- Styling ---
        self._apply_styles() # Apply QSS

        # Initial state based on session
        if self.session_manager.spec and len(self.session_manager.spec.x) > 0:
            self.add_session(session, initial_load=True)
            self._update_session_collapse_button_icon(True) # Set initial icon
            self.session_collapse_button.setVisible(True) # ** Show button **
            self._update_plot_controls_collapse_button_icon(True) # Update right button
            self.plot_controls_collapse_button.setVisible(True)
            self.plot_controls_dock.setVisible(True) # Make sure right dock is visible
        else:
            self.central_stack.setCurrentIndex(1) # Show empty view
            self.session_dock.setVisible(False) # Hide dock initially
            self.plot_controls_dock.setVisible(False) # Hide dock initially
            self._update_session_collapse_button_icon(False) # Set initial icon
            self.session_collapse_button.setVisible(False) # ** Hide button initially **
            self._update_plot_controls_collapse_button_icon(False)
            self.plot_controls_collapse_button.setVisible(False)

        # Connect dock visibility change AFTER button exists
        self.session_dock.visibilityChanged.connect(self._update_session_collapse_button_icon)
        self.plot_controls_dock.visibilityChanged.connect(self._update_plot_controls_collapse_button_icon) # Connect right dock

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

    def _setup_session_panel(self):
        self.session_dock = QDockWidget("Sessions", self)
        self.session_dock.setAllowedAreas(Qt.LeftDockWidgetArea)
        self.session_dock.setFeatures(QDockWidget.DockWidgetMovable) # Only movable
        self.session_dock.setTitleBarWidget(QWidget()) # Hide title bar

        container_widget = QWidget()
        container_widget.setFixedWidth(200)
        container_layout = QVBoxLayout(container_widget)
        container_layout.setContentsMargins(0,0,0,0)
        container_layout.setSpacing(0)

        self.session_list_view = QListView()
        self.session_list_view.setModel(self.session_model)
        font = self.session_list_view.font()
        font.setPointSize(14) # Set desired point size
        self.session_list_view.setFont(font)

        container_layout.addWidget(self.session_list_view)
        self.session_dock.setWidget(container_widget)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.session_dock)
        self.session_list_view.clicked.connect(self._on_session_switched)

        # Give IDs for styling
        container_widget.setObjectName("SessionContainer")
        self.session_list_view.setObjectName("SessionListView")

    def _setup_plot_controls_panel(self):
        """Creates the right-hand dock widget for plot toggles."""
        self.plot_controls_dock = QDockWidget("Plot Controls", self)
        self.plot_controls_dock.setAllowedAreas(Qt.RightDockWidgetArea)
        # **Hide title bar and limit features for consistency**
        self.plot_controls_dock.setFeatures(QDockWidget.DockWidgetMovable)
        self.plot_controls_dock.setTitleBarWidget(QWidget())

        controls_container = QWidget()
        controls_container.setFixedWidth(200)
        controls_layout = QVBoxLayout(controls_container)
        controls_layout.setContentsMargins(10, 10, 10, 10)
        controls_layout.setSpacing(8)

        # --- Checkboxes (unchanged) ---
        self.error_checkbox = QCheckBox("Show 1σ Error"); self.error_checkbox.setChecked(True)
        self.error_checkbox.stateChanged.connect(self._trigger_replot)
        controls_layout.addWidget(self.error_checkbox)

        self.continuum_checkbox = QCheckBox("Show Continuum"); self.continuum_checkbox.setChecked(False)
        self.continuum_checkbox.stateChanged.connect(self._trigger_replot)
        controls_layout.addWidget(self.continuum_checkbox)
        # ... (add more later) ...

        spacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        controls_layout.addItem(spacer)
        self.plot_controls_dock.setWidget(controls_container)
        self.addDockWidget(Qt.RightDockWidgetArea, self.plot_controls_dock)

        controls_container.setObjectName("PlotControlsContainer") # For styling
        self.error_checkbox.setObjectName("PlotControlCheckbox")
        self.continuum_checkbox.setObjectName("PlotControlCheckbox")

    def _trigger_replot(self):
        """Slot to be called when any plot toggle checkbox changes state."""
        if self.plot_viewer: # Check if plot viewer exists
            # Call the plot function, which will now read checkbox states
            self.plot_viewer.plot_spectrum()

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
        open_action = QAction("&Open Spectrum...", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self._on_open_spectrum)
        file_menu.addAction(open_action)

        # The QDockWidget automatically provides an action to toggle its visibility.
        view_menu.addAction(self.session_dock.toggleViewAction())

        # Action to toggle the NEW plot controls dock
        view_menu.addAction(self.plot_controls_dock.toggleViewAction())

        # --- Add QActions for V2 Recipes ---
        
        # RECIPES FOR 'EDIT' MENU (x_convert, y_convert)
        
        # x_convert Action
        x_convert_action = QAction("Convert &X Axis...", self)
        x_convert_action.triggered.connect(self._launch_x_convert_dialog)
        edit_menu.addAction(x_convert_action)
        
        # y_convert Action
        y_convert_action = QAction("Convert &Y Axis...", self)
        y_convert_action.triggered.connect(self._launch_y_convert_dialog)
        edit_menu.addAction(y_convert_action)
        
        # RECIPES FOR 'FLUX' MENU (rebin)
        
        # rebin Action
        rebin_action = QAction("&Rebin Spectrum...", self)
        rebin_action.triggered.connect(self._launch_rebin_dialog)
        flux_menu.addAction(rebin_action)
        
        # ... (Other menu item actions will be added here later) ...

    def _apply_styles(self):
        """Applies QSS styling to the relevant widgets."""

        # Try to get a theme-aware background color (might just be standard window bg)
        try:
            palette = QApplication.palette()
            # Use Button role for a slightly darker/more distinct sidebar usually
            sidebar_bg = palette.color(palette.ColorRole.Button).name()
            item_selected_bg = palette.color(palette.ColorRole.Highlight).name()
            item_selected_text = palette.color(palette.ColorRole.HighlightedText).name()
            item_hover_bg = palette.color(palette.ColorRole.Midlight).name() # Or Mid
            #window_bg = palette.color(palette.ColorRole.Window).name() # Get standard window bg
            text_color = palette.color(palette.ColorRole.WindowText).name() # Standard text
            button_bg = palette.color(palette.ColorRole.Button).name() # Use Button bg
            button_fg = palette.color(palette.ColorRole.ButtonText).name() # Use Button text for icon
            button_hover_bg = palette.color(palette.ColorRole.Midlight).name()
            border_color = palette.color(palette.ColorRole.Mid).name()
        except:
            # Fallback colors if palette query fails
            sidebar_bg = "#E0E0E0" # Light gray fallback
            item_selected_bg = "#808080" # Mid-gray selected
            item_selected_text = "#FFFFFF" # White text
            item_hover_bg = "#D0D0D0"
            border_color = "#B0B0B0"


        qss = f"""
            /* Left Sidebar Container */
            QWidget#SessionContainer {{
                background-color: {sidebar_bg};
                border: none; /* No internal border */
            }}
            /* Left Sidebar List View */
            QListView#SessionListView {{
                background-color: transparent;
                border: none;
            }}
            QListView#SessionListView::item {{
                padding: 6px 10px;
                border: none;
                background-color: transparent;
                /* Font set directly on widget */
            }}
            QListView#SessionListView::item:selected {{ background-color: {item_selected_bg}; color: {item_selected_text}; }}
            QListView#SessionListView::item:hover {{ background-color: {item_hover_bg}; }}

            /* Right Sidebar Container */
            QWidget#PlotControlsContainer {{
                background-color: {sidebar_bg}; /* ** Match left sidebar bg ** */
                border: none; /* No internal border */
            }}
            /* Right Sidebar Checkboxes */
            QCheckBox#PlotControlCheckbox {{
                color: {text_color};
                spacing: 5px;
                padding: 4px 0px; /* Add some vertical padding */
            }}

            /* General Style for BOTH Collapse Buttons */
            QPushButton#CollapseButton {{
                background-color: {button_bg}; /* ** Use explicit button background ** */
                color: {button_fg}; /* ** Explicitly set icon/text color ** */
                border: 1px solid {border_color}; /* Add a border for definition */
                border-radius: 3px; /* Slightly rounded */
                margin: 0px;
                padding: 0px;
                /* Ensure icon size fits */
                icon-size: 16px; /* Adjust if needed */
            }}
            QPushButton#CollapseButton:hover {{
                 background-color: {button_hover_bg};
            }}

            /* Style Dock Separators (Optional) */
            QMainWindow::separator {{
                background: {palette.color(palette.ColorRole.Mid).name() if 'palette' in locals() else '#C0C0C0'};
                width: 1px; /* Vertical */
                height: 1px; /* Horizontal */
            }}
            QMainWindow::separator:hover {{
                background: {palette.color(palette.ColorRole.Highlight).name() if 'palette' in locals() else '#808080'};
            }}
        """
        self.setStyleSheet(qss) # Apply to main window, cascade down
        self.session_collapse_button.setObjectName("CollapseButton") # Set ID for QSS
        self.plot_controls_collapse_button.setObjectName("CollapseButton")

    def _toggle_session_dock(self):
        is_visible = self.session_dock.isVisible()
        self.session_dock.setVisible(not is_visible)
        # Icon update handled by visibilityChanged signal

    def _toggle_plot_controls_dock(self): # <<< NEW METHOD
        """Toggles the visibility of the plot controls dock widget."""
        is_visible = self.plot_controls_dock.isVisible()
        self.plot_controls_dock.setVisible(not is_visible)

    def _update_session_collapse_button_icon(self, visible):
        """Updates the collapse button icon based on dock visibility."""
        if visible:
            # Use standard Pixmap for left arrow
            icon = self.style().standardIcon(QStyle.SP_ArrowLeft)
            self.session_collapse_button.setToolTip("Collapse Session Panel")
        else:
            # Use standard Pixmap for right arrow
            icon = self.style().standardIcon(QStyle.SP_ArrowRight)
            self.session_collapse_button.setToolTip("Expand Session Panel")
        self.session_collapse_button.setIcon(icon)

    def _update_plot_controls_collapse_button_icon(self, visible): # <<< NEW METHOD
        """Updates the right collapse button icon."""
        if visible:
            # Right sidebar collapses TO the right, so arrow points right
            icon = self.style().standardIcon(QStyle.SP_ArrowRight)
            self.plot_controls_collapse_button.setToolTip("Collapse Plot Controls")
        else:
            # Right sidebar expands FROM the right, so arrow points left
            icon = self.style().standardIcon(QStyle.SP_ArrowLeft)
            self.plot_controls_collapse_button.setToolTip("Expand Plot Controls")
        self.plot_controls_collapse_button.setIcon(icon)

    def update_session(self, new_session, set_current_list_item=False):
        """Swaps the central session manager object and updates the view."""
        self.session_manager = new_session
        is_valid_session = new_session and new_session.spec and len(new_session.spec.x) > 0

        if is_valid_session:
            if self.central_stack.currentIndex() != 0: # If switching from empty
                self.resize(1400, 900)
                # Recenter after resize
                screen_geometry = QApplication.primaryScreen().geometry()
                x = (screen_geometry.width() - self.width()) // 2
                y = (screen_geometry.height() - self.height()) // 2
                self.move(x, y)
                # Ensure dock is visible when a valid session is active
                if not self.session_dock.isVisible():
                    self.session_dock.setVisible(True)

            if not self.session_collapse_button.isVisible():
                 self.session_collapse_button.setVisible(True)
            if not self.plot_controls_collapse_button.isVisible():
                 self.plot_controls_collapse_button.setVisible(True)
            if not self.session_dock.isVisible():
                self.session_dock.setVisible(True)
            if not self.plot_controls_dock.isVisible():
                 self.plot_controls_dock.setVisible(True)

            self.central_stack.setCurrentIndex(0) # Switch to plot view
            self.plot_viewer.update_plot(new_session)

            if set_current_list_item:
                 try:
                    idx = self.active_sessions.index(new_session)
                    q_model_index = self.session_model.index(idx)
                    # Use selection model for better control
                    selection_model = self.session_list_view.selectionModel()
                    selection_flag = QItemSelectionModel.SelectionFlag.ClearAndSelect
                    selection_model.setCurrentIndex(q_model_index, selection_flag)
                 except (ValueError, IndexError):
                     logging.warning("Could not find new session in list view to select.")

        else: # Handle empty/invalid session
            logging.info("Updating to an empty session view.")
            self.central_stack.setCurrentIndex(1)

            # ** Hide Plot Controls dock for empty view **
            if self.plot_controls_dock.isVisible():
                 self.plot_controls_dock.setVisible(False)

            # Hide BOTH collapse buttons and docks for empty view
            if self.session_collapse_button.isVisible():
                 self.session_collapse_button.setVisible(False)
            if self.plot_controls_collapse_button.isVisible():
                 self.plot_controls_collapse_button.setVisible(False)
            if self.session_dock.isVisible():
                 self.session_dock.setVisible(False)
            if self.plot_controls_dock.isVisible():
                 self.plot_controls_dock.setVisible(False)

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

    def _launch_x_convert_dialog(self):
        # Placeholder function called by the QAction
        print("Launching X Convert Dialog...")
        
    def _launch_y_convert_dialog(self):
        # Placeholder function called by the QAction
        print("Launching Y Convert Dialog...")
        
    def _launch_rebin_dialog(self):
        # Placeholder function called by the QAction
        print("Launching Rebin Dialog...")
        

    def add_session(self, new_session, initial_load=False):
        """Adds a new session object and sets it as active."""

        # Make sure the session object is valid
        if new_session is None or isinstance(new_session, int):
             logging.error("add_session received invalid session object.")
             return

        self.active_sessions.append(new_session)

        # Set session name if not already set during loading
        if not hasattr(new_session, 'name') or not new_session.name:
             # Fallback name based on file path if needed
             new_session.name = os.path.splitext(os.path.basename(getattr(new_session,'path', f"Session_{len(self.active_sessions)}")))[0]


        # Update the model list after modification
        self.session_model.setStringList([s.name for s in self.active_sessions])

        # Set the new session as active and redraw
        self.update_session(new_session, set_current_list_item=True)

    def _on_session_switched(self, index):
        """Switches the main window's active session when the list view is clicked."""
        try:
            session_index = index.row()
            if 0 <= session_index < len(self.active_sessions):
                new_active_session = self.active_sessions[session_index]
                self.update_session(new_active_session)
            else:
                 logging.warning(f"Invalid index {session_index} clicked in session list.")
        except Exception as e:
            logging.error(f"Error switching session: {e}")
