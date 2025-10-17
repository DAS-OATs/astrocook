# astrocook/v2/gui/main_window.py

import logging
import os
from PySide6.QtCore import Qt, QStringListModel
from PySide6.QtGui import QAction, QIcon
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QMenuBar, QStackedWidget, QLabel, QDockWidget, QListView, QFileDialog
)

from .pyside_plot import SpectrumPlotWidget
from ..session import load_session_from_file

class MainWindowV2(QMainWindow):
    def __init__(self, session): # Accept the loaded session object
        super().__init__()
        
        # 1. Initialize the internal list to empty.
        self.active_sessions = [] 
        
        # 2. Set the session manager (the first session object)
        self.session_manager = session 
        
        # 3. Create the list model for QListView (must be defined early)
        self.session_model = QStringListModel()

        self.setGeometry(100, 100, 1400, 900)
        self.central_stack = QStackedWidget()
        self.setCentralWidget(self.central_stack)
        
        # Setup plot view (Central widget) and empty view
        self._setup_plot_view() 
        self._setup_empty_view()
        
        # Setup Dockable Session Panel
        self._setup_session_panel()
        
        self._create_menubar()
        
        if self.session_manager.spec and len(self.session_manager.spec.x) > 0:
            # If a session was successfully loaded initially (e.g., from command line), add it.
            self.add_session(session, initial_load=True)
        else:
            # If the session is an empty placeholder, explicitly set the view to empty.
            self.central_stack.setCurrentIndex(1)
    
    def _setup_plot_view(self):
        """Sets up the primary spectrum plot area."""
        # NOTE: SpectrumPlotWidget is the refactored widget from your prototype.
        # It handles the plot, toolbar, and toggles for the *active* session.
        self.plot_viewer = SpectrumPlotWidget(self.session_manager)
        self.central_stack.addWidget(self.plot_viewer)

    def _setup_empty_view(self):
        """Sets up the view when no spectrum is loaded."""
        empty_widget = QWidget()
        layout = QVBoxLayout(empty_widget)
        
        # Use a QLabel to provide clear instructions
        label = QLabel("Welcome to Astrocook v2.\n\n"
                       "Please use the 'File' menu to load a spectrum (.fits, .acs, etc.).")
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setFont(label.font()) # Optional: make it larger
        
        layout.addWidget(label)
        self.central_stack.addWidget(empty_widget)

    def _setup_session_panel(self):
        self.session_dock = QDockWidget("Sessions", self)
        self.session_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        
        # NOTE: QListView is typically vertical. For a horizontal bar, we use a QToolBar or a QListView with a horizontal flow delegate.
        # For simplicity, we use QListView but dock it horizontally.
        self.session_list_view = QListView()
        
        # Use QModel to link Python list to QListView
        self.session_list_view.setModel(self.session_model)
        
        # Set the dock widget's content
        self.session_dock.setWidget(self.session_list_view)
        
        # Add the dock widget to the main window
        self.addDockWidget(Qt.LeftDockWidgetArea, self.session_dock)
        
        # Connect signal for session switching (when a list item is clicked)
        self.session_list_view.clicked.connect(self._on_session_switched)

    def _on_session_switched(self, index):
        """Switches the main window's active session when the list view is clicked."""
        session_index = index.row()
        new_active_session = self.active_sessions[session_index]
        self.update_session(new_active_session)


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

    def _on_open_spectrum(self):
        """Launches the file dialog and initiates V2 loading."""
        file_name, _ = QFileDialog.getOpenFileName(
            self, 
            "Open Spectrum File", 
            os.getcwd(), # Start in current working directory
            "Spectrum Files (*.fits *.acs *.dat *.txt *.json);;All Files (*)"
        )
        
        if file_name:
            # We assume a fixed format name for now, as V1 auto-detection is complex
            format_name = 'generic_spectrum' 
            
            try:
                # CRITICAL: Call the utility function directly
                new_session = load_session_from_file(
                    file_path=file_name, 
                    format_name=format_name,
                    gui_context=self.session_manager._gui # Assume GUI placeholder is accessible
                )
                
                # Add the new session to the application state
                self.add_session(new_session)
                
            except Exception as e:
                logging.error(f"Failed to load file via V2 adapter: {e}")

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
        
        self.active_sessions.append(new_session)
        
        # FIX: Ensure session name is set from the file path
        if hasattr(new_session, 'path') and new_session.path:
            new_session.name = os.path.basename(new_session.path)
        
        # CRITICAL FIX: Always update the model list after modification
        self.session_model.setStringList([s.name for s in self.active_sessions])
        
        # Set the new session as active and redraw
        self.update_session(new_session)


    def update_session(self, new_session):
        """Swaps the central session manager object and updates the view."""
        
        # FIX: The update logic must handle replacing the current session manager instance
        
        # 1. Update the session manager reference
        self.session_manager = new_session

        # 2. Find the index of the old session if a replacement occurred (not necessary here, as we only swap pointers)
        
        # 3. Update the central view
        if new_session.spec and len(new_session.spec.x) > 0:
            self.central_stack.setCurrentIndex(0) 
            self.plot_viewer.update_plot(new_session)
        else:
            self.central_stack.setCurrentIndex(1)

# NOTE: This structure requires a central Session Manager to be passed during initialization.