# astrocook/v2/gui/main_window.py

import logging
import os
from PySide6.QtCore import ( # <<< Modify this import
    Qt, QStringListModel, QSize,
    QItemSelectionModel, # <<< Add QItemSelectionModel
    QPropertyAnimation, QEasingCurve, QRect, QPoint,
    QParallelAnimationGroup
)
from PySide6.QtGui import QAction, QIcon, QPalette, QFont # Import QPalette
from PySide6.QtWidgets import (
    QCheckBox, 
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QStackedWidget, QLabel, QDockWidget, QListView, QFileDialog, QApplication, 
    QPushButton, QSizePolicy, QSpacerItem, QStyle
)

from .pyside_plot import SpectrumPlotWidget
from ..session import load_session_from_file

# --- Constants for Sidebar Widths ---
LEFT_SIDEBAR_WIDTH = 300
RIGHT_SIDEBAR_WIDTH = 300
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
        self._update_ui_state(is_initial_session_valid, is_startup=True)

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
        """Creates the right sidebar widget as a child of the main window."""
        self.right_sidebar_widget = QWidget(self) # ** Parent is main window **
        sidebar_layout = QVBoxLayout(self.right_sidebar_widget)
        sidebar_layout.setContentsMargins(10, 10, 10, 10); sidebar_layout.setSpacing(8)

        # --- Checkboxes ---
        self.error_checkbox = QCheckBox("Show 1σ Error"); # ... setup ...
        self.continuum_checkbox = QCheckBox("Show Continuum"); # ... setup ...
        self.error_checkbox.setChecked(True); self.error_checkbox.stateChanged.connect(self._trigger_replot)
        self.continuum_checkbox.setChecked(False); self.continuum_checkbox.stateChanged.connect(self._trigger_replot)
        sidebar_layout.addWidget(self.error_checkbox)
        sidebar_layout.addWidget(self.continuum_checkbox)
        # ... (spacer) ...
        spacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        sidebar_layout.addItem(spacer)

        self.right_sidebar_widget.setObjectName("PlotControlsContainer")
        self.error_checkbox.setObjectName("PlotControlCheckbox")
        self.continuum_checkbox.setObjectName("PlotControlCheckbox")

        # Initial geometry set later, start hidden/collapsed
        self.right_sidebar_widget.resize(0, self.height()) # Set initial size for geometry
        self.right_sidebar_widget.move(self.width(), self.right_sidebar_widget.y()) # Also set initial X pos offscreen
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

        # --- View Menu Actions ---
        toggle_left_action = QAction("Toggle Session Panel", self)
        toggle_left_action.triggered.connect(self._toggle_left_sidebar)
        view_menu.addAction(toggle_left_action)

        toggle_right_action = QAction("Toggle Plot Controls", self)
        toggle_right_action.triggered.connect(self._toggle_right_sidebar)
        view_menu.addAction(toggle_right_action)

        self.toggle_left_action = toggle_left_action
        self.toggle_right_action = toggle_right_action

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
        # ... (Get palette colors) ...
        # ... (Define sidebar_bg, item_selected_bg, etc.) ...
        try: # Added try-except for palette
            palette = QApplication.palette()
            # Get RGB components from the desired background color
            sidebar_base_color = palette.color(palette.ColorRole.Button)
            sidebar_r = sidebar_base_color.red()
            sidebar_g = sidebar_base_color.green()
            sidebar_b = sidebar_base_color.blue()
            sidebar_a = 0.5 # Alpha value (0.0 to 1.0), e.g., 90% opaque

            sidebar_bg = palette.color(palette.ColorRole.Button).name()
            item_selected_bg = palette.color(palette.ColorRole.Highlight).name()
            item_selected_text = palette.color(palette.ColorRole.HighlightedText).name()
            #item_hover_bg = palette.color(palette.ColorRole.Highlight).name()
            border_color = palette.color(palette.ColorRole.Mid).name()
            button_fg = palette.color(palette.ColorRole.ButtonText).name()
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
            
            /* Ensure Checkboxes are NOT transparent */
            QCheckBox#PlotControlCheckbox {{
                color: {button_fg}; spacing: 5px; padding: 4px 0px;
                background-color: transparent; /* Let container show through */
                padding-left: 15px;
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


    def update_session(self, new_session, set_current_list_item=False):
        """Swaps the central session manager object and updates the view."""
        self.session_manager = new_session
        is_valid = bool(new_session and new_session.spec and len(new_session.spec.x) > 0)

        # This call now handles all UI state changes correctly
        self._update_ui_state(is_valid, is_startup=(len(self.active_sessions) <= 1) ) # Pass startup hint

        if is_valid:
            self.plot_viewer.update_plot(new_session) # Update plot content
            if set_current_list_item:
                # ... (list item selection logic) ...
                try:
                    idx = self.active_sessions.index(new_session)
                    q_model_index = self.session_model.index(idx)
                    selection_model = self.session_list_view.selectionModel()
                    selection_flag = QItemSelectionModel.SelectionFlag.ClearAndSelect
                    selection_model.setCurrentIndex(q_model_index, selection_flag)
                except (ValueError, IndexError, AttributeError) as e:
                    logging.warning(f"Could not select session in list view: {e}")
        else:
            logging.info("Updating to an empty session view.")

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
