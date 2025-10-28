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
LEFT_SIDEBAR_WIDTH = 250
RIGHT_SIDEBAR_WIDTH = 200
ANIMATION_DURATION = 150 # ** Speed up animation **
BUTTON_WIDTH = 20
class MainWindowV2(QMainWindow):
    def __init__(self, session):
        super().__init__()
        self.active_sessions = []
        self.session_manager = session
        self.session_model = QStringListModel()

        # Animation objects (initialize early)
        self.left_sidebar_animation = None
        self.right_sidebar_animation = None

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
        content_y = menubar_height
        content_height = window_height - content_y

        # --- Left Sidebar & Button ---
        left_visible = self.left_sidebar_widget.isVisible()
        left_width = LEFT_SIDEBAR_WIDTH if left_visible else 0
        # **Button at absolute left edge**
        self.session_collapse_button.setGeometry(0, content_y, BUTTON_WIDTH, content_height)
        # **Sidebar positioned next to the button**
        self.left_sidebar_widget.setGeometry(BUTTON_WIDTH, content_y, left_width, content_height)

        # --- Right Sidebar & Button ---
        right_visible = self.right_sidebar_widget.isVisible()
        right_width = RIGHT_SIDEBAR_WIDTH if right_visible else 0
        # **Button at absolute right edge**
        self.plot_controls_collapse_button.setGeometry(window_width - BUTTON_WIDTH, content_y, BUTTON_WIDTH, content_height)
        # **Sidebar positioned next (left) to the button**
        self.right_sidebar_widget.setGeometry(window_width - BUTTON_WIDTH - right_width, content_y, right_width, content_height)

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
        font = self.session_list_view.font(); font.setPointSize(10); self.session_list_view.setFont(font)
        self.session_list_view.clicked.connect(self._on_session_switched)

        sidebar_layout.addWidget(self.session_list_view)
        self.left_sidebar_widget.setObjectName("SessionContainer")
        self.session_list_view.setObjectName("SessionListView")

        # Initial geometry set later, start hidden/collapsed
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
        self.right_sidebar_widget.setVisible(False)

    def _setup_collapse_buttons(self):
        """Creates and positions the collapse buttons as children of main window."""
        # Left Button
        self.session_collapse_button = QPushButton(self) # ** Parent is main window **
        self.session_collapse_button.setToolTip("Collapse/Expand Session Panel")
        self.session_collapse_button.setFixedSize(QSize(BUTTON_WIDTH, 40)) # Adjust height?
        self.session_collapse_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding) # Fill height
        self.session_collapse_button.clicked.connect(self._toggle_left_sidebar)
        self.session_collapse_button.setObjectName("CollapseButton")
        self.session_collapse_button.setVisible(False) # Start hidden

        # Right Button
        self.plot_controls_collapse_button = QPushButton(self) # ** Parent is main window **
        self.plot_controls_collapse_button.setToolTip("Collapse/Expand Plot Controls")
        self.plot_controls_collapse_button.setFixedSize(QSize(BUTTON_WIDTH, 40)) # Adjust height?
        self.plot_controls_collapse_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding) # Fill height
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
            sidebar_a = 0.3 # Alpha value (0.0 to 1.0), e.g., 90% opaque

            sidebar_bg = palette.color(palette.ColorRole.Button).name()
            item_selected_bg = palette.color(palette.ColorRole.Highlight).name()
            item_selected_text = palette.color(palette.ColorRole.HighlightedText).name()
            item_hover_bg = palette.color(palette.ColorRole.Midlight).name()
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
            /* Add side borders for visual separation */
            QWidget#SessionContainer {{ border-right: 1px solid {border_color}; }}
            QWidget#PlotControlsContainer {{ border-left: 1px solid {border_color}; }}

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
            QListView#SessionListView::item:hover {{ background-color: {item_hover_bg}; }}

            /* Ensure Checkboxes are NOT transparent */
            QCheckBox#PlotControlCheckbox {{
                color: {button_fg}; spacing: 5px; padding: 4px 0px;
                background-color: transparent; /* Let container show through */
            }}

            /* Collapse Buttons: Blend, NO borders */
            QPushButton#CollapseButton {{
                background-color: rgba({sidebar_r}, {sidebar_g}, {sidebar_b}, {sidebar_a});
                color: {button_fg}; /* Icon color */
                border: none; /* ** Remove all borders ** */
                margin: 0px; padding: 0px;
                icon-size: 16px;
                /* Maybe add slight rounding */
                border-radius: 3px;
            }}
            QPushButton#CollapseButton:hover {{
                 background-color: {item_hover_bg};
            }}
        """
        self.setStyleSheet(qss)
        # Set object names after creation if not done elsewhere

    # --- ** NEW Toggle & Animation Methods ** ---

    def _toggle_left_sidebar(self):
        """Animates the left sidebar overlay."""
        is_currently_visible = self.left_sidebar_widget.isVisible()
        start_geometry = QRect(self.left_sidebar_widget.geometry())
        end_geometry = QRect(start_geometry)

        if is_currently_visible: # Closing
            end_geometry.setWidth(0) # Animate width to 0
            target_visible = False
            # Button stays at x=0
        else: # Opening
            start_geometry.setWidth(0) # Start from width 0
            end_geometry.setWidth(LEFT_SIDEBAR_WIDTH)
            target_visible = True
            self.left_sidebar_widget.setVisible(True) # Show before animation

        # Animate sidebar width
        self.left_sidebar_animation = QPropertyAnimation(self.left_sidebar_widget, b"geometry") # Animate geometry
        self.left_sidebar_animation.setDuration(ANIMATION_DURATION)
        self.left_sidebar_animation.setStartValue(start_geometry)
        self.left_sidebar_animation.setEndValue(end_geometry)
        self.left_sidebar_animation.setEasingCurve(QEasingCurve.InOutQuad)
        if not target_visible: # Hide after closing animation
            self.left_sidebar_animation.finished.connect(lambda: self.left_sidebar_widget.setVisible(False))

        self.left_sidebar_animation.start()
        self._update_sidebar_button_icon(self.session_collapse_button, target_visible, is_left=True)
        # Button position doesn't change, no need to animate it

    def _toggle_right_sidebar(self):
        """Animates the right sidebar overlay."""
        is_currently_visible = self.right_sidebar_widget.isVisible()
        start_geometry = QRect(self.right_sidebar_widget.geometry())
        end_geometry = QRect(start_geometry)
        window_width = self.width()

        if is_currently_visible: # Closing
            end_geometry.setX(window_width - BUTTON_WIDTH) # Move X to button
            end_geometry.setWidth(0) # Animate width to 0
            target_visible = False
        else: # Opening
            start_geometry.setX(window_width - BUTTON_WIDTH) # Start at button
            start_geometry.setWidth(0)
            end_geometry.setX(window_width - BUTTON_WIDTH - RIGHT_SIDEBAR_WIDTH) # Move X left
            end_geometry.setWidth(RIGHT_SIDEBAR_WIDTH)
            target_visible = True
            self.right_sidebar_widget.setVisible(True) # Show before animation

        # Animate sidebar geometry (position and width)
        self.right_sidebar_animation = QPropertyAnimation(self.right_sidebar_widget, b"geometry")
        self.right_sidebar_animation.setDuration(ANIMATION_DURATION)
        self.right_sidebar_animation.setStartValue(start_geometry)
        self.right_sidebar_animation.setEndValue(end_geometry)
        self.right_sidebar_animation.setEasingCurve(QEasingCurve.InOutQuad)
        if not target_visible: # Hide after closing
            self.right_sidebar_animation.finished.connect(lambda: self.right_sidebar_widget.setVisible(False))

        self.right_sidebar_animation.start()
        self._update_sidebar_button_icon(self.plot_controls_collapse_button, target_visible, is_left=False)
        # Button position doesn't change, no need to animate it 

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

            # Update icons based on current *actual* visibility state of sidebars
            self._update_sidebar_button_icon(self.session_collapse_button, self.left_sidebar_widget.isVisible(), is_left=True)
            self._update_sidebar_button_icon(self.plot_controls_collapse_button, self.right_sidebar_widget.isVisible(), is_left=False)

            # Ensure sidebars are potentially visible (geometry/animation handles appearance)
            # No setVisible calls needed here, _reposition handles geometry         

        else:
            # --- State when NO valid session is loaded ---
            if self.central_stack.currentIndex() != 1: self.central_stack.setCurrentIndex(1)

            # Hide buttons
            self.session_collapse_button.setVisible(False)
            self.plot_controls_collapse_button.setVisible(False)
            # Hide sidebars directly (no animation needed when invalid)
            if self.left_sidebar_widget.isVisible(): self.left_sidebar_widget.setVisible(False)
            if self.right_sidebar_widget.isVisible(): self.right_sidebar_widget.setVisible(False)

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
