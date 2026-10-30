import logging
import os
from PySide6.QtCore import (
    Qt, QThreadPool
)
from PySide6.QtGui import (
    QIcon, QPixmap
)
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QDialog, QMainWindow, QMessageBox, QProgressDialog, QStackedWidget, QTextEdit,
)
from typing import List, Optional, Union

from .continuum_controls import ContinuumControlsMixin
from .debug_utils import GLOBAL_PLOTTER
from .identification_viewer_dialog import IdentificationViewerDialog
from .log_scripter_dialog import LogScripterDialog
from .main_menus import MainMenusMixin
from .main_panels import MainPanelsMixin
from .main_recipes import MainRecipesMixin
from .main_scripts import MainScriptsMixin
from .main_sessions import MainSessionsMixin, SessionListModel
from .script_dialog import ScriptDialog
from astrocook.core.debug_hooks import set_debug_plot_handler
from astrocook.core.session_manager import SessionHistory
from astrocook.core.session import SessionV2, LogManager
from astrocook.core.utils import resource_path # Import recipe helpers

# --- *** RECIPE LOOKUP MAP *** ---
from astrocook.recipes.absorbers import ABSORBERS_RECIPES_SCHEMAS
from astrocook.recipes.continuum import CONTINUUM_RECIPES_SCHEMAS
from astrocook.recipes.edit import EDIT_RECIPES_SCHEMAS
from astrocook.recipes.file import FILE_RECIPES_SCHEMAS
from astrocook.recipes.flux import FLUX_RECIPES_SCHEMAS
RECIPE_CATEGORY_MAP = {}
for name in ABSORBERS_RECIPES_SCHEMAS:
    RECIPE_CATEGORY_MAP[name] = 'absorbers'
for name in CONTINUUM_RECIPES_SCHEMAS: 
    RECIPE_CATEGORY_MAP[name] = 'continuum'
for name in EDIT_RECIPES_SCHEMAS:
    RECIPE_CATEGORY_MAP[name] = 'edit'
for name in FILE_RECIPES_SCHEMAS:
    RECIPE_CATEGORY_MAP[name] = 'file'
for name in FLUX_RECIPES_SCHEMAS:
    RECIPE_CATEGORY_MAP[name] = 'flux'
# (Add other recipe categories here as they are created)

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


class MainWindowV2(MainMenusMixin, MainPanelsMixin, ContinuumControlsMixin, MainSessionsMixin,
                   MainRecipesMixin, MainScriptsMixin, QMainWindow):
    def __init__(self, initial_sessions: Union[SessionV2, List[SessionV2], None], initial_log_object: Optional[LogManager]):
        super().__init__()
        # Forces the app name to be "Astrocook" instead of "Python"
        app = QApplication.instance()
        if app:
            app.setApplicationName("Astrocook")
            app.setApplicationDisplayName("Astrocook")

        self.session_histories: List[SessionHistory] = [] # List of history managers
        self.active_history: Optional[SessionHistory] = None # Reference to the selected manager
        self.session_model = SessionListModel()
        
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
            icon_path = resource_path(os.path.join("assets", "logo_3d_HR.png"))
            self.setWindowIcon(QIcon(icon_path))
        except Exception as e:
            logging.warning(f"Could not set window icon: {e}")

        # --- ** Central Widget is NOW the Stack ** ---
        self.central_stack = QStackedWidget()
        
        # --- ** Sidebars are Children of the MainWindow, floating above ** ---
        self._setup_left_sidebar()  # Creates self.left_sidebar_widget
        self._setup_right_sidebar() # Creates self.right_sidebar_widget
        self._setup_collapse_buttons()
        self._sessions_fader.add_handle(self.session_collapse_button)
        self._plot_controls_fader.add_handle(self.plot_controls_collapse_button)

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
                
                log_object = self._log_for_new_session(session)
                new_history = SessionHistory(session, log_object)
                
                # Link real GUI
                session._gui = self
                
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
        # The debug plots requested by the core (from worker threads) reach
        # the main thread through the signal (F6)
        set_debug_plot_handler(GLOBAL_PLOTTER.plot_requested.emit)
        GLOBAL_PLOTTER.plot_requested.connect(self._on_debug_plot_requested)

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
            /* Floating cards, as tall as their contents and flush with the
               window edge: rounded on the inner side, no border */
            QWidget#SessionContainer {{
                background-color: rgba({sidebar_r}, {sidebar_g}, {sidebar_b}, 0.95);
                border: none;
                border-top-right-radius: 8px;
                border-bottom-right-radius: 8px;
            }}
            QWidget#PlotControlsContainer {{
                background-color: rgba({sidebar_r}, {sidebar_g}, {sidebar_b}, 0.95);
                border: none;
                border-top-left-radius: 8px;
                border-bottom-left-radius: 8px;
            }}
            
            /* Make ListView background fully transparent */
            QListView#SessionListView {{ background-color: transparent; border: none; }}
            /* Ensure list items are NOT transparent */
            QListView#SessionListView::item {{
                padding: 10px 10px; border: none;
                background-color: transparent; /* Inherit sidebar bg initially */
                color: {button_fg}; /* Use button text color for readability */
             }}
            /* Selected/Hover items should be opaque */
            QListView#SessionListView::item:selected {{ background-color: {item_selected_bg}; color: {item_selected_text}; border-radius: 6px; }}

            /* Right Sidebar Checkboxes & Labels */
            QWidget#PlotControlsContainer QLabel {{ /* Style section labels */
                margin-bottom: 4px; /* Space below label */
            }}
            QCheckBox#PlotControlCheckbox {{
                color: {text_color}; spacing: 5px; padding: 2px 0px; /* Reduced padding */
                padding-left: 15px; /* Indentation */
                /* No background-color here: with the non-native styles
                   (Fusion: Linux, Windows) it also clears the frame of the
                   indicator, and unchecked boxes become invisible */
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
            /* Collapse tabs: handles with the background of the cards,
               rounded on the side away from their card */
            QPushButton#CollapseButton {{
                background-color: rgba({sidebar_r}, {sidebar_g}, {sidebar_b}, 0.95);
                color: {button_fg}; /* Icon color */
                border: none;
                margin: 0px; padding: 0px;
                icon-size: 18px;
            }}
            QPushButton#CollapseButton[side="left"] {{
                border-top-right-radius: 6px; border-bottom-right-radius: 6px;
            }}
            QPushButton#CollapseButton[side="right"] {{
                border-top-left-radius: 6px; border-bottom-left-radius: 6px;
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
                             checkbox_text=None,
                             detailed_text=None):
        target_parent = parent if parent else self
        
        msg_box = QMessageBox(target_parent)
        msg_box.setWindowTitle(title)

        # Convert standard \n line breaks to HTML <br> tags
        formatted_body = text.replace('\n', '<br>')
        
        formatted_text = (
            f"<p style='font-size: 13pt; font-weight: bold; margin: 0px;'>{header}</p>"
            f"<p style='font-size: 12pt; font-weight: normal; margin-top: 8px;'>{formatted_body}</p>"
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
        
        msg_box.setStyleSheet("""
            QLabel { padding-top: 10px; }
            QCheckBox { margin-top: 10px; font-size: 12pt; }
        """)
        
        # --- [NEW] Permanently visible scrollable box ---
        if detailed_text:
            details_box = QTextEdit(msg_box)
            details_box.setReadOnly(True)
            details_box.setPlainText(detailed_text)
            details_box.setMinimumWidth(500)
            details_box.setMinimumHeight(200)
            details_box.setMaximumHeight(400)
            # Add some spacing above it and a slightly darker background to inset it
            #details_box.setStyleSheet("""
            #    QTextEdit {
            #        font-family: monospace; 
            #        margin-top: 15px; 
            #        background-color: #f4f5f7;
            #        border: 1px solid #dcdcdc;
            #        border-radius: 5px;
            #        padding: 5px;
            #    }
            #    QScrollBar:vertical {
            #        border: none;
            #        background: transparent;
            #        width: 10px;
            #        margin: 2px;
            #    }
            #    QScrollBar::handle:vertical {
            #        background: #c1c1c1;
            #        min-height: 30px;
            #        border-radius: 4px;
            #    }
            #    QScrollBar::handle:vertical:hover {
            #        background: #a8a8a8;
            #    }
            #    QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
            #        height: 0px; /* Hides the ugly up/down arrow buttons */
            #    }
            #    QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {
            #        background: none; /* Removes the dark background track */
            #    }
            #""")

            # Set the monospace font natively without invoking the CSS engine
            font = details_box.font()
            font.setFamily("monospace")
            details_box.setFont(font)

            # Inject it into the QMessageBox's internal grid layout
            grid = msg_box.layout()
            # Add it at the bottom, spanning all columns
            grid.addWidget(details_box, grid.rowCount(), 0, 1, grid.columnCount())
            
        # Add Checkbox if requested
        cb = None
        if checkbox_text:
            cb = QCheckBox(checkbox_text)
            msg_box.setCheckBox(cb)
        
        ret = msg_box.exec()
        
        if cb:
            return ret, cb.isChecked()
        return ret
     

    # --- ** NEW Toggle & Animation Methods ** ---

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
