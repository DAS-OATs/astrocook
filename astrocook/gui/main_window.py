from copy import deepcopy
import logging
import numpy as np
import os
from PySide6.QtCore import (
    Qt, QItemSelectionModel,
    QLocale,
    QPoint,
    QStringListModel, QThreadPool, QTimer, QSettings
)
from PySide6.QtGui import (
    QAction, QDoubleValidator, QIcon, QPixmap
)
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QDialog, QDialogButtonBox, QFileDialog, QMainWindow, QVBoxLayout, QFormLayout, QLabel, QLineEdit, QMenu, QMessageBox, QProgressDialog, QSizePolicy, QStackedWidget, QTextEdit,
)
import re
from typing import Any, Dict, List, Optional, Union

from .continuum_controls import ContinuumControlsMixin
from .debug_utils import GLOBAL_PLOTTER
from .identification_viewer_dialog import IdentificationViewerDialog
from .log_scripter_dialog import LogScripterDialog
from .main_menus import MainMenusMixin
from .main_panels import MainPanelsMixin
from .qt_workers import RecipeWorker, ScriptWorker
from .recipe_dialog import RecipeDialog
from .recipe_progress import (LogSignalBridge, QtLogHandler, RecipeProgressDialog,
                              WarningCaptureHandler)
from .script_dialog import ScriptDialog
from astrocook.core.debug_hooks import set_debug_plot_handler
from astrocook.core.resolution import resolution_summary
from astrocook.core.session_manager import SessionHistory
from astrocook.core.session import SessionV2, load_session_from_file, LogManager
from astrocook.core.structures import HistoryLogV2
from astrocook.core.utils import is_branching_recipe, resource_path # Import recipe helpers

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


class MainWindowV2(MainMenusMixin, MainPanelsMixin, ContinuumControlsMixin, QMainWindow):
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
            logging.debug("Branching: Creating new SessionHistory from current state.")
            # 1. Get the source log and copy it (HistoryLogV2 is safe to deepcopy)
            source_log_manager = target_history.log_manager
            try:
                new_log_copy = deepcopy(source_log_manager)
            except Exception as e:
                logging.error(f"Standard deepcopy failed on log manager: {e}")
                new_log_copy = None

            if new_log_copy is None:
                logging.error("Failed to copy the session log for branching, aborting.")
                QMessageBox.critical(self, "Branching Error", "Could not copy session log state.")
                return

            # 2. Update the new state from the recipe to point to the *new* log
            new_session.log_manager = new_log_copy


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
            logging.debug("Linear update: Adding state to active SessionHistory.")
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
    
    def _log_for_new_session(self, session):
        """Return the log for a new session history.

        A V2 log restored from an ``.acs2`` archive is kept; any other session
        (including a duplicate, whose log belongs to another history) starts
        with a new, empty log.
        """
        log = getattr(session, 'log_manager', None)
        histories = getattr(self, 'session_histories', [])
        if isinstance(log, HistoryLogV2) and not any(
                h.log_manager is log for h in histories):
            return log
        return HistoryLogV2()

    def _add_session_internal(self, new_session, is_initial=False):
        """Creates and adds a new SessionHistory manager."""
        if new_session is None or isinstance(new_session, int): return

        log_object = self._log_for_new_session(new_session)
        
        # 1. Always create a new history for a newly loaded session
        new_history = SessionHistory(new_session, log_object)
        
        # Overwrite the mock GUI context with the real one
        new_session._gui = self
        
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

        # --- Safely uncheck Normalize F if no continuum exists ---
        if is_valid and session_state_to_show and session_state_to_show.spec:
            if session_state_to_show.spec.cont is None:
                if self.norm_y_checkbox.isChecked():
                    self.norm_y_checkbox.blockSignals(True)
                    self.norm_y_checkbox.setChecked(False)
                    if hasattr(self, 'lbl_flux'):
                        self.lbl_flux.setText("F:")
                    self.norm_y_checkbox.blockSignals(False)

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

    def _on_session_switched(self, index):
        """Sets the active SessionHistory manager when the list view is clicked."""
        try:
            new_history_list_index = index.row()
            if 0 <= new_history_list_index < len(self.session_histories):
                new_active_history = self.session_histories[new_history_list_index]
                if self.active_history != new_active_history:
                    # Only switch if user confirms discarding edits
                    if not self._confirm_discard_continuum_edits():
                        # Re-select the old item in the list visually
                        self._update_view_for_session(self.active_history.current_state, set_current_list_item=True)
                        return

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
        # [Check for pending edits BEFORE opening the file browser
        if not self._confirm_discard_continuum_edits():
            return
        
        # --- Retrieve the last used directory ---
        settings = QSettings("Astrocook", "AstrocookV2")
        last_dir = settings.value("last_working_dir", os.getcwd())

        file_names, _ = QFileDialog.getOpenFileNames(
            self, 
            "Open Spectrum File(s)", 
            last_dir,
            "All Supported (*.acs *.acs2 *.fits *.txt *.dat);;Astrocook Sessions (*.acs *.acs2);;FITS Files (*.fits);;Text Files (*.txt *.dat);;All Files (*)"
        )
        
        if file_names:
            # --- Save the directory for next time ---
            settings.setValue("last_working_dir", os.path.dirname(file_names[0]))
            
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

            stitch_eq_action = QAction(f"Equalize and Stitch {len(indexes)} Arms...", self)
            stitch_eq_action.setToolTip("Equalize flux levels and stitch arms at specific wavelengths.")
            stitch_eq_action.triggered.connect(
                lambda: self._on_equalized_stitch_dialog_requested(primary_history, others_names)
            )
            menu.addAction(stitch_eq_action)

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
        info_action = QAction("View Info", self)
        info_action.triggered.connect(
            # Use lambda to pass the specific history item
            lambda checked=False, item=history_item: self._on_session_info(item)
        )
        menu.addAction(info_action)

        # --- View Log action ---
        view_action = QAction("View Log", self)
        view_action.triggered.connect(
            lambda: self._launch_log_scripter(history_item)
        )
        menu.addAction(view_action)

        view_ids_action = QAction("View Identifications", self)
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
        save_action = QAction("Save as...", self)
        save_action.triggered.connect(
            lambda: self._on_save_session_context(history_item)
        )
        menu.addAction(save_action)

        # --- Close Session Action ---
        close_action = QAction("Close", self)
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

    def _on_equalized_stitch_dialog_requested(self, primary_history, others_names: List[str]):
        """
        Callback for the Equalize and Stitch context menu action.
        Pre-fills the 'other_sessions' in the dialog.
        """
        if self.active_history != primary_history:
            self.active_history = primary_history
            self._update_view_for_session(primary_history.current_state, set_current_list_item=True)

        # N total sessions selected means we need exactly N-1 cuts.
        # Since 'others_names' excludes the primary session, its length IS exactly N-1!
        n_minus_one = len(others_names)

        # Pluralization check
        val_word = "value" if n_minus_one == 1 else "values"

        # 1. Dynamically inject the exact number into the Recipe Schema in memory
        from astrocook.recipes.edit import EDIT_RECIPES_SCHEMAS
        params = EDIT_RECIPES_SCHEMAS["equalize_and_stitch"]["params"]
        
        for p in params:
            if p["name"] == "stitch_wavelengths":
                p["doc"] = f"Comma-separated cut-offs (requires {n_minus_one} {val_word})."
            elif p["name"] == "equalize_ranges":
                p["doc"] = f"Overlap ranges 'min-max' (requires {n_minus_one} {val_word}) or 'auto'."
            elif p["name"] == "manual_factors":
                p["doc"] = f"Fixed multipliers (requires {n_minus_one} {val_word}) or 'auto'."

        others_str = ", ".join(others_names)
        
        # 3. Store the hidden string so we can inject it later
        self._hidden_stitch_sessions = others_str

        self._launch_recipe_dialog(
            "edit", "equalize_and_stitch"
        )

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

        # Resolution: median R of the resol column (#526)
        median_r, variable_r = resolution_summary(spec)
        current_resol_str = f"{median_r:.0f}" if median_r > 0 else "undefined"
        self._info_resol_shown = current_resol_str
            
        self.info_resol_edit = _expand_field(QLineEdit(current_resol_str))
        self.info_resol_edit.setPlaceholderText("Resolving Power R")
        if variable_r:
            self.info_resol_edit.setToolTip(
                "The resolution varies with wavelength: this is its median. "
                "A new value replaces it only if you confirm.")
        form_layout.addRow("Resolution (R, median):" if variable_r else "Resolution (R):",
                           self.info_resol_edit)
        
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

        chi2_html = "<tr><td>Red. χ² (Global):</td><td>N/A (No model)</td></tr>"
        
        if spec.model is not None and len(spec.model.value) > 0:
            try:
                x_vals = spec.x.value
                y_vals = spec.y.value
                dy_vals = spec.dy.value
                mod_vals = spec.model.value
                
                # Filter out NaNs and division-by-zero
                valid = (dy_vals > 0) & np.isfinite(y_vals) & np.isfinite(mod_vals) & np.isfinite(dy_vals)
                
                if np.any(valid):
                    x_v = x_vals[valid]
                    y_v = y_vals[valid]
                    dy_v = dy_vals[valid]
                    mod_v = mod_vals[valid]
                    
                    # Calculate pixel-by-pixel chi-squared
                    chi2_array = ((y_v - mod_v) / dy_v) ** 2
                    
                    z_em = spec._data.z_em
                    
                    if z_em > 0:
                        # Split by Ly-alpha (121.567 nm)
                        lya_obs = 121.567 * (1.0 + z_em)
                        
                        forest_mask = x_v < lya_obs
                        red_mask = x_v >= lya_obs
                        
                        c_forest_str = "N/A"
                        c_red_str = "N/A"
                        
                        if np.any(forest_mask):
                            c_forest = np.sum(chi2_array[forest_mask]) / np.sum(forest_mask)
                            c_forest_str = f"{c_forest:.2f}"
                            
                        if np.any(red_mask):
                            c_red = np.sum(chi2_array[red_mask]) / np.sum(red_mask)
                            c_red_str = f"{c_red:.2f}"
                            
                        chi2_html = f"<tr><td>Red. χ² (forest / red side):</td><td>{c_forest_str} / {c_red_str}</td></tr>"
                    else:
                        # Fallback if z_em is not set
                        c_glob = np.sum(chi2_array) / len(chi2_array)
                        chi2_html = f"<tr><td>Red. χ² (global):</td><td>{c_glob:.2f}</td></tr>"
            except Exception as e:
                logging.warning(f"Could not calculate chi-squared for info panel: {e}")

        stats_html = (
            "<p style='margin-bottom: 5px;'><b>Other info:</b></p>"
            "<style>td { padding-left: 15px; padding-top: 5px; }</style>"
            "<table width='100%'>"
            f"<tr><td>Wavelength Range:</td><td>{x_range}</td></tr>"
            f"<tr><td>Median Step:</td><td>{step_str}</td></tr>" # Added Row
            f"<tr><td>Median SNR:</td><td>{snr_str}</td></tr>"
            f"{chi2_html}"  # <--- [NEW] Inject the Chi-Squared row here
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
        if raw_resol == getattr(self, '_info_resol_shown', None):
            new_resol = '_current_'  # not edited: leave the resolution as it is
        elif raw_resol.lower() == "undefined" or not raw_resol:
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
        
        # --- Sync with the last used directory ---
        settings = QSettings("Astrocook", "AstrocookV2")
        last_dir = settings.value("last_working_dir", os.getcwd())
        default_path = os.path.join(last_dir, session_to_save_final.name + ".acs2")

        file_name, selected_filter = QFileDialog.getSaveFileName(
            self,
            "Save Astrocook Session",
            default_path,
            "Astrocook V2 Session (*.acs2)"
        )

        if file_name:
            # --- Save the directory for next time ---
            settings.setValue("last_working_dir", os.path.dirname(file_name))

            logging.info(f"Saving session to: {file_name}")
            try:
                # Pass *both* sessions to the save method
                session_to_save_final.save(
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
        if history_to_close is None or isinstance(history_to_close, bool):
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
            # The resol column, the single source (#526)
            has_resol = resolution_summary(spec)[0] > 0
            
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
                
                median_r = resolution_summary(spec)[0]
                current_params['resol'] = f"{median_r:.0f}" if median_r > 0 else "0"
            
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

        # --- 1. THE INTERCEPT (MUST BE AT THE TOP) ---
        if recipe_name == "export_ascii":
            # 1. CLEAN PARSING OF TARGETS
            raw_targets = params.get('targets', '')
            target_list = [t.strip() for t in raw_targets.split(',') if t.strip()]
            
            # Authoritative identification
            # A 'spec' target is anything in the spectrum table (excluding 'systems')
            spec_colnames = self.active_history.current_state.spec.t.colnames
            has_spec = any(t in spec_colnames or t == "spec" for t in target_list if t != "systems")
            
            # A 'systems' target is the literal 'systems' or any column in the systems table
            syst_colnames = self.active_history.current_state.systs.t.colnames
            has_systems = "systems" in target_list or any(t in syst_colnames for t in target_list)
        
            # 2. SMART FILENAME SUGGESTION (REFINED)
            session_name = self.active_history.current_state.name
            if has_spec and not has_systems:
                suggested_name = f"{session_name}_spec.csv"
            elif has_systems and not has_spec:
                suggested_name = f"{session_name}_systems.csv"
            else:
                # Mixed export: Suggested name is just the prefix (no .csv)
                suggested_name = session_name
            
            full_path, _ = QFileDialog.getSaveFileName(
                self, "Export ASCII Table", suggested_name, 
                "CSV Files (*.csv);;Text Files (*.txt);;All Files (*)"
            )
            
            if not full_path: return
        
            import os
            folder = os.path.dirname(full_path)
            # Extract prefix and strip redundant suffixes if user clicked an existing file
            filename_prefix = os.path.splitext(os.path.basename(full_path))[0]
            filename_prefix = filename_prefix.replace('_spec', '').replace('_systems', '')
        
            # 3. OVERWRITE SAFETY CHECK
            conflicts = []
            if has_spec:
                spec_path = os.path.join(folder, f"{filename_prefix}_spec.csv")
                if os.path.exists(spec_path): conflicts.append(os.path.basename(spec_path))
            if has_systems:
                syst_path = os.path.join(folder, f"{filename_prefix}_systems.csv")
                if os.path.exists(syst_path): conflicts.append(os.path.basename(syst_path))
        
            if conflicts:
                conflict_list = "\n".join([f"• {c}" for c in conflicts])
                ret = self._show_custom_message(
                    title="Confirm Overwrite",
                    header="Files already exist",
                    text=f"The following files will be overwritten:\n\n{conflict_list}\n\nDo you want to proceed?",
                    buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    default_btn=QMessageBox.StandardButton.No
                )
                if ret != QMessageBox.StandardButton.Yes: return
        
            # 4. EXECUTION
            try:
                saved = self.active_history.current_state.file.export_ascii(
                    targets=raw_targets, 
                    path=folder, 
                    prefix=filename_prefix
                )
                
                if saved:
                    files_created = []
                    if has_spec: files_created.append(f"• {filename_prefix}_spec.csv")
                    if has_systems: files_created.append(f"• {filename_prefix}_systems.csv")
                    
                    self._show_custom_message(
                        title="Export Complete",
                        header="Data exported successfully",
                        text=f"Saved in {folder}:\n\n" + "\n".join(files_created)
                    )
            except Exception as e:
                logging.error(f"Export failed: {e}")
            return

        # --- [NEW] GUI-Thread Pre-checks ---
        if recipe_name == "set_properties":
            new_resol = str(params.get('resol', '_current_'))
            # If the user is actively changing the resolution
            if new_resol != '_current_':
                # Ask only if the resolution varies with wavelength: a flat
                # or missing one is always replaced (#526)
                if (self.active_history
                        and resolution_summary(self.active_history.current_state.spec)[1]):
                    from PySide6.QtWidgets import QMessageBox, QApplication
                    
                    # Find the currently active floating window/dialog
                    active_win = QApplication.activeWindow()
                    if active_win is None:
                        active_win = self
                        
                    reply = self._show_custom_message(
                        title="Overwrite Resolution Column?",
                        header="The resolution varies with wavelength.",
                        text="This spectrum has a resolution that varies with wavelength (the 'resol' column).\n\nDo you want to replace it with a single value?",
                        buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                        default_btn=QMessageBox.StandardButton.No,
                        parent=active_win  # <--- [FIX] Parent to the active floating window!
                    )
                    # Inject the user's choice into the params before sending to the worker
                    params['overwrite_resol_col'] = (reply == QMessageBox.StandardButton.Yes)

        if recipe_name == "coadd":
            names_str = params.get('session_names', '')
            names_list = [n.strip() for n in names_str.split(',') if n.strip()]
            
            if names_list and hasattr(self, 'session_histories'):
                # Find the actual Session objects requested
                specs_to_coadd = []
                for hist in self.session_histories:
                    if hist.display_name in names_list:
                        specs_to_coadd.append(hist.current_state.spec)
                
                # Count how many have defined resolutions
                if specs_to_coadd:
                    resol_count = 0
                    for s in specs_to_coadd:
                        if resolution_summary(s)[0] > 0:  # #526
                            resol_count += 1
                            
                    # Mismatch detected!
                    if 0 < resol_count < len(specs_to_coadd):
                        from PySide6.QtWidgets import QMessageBox, QApplication
                        active_win = QApplication.activeWindow()
                        if active_win is None: active_win = self
                            
                        reply = self._show_custom_message(
                            title="Resolution Mismatch",
                            header="Missing Resolution Data",
                            text=f"Only {resol_count} out of {len(specs_to_coadd)} selected sessions have a defined resolution.\n\nMixing known and unknown resolutions means the final LSF cannot be calculated accurately. The co-added spectrum will NOT have a resolution column.\n\nDo you want to proceed anyway?",
                            buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                            default_btn=QMessageBox.StandardButton.No,
                            parent=active_win
                        )
                        if reply == QMessageBox.StandardButton.No:
                            return  # Abort the recipe dispatch!

        # --- 2. STANDARD RECIPE PREPARATION ---

        # [FIX] Re-inject the hidden session names for equalize_and_stitch
        if recipe_name == "equalize_and_stitch" and hasattr(self, '_hidden_stitch_sessions'):
            params['other_sessions'] = self._hidden_stitch_sessions
            del self._hidden_stitch_sessions
        
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
            'refit_all', 'metals_auto', 'add_component', 'fit_component', 'fit_bayesian'
        }
        
        is_long_running = recipe_name in long_running_recipes

        if recipe_name == 'import_systems':
            # Check if refit is requested (handles both boolean True and string 'True')
            if str(params.get('refit', False)) == 'True':
                is_long_running = True

        if is_long_running:
            # A. Create Detailed Dialog
            self.progress_dialog = RecipeProgressDialog(f"{recipe_name}", self)

            # Determine if we should show an indeterminate bar.
            # 1. add_component and fit_component are always indeterminate.
            # 2. refit_all is indeterminate ONLY if selected_uuids is provided (Consolidated mode).
            is_indeterminate = recipe_name in {'add_component', 'fit_component'}
            
            if recipe_name == 'refit_all' and params.get('selected_uuids') is not None:
                is_indeterminate = True

            if is_indeterminate:
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

                # --- 2. LOGGING ---
                # A linear recipe is logged once on the parent. A branching
                # recipe is logged in each new session it opens, after the
                # branch has copied the parent log: the parent session does
                # not change, so its log does not either (F41, #575)
                branching = is_branching_recipe(recipe_name)
                params_to_log = params.copy()
                params_to_log.pop('alias_map', None)
                try:
                    if (i == 0 and not branching
                            and isinstance(original_parent_history.log_manager, HistoryLogV2)):
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

                # set_properties no longer writes the resolution into the
                # components (#526, F47): they follow the spectrum unless the
                # user set their own value
                original_history_index = self.session_histories.index(original_parent_history)
                
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

                branch_history = self.active_history
                try:
                    if (branching and branch_history is not original_parent_history
                            and isinstance(branch_history.log_manager, HistoryLogV2)):
                        branch_history.log_manager.add_entry(recipe_name, params_to_log)
                        if (self.log_scripter_dialog and self.log_scripter_dialog.isVisible()
                                and self.log_scripter_dialog.log_object is branch_history.log_manager):
                            self.log_scripter_dialog.refresh()
                except Exception: pass

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

        # [FIX] Clear the hidden string if the user cancelled the dialog
        if hasattr(self, '_hidden_stitch_sessions'):
            del self._hidden_stitch_sessions

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
            
        # 2. Build the scrollable details text (Traceback + Warnings)
        details_text = ""
        
        if trace:
            details_text += "Traceback:\n" + "-"*40 + "\n"
            # Truncate if massive to prevent Qt from lagging
            if len(trace) > 5000:
                details_text += trace[:5000] + "\n... (traceback truncated)"
            else:
                details_text += trace
            details_text += "\n\n"

        # Check for captured warnings and add them to the details box
        if hasattr(self, '_captured_warnings') and self._captured_warnings:
            unique_warnings = list(dict.fromkeys(self._captured_warnings))
            if unique_warnings:
                details_text += "Previous Warnings:\n" + "-"*40 + "\n"
                details_text += "\n".join(unique_warnings)

        if trace:
            # --- CASE A: System Crash (Show Traceback) ---
            self._show_custom_message(
                title="Internal Error",
                header=title, # e.g. "Recipe Error: fit_continuum"
                text=message,
                detailed_text=details_text if details_text else None, # <--- Both go here now
                buttons=QMessageBox.StandardButton.Ok,
                icon_name="icon_3d_HR.png" 
            )

        else:
            # --- CASE B: User Error (Offer Retry) ---
            ret = self._show_custom_message(
                title="Recipe Failed",
                header=title,   
                text=message,   
                detailed_text=details_text if details_text else None, # <--- Both go here now
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
        expression = f"(λ > {xmin:.4f}) & (λ < {xmax:.4f})"
        
        logging.info(f"Splitting from view: {expression}")

        # 3. Launch the dialog with these params pre-filled
        self._launch_recipe_dialog(
            "edit", 
            "split", 
            initial_params={"expression": expression}
        )

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
