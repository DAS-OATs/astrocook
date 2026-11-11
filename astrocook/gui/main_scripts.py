"""The scripts of the main window: the log scripter (the log of a session,
edited and run again), and a script run from the original data in place
of the history of the session.

Moved out of ``main_window.py`` in phase 5 (split of the GUI), as a mixin
of ``MainWindowV2``.
"""
import logging
from typing import TYPE_CHECKING, cast

from PySide6.QtCore import QTimer, Qt
from PySide6.QtWidgets import QApplication, QMessageBox, QProgressDialog

from astrocook.session.history import SessionHistory

from .log_scripter_dialog import LogScripterDialog
from .qt_workers import ScriptWorker

from astrocook.gui.errors import QT_ERRORS
from astrocook.gui._window_host import HostBase

if TYPE_CHECKING:
    from astrocook.gui.main_window import MainWindowV2
logger = logging.getLogger(__name__)


class MainScriptsMixin(HostBase):
    """The scripts of ``MainWindowV2``.

    Uses the window's ``active_history``, ``session_histories``, the session
    list's model, ``thread_pool``, ``progress_dialog``,
    ``log_scripter_dialog``, ``recipe_category_map``, ``plot_viewer`` and
    ``_safely_close_progress_dialog`` (``MainRecipesMixin``).
    """

    def _launch_log_scripter(self, history_object: 'SessionHistory'):
        """
        Creates and shows the *single, persistent* LogScripterDialog,
        or activates it if it already exists.
        """
        log_object_to_view = history_object.log_manager
        session_name = history_object.display_name

        if self.log_scripter_dialog:
            self.log_scripter_dialog.set_log_object(log_object_to_view) 
            self.log_scripter_dialog.setWindowTitle(f"Session Log: {session_name}")
            self.log_scripter_dialog.raise_()
            self.log_scripter_dialog.activateWindow()
            logger.debug(f"Activating existing log scripter for {session_name}")
            return

        dialog = LogScripterDialog(
            log_object_to_view, 
            session_name, 
            cast('MainWindowV2', self), 
            self.recipe_category_map # Pass the map
        )
        self.log_scripter_dialog = dialog
        dialog.finished.connect(
            lambda: self._on_log_scripter_closed() 
        )
        dialog.show()
        logger.debug(f"Creating new single log scripter for {session_name}")

    def _on_log_scripter_closed(self):
        """
        Slot called when the single LogScripterDialog is closed.
        Cleans up the reference.
        """
        self.log_scripter_dialog = None
        logger.debug("Single log scripter closed, reference removed.")
        
        try:
            self.session_collapse_button.raise_()
            self.plot_controls_collapse_button.raise_()
        except QT_ERRORS as e:
            logger.debug(f"Side buttons not raised: {e}")

    def run_script(self, script_text: str):
        """
        Launches the ScriptWorker on a background thread to run the script.
        """
        
        if not self.active_history:
            QMessageBox.warning(self, "No Session", "Cannot run script: No active session.")
            return
        if self._refuse_if_followed("Running a script"):
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
                self.plot_viewer.toolbar.deactivate_tools()
            except QT_ERRORS as e:
                logger.debug(f"Toolbar not reset: {e}")
        # --- *** END FIX 1 *** ---

        logger.info("Running the script.")

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
        worker = ScriptWorker(script_text, initial_state, cast('MainWindowV2', self))
        
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
            
        logger.error(f"The script stopped at line {line_num} ({line}): {error_msg}")
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
        
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        # Schedule the heavy GUI update
        QTimer.singleShot(10, lambda: self._process_script_result(new_history))

    def _process_script_result(self, new_history: SessionHistory):
        """The actual GUI update for 'Run All', called by a QTimer."""
        # 1. Swap the old history with the new one
        try:
            try:
                old_history_index = self.session_histories.index(cast(SessionHistory, self.active_history))
                self.session_histories[old_history_index] = new_history
                self.active_history = new_history
                self.session_model.setData(self.session_model.index(old_history_index), new_history.display_name)
            except ValueError:
                logger.warning("The session of the script was not found: the result is added as a new session.")
                self.session_histories.append(new_history)
                self.active_history = new_history
                self.session_model.setStringList([h.display_name for h in self.session_histories])

            # 2. Refresh everything
            logger.info("Script done.")
            self._update_view_for_session(self.active_history.current_state, set_current_list_item=True)
            self._update_undo_redo_actions()
            if self.log_scripter_dialog:
                self.log_scripter_dialog.set_log_object(new_history.log_manager)
        finally:
            # --- *** ALWAYS restore the cursor *** ---
            QApplication.restoreOverrideCursor()
