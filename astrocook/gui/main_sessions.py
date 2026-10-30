"""The sessions of the main window: opening them (file dialog, path, drag
and drop), the session list (switching, renaming, its context menu), the
view of a session, a recipe result added to the history or as a new
branch, duplicating, stitching and co-adding, the session inspector,
saving and closing.

Moved out of ``main_window.py`` in phase 5 (split of the GUI), as a mixin
of ``MainWindowV2``.
"""
import logging
import os
from copy import deepcopy
from typing import List, Optional

import numpy as np
from PySide6.QtCore import QItemSelectionModel, QLocale, QPoint, QSettings, QStringListModel, Qt
from PySide6.QtGui import QAction, QDoubleValidator, QPalette
from PySide6.QtWidgets import (QApplication, QDialog, QDialogButtonBox, QFileDialog,
                               QFormLayout, QLabel, QLineEdit, QMenu, QMessageBox,
                               QSizePolicy, QVBoxLayout)

from astrocook.core.resolution import resolution_summary
from astrocook.core.session import SessionV2, load_session_from_file
from astrocook.core.session_manager import SessionHistory
from astrocook.core.structures import HistoryLogV2

from .theme import FOCUS_BORDER, palette_hex


class SessionListModel(QStringListModel):
    """The names of the sessions, with the swatch of each other session
    shown in the plot next to its name (the plot has no legend, F86)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._swatches = {}  # row -> QIcon
        self._blank = None

    def set_swatches(self, swatches: dict, blank):
        """``swatches``: row -> icon of the sessions shown; with some, the
        other rows get ``blank`` so that the names stay aligned."""
        if swatches == self._swatches:
            return
        self._swatches, self._blank = swatches, blank
        if self.rowCount():
            # Only the decoration: the names are untouched (dataChanged on
            # the names renames the sessions, _on_session_name_changed)
            self.dataChanged.emit(self.index(0), self.index(self.rowCount() - 1),
                                  [Qt.ItemDataRole.DecorationRole])

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DecorationRole:
            if not self._swatches:
                return None
            return self._swatches.get(index.row(), self._blank)
        return super().data(index, role)


class MainSessionsMixin:
    """The sessions of ``MainWindowV2``.

    Uses the window's ``session_histories``, ``active_history``, the
    session list and its model, ``plot_viewer``, the plot controls (the
    auxiliary and SNR columns, Normalize F), the recipe machinery
    (``_on_recipe_requested``, ``_launch_recipe_dialog``), the log
    scripter and the identification viewer, and ``_update_ui_state``.
    """

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
        dialog.setWindowTitle(f"Session Info: {history_item.display_name}")
        dialog.setMinimumWidth(400)
        
        # A real f-string: it was a plain string with the placeholders of
        # one, not valid CSS (findings register F87)
        dialog.setStyleSheet(f"""
            QLineEdit {{
                padding: 3px;
                border-radius: 5px;
                background-color: {palette_hex(QPalette.ColorRole.Base)};
                color: {palette_hex(QPalette.ColorRole.Text)};
            }}
            QLineEdit:focus {{
                {FOCUS_BORDER}
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
        form_layout.addRow("Session name:", self.info_name_edit)
        
        # Object
        current_obj = meta.get('OBJECT', '')
        obj_display = str(current_obj) if current_obj else "undefined"
        self.info_object_edit = _expand_field(QLineEdit(obj_display))
        form_layout.addRow("Object:", self.info_object_edit)
        
        # z_em
        self.info_zem_edit = _expand_field(QLineEdit(f"{spec._data.z_em:.5f}"))
        dbl_validator = QDoubleValidator(); dbl_validator.setLocale(QLocale.C)
        self.info_zem_edit.setValidator(dbl_validator)
        form_layout.addRow("Emission redshift (z_em):", self.info_zem_edit)

        # Resolution: median R of the resol column (#526)
        median_r, variable_r = resolution_summary(spec)
        current_resol_str = f"{median_r:.0f}" if median_r > 0 else "undefined"
        self._info_resol_shown = current_resol_str
            
        self.info_resol_edit = _expand_field(QLineEdit(current_resol_str))
        self.info_resol_edit.setPlaceholderText("e.g. 50000")
        if variable_r:
            self.info_resol_edit.setToolTip(
                "The resolution varies with wavelength: this is its median. "
                "A new value replaces it only if you confirm.")
        form_layout.addRow("Resolving power (R, median):" if variable_r else "Resolving power (R):",
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
                        # In velocity first: the step of a logarithmic grid
                        step_str = f"{dv:.2f} km/s (logarithmic grid, {med_dx:.4f} {spec.x.unit})"
                    else:
                        step_str = f"{med_dx:.4f} {spec.x.unit} (variable)"

        chi2_html = "<tr><td>Reduced χ²:</td><td>N/A (no model)</td></tr>"
        
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
                            
                        chi2_html = f"<tr><td>Reduced χ² (forest / red side):</td><td>{c_forest_str} / {c_red_str}</td></tr>"
                    else:
                        # Fallback if z_em is not set
                        c_glob = np.sum(chi2_array) / len(chi2_array)
                        chi2_html = f"<tr><td>Reduced χ²:</td><td>{c_glob:.2f}</td></tr>"
            except Exception as e:
                logging.warning(f"Could not calculate chi-squared for info panel: {e}")

        stats_html = (
            "<p style='margin-bottom: 5px;'><b>Statistics:</b></p>"
            "<style>td { padding-left: 15px; padding-top: 5px; }</style>"
            "<table width='100%'>"
            f"<tr><td>Wavelength range:</td><td>{x_range}</td></tr>"
            f"<tr><td>Median step:</td><td>{step_str}</td></tr>" # Added Row
            f"<tr><td>Median S/N:</td><td>{snr_str}</td></tr>"
            f"{chi2_html}"  # <--- [NEW] Inject the Chi-Squared row here
            f"<tr><td>Pixels:</td><td>{len(spec.x)}</td></tr>"
            f"<tr><td>Components:</td><td>{len(systs.components)}</td></tr>"
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

    def _on_session_name_changed(self, index_top_left, index_bottom_right, roles=()):
        """ 
        Slot called when the user renames a session directly in the sidebar list.
        """
        if roles and not ({Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.EditRole} & set(roles)):
            return  # Only the swatches changed (SessionListModel)
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
