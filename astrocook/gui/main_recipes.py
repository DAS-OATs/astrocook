"""The recipes run from the main window: the requirements (z_em,
resolution) and the recipe resumed after Set Properties, the recipe
dialog, the checks before a recipe runs (export of ASCII files, overwrite
of a variable resolution, co-add with mixed resolutions), the progress
dialog and the worker, the result (log, warnings, cancellation), errors
and Stop.

Moved out of ``main_window.py`` in phase 5 (split of the GUI), as a mixin
of ``MainWindowV2``.
"""
import logging
import os
from typing import Any, Dict

from PySide6.QtCore import QTimer, Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (QApplication, QCheckBox, QDialog, QFileDialog, QLineEdit,
                               QMessageBox, QProgressDialog)

from astrocook.core.resolution import resolution_summary
from astrocook.core.structures import HistoryLogV2
from astrocook.core.utils import is_branching_recipe, resource_path

from .qt_workers import RecipeWorker
from .recipe_dialog import RecipeDialog
from .recipe_progress import (LogSignalBridge, QtLogHandler, RecipeProgressDialog,
                              WarningCaptureHandler)


class MainRecipesMixin:
    """The recipes of ``MainWindowV2``.

    Uses the window's ``active_history``, ``session_histories``,
    ``thread_pool``, ``progress_dialog``, ``active_recipe_dialog``, the
    pending recipe and the refit warning flag (set in ``__init__``),
    ``_show_custom_message``, ``plot_viewer`` and the session machinery
    (``update_gui_session_state``, ``_update_view_for_session``).
    """

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

            # The files the export writes, as the recipe decides them (F82:
            # decided here from the internal names of the columns only, an
            # existing spectrum file could be overwritten without asking)
            outputs = self.active_history.current_state.file.export_ascii_outputs(raw_targets)
            has_spec = 'spec' in outputs
            has_systems = 'systems' in outputs
        
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
