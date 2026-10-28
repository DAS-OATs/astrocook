"""The controls of the continuum editor in the plot controls: Start/Save
(the edited continuum saved as a new state of the session), Cancel and
Reset, and the spacing of the knots, on a slider with a log scale (100 to
20000 km/s) or typed.

Moved out of ``main_window.py`` in phase 5 (split of the GUI), as a mixin
of ``MainWindowV2``; the editor on the plot is ``ContinuumEditorMixin``
(``continuum_editor.py``).
"""
import logging

import numpy as np
from PySide6.QtWidgets import QMessageBox

from astrocook.recipes.continuum import CONTINUUM_RECIPES_SCHEMAS


class ContinuumControlsMixin:
    """The controls of the continuum editor of ``MainWindowV2``.

    Uses the window's ``plot_viewer`` (the editor itself), ``active_history``,
    the widgets of the Edit Continuum section and of the redshift cursor
    (made by ``MainPanelsMixin``), ``_show_custom_message``,
    ``_launch_recipe_dialog``, ``_ask_to_renormalize_model`` and the session
    list.
    """

    def _on_cancel_continuum_clicked(self):
        """Slot for the dedicated Cancel button in the sidebar."""
        # Re-use the confirmation logic. 
        # If the user clicks 'Cancel' in the dialog, the editor stays open.
        # If they click 'Save' or 'Discard', the editor closes.
        self._confirm_discard_continuum_edits()

    def _confirm_discard_continuum_edits(self) -> bool:
        """
        Checks if continuum editing is active and asks user to Save, Discard, or Stay.
        Returns True if the workflow (Open/Switch) can proceed.
        Returns False if the user wants to stay and continue editing.
        """
        if hasattr(self, 'plot_viewer') and getattr(self.plot_viewer, '_edit_mode_active', False):
            ret = self._show_custom_message(
                title="Continuum Editor Active",
                header="Unsaved Continuum Edits",
                text="You are currently editing the continuum. What would you like to do with your changes?",
                buttons=QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
                default_btn=QMessageBox.StandardButton.Save
            )
            
            if ret == QMessageBox.StandardButton.Save:
                # Trigger the same logic as the "Save" (formerly toggle) button
                self._toggle_continuum_editor() 
                return True
            elif ret == QMessageBox.StandardButton.Discard:
                # Discard and clean up visuals
                self.plot_viewer.cancel_continuum_edit()
                self._reset_continuum_ui()
                return True
            else: 
                # User clicked 'Cancel' (meaning 'Cancel the switch/open, I want to stay')
                return False
        return True

    def _reset_continuum_ui(self):
        """Helper to reset the continuum editor panel to its default visual state."""
        self.edit_cont_button.setText("Start")
        self.edit_cont_button.setStyleSheet("") 
        self.cancel_cont_button.setEnabled(False)
        self.stride_slider.setEnabled(False)
        self.stride_entry.setEnabled(False)
        self.reset_cont_button.setEnabled(False)
        self.edit_hints_label.setVisible(False)
        self.stride_entry.clear()
        if hasattr(self, '_default_continuum_stride'):
            del self._default_continuum_stride

    def _toggle_continuum_editor(self):
        # Safety check
        if not hasattr(self, 'plot_viewer'): return

        is_editing = getattr(self.plot_viewer, '_edit_mode_active', False)

        if not is_editing:
            # --- STARTING ---
            if not self.active_history: return

            # --- NEW: Force cursor off before taking background snapshot ---
            if self.cursor_show_checkbox.isChecked():
                self.cursor_show_checkbox.blockSignals(True)
                self.cursor_show_checkbox.setChecked(False)
                self.cursor_show_checkbox.blockSignals(False)
                
                # Erase cursor lines from Matplotlib so they aren't captured
                if hasattr(self, 'plot_viewer'):
                    for artist in self.plot_viewer.canvas.cursor_artists:
                        try: artist.remove()
                        except: pass
                    self.plot_viewer.canvas.cursor_artists = []

            self.edit_hints_label.setVisible(True)

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
            self.cancel_cont_button.setEnabled(True)
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
            self._reset_continuum_ui()

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
                            # The view for the new state, as a switch of
                            # session shows it (F74: this called a method
                            # that does not exist)
                            self._update_view_for_session(new_session, set_current_list_item=True)
                            
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
        default_stride = getattr(self, '_default_continuum_stride', 1500.0)
        
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
