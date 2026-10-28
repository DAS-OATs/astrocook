"""The continuum editor of the main plot: knots spaced in velocity on a draft
curve frozen on the continuum; a moved knot patches the curve around it
(PCHIP between its neighbours); the stride of the knots, reset, save and
cancel.

Moved out of ``pyside_plot.py`` in phase 5 (split of the GUI), as a mixin
of ``SpectrumPlotWidget``, which drives it from its mouse handlers.
"""
import logging
import numpy as np


class ContinuumEditorMixin:
    """The continuum editor of ``SpectrumPlotWidget``.

    Uses of the widget: ``main_window``, ``canvas``, ``toolbar``,
    ``plot_spectrum``, and the state of the editor it initializes
    (``_edit_mode_active``, ``_knots_x``, ``_knots_y``, ``_is_anchor``,
    ``_active_knot_idx``, ``knot_artist``, ``spline_artist``).
    """

    def _generate_velocity_knots(self, x_full, stride_kms):
        """
        Helper: Generates indices for knots spaced by `stride_kms` in velocity.
        """
        # Constants
        c_kms = 299792.458
        
        # 1. Convert wavelength grid to relative velocity (starting at 0)
        # v = c * ln(lambda / lambda_0)
        # We assume x_full is in Angstroms/nm (linear wavelength)
        v_rel = c_kms * np.log(x_full / x_full[0])
        
        # 2. Create the target velocity grid
        # From 0 to max_velocity, step = stride_kms
        target_v = np.arange(0, v_rel[-1], stride_kms)
        
        # 3. Find the indices in the real data that are closest to these velocities
        # searchsorted finds the insertion points to maintain order
        indices = np.searchsorted(v_rel, target_v)
        
        # 4. Ensure the last point is included
        if indices[-1] != len(x_full) - 1:
            indices = np.append(indices, len(x_full) - 1)
            
        # 5. Sanity check: Ensure unique indices (in case stride is smaller than 1 pixel)
        indices = np.unique(indices)
        
        return indices

    def start_continuum_edit(self, initial_stride=1500):
        if not self.main_window.active_history: return
        
        if hasattr(self, 'toolbar') and self.toolbar:
             if self.toolbar.mode == 'pan/zoom': self.toolbar.pan() 
             elif self.toolbar.mode == 'zoom rect': self.toolbar.zoom()

        spec = self.main_window.active_history.current_state.spec
        x_full = spec.x.value
        y_source = spec.cont.value if spec.cont is not None else spec.y.value
            
        indices = self._generate_velocity_knots(x_full, float(initial_stride))
            
        self._knots_x = x_full[indices]
        self._knots_y = y_source[indices]
        self._is_anchor = np.zeros(len(self._knots_x), dtype=bool) 
        
        # --- Set master curve EXACTLY to the existing continuum ---
        self._frozen_spline_y = np.copy(y_source)
        
        self._remove_editor_artists()

        self.knot_artist, = self.canvas.axes.plot(
            self._knots_x, self._knots_y, 
            'o', color='black', markersize=2, alpha=0.5, 
            label='Knots', zorder=10, animated=True
        )
        
        # Notice we plot the frozen array here, not an empty list
        self.spline_artist, = self.canvas.axes.plot(
            x_full, self._frozen_spline_y, 
            '-', color='black', lw=3, alpha=0.3, 
            label='Draft Cont.', zorder=9, animated=True
        )

        self._edit_mode_active = True
        self._draw_cid = self.canvas.mpl_connect('draw_event', self.on_background_draw)
        self.canvas.draw()

    def stop_continuum_edit(self, save=False):
        """Cleans up the editor."""
        self._edit_mode_active = False
        self._active_knot_idx = None
        
        if hasattr(self, '_draw_cid'):
            self.canvas.mpl_disconnect(self._draw_cid)
            del self._draw_cid
        
        self.background = None

        final_x, final_y = None, None
        if save:
            # --- NEW: PASS THE DENSE FROZEN CURVE TO THE RECIPE ---
            # By passing every pixel, the backend PCHIP will exactly recreate your frozen curve.
            spec = self.main_window.active_history.current_state.spec
            final_x = spec.x.value.tolist()
            final_y = self._frozen_spline_y.tolist()

        self._remove_editor_artists()
        self.canvas.draw()
        
        return final_x, final_y
    
    def cancel_continuum_edit(self):
        """
        Forcefully stops the continuum editor and discards all pending visual edits.
        Cleans up artists and background buffers.
        """
        if not self._edit_mode_active:
            return

        logging.info("Continuum edit cancelled.")
        
        # 1. Stop the editor state without saving data
        self._edit_mode_active = False
        self._active_knot_idx = None
        
        # 2. Disconnect the background draw listener
        if hasattr(self, '_draw_cid'):
            self.canvas.mpl_disconnect(self._draw_cid)
            del self._draw_cid
        
        # 3. Wipe visual buffers
        self.canvas.background = None
        self._remove_editor_artists()
        
        # 4. Trigger a clean redraw of the spectrum (using the session's 'cont' column)
        session_state = None
        if self.main_window.active_history:
            session_state = self.main_window.active_history.current_state
        self.plot_spectrum(session_state=session_state)

    def _apply_local_patch(self, center_x):
        """Recomputes the spline ONLY between i-1 and i+1, and pastes it into the frozen curve."""
        x_clean, y_clean = self._get_sanitized_knots()
        spec = self.main_window.active_history.current_state.spec
        x_full = spec.x.value
        
        if len(x_clean) < 4:
            # Fallback for extreme edge cases
            from scipy.interpolate import PchipInterpolator
            cs = PchipInterpolator(x_clean, y_clean)
            self._frozen_spline_y = cs(x_full)
        else:
            # 1. Find the index of the knot closest to the action
            idx = np.argmin(np.abs(x_clean - center_x))
            
            # 2. Define the neighborhood (i-1 to i+1)
            start_idx = max(0, idx - 1)
            end_idx = min(len(x_clean) - 1, idx + 1)
            
            local_x = x_clean[start_idx : end_idx + 1]
            local_y = y_clean[start_idx : end_idx + 1]
            
            # 3. Create a purely local PCHIP
            from scipy.interpolate import PchipInterpolator
            cs_local = PchipInterpolator(local_x, local_y)
            
            # 4. Paste the patch into the master frozen array
            mask = (x_full >= local_x[0]) & (x_full <= local_x[-1])
            self._frozen_spline_y[mask] = cs_local(x_full[mask])
            
        # Update the visual green line
        self.spline_artist.set_data(x_full, self._frozen_spline_y)

    def on_background_draw(self, event):
        """
        Called on every Matplotlib draw (Pan, Zoom, Resize).
        We capture the 'clean' background (Grid+Spectrum) BEFORE painting our artists.
        """
        if not self._edit_mode_active: return
        
        # 1. Capture the CLEAN background now!
        # Since 'animated=True' artists are skipped by the main draw, 
        # the canvas currently contains only the Grid and Spectrum.
        # We save this state immediately.
        bbox = self.canvas.axes.bbox
        self.background = self.canvas.copy_from_bbox(bbox)
        
        # 2. Now paint our animated artists on top
        # This makes them visible to the user, but doesn't 'pollute' 
        # the self.background buffer we just saved.
        if self.knot_artist:
            self.knot_artist.draw(event.renderer)
        if self.spline_artist:
            self.spline_artist.draw(event.renderer)

    def _remove_editor_artists(self):
        """Helper to cleanly remove artists from both Canvas and Matplotlib."""
        # 1. Remove from the Blitting List (Canvas Loop)
        if hasattr(self, 'knot_artist') and self.knot_artist in self.canvas.cursor_artists:
            self.canvas.cursor_artists.remove(self.knot_artist)
        if hasattr(self, 'spline_artist') and self.spline_artist in self.canvas.cursor_artists:
            self.canvas.cursor_artists.remove(self.spline_artist)

        # 2. Remove from the Plot (Matplotlib Axes)
        if hasattr(self, 'knot_artist') and self.knot_artist:
            try: self.knot_artist.remove()
            except: pass # Already removed or not attached
            self.knot_artist = None

        if hasattr(self, 'spline_artist') and self.spline_artist:
            try: self.spline_artist.remove()
            except: pass
            self.spline_artist = None

    def begin_stride_interaction(self):
        """Called when user starts clicking the slider. Caches the current high-fidelity knots."""
        if not self._edit_mode_active: return
        # Save a clean copy of the knots as they were BEFORE resampling started
        self._stride_backup_x, self._stride_backup_y = self._get_sanitized_knots()

    def end_stride_interaction(self):
        """Called when user releases the slider."""
        # Clear the backup so subsequent manual edits are respected
        self._stride_backup_x = None
        self._stride_backup_y = None

    def update_continuum_stride(self, new_stride):
        """
        Resamples visual knots along the purely frozen master curve.
        """
        if not self._edit_mode_active: return
        
        # 1. Keep anchors
        x_live, _, a_live = self._get_sanitized_knots_with_anchors()
        x_anchors = x_live[a_live]
        
        # 2. Generate new floating X grid
        spec = self.main_window.active_history.current_state.spec
        x_full = spec.x.value
        new_indices = self._generate_velocity_knots(x_full, float(new_stride))
        x_new_floating = x_full[new_indices]
        
        # 3. Merge X
        x_merged = np.unique(np.concatenate([x_anchors, x_new_floating]))
        
        # 4. SAMPLE DIRECTLY FROM THE FROZEN CURVE
        y_merged = np.interp(x_merged, x_full, self._frozen_spline_y)
        
        # 5. Restore flags
        self._knots_x = x_merged
        self._knots_y = y_merged
        self._is_anchor = np.isin(x_merged, x_anchors)
        
        self.knot_artist.set_data(self._knots_x, self._knots_y)
        
        # Notice we DO NOT call _apply_local_patch here. The green line stays perfectly frozen.
        self.canvas.draw_idle()

    def _get_sanitized_knots_with_anchors(self):
        if len(self._knots_x) == 0: return np.array([]), np.array([]), np.array([])
        
        sorter = np.argsort(self._knots_x)
        x_s = self._knots_x[sorter]
        y_s = self._knots_y[sorter]
        a_s = self._is_anchor[sorter]
        
        x_u, unique_indices = np.unique(x_s, return_index=True)
        y_u = y_s[unique_indices]
        a_u = a_s[unique_indices]
        
        return x_u, y_u, a_u

    def _get_sanitized_knots(self):
        """Returns sorted, unique copies of knots for spline calculation."""
        x_u, y_u, _ = self._get_sanitized_knots_with_anchors()
        return x_u, y_u
    
    def reset_continuum_to_original(self, stride_kms):
        """
        Discards all manual edits and re-generates the frozen master curve 
        from the underlying session data (spec.cont or spec.y).
        """
        if not self._edit_mode_active or not self.main_window.active_history: 
            return

        logging.info(f"Resetting continuum to original data with stride {stride_kms:.0f} km/s.")

        # 1. Clear any 'dragging' backups
        self._stride_backup_x = None
        self._stride_backup_y = None

        # 2. Retrieve the SOURCE data (Original State)
        spec = self.main_window.active_history.current_state.spec
        x_full = spec.x.value
        if spec.cont is not None:
            y_source = spec.cont.value
        else:
            y_source = spec.y.value
            
        # 3. Generate fresh indices based on the requested stride
        indices = self._generate_velocity_knots(x_full, float(stride_kms))
        
        # 4. Overwrite current knots and reset anchors
        self._knots_x = x_full[indices]
        self._knots_y = y_source[indices]
        self._is_anchor = np.zeros(len(self._knots_x), dtype=bool) 
        
        # --- FIX: Regenerate frozen master curve EXACTLY from source ---
        self._frozen_spline_y = np.copy(y_source)
        
        # 6. Update Visuals
        self.knot_artist.set_data(self._knots_x, self._knots_y)
        self.spline_artist.set_data(x_full, self._frozen_spline_y)
        
        # Force redraw
        self.canvas.draw_idle()

    def refresh_edit_layer(self):
        """
        Refreshes the continuum editing visuals (knots and draft spline).
        Called by MainWindow to ensure the editor layer stays in sync with 
        the background plot during updates or replots.
        """
        if not self._edit_mode_active:
            return

        try:
            # 1. Ensure knots and spline data are up to date with the internal state
            if self.knot_artist:
                self.knot_artist.set_data(self._knots_x, self._knots_y)
            
            if self.spline_artist:
                # The draft spline is already stored in self._frozen_spline_y
                # We ensure the artist matches the full current x-grid
                spec = self.main_window.active_history.current_state.spec
                self.spline_artist.set_data(spec.x.value, self._frozen_spline_y)

            # 2. Request an animation frame update (this uses blitting if background is valid)
            self.canvas.draw_animated()
            
        except Exception as e:
            logging.debug(f"Failed to refresh edit layer: {e}")
