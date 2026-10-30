"""The mouse on the main plot: the redshift cursor, the tooltip after a
pause (component, region of an identification, data near the mouse), the
menu of the right click, the knots of the continuum editor moved, added
and removed with the mouse, and the data inspector opened from the menu.

Moved out of ``pyside_plot.py`` in phase 5 (split of the GUI), as a mixin
of ``SpectrumPlotWidget``; the canvas (``plot_canvas.py``) passes it the
events of matplotlib.
"""
import astropy.units as au
import json
import logging
import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QCursor
from PySide6.QtWidgets import QApplication, QMenu, QToolTip
from typing import Optional, TYPE_CHECKING

from astrocook.core.atomic_data import STANDARD_MULTIPLETS, trans_parse, xem_d

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2

from .data_inspector import DataInspectorDialog
from .main_recipes import region_expression
from .plot_style import SELECTION_STYLE, z_convert_inverse


class PlotMouseMixin:
    """The mouse handlers and the tooltip of ``SpectrumPlotWidget``.

    Uses the widget's ``main_window``, ``canvas``, ``toolbar``,
    ``_tick_visuals`` (the ticks of the components of the last plot), the
    state of the tooltip (``_tip_*``), of the region selection and of the
    continuum editor (``ContinuumEditorMixin``).
    """

    def toggle_region_selector(self):
        """ Enables/disables the region selection mode (the Select button). """
        self.set_region_selection(not self._is_selecting_region)

    def set_region_selection(self, active: bool):
        """Switch the region selection mode on or off (F73).

        On: Pan and Zoom off (the three tools exclude each other), a cross
        cursor; a drag along the x axis selects a region (``on_press``,
        the canvas, ``on_release``)."""
        if active == self._is_selecting_region:
            return
        if active:
            # Before the flag: switching a tool off calls the toolbar's pan
            # or zoom, which would switch the selection off again
            self.toolbar.deactivate_tools()
            self._is_selecting_region = True
            logging.info("Region selection mode ACTIVE.")
            self.setCursor(Qt.CrossCursor)
        else:
            self._is_selecting_region = False
            logging.info("Region selection mode DEACTIVATED.")
            self.setCursor(Qt.ArrowCursor)
            # Clean up any half-drawn selection
            self.canvas.cleanup_selection()
        self.toolbar.sync_select_action()

    def _region_menu(self, xmin: float, xmax: float) -> QMenu:
        """The actions on a selected region (x in nm, the unit of the data):
        each recipe opens its dialog with the region filled in."""
        mw = self.main_window
        expression = region_expression(xmin, xmax)
        menu = QMenu(self)
        menu.addSection(f"Region {xmin:.2f}\u2013{xmax:.2f} nm")
        for text, tip, slot in (
                ("Split region", "A new session with the region only (Edit > Split)",
                 lambda: mw.launch_split_from_region(xmin, xmax)),
                ("Mask region...", "Mask the region in a column (Edit > Mask Expression)",
                 lambda: mw._launch_recipe_dialog(
                     "edit", "mask_expression", initial_params={"expression": expression})),
                ("Rebin region...", "Rebin the region only (Flux > Rebin)",
                 lambda: mw._launch_recipe_dialog(
                     "flux", "rebin", initial_params={"xstart": f"{xmin:.4f}",
                                                      "xend": f"{xmax:.4f}"})),
                ("Zoom to region", "Show the region on the whole plot",
                 lambda: self.canvas.axes.set_xlim(xmin, xmax))):
            action = QAction(text, menu)
            action.setToolTip(tip)
            action.triggered.connect(slot)
            menu.addAction(action)
        menu.setToolTipsVisible(True)
        return menu

    def get_cursor_line_positions_at_z(self, session_state: 'SessionV2', z_cursor):
        """Helper to get theoretical X positions for cursor lines at a specific Z."""
        positions = []
        spec = session_state.spec
        main_window = self.main_window # type: MainWindowV2
        if not spec: return []
        try:
            # V1/V2 Fallback Logic for Transitions
            series_str = main_window.cursor_series_input.text()
            transitions = []
            try:
                # 1. Try V2-native lookup
                transitions = STANDARD_MULTIPLETS[series_str]
            except KeyError:
                # 2. Try the V1 series syntax
                try:
                    transitions = trans_parse(series_str)
                except Exception:
                    pass
                        
                # 3. Fallback for single lines
                if not transitions and series_str in xem_d:
                    transitions = [series_str]


            for t in transitions:
                if t in xem_d:
                    xem_nm = xem_d[t].to_value(au.nm)
                    x_plot_nm = (1 + z_cursor) * xem_nm # This is the data coordinate
                    positions.append(x_plot_nm)
        except (ValueError, AttributeError, ImportError) as e:
            logging.warning(f"Could not calculate cursor positions: {e}")
            return []
        return positions
    
    def calculate_z_from_x(self, mouse_x_nm: float) -> Optional[float]:
        """
        Helper to calculate redshift from a wavelength (nm) based on
        the currently active cursor series.
        """
        if not self.main_window:
            return None
        try:
            # V1/V2 Fallback Logic for Transitions
            series_str = self.main_window.cursor_series_input.text()
            transitions = []
            try:
                # 1. Try V2-native lookup
                transitions = STANDARD_MULTIPLETS[series_str]
            except KeyError:
                # 2. Try the V1 series syntax
                try:
                    transitions = trans_parse(series_str)
                except Exception:
                    pass

                # 3. Fallback for single lines
                if not transitions and series_str in xem_d:
                    transitions = [series_str]

            if not transitions: return None
            
            # Find first valid transition to use as reference
            ref_xem_nm = None
            for t in transitions:
                if t in xem_d:
                    ref_xem_nm = xem_d[t].to_value(au.nm)
                    break
            if ref_xem_nm is None: return None

            return z_convert_inverse(mouse_x_nm, ref_xem_nm)
        except Exception as e:
            logging.error(f"Error calculating Z from X: {e}")
            return None

    def show_data_tooltip(self, event):
        """
        Slow handler called by timer (hover stop).
        Calculates Component, Region, and Data info and displays the unified tooltip.
        """
        # 0. CHECK MODIFIERS: Disable tooltip if Ctrl/Cmd is held
        modifiers = QApplication.keyboardModifiers()
        if (modifiers & Qt.ControlModifier):
            return
        
        if not event or not event.xdata or not event.ydata:
            return
        
        # Don't show data tooltip if redshift cursor is active
        if self.main_window.cursor_show_checkbox.isChecked():
            return
            
        session = None
        if self.main_window and self.main_window.active_history:
            session = self.main_window.active_history.current_state
        if not session or not session.spec:
            return
            
        spec = session.spec
        x_full = spec.x.value

        # --- 1. Detect Component ---
        self._tip_component_html = None
        found_comp = None
        
        if event.inaxes:
            xlim = event.inaxes.get_xlim()
            tol = (xlim[1] - xlim[0]) * 0.01 
            
            for ax, x_pos, comp in self._tick_visuals:
                if ax == event.inaxes and abs(event.xdata - x_pos) < tol:
                    found_comp = comp
                    break
        
        if found_comp:
            c_chi2 = f"{found_comp.chi2:.2f}" if found_comp.chi2 is not None else "N/A"
            c_resol = f"{found_comp.resol:.0f}" if found_comp.resol is not None else "spectrum"
            
            self._tip_component_html = (
                f"<div style='font-size: 13px'>"
                f"<b>ID:</b> {found_comp.id} | <b>{found_comp.series}</b><br>"
                f"<b>z:</b> {found_comp.z:.5f}<br>"
                f"<b>\u03c7\u00b2:</b> {c_chi2} | <b>R:</b> {c_resol}"
                f"</div>"
            )

        # --- 2. Find Nearest Data Index ---
        idx = np.searchsorted(x_full, event.xdata)
        if idx > 0 and (idx == len(x_full) or 
                        np.abs(event.xdata - x_full[idx-1]) < np.abs(event.xdata - x_full[idx])):
            idx -= 1
        
        if idx < 0 or idx >= len(x_full):
            if self._tip_component_html: self._update_tooltip_display()
            return
        
        # --- 3. Region Identification ---
        self._tip_region_html = None
        if 'abs_ids' in spec._data.aux_cols:
            try:
                col_obj = spec._data.aux_cols['abs_ids']
                ids_array = col_obj.values if hasattr(col_obj, 'values') else col_obj
                val = ids_array[idx]
                
                # --- We only print if it looks like a valid region to avoid spamming the console on empty space ---
                is_valid_pixel = not (np.ma.is_masked(val) or (isinstance(val, (float, np.floating)) and np.isnan(val)) or val == 0)
                
                if is_valid_pixel:
                    region_id_val = int(val)                    
                    ident_labels = spec._data.meta.get('region_identifications')
                    
                    if ident_labels:
                        while isinstance(ident_labels, str):
                            try:
                                ident_labels = json.loads(ident_labels)
                            except json.JSONDecodeError:
                                break 
                        
                        if isinstance(ident_labels, dict):
                            
                            # Try both lookups
                            res_int = ident_labels.get(region_id_val)
                            res_str = ident_labels.get(str(region_id_val))
                            region_data = res_int or res_str
                            
                            if region_data:
                                if isinstance(region_data, list):
                                    labels = []
                                    for entry in region_data:
                                        if isinstance(entry, list) and len(entry) >= 2:
                                            name, score = entry[0], entry[1]
                                            labels.append(f"{name} ({score:.2f})" if score > 0 else name)
                                        else:
                                            labels.append(str(entry))
                                            
                                    label_str = ", ".join(labels)
                                    if label_str:
                                        self._tip_region_html = f"<tr><td style='padding-right:8px'><b>Likely:</b></td><td>{label_str}</td></tr>"
                else:
                    region_id_val = 0

            except Exception:
                pass

        # --- 4. Data Points (Zoom Check) ---
        self._tip_data_rows = []
        
        xlim = self.canvas.axes.get_xlim()
        idx_start_vis = np.searchsorted(x_full, xlim[0], side='left')
        idx_end_vis = np.searchsorted(x_full, xlim[1], side='right')
        n_visible = idx_end_vis - idx_start_vis
        TOOLTIP_MAX_POINTS = 20000 
        
        if n_visible <= TOOLTIP_MAX_POINTS:
            # Determine Normalization State
            is_norm = self.main_window.norm_y_checkbox.isChecked()
            sym_flux = "\u0192" if is_norm else "F"
            sym_dflux = "d" + sym_flux
            
            # [NEW] Calculate dLambda (Pixel Size)
            d_lambda_val = None
            if spec.xmin is not None and spec.xmax is not None:
                try:
                    # Access values safely assuming they are Columns/Quantities
                    v_xmin = spec.xmin.value[idx] if hasattr(spec.xmin, 'value') else spec.xmin[idx]
                    v_xmax = spec.xmax.value[idx] if hasattr(spec.xmax, 'value') else spec.xmax[idx]
                    d_lambda_val = v_xmax - v_xmin
                except Exception:
                    pass

            # Pre-calculate Continuum Value
            norm_val_at_idx = 1.0
            has_norm = False
            # [CHANGE] Access .norm property
            norm_q = spec.norm
            if norm_q is not None:
                norm_val_at_idx = norm_q.value[idx]
                has_norm = True
                
            if norm_val_at_idx == 0: norm_val_at_idx = np.nan

            # Helper: Calculate PLOTTED value and DISPLAY value
            def get_vals(col_type):
                if col_type == 'y':
                    raw = spec.y.value[idx]
                    if is_norm: return (raw / norm_val_at_idx), (raw / norm_val_at_idx), ""
                    return raw, raw, str(spec.y.unit)
                
                elif col_type == 'dy':
                    raw = spec.dy.value[idx]
                    if is_norm: return (raw / norm_val_at_idx), (raw / norm_val_at_idx), ""
                    return raw, raw, str(spec.dy.unit)

                elif col_type == 'cont':
                    if not has_norm: return None, None, ""
                    if is_norm: return 1.0, 1.0, ""
                    return norm_val_at_idx, norm_val_at_idx, str(spec.y.unit)

                elif col_type == 'model':
                    if spec.model is None: return None, None, ""
                    raw = spec.model.value[idx]
                    if is_norm: return (raw / norm_val_at_idx), (raw / norm_val_at_idx), ""
                    return raw, raw, str(spec.y.unit)
                return None, None, ""

            # Check Proximity
            THRESHOLD_PIX = 40
            def check_dist(y_val_plotted):
                if y_val_plotted is None or not np.isfinite(y_val_plotted): return False
                data_pos = self.canvas.axes.transData.transform((x_full[idx], y_val_plotted))
                mouse_pos = (event.x, event.y)
                return np.hypot(data_pos[0] - mouse_pos[0], data_pos[1] - mouse_pos[1]) < THRESHOLD_PIX

            # Helper to add Lambda and dLambda rows
            def _add_lambda_rows():
                self._tip_data_rows.append(("\u03bb", f"{x_full[idx]:.4f}", f"{spec.x.unit}"))
                if d_lambda_val is not None:
                    self._tip_data_rows.append(("d\u03bb", f"{d_lambda_val:.4f}", f"{spec.x.unit}"))

            # A. Check Flux
            y_plot, y_disp, y_unit = get_vals('y')
            if check_dist(y_plot):
                _add_lambda_rows()
                self._tip_data_rows.append((sym_flux, f"{y_disp:.4e}", y_unit))
                
                # Add error
                _, dy_disp, _ = get_vals('dy')
                self._tip_data_rows.append((sym_dflux, f"{dy_disp:.4e}", y_unit))

            # B. Check Continuum
            if self.main_window.continuum_checkbox.isChecked():
                c_plot, c_disp, c_unit = get_vals('cont')
                if check_dist(c_plot):
                    # Add lambda if not already added
                    if not self._tip_data_rows:
                        _add_lambda_rows()
                    self._tip_data_rows.append(("Cont", f"{c_disp:.4e}", c_unit))

            # C. Check Model
            if self.main_window.model_checkbox.isChecked():
                m_plot, m_disp, m_unit = get_vals('model')
                if check_dist(m_plot):
                    if not self._tip_data_rows:
                        _add_lambda_rows()
                    self._tip_data_rows.append(("Mod", f"{m_disp:.4e}", m_unit))

        # 5. Refresh Display
        self._update_tooltip_display()

    def _update_tooltip_display(self):
        """
        Combines component, region, and data info into a single HTML tooltip.
        Order: Data Points -> Region Info -> Component Info
        """
        html_parts = []
        has_previous_content = False

        # --- 1. Data Points (Top) ---
        if self._tip_data_rows:
            rows_html = ""
            for label, val_str, unit_str in self._tip_data_rows:
                # Format: <b>Symbol:</b> Value Unit
                rows_html += f"<tr><td style='padding-right:8px'><b>{label}:</b></td><td>{val_str} {unit_str}</td></tr>"
            
            html_parts.append(f"<table style='font-size: 13px'>{rows_html}</table>")
            has_previous_content = True

        # --- 2. Region Info (Middle) ---
        if self._tip_region_html:
            # Add separator if Data Points were shown
            if has_previous_content:
                html_parts.append("<hr style='margin: 4px 0px;'>")
            
            # _tip_region_html is a pre-formatted <tr>...</tr> string. 
            # We wrap it in a fresh table to ensure alignment and structure.
            html_parts.append(f"<table style='font-size: 13px'>{self._tip_region_html}</table>")
            has_previous_content = True

        # --- 3. Component Info (Bottom) ---
        if self._tip_component_html:
            # Add separator if Data or Region info was shown
            if has_previous_content:
                html_parts.append("<hr style='margin: 4px 0px;'>")
            
            # _tip_component_html is a pre-formatted <div>
            html_parts.append(self._tip_component_html)

        # --- Finalize ---
        if not html_parts:
            QToolTip.hideText()
            return

        content = "".join(html_parts)
        # The CSS in __init__ now handles the border and rounding.
        QToolTip.showText(QCursor.pos(), content, self)

    def on_mouse_move(self, event):
        if not event.inaxes: return

        # Gate: Ignore if toolbar is active
        if hasattr(self, 'toolbar') and self.toolbar.mode != '':
            return

        # --- CONTINUUM EDITING LOGIC ---
        if self._edit_mode_active and self._active_knot_idx is not None:
            new_x = event.xdata
            new_y = event.ydata
            
            self._knots_x[self._active_knot_idx] = new_x
            self._knots_y[self._active_knot_idx] = new_y
            self._is_anchor[self._active_knot_idx] = True # NEW: Mark moved knots as anchors
            
            # 2. Update the Artist Objects
            self.knot_artist.set_data(self._knots_x, self._knots_y)
            self._apply_local_patch(new_x)
            
            # 3. BLITTING (The "Live" Update)
            if self.background:
                self.canvas.restore_region(self.background)
                self.canvas.axes.draw_artist(self.knot_artist)
                self.canvas.axes.draw_artist(self.spline_artist)
                self.canvas.blit(self.canvas.axes.bbox)
                
            return # Stop processing
            
    def on_press(self, event):
        if not event.inaxes: return
        
        # --- 0. CHECK MODIFIERS ---
        modifiers = QApplication.keyboardModifiers()
        is_ctrl_held = bool(modifiers & Qt.ControlModifier)
        is_shift_held = bool(modifiers & Qt.ShiftModifier)

        # Check if a toolbar mode (Zoom/Pan) is currently active
        is_tool_active = bool(self.toolbar.mode)

        # --- 1. CONTINUUM EDITING LOGIC (Highest Priority) ---
        if self._edit_mode_active and event.inaxes == self.canvas.axes:
            click_x, click_y = event.xdata, event.ydata
            
            # Distance check for existing knots
            knot_pixels = self.canvas.axes.transData.transform(np.column_stack([self._knots_x, self._knots_y]))
            click_pixel = self.canvas.axes.transData.transform((click_x, click_y))
            
            nearest_idx = None
            min_dist = float('inf')
            if len(knot_pixels) > 0:
                dists = np.hypot(knot_pixels[:,0] - click_pixel[0], knot_pixels[:,1] - click_pixel[1])
                nearest_idx = np.argmin(dists)
                min_dist = dists[nearest_idx]
                
            THRESHOLD = 10 # pixels
            
            # A. SHIFT + LEFT CLICK: Add Knot
            if event.button == 1 and is_shift_held:
                spec = self.main_window.active_history.current_state.spec
                # NEW: Snap the new knot exactly to the frozen curve
                new_y = np.interp(click_x, spec.x.value, self._frozen_spline_y) 
                    
                self._knots_x = np.append(self._knots_x, click_x)
                self._knots_y = np.append(self._knots_y, new_y)
                self._is_anchor = np.append(self._is_anchor, True) 
                
                sorter = np.argsort(self._knots_x)
                self._knots_x = self._knots_x[sorter]
                self._knots_y = self._knots_y[sorter]
                self._is_anchor = self._is_anchor[sorter]
                
                self.knot_artist.set_data(self._knots_x, self._knots_y)
                # NEW: Patch locally around the new knot
                self._apply_local_patch(click_x) 
                
                self.canvas.draw_idle()
                return True 

            # B. RIGHT CLICK (on a knot): Remove Knot
            elif event.button == 3 and min_dist < THRESHOLD:
                if len(self._knots_x) > 2:
                    removed_x = self._knots_x[nearest_idx] # Remember where it was
                    self._knots_x = np.delete(self._knots_x, nearest_idx)
                    self._knots_y = np.delete(self._knots_y, nearest_idx)
                    self._is_anchor = np.delete(self._is_anchor, nearest_idx)
                    
                    self.knot_artist.set_data(self._knots_x, self._knots_y)
                    # NEW: Patch locally across the gap we just created
                    self._apply_local_patch(removed_x) 
                    
                    self.canvas.draw_idle()
                return True
                
            # C. LEFT CLICK (on a knot, NO shift): Grab Knot to move
            elif event.button == 1 and min_dist < THRESHOLD and not is_shift_held:
                self._active_knot_idx = nearest_idx
                # Pause toolbar if active
                if hasattr(self, 'toolbar') and self.toolbar.mode != '':
                    self._paused_tool_mode = self.toolbar.mode
                    if self.toolbar.mode == 'pan/zoom': self.toolbar.pan()
                    elif self.toolbar.mode == 'zoom rect': self.toolbar.zoom()
                return True

        # --- 2. EXISTING RIGHT CLICK CONTEXT MENU ---
        # Condition: Right Click (Button 3) AND (Ctrl held OR No tool active)
        if event.button == 3 and (is_ctrl_held or not is_tool_active):
            
            # A. Pause Toolbar (if active) to prevent artifacts
            paused_mode = None
            if hasattr(self, 'toolbar') and self.toolbar.mode != '':
                paused_mode = self.toolbar.mode
                if self.toolbar.mode == 'pan/zoom': 
                    self.toolbar.pan() # Toggle off
                elif self.toolbar.mode == 'zoom rect': 
                    self.toolbar.zoom() # Toggle off

            # B. Prepare Menu
            found_comp = None
            xlim = event.inaxes.get_xlim(); tol = (xlim[1] - xlim[0]) * 0.015 
            for ax, x_pos, comp in self._tick_visuals:
                if ax == event.inaxes and abs(event.xdata - x_pos) < tol: found_comp = comp; break
            
            menu = QMenu(self)
            
            # --- Menu Items ---
            if found_comp:
                act_inspect = QAction(f"Inspect {found_comp.series} (its group)", menu)
                if self.main_window: act_inspect.triggered.connect(lambda checked=False, u=found_comp.uuid: self.main_window.open_inspector_on_component(u))
                menu.addAction(act_inspect)
                menu.addSeparator()

            z_click = self.calculate_z_from_x(event.xdata)
            if z_click is not None and self.main_window and self.main_window.cursor_show_checkbox.isChecked():
                series_str = self.main_window.cursor_series_input.text()
                act_add = QAction(f"Add {series_str} at z={z_click:.5f}", menu)
                act_add.triggered.connect(lambda checked=False, z=z_click, s=series_str: self.main_window._on_recipe_requested("absorbers", "add_component", {'series': s, 'z': z}, {}))
                menu.addAction(act_add)
                
                act_zem = QAction(f"Set Emission Redshift to z={z_click:.5f}", menu)
                act_zem.triggered.connect(lambda checked=False, z=z_click: self.main_window._on_recipe_requested("edit", "set_properties", {"z_em": str(z)}, {}))
                menu.addAction(act_zem)

            # --- NEW: Split Action ---
            menu.addSeparator()
            act_split = QAction("Split from Current Zoom", menu)
            act_split.setToolTip("Create a new session containing only the visible spectral region.")
            if self.main_window:
                act_split.triggered.connect(lambda: self.main_window.launch_split_from_current_view())
            menu.addAction(act_split)

            # --- [NEW] Inspect Action ---
            menu.addSeparator()
            
            # Find nearest index for the label
            if self.main_window and self.main_window.active_history:
                session = self.main_window.active_history.current_state
                if session and session.spec:
                    x_full = session.spec.x.value
                    idx_click = np.searchsorted(x_full, event.xdata)
                    idx_click = np.clip(idx_click, 0, len(x_full)-1)
                    
                    x_val_click = x_full[idx_click]
                    x_unit = str(session.spec.x.unit)
                    
                    act_inspect_data = QAction(f"Inspect Data (near {x_val_click:.4f} {x_unit})", menu)
                    act_inspect_data.triggered.connect(
                        lambda checked=False, s=session, i=idx_click: self._open_data_inspector(s, i)
                    )
                    menu.addAction(act_inspect_data)
            # ---------------------------

            # C. Execute Menu
            if not menu.isEmpty(): 
                menu.exec(QCursor.pos())

            # D. Restore Toolbar immediately
            if paused_mode:
                if paused_mode == 'pan/zoom': self.toolbar.pan()
                elif paused_mode == 'zoom rect': self.toolbar.zoom()

            return True # Consume event

        # --- 3. Region Selection (Button 1) ---
        if event.button == 1 and event.inaxes == self.canvas.axes and self._is_selecting_region:
            self.canvas.cleanup_selection()
            self.canvas.selection_start_x = event.xdata
            self.canvas.selection_artist = self.canvas.axes.axvspan(
                event.xdata, event.xdata, **SELECTION_STYLE)
            self.canvas.draw_idle() 
            return True

        return False
    
    def on_release(self, event):
        # 1. Drop the knot
        if self._active_knot_idx is not None:
            self._active_knot_idx = None
            
        # --- [START NEW BLOCK] ---
        # 2. Restore Toolbar if we paused it
        if hasattr(self, '_paused_tool_mode') and self._paused_tool_mode:
            if self._paused_tool_mode == 'pan/zoom': 
                self.toolbar.pan()
            elif self._paused_tool_mode == 'zoom rect': 
                self.toolbar.zoom()
            self._paused_tool_mode = None

        # 3. A region selected: its span stays while the menu of the
        #    actions is open; the menu closed (an action or Esc), it goes
        if self._is_selecting_region and self.canvas.selection_start_x is not None:
            start = self.canvas.selection_start_x
            end = event.xdata if event.inaxes == self.canvas.axes else None
            self.canvas.selection_start_x = None
            if start is None or end is None or start == end:
                self.canvas.cleanup_selection()
                return
            xmin, xmax = sorted([start, end])
            if self.canvas.selection_artist is not None:
                self.canvas.selection_artist.remove()
            self.canvas.selection_artist = self.canvas.axes.axvspan(
                xmin, xmax, **SELECTION_STYLE)
            self.canvas.draw_idle()
            logging.info(f"Region selected: {xmin:.4f} to {xmax:.4f}")
            if self.main_window:
                self._region_menu(xmin, xmax).exec(QCursor.pos())
            self.canvas.cleanup_selection()

    def _open_data_inspector(self, session, center_idx):
        """Launches the Data Inspector Dialog (Non-blocking)."""
        # 1. Close existing dialog if it exists
        # This prevents opening 50 windows if the user keeps clicking
        if hasattr(self, '_data_inspector_dialog') and self._data_inspector_dialog is not None:
            try:
                self._data_inspector_dialog.close()
            except Exception:
                pass
        
        # 2. Create and Store Reference
        # We save it to 'self' so Python doesn't garbage-collect it immediately
        self._data_inspector_dialog = DataInspectorDialog(session, center_idx, window=15, parent=self)
        
        # 3. Clean up reference when closed
        # This ensures we don't hold onto memory unnecessarily
        self._data_inspector_dialog.finished.connect(
            lambda: setattr(self, '_data_inspector_dialog', None)
        )
        
        # 4. Show non-modally
        self._data_inspector_dialog.show()
