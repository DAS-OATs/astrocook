# astrocook/v2/gui/pyside_plot.py

from matplotlib import pylab
import astropy.units as au
import logging
import matplotlib
matplotlib.use('QtAgg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.style as mplstyle
import numpy as np
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QApplication, QWidget, QVBoxLayout
import scienceplots
from typing import Optional, TYPE_CHECKING

# Use TYPE_CHECKING to avoid circular import errors at runtime
if TYPE_CHECKING:
    from ..session import SessionV2 

try:
    from ...v1.functions import trans_parse, x_convert # Need x_convert for systs
    from ...v1.vars import xem_d
    V1_FUNCTIONS_AVAILABLE = True
except ImportError:
    V1_FUNCTIONS_AVAILABLE = False
    logging.error("V1 functions/vars not found, cursor/system plotting may fail.")
    
# The x_plot it receives is *always* in the data's native units (nm)
def z_convert_inverse(x_plot_nm, xem_nm): # <<< Removed zem_spec, x_unit
    """Converts plot X coordinate (in data units, nm) back to redshift."""
    if not V1_FUNCTIONS_AVAILABLE: return None
    try:
        z = (x_plot_nm / xem_nm) - 1.0 # Simple calculation
        return z
    except Exception as e:
        logging.debug(f"Inverse z conversion failed: {e}")
        return None

# --- Helper to get colors ---
def get_color_cycle(n=10, fallback_cmap='viridis'):
    # ... (function definition as before) ...
    try:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        return [colors[i % len(colors)] for i in range(n)]
    except Exception as e:
        logging.warning(f"Could not get color cycle: {e}. Using fallback.")
        cmap = plt.get_cmap(fallback_cmap); return [cmap(i / n) for i in range(n)]
# ----------------------------


class MatplotlibCanvas(FigureCanvasQTAgg):
    """A custom Matplotlib canvas enabling blitting for cursor dragging."""
    def __init__(self, parent=None, width=5, height=4, dpi=100, plot_widget=None): # Added plot_widget ref
        try:
            plt.style.use(['science', 'fast'])
            #plt.style.use('fast')
        except Exception as e:
            logging.warning(f"Could not apply 'scienceplots': {e}. Using default.")
            plt.style.use('fast')

        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        fig.tight_layout()

        super(MatplotlibCanvas, self).__init__(fig)

        # We rely on the QtAgg backend's native speed for now.
        if parent:
            self.setParent(parent)

        self.plot_widget = plot_widget # Store reference

        # Blitting Attributes
        self.background = None # Store the clean background pixels
        self.cursor_artists = [] # Store references to the cursor axvline objects
        self.draw_event_cid = None # Store connection ID for draw_event

        # Connect Matplotlib events
        self.mpl_connect('motion_notify_event', self.on_motion)
        self.xlim_cid = self.axes.callbacks.connect('xlim_changed', self.on_lim_changed)
        self.ylim_cid = self.axes.callbacks.connect('ylim_changed', self.on_lim_changed)


    def on_lim_changed(self, axes):
        """Callback: Invalidate background and schedule a full redraw."""
        # Only act if a background existed (prevents initial draw loops)
        if self.background is not None:
            logging.debug(f"Limits changed ({axes.get_label()}): Invalidating background & scheduling plot_spectrum.")
            self.background = None # Invalidate background cache

            # Schedule plot_spectrum to run soon via the event loop
            if self.plot_widget:
                QTimer.singleShot(0, self.plot_widget.plot_spectrum) # <<< Schedule call

    def _capture_background(self, event):
        """Callback for draw_event to capture the background for blitting."""
        capture_successful = False
        try:
            if self.background is None and self.figure.get_size_inches().any():
                if self.axes.bbox.width > 0 and self.axes.bbox.height > 0:
                    self.background = self.copy_from_bbox(self.axes.bbox)
                    capture_successful = True
        except Exception as e:
            logging.error(f"Failed during background copy: {e}")
            self.background = None
        finally:
            if self.draw_event_cid is not None:
                try: self.mpl_disconnect(self.draw_event_cid)
                except Exception: pass
                self.draw_event_cid = None
                 
    def draw_animated(self):
        if self.background is None:
            # This can happen in the small gap between on_lim_changed
            # and the full redraw completing. Calling draw_idle()
            # is a safe, lightweight way to request the redraw.
            self.draw_idle()
            return

        # Proceed with blitting
        self.restore_region(self.background)
        drawn_artists = []
        for artist in self.cursor_artists:
            try: self.axes.draw_artist(artist); drawn_artists.append(artist)
            except Exception as e: logging.error(f"Error drawing anim artist {artist}: {e}")
        if drawn_artists:
            try: self.blit(self.axes.bbox)
            except Exception as e: logging.error(f"Blitting failed: {e}"); self.draw_idle()


    def on_press(self, event):
        """Initiates cursor drag OR handles other clicks."""
        # ... (Safety check for plot_widget/main_window) ...
        if not self.plot_widget or not hasattr(self.plot_widget, 'main_window') or not self.plot_widget.main_window: return

        main_window = self.plot_widget.main_window

        if event.button == 3 and event.inaxes == self.axes:
             logging.debug("Right-click detected.")
             # Add context menu logic here if needed later


    def on_motion(self, event):
        """Handles mouse motion: Updates cursor position if active."""
        # Safety checks
        if (not self.plot_widget or not hasattr(self.plot_widget, 'main_window')
            or not self.plot_widget.main_window):
            # logging.debug("on_motion: plot_widget or main_window invalid.")
            return

        main_window = self.plot_widget.main_window # type: MainWindowV2
        
        # Check if inside axes and cursor checkbox is checked
        if (event.inaxes == self.axes and event.xdata is not None
            and main_window.cursor_show_checkbox.isChecked()):

            # --- Always update on motion if cursor is active ---
            mouse_x = event.xdata
            spec = self.plot_widget.session_state.spec
            if not spec: return

            try:
                # Calculate new_z based on mouse_x
                series_str = main_window.cursor_series_input.text(); transitions = trans_parse(series_str)
                if not transitions or not V1_FUNCTIONS_AVAILABLE: return
                ref_transition = None; ref_xem_nm = None
                for t in transitions:
                    if t in xem_d: ref_transition = t; ref_xem_nm = xem_d[t].to_value(au.nm); break
                if ref_xem_nm is None: return
                zem_spec = getattr(spec, '_zem', 0.0); x_unit = spec.x.unit
                new_z = z_convert_inverse(mouse_x, ref_xem_nm)

                if new_z is not None:
                    # Update QLineEdit (no need for blockSignals if not dragging)
                    current_text_z = None
                    try: current_text_z = float(main_window.cursor_z_input.text())
                    except ValueError: pass
                    # Update text only if significantly different to avoid excessive signals
                    if current_text_z is None or not np.isclose(new_z, current_text_z, atol=1e-7):
                        main_window.cursor_z_input.setText(f"{new_z:.7f}")

                    # Update Artist Positions & Redraw using Blit
                    new_x_positions = self.plot_widget.get_cursor_line_positions_at_z(new_z)

                    # Ensure artists exist and match count
                    if not self.cursor_artists or len(new_x_positions) != len(self.cursor_artists):
                        # Force a full redraw if artists are missing/mismatched
                        logging.warning("Cursor artists missing or mismatched on motion. Forcing full redraw.")
                        self.plot_widget.plot_spectrum() # Recreate artists
                    else:
                        # Update existing artists
                        for artist, x_pos in zip(self.cursor_artists, new_x_positions):
                            artist.set_xdata([x_pos])
                        self.draw_animated() # Use blitting redraw

            except Exception as e:
                logging.error(f"Error during cursor hover update: {e}", exc_info=True)
        # else:
            # Optional: Clear coordinate display if outside axes or cursor off?
            # Toolbar should handle clearing coordinate display by default.
            # pass

    def on_release(self, event):
        """Handles mouse button release events to stop cursor drag."""
        # ... (Safety checks) ...
        if not self.plot_widget: return

        pass

class SpectrumPlotWidget(QWidget):
    """
    Refactored widget that contains the plot canvas, toolbar, and toggles.
    This replaces the content that was previously inside SpectrumViewerPySide.
    """
    def __init__(self, initial_session_state: Optional['SessionV2'], main_window_ref):
        super().__init__()
        self.main_window = main_window_ref
        
        # Main Vertical Layout for the plot area
        self.main_layout = QVBoxLayout(self)

        # 1. Canvas
        self.canvas = MatplotlibCanvas(parent=self, plot_widget=self)
        self.main_layout.addWidget(self.canvas, 1)

        # 2. Toolbar (NO Checkboxes here anymore)
        self.toolbar = AstrocookToolbar(self.canvas, self) # Pass self as parent

        # Add ONLY toolbar to main layout
        self.main_layout.addWidget(self.toolbar, 0) # Stretch 0

        # Draw initial plot (will read checkbox state from main window)
        self.plot_spectrum(session_state=initial_session_state, initial_draw=True)

    # --- Refactored Plotting Method (Requires an update) ---
    def plot_spectrum(self, session_state: Optional['SessionV2'], initial_draw=False):
        """
        Retrieves data from the immutable V2 Session and plots it.
        """
        if session_state is None:
            logging.debug("plot_spectrum called with no session manager. Clearing plot.")
            ax = self.canvas.axes
            ax.clear()
            ax.set_title("No Spectrum Data Loaded")
            # Ensure background capture logic is handled if needed after clear
            self.canvas.background = None # Reset background
            # Connect draw_event for potential background capture if needed later
            if self.canvas.draw_event_cid is not None:
                try: self.canvas.mpl_disconnect(self.canvas.draw_event_cid)
                except Exception: pass
            self.canvas.draw_event_cid = self.canvas.mpl_connect('draw_event', self.canvas._capture_background)
            # Trigger draw
            if initial_draw: self.canvas.draw()
            else: self.canvas.draw_idle()
            return # Exit early

        spec = session_state.spec
        systs = session_state.systs # Get systs object

        ax = self.canvas.axes

        # ** Store current limits BEFORE clearing **
        # ** Make sure to store current xlim/ylim BEFORE ax.clear() **
        previous_xlim = ax.get_xlim()
        previous_ylim = ax.get_ylim()
        was_zoomed = False
        try: 
            if previous_xlim != (0.0, 1.0) or previous_ylim != (0.0, 1.0):
                if all(isinstance(v,(int,float)) for v in previous_xlim+previous_ylim) and \
                    previous_xlim[0] != previous_xlim[1] and previous_ylim[0] != previous_ylim[1]:
                    was_zoomed = True
        except Exception: pass # Ignore if axes not ready
        ax.clear() # Now clear axes

        # ** Reconnect limit change callbacks here **
        # Disconnect previous first, just in case
        if self.canvas.xlim_cid is not None:
            try: self.canvas.axes.callbacks.disconnect(self.canvas.xlim_cid)
            except Exception: pass
        if self.canvas.ylim_cid is not None:
            try: self.canvas.axes.callbacks.disconnect(self.canvas.ylim_cid)
            except Exception: pass
        # Reconnect
        try:
            self.canvas.xlim_cid = self.canvas.axes.callbacks.connect('xlim_changed', self.canvas.on_lim_changed)
            self.canvas.ylim_cid = self.canvas.axes.callbacks.connect('ylim_changed', self.canvas.on_lim_changed)
            logging.debug("Reconnected limit change callbacks in plot_spectrum.")
        except Exception as e:
            logging.error(f"Failed to reconnect limit callbacks: {e}")
        # ----------------------------------------

        # ** Determine Legend Location **
        legend_loc = 'best' # Default
        if not initial_draw:
            try:
                # Check width as a proxy for visibility/expanded state
                left_sidebar_open = self.main_window.left_sidebar_widget.width() > 1
                right_sidebar_open = self.main_window.right_sidebar_widget.width() > 1

                if left_sidebar_open and not right_sidebar_open:
                    legend_loc = 'upper right'
                elif not left_sidebar_open and right_sidebar_open:
                    legend_loc = 'upper left'
                elif left_sidebar_open and right_sidebar_open:
                    # Maybe center top if both are open? Or stick to one side.
                    legend_loc = 'upper center'
                # else: use 'best' if both closed
                logging.debug(f"Legend location set to: {legend_loc} (L:{left_sidebar_open}, R:{right_sidebar_open})")
            except Exception as e:
                logging.warning(f"Could not determine sidebar state for legend: {e}")
                legend_loc = 'best' # Fallback

        # ** Reset background and clear artists for blitting **
        self.canvas.background = None
        self.canvas.cursor_artists = []
        
        plot_occurred = False
        if spec and len(spec.x) > 0:
            plot_occurred = True

            # --- ** Data Slicing Logic ** ---
            # Get full data arrays
            full_x_data = spec.x.value
            full_y_data = spec.y.value
            full_dy_data = spec.dy.value if spec.dy is not None else np.zeros_like(full_y_data)
            full_cont_data = spec.cont   # Returns np.ndarray or None
            full_model_data = spec.model # Returns np.ndarray or None
            
            x_unit_str = str(spec.x.unit); y_unit_str = str(spec.y.unit)

            # --- Check View Toggles ---
            is_norm_y = self.main_window.norm_y_checkbox.isChecked()
            is_log_x = self.main_window.log_x_checkbox.isChecked()
            is_log_y = self.main_window.log_y_checkbox.isChecked()
            selected_x_unit = self.main_window.x_unit_combo.currentText()

            # --- ** Get CURRENT target Axes Limits ** ---
            # Get limits AFTER ax.clear() but BEFORE plotting.
            # These might be autoscaled (inf/-inf initially) or set by zoom/pan.
            target_xlim = ax.get_xlim()
            target_ylim = ax.get_ylim()
            logging.debug(f"Target limits for slicing: X={target_xlim}, Y={target_ylim}")
            # ------------------------------------------

            data_slice = slice(None) # Default: use full array
            

            # If zoomed in, calculate the slice
            if was_zoomed and not initial_draw:
                try:
                    # Find indices for the visible range
                    # Use searchsorted for fast lookup on sorted x data
                    idx_start = np.searchsorted(full_x_data, previous_xlim[0], side='left')
                    idx_end = np.searchsorted(full_x_data, previous_xlim[1], side='right')

                    # Add a small buffer (e.g., 50 points) to each side for cleaner edges
                    buffer = 50 
                    idx_start = max(0, idx_start - buffer)
                    idx_end = min(len(full_x_data), idx_end + buffer)

                    if idx_end > idx_start: # Ensure slice is valid
                        data_slice = slice(idx_start, idx_end)
                        logging.debug(f"Plotting sliced data: {idx_end - idx_start} points (was {len(full_x_data)})")
                    else:
                        logging.debug("Zoom slice is empty, plotting full data.")
                except Exception as e:
                    logging.warning(f"Failed to calculate plot slice: {e}")

            # Apply the slice
            x_data = full_x_data[data_slice]
            y_data = full_y_data[data_slice]
            dy_data = full_dy_data[data_slice]
            cont_data = full_cont_data[data_slice] if full_cont_data is not None else None
            model_data = full_model_data[data_slice] if full_model_data is not None else None

            if is_norm_y and cont_data is not None:
                # Avoid division by zero
                y_data = np.divide(y_data, cont_data, out=np.full_like(y_data, np.nan), where=cont_data!=0)
                dy_data = np.divide(dy_data, cont_data, out=np.full_like(dy_data, np.nan), where=cont_data!=0)
                if model_data is not None:
                    model_data = np.divide(model_data, cont_data, out=np.full_like(model_data, np.nan), where=cont_data!=0)
                cont_data = np.ones_like(cont_data) # Normalized continuum is 1

            x_unit = str(spec.x.unit)
            y_unit = str(spec.y.unit)
                    
            colors = get_color_cycle(5)

            # 1. Plot Main Flux (Uses Matplotlib style defaults for color)
            ax.step(x_data, y_data, where='mid', label="Spectrum", lw=0.5, color=colors[0], rasterized=True)

            # --- Check state from main_window reference ---
            # 2. Plot Error Shading (Conditional)
            if self.main_window.error_checkbox.isChecked(): # <<< Check main window's checkbox
                ax.fill_between(
                    x_data, y_data - dy_data, y_data + dy_data,
                    step='mid', color='#aaaaaa', alpha=0.5,
                    label='1-sigma error', rasterized=True
                )

            # 3. Plot Continuum (Conditional)
            if self.main_window.continuum_checkbox.isChecked(): # <<< Check main window's checkbox
                if cont_data is not None:
                    ax.plot(
                        x_data, cont_data.value, linestyle='--',
                        color='black', lw=0.8, label='Continuum', rasterized=True
                    )
                else:
                    logging.warning("Continuum requested but 'cont' not found.")

            # 4. Plot Model
            if self.main_window.model_checkbox.isChecked():
                if model_data is not None:
                    ax.plot(x_data, model_data.value, ls='-', color=colors[1], lw=0.8, label='Absorption model', rasterized=True)
                # No warning needed

            # 5. Plot Systems
            if self.main_window.systems_checkbox.isChecked() and V1_FUNCTIONS_AVAILABLE and systs and systs.components:
                # Get current plot limits to only draw visible lines
                xlim = ax.get_xlim()
                added_hi_label = False
                added_metal_label = False

                # Define Y position for markers (e.g., slightly below top, in axis coords)
                marker_y_axis_coord = 0.05 # 5% from the bottom
                
                trans = ax.get_xaxis_transform()

                # ** Create lists to batch coordinates **
                hi_lines_x = []
                metal_lines_x = []
                zem = getattr(spec, '_zem', 0.0) # Get zem once

                for comp in systs.components:
                    z = comp.z
                    # Use V1 function to get transitions for the series
                    transitions = trans_parse(comp.series)
                    for t in transitions:
                        if t in xem_d:
                            try:
                                x_plot = x_convert((1 + z) * xem_d[t].to_value(au.nm) * au.nm, zem=zem, xunit=spec.x.unit).value
                            except Exception: continue

                            # Check xlim *before* appending to keep lists small
                            if xlim[0] <= x_plot <= xlim[1]:
                                if t.startswith('Ly_'):
                                    hi_lines_x.append(x_plot)
                                else:
                                    metal_lines_x.append(x_plot)
                        else:
                            logging.warning(f"Transition '{t}' for series '{comp.series}' not found in xem_d.")

                # ** Make only TWO plot calls, outside the loop **
                if hi_lines_x:
                    ax.plot(hi_lines_x, [marker_y_axis_coord] * len(hi_lines_x),
                            marker='|', markersize=12, linestyle='None',
                            color='red', alpha=0.9, label="HI Systems",
                            transform=trans)
                if metal_lines_x:
                    ax.plot(metal_lines_x, [marker_y_axis_coord] * len(metal_lines_x),
                            marker='|', markersize=12, linestyle='None',
                            color='darkgray', alpha=0.9, label="Metal Systems",
                            transform=trans)

            # 6. Plot Redshift Cursor
            if self.main_window.cursor_show_checkbox.isChecked() and V1_FUNCTIONS_AVAILABLE:
                try:
                    z_cursor = float(self.main_window.cursor_z_input.text())
                    # Calculate positions using helper
                    cursor_x_positions = self.get_cursor_line_positions_at_z(z_cursor)
                    added_cursor_label = False
                    colors = get_color_cycle(5) # Get colors

                    for x_plot in cursor_x_positions:
                        # ** REMOVED check against xlim_plot **
                        # if xlim_plot[0] <= x_plot <= xlim_plot[1]:
                        label = f"Cursor (z={z_cursor:.4f})" if not added_cursor_label else None
                        line = self.canvas.axes.axvline(x_plot, ls='--', color=colors[3],
                                                         alpha=0.8, lw=1.0, label=label,
                                                         animated=True) # <<< ANIMATED
                        self.canvas.cursor_artists.append(line)
                        added_cursor_label = True # Only label first line potentially visible
                        # *************************************
                    logging.debug(f"plot_spectrum: Added {len(self.canvas.cursor_artists)} cursor artists.")
                except (ValueError, AttributeError) as e:
                    logging.warning(f"Could not draw cursor lines: {e}")
                except Exception as e:
                    logging.error(f"Unexpected error drawing cursor: {e}", exc_info=True)  

            # 4. Final Touches
            if ax.has_data(): # Only add legend if something was plotted
                ax.legend(loc=legend_loc, markerscale=0.2) # Smaller legend font
            ax.set_xlabel(f"Wavelength ({x_unit})")
            ax.set_ylabel(f"Flux ({y_unit})")
            ax.set_title(f"Spectrum: {session_state.name}")
            ax.grid(True, linestyle=':')

            # 1. Scales (Log/Linear)
            ax.set_xscale('log' if is_log_x else 'linear')
            ax.set_yscale('log' if is_log_y else 'linear')

            # 2. X-Axis Formatting
            xlabel = "Wavelength (nm)" # Default
            if selected_x_unit == 'Angstrom':
                ax.xaxis.set_major_formatter(plt_ticker.FuncFormatter(lambda x, pos: f"{x * 10:g}"))
                xlabel = "Wavelength (Angstrom)"
            elif selected_x_unit == 'um':
                ax.xaxis.set_major_formatter(plt_ticker.FuncFormatter(lambda x, pos: f"{x / 1000:g}"))
                xlabel = "Wavelength (μm)"
            # else: no formatter needed for nm
            ax.set_xlabel(xlabel)

            # 3. Y-Axis Label
            ax.set_ylabel("Normalized Flux" if is_norm_y else f"Flux ({y_unit_str})")
            
            # 4. Set Limits
            # X-Axis: Restore zoom if needed
            if was_zoomed and target_xlim is None: ax.set_xlim(previous_xlim)
            elif target_xlim is not None: ax.set_xlim(target_xlim)
            # else: autoscale (Matplotlib default)
            
            # Y-Axis: Handle normalization OR restore zoom
            if is_norm_y:
                ax.set_ylim(-0.3, 1.3) # Apply fixed limits
            elif target_ylim is not None: ax.set_ylim(target_ylim)
            elif was_zoomed_y: ax.set_ylim(previous_ylim)
            # else: autoscale (Matplotlib default)

        if not plot_occurred:
            ax = self.canvas.axes # Ensure ax is defined
            ax.clear() # Clear axes even if no data
            ax.set_title("No Spectrum Data Loaded") 

        # Disconnect any previous connection first
        if self.canvas.draw_event_cid is not None:
            try: self.canvas.mpl_disconnect(self.canvas.draw_event_cid)
            except Exception: pass # Ignore if already disconnected
        self.canvas.draw_event_cid = self.canvas.mpl_connect('draw_event', self.canvas._capture_background)

        # Redraw the canvas
        if initial_draw:
            self.canvas.draw()
        else:
            self.canvas.draw_idle()

        # NOTE: You must update the self.session reference here if the session changes.
        # This function will be called by update_plot().

    def get_cursor_line_positions_at_z(self, z_cursor):
        """Helper to get theoretical X positions for cursor lines at a specific Z."""
        positions = []
        spec = session_state.spec
        main_window = self.main_window # type: MainWindowV2
        if not spec or not V1_FUNCTIONS_AVAILABLE: return []
        try:
            series_str = main_window.cursor_series_input.text()
            transitions = trans_parse(series_str)
            zem = getattr(spec, '_zem', 0.0) # Get spec RF z
            x_unit = spec.x.unit # Get spec current x unit

            for t in transitions:
                if t in xem_d:
                    xem_nm = xem_d[t].to_value(au.nm)
                    x_plot_nm = (1 + z_cursor) * xem_nm # This is the data coordinate
                    positions.append(x_plot_nm)
        except (ValueError, AttributeError, ImportError) as e:
            logging.warning(f"Could not calculate cursor positions: {e}")
            return []
        return positions
    
    def update_plot(self, new_session_state: Optional['SessionV2']):
        """Called by MainWindowV2 to swap the immutable session object and redraw."""
        self.plot_spectrum(session_state=new_session_state, initial_draw=False)

# NOTE: The rest of the plotting logic (plot_spectrum content, resizeEvent, etc.) 
# must be placed into this class, and the old SpectrumViewerPySide class should be deleted.

class AstrocookToolbar(NavigationToolbar):
    """
    Custom Matplotlib toolbar that inherits from NavigationToolbar2QT 
    to remove unwanted buttons and simplify the interface.
    """
    # Define the order of toolbar items to keep (based on V1 functionality)
    toolitems = [
        # Label, Tooltip, Icon Filename, Tool ID
        ('Home', 'Reset original view', 'home', 'home'),
        ('Back', 'Back to previous view', 'back', 'back'),
        ('Forward', 'Forward to next view', 'forward', 'forward'),
        (None, None, None, None), # Separator placeholder
        
        # FIX: Use Matplotlib's internal symbolic names for icons (third element).
        # We rely on the base class finding 'hand' and 'magnifying_glass' equivalent icons.
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'), 
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
        
        (None, None, None, None), # Separator placeholder
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
    ]

    def __init__(self, canvas, parent):
        # 1. Base Class Initialization: This CREATES and stores the C++ QAction objects.
        # This MUST run first.
        super().__init__(canvas, parent, coordinates=True) 
        self.viewer_parent = parent

        # 2. Final cleanup (Remove unwanted actions created by the base class)
        actions_to_remove = ['Customize', 'Subplots', 'Edit axis, curve and image parameters']
        
        for action in self.actions():
            if action.text() in actions_to_remove or action.iconText() in actions_to_remove:
                self.removeAction(action)


# --- Utility Function to Launch the Viewer for Testing ---

def launch_viewer(session: 'SessionV2'):
    """Utility to run the viewer application."""
    import sys
    
    # Check if a Qt application instance is already running (crucial for integration)
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
        
    viewer = SpectrumViewerPySide(session)
    viewer.show()
    
    # If the app was newly created, start the event loop
    if app and QApplication.instance() == app:
        sys.exit(app.exec())