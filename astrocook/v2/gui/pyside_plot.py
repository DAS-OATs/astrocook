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
import matplotlib.ticker as plt_ticker
import numpy as np
from PySide6.QtCore import QPoint, Qt, QTimer
from PySide6.QtGui import QAction, QCursor, QIcon
from PySide6.QtWidgets import QApplication, QMenu, QStyle, QToolTip, QVBoxLayout, QWidget
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
def get_color_cycle(n=10, cmap=None, fallback_cmap='tab20'):
    # ... (function definition as before) ...
    if cmap is not None:
        try:
            cmap_obj = plt.get_cmap(cmap)
            return [cmap_obj(i / n) for i in range(n)]
        except Exception as e:
            logging.warning(f"Could not get specified colormap '{cmap}': {e}. Using fallback.")
    try:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        return [colors[i % len(colors)] for i in range(n)]
    except Exception as e:
        logging.warning(f"Could not get color cycle: {e}. Using fallback.")
        cmap = plt.get_cmap(fallback_cmap); return [cmap(i / n) for i in range(n)]
# ----------------------------

def decimate_y_min_max(y: np.ndarray, factor: int) -> np.ndarray:
    """
    Performs a min-max decimation on a Y-axis array.
    Returns an array of [y_min1, y_max1, y_min2, y_max2, ...]
    """
    if y is None:
        return None
    try:
        # 1. Trim array to be divisible by factor
        size = (len(y) // factor) * factor
        if size == 0:
            return y # Not enough data to decimate, return original
            
        y_trim = y[:size]
        
        # 2. Reshape
        y_reshaped = y_trim.reshape(-1, factor)
        
        # 3. Get min/max for y
        y_mins = y_reshaped.min(axis=1)
        y_maxs = y_reshaped.max(axis=1)
        
        # 4. Interleave
        return np.stack((y_mins, y_maxs), axis=1).ravel()
    except Exception as e:
        logging.warning(f"y-axis decimation failed: {e}")
        return y # Return original on failure
        
def decimate_x_for_min_max(x: np.ndarray, factor: int) -> np.ndarray:
    """
    Decimates an X-axis to match a min-max-decimated Y-axis.
    Returns an array of [x_start1, x_end1, x_start2, x_end2, ...]
    """
    try:
        # 1. Trim array to be divisible by factor
        size = (len(x) // factor) * factor
        if size == 0:
            return x # Not enough data to decimate, return original

        x_trim = x[:size]
        
        # 2. Reshape
        x_reshaped = x_trim.reshape(-1, factor)
        
        # 3. Get start/end for x
        x_starts = x_reshaped[:, 0]  # Get the start x of each bin
        x_ends = x_reshaped[:, -1] # Get the end x of each bin
        
        # 4. Interleave
        return np.stack((x_starts, x_ends), axis=1).ravel()
    except Exception as e:
        logging.warning(f"x-axis decimation failed: {e}")
        return x # Return original on failure

class MatplotlibCanvas(FigureCanvasQTAgg):
    """A custom Matplotlib canvas enabling blitting for cursor dragging."""
    def __init__(self, parent=None, width=5, height=4, dpi=100, plot_widget=None): # Added plot_widget ref
        try:
            #plt.style.use(['science', 'fast'])
            plt.style.use('fast')
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

        # Hover Timer for Tooltips
        self.hover_timer = QTimer(self)
        self.hover_timer.setSingleShot(True)
        self.hover_timer.setInterval(500) # 500ms delay before showing tooltip
        self.hover_timer.timeout.connect(self._on_hover_timeout)
        self._last_mouse_event = None

        # Region Selection State
        self.selection_artist = None # Will hold the axvspan
        self.selection_start_x = None

        # Connect Matplotlib events
        self.mpl_connect('motion_notify_event', self.on_motion)
        self.mpl_connect('button_press_event', self.on_press)
        self.mpl_connect('button_release_event', self.on_release)
        self.xlim_cid = self.axes.callbacks.connect('xlim_changed', self.on_lim_changed)
        self.ylim_cid = self.axes.callbacks.connect('ylim_changed', self.on_lim_changed)


    def on_lim_changed(self, axes):
        """Callback: Invalidate background and schedule a full redraw."""
        if self.background is not None and self.plot_widget: # Check for plot_widget
            logging.debug(f"Limits changed ({axes.get_label()}): Invalidating background & scheduling plot_spectrum.")
            self.background = None # Invalidate background cache

            # Get the current state to pass to the scheduled plot call
            session_state = None
            if self.plot_widget.main_window and self.plot_widget.main_window.active_history:
                session_state = self.plot_widget.main_window.active_history.current_state

            # Schedule plot_spectrum to run soon via the event loop
            # ** Use lambda to pass the session_state argument **
            QTimer.singleShot(0, lambda: self.plot_widget.plot_spectrum(session_state=session_state))

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


    def on_motion(self, event):
        """Handles mouse motion: Updates cursor position if active."""
        # 1. Always restart the hover timer on any motion
        # Immediately hide any existing tooltip when the mouse moves.
        # It will reappear only when the user stops again (timer timeout).
        QToolTip.hideText()
        
        self._last_mouse_event = event
        self.hover_timer.start()

        # Safety checks
        if (not self.plot_widget or not hasattr(self.plot_widget, 'main_window')
            or not self.plot_widget.main_window):
            # logging.debug("on_motion: plot_widget or main_window invalid.")
            return

        main_window = self.plot_widget.main_window # type: MainWindowV2

        # Update Selection Span while dragging
        if (self.plot_widget and self.plot_widget._is_selecting_region and
            self.selection_start_x is not None and event.inaxes == self.axes):
            
            # Update the span's extent
            # axvspan is a Polygon, we convert it to a simple Rect for easy updating if needed,
            # but .set_xy is the standard way to update a Polygon.
            # Actually, for axvspan it's easier to just update the vertices.
            # A simpler way for axvspan: it works by setting xmin/xmax.
            # BUT axvspan doesn't have simple set_xlim methods.
            # Re-drawing might be easiest, or using a Rectangle instead.
            
            # Let's use a simpler approach: remove old, draw new (fast enough for this)
            if self.selection_artist:
                self.selection_artist.remove()
            
            self.selection_artist = self.axes.axvspan(
                self.selection_start_x, event.xdata, color='orange', alpha=0.3
            )
            self.draw_idle()
            return # Don't do cursor updates while selecting
        
        # Check if inside axes and cursor checkbox is checked
        if (event.inaxes == self.axes and event.xdata is not None
            and main_window.cursor_show_checkbox.isChecked()):

            session_state = None
            if main_window.active_history:
                session_state = main_window.active_history.current_state

            if not session_state or not session_state.spec:
                logging.debug(
                    "on_motion: No active session state or spec. Skipping.")
                return  # Can't do conversions without a spec
            
            # --- Always update on motion if cursor is active ---
            mouse_x = event.xdata
            spec = session_state.spec
            if not spec: return

            try:
                new_z = self.plot_widget.calculate_z_from_x(event.xdata)

                if new_z is not None:
                    # Update QLineEdit (no need for blockSignals if not dragging)
                    current_text_z = None
                    try: current_text_z = float(main_window.cursor_z_input.text())
                    except ValueError: pass
                    # Update text only if significantly different to avoid excessive signals
                    if current_text_z is None or not np.isclose(new_z, current_text_z, atol=1e-7):
                        main_window.cursor_z_input.setText(f"{new_z:.7f}")

                    # Update Artist Positions & Redraw using Blit
                    new_x_positions = self.plot_widget.get_cursor_line_positions_at_z(session_state, new_z)

                    # Ensure artists exist and match count
                    if not self.cursor_artists or len(new_x_positions) != len(self.cursor_artists):
                        # Force a full redraw if artists are missing/mismatched
                        logging.warning("Cursor artists missing or mismatched on motion. Forcing full redraw.")
                        self.plot_widget.plot_spectrum(session_state=session_state) # Recreate artists
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

    def _on_hover_timeout(self):
        """Called when the mouse has stopped moving for a while."""
        if self.plot_widget and self._last_mouse_event and self._last_mouse_event.inaxes == self.axes:
            self.plot_widget.show_data_tooltip(self._last_mouse_event)

    def on_press(self, event):
        """
        Handle mouse clicks. Right-click (Button 3) opens a context menu
        with actions relevant to the clicked position.
        """
        if not self.plot_widget or not self.plot_widget.main_window: return
        main_win = self.plot_widget.main_window

        if (event.button == 1 and event.inaxes == self.axes and 
            self.plot_widget._is_selecting_region):
            
            self.selection_start_x = event.xdata
            # Create the visual artist (invisible at first, 0 width)
            self.selection_artist = self.axes.axvspan(
                event.xdata, event.xdata, color='orange', alpha=0.3
            )
            self.draw_idle() # Force a draw to show it starts
            return # Don't do other click actions

        # Only handle Right-Click (Button 3) inside the axes
        if event.button == 3 and event.inaxes == self.axes:
            
            menu = QMenu(self)
            
            # --- Feature 1: Set z_em (only if cursor is active) ---
            if main_win.cursor_show_checkbox.isChecked():
                new_z = self.plot_widget.calculate_z_from_x(event.xdata)
                if new_z is not None:
                    # Create the action with a clear label
                    z_action = QAction(f"Set emission redshift to z={new_z:.5f}", menu)
                    z_action.triggered.connect(lambda: main_win._on_recipe_requested(
                        category="edit",
                        recipe_name="set_properties",
                        params={"z_em": str(new_z)},
                        alias_map={}
                    ))
                    menu.addAction(z_action)

            # --- (Future features will go here: Fit Line, Mask Region, etc.) ---

            # Show the menu if it has any actions
            if not menu.isEmpty():
                # Use QCursor.pos() to show menu at the global mouse coordinates
                menu.exec(QCursor.pos())

    def on_release(self, event):
        """Handle mouse button release events."""
        # Check if we were selecting a region
        if self.selection_start_x is not None:
            start = self.selection_start_x
            end = event.xdata
            
            # 1. Cleanup the visual
            self.cleanup_selection()
            
            # 2. Validate the selection
            if start is None or end is None or start == end:
                 return # Just a click, not a drag
            
            # Ensure min/max order
            xmin, xmax = sorted([start, end])
            
            # 3. Trigger the action in the main window
            if self.plot_widget and self.plot_widget.main_window:
                logging.info(f"Region selected: {xmin:.4f} to {xmax:.4f}")
                # Call a new helper in main_window to launch split
                self.plot_widget.main_window.launch_split_from_region(xmin, xmax)
                
            # 4. Optional: Auto-exit select mode?
            # self.plot_widget.toggle_region_selector() 

    def cleanup_selection(self):
        """Helper to remove selection artist and reset state."""
        self.selection_start_x = None
        if self.selection_artist:
            try:
                self.selection_artist.remove()
            except Exception: pass # Might already be gone
            self.selection_artist = None
            self.draw_idle()

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

        # ** Store the last-drawn normalization state **
        self._last_draw_norm_state = False # Default to non-normalized
        self._last_draw_snr_state = False

        self._cached_data_id = None # Tracks the ID of the cached data

        # --- NEW: Cache for RAW plot-ready (decimated) arrays ---
        self._cached_x_plot_raw = None
        self._cached_y_plot_raw = None
        self._cached_dy_plot_raw = None
        self._cached_cont_plot_raw = None
        self._cached_model_plot_raw = None

        # Cache for NORMALIZED plot-ready (decimated) arrays
        self._cached_y_plot_norm = None
        self._cached_dy_plot_norm = None
        self._cached_cont_plot_norm = None
        self._cached_model_plot_norm = None

        # Cache for SNR plot-ready (decimated) arrays
        self._cached_y_plot_snr = None
        self._cached_dy_plot_snr = None
        self._cached_model_plot_snr = None # (e.g., model / dy)
        self._cached_cont_plot_snr = None # (e.g., cont / dy)
        self._cached_snr_error_col = "" # Tracks which error col was used

        # Selection Mode State
        self._is_selecting_region = False

        # Add ONLY toolbar to main layout
        self.main_layout.addWidget(self.toolbar, 0) # Stretch 0

        # Draw initial plot (will read checkbox state from main window)
        self.plot_spectrum(session_state=initial_session_state, initial_draw=True)

    # --- Refactored Plotting Method (Requires an update) ---
    def plot_spectrum(self, session_state: Optional['SessionV2'], initial_draw=False, force_autoscale: bool = False):
        """
        Retrieves data from the immutable V2 Session and plots it.
        """
        # Get the unique ID of the immutable data core
        current_data_id = None
        if session_state and session_state.spec:
            # The ._data object is the frozen SpectrumDataV2
            current_data_id = id(session_state.spec._data) 

        # --- Get the selected SNR error column ---
        selected_snr_col = self.main_window.snr_col_combo.currentText()

        # If data changed OR the SNR error col changed, invalidate caches
        if (current_data_id != self._cached_data_id or 
            force_autoscale or
            (self.main_window.snr_checkbox.isChecked() and 
             self._cached_snr_error_col != selected_snr_col)
           ):
            
            if force_autoscale:
                logging.debug("Home button pressed, invalidating plot caches.")
            elif current_data_id != self._cached_data_id:
                logging.debug("Data change detected, invalidating ALL plot caches.")
            else:
                logging.debug("SNR error column changed, invalidating plot caches.")

            self._cached_data_id = current_data_id
            # --- Invalidate ALL caches ---
            self._cached_x_plot_raw = None
            self._cached_y_plot_raw = None
            self._cached_dy_plot_raw = None
            self._cached_cont_plot_raw = None
            self._cached_model_plot_raw = None
            self._cached_y_plot_norm = None
            self._cached_dy_plot_norm = None
            self._cached_cont_plot_norm = None
            self._cached_model_plot_norm = None
            self._cached_y_plot_snr = None
            self._cached_dy_plot_snr = None
            self._cached_cont_plot_snr = None
            self._cached_model_plot_snr = None
            self._cached_snr_error_col = "" # Clear the cached col name

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
        #if previous_ylim == (-0.3, 1.3): previous_ylim = (0.0, 1.0) # Normalize special case
        if force_autoscale:
            was_zoomed = False
            logging.debug("Home button pressed: Forcing autoscale.")
        else:
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
            full_cont_data = spec.cont.value if hasattr(spec.cont, 'value') else spec.cont
            full_model_data = spec.model.value if hasattr(spec.model, 'value') else spec.model

            # --- Check View Toggles ---
            is_norm_y = self.main_window.norm_y_checkbox.isChecked()
            is_snr = self.main_window.snr_checkbox.isChecked() # <-- NEW
            is_log_x = self.main_window.log_x_checkbox.isChecked()
            is_log_y = self.main_window.log_y_checkbox.isChecked()
            selected_x_unit = self.main_window.x_unit_combo.currentText()

            # ** Check if normalization state CHANGED since last draw **
            norm_state_changed = (is_norm_y != self._last_draw_norm_state) or \
                             (is_snr != self._last_draw_snr_state)
            if norm_state_changed:
                logging.debug(f"View state changed (Norm: {is_norm_y}, SNR: {is_snr})")
            # ** Update the stored state for the next draw **
            self._last_draw_norm_state = is_norm_y
            self._last_draw_snr_state = is_snr
        
            data_slice = slice(None) # Default: use full array

            # Define decimation parameters
            DECIMATION_THRESHOLD = 2000000 # Only decimate if plot has > 20k points
            DECIMATION_FACTOR = 10       # Plot 2 points for every 10
            
            use_decimation = False
            
            if not was_zoomed and len(full_x_data) > DECIMATION_THRESHOLD and not force_autoscale:
                use_decimation = True
                logging.debug(f"Using Min-Max decimation (factor {DECIMATION_FACTOR}) for full-view plot.")

            if is_norm_y:
                # 1a. We want NORMALIZED data. Check cache.
                if self._cached_y_plot_norm is not None:
                    logging.debug("Using cached NORMALIZED plot data.")
                    x_data = self._cached_x_plot_raw
                    y_data = self._cached_y_plot_norm
                    dy_data = self._cached_dy_plot_norm
                    cont_data = self._cached_cont_plot_norm
                    model_data = self._cached_model_plot_norm
                
                else:
                    logging.debug("Cache miss for normalized data. Computing...")
                    # 2a. Need to COMPUTE. First, get RAW data.
                    if self._cached_x_plot_raw is not None:
                        logging.debug("...using cached RAW plot data.")
                        x_data_raw = self._cached_x_plot_raw
                        y_data_raw = self._cached_y_plot_raw
                        dy_data_raw = self._cached_dy_plot_raw
                        cont_data_raw = self._cached_cont_plot_raw
                        model_data_raw = self._cached_model_plot_raw
                    else:
                        logging.debug("...RAW plot data cache is also empty. Computing...")
                        # 3a. TOTAL cache miss. Compute RAW data.
                        if use_decimation:
                            x_data_raw = decimate_x_for_min_max(full_x_data, DECIMATION_FACTOR)
                            y_data_raw = decimate_y_min_max(full_y_data, DECIMATION_FACTOR)
                            dy_data_raw = decimate_y_min_max(full_dy_data, DECIMATION_FACTOR)
                            cont_data_raw = decimate_y_min_max(full_cont_data, DECIMATION_FACTOR)
                            model_data_raw = decimate_y_min_max(full_model_data, DECIMATION_FACTOR)
                        else: # We are zoomed, forcing home, or data is small
                            x_data_raw = full_x_data
                            y_data_raw = full_y_data
                            dy_data_raw = full_dy_data
                            cont_data_raw = full_cont_data
                            model_data_raw = full_model_data
                        
                        # Save to RAW cache
                        self._cached_x_plot_raw = x_data_raw
                        self._cached_y_plot_raw = y_data_raw
                        self._cached_dy_plot_raw = dy_data_raw
                        self._cached_cont_plot_raw = cont_data_raw
                        self._cached_model_plot_raw = model_data_raw

                    # 4a. NOW, compute and cache the NORMALIZED data
                    if cont_data_raw is not None:
                        x_data = x_data_raw
                        cont_val = cont_data_raw
                        y_data = np.divide(y_data_raw, cont_val, out=np.full_like(y_data_raw, np.nan), where=cont_val!=0)
                        dy_data = np.divide(dy_data_raw, cont_val, out=np.full_like(dy_data_raw, np.nan), where=cont_val!=0)
                        if model_data_raw is not None:
                            model_data = np.divide(model_data_raw, cont_val, out=np.full_like(model_data_raw, np.nan), where=cont_val!=0)
                        cont_data = np.ones_like(cont_data_raw)
                        
                        # Save to NORMALIZED cache
                        self._cached_x_plot_norm = x_data
                        self._cached_y_plot_norm = y_data
                        self._cached_dy_plot_norm = dy_data
                        self._cached_cont_plot_norm = cont_data
                        self._cached_model_plot_norm = model_data
                    else:
                        logging.warning("Normalize checked, but continuum data is missing.")
                        x_data, y_data, dy_data, cont_data, model_data = \
                            x_data_raw, y_data_raw, dy_data_raw, cont_data_raw, model_data_raw

            elif is_snr:
                # --- *** NEW BLOCK FOR SNR *** ---
                if self._cached_y_plot_snr is not None:
                    logging.debug("Using cached SNR plot data.")
                    x_data = self._cached_x_plot_raw
                    y_data = self._cached_y_plot_snr
                    dy_data = self._cached_dy_plot_snr
                    cont_data = self._cached_cont_plot_snr
                    model_data = self._cached_model_plot_snr
                
                else:
                    logging.debug("Cache miss for SNR data. Computing...")
                    # 1. Get RAW plot data (from cache or compute)
                    if self._cached_x_plot_raw is not None:
                        # ... (read from _cached_x_plot_raw, etc.) ...
                        x_data_raw = self._cached_x_plot_raw
                        y_data_raw = self._cached_y_plot_raw
                        dy_data_raw = self._cached_dy_plot_raw
                        cont_data_raw = self._cached_cont_plot_raw
                        model_data_raw = self._cached_model_plot_raw
                    else:
                        # ... (compute raw decimated data) ...
                        if use_decimation:
                            # ... (decimate_x_for_min_max, decimate_y_min_max) ...
                            x_data_raw = decimate_x_for_min_max(full_x_data, DECIMATION_FACTOR)
                            y_data_raw = decimate_y_min_max(full_y_data, DECIMATION_FACTOR)
                            dy_data_raw = decimate_y_min_max(full_dy_data, DECIMATION_FACTOR)
                            cont_data_raw = decimate_y_min_max(full_cont_data, DECIMATION_FACTOR)
                            model_data_raw = decimate_y_min_max(full_model_data, DECIMATION_FACTOR)
                        else:
                            # ... (get from data_slice) ...
                            x_data_raw = full_x_data[data_slice]
                            y_data_raw = full_y_data[data_slice]
                            dy_data_raw = full_dy_data[data_slice]
                            cont_data_raw = full_cont_data[data_slice] if full_cont_data is not None else None
                            model_data_raw = full_model_data[data_slice] if full_model_data is not None else None
                        
                        # Save to RAW cache
                        self._cached_x_plot_raw = x_data_raw
                        self._cached_y_plot_raw = y_data_raw
                        self._cached_dy_plot_raw = dy_data_raw
                        self._cached_cont_plot_raw = cont_data_raw
                        self._cached_model_plot_raw = model_data_raw

                    # 2. Get the ERROR array to use
                    error_col_name = self.main_window.snr_col_combo.currentText()
                    full_error_data = None
                    try:
                        if error_col_name == 'dy':
                            full_error_data = full_dy_data
                        else:
                            # This path is for 'rms', 'cont_err', etc.
                            full_error_data = spec.get_column(error_col_name).value
                    except Exception as e:
                        logging.error(f"Could not get SNR error column '{error_col_name}': {e}. Falling back to 'dy'.")
                        full_error_data = full_dy_data # Fallback
                        self.main_window.snr_col_combo.setCurrentIndex(0) # Reset combo
                    
                    # 3. Decimate/slice the error array
                    if use_decimation:
                        error_data_raw = decimate_y_min_max(full_error_data, DECIMATION_FACTOR)
                    else:
                        error_data_raw = full_error_data[data_slice]
                    
                    # 4. NOW, compute and cache the SNR data
                    y_data = np.divide(y_data_raw, error_data_raw, out=np.full_like(y_data_raw, np.nan), where=error_data_raw!=0)
                    dy_data = np.ones_like(y_data) # Error on SNR is (dy/dy) = 1.0
                    
                    cont_data = None
                    if cont_data_raw is not None:
                        cont_data = np.divide(cont_data_raw, error_data_raw, out=np.full_like(cont_data_raw, np.nan), where=error_data_raw!=0)
                    
                    model_data = None
                    if model_data_raw is not None:
                        model_data = np.divide(model_data_raw, error_data_raw, out=np.full_like(model_data_raw, np.nan), where=error_data_raw!=0)

                    # Save to SNR cache
                    x_data = x_data_raw # X is always the same
                    self._cached_y_plot_snr = y_data
                    self._cached_dy_plot_snr = dy_data
                    self._cached_cont_plot_snr = cont_data
                    self._cached_model_plot_snr = model_data
                    self._cached_snr_error_col = error_col_name # Store the col used

            else:
                # 1b. We want RAW data. Check cache.
                if self._cached_x_plot_raw is not None:
                    logging.debug("Using cached RAW plot data.")
                    x_data = self._cached_x_plot_raw
                    y_data = self._cached_y_plot_raw
                    dy_data = self._cached_dy_plot_raw
                    cont_data = self._cached_cont_plot_raw
                    model_data = self._cached_model_plot_raw
                
                else:
                    logging.debug("Cache miss for raw plot data. Computing...")
                    # 2b. Compute RAW data.
                    if use_decimation:
                        x_data = decimate_x_for_min_max(full_x_data, DECIMATION_FACTOR)
                        y_data = decimate_y_min_max(full_y_data, DECIMATION_FACTOR)
                        dy_data = decimate_y_min_max(full_dy_data, DECIMATION_FACTOR)
                        cont_data = decimate_y_min_max(full_cont_data, DECIMATION_FACTOR)
                        model_data = decimate_y_min_max(full_model_data, DECIMATION_FACTOR)
                    else: # We are zoomed, forcing home, or data is small
                        x_data = full_x_data
                        y_data = full_y_data
                        dy_data = full_dy_data
                        cont_data = full_cont_data
                        model_data = full_model_data
                    
                    # Save to RAW cache
                    self._cached_x_plot_raw = x_data
                    self._cached_y_plot_raw = y_data
                    self._cached_dy_plot_raw = dy_data
                    self._cached_cont_plot_raw = cont_data
                    self._cached_model_plot_raw = model_data

            # --- 2. SLICING BLOCK (NEW) ---
            # At this point, x_data, y_data, etc., are ALL full-range.
            # Now we apply the zoom-slice if necessary.
            
            final_slice = slice(None) # Default to no extra slicing

            if was_zoomed and not force_autoscale:
                try:
                    # Find indices for the visible range *on the plot data*
                    idx_start = np.searchsorted(x_data, previous_xlim[0], side='left')
                    idx_end = np.searchsorted(x_data, previous_xlim[1], side='right')

                    # Add a buffer (in *indices*, not data units)
                    # 10 bins for decimated, 50 for full
                    buffer = (DECIMATION_FACTOR * 2) if use_decimation else 50
                    idx_start = max(0, idx_start - buffer)
                    idx_end = min(len(x_data), idx_end + buffer)

                    if idx_end > idx_start:
                        logging.debug(f"Applying zoom-slice to plot data: {idx_end - idx_start} points")
                        
                        final_slice = slice(idx_start, idx_end)

                        # Apply the slice to all arrays
                        x_data = x_data[final_slice]
                        y_data = y_data[final_slice]
                        dy_data = dy_data[final_slice]
                        if cont_data is not None:
                            cont_data = cont_data[final_slice]
                        if model_data is not None:
                            model_data = model_data[final_slice]
                    else:
                        logging.debug("Zoom slice is empty, plotting full data.")
                except Exception as e:
                    logging.warning(f"Failed to apply zoom slice: {e}")

            x_unit = str(spec.x.unit)
                    
            colors = get_color_cycle(5, cmap='tab20') # Get colors

            # --- Check state from main_window reference ---
            # 2. Plot Error Shading (Conditional)
            if self.main_window.error_checkbox.isChecked(): # <<< Check main window's checkbox
                #ax.fill_between(
                #    x_data, y_data - dy_data, y_data + dy_data,
                #    step='mid', color='#aaaaaa', alpha=0.5,
                #    label='1-sigma error', rasterized=True
                #)
                if dy_data is not None:
                    label = '1-sigma error'
                    if use_decimation:
                        # Use plot() for decimated data (fast)
                        ax.plot(x_data, y_data - dy_data, lw=0.3, color='#aaaaaa', rasterized=True, label=label)
                        ax.plot(x_data, y_data + dy_data, lw=0.3, color='#aaaaaa', rasterized=True)
                    else:
                        # Use step() for non-decimated data (correct)
                        ax.step(x_data, y_data - dy_data, where='mid', lw=0.3, color='#aaaaaa', rasterized=True, label=label)
                        ax.step(x_data, y_data + dy_data, where='mid', lw=0.3, color='#aaaaaa', rasterized=True)
                else:
                    logging.warning("Skipping error plot because dy_data is None (this may be a cache bug).")

            # 1. Plot Main Flux (Uses Matplotlib style defaults for color)
            if use_decimation:
                ax.plot(x_data, y_data, label="Spectrum (decimated)", lw=0.5, color=colors[0], rasterized=True)
            else:
                ax.step(x_data, y_data, where='mid', label="Spectrum", lw=0.5, color=colors[0], rasterized=True)

            # 3. Plot Continuum (Conditional)
            if self.main_window.continuum_checkbox.isChecked(): # <<< Check main window's checkbox
                if cont_data is not None:
                    ax.plot(
                        x_data, cont_data, linestyle='--',
                        color='black', lw=0.8, label='Continuum', rasterized=True
                    )
                else:
                    logging.warning("Continuum requested but 'cont' not found.")

            # 4. Plot Model
            if self.main_window.model_checkbox.isChecked():
                if model_data is not None:
                    ax.plot(x_data, model_data, ls='-', color=colors[1], lw=0.8, label='Absorption model', rasterized=True)
                # No warning needed

            # 4b. Plot Auxiliary Column (Dynamic)
            aux_col_to_plot = self.main_window.aux_col_combo.currentText()
            if aux_col_to_plot and aux_col_to_plot != "None":
                try:
                    # Get the *full* original aux column data
                    full_aux_col_data = spec.get_column(aux_col_to_plot)
                    if full_aux_col_data is None:
                        raise ValueError(f"Column '{aux_col_to_plot}' returned None.")
                    
                    full_aux_col_data = full_aux_col_data.value # Get numpy array
                    
                    # 1. Apply Decimation (if it was applied to x_data)
                    if use_decimation:
                        aux_data = decimate_y_min_max(full_aux_col_data, DECIMATION_FACTOR)
                    else:
                        aux_data = full_aux_col_data
                        
                    # 2. Apply Zoom Slicing (must match x_data's final state)
                    aux_data = aux_data[final_slice]
                    
                    # Plot based on data type
                    if aux_data.dtype.kind in 'fi': # Float or Int
                        if use_decimation:
                            ax.plot(x_data, aux_data, label=aux_col_to_plot, 
                                    linestyle=':', lw=1.2, color='purple', 
                                    rasterized=True)
                        else:
                            ax.step(x_data, aux_data, where='mid', label=aux_col_to_plot, 
                                    linestyle=':', lw=1.2, color='purple', 
                                    rasterized=True)
                                    
                    elif aux_data.dtype.kind == 'b': # Boolean mask
                        trans = ax.get_xaxis_transform()
                        ax.fill_between(x_data, 0.95, 1.0, 
                                        where=aux_data.astype(bool), 
                                        color='cyan', alpha=0.5, 
                                        label=f'Mask: {aux_col_to_plot}', 
                                        transform=trans,
                                        step='mid' if not use_decimation else None,
                                        rasterized=True)

                except Exception as e:
                    logging.warning(f"Could not plot auxiliary column '{aux_col_to_plot}': {e}")

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
                    cursor_x_positions = self.get_cursor_line_positions_at_z(session_state, z_cursor)
                    added_cursor_label = False
                    
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
            ax.set_title(f"Spectrum: {session_state.name}")
            
            # 1. Scales (Log/Linear)
            ax.set_xscale('log' if is_log_x else 'linear')
            ax.set_yscale('log' if is_log_y else 'linear')

            ax.grid(True, linestyle=':')
            
            # 2. X-Axis Formatting
            xlabel = "Wavelength (nm)" # Default
            if selected_x_unit == 'Angstrom':
                ax.xaxis.set_major_formatter(plt_ticker.FuncFormatter(lambda x, pos: f"{x * 10:g}"))
                xlabel = "Wavelength (Angstrom)"
            elif selected_x_unit == 'micron':
                ax.xaxis.set_major_formatter(plt_ticker.FuncFormatter(lambda x, pos: f"{x / 1000:g}"))
                xlabel = "Wavelength (micron)"
            # else: no formatter needed for nm
            ax.set_xlabel(xlabel)
            
            # 3. Y-Axis Label
            ax.set_ylabel("Normalized Flux" if is_norm_y else f"Flux (arbitrary units)")
            
            # --- X Limits ---
            # Priority: Autoscale if requested (e.g., by unit change)
            if was_zoomed and not force_autoscale:
                logging.debug(f"Restoring previous X zoom: {previous_xlim}")
                ax.set_xlim(previous_xlim)
            # Else (initial draw or zoomed out), let Matplotlib autoscale
            else:
                logging.debug("Allowing default X autoscale.")
                # No command needed, autoscale is default after clear()

            # --- Y Limits (Your simplified logic) ---
            if norm_state_changed:
                if is_norm_y:
                    logging.debug("Normalization toggled ON: Setting default Y limits.")
                    ax.set_ylim(-0.3, 1.3)
                else:
                    logging.debug("Normalization toggled OFF: Autoscaling Y.")
                    ax.autoscale(enable=True, axis='y')
            elif was_zoomed and not force_autoscale:
                logging.debug(f"Restoring previous Y zoom: {previous_ylim}")
                ax.set_ylim(previous_ylim)
            # Else (initial draw or zoomed out), let Matplotlib autoscale
            else:
                logging.debug("Allowing default Y autoscale.")

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

        # Update the limit boxes in the main window *after* plot is drawn
        if self.main_window:
            # Schedule this to run just after the draw completes
            QTimer.singleShot(0, self.main_window._update_limit_boxes_from_plot)
    
    def toggle_region_selector(self):
        """ Enables/disables the region selection mode. """
        # 1. Toggle state
        self._is_selecting_region = not self._is_selecting_region
        
        if self._is_selecting_region:
            logging.info("Region selection mode ACTIVE.")
            # Disable standard Matplotlib tools to avoid conflicts
            if self.toolbar.mode != '':
                self.toolbar.zoom() # Toggles off if active
                self.toolbar.pan()  # Toggles off if active
            
            # Change cursor to indicate mode
            self.setCursor(Qt.CrossCursor)
        else:
            logging.info("Region selection mode DEACTIVATED.")
            self.setCursor(Qt.ArrowCursor)
            # Clean up any half-drawn selection
            self.canvas.cleanup_selection()

    def get_cursor_line_positions_at_z(self, session_state: 'SessionV2', z_cursor):
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
    
    def calculate_z_from_x(self, mouse_x_nm: float) -> Optional[float]:
        """
        Helper to calculate redshift from a wavelength (nm) based on
        the currently active cursor series.
        """
        if not self.main_window or not V1_FUNCTIONS_AVAILABLE:
            return None
        try:
            series_str = self.main_window.cursor_series_input.text()
            transitions = trans_parse(series_str)
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
        Finds the data point nearest to the mouse event and shows a tooltip
        if it's close enough.
        """
        if not event or not event.xdata or not event.ydata:
            return
            
        # 1. Get current spectrum
        session = None
        if self.main_window and self.main_window.active_history:
             session = self.main_window.active_history.current_state
        
        if not session or not session.spec:
            return
            
        spec = session.spec
        x_full = spec.x.value
        y_full = spec.y.value

        # Calculate how many points are currently visible.
        xlim = self.canvas.axes.get_xlim()
        
        # Use searchsorted to quickly find the visible range indices
        idx_start_vis = np.searchsorted(x_full, xlim[0], side='left')
        idx_end_vis = np.searchsorted(x_full, xlim[1], side='right')
        n_visible = idx_end_vis - idx_start_vis
        
        # Threshold: if more than 5000 points are visible, it's too dense.
        # (You can adjust this number based on your preference)
        TOOLTIP_MAX_POINTS = 5000 
        if n_visible > TOOLTIP_MAX_POINTS:
            return

        # 2. Find nearest index in X (fast search on sorted array)
        # np.searchsorted finds the insertion point to maintain order
        idx = np.searchsorted(x_full, event.xdata)
        
        # Check the neighboring point to find the *true* nearest
        if idx > 0 and (idx == len(x_full) or 
                        np.abs(event.xdata - x_full[idx-1]) < np.abs(event.xdata - x_full[idx])):
            idx -= 1
            
        # Safety check for empty/out-of-bounds
        if idx < 0 or idx >= len(x_full):
            return

        # 3. Check proximity in PIXELS (crucial for usability)
        # Convert the data point back to screen coordinates
        data_point_screen = self.canvas.axes.transData.transform((x_full[idx], y_full[idx]))
        mouse_screen = (event.x, event.y)
        
        # Calculate euclidean distance in pixels
        dist_pix = np.hypot(data_point_screen[0] - mouse_screen[0], 
                            data_point_screen[1] - mouse_screen[1])
        
        THRESHOLD_PIX = 40 # Slightly tighter threshold for better separation
        
        # Helper to check distance
        def check_dist(y_array):
            if y_array is None: return False
            data_pos = self.canvas.axes.transData.transform((x_full[idx], y_array[idx]))
            mouse_pos = (event.x, event.y)
            return np.hypot(data_pos[0] - mouse_pos[0], data_pos[1] - mouse_pos[1]) < THRESHOLD_PIX

        # List of (label, value, unit) tuples to display
        rows_to_show = []

        # 1. Check Flux
        if check_dist(spec.y.value):
            rows_to_show.append(("y", spec.y.value[idx], spec.y.unit))
            rows_to_show.append(("dy", spec.dy.value[idx], spec.dy.unit))
            
        # 2. Check Continuum (if visible AND exists)
        if self.main_window.continuum_checkbox.isChecked():
            cont_col = spec._data.aux_cols.get('cont')
            if cont_col and check_dist(cont_col.values):
                rows_to_show.append(("cont", cont_col.values[idx], cont_col.unit))

        # 3. Check Model (if visible AND exists)
        if self.main_window.model_checkbox.isChecked():
            model_col = spec._data.aux_cols.get('model')
            if model_col and check_dist(model_col.values):
                rows_to_show.append(("model", model_col.values[idx], model_col.unit))

        # 4. Check Dynamic Aux Column
        selected_aux = self.main_window.aux_col_combo.currentText()
        if selected_aux and selected_aux != "None":
             # Avoid double-showing if it's already covered by the checkboxes
             is_redundant = (selected_aux == 'cont' and self.main_window.continuum_checkbox.isChecked()) or \
                            (selected_aux == 'model' and self.main_window.model_checkbox.isChecked())
             
             if not is_redundant:
                 aux_col = spec._data.aux_cols.get(selected_aux)
                 # We only show numerical columns in the tooltip for now
                 if aux_col and aux_col.values.dtype.kind in 'fiu' and check_dist(aux_col.values):
                      rows_to_show.append((selected_aux, aux_col.values[idx], aux_col.unit))
        
        # --- If no lines are near, don't show anything ---
        if not rows_to_show:
            return

        # 4. Build the HTML table
        try:
            tip_html = "<table style='font-size: 13px'>"
            
            # Always show X at the top if we are showing anything
            tip_html += f"<tr><td style='padding-right:10px'><b>x:</b></td><td>{x_full[idx]:.4f} {spec.x.unit}</td></tr>"
            
            # Add the dynamic rows
            for label, val, unit in rows_to_show:
                 tip_html += f"<tr><td style='padding-right:10px'><b>{label}:</b></td><td>{val:.4e} {unit}</td></tr>"

            tip_html += "</table>"

            QToolTip.showText(QCursor.pos(), tip_html, self)
                    
        except Exception as e:
            logging.warning(f"Error building tooltip: {e}")

    def update_plot(self, new_session_state: Optional['SessionV2'], force_autoscale: bool = False):
        """Called by MainWindowV2 to swap the immutable session object and redraw."""
        self.plot_spectrum(session_state=new_session_state, initial_draw=False, force_autoscale=force_autoscale)

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
        ('Select', 'Select region', 'toggle_select_mode', 'toggle_select_mode'),

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

        for action in self.actions():
            if action.text() == 'Select':
                # Use a standard Qt icon. 
                # SP_FileDialogDetailedView often looks like a list/selection.
                # You can try others like SP_ArrowRight or SP_ToolBarHorizontalExtensionButton
                icon = QApplication.style().standardIcon(QStyle.SP_FileDialogDetailedView)
                action.setIcon(icon)
                break
            
    def toggle_select_mode(self):
        """Toggles the region selection mode on/off."""
        if not self.viewer_parent: return
        
        # Toggle the mode in the parent widget
        self.viewer_parent.toggle_region_selector()
        
        # Visual feedback: "press" or "unpress" the button based on new state
        # (We need to find our action in the toolbar list)
        for action in self.actions():
             if action.text() == 'Select':
                 action.setChecked(self.viewer_parent._is_selecting_region)
                 break

    def home(self):
        """ Overrides the default 'home' to reset to the *current*
            data state, respecting normalization.
        """
        logging.debug("Custom 'Home' button pressed.")
        if self.viewer_parent and self.viewer_parent.main_window:
            session_state = None
            if self.viewer_parent.main_window.active_history:
                session_state = self.viewer_parent.main_window.active_history.current_state
            
            # Call our plot function with the new flag
            self.viewer_parent.plot_spectrum(
                session_state=session_state,
                force_autoscale=True
            )
        else:
            # Fallback to original behavior
            super().home()


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