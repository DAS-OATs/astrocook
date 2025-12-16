import traceback
from matplotlib import pylab
import astropy.constants as const
import astropy.units as au
import json
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
import qtawesome as qta
import scienceplots
from scipy.interpolate import CubicSpline
from typing import Optional, TYPE_CHECKING

from astrocook.core.atomic_data import STANDARD_MULTIPLETS, xem_d, is_hydrogen_line
from astrocook.core.photometry import STANDARD_FILTERS, get_filter_transmission

# Use TYPE_CHECKING to avoid circular import errors at runtime
if TYPE_CHECKING:
    from astrocook.core.session import SessionV2 

try:
    from astrocook.legacy.functions import trans_parse, x_convert
    V1_FUNCTIONS_AVAILABLE = True
except ImportError:
    V1_FUNCTIONS_AVAILABLE = False
    trans_parse = lambda s: [s] # Dummy function
    logging.error("V1 'trans_parse' not found. Displaying V1-style systems may fail.")

# Keys must match exactly what is in your legacy.vars.xem_d dictionary
STRONG_EMISSION_LINES = {
    'Ly_a', 'Ly_b', 'NV_1238', 'SiIV_1393', 'CIV_1548', 
    'CIII_1908', 'MgII_2796', 'H_b', 'OIII_5008', 'H_a'
}

# --- Centralized Plot Style Configuration ---
PLOT_STYLE = {
    'flux': {
        'color_idx': 0,       # Index in color cycle (tab20)
        'lw': 0.5,
        'alpha': 1.0,
        'label': 'Flux'
    },
    'model': {
        'color_idx': 1,       # Index 2 is usually Orange/Red in tab20
        'ls': '-',
        'lw': 0.8,
        'alpha': 1.0,
        'label': 'Model'
    },
    'error': {
        'color': '#aaaaaa',   # Grey shading
        'lw': 0.3,
        'alpha': 1.0
    },
    'continuum': {
        'color': 'black',
        'ls': '--',
        'lw': 0.8
    },
    'cursor': {
        'color_idx': 3,
        'ls': '--',
        'lw': 1.0,
        'alpha': 0.8
    },
    'zero_line': {
        'color': 'gray',
        'ls': '-',
        'lw': 0.5,
        'alpha': 1.0
    },
    'center_line': {
        'color': 'gray',
        'ls': '-',
        'lw': 0.5,
        'alpha': 1.0
    }
}

def get_style_color(key: str, cycle: list):
    """Helper to get color from style dict or cycle."""
    if 'color' in PLOT_STYLE[key]:
        return PLOT_STYLE[key]['color']
    else:
        idx = PLOT_STYLE[key].get('color_idx', 0)
        return cycle[idx % len(cycle)]

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
        cmap = plt.get_cmap(fallback_cmap)
        return [cmap(i / n) for i in range(n)]

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

        fig = Figure(figsize=(width, height), dpi=dpi, constrained_layout=True)
        self.axes = fig.add_subplot(111)
        #fig.tight_layout()

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
                self.selection_start_x, event.xdata, color='yellow', alpha=0.3
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
        self.plot_widget.on_mouse_move(event)


    def _on_hover_timeout(self):
        """Called when the mouse has stopped moving for a while."""
        if self.plot_widget and self._last_mouse_event and self._last_mouse_event.inaxes == self.axes:
            self.plot_widget.show_data_tooltip(self._last_mouse_event)

    def on_press(self, event):
        """Delegate press to widget."""
        if self.plot_widget:
            self.plot_widget.on_press(event)

    def on_release(self, event):
        """Delegate release to widget."""
        if self.plot_widget:
            self.plot_widget.on_release(event)

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
        
        # This style applies to any QToolTip triggered from this widget
        self.setStyleSheet("""
            QToolTip {
                border: 1px solid gray;
                border-radius: 5px;     /* The rounded corners */
                background-color: #F2F8F8;
                color: black;
                padding: 4px;           /* internal padding */
                background-clip: border-box;
            }
        """)

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
        self._cached_norm_plot_raw = None
        self._cached_model_plot_raw = None

        # Cache for NORMALIZED plot-ready (decimated) arrays
        self._cached_y_plot_norm = None
        self._cached_dy_plot_norm = None
        self._cached_norm_plot_norm = None
        self._cached_model_plot_norm = None

        # Cache for SNR plot-ready (decimated) arrays
        self._cached_y_plot_snr = None
        self._cached_dy_plot_snr = None
        self._cached_model_plot_snr = None # (e.g., model / dy)
        self._cached_norm_plot_snr = None # (e.g., norm / dy)
        self._cached_snr_error_col = "" # Tracks which error col was used

        # List to store active iso-mag curves
        self._show_isomag_grid = False
        self.isomag_artists = []

        # Selection Mode State
        self._is_selecting_region = False
        self._selection_start_x = None
        
        # Store the original Matplotlib coordinate formatter ONCE
        self._default_format_coord = self.canvas.axes.format_coord

        # Tooltip state
        self._tip_component_html = None
        self._tip_data_rows = None # List of (label, value, unit) tuples
        self._tip_region_html = None

        # Add ONLY toolbar to main layout
        self.main_layout.addWidget(self.toolbar, 0) # Stretch 0

        # --- Continuum Editor State ---
        self._edit_mode_active = False
        self._knots_x = None  # np.array
        self._knots_y = None  # np.array
        self._active_knot_idx = None # Index of knot being dragged
        
        self.knot_artist = None  # Matplotlib Line2D (markers)
        self.spline_artist = None # Matplotlib Line2D (smooth line)

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

        selected_snr_col = self.main_window.snr_col_combo.currentText()
        
        # [CHANGE] Map "dF" display name back to internal "dy" name
        if selected_snr_col == "dF":
            selected_snr_col = "dy"

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
            self._cached_norm_plot_raw = None
            self._cached_model_plot_raw = None
            self._cached_y_plot_norm = None
            self._cached_dy_plot_norm = None
            self._cached_norm_plot_norm = None
            self._cached_model_plot_norm = None
            self._cached_y_plot_snr = None
            self._cached_dy_plot_snr = None
            self._cached_norm_plot_snr = None
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
        
        self._tick_visuals = []

        plot_occurred = False
        if spec and len(spec.x) > 0:
            plot_occurred = True

            # --- ** Data Slicing Logic ** ---
            # Get full data arrays
            full_x_data = spec.x.value
            full_y_data = spec.y.value
            full_dx_data = spec.xmax.value-spec.xmin.value if (spec.xmin, spec.xmax) is not (None, None) else np.zeros_like(full_x_data)
            full_dy_data = spec.dy.value if spec.dy is not None else np.zeros_like(full_y_data)
            full_norm_data = spec.norm.value if hasattr(spec.norm, 'value') else spec.norm
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
                    dx_data = self._cached_dx_plot_norm
                    dy_data = self._cached_dy_plot_norm
                    norm_data = self._cached_norm_plot_norm
                    model_data = self._cached_model_plot_norm
                
                else:
                    logging.debug("Cache miss for normalized data. Computing...")
                    # 2a. Need to COMPUTE. First, get RAW data.
                    if self._cached_x_plot_raw is not None:
                        logging.debug("...using cached RAW plot data.")
                        x_data_raw = self._cached_x_plot_raw
                        y_data_raw = self._cached_y_plot_raw
                        dx_data_raw = self._cached_dx_plot_raw
                        dy_data_raw = self._cached_dy_plot_raw
                        norm_data_raw = self._cached_norm_plot_raw
                        model_data_raw = self._cached_model_plot_raw
                    else:
                        logging.debug("...RAW plot data cache is also empty. Computing...")
                        # 3a. TOTAL cache miss. Compute RAW data.
                        if use_decimation:
                            x_data_raw = decimate_x_for_min_max(full_x_data, DECIMATION_FACTOR)
                            y_data_raw = decimate_y_min_max(full_y_data, DECIMATION_FACTOR)
                            dx_data_raw = decimate_y_min_max(full_dx_data, DECIMATION_FACTOR)
                            dy_data_raw = decimate_y_min_max(full_dy_data, DECIMATION_FACTOR)
                            norm_data_raw = decimate_y_min_max(full_norm_data, DECIMATION_FACTOR)
                            model_data_raw = decimate_y_min_max(full_model_data, DECIMATION_FACTOR)
                        else: # We are zoomed, forcing home, or data is small
                            x_data_raw = full_x_data
                            y_data_raw = full_y_data
                            dx_data_raw = full_dx_data
                            dy_data_raw = full_dy_data
                            norm_data_raw = full_norm_data
                            model_data_raw = full_model_data
                        
                        # Save to RAW cache
                        self._cached_x_plot_raw = x_data_raw
                        self._cached_y_plot_raw = y_data_raw
                        self._cached_dx_plot_raw = dx_data_raw
                        self._cached_dy_plot_raw = dy_data_raw
                        self._cached_norm_plot_raw = norm_data_raw
                        self._cached_model_plot_raw = model_data_raw

                    # 4a. NOW, compute and cache the NORMALIZED data
                    if norm_data_raw is not None:
                        x_data = x_data_raw
                        norm_val = norm_data_raw
                        y_data = np.divide(y_data_raw, norm_val, out=np.full_like(y_data_raw, np.nan), where=norm_val!=0)
                        dx_data = np.divide(dx_data_raw, norm_val, out=np.full_like(dx_data_raw, np.nan), where=norm_val!=0)
                        dy_data = np.divide(dy_data_raw, norm_val, out=np.full_like(dy_data_raw, np.nan), where=norm_val!=0)
                        model_data = None
                        if model_data_raw is not None:
                            model_data = np.divide(model_data_raw, norm_val, out=np.full_like(model_data_raw, np.nan), where=norm_val!=0)
                        norm_data = np.ones_like(norm_data_raw)
                        
                        # Save to NORMALIZED cache
                        self._cached_x_plot_norm = x_data
                        self._cached_y_plot_norm = y_data
                        self._cached_dx_plot_norm = dx_data
                        self._cached_dy_plot_norm = dy_data
                        self._cached_norm_plot_norm = norm_data
                        self._cached_model_plot_norm = model_data
                    else:
                        logging.warning("Normalize checked, but continuum data is missing.")
                        x_data, y_data, dx_data, dy_data, norm_data, model_data = \
                            x_data_raw, y_data_raw, dx_data_raw, dy_data_raw, norm_data_raw, model_data_raw

            elif is_snr:
                # --- *** NEW BLOCK FOR SNR *** ---
                if self._cached_y_plot_snr is not None:
                    logging.debug("Using cached SNR plot data.")
                    x_data = self._cached_x_plot_raw
                    y_data = self._cached_y_plot_snr
                    dx_data = self._cached_dx_plot_snr
                    dy_data = self._cached_dy_plot_snr
                    norm_data = self._cached_norm_plot_snr
                    model_data = self._cached_model_plot_snr
                
                else:
                    logging.debug("Cache miss for SNR data. Computing...")
                    # 1. Get RAW plot data (from cache or compute)
                    if self._cached_x_plot_raw is not None:
                        # ... (read from _cached_x_plot_raw, etc.) ...
                        x_data_raw = self._cached_x_plot_raw
                        y_data_raw = self._cached_y_plot_raw
                        dx_data_raw = self._cached_dx_plot_raw
                        dy_data_raw = self._cached_dy_plot_raw
                        norm_data_raw = self._cached_norm_plot_raw
                        model_data_raw = self._cached_model_plot_raw
                    else:
                        # ... (compute raw decimated data) ...
                        if use_decimation:
                            # ... (decimate_x_for_min_max, decimate_y_min_max) ...
                            x_data_raw = decimate_x_for_min_max(full_x_data, DECIMATION_FACTOR)
                            y_data_raw = decimate_y_min_max(full_y_data, DECIMATION_FACTOR)
                            dx_data_raw = decimate_y_min_max(full_dx_data, DECIMATION_FACTOR)
                            dy_data_raw = decimate_y_min_max(full_dy_data, DECIMATION_FACTOR)
                            norm_data_raw = decimate_y_min_max(full_norm_data, DECIMATION_FACTOR)
                            model_data_raw = decimate_y_min_max(full_model_data, DECIMATION_FACTOR)
                        else:
                            # ... (get from data_slice) ...
                            x_data_raw = full_x_data[data_slice]
                            y_data_raw = full_y_data[data_slice]
                            dx_data_raw = full_dx_data[data_slice]
                            dy_data_raw = full_dy_data[data_slice]
                            norm_data_raw = full_norm_data[data_slice] if full_norm_data is not None else None
                            model_data_raw = full_model_data[data_slice] if full_model_data is not None else None
                        
                        # Save to RAW cache
                        self._cached_x_plot_raw = x_data_raw
                        self._cached_y_plot_raw = y_data_raw
                        self._cached_dx_plot_raw = dx_data_raw
                        self._cached_dy_plot_raw = dy_data_raw
                        self._cached_norm_plot_raw = norm_data_raw
                        self._cached_model_plot_raw = model_data_raw

                    # 2. Get the ERROR array to use
                    error_col_name = self.main_window.snr_col_combo.currentText()
                    full_error_data = None
                    try:
                        if error_col_name == 'dy':
                            full_error_data = full_dy_data
                        else:
                            # This path is for 'running_std', 'norm_err', etc.
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
                    
                    norm_data = None
                    if norm_data_raw is not None:
                        norm_data = np.divide(norm_data_raw, error_data_raw, out=np.full_like(norm_data_raw, np.nan), where=error_data_raw!=0)
                    
                    model_data = None
                    if model_data_raw is not None:
                        model_data = np.divide(model_data_raw, error_data_raw, out=np.full_like(model_data_raw, np.nan), where=error_data_raw!=0)

                    # Save to SNR cache
                    x_data = x_data_raw # X is always the same
                    dx_data = dx_data_raw # X is always the same
                    self._cached_y_plot_snr = y_data
                    self._cached_dx_plot_snr = dx_data
                    self._cached_dy_plot_snr = dy_data
                    self._cached_norm_plot_snr = norm_data
                    self._cached_model_plot_snr = model_data
                    self._cached_snr_error_col = error_col_name # Store the col used

            else:
                # 1b. We want RAW data. Check cache.
                if self._cached_x_plot_raw is not None:
                    logging.debug("Using cached RAW plot data.")
                    x_data = self._cached_x_plot_raw
                    y_data = self._cached_y_plot_raw
                    dx_data = self._cached_dx_plot_raw
                    dy_data = self._cached_dy_plot_raw
                    norm_data = self._cached_norm_plot_raw
                    model_data = self._cached_model_plot_raw
                
                else:
                    logging.debug("Cache miss for raw plot data. Computing...")
                    # 2b. Compute RAW data.
                    if use_decimation:
                        x_data = decimate_x_for_min_max(full_x_data, DECIMATION_FACTOR)
                        y_data = decimate_y_min_max(full_y_data, DECIMATION_FACTOR)
                        dx_data = decimate_y_min_max(full_dx_data, DECIMATION_FACTOR)
                        dy_data = decimate_y_min_max(full_dy_data, DECIMATION_FACTOR)
                        norm_data = decimate_y_min_max(full_norm_data, DECIMATION_FACTOR)
                        model_data = decimate_y_min_max(full_model_data, DECIMATION_FACTOR)
                    else: # We are zoomed, forcing home, or data is small
                        x_data = full_x_data
                        y_data = full_y_data
                        dx_data = full_dx_data
                        dy_data = full_dy_data
                        norm_data = full_norm_data
                        model_data = full_model_data
                    
                    # Save to RAW cache
                    self._cached_x_plot_raw = x_data
                    self._cached_y_plot_raw = y_data
                    self._cached_dx_plot_raw = dx_data
                    self._cached_dy_plot_raw = dy_data
                    self._cached_norm_plot_raw = norm_data
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
                        dx_data = dx_data[final_slice]
                        dy_data = dy_data[final_slice]
                        if norm_data is not None:
                            norm_data = norm_data[final_slice]
                        if model_data is not None:
                            model_data = model_data[final_slice]
                    else:
                        logging.debug("Zoom slice is empty, plotting full data.")
                except Exception as e:
                    logging.warning(f"Failed to apply zoom slice: {e}")

            x_unit = str(spec.x.unit)
                    
            colors = get_color_cycle(5, cmap='tab20') # Get colors
            #colors[0] = "#296bff"

            plot_style = spec.meta.get('plot_style', 'step')
            
            is_scatter = False
            n_points = len(x_data)
            if n_points < 20000:
                is_scatter = (plot_style == 'scatter')
                scatter_ms = 4.0
                if n_points < 600:
                    scatter_ms = 10.0
                elif n_points < 2000:
                    scatter_ms = 8.0
                elif n_points < 6000:
                    scatter_ms = 6.0

            # --- Check state from main_window reference ---
            # 2. Plot Error Shading (Conditional)
            if self.main_window.error_checkbox.isChecked(): # <<< Check main window's checkbox
                #ax.fill_between(
                #    x_data, y_data - dy_data, y_data + dy_data,
                #    step='mid', color='#aaaaaa', alpha=0.5,
                #    label='1-sigma error', rasterized=True
                #)
                if dx_data is not None and dy_data is not None:
                    label = '1-sigma error'
                    if is_scatter:
                        # fmt='none' ensures no connecting lines or extra markers
                        ax.errorbar(x_data, y_data, xerr=dx_data, 
                                    fmt='none', ecolor='#aaaaaa', elinewidth=0.8, 
                                    alpha=0.5, rasterized=True, zorder=0.9, label='bin size')
                        ax.errorbar(x_data, y_data, yerr=dy_data, 
                                    fmt='none', ecolor='#aaaaaa', elinewidth=0.8, 
                                    alpha=0.5, rasterized=True, zorder=0.9, label=label)
                    elif use_decimation:
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
            if is_scatter:
                # Use a fast marker plot for echelle data
                # 'ls="None"' ensures no connecting lines (sawtooth removal)
                # 'marker="."' is a pixel-like dot
                ax.plot(x_data, y_data, label="Spectrum", 
                        linestyle='None', marker='.', markersize=scatter_ms, markeredgewidth=0,
                        color=colors[0], alpha=0.7, rasterized=True)
            elif use_decimation:
                ax.plot(x_data, y_data, label="Spectrum (decimated)", lw=0.5, color=colors[0], rasterized=True)
            else:
                ax.step(x_data, y_data, where='mid', label="Spectrum", lw=0.5, color=colors[0], rasterized=True)

            # 3. Plot Continuum (Conditional)
            if self.main_window.continuum_checkbox.isChecked(): # <<< Check main window's checkbox
                if norm_data is not None:
                    ax.plot(
                        x_data, norm_data, linestyle='--',
                        color='black', lw=0.8, label='Normalization', rasterized=True
                    )
                else:
                    logging.warning("Normalization requested but spec.norm is None.")

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
                    
                    if aux_col_to_plot in ['abs_mask', 'abs_ids']:
                        # Create a boolean mask where regions exist
                        region_mask = (aux_data > 0)
                        
                        # Use a blended transform to shade the full vertical extent
                        # X in data coordinates, Y in axes coordinates (0 to 1)
                        trans = ax.get_xaxis_transform()
                        
                        # Fill with a subtle color (e.g., teal/cyan with low alpha)
                        ax.fill_between(x_data, 0.05, 0.95, # 0 to 1 in axes coords covers full height
                                        where=region_mask,
                                        color='teal', alpha=0.05, 
                                        transform=trans,
                                        label='Detected Regions',
                                        rasterized=True)

                    # Plot based on data type
                    elif aux_data.dtype.kind in 'fiu': # Float or Int
                        if is_scatter:
                             ax.plot(x_data, aux_data, label=aux_col_to_plot, 
                                     linestyle='None', marker=scatter_ms, markersize=1, 
                                     color='purple', rasterized=True)
                        elif use_decimation:
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
                    
                    # V1/V2 Fallback Logic for Transitions
                    transitions = []
                    try:
                        # 1. Try V2-native lookup (e.g., "CIV")
                        transitions = STANDARD_MULTIPLETS[comp.series]
                    except KeyError:
                        if V1_FUNCTIONS_AVAILABLE:
                            # 2. Try V1 trans_parse (e.g., "CIV_1548,CIV_1550")
                            try:
                                transitions = trans_parse(comp.series)
                            except Exception:
                                pass # Will be caught by logging below
                        
                        # 3. Fallback for single lines (e.g., "Ly_a")
                        if not transitions and comp.series in xem_d:
                            transitions = [comp.series]
                    
                    if not transitions:
                        logging.warning(f"Could not parse series '{comp.series}' for system plot.")

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

                                self._tick_visuals.append((ax, x_plot, comp))
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
            ax.set_title(f"Spectrum: {session_state.name}")
            
            # 1. Scales (Log/Linear)
            ax.set_xscale('log' if is_log_x else 'linear')
            ax.set_yscale('log' if is_log_y else 'linear')

            ax.grid(True, linestyle=':')
            
            # 2. X-Axis Formatting
            #xlabel = "Wavelength (nm)" # Default
            xaxis_formatter = None  # Store the formatter to re-use it

            if selected_x_unit == 'Angstrom':
                xaxis_formatter = plt_ticker.FuncFormatter(lambda x, pos: f"{x * 10:g}")
                ax.xaxis.set_major_formatter(xaxis_formatter)
                xlabel = "Wavelength (Angstrom)"
            elif selected_x_unit == 'micron':
                xaxis_formatter = plt_ticker.FuncFormatter(lambda x, pos: f"{x / 1000:g}")
                ax.xaxis.set_major_formatter(xaxis_formatter)
                xlabel = "Wavelength (micron)"
            # else: no formatter needed for nm
            unit_str = au.Unit(selected_x_unit).to_string('latex_inline')
            xlabel = f"Wavelength ({unit_str})"
            ax.set_xlabel(xlabel)
            
            # Secondary Axis & Hover

            z_em = spec._data.z_em
            
            # 1. Define the custom hover-text formatter
            # Rimpiazziamo completamente il formatter di default per avere il formato:
            # Wave=... nm, Flux=... erg/..., Rest=... nm
            
            # Symbols
            SYM_LAMBDA = "\u03bb"      # λ
            SYM_LAMBDA_R = "\u03bb\u1d63" # λᵣ (lambda sub r)
            SYM_FLUX = "F"        # 
            SYM_FLUX_NORM = "\u0192"        # ƒ (hook f for flux/function)
            SYM_ANGSTROM = "\u00c5"    # Å
            SYM_MICRON = "\u00b5m"     # µm

            # Helper per le unità X
            x_label_unit = "nm"
            if selected_x_unit == 'Angstrom': x_label_unit = SYM_ANGSTROM
            elif selected_x_unit == 'micron': x_label_unit = SYM_MICRON
            if is_norm_y: y_label_symbol = SYM_FLUX_NORM
            else: y_label_symbol = SYM_FLUX
            
            # Helper per le unità Y (pulizia stringa Astropy)
            y_label_unit = ""
            if is_norm_y:
                y_label_unit = "" 
            elif is_snr:
                y_label_unit = "" 
            elif spec.y.unit:
                # Convertiamo l'unità astropy in stringa
                # Sostituiamo le notazioni "brutte" con caratteri Unicode se possibile
                unit_str = f"{spec.y.unit:s}"
                unit_str = unit_str.replace("Angstrom", SYM_ANGSTROM)
                unit_str = unit_str.replace("**", "^") # Es. cm^2
                y_label_unit = unit_str

            def new_format_coord(x, y):
                # x è in nm (unità nativa dei dati)
                # y è il valore visualizzato (flusso, norm o snr)
                
                # 1. Calcolo X Visualizzato
                x_disp = x
                if selected_x_unit == 'Angstrom': x_disp = x * 10.0
                elif selected_x_unit == 'micron': x_disp = x / 1000.0
                
                # 2. Costruzione Stringa Base
                if is_norm_y or is_snr:
                    y_fmt = f"{y:.3f}"
                else:
                    y_fmt = f"{y:.3e}"
                
                # Formato: λ=1215.67 Å, ƒ=1.23e-17 ...
                coord_str = f"{SYM_LAMBDA} = {x_disp:.4f} {x_label_unit}, {y_label_symbol} = {y_fmt} {y_label_unit}"
                
                # 3. Aggiunta Rest Frame
                if z_em > 0.0 and spec.x.unit.is_equivalent(au.nm):
                    x_rf_nm = x / (1 + z_em)
                    
                    if selected_x_unit == 'Angstrom': x_rf_disp = x_rf_nm * 10.0
                    elif selected_x_unit == 'micron': x_rf_disp = x_rf_nm / 1000.0
                    else: x_rf_disp = x_rf_nm
                        
                    coord_str += f", {SYM_LAMBDA_R} = {x_rf_disp:.4f} {x_label_unit}"
                
                return coord_str

            ax.format_coord = new_format_coord

            # 2. Add Secondary Rest-Frame Axis (if applicable)
            if z_em > 0.0 and spec.x.unit.is_equivalent(au.nm):
                
                def obs_to_rf(x): return x / (1 + z_em)
                def rf_to_obs(x): return x * (1 + z_em)
                
                sec_ax = ax.secondary_xaxis('top', functions=(obs_to_rf, rf_to_obs))
                sec_ax.set_xlabel(f"Rest-frame ({SYM_LAMBDA_R}, {x_label_unit})")
                
                # Apply the unit formatter (or reset to default scalar if None)
                if xaxis_formatter:
                    sec_ax.xaxis.set_major_formatter(xaxis_formatter)
                else:
                    # Explicitly reset to default if we switched back to 'nm'
                    sec_ax.xaxis.set_major_locator(plt_ticker.AutoLocator())

            # 3. Emission Line Labels
            if z_em > 0.0 and V1_FUNCTIONS_AVAILABLE and spec.x.unit.is_equivalent(au.nm):
                # Get current view limits (native units, usually nm)
                # We use the limits *we are about to set* if we know them, 
                # otherwise we rely on what Matplotlib will autoscale to.
                # For simplicity, let's use the full data range if autoscaling,
                # or the zoomed range if known.
                if force_autoscale:
                     view_xmin, view_xmax = full_x_data[0], full_x_data[-1]
                elif was_zoomed:
                     view_xmin, view_xmax = previous_xlim
                else:
                     view_xmin, view_xmax = full_x_data[0], full_x_data[-1]

                # Transformation: X is data coords, Y is axes coords (0=bottom, 1=top)
                trans = ax.get_xaxis_transform()
                
                strong_only = self.main_window.strong_lines_checkbox.isChecked()

                # Loop through standard emission lines
                for label, rest_wave_q in xem_d.items():
                    if strong_only and label not in STRONG_EMISSION_LINES:
                        continue

                    if label not in STRONG_EMISSION_LINES:
                        color = 'gray'
                        alpha_text = 0.5
                        alpha_line = 0.15
                    else:
                        color = 'black'
                        alpha_text = 1.0
                        alpha_line = 0.3

                    # Calculate observed wavelength in native units (nm)
                    try:
                        rest_nm = rest_wave_q.to_value(au.nm)
                        obs_nm = rest_nm * (1.0 + z_em)
                        
                        # Only plot if within the current view
                        if view_xmin <= obs_nm <= view_xmax:
                            ax.text(obs_nm, 0.95, label, 
                                    transform=trans,
                                    rotation=90,
                                    verticalalignment='top',
                                    horizontalalignment='center',
                                    color=color, alpha=alpha_text,
                                    fontsize=8, clip_on=True)
                            # Optional: Add a faint vertical line
                            ax.axvline(obs_nm, color=color, alpha=alpha_line, lw=1.0, ls=':')
                            
                    except Exception:
                        continue # Skip lines with incompatible units (rare)

            # 3. X/Y-Axis Label
            ax.set_xlabel(f"Wavelength ({SYM_LAMBDA}, {x_label_unit})")
            
            if is_norm_y:
                y_label_text = f"Normalized Flux ({SYM_FLUX_NORM})"
            elif is_snr:
                y_label_text = f"Signal-to-Noise Ratio ({SYM_FLUX} / {selected_snr_col})"
            else:
                # Check if we have a real physical unit
                if spec.y.unit is None or spec.y.unit == au.dimensionless_unscaled:
                     y_label_text = f"Flux ({SYM_FLUX}, arbitrary units)"
                else:
                     # Use Astropy's pretty formatter for the unit string
                     unit_str = spec.y.unit.to_string('latex_inline')
                     y_label_text = f"Flux ({SYM_FLUX}, {unit_str})"
                
            ax.set_ylabel(y_label_text)
            
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

        # Draw the iso-mag grid last. It will self-manage its impact on limits.
        self._draw_isomag_grid()

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
            # V1/V2 Fallback Logic for Transitions
            series_str = main_window.cursor_series_input.text()
            transitions = []
            try:
                # 1. Try V2-native lookup
                transitions = STANDARD_MULTIPLETS[series_str]
            except KeyError:
                if V1_FUNCTIONS_AVAILABLE:
                    # 2. Try V1 trans_parse
                    try:
                        transitions = trans_parse(series_str)
                    except Exception:
                        pass
                        
                # 3. Fallback for single lines
                if not transitions and series_str in xem_d:
                    transitions = [series_str]

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
            # V1/V2 Fallback Logic for Transitions
            series_str = self.main_window.cursor_series_input.text()
            transitions = []
            try:
                # 1. Try V2-native lookup
                transitions = STANDARD_MULTIPLETS[series_str]
            except KeyError:
                if V1_FUNCTIONS_AVAILABLE:
                    # 2. Try V1 trans_parse
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
        
    def _draw_isomag_grid(self):
        """
        Helper to draw the iso-magnitude grid.
        Crucially, it saves and restores axis limits so these lines
        do NOT affect autoscaling or the 'Home' view.
        """
        # 1. Clear old artists first
        for artist in self.isomag_artists:
            try: artist.remove()
            except: pass
        self.isomag_artists = []

        if not self._show_isomag_grid:
            return

        # 2. Get data and perform checks
        if not self.main_window.active_history: return
        spec = self.main_window.active_history.current_state.spec
        if not spec: return

        standard_flux = au.erg / (au.cm**2 * au.s * au.Angstrom)
        if spec.y.unit is None or not spec.y.unit.is_equivalent(standard_flux):
             return

        # 3. SAVE CURRENT LIMITS (The secret to acting like a passive grid)
        current_xlim = self.canvas.axes.get_xlim()
        current_ylim = self.canvas.axes.get_ylim()

        # 4. Calculate and Plot
        try:
            x_aa = spec.x.to(au.Angstrom).value
            x_plot = spec.x.value
            current_y_unit = spec.y.unit
            
            c_aa_s = const.c.to(au.Angstrom/au.s).value
            jy_to_cgs_fnu = (1.0 * au.Jy).to(au.erg / (au.cm**2 * au.s * au.Hz)).value
            
            grid_mags = [15, 16, 17, 18, 19, 20, 21]
            
            x_label_pos = min(current_xlim[1], x_plot[-1])

            for mag in grid_mags:
                # Calculate curve
                f_nu_val = 3631.0 * (10**(-0.4 * mag)) * jy_to_cgs_fnu
                f_lambda_val = f_nu_val * c_aa_s / (x_aa**2)
                y_plot_iso = (f_lambda_val * standard_flux).to(current_y_unit).value

                # Plot curve
                artist, = self.canvas.axes.plot(x_plot, y_plot_iso, 
                                                ls=':', lw=1.0, color='orange', alpha=0.4,
                                                zorder=0.1) 
                self.isomag_artists.append(artist)
                
                # --- NEW: Smart Labeling ---
                # Interpolate to find Y at our chosen X position
                # We can safely use the last data point as 'right' fill value because
                # we ensured x_label_pos <= x_plot[-1]
                y_label_pos = np.interp(x_label_pos, x_plot, y_plot_iso)
                
                # Only draw label if it's vertically within the current view
                if current_ylim[0] <= y_label_pos <= current_ylim[1]:
                    # Use annotate to place text exactly at (x_label_pos, y_label_pos)
                    # and then offset it slightly to the right (10 points) so it doesn't overlap.
                    lbl = self.canvas.axes.annotate(
                        f"m={mag}", 
                        xy=(x_label_pos, y_label_pos), 
                        xytext=(5, 0),                 # 5 points to the right
                        textcoords="offset points", 
                        color='orange', fontsize=9, alpha=0.8,
                        verticalalignment='center',
                        clip_on=False # Allow it to spill into the margin if needed
                    )
                    self.isomag_artists.append(lbl)
                # ----------------------------------------------------

        except Exception as e:
            logging.error(f"Error drawing iso-mag grid: {e}")

        # 5. RESTORE LIMITS
        self.canvas.axes.set_xlim(current_xlim)
        self.canvas.axes.set_ylim(current_ylim)   

    def toggle_isomag_grid(self, visible: bool):
        """ Toggles the state and triggers a redraw. """
        self._show_isomag_grid = visible
        # We call plot_spectrum to ensure everything (including limits) 
        # is in a consistent state.
        self.plot_spectrum(
            session_state=self.main_window.active_history.current_state,
            force_autoscale=False 
        )
    
    def show_data_tooltip(self, event):
        """
        Slow handler called by timer (hover stop).
        Calculates Component, Region, and Data info and displays the unified tooltip.
        """
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
            c_resol = f"{found_comp.resol:.0f}" if found_comp.resol is not None else "N/A"
            
            self._tip_component_html = (
                f"<div style='font-size: 13px'>"
                f"<b>ID:</b> {found_comp.id} | <b>{found_comp.series}</b><br>"
                f"<b>z:</b> {found_comp.z:.5f}<br>"
                f"<b>{u'\u03c7'}\u00b2:</b> {c_chi2} | <b>R:</b> {c_resol}"
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
                region_id_val = int(spec._data.aux_cols['abs_ids'].values[idx])
                if region_id_val != 0:
                    ident_labels_json = spec.meta.get('region_identifications')
                    if ident_labels_json:
                        ident_labels_raw = json.loads(ident_labels_json)
                        ident_labels_dict = {int(k): v for k, v in ident_labels_raw.items()}
                        
                        if region_id_val in ident_labels_dict:
                            ids_for_region = ident_labels_dict[region_id_val]
                            labels = [f"{name} ({score:.2f})" if score > 0 else name for name, score in ids_for_region]
                            label_str = ", ".join(labels)
                            if label_str:
                                self._tip_region_html = f"<tr><td style='padding-right:8px'><b>ID:</b></td><td>{label_str}</td></tr>"
            except Exception: pass

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

    def update_plot(self, new_session_state: Optional['SessionV2'], force_autoscale: bool = False):
        """Called by MainWindowV2 to swap the immutable session object and redraw."""
        self.plot_spectrum(session_state=new_session_state, initial_draw=False, force_autoscale=force_autoscale)

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
        """
        Handler called on mouse move.
        We do nothing here to prevent flickering. All tooltip logic is deferred 
        to show_data_tooltip (called when the mouse stops).
        """
        if self._edit_mode_active and self._active_knot_idx is not None and event.inaxes:
            # Vertical Only: Update Y, keep X original
            self._knots_y[self._active_knot_idx] = event.ydata
            
            # Update visuals
            self.update_spline_preview()
            
            # Reuse the canvas's built-in animator. 
            # It restores the background -> draws all cursor_artists (including ours) -> blits.
            self.canvas.draw_animated()
            
    # Add on_press for Context Menu
    def on_press(self, event):
        if event.button == 1 and event.inaxes == self.canvas.axes and self._is_selecting_region:
            self._selection_start_x = event.xdata
            self.canvas.selection_artist = self.canvas.axes.axvspan(event.xdata, event.xdata, color='orange', alpha=0.3)
            self.canvas.draw_idle() 
            return True

        if event.button == 3 and event.inaxes:
            found_comp = None
            xlim = event.inaxes.get_xlim(); tol = (xlim[1] - xlim[0]) * 0.015 
            for ax, x_pos, comp in self._tick_visuals:
                if ax == event.inaxes and abs(event.xdata - x_pos) < tol: found_comp = comp; break
            
            menu = QMenu(self)
            if found_comp:
                act_inspect = QAction(f"Inspect {found_comp.series} (Group View)", menu)
                if self.main_window: act_inspect.triggered.connect(lambda checked=False, u=found_comp.uuid: self.main_window.open_inspector_on_component(u))
                menu.addAction(act_inspect)
                menu.addSeparator()

            z_click = self.calculate_z_from_x(event.xdata)
            if z_click is not None and self.main_window and self.main_window.cursor_show_checkbox.isChecked():
                series_str = self.main_window.cursor_series_input.text()
                act_add = QAction(f"Add {series_str} at z={z_click:.5f}", menu)
                act_add.triggered.connect(lambda checked=False, z=z_click, s=series_str: self.main_window._on_recipe_requested("absorbers", "add_component", {'series': s, 'z': z}, {}))
                menu.addAction(act_add)
                
                act_zem = QAction(f"Set emission redshift to z={z_click:.5f}", menu)
                act_zem.triggered.connect(lambda checked=False, z=z_click: self.main_window._on_recipe_requested("edit", "set_properties", {"z_em": str(z)}, {}))
                menu.addAction(act_zem)

            if not menu.isEmpty(): menu.exec(QCursor.pos()); return True

        if self._edit_mode_active and event.inaxes == self.canvas.axes:
            click_x, click_y = event.xdata, event.ydata
            
            # Calculate distance in pixels
            knot_pixels = self.canvas.axes.transData.transform(np.column_stack([self._knots_x, self._knots_y]))
            click_pixel = self.canvas.axes.transData.transform((click_x, click_y))
            dists = np.hypot(knot_pixels[:,0] - click_pixel[0], knot_pixels[:,1] - click_pixel[1])
            nearest_idx = np.argmin(dists)
            min_dist = dists[nearest_idx]
            
            THRESHOLD = 10 # pixels
            
            # RIGHT CLICK: Precise Add/Remove
            if event.button == 3: 
                if min_dist < THRESHOLD:
                    # Remove Knot (Delete) if clicked ON a knot
                    if len(self._knots_x) > 2:
                        self._knots_x = np.delete(self._knots_x, nearest_idx)
                        self._knots_y = np.delete(self._knots_y, nearest_idx)
                        self.update_spline_preview()
                        self.canvas.draw_animated() # Instant update
                else:
                    # Add Knot (Create) if clicked in EMPTY space
                    # We insert it sorted by X to maintain spline logic
                    insert_idx = np.searchsorted(self._knots_x, click_x)
                    self._knots_x = np.insert(self._knots_x, insert_idx, click_x)
                    self._knots_y = np.insert(self._knots_y, insert_idx, click_y)
                    self.update_spline_preview()
                    self.canvas.draw_animated() # Instant update
                return True 

            # LEFT CLICK: "Magnet" Grab
            if event.button == 1:
                # --- FIX 2: No Threshold ---
                # We grab 'nearest_idx' regardless of distance.
                # The entire vertical strip belongs to this node.
                self._active_knot_idx = nearest_idx
                return True

        return False
    
    def on_release(self, event):
        # 1. Handle Continuum Editor Release
        if self._edit_mode_active and self._active_knot_idx is not None:
            logging.debug(f"Released knot {self._active_knot_idx}")
            self._active_knot_idx = None  # <--- This "drops" the knot
            return  # Stop processing further (e.g. region selection)

        # 2. Handle Existing Region Selection Logic
        if self._selection_start_x is not None:
            start = self._selection_start_x
            end = event.xdata
            self.canvas.cleanup_selection()
            self._selection_start_x = None
            
            if start is None or end is None or start == end: 
                return
            
            xmin, xmax = sorted([start, end])
            if self.main_window:
                logging.info(f"Region selected: {xmin:.4f} to {xmax:.4f}")
                self.main_window.launch_split_from_region(xmin, xmax)

    def start_continuum_edit(self):
        """Initializes the continuum editor with decimated knots."""
        if not self.main_window.active_history: return
        
        # 1. PREPARE THE CANVAS (The "Clean Slate" Step)
        # We process toolbar changes and force a full draw FIRST.
        # This flushes any pending events and ensures the background is captured
        # BEFORE we introduce our fragile animated artists.
        if hasattr(self, 'toolbar') and self.toolbar:
            if self.toolbar.mode == 'pan/zoom':
                self.toolbar.pan() 
            elif self.toolbar.mode == 'zoom rect':
                self.toolbar.zoom()

        # Force the full redraw now. 
        # If this method clears 'cursor_artists', our new artists won't be there yet to get killed.
        self.canvas.draw() 

        # 2. DATA PREPARATION
        spec = self.main_window.active_history.current_state.spec
        x_full = spec.x.value
        if spec.cont is not None:
            y_source = spec.cont.value
        else:
            y_source = spec.y.value
            
        stride = 500
        indices = np.arange(0, len(x_full), stride, dtype=int)
        if indices[-1] != len(x_full) - 1:
            indices = np.append(indices, len(x_full) - 1)
            
        self._knots_x = x_full[indices]
        self._knots_y = y_source[indices]
        
        # 3. NUCLEAR CLEANUP (Just in case)
        self._remove_editor_artists()

        # 4. CREATE ARTISTS (Fresh)
        # Note: We create them AFTER the self.canvas.draw() above.
        self.knot_artist, = self.canvas.axes.plot(
            self._knots_x, self._knots_y, 
            'o', color='black', markersize=4, alpha=0.5, 
            label='Knots', zorder=10, animated=True
        )
        
        self.spline_artist, = self.canvas.axes.plot(
            [], [], 
            '-', color='black', lw=3, alpha=0.3, 
            label='Draft Cont.', zorder=9, animated=True
        )

        self._edit_mode_active = True
        self.update_spline_preview()

        # 5. REGISTRATION & MANUAL PAINT
        # We append them now. They missed the main draw() cycle (good!), 
        # so we must paint them manually.
        if self.knot_artist not in self.canvas.cursor_artists:
            self.canvas.cursor_artists.append(self.knot_artist)
        if self.spline_artist not in self.canvas.cursor_artists:
            self.canvas.cursor_artists.append(self.spline_artist)

        # Paint them on top of the already-drawn background
        self.canvas.draw_animated()

        # 1. Force Qt focus to the canvas widget
        self.canvas.setFocus()
        
        # 2. Force a Qt-level repaint (not just Matplotlib)
        # This tells the window manager "this area is dirty, please allow updates"
        self.canvas.update()

    def stop_continuum_edit(self, save=False):
        import traceback
        print("STOP called by:")
        traceback.print_stack()
        """Cleans up the editor."""
        self._edit_mode_active = False
        self._active_knot_idx = None
        
        # Capture data before destruction if saving
        final_x, final_y = None, None
        if save:
            final_x, final_y = self._knots_x, self._knots_y

        # --- FIX: FULL DESTRUCTION ---
        # Don't just hide them. Remove them. 
        # This ensures the next 'start' begins with a clean slate.
        self._remove_editor_artists()
            
        self.canvas.draw()
        
        return final_x, final_y

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

    def update_spline_preview(self):
        """Recalculates the spline and updates the artist."""
        if self._knots_x is None or len(self._knots_x) < 2: return
        
        try:
            # Sort is required for Spline logic if we added knots out of order
            sort_idx = np.argsort(self._knots_x)
            self._knots_x = self._knots_x[sort_idx]
            self._knots_y = self._knots_y[sort_idx]
            
            # Create Spline
            cs = CubicSpline(self._knots_x, self._knots_y)
            
            # Evaluate on visible range (or full range)
            # For speed, we can evaluate on the same x-grid as the plot data
            # cached in self._cached_x_plot_raw, or generate a linspace
            xlim = self.canvas.axes.get_xlim()
            x_eval = np.linspace(xlim[0], xlim[1], 1000)
            y_eval = cs(x_eval)
            
            self.spline_artist.set_data(x_eval, y_eval)
            self.knot_artist.set_data(self._knots_x, self._knots_y)
            
        except Exception as e:
            logging.error(f"Spline update failed: {e}")

    def get_knots(self):
        return self._knots_x, self._knots_y

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
                icon = qta.icon('mdi.selection-drag', color='black')
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