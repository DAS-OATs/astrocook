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
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QCheckBox, QHBoxLayout, QStyle
)
from PySide6.QtGui import QAction, QIcon
import scienceplots
from typing import TYPE_CHECKING

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
    
# --- Helper function for inverse Z calculation ---
def z_convert_inverse(x_plot, xem_nm, zem_spec, x_unit):
    """Converts plot X coordinate back to redshift for a given line."""
    if not V1_FUNCTIONS_AVAILABLE or not isinstance(x_unit, au.UnitBase):
        logging.warning("Cannot perform inverse Z conversion: V1 functions/units unavailable.")
        return None
    try:
        # Convert plot x back to observed wavelength in nm
        x_obs_nm = x_convert(x_plot * x_unit, zem=zem_spec, xunit=au.nm).to_value(au.nm)
        # Calculate redshift
        z = (x_obs_nm / xem_nm) - 1.0
        return z
    except Exception as e:
        logging.debug(f"Inverse z conversion failed for x={x_plot} ({x_unit}): {e}")
        return None
# -----------------------------------------------

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
        except Exception as e:
            logging.warning(f"Could not apply 'scienceplots': {e}. Using default.")
            plt.style.use('fast')

        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        # Apply basic style improvements immediately:
        #self.axes.tick_params(direction='in', top=True, right=True)
        #self.axes.grid(True, which='major', linestyle=':', alpha=0.6)
        
        # Apply the tight layout after setting up the axes
        fig.tight_layout()

        super(MatplotlibCanvas, self).__init__(fig)

        # We rely on the QtAgg backend's native speed for now.
        if parent:
            self.setParent(parent)

        self.plot_widget = plot_widget # Store reference

        # NOTE: Blitting is usually handled automatically by NavigationToolbar2QT,
        # but performance issues can arise if the canvas draw method is overriding defaults.
        # Check that you are not suppressing the backend's default interactive rendering.
        # If the slowness persists after checking the basics, 
        # the simplest general speedup is using 'fast' style:
        
        #from matplotlib import rcParams
        #rcParams['path.simplify_threshold'] = 1.0 # Simplify plotting paths
        
        # --- Blitting Attributes ---
        self.background = None # Store the clean background pixels
        self.cursor_artists = [] # Store references to the cursor axvline objects
        self.draw_event_cid = None # Store connection ID for draw_event
        self._needs_full_redraw = False
        # ---------------------------

        # --- Dragging Attributes ---
        # self.dragging_cursor = False
        # self.active_cursor_line = None # Optional: Store the specific line being dragged
        # ---------------------------

        # Connect Matplotlib events
        # self.mpl_connect('button_press_event', self.on_press)
        # self.mpl_connect('button_release_event', self.on_release)
        # self.mpl_connect('motion_notify_event', self.on_motion)
        # ... (Connect motion, limit events) ...
        self.mpl_connect('motion_notify_event', self.on_motion)
        self.xlim_cid = self.axes.callbacks.connect('xlim_changed', self.on_lim_changed)
        self.ylim_cid = self.axes.callbacks.connect('ylim_changed', self.on_lim_changed)

        #self._lim_changed = False

    def on_lim_changed(self, axes):
        """Callback: Invalidate background and flag for full redraw,
           ONLY if a stable background already exists."""

        # ** CRITICAL CHECK: **
        # Only invalidate if a background has already been captured.
        # This prevents on_lim_changed from firing during the
        # *initial* plot_spectrum call (when background is None).
        if self.background is not None:
            logging.debug(f"Limits changed ({axes.get_label()}) *after* init: Invalidating background & flagging redraw.")
            self.background = None # Invalidate background
            self._needs_full_redraw = True # Flag that the *next* motion event needs to trigger a full plot
        else:
            logging.debug(f"Limits changed ({axes.get_label()}) *during* init/redraw. Ignoring.")
            # Do nothing - a redraw is already in progress which will capture the background.

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


        # ** 1. Check if a full redraw is flagged **
        if self._needs_full_redraw:
            self._needs_full_redraw = False # Reset flag
            self.plot_widget.plot_spectrum() # Trigger full redraw
            return # Stop here, don't try to blit this frame
        
        # Check if inside axes and cursor checkbox is checked
        if (event.inaxes == self.axes and event.xdata is not None
            and main_window.cursor_show_checkbox.isChecked()):

            # --- Always update on motion if cursor is active ---
            mouse_x = event.xdata
            spec = self.plot_widget.session_manager.spec
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
                new_z = z_convert_inverse(mouse_x, ref_xem_nm, zem_spec, x_unit)

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
    def __init__(self, session_manager: 'SessionV2', main_window_ref):
        super().__init__()
        self.session_manager = session_manager
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
        self.plot_spectrum(initial_draw=True)

    # --- Refactored Plotting Method (Requires an update) ---
    def plot_spectrum(self, initial_draw=False):
        """
        Retrieves data from the immutable V2 Session and plots it.
        """
        spec = self.session_manager.spec
        systs = self.session_manager.systs # Get systs object

        ax = self.canvas.axes

        # ** Store current limits BEFORE clearing **
        # ** Make sure to store current xlim/ylim BEFORE ax.clear() **
        was_zoomed = False; current_xlim = (0.0, 1.0); current_ylim = (0.0, 1.0) # Initialize
        try: # Get current limits robustly
             current_xlim = ax.get_xlim()
             current_ylim = ax.get_ylim()
             if current_xlim != (0.0, 1.0) or current_ylim != (0.0, 1.0):
                 if all(isinstance(v,(int,float)) for v in current_xlim+current_ylim) and \
                    current_xlim[0]!=current_xlim[1] and current_ylim[0]!=current_ylim[1]:
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
        # -----------------------------

        # ** Reset background and clear artists for blitting **
        self.canvas.background = None
        self.canvas._needs_full_redraw = False # Reset flag on *every* full draw
        self.canvas.cursor_artists = []
        # ----------------------------------------------------
        
        if spec and len(spec.x) > 0:
            
            # --- V2 Data Access: Retrieve data and check for safety ---
            x_data = spec.x.value
            y_data = spec.y.value
            # Safely handle missing error data by defaulting to zeros
            dy_data = spec.dy.value if spec.dy is not None else np.zeros_like(y_data)
            x_unit = str(spec.x.unit)
            y_unit = str(spec.y.unit)
                    
            colors = get_color_cycle(5)

            # 1. Plot Main Flux (Uses Matplotlib style defaults for color)
            ax.step(x_data, y_data, where='mid', label="Spectrum", lw=0.5, color=colors[0])

            # --- Check state from main_window reference ---
            # 2. Plot Error Shading (Conditional)
            if self.main_window.error_checkbox.isChecked(): # <<< Check main window's checkbox
                ax.fill_between(
                    x_data, y_data - dy_data, y_data + dy_data,
                    step='mid', color='#aaaaaa', alpha=0.5,
                    label='1-sigma error'
                )

            # 3. Plot Continuum (Conditional)
            if self.main_window.continuum_checkbox.isChecked(): # <<< Check main window's checkbox
                cont_data_q = spec.get_column('cont')
                if cont_data_q is not None:
                    ax.plot(
                        x_data, cont_data_q.value, linestyle='--',
                        color='black', lw=0.8, label='Continuum'
                    )
                else:
                    logging.warning("Continuum requested but 'cont' not found.")

            # 4. Plot Model
            if self.main_window.model_checkbox.isChecked():
                model_data_q = spec.get_column('model')
                if model_data_q is not None:
                    ax.plot(x_data, model_data_q.value, ls='-', color=colors[1], lw=0.8, label='Model')
                # No warning needed

            # 5. Plot Systems
            if self.main_window.systems_checkbox.isChecked() and V1_FUNCTIONS_AVAILABLE and systs and systs.components:
                # Get current plot limits to only draw visible lines
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                added_hi_label = False
                added_metal_label = False

                # Define Y position for markers (e.g., slightly below top, in axis coords)
                marker_y_axis_coord = 0.05 # 5% from the bottom

                for comp in systs.components:
                    z = comp.z
                    # Use V1 function to get transitions for the series
                    transitions = trans_parse(comp.series)
                    for t in transitions:
                        if t in xem_d:
                            # Calculate observed wavelength in current spectrum units
                            xem_nm = xem_d[t].to_value(au.nm)
                            x_obs = (1 + z) * xem_nm * au.nm
                            # Convert to plot's x-unit (might be velocity)
                            try:
                                # Need zem for velocity conversion if spec._xunit is velocity
                                zem = getattr(spec, '_zem', 0.0) # Get emission redshift if available
                                x_plot = x_convert(x_obs, zem=zem, xunit=spec.x.unit).value
                            except Exception as e:
                                # Fallback if conversion fails (e.g., incompatible units)
                                logging.debug(f"Could not convert syst wavelength {x_obs} to plot unit {spec.x.unit}: {e}")
                                continue

                            if xlim[0] <= x_plot <= xlim[1]:
                                # --- Color and Label Logic ---
                                is_hi = t.startswith('Ly_') # Simple check for Lyman series
                                color = colors[2] if is_hi else colors[3] # HI = red, Metals = blue
                                label = None
                                if is_hi and not added_hi_label:
                                    label = "HI Systems"
                                    added_hi_label = True
                                elif not is_hi and not added_metal_label:
                                    label = "Metal Systems"
                                    added_metal_label = True
                                # ---------------------------

                                # Convert marker Y position from axis to data coordinates
                                # Use blended transform: data for X, axes for Y
                                trans = ax.get_xaxis_transform() # Shortcut for blended transform
                                marker_y_data_coord = marker_y_axis_coord # For plot in axis coords

                                # Use ax.plot with marker style
                                ax.plot(x_plot, marker_y_data_coord,
                                        marker='|', markersize=24, linestyle='None', # Use '|' marker
                                        color=color, alpha=0.5, label=label,
                                        transform=trans) # Apply transform for Y axis coord
                        else:
                            logging.warning(f"Transition '{t}' for series '{comp.series}' not found in xem_d.")


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
            ax.set_title(f"Spectrum: {self.session_manager.name}")
            ax.grid(True, linestyle=':')
            
            # ** Restore Zoom/Pan if applicable **
            if was_zoomed and not initial_draw:
                logging.debug(f"Restoring zoom: Xlim={current_xlim}, Ylim={current_ylim}")
                # Check if limits are still somewhat valid (optional, prevents extreme zooms on different data)
                # For now, just restore them
                ax.set_xlim(current_xlim)
                ax.set_ylim(current_ylim)
            else:
                logging.debug("Not restoring zoom (initial draw or wasn't zoomed).")

        else:
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
        spec = self.session_manager.spec
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
                    x_obs = (1 + z_cursor) * xem_nm * au.nm
                    try:
                        # Convert observed nm to current plot unit
                        x_plot = x_convert(x_obs, zem=zem, xunit=x_unit).value
                        positions.append(x_plot)
                    except Exception as e:
                        logging.debug(f"Cursor pos calculation failed for {t}: {e}")
                        continue # Skip this transition if conversion fails
        except (ValueError, AttributeError, ImportError) as e:
            logging.warning(f"Could not calculate cursor positions: {e}")
            return []
        return positions
    
    def update_plot(self, new_session):
        """Called by MainWindowV2 to swap the immutable session object and redraw."""
        self.session_manager = new_session
        self.plot_spectrum(initial_draw=False)

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