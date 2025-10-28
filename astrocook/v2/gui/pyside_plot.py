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
class MatplotlibCanvas(FigureCanvasQTAgg):
    """A custom Matplotlib canvas for embedding in a PySide widget."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):

        # --- Apply Matplotlib Style ---
        try:
            plt.style.use(['science', 'nature','fast'])
            params = {'legend.fontsize': 11,
                      'axes.labelsize': 11,
                      'axes.titlesize': 12,
                      'xtick.labelsize': 10,
                      'ytick.labelsize': 10}
            plt.rcParams.update(params)
        except Exception as e:
            logging.warning(f"Could not apply 'scienceplots' style: {e}. Using default.")
            plt.style.use('fast') # Fallback to fast

        # Create a Matplotlib Figure object
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        # Apply basic style improvements immediately:
        #self.axes.tick_params(direction='in', top=True, right=True)
        #self.axes.grid(True, which='major', linestyle=':', alpha=0.6)
        
        # Apply the tight layout after setting up the axes
        fig.tight_layout()

        super(MatplotlibCanvas, self).__init__(fig)
        
        # NOTE: Blitting is usually handled automatically by NavigationToolbar2QT,
        # but performance issues can arise if the canvas draw method is overriding defaults.
        # Check that you are not suppressing the backend's default interactive rendering.
        # If the slowness persists after checking the basics, 
        # the simplest general speedup is using 'fast' style:
        
        #from matplotlib import rcParams
        #rcParams['path.simplify_threshold'] = 1.0 # Simplify plotting paths
        
        # We rely on the QtAgg backend's native speed for now.
        if parent:
            self.setParent(parent)

    def _on_draw_event(self, event):
        """Captures the canvas background after the first full draw."""
        # This hook ensures we only try to capture the background once, after axes are ready.
        if not self.blit_initialized and self.figure.get_size_inches().any():
            self.background = self.copy_from_bbox(self.figure.bbox)
            self.blit_initialized = True


def get_color_cycle(n=10, fallback_cmap='viridis'):
    """Gets the current axes.prop_cycle colors, with a fallback."""
    try:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        # Ensure we have enough colors, cycle if needed
        return [colors[i % len(colors)] for i in range(n)]
    except Exception as e:
        logging.warning(f"Could not get color cycle from rcParams: {e}. Using fallback cmap.")
        # Fallback using a colormap
        cmap = plt.get_cmap(fallback_cmap)
        return [cmap(i / n) for i in range(n)]
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
        self.canvas = MatplotlibCanvas(self)
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
        was_zoomed = False
        current_xlim = ax.get_xlim()
        current_ylim = ax.get_ylim()
        # Check if axes were not in default state (0,1) or uninitialized
        if current_xlim != (0.0, 1.0) or current_ylim != (0.0, 1.0):
             # More robust check: are limits valid numbers and different?
             if all(isinstance(v, (int, float)) for v in current_xlim + current_ylim) and \
                current_xlim[0] != current_xlim[1] and current_ylim[0] != current_ylim[1]:
                 was_zoomed = True
                 logging.debug(f"Plot was zoomed. Xlim={current_xlim}, Ylim={current_ylim}")

        ax.clear() # Clear after getting limits

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
                    series_str = self.main_window.cursor_series_input.text()
                    z_cursor = float(self.main_window.cursor_z_input.text())
                    transitions = trans_parse(series_str)
                    xlim = ax.get_xlim() # Get current limits
                    added_cursor_label = False
                    for t in transitions:
                        if t in xem_d:
                            xem_nm = xem_d[t].to_value(au.nm)
                            x_obs = (1 + z_cursor) * xem_nm * au.nm
                            # Convert to plot's x-unit
                            try:
                                zem = getattr(spec, '_zem', 0.0)
                                x_plot = x_convert(x_obs, zem=zem, xunit=spec.x.unit).value
                            except Exception: continue # Skip if conversion fails

                            if xlim[0] <= x_plot <= xlim[1]:
                                label = f"Cursor (z={z_cursor:.4f})" if not added_cursor_label else None
                                ax.axvline(x_plot, ls='--', color='blue', alpha=0.8, lw=1.0, label=label)
                                added_cursor_label = True
                        else:
                             logging.warning(f"Cursor transition '{t}' for series '{series_str}' not found in xem_d.")

                except ValueError:
                    logging.warning("Invalid redshift value entered for cursor.")
                except Exception as e:
                    logging.error(f"Error drawing redshift cursor: {e}")
            # -----------------------------------------------

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

        # Redraw the canvas
        if initial_draw:
            self.canvas.draw()
        else:
            self.canvas.draw_idle()

        # NOTE: You must update the self.session reference here if the session changes.
        # This function will be called by update_plot().

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