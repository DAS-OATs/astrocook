# astrocook/v2/gui/pyside_plot.py

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

class MatplotlibCanvas(FigureCanvasQTAgg):
    """A custom Matplotlib canvas for embedding in a PySide widget."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):

        # --- Apply Matplotlib Style ---
        # Choose one style you like:
        #mplstyle.use('ggplot') 
        mplstyle.use('science') 
        #mplstyle.use('default') 
        mplstyle.use('fast') # Leave 'fast' to improve performance
        # ----------------------------

        # Create a Matplotlib Figure object
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        # Apply basic style improvements immediately:
        #self.axes.tick_params(direction='in', top=True, right=True)
        #self.axes.grid(True, which='major', linestyle=':', alpha=0.6)
        
        # Apply the tight layout after setting up the axes
        #fig.tight_layout()

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

class SpectrumPlotWidget(QWidget):
    """
    Refactored widget that contains the plot canvas, toolbar, and toggles.
    This replaces the content that was previously inside SpectrumViewerPySide.
    """
    def __init__(self, session_manager: 'SessionV2'):
        super().__init__()
        self.session_manager = session_manager
        
        # Main Vertical Layout for the plot area
        self.main_layout = QVBoxLayout(self)

        # 1. Initialize Canvas (Draw object is created)
        self.canvas = MatplotlibCanvas(self)
        self.main_layout.addWidget(self.canvas, 1) # Stretch factor 1

        # 2. Controls Container Setup (Stretch factor 0)
        self.controls_container = QWidget()
        self.controls_layout = QHBoxLayout(self.controls_container)
        self.controls_layout.setContentsMargins(0, 0, 0, 0)
        
        # 3. Initialize Toolbar (Needs canvas)
        # NOTE: We assume AstrocookToolbar is defined and functional.
        self.toolbar = AstrocookToolbar(self.canvas, self)
        self.controls_layout.addWidget(self.toolbar)

        # 4. Add Stretch Spacer (Pushes controls right)
        self.controls_layout.addStretch(1) 

        # 5. Initialize and Place Toggles (Fixed position)
        self.error_checkbox = QCheckBox("Show 1σ Error")
        self.error_checkbox.setChecked(True)
        self.error_checkbox.stateChanged.connect(self.plot_spectrum)
        
        self.continuum_checkbox = QCheckBox("Show Continuum")
        self.continuum_checkbox.setChecked(False)
        self.continuum_checkbox.stateChanged.connect(self.plot_spectrum)
        
        self.controls_layout.addWidget(self.error_checkbox)
        self.controls_layout.addWidget(self.continuum_checkbox)
        
        # 6. Add the combined Controls Container to the main layout
        self.main_layout.addWidget(self.controls_container, 0) # Stretch factor 0
        
        # Draw initial plot
        self.plot_spectrum(initial_draw=True)

    # --- Refactored Plotting Method (Requires an update) ---
    def plot_spectrum(self, initial_draw=False):
        """
        Retrieves data from the immutable V2 Session and plots it.
        """
        spec = self.session_manager.spec
        
        if spec and len(spec.x) > 0:
            
            # --- V2 Data Access: Retrieve data and check for safety ---
            x_data = spec.x.value
            y_data = spec.y.value
            # Safely handle missing error data by defaulting to zeros
            dy_data = spec.dy.value if spec.dy is not None else np.zeros_like(y_data)
            x_unit = str(spec.x.unit)
            y_unit = str(spec.y.unit)
            
            ax = self.canvas.axes
            ax.clear()
            
            # 1. Plot Main Flux (Uses Matplotlib style defaults for color)
            ax.step(x_data, y_data, where='mid', label="Spectrum", lw=1.0)

            # 2. Plot Error Shading (Conditional)
            if hasattr(self, 'error_checkbox') and self.error_checkbox.isChecked():
                ax.fill_between(
                    x_data, 
                    y_data - dy_data, 
                    y_data + dy_data, 
                    step='mid', 
                    color='#aaaaaa',  # Muted gray for background uncertainty
                    alpha=0.5, 
                    label='1 sigma uncertainty'
                )

            # 3. Plot Continuum (Conditional and V2 API check)
            if hasattr(self, 'continuum_checkbox') and self.continuum_checkbox.isChecked():
                cont_data_q = spec.get_column('cont')
                if cont_data_q is not None:
                    ax.plot(
                        x_data, 
                        cont_data_q.value, 
                        linestyle='--', 
                        color='red', 
                        label='Continuum'
                    )
                else:
                    logging.warning("Continuum requested but 'cont' column not found in V2 spectrum.")

            # 4. Final Touches
            ax.legend(loc='best')
            ax.set_xlabel(f"Wavelength ({x_unit})")
            ax.set_ylabel(f"Flux ({y_unit})")
            ax.set_title(f"Spectrum: {self.session_manager.name}")
            ax.grid(True, linestyle=':')
            
        else:
            self.canvas.axes.set_title("No Spectrum Data Loaded")
        
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
        super().__init__(canvas, parent, coordinates=False) 
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