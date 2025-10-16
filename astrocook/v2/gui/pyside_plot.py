# astrocook/v2/gui/pyside_plot.py

from PySide6.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QCheckBox, QGroupBox, QHBoxLayout
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np
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
        # plt.style.use('ggplot') 
        plt.style.use('science') 
        # plt.style.use('default') 
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
        
        from matplotlib import rcParams
        rcParams['path.simplify_threshold'] = 1.0 # Simplify plotting paths
        
        # We rely on the QtAgg backend's native speed for now.
        if parent:
            self.setParent(parent)

class SpectrumViewerPySide(QMainWindow):
    """
    Embryo for the V2 PySide spectrum viewer.
    Replaces the primary plotting logic of v1/gui_graph.py.
    """
    def __init__(self, session: 'SessionV2'):
        super().__init__()
        self.setWindowTitle("Astrocook V2 PySide Plot (Prototype)")
        self.session = session
        
        # 1. Set Initial Size (Aesthetics)
        self.resize(1200, 800) 

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        
        # Main Vertical Layout (Holds everything)
        self.main_layout = QVBoxLayout(self.central_widget)

        # 2. Initialize and Place the Matplotlib Canvas (MUST BE FIRST)
        self.canvas = MatplotlibCanvas(self.central_widget)
        self.main_layout.addWidget(self.canvas, 1)
        
        # 3. Create Horizontal Container for Toolbar and Toggles
        self.controls_container = QWidget()
        self.controls_layout = QHBoxLayout(self.controls_container)
        self.controls_layout.setContentsMargins(0, 0, 0, 0) 
        
        # 4. Initialize and Place the Toolbar (Needs self.canvas)
        # This is the line that previously crashed due to timing.
        self.toolbar = NavigationToolbar(self.canvas, self) 
        self.controls_layout.addWidget(self.toolbar)

        self.controls_layout.addStretch(1)

        # 5. Initialize and Place Toggles (Aesthetics)
        self.error_checkbox = QCheckBox("Show 1σ Error")
        self.error_checkbox.setChecked(True)
        self.error_checkbox.stateChanged.connect(self.plot_spectrum)
        
        self.continuum_checkbox = QCheckBox("Show Continuum")
        self.continuum_checkbox.setChecked(False)
        self.continuum_checkbox.stateChanged.connect(self.plot_spectrum)
        
        self.controls_layout.addWidget(self.error_checkbox)
        self.controls_layout.addWidget(self.continuum_checkbox)
        
        # 6. Add the combined Controls Container below the canvas
        self.main_layout.addWidget(self.controls_container, 0)
        
        # 7. Draw the initial spectrum
        self.plot_spectrum()

    def plot_spectrum(self):
        """
        Retrieves data from the immutable V2 Session and plots it.
        This validates the V2 data access (SessionV2 -> SpectrumV2 -> Quantity).
        """
        spec = self.session.spec
        
        if spec and len(spec.x) > 0:
            # V2 Data Access: Use the clean, immutable properties
            x_data = spec.x.value
            y_data = spec.y.value
            dy_data = spec.dy.value if spec.dy is not None else np.zeros_like(y_data)
            x_unit = str(spec.x.unit)
            y_unit = str(spec.y.unit)
            
            ax = self.canvas.axes
            ax.clear()
            
            # Use step plot for spectrum data integrity check
            ax.step(x_data, y_data, where='mid', label="Spectrum")

            # 2. Plot Error (Conditional)
            if hasattr(self, 'error_checkbox') and self.error_checkbox.isChecked():
                ax.fill_between(
                    x_data, 
                    y_data - dy_data, 
                    y_data + dy_data, 
                    step='mid', 
                    color='#aaaaaa', 
                    alpha=0.5, 
                    label='1 sigma uncertainty'
                )

            # 3. Plot Continuum (Conditional - Requires V2 data access)
            if hasattr(self, 'continuum_checkbox') and self.continuum_checkbox.isChecked():
                # Use the V2 API: Check if 'cont' column exists and retrieve it
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

            ax.legend(loc='best')
            ax.set_xlabel(f"Wavelength ({x_unit})")
            ax.set_ylabel(f"Flux ({y_unit})")
            ax.set_title(f"Spectrum: {self.session.name}")
            ax.grid(True, linestyle=':')
            self.canvas.draw()
        else:
            self.canvas.axes.set_title("No Spectrum Data Loaded")
            self.canvas.draw()

    def resizeEvent(self, event):
        """Overrides Qt resize event to redraw and apply tight layout."""
        super().resizeEvent(event)
        
        # Ensure the canvas and figure exist before attempting to draw
        if hasattr(self, 'canvas') and self.canvas.figure:
            # Tell Matplotlib to fit everything within the new bounds
            self.canvas.figure.tight_layout()
            
            # Request a redraw
            self.canvas.draw_idle()

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