from PySide6.QtWidgets import QMainWindow, QVBoxLayout, QWidget, QScrollArea
from PySide6.QtCore import Qt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import corner
import numpy as np

class CornerPlotWindow(QMainWindow):
    """
    A window for displaying corner plots of Bayesian sampling results.
    Features a scrollable area for large parameter sets and a 10pt border.
    """
    def __init__(self, samples: np.ndarray, labels: list, title: str = "Corner Plot", parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(900, 900)
        
        # Central widget container
        container = QWidget()
        self.setCentralWidget(container)
        main_layout = QVBoxLayout(container)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        
        # Scroll Area for the plot
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setStyleSheet("QScrollArea { border: none; }")
        
        # Content widget with 10pt margins
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)
        content_layout.setContentsMargins(10, 10, 10, 10)
        
        # Matplotlib Figure
        # We estimate size based on number of parameters
        num_params = len(labels)
        fig_size = max(8, num_params * 1.5)
        self.figure = Figure(figsize=(fig_size, fig_size), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        content_layout.addWidget(self.canvas)
        
        self.scroll_area.setWidget(content_widget)
        main_layout.addWidget(self.scroll_area)
        
        # Navigation Toolbar at the bottom
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setStyleSheet("background-color: white; border-top: 1px solid #ccc;")
        main_layout.addWidget(self.toolbar)
        
        # Generate Corner Plot
        self.plot_corner(samples, labels)
        
    def plot_corner(self, samples, labels):
        self.figure.clear()
        
        # Determine label and title font sizes based on grid density
        num_params = len(labels)
        fontsize = max(6, 12 - num_params // 2)
        
        corner.corner(
            samples, 
            labels=labels, 
            fig=self.figure, 
            show_titles=True, 
            quantiles=[0.16, 0.5, 0.84],
            title_fmt='.3f',
            label_kwargs={'fontsize': fontsize},
            title_kwargs={'fontsize': fontsize}
        )
        
        self.canvas.draw()
