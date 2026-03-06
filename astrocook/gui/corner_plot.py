from PySide6.QtWidgets import QMainWindow, QVBoxLayout, QHBoxLayout, QWidget, QScrollBar, QLabel, QSizePolicy
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from .pyside_plot import AstrocookToolbar 
from matplotlib.figure import Figure
import corner
import numpy as np

class CornerPlotWindow(QMainWindow):
    """
    A window for displaying corner plots of Bayesian sampling results.
    Features paged scrolling (horizontal and vertical) for large parameter sets.
    Matches layout and styling of System Inspector.
    """
    def __init__(self, samples: np.ndarray, labels: list, title: str = "Corner Plot", parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(1000, 1000)
        
        # 1. Parameter Re-ordering (z, then logN, then b, then others)
        param_order = []
        for p_type in ['z', 'logN', 'b']:
            for i, label in enumerate(labels):
                if label.startswith(f"{p_type}_") or label == p_type:
                    param_order.append(i)
        
        remaining = set(range(len(labels))) - set(param_order)
        param_order.extend(sorted(list(remaining)))
        
        self.samples = samples[:, param_order]
        self.labels = [labels[i] for i in param_order]
        
        self.PAGE_SIZE = 6 
        
        # --- Layout Matching System Inspector ---
        container = QWidget()
        self.setCentralWidget(container)
        
        # Main layout of the window (Margins 12, 12, 12, 12 as per user request/standard)
        self.main_layout = QVBoxLayout(container)
        self.main_layout.setContentsMargins(12, 12, 12, 12)
        self.main_layout.setSpacing(0)
        
        # Center area (HBox) containing Plot Stack and Vertical Scrollbar
        self.center_hbox = QHBoxLayout()
        self.center_hbox.setContentsMargins(0, 0, 0, 0)
        self.center_hbox.setSpacing(2) # Spacing between plot and scrollbar
        
        # Left Stack (VBox) containing Canvas, Horizontal Scrollbar, and Toolbar
        self.left_stack = QVBoxLayout()
        self.left_stack.setContentsMargins(0, 0, 0, 0)
        self.left_stack.setSpacing(5) # Spacing between canvas/scrollbar/toolbar
        
        self.figure = Figure(figsize=(8, 8), dpi=100)
        self.figure.patch.set_facecolor('#ffffff') 
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.left_stack.addWidget(self.canvas, stretch=1)
        
        self.scroll_h = QScrollBar(Qt.Horizontal)
        self.scroll_h.valueChanged.connect(self.update_plot)
        self.left_stack.addWidget(self.scroll_h)
        
        self.toolbar = AstrocookToolbar(self.canvas, self)
        # Match System Inspector: no special container, just added to the layout
        # We ensure it's transparent and flat
        self.toolbar.setStyleSheet("border: none; background: transparent;")
        self.left_stack.addWidget(self.toolbar)
        
        self.center_hbox.addLayout(self.left_stack, stretch=1)
        
        # Vertical Scrollbar on the right, spanning the whole left stack
        self.scroll_v = QScrollBar(Qt.Vertical)
        self.scroll_v.valueChanged.connect(self.update_plot)
        self.center_hbox.addWidget(self.scroll_v)
        
        self.main_layout.addLayout(self.center_hbox, stretch=1)
        
        num_params = len(self.labels)
        max_idx = max(0, num_params - self.PAGE_SIZE)
        self.scroll_h.setRange(0, max_idx)
        self.scroll_v.setRange(0, max_idx)
        
        self.update_plot()
        
    def update_plot(self):
        self.figure.set_layout_engine(None)
        self.figure.clear()
        
        x_start = self.scroll_h.value()
        y_start = self.scroll_v.value()
        
        x_end = min(x_start + self.PAGE_SIZE, len(self.labels))
        y_end = min(y_start + self.PAGE_SIZE, len(self.labels))
        
        x_labels = self.labels[x_start:x_end]
        y_labels = self.labels[y_start:y_end]
        
        color_blue = '#1f77b4'   
        color_orange = '#ff7f0e' 
        
        num_x = len(x_labels)
        num_y = len(y_labels)
        
        axs = self.figure.subplots(num_y, num_x, sharex='col', sharey=False, squeeze=False)
        self.figure.subplots_adjust(hspace=0, wspace=0)
        
        for i in range(num_y): 
            idx_y = y_start + i
            master_y_ax = None
            for j in range(num_x):
                idx_x = x_start + j
                if idx_y > idx_x:
                    master_y_ax = axs[i, j]
                    break

            for j in range(num_x): 
                idx_x = x_start + j
                ax = axs[i, j]
                ax.set_facecolor('#ffffff') 
                # [NEW] Enforce square aspect ratio for boxes
                ax.set_box_aspect(1)
                
                if idx_y < idx_x:
                    ax.set_visible(False)
                    continue

                if idx_x == idx_y:
                    # Diagonal: Histogram
                    ax.get_yaxis().set_visible(False)
                    data = self.samples[:, idx_x]
                    data = data[np.isfinite(data)]
                    
                    ax.hist(data, bins=30, color=color_blue, 
                            histtype='step', density=True, lw=1.5)
                    
                    qs = np.percentile(data, [16, 50, 84])
                    for q in qs:
                        ax.axvline(q, color=color_orange, ls='--', lw=1, alpha=0.5)
                    
                    ax.set_title(self.labels[idx_x], fontsize=10, fontweight='bold', pad=10, color='#333333')
                    ax.text(0.5, 0.85, f"{med:.3f} +{std_p:.3f} -{std_m:.3f}", 
                            transform=ax.transAxes, ha='center', va='top', fontsize=8, color=color_orange) if 'med' in locals() else None
                    # Recalculate med for title text if needed
                    med = qs[1]
                    std_p = qs[2] - qs[1]
                    std_m = qs[1] - qs[0]
                    ax.text(0.5, 0.85, f"{med:.3f}\n+{std_p:.3f} -{std_m:.3f}", 
                            transform=ax.transAxes, ha='center', va='top', fontsize=7, color=color_orange)

                else: 
                    if master_y_ax and ax != master_y_ax:
                        ax.sharey(master_y_ax)

                    corner.hist2d(
                        self.samples[:, idx_x], 
                        self.samples[:, idx_y], 
                        ax=ax,
                        plot_datapoints=False,
                        plot_density=True,
                        no_fill_contours=False,
                        color=color_blue
                    )
                    
                    med_x = np.median(self.samples[:, idx_x])
                    med_y = np.median(self.samples[:, idx_y])
                    ax.axvline(med_x, color=color_orange, linestyle='--', alpha=0.6, lw=1)
                    ax.axhline(med_y, color=color_orange, linestyle='--', alpha=0.6, lw=1)
                
                ax.tick_params(axis='both', which='major', labelsize=7, 
                              labelleft=(j==0 and idx_y > idx_x), 
                              labelbottom=(i==num_y-1))
                
                if idx_y != idx_x:
                     ax.text(0.95, 0.05, f"{x_labels[j]} vs {y_labels[i]}", 
                            transform=ax.transAxes, ha='right', va='bottom', fontsize=5, 
                            alpha=0.4, color='#666666')

        self.canvas.draw()
