import logging
import numpy as np
from typing import List, Optional
from PySide6.QtCore import Qt, QAbstractTableModel, QModelIndex, Signal
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QTableView, 
    QPushButton, QLabel, QHeaderView, QAbstractItemView, QFrame,
    QScrollArea, QSizePolicy, QCheckBox
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
import matplotlib.ticker as ticker

from ..structures import ComponentDataV2
from ..atomic_data import STANDARD_MULTIPLETS, xem_d

# --- 1. Table Configuration ---
COL_MAP = {
    "ID": "id",
    "Transitions": "series", 
    "z": "z", 
    "logN": "logN", 
    "b": "b", 
    "btur": "btur", 
    "dz": "dz", 
    "dlogN": "dlogN", 
    "db": "db"
}
COLUMNS = list(COL_MAP.keys())

class SystemTableModel(QAbstractTableModel):
    """ Read-only model for ComponentDataV2. """
    def __init__(self, components: List[ComponentDataV2] = None):
        super().__init__()
        self._components = components if components else []

    def update_data(self, components: List[ComponentDataV2]):
        self.beginResetModel()
        self._components = components
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return len(self._components)

    def columnCount(self, parent=QModelIndex()):
        return len(COLUMNS)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        
        comp = self._components[index.row()]
        col_label = COLUMNS[index.column()]
        attr_name = COL_MAP[col_label] 

        if role == Qt.DisplayRole:
            val = getattr(comp, attr_name, None)
            if isinstance(val, float):
                if attr_name == 'z': return f"{val:.6f}"
                elif attr_name in ['logN', 'dlogN']: return f"{val:.3f}"
                else: return f"{val:.2f}"
            return str(val) if val is not None else ""
            
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal: return COLUMNS[section]
            else: return str(section + 1)
        return None
    
    def get_component_at(self, row: int) -> Optional[ComponentDataV2]:
        if 0 <= row < len(self._components):
            return self._components[row]
        return None

class VelocityPlotWidget(QWidget):
    """
    A unified Matplotlib widget that handles multiple stacked panels
    with shared X-axis, synchronization, and cursor tracking.
    """
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0,0,0,0)
        self.layout.setSpacing(0)

        # 1. Create Figure and Canvas
        self.fig = Figure(figsize=(5, 8), dpi=100) # Default tall size
        self.canvas = FigureCanvasQTAgg(self.fig)
        
        # 2. Add Toolbar
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.layout.addWidget(self.toolbar)
        
        # 3. Scroll Area for the Canvas (to handle many panels)
        self.scroll = QScrollArea()
        self.scroll.setWidget(self.canvas)
        self.scroll.setWidgetResizable(True)
        self.layout.addWidget(self.scroll)

        # Internal State
        self.axes = []
        self.cursor_lines = [] # Vertical lines for cursor
        self.bg_cache = None   # Background for blitting
        
        # Connect events for cursor
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('draw_event', self.on_draw)

    def plot_system(self, session, component: ComponentDataV2):
        """
        Main plotting routine. Calculates velocity, normalizes flux, 
        and sets up the subplot grid.
        """
        self.fig.clear()
        self.axes = []
        self.cursor_lines = []
        self.bg_cache = None

        if not session or not component:
            self.canvas.draw()
            return

        # 1. Determine Transitions
        series_name = component.series
        if series_name in STANDARD_MULTIPLETS:
            transitions = STANDARD_MULTIPLETS[series_name]
        else:
            # Check for single line fallback
            found = False
            for group, lines in STANDARD_MULTIPLETS.items():
                if series_name in lines:
                    transitions = lines
                    found = True
                    break
            if not found:
                transitions = [series_name]

        num_plots = len(transitions)
        if num_plots == 0: return

        # 2. Resize Canvas Height dynamically based on number of panels
        # 150px per panel minimum
        total_height = max(600, num_plots * 150) 
        self.canvas.setMinimumHeight(total_height)
        
        # 3. Create Subplots with Shared X and No Gap
        # sharex=True enables synchronized zooming!
        axs = self.fig.subplots(nrows=num_plots, ncols=1, sharex=True, sharey=True)
        if num_plots == 1: axs = [axs] # Ensure list
        
        # Remove vertical space between plots
        self.fig.subplots_adjust(hspace=0, left=0.1, right=0.95, top=0.95, bottom=0.05)

        # 4. Prepare Data
        # Check for Continuum
        has_cont = session.spec.cont is not None
        y_full = session.spec.y.value
        if has_cont:
            cont_full = session.spec.cont.value
            # Normalize (avoid div by zero)
            y_norm = np.divide(y_full, cont_full, out=np.ones_like(y_full), where=cont_full!=0)
        else:
            y_norm = y_full # Fallback to raw if no cont
            
        x_full = session.spec.x.value
        c_kms = 299792.458
        window_kms = 600 # Default view window

        # 5. Plot Loop
        for i, (ax, trans_name) in enumerate(zip(axs, transitions)):
            self.axes.append(ax)
            
            if trans_name not in xem_d:
                ax.text(0.5, 0.5, f"Unknown: {trans_name}", ha='center', transform=ax.transAxes)
                continue

            try:
                lambda_0 = xem_d[trans_name].to(session.spec.x.unit).value
            except: continue

            lambda_obs = lambda_0 * (1.0 + component.z)
            
            # Convert X to Velocity relative to this line
            # Optimization: Only process data roughly in the window to save memory/time? 
            # Actually, for zoom to work, we usually plot more. 
            # Let's plot a generous chunk (e.g. +/- 20000 km/s) to allow panning.
            
            v_full = c_kms * (x_full - lambda_obs) / lambda_obs
            
            # Limit plotting range to avoid rendering millions of points
            mask = (v_full > -5000) & (v_full < 5000)
            v_plot = v_full[mask]
            y_plot = y_norm[mask]
            
            # --- Draw Data (Step style like main plot) ---
            ax.step(v_plot, y_plot, where='mid', color='black', lw=0.8)
            
            # --- Draw Model ---
            if session.spec.model is not None:
                mod_full = session.spec.model.value
                if has_cont:
                    mod_plot = np.divide(mod_full[mask], cont_full[mask], out=np.ones_like(y_plot), where=cont_full[mask]!=0)
                else:
                    mod_plot = mod_full[mask]
                ax.plot(v_plot, mod_plot, color='red', lw=1.2, alpha=0.8)

            # --- Decorations ---
            ax.axvline(0, color='blue', ls='--', alpha=0.4) # Center
            ax.axhline(1.0, color='green', ls=':', alpha=0.5, lw=0.8) # Continuum level
            ax.axhline(0.0, color='gray', lw=0.5) # Zero level

            # Add Label inside the plot (Top Right)
            ax.text(0.98, 0.85, trans_name, transform=ax.transAxes, 
                    ha='right', fontsize=10, fontweight='bold',
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

            # Grid
            ax.grid(True, linestyle=':', alpha=0.6)
            
            # Initialize Cursor Line (invisible at start)
            line = ax.axvline(0, color='orange', ls='-', alpha=0.7, lw=1.0, visible=False)
            self.cursor_lines.append(line)

        # 6. Final Formatting
        axs[-1].set_xlabel("Velocity (km/s)")
        axs[len(axs)//2].set_ylabel("Normalized Flux") # Middle label
        
        # Set default view
        axs[0].set_xlim(-window_kms, window_kms)
        axs[0].set_ylim(-0.1, 1.2) # Standard normalized view

        self.canvas.draw()

    # --- Cursor Logic (Blitting) ---
    def on_draw(self, event):
        """Capture background for fast cursor updates."""
        if self.axes:
            self.bg_cache = self.canvas.copy_from_bbox(self.fig.bbox)
            for line in self.cursor_lines:
                line.set_visible(False)

    def on_mouse_move(self, event):
        """Update crosshair cursor across all panels."""
        if not self.axes or not event.inaxes: return
        
        # 1. Restore background
        if self.bg_cache:
            self.canvas.restore_region(self.bg_cache)
            
        # 2. Move all lines to mouse X (Velocity)
        v_mouse = event.xdata
        for line in self.cursor_lines:
            line.set_xdata([v_mouse])
            line.set_visible(True)
            line.axes.draw_artist(line)
            
        # 3. Blit
        self.canvas.blit(self.fig.bbox)


class SystemInspector(QWidget):
    """
    Floating window for inspecting systems.
    Layout: Table (Left) | Velocity Plots (Right)
    """
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.current_session = None
        
        self.setWindowTitle("System Inspector")
        self.resize(1200, 700) # Larger default size
        self.setWindowFlags(Qt.Window)
        
        self._setup_ui()
        
    def _setup_ui(self):
        main_layout = QVBoxLayout(self)
        
        # --- Top Toolbar ---
        toolbar = QHBoxLayout()
        self.btn_refit = QPushButton("Refit Selected")
        self.btn_delete = QPushButton("Delete")
        self.btn_freeze = QPushButton("Freeze") 
        
        # Placeholders
        self.btn_refit.setEnabled(False)
        self.btn_delete.setEnabled(False)
        self.btn_freeze.setEnabled(False)
        
        toolbar.addWidget(self.btn_refit)
        toolbar.addWidget(self.btn_freeze)
        toolbar.addWidget(self.btn_delete)
        toolbar.addStretch()
        main_layout.addLayout(toolbar)
        
        # --- Main Splitter (Table Left, Plot Right) ---
        splitter = QSplitter(Qt.Horizontal)
        
        # 1. Left: Table View
        self.table_view = QTableView()
        self.table_model = SystemTableModel()
        self.table_view.setModel(self.table_model)
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_view.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table_view.setAlternatingRowColors(True)
        
        # Styling columns
        header = self.table_view.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeToContents) # Compact columns
        
        # Connect selection
        self.table_view.selectionModel().selectionChanged.connect(self._on_selection_changed)
        
        splitter.addWidget(self.table_view)
        
        # 2. Right: Velocity Plot Widget
        self.vel_plot = VelocityPlotWidget()
        splitter.addWidget(self.vel_plot)
        
        # Set splitter ratio (Table 40%, Plot 60% for better visibility)
        splitter.setStretchFactor(0, 4)
        splitter.setStretchFactor(1, 6)
        
        main_layout.addWidget(splitter)

    def set_session(self, session):
        self.current_session = session
        if session and session.systs:
            self.table_model.update_data(session.systs.components)
        else:
            self.table_model.update_data([])
            
        # Clear plot if data reset
        # (Optional: Keep last plot or clear? Clearing is safer)
        if self.table_model.rowCount() == 0:
            self.vel_plot.fig.clear()
            self.vel_plot.canvas.draw()

    def _on_selection_changed(self, selected, deselected):
        indexes = self.table_view.selectionModel().selectedRows()
        if not indexes: return
            
        row = indexes[0].row()
        comp = self.table_model.get_component_at(row)
        
        if comp and self.current_session:
            self.vel_plot.plot_system(self.current_session, comp)
            
            # Enable buttons
            self.btn_refit.setEnabled(True)
            self.btn_delete.setEnabled(True)