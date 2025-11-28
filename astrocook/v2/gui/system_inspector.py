import logging
import numpy as np
import astropy.units as au
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.ticker as ticker
from PySide6.QtCore import Qt, QAbstractTableModel, QModelIndex, Signal, QItemSelectionModel
from PySide6.QtGui import QAction, QCursor
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QTableView, 
    QPushButton, QHeaderView, QAbstractItemView, QMenu,
    QScrollBar, QSizePolicy, QMessageBox, QLabel, QLineEdit
)
from typing import List, Optional

from ..structures import ComponentDataV2
from ..atomic_data import STANDARD_MULTIPLETS, xem_d
from .pyside_plot import AstrocookToolbar, get_color_cycle, PLOT_STYLE, get_style_color

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
EDITABLE_COLS = {"Transitions", "z", "logN", "b", "btur"}

class SystemTableModel(QAbstractTableModel):
    data_changed_request = Signal(str, dict)

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

    def flags(self, index):
        base_flags = super().flags(index)
        col_label = COLUMNS[index.column()]
        if col_label in EDITABLE_COLS:
            return base_flags | Qt.ItemIsEditable
        return base_flags

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid(): return None
        comp = self._components[index.row()]
        col_label = COLUMNS[index.column()]
        attr_name = COL_MAP[col_label] 

        if role == Qt.DisplayRole or role == Qt.EditRole:
            val = getattr(comp, attr_name, None)
            if role == Qt.DisplayRole and isinstance(val, float):
                if attr_name == 'z': return f"{val:.6f}"
                elif attr_name in ['logN', 'dlogN']: return f"{val:.3f}"
                else: return f"{val:.2f}"
            return str(val) if val is not None else ""
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter
        return None

    def setData(self, index, value, role=Qt.EditRole):
        if not index.isValid() or role != Qt.EditRole: return False
        col_label = COLUMNS[index.column()]
        attr_name = COL_MAP[col_label]
        comp = self._components[index.row()]
        try:
            if attr_name in ['z', 'logN', 'b', 'btur']:
                float_val = float(value)
                params = {attr_name: float_val}
            else:
                params = {attr_name: str(value)}
            params['uuid'] = comp.uuid
            self.data_changed_request.emit("update_component", params)
            return True 
        except ValueError:
            return False

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal: return COLUMNS[section]
            else: return str(section + 1)
        return None
    
    def get_component_at(self, row: int) -> Optional[ComponentDataV2]:
        if 0 <= row < len(self._components):
            return self._components[row]
        return None

# --- 2. The Paged Plot Widget ---

class VelocityPlotWidget(QWidget):
    """
    A widget that displays a fixed number of panels (PAGE_SIZE).
    A QScrollBar on the right allows the user to 'page' through the list of transitions.
    This avoids all issues with QScrollArea and Matplotlib resizing.
    """
    PAGE_SIZE = 3 # Number of panels visible at once

    def __init__(self, inspector_parent): 
        super().__init__()
        self.inspector = inspector_parent
        
        # Main Horizontal Layout (Canvas | ScrollBar)
        self.main_layout = QHBoxLayout(self)
        self.main_layout.setContentsMargins(0, 0, 0, 0)
        self.main_layout.setSpacing(2)

        # 1. Left Side: Canvas + Toolbar (Vertical)
        left_container = QWidget()
        self.left_layout = QVBoxLayout(left_container)
        self.left_layout.setContentsMargins(0, 0, 0, 0)
        self.left_layout.setSpacing(0)

        # Figure
        self.fig = Figure(figsize=(5, 6), dpi=100) # Fixed reasonable size
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # Toolbar
        self.toolbar = AstrocookToolbar(self.canvas, self)
        for action in self.toolbar.actions():
            if action.text() == 'Select': self.toolbar.removeAction(action); break
        
        self.left_layout.addWidget(self.canvas)
        self.left_layout.addWidget(self.toolbar)

        # 2. Right Side: ScrollBar
        self.scrollbar = QScrollBar(Qt.Vertical)
        self.scrollbar.valueChanged.connect(self._on_scroll)
        self.scrollbar.setRange(0, 0)
        self.scrollbar.setPageStep(1)

        # Assemble
        self.main_layout.addWidget(left_container)
        self.main_layout.addWidget(self.scrollbar)

        # Data State
        self._all_transitions = [] # Full list of transitions to show
        self._panel_map = {} # Maps Axes object -> Transition Name
        self.cursor_lines = []
        
        self._current_session = None
        self._current_component = None
        
        # Connect Events
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('button_press_event', self.on_press)

    @property
    def main_window(self):
        return self.inspector.main_window

    def toggle_region_selector(self): pass # No-op for compatibility

    def plot_spectrum(self, session_state=None, force_autoscale=False):
        """Standard entry point for plot updates."""
        # Use cached state if args are missing
        session = session_state if session_state else self._current_session
        comp = self._current_component
        
        if session and comp:
            self.plot_system(session, comp)
        else:
            self.fig.clear()
            self.canvas.draw()

    def plot_system(self, session, component: ComponentDataV2):
        """Prepare the list of transitions and reset the view."""
        self._current_session = session
        self._current_component = component
        
        if not session or not component:
            self._all_transitions = []
            self._update_plot()
            return

        # 1. Build Transition List
        series_name = component.series
        base_transitions = STANDARD_MULTIPLETS.get(series_name, [series_name]) if series_name else []
        
        extra_transitions = []
        for item in self.inspector.active_transitions:
            if item in STANDARD_MULTIPLETS: extra_transitions.extend(STANDARD_MULTIPLETS[item])
            elif item in xem_d: extra_transitions.append(item)

        # Combine and Deduplicate
        self._all_transitions = list(base_transitions)
        for t in extra_transitions:
            if t not in self._all_transitions: self._all_transitions.append(t)

        # 2. Configure Scrollbar
        total_items = len(self._all_transitions)
        if total_items <= self.PAGE_SIZE:
            self.scrollbar.setRange(0, 0) # Disable scrolling
            self.scrollbar.setEnabled(False)
        else:
            self.scrollbar.setRange(0, total_items - self.PAGE_SIZE)
            self.scrollbar.setEnabled(True)
            self.scrollbar.setValue(0) # Reset to top

        # 3. Draw
        self._update_plot()

    def _on_scroll(self, value):
        """Slot called when scrollbar moves. Simply redraws the page."""
        self._update_plot()

    def _update_plot(self):
        """The core rendering logic. Draws a 'page' of N panels."""
        self.fig.clear()
        self.axes = []
        self.cursor_lines = []
        self._panel_map = {}

        if not self._all_transitions or not self._current_session:
            self.canvas.draw()
            return

        # 1. Determine Window
        start_idx = self.scrollbar.value()
        # Ensure we don't go out of bounds
        end_idx = min(start_idx + self.PAGE_SIZE, len(self._all_transitions))
        visible_transitions = self._all_transitions[start_idx:end_idx]
        
        num_plots = len(visible_transitions)
        if num_plots == 0: return

        # 2. Create Subplots
        # constrained_layout=True works perfectly here because the Figure size is fixed to the widget
        self.fig.set_layout_engine(layout='constrained')
        axs = self.fig.subplots(nrows=num_plots, ncols=1, sharex=True, sharey=True)
        if num_plots == 1: axs = [axs]
        else: axs = axs.flatten()

        # 3. Plot Data
        session = self._current_session
        comp = self._current_component
        
        # Pre-fetch data once
        has_cont = session.spec.cont is not None
        y_full = session.spec.y.value
        dy_full = session.spec.dy.value if session.spec.dy is not None else np.zeros_like(y_full)
        if has_cont:
            cont_full = session.spec.cont.value
            y_norm = np.divide(y_full, cont_full, out=np.ones_like(y_full), where=cont_full!=0)
            dy_norm = np.divide(dy_full, cont_full, out=np.zeros_like(dy_full), where=cont_full!=0)
        else:
            y_norm = y_full; dy_norm = dy_full
        
        x_full = session.spec.x.value
        c_kms = 299792.458
        window_kms = 600
        z_sys = comp.z

        colors = get_color_cycle(5, cmap='tab20')

        for i, trans_name in enumerate(visible_transitions):
            ax = axs[i]
            self.axes.append(ax)
            self._panel_map[ax] = trans_name # For mouse click logic

            # Formatter
            def make_fmt(z):
                def fmt(x, y):
                    val_z = (1 + z) * (1 + x / c_kms) - 1
                    return f"v={x:.1f}, y={y:.2f}, z={val_z:.5f}"
                return fmt
            ax.format_coord = make_fmt(z_sys)

            if trans_name not in xem_d:
                ax.text(0.5, 0.5, f"Unknown: {trans_name}", ha='center'); continue

            try: lam_0 = xem_d[trans_name].to(session.spec.x.unit).value
            except: ax.text(0.5, 0.5, f"Unit Error", ha='center'); continue

            lam_obs = lam_0 * (1.0 + z_sys)
            v = c_kms * (x_full - lam_obs) / lam_obs
            
            # Mask data to window
            mask = (v > -window_kms - 200) & (v < window_kms + 200)
            v_p = v[mask]; y_p = y_norm[mask]; dy_p = dy_norm[mask]

            ax.step(v_p, y_p-dy_p, where='mid', color=PLOT_STYLE['error']['color'], lw=0.5, alpha=0.5)
            ax.step(v_p, y_p+dy_p, where='mid', color=PLOT_STYLE['error']['color'], lw=0.5, alpha=0.5)
            ax.step(v_p, y_p, where='mid', color=get_style_color('flux', colors), lw=0.8)

            if session.spec.model is not None:
                mod = session.spec.model.value
                mod_p = np.divide(mod[mask], cont_full[mask], out=np.ones_like(y_p), where=cont_full[mask]!=0) if has_cont else mod[mask]
                ax.plot(v_p, mod_p, color=get_style_color('model', colors), lw=1.0)

            ax.axvline(0, color='gray', ls='--', lw=0.8)
            ax.axhline(1.0, color='green', ls=':', alpha=0.5)
            ax.axhline(0.0, color='gray', lw=0.5)
            ax.text(0.98, 0.85, trans_name, transform=ax.transAxes, ha='right', fontweight='bold', 
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            
            ax.set_xlim(-window_kms, window_kms)
            ax.set_ylim(-0.2, 1.4)
            
            # Cursor
            l = ax.axvline(0, color=get_style_color('cursor', colors), ls='--', lw=1, visible=False)
            self.cursor_lines.append(l)

        self.fig.supxlabel("Velocity (km/s)")
        self.fig.supylabel("Normalized Flux")
        
        self.canvas.draw()

    def on_mouse_move(self, event):
        if not self.axes or not event.inaxes: return
        for l in self.cursor_lines:
            l.set_xdata([event.xdata]); l.set_visible(True)
        self.canvas.draw_idle()

    def on_press(self, event):
        if event.button == 3 and event.inaxes in self.axes:
            ax = event.inaxes
            trans = self._panel_map.get(ax)
            if trans and self.inspector.main_window:
                v = event.xdata
                z_new = (1 + self._current_component.z) * (1 + v / 299792.458) - 1
                
                # Identify series
                series = trans
                for g, l in STANDARD_MULTIPLETS.items():
                    if trans in l: series = g; break
                
                menu = QMenu(self)
                act = QAction(f"Add {series} at z={z_new:.5f}", menu)
                act.triggered.connect(lambda: self.inspector.main_window._on_recipe_requested(
                    "absorbers", "add_component", {'series': series, 'z': z_new}, {}))
                menu.addAction(act)
                menu.exec(QCursor.pos())

# --- 3. Inspector Window (Unchanged logic) ---

class SystemInspector(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.current_session = None
        self.active_transitions: List[str] = [] 
        self.setWindowTitle("System Inspector")
        self.resize(1200, 800) 
        self._setup_ui()
        
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        
        # Left
        self.table_view = QTableView()
        self.table_model = SystemTableModel()
        self.table_view.setModel(self.table_model)
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_view.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table_view.setAlternatingRowColors(True)
        self.table_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table_view.selectionModel().selectionChanged.connect(self._on_selection_changed)
        self.table_view.customContextMenuRequested.connect(self._on_context_menu)
        self.table_model.data_changed_request.connect(self._forward)
        splitter.addWidget(self.table_view)
        
        # Right
        right = QWidget()
        r_layout = QVBoxLayout(right)
        r_layout.setContentsMargins(0,0,0,0)
        
        c_layout = QHBoxLayout()
        c_layout.setContentsMargins(5,5,5,0)
        c_layout.addWidget(QLabel("Shown transitions:"))
        self.trans_in = QLineEdit()
        self.trans_in.setPlaceholderText("e.g. CIV, SiIV")
        self.trans_in.returnPressed.connect(self._apply)
        
        pal = QApplication.palette()
        self.trans_in.setStyleSheet(f"QLineEdit {{ padding: 4px; border-radius: 5px; background: {pal.color(pal.ColorRole.Base).name()}; color: {pal.color(pal.ColorRole.Text).name()}; }}")
        
        c_layout.addWidget(self.trans_in)
        r_layout.addLayout(c_layout)
        
        self.vel_plot = VelocityPlotWidget(self)
        r_layout.addWidget(self.vel_plot, 1) 
        
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 4); splitter.setStretchFactor(1, 6)
        layout.addWidget(splitter)

    def set_session(self, session):
        sel = self.table_view.selectionModel().selectedRows()
        prev_uuid = self.table_model.get_component_at(sel[0].row()).uuid if sel else None

        self.current_session = session
        self.table_model.update_data(session.systs.components if session and session.systs else [])
        
        if prev_uuid:
            for r in range(self.table_model.rowCount()):
                if self.table_model.get_component_at(r).uuid == prev_uuid:
                    self.table_view.selectRow(r)
                    # Note: plot_system called via selection signal
                    return
        
        if self.table_model.rowCount() > 0:
            self.table_view.selectRow(0)
        else:
            self.vel_plot.plot_spectrum(None)

    def _on_selection_changed(self, sel, desel):
        if not sel.indexes(): return
        row = sel.indexes()[0].row()
        comp = self.table_model.get_component_at(row)
        if comp and self.current_session:
            cur = self.trans_in.text()
            toks = [t.strip() for t in cur.split(',') if t.strip()]
            if comp.series and comp.series not in toks: toks.append(comp.series)
            new_txt = ", ".join(toks)
            self.trans_in.setText(new_txt)
            self._parse(new_txt)
            self.vel_plot.plot_system(self.current_session, comp)

    def _apply(self):
        self._parse(self.trans_in.text())
        self.vel_plot.plot_spectrum()

    def _parse(self, text):
        self.active_transitions = [t.strip() for t in text.split(',') if t.strip()]

    def _forward(self, name, params):
        if self.main_window: self.main_window._on_recipe_requested("absorbers", name, params, {})

    def _on_context_menu(self, pos):
        if not self.table_view.indexAt(pos).isValid(): return
        m = QMenu(self)
        m.addAction("Refit Selected", self._refit)
        m.addAction("Delete", self._delete)
        m.exec(self.table_view.mapToGlobal(pos))

    def _delete(self):
        row = self.table_view.currentIndex().row()
        comp = self.table_model.get_component_at(row)
        if comp and self.main_window:
            if QMessageBox.question(self, "Delete", f"Delete {comp.series}?", QMessageBox.Yes|QMessageBox.No) == QMessageBox.Yes:
                self.main_window._on_recipe_requested("absorbers", "delete_component", {"uuid": comp.uuid}, {})

    def _refit(self):
        row = self.table_view.currentIndex().row()
        comp = self.table_model.get_component_at(row)
        if comp and self.main_window:
            self.main_window._on_recipe_requested("absorbers", "fit_component", {"uuid": comp.uuid}, {})