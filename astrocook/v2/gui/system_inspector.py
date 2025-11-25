import logging
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.ticker as ticker
from PySide6.QtCore import Qt, QAbstractTableModel, QModelIndex, Signal, QItemSelectionModel
from PySide6.QtGui import QAction, QCursor
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QTableView, 
    QPushButton, QHeaderView, QAbstractItemView, QMenu,
    QScrollArea, QSizePolicy, QMessageBox
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

class VelocityPlotWidget(QWidget):
    def __init__(self, inspector_parent): 
        super().__init__()
        self.inspector = inspector_parent
        
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0,0,0,0)
        self.layout.setSpacing(9) 

        self.fig = Figure(figsize=(5, 8), dpi=100, constrained_layout=True) 
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True) 
        
        self.content_widget = QWidget()
        self.content_layout = QVBoxLayout(self.content_widget)
        self.content_layout.setContentsMargins(0,0,0,0)
        self.content_layout.addWidget(self.canvas)
        
        self.scroll.setWidget(self.content_widget)
        self.layout.addWidget(self.scroll, 1) 

        self.toolbar = AstrocookToolbar(self.canvas, self)
        for action in self.toolbar.actions():
            if action.text() == 'Select':
                self.toolbar.removeAction(action); break
        self.layout.addWidget(self.toolbar, 0) 

        # Internal State
        self.axes = []
        self.cursor_lines = [] 
        self.bg_cache = None
        self._current_z_sys = 0.0
        self._panel_map = {} 
        self._last_xlim = None
        self._last_ylim = None
        
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('draw_event', self.on_draw)
        self.canvas.mpl_connect('button_press_event', self.on_press)

    @property
    def main_window(self):
        return self.inspector.main_window

    def toggle_region_selector(self): pass
        
    def plot_spectrum(self, session_state=None, force_autoscale=False):
        if force_autoscale:
            self._last_xlim = None
            self._last_ylim = None
        if hasattr(self, '_last_session') and hasattr(self, '_last_component'):
            self.plot_system(self._last_session, self._last_component)

    def plot_system(self, session, component: ComponentDataV2):
        """
        Main plotting routine.
        """
        self._last_session = session
        self._last_component = component

        self._current_z_sys = component.z if component else 0.0

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

        # 2. Resize Canvas Height
        total_height = max(600, num_plots * 250) 
        self.content_widget.setMinimumHeight(total_height)
        self.canvas.setMinimumHeight(total_height)
        
        # 3. Create Subplots
        axs = self.fig.subplots(nrows=num_plots, ncols=1, sharex=True, sharey=True)
        if num_plots == 1: axs = [axs] 
        
        # 4. Prepare Data
        has_cont = session.spec.cont is not None
        y_full = session.spec.y.value
        dy_full = session.spec.dy.value if session.spec.dy is not None else np.zeros_like(y_full)
        
        if has_cont:
            cont_full = session.spec.cont.value
            y_norm = np.divide(y_full, cont_full, out=np.ones_like(y_full), where=cont_full!=0)
            dy_norm = np.divide(dy_full, cont_full, out=np.zeros_like(dy_full), where=cont_full!=0)
        else:
            y_norm = y_full 
            dy_norm = dy_full
            
        x_full = session.spec.x.value
        c_kms = 299792.458
        window_kms = 200 

        # --- GET COLORS ---
        colors = get_color_cycle(5, cmap='tab20')
        
        col_flux = get_style_color('flux', colors)
        col_model = get_style_color('model', colors)
        col_err = PLOT_STYLE['error']['color']
        col_cursor = get_style_color('cursor', colors)
        alpha_err = PLOT_STYLE['error']['alpha']
        
        style_flux = PLOT_STYLE['flux']
        style_error = PLOT_STYLE['error']
        style_model = PLOT_STYLE['model']

        # 5. Plot Loop
        for i, (ax, trans_name) in enumerate(zip(axs, transitions)):
            self.axes.append(ax)
            self._panel_map[ax] = trans_name

            # This closure captures 'self._current_z_sys'
            def make_format_coord(current_z_sys):
                def format_coord(x, y):
                    # x is velocity in km/s
                    # Formula: z_local = (1 + z_sys) * (1 + v/c) - 1
                    z_val = (1 + current_z_sys) * (1 + x / c_kms) - 1
                    return f"v={x:.1f} km/s, y={y:.2f}, z={z_val:.5f}"
                return format_coord
            
            ax.format_coord = make_format_coord(self._current_z_sys)
            
            if trans_name not in xem_d:
                ax.text(0.5, 0.5, f"Unknown: {trans_name}", ha='center', transform=ax.transAxes)
                continue

            try:
                lambda_0 = xem_d[trans_name].to(session.spec.x.unit).value
            except: continue

            lambda_obs = lambda_0 * (1.0 + component.z)
            v_full = c_kms * (x_full - lambda_obs) / lambda_obs
            
            mask = (v_full > -5000) & (v_full < 5000)
            v_plot = v_full[mask]
            y_plot = y_norm[mask]
            dy_plot = dy_norm[mask]
            
            # --- Draw Error Shading ---
            ax.step(v_plot, y_plot - dy_plot, where='mid', 
                    color=col_err, lw=style_error['lw'], alpha=alpha_err, label='Error')
            ax.step(v_plot, y_plot + dy_plot, where='mid', 
                    color=col_err, lw=style_error['lw'], alpha=alpha_err)

            # --- Draw Flux (Step Style) ---
            ax.step(v_plot, y_plot, where='mid', 
                    color=col_flux, lw=style_flux['lw'], 
                    label=style_flux['label'])
            
            # --- Draw Model ---
            if session.spec.model is not None:
                mod_full = session.spec.model.value
                if has_cont:
                    mod_plot = np.divide(mod_full[mask], cont_full[mask], out=np.ones_like(y_plot), where=cont_full[mask]!=0)
                else:
                    mod_plot = mod_full[mask]
                
                ax.plot(v_plot, mod_plot, 
                        color=col_model, lw=style_model['lw'], alpha=style_model['alpha'],
                        label=style_model['label'])

            # --- Decorations ---
            # Center Line
            ax.axvline(0, color=PLOT_STYLE['center_line']['color'], 
                       ls=PLOT_STYLE['center_line']['ls'], lw=PLOT_STYLE['center_line']['lw'], alpha=PLOT_STYLE['center_line']['alpha']) 
            
            # Continuum (1.0) and Zero (0.0)
            # (Assuming normalized plot)
            ax.axhline(1.0, color='green', ls=':', alpha=0.5, lw=0.8) 
            ax.axhline(0.0, color=PLOT_STYLE['zero_line']['color'], 
                       ls=PLOT_STYLE['zero_line']['ls'], lw=PLOT_STYLE['zero_line']['lw'], alpha=PLOT_STYLE['zero_line']['alpha']) 

            # Label
            ax.text(0.98, 0.85, trans_name, transform=ax.transAxes, 
                    ha='right', fontsize=10, fontweight='bold',
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

            ax.grid(True, linestyle=':', alpha=0.6)
            
            # Cursor Line
            line = ax.axvline(0, color=col_cursor, 
                              ls=PLOT_STYLE['cursor']['ls'], 
                              lw=PLOT_STYLE['cursor']['lw'], 
                              alpha=PLOT_STYLE['cursor']['alpha'], 
                              visible=False)
            self.cursor_lines.append(line)

        # 6. Final Formatting
        self.fig.supxlabel("Velocity (km/s)")
        self.fig.supylabel("Normalized Flux")
        axs[0].set_xlim(-window_kms, window_kms)
        axs[0].set_ylim(-0.1, 1.2) 
        self.canvas.draw()
    
    def on_draw(self, event):
        if self.axes:
            self.bg_cache = self.canvas.copy_from_bbox(self.fig.bbox)
            for line in self.cursor_lines: line.set_visible(False)

    def on_mouse_move(self, event):
        if not self.axes or not event.inaxes: return
        if self.bg_cache: self.canvas.restore_region(self.bg_cache)
        v_mouse = event.xdata
        for line in self.cursor_lines:
            line.set_xdata([v_mouse])
            line.set_visible(True)
            line.axes.draw_artist(line)
        self.canvas.blit(self.fig.bbox)

    def on_press(self, event):
        if event.button == 3 and event.inaxes in self.axes:
            ax = event.inaxes
            trans_name = self._panel_map.get(ax)
            
            if trans_name and self.inspector.main_window:
                v_click = event.xdata
                c_kms = 299792.458
                z_new = (1 + self._current_z_sys) * (1 + v_click / c_kms) - 1
                
                # --- FIX: Detect Multiplet for New Component ---
                # If clicked transition is part of a multiplet, add the MULTIPLET series name
                series_to_add = trans_name
                for group, lines in STANDARD_MULTIPLETS.items():
                    if trans_name in lines:
                        series_to_add = group # Use "CIV" instead of "CIV_1548"
                        break
                # ------------------------------------------------
                
                menu = QMenu(self)
                add_act = QAction(f"Add {series_to_add} at z={z_new:.5f}", menu)
                add_act.triggered.connect(lambda checked=False: self.inspector.main_window._on_recipe_requested(
                    "absorbers", "add_component", 
                    {'series': series_to_add, 'z': z_new}, {}
                ))
                menu.addAction(add_act)
                menu.exec(QCursor.pos())


class SystemInspector(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.current_session = None
        self.setWindowTitle("System Inspector")
        self.resize(1200, 800) 
        self.setWindowFlags(Qt.Window)
        self._setup_ui()
        
    def _setup_ui(self):
        main_layout = QVBoxLayout(self)
        
        splitter = QSplitter(Qt.Horizontal)
        
        self.table_view = QTableView()
        self.table_model = SystemTableModel()
        self.table_view.setModel(self.table_model)
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_view.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table_view.setAlternatingRowColors(True)
        self.table_view.setContextMenuPolicy(Qt.CustomContextMenu)
        
        header = self.table_view.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeToContents)
        
        self.table_view.selectionModel().selectionChanged.connect(self._on_selection_changed)
        self.table_view.customContextMenuRequested.connect(self._on_context_menu)
        self.table_model.data_changed_request.connect(self._forward_recipe_request)
        
        splitter.addWidget(self.table_view)
        
        self.vel_plot = VelocityPlotWidget(self)
        splitter.addWidget(self.vel_plot)
        
        splitter.setStretchFactor(0, 5)
        splitter.setStretchFactor(1, 5)
        main_layout.addWidget(splitter)

    def set_session(self, session):
        # 1. Store previous UUID to restore selection
        previous_uuid = None
        indexes = self.table_view.selectionModel().selectedRows()
        if indexes:
            row = indexes[0].row()
            comp = self.table_model.get_component_at(row)
            if comp: previous_uuid = comp.uuid

        self.current_session = session
        if session and session.systs:
            self.table_model.update_data(session.systs.components)
        else:
            self.table_model.update_data([])
            
        # 2. Try to restore selection
        if previous_uuid:
            found = False
            for r in range(self.table_model.rowCount()):
                c = self.table_model.get_component_at(r)
                if c.uuid == previous_uuid:
                    idx = self.table_model.index(r, 0)
                    # FIX: Force the table to acknowledge selection change
                    self.table_view.selectionModel().setCurrentIndex(idx, QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)
                    # Explicitly update plot (in case signal didn't fire if row # is same)
                    self.vel_plot.plot_system(session, c)
                    found = True
                    break
            # Fallback: If not found, select the last one? (Good for "Add" workflow)
            if not found and self.table_model.rowCount() > 0:
                 idx = self.table_model.index(self.table_model.rowCount()-1, 0)
                 self.table_view.selectionModel().setCurrentIndex(idx, QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)
                 
        elif self.table_model.rowCount() == 0:
            self.vel_plot.fig.clear()
            self.vel_plot.canvas.draw()

    def _on_selection_changed(self, selected, deselected):
        indexes = self.table_view.selectionModel().selectedRows()
        if not indexes: return
        row = indexes[0].row()
        comp = self.table_model.get_component_at(row)
        if comp and self.current_session:
            self.vel_plot.plot_system(self.current_session, comp)

    def _forward_recipe_request(self, recipe_name, params):
        if self.main_window:
            self.main_window._on_recipe_requested("absorbers", recipe_name, params, {})

    def _on_context_menu(self, pos):
        index = self.table_view.indexAt(pos)
        if not index.isValid(): return
        menu = QMenu(self)
        refit_act = QAction("Refit Selected", self)
        delete_act = QAction("Delete", self)
        refit_act.triggered.connect(self._on_refit_clicked)
        delete_act.triggered.connect(self._on_delete_clicked)
        menu.addAction(refit_act)
        menu.addSeparator()
        menu.addAction(delete_act)
        menu.exec(self.table_view.mapToGlobal(pos))

    def _on_delete_clicked(self):
        indexes = self.table_view.selectionModel().selectedRows()
        if not indexes: return
        row = indexes[0].row()
        comp = self.table_model.get_component_at(row)
        if comp and self.main_window:
            confirm = QMessageBox.question(self, "Delete Component", 
                                           f"Delete component {comp.series} at z={comp.z:.4f}?",
                                           QMessageBox.Yes | QMessageBox.No)
            if confirm == QMessageBox.Yes:
                self.main_window._on_recipe_requested(
                    "absorbers", "delete_component", {"uuid": comp.uuid}, {})

    def _on_refit_clicked(self):
        indexes = self.table_view.selectionModel().selectedRows()
        if not indexes: return
        row = indexes[0].row()
        comp = self.table_model.get_component_at(row)
        if comp and self.main_window:
            self.main_window._on_recipe_requested(
                "absorbers", "fit_component", {"uuid": comp.uuid}, {})