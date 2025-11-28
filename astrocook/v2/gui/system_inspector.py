import logging
import numpy as np
import astropy.units as au
from scipy.special import wofz 
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.ticker as ticker
from PySide6.QtCore import Qt, QAbstractTableModel, QModelIndex, Signal, QItemSelectionModel, QEvent
from PySide6.QtGui import QAction, QCursor, QResizeEvent, QDoubleValidator
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QTableView, 
    QPushButton, QHeaderView, QAbstractItemView, QMenu,
    QScrollArea, QSizePolicy, QMessageBox, QLabel, QLineEdit, QLayout, QFrame, QScrollBar
)
from typing import List, Optional

from ..structures import ComponentDataV2
from ..atomic_data import STANDARD_MULTIPLETS, xem_d, ATOM_DATA
from .pyside_plot import AstrocookToolbar, get_color_cycle, PLOT_STYLE, get_style_color

# --- 1. Table Configuration ---
COL_MAP = {
    "ID": "id", "Transitions": "series", "z": "z", "logN": "logN", 
    "b": "b", "btur": "btur", "dz": "dz", "dlogN": "dlogN", "db": "db"
}
COLUMNS = list(COL_MAP.keys())
EDITABLE_COLS = {"Transitions", "z", "logN", "b", "btur"}

class SystemTableModel(QAbstractTableModel):
    data_changed_request = Signal(str, dict)
    def __init__(self, components: List[ComponentDataV2] = None):
        super().__init__()
        self._components = components if components else []
    def update_data(self, components: List[ComponentDataV2]):
        self.beginResetModel(); self._components = components; self.endResetModel()
    def rowCount(self, parent=QModelIndex()): return len(self._components)
    def columnCount(self, parent=QModelIndex()): return len(COLUMNS)
    def flags(self, index):
        return super().flags(index) | Qt.ItemIsEditable if COLUMNS[index.column()] in EDITABLE_COLS else super().flags(index)
    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid(): return None
        comp = self._components[index.row()]
        attr = COL_MAP[COLUMNS[index.column()]]
        if role in (Qt.DisplayRole, Qt.EditRole):
            val = getattr(comp, attr, None)
            if role == Qt.DisplayRole and isinstance(val, float):
                return f"{val:.6f}" if attr == 'z' else f"{val:.3f}" if attr in ['logN', 'dlogN'] else f"{val:.2f}"
            return str(val) if val is not None else ""
        return Qt.AlignCenter if role == Qt.TextAlignmentRole else None
    def setData(self, index, value, role=Qt.EditRole):
        if not index.isValid() or role != Qt.EditRole: return False
        attr = COL_MAP[COLUMNS[index.column()]]
        try:
            val = float(value) if attr in ['z', 'logN', 'b', 'btur'] else str(value)
            self.data_changed_request.emit("update_component", {attr: val, 'uuid': self._components[index.row()].uuid})
            return True
        except ValueError: return False
    def headerData(self, section, orientation, role):
        return COLUMNS[section] if role == Qt.DisplayRole and orientation == Qt.Horizontal else None
    def get_component_at(self, row: int): return self._components[row] if 0 <= row < len(self._components) else None

# --- Helper for Voigt Profile ---
def calc_voigt_profile(wave_grid_ang, lambda_0, f_val, gamma, z, N, b_kms):
    """Pure logic to generate a Voigt profile for visualization."""
    lambda_c = lambda_0 * (1.0 + z)
    c_kms = 2.99792458e5
    b_safe = max(b_kms, 0.1)
    dop_width_ang = (b_safe / c_kms) * lambda_c
    
    x = (wave_grid_ang - lambda_c) / dop_width_ang
    c_ang_s = 2.99792458e18
    a = (gamma * lambda_c**2) / (4.0 * np.pi * c_ang_s * dop_width_ang)
    H_ax = wofz(x + 1j * a).real
    
    tau_factor = 1.4974e-15
    tau = tau_factor * N * f_val * lambda_0 / b_safe * H_ax
    return np.exp(-tau)

# --- 2. The Paged Plot Widget ---

class VelocityPlotWidget(QWidget):
    def __init__(self, inspector_parent): 
        super().__init__()
        self.inspector = inspector_parent
        self.PAGE_SIZE = 3 
        self._xlim = (-300, 300)
        
        self.main_layout = QHBoxLayout(self)
        self.main_layout.setContentsMargins(0, 0, 0, 0)
        self.main_layout.setSpacing(2)

        # Left Side
        left_container = QWidget()
        self.left_layout = QVBoxLayout(left_container)
        self.left_layout.setContentsMargins(0, 0, 0, 0)
        self.left_layout.setSpacing(0)

        self.fig = Figure(figsize=(5, 6), dpi=100)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        self.toolbar = AstrocookToolbar(self.canvas, self)
        for action in self.toolbar.actions():
            if action.text() == 'Select': self.toolbar.removeAction(action); break
        
        self.left_layout.addWidget(self.canvas)
        self.left_layout.addWidget(self.toolbar)

        # Right Side
        self.scrollbar = QScrollBar(Qt.Vertical)
        self.scrollbar.valueChanged.connect(self._on_scroll)
        self.scrollbar.setRange(0, 0)

        self.main_layout.addWidget(left_container)
        self.main_layout.addWidget(self.scrollbar)

        # State
        self._all_transitions = [] 
        self._panel_map = {} 
        self.cursor_lines = []
        
        self._current_session = None
        self._current_component = None
        self._selected_components = []
        
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('draw_event', self.on_draw_update_limits)

    @property
    def main_window(self):
        return self.inspector.main_window

    def toggle_region_selector(self): pass

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.PAGE_SIZE = max(1, event.size().height() // 200)
        self._update_scrollbar_range()
        self._update_plot()

    def plot_spectrum(self, session_state=None, force_autoscale=False):
        session = session_state if session_state else self._current_session
        comp = self._current_component
        if session and comp:
            self.plot_system(session, comp, self._selected_components)
        else:
            self.fig.clear()
            self.canvas.draw()

    def set_velocity_limits(self, v_min, v_max):
        self._xlim = (v_min, v_max)
        self._update_plot()

    def plot_system(self, session, component: ComponentDataV2, selected_components: List[ComponentDataV2] = None):
        self._current_session = session
        self._current_component = component
        self._selected_components = selected_components if selected_components else [component]
        
        if not session or not component:
            self._all_transitions = []
            self._update_plot()
            return

        # 1. Transitions
        expanded = []
        for item in self.inspector.active_transitions:
            if item in STANDARD_MULTIPLETS: expanded.extend(STANDARD_MULTIPLETS[item])
            elif item in xem_d: expanded.append(item)
        
        seen = set()
        self._all_transitions = [x for x in expanded if not (x in seen or seen.add(x))]

        self._update_scrollbar_range()
        self._update_plot()

    def _update_scrollbar_range(self):
        total = len(self._all_transitions)
        if total <= self.PAGE_SIZE:
            self.scrollbar.setRange(0, 0); self.scrollbar.setEnabled(False)
        else:
            self.scrollbar.setRange(0, total - self.PAGE_SIZE); self.scrollbar.setEnabled(True)

    def _on_scroll(self, value): self._update_plot()

    def on_draw_update_limits(self, event):
        if not self.axes: return
        xlim = self.axes[0].get_xlim()
        if not np.allclose(xlim, self._xlim, atol=0.1):
            self._xlim = xlim
            self.inspector.update_limit_boxes(xlim[0], xlim[1])

    def _update_plot(self):
        self.fig.clear()
        self.axes = []
        self.cursor_lines = []
        self._panel_map = {}

        if not self._all_transitions or not self._current_session:
            self.canvas.draw(); return

        start = self.scrollbar.value() if self.scrollbar.isEnabled() else 0
        end = min(start + self.PAGE_SIZE, len(self._all_transitions))
        visible_trans = self._all_transitions[start:end]
        
        num_plots = len(visible_trans)
        if num_plots == 0: self.canvas.draw(); return

        self.fig.set_layout_engine(layout='constrained')
        axs = self.fig.subplots(nrows=num_plots, ncols=1, sharex=True, sharey=True)
        if num_plots == 1: axs = [axs]
        else: axs = axs.flatten()

        session = self._current_session
        comp = self._current_component
        all_comps = session.systs.components
        
        x_full = session.spec.x.value
        c_kms = 299792.458
        z_sys = comp.z
        
        has_cont = session.spec.cont is not None
        y = session.spec.y.value
        dy = session.spec.dy.value if session.spec.dy is not None else np.zeros_like(y)
        cont = session.spec.cont.value if has_cont else np.ones_like(y)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            y_norm = np.divide(y, cont, where=cont!=0)
            dy_norm = np.divide(dy, cont, where=cont!=0)
        
        colors = get_color_cycle(5, cmap='tab20')

        for i, trans_name in enumerate(visible_trans):
            ax = axs[i]
            self.axes.append(ax)
            self._panel_map[ax] = trans_name

            # 1. Formatters
            def make_fmt(z):
                def fmt(x, y):
                    val_z = (1 + z) * (1 + x / c_kms) - 1
                    return f"v={x:.1f}, y={y:.2f}, z={val_z:.5f}"
                return fmt
            ax.format_coord = make_fmt(z_sys)

            # 2. Physics check
            if trans_name not in xem_d:
                ax.text(0.5, 0.5, f"Unknown: {trans_name}", ha='center'); continue
            try: 
                lam_0 = xem_d[trans_name].to(session.spec.x.unit).value
            except: 
                ax.text(0.5, 0.5, f"Unit Error", ha='center'); continue

            lam_obs = lam_0 * (1.0 + z_sys)
            v = c_kms * (x_full - lam_obs) / lam_obs
            
            # Mask
            v_min, v_max = self._xlim
            mask = (v > v_min - 500) & (v < v_max + 500)
            v_p = v[mask]; y_p = y_norm[mask]; dy_p = dy_norm[mask]
            
            # 3. Plot Data
            ax.step(v_p, y_p-dy_p, where='mid', color=PLOT_STYLE['error']['color'], lw=0.5, alpha=0.5)
            ax.step(v_p, y_p+dy_p, where='mid', color=PLOT_STYLE['error']['color'], lw=0.5, alpha=0.5)
            ax.step(v_p, y_p, where='mid', color=get_style_color('flux', colors), lw=0.8)

            # 4. Plot Global Model
            if session.spec.model is not None:
                mod = session.spec.model.value
                with np.errstate(divide='ignore', invalid='ignore'):
                    mod_p = np.divide(mod[mask], cont[mask], where=cont[mask]!=0)
                ax.plot(v_p, mod_p, color=get_style_color('model', colors), lw=1.0)

            # --- 5. Plot Individual Selected Components (Corrected) ---
            atom_info = ATOM_DATA.get(trans_name)
            if atom_info and self._selected_components:
                x_ang_p = session.spec.x.to(au.Angstrom).value[mask]
                
                for sel_c in self._selected_components:
                    # FIX: Check if this component actually corresponds to this panel's transition
                    is_match = False
                    if sel_c.series == trans_name:
                        is_match = True
                    elif sel_c.series in STANDARD_MULTIPLETS and trans_name in STANDARD_MULTIPLETS[sel_c.series]:
                        is_match = True
                    
                    if not is_match: continue # Skip if this component isn't relevant to this panel

                    prof = calc_voigt_profile(
                        x_ang_p, atom_info['wave'], atom_info['f'], atom_info['gamma'],
                        sel_c.z, 10**sel_c.logN, sel_c.b
                    )
                    ax.plot(v_p, prof, color='orange', ls='--', lw=1.2, alpha=0.9)

            # --- 6. Plot System Markers ---
            tick_ymin, tick_ymax = 0.02, 0.08 
            trans = ax.get_xaxis_transform()
            
            for other_c in all_comps:
                v_shift = c_kms * (other_c.z - z_sys) / (1.0 + z_sys)
                if v_min <= v_shift <= v_max:
                    is_h = other_c.series.startswith('Ly') or other_c.series.startswith('H')
                    col = 'red' if is_h else 'gray'
                    alpha = 0.8 if is_h else 0.5
                    lw = 1.5 if is_h else 1.0
                    zorder = 5 if is_h else 4
                    ax.plot([v_shift, v_shift], [tick_ymin, tick_ymax], 
                            transform=trans, color=col, lw=lw, alpha=alpha, zorder=zorder)

            # 7. Decorations
            ax.axvline(0, color='gray', ls='--', lw=0.8)
            ax.axhline(1.0, color='green', ls=':', alpha=0.5)
            ax.axhline(0.0, color='gray', lw=0.5)
            ax.text(0.98, 0.85, trans_name, transform=ax.transAxes, ha='right', fontweight='bold', 
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            
            ax.set_xlim(self._xlim)
            ax.set_ylim(-0.2, 1.4)
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
                
                series = trans
                for g, l in STANDARD_MULTIPLETS.items():
                    if trans in l: series = g; break
                
                menu = QMenu(self)
                act = QAction(f"Add {series} at z={z_new:.5f}", menu)
                act.triggered.connect(lambda: self.inspector.main_window._on_recipe_requested(
                    "absorbers", "add_component", {'series': series, 'z': z_new}, {}))
                menu.addAction(act)
                menu.exec(QCursor.pos())

# --- 3. The Inspector Window ---

class SystemInspector(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.current_session = None
        self.active_transitions: List[str] = [] 
        self.setWindowTitle("System Inspector")
        self.resize(1200, 800) 
        self.setWindowFlags(Qt.Window)
        self._setup_ui()
        
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        
        # Table
        self.table_view = QTableView()
        self.table_model = SystemTableModel()
        self.table_view.setModel(self.table_model)
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_view.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.table_view.setAlternatingRowColors(True)
        self.table_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table_view.selectionModel().selectionChanged.connect(self._on_selection_changed)
        self.table_view.customContextMenuRequested.connect(self._on_context_menu)
        self.table_model.data_changed_request.connect(self._forward_recipe_request)
        splitter.addWidget(self.table_view)
        
        # Right Side
        right = QWidget()
        r_layout = QVBoxLayout(right)
        r_layout.setContentsMargins(0,0,0,0)
        
        c_layout = QHBoxLayout()
        c_layout.setContentsMargins(5,5,5,0)
        
        c_layout.addWidget(QLabel("Shown transitions:"))
        self.trans_in = QLineEdit()
        self.trans_in.setPlaceholderText("e.g. CIV, SiIV")
        self.trans_in.returnPressed.connect(self._apply)
        c_layout.addWidget(self.trans_in)
        
        c_layout.addSpacing(15)
        c_layout.addWidget(QLabel("V Min:"))
        self.vmin_in = QLineEdit("-300")
        self.vmin_in.setFixedWidth(60)
        self.vmin_in.setValidator(QDoubleValidator())
        self.vmin_in.returnPressed.connect(self._apply_limits)
        c_layout.addWidget(self.vmin_in)
        
        c_layout.addWidget(QLabel("V Max:"))
        self.vmax_in = QLineEdit("300")
        self.vmax_in.setFixedWidth(60)
        self.vmax_in.setValidator(QDoubleValidator())
        self.vmax_in.returnPressed.connect(self._apply_limits)
        c_layout.addWidget(self.vmax_in)
        
        pal = QApplication.palette()
        style = f"QLineEdit {{ padding: 3px; border-radius: 4px; background: {pal.color(pal.ColorRole.Base).name()}; color: {pal.color(pal.ColorRole.Text).name()}; }}"
        self.trans_in.setStyleSheet(style)
        self.vmin_in.setStyleSheet(style)
        self.vmax_in.setStyleSheet(style)
        
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
                    return
        if self.table_model.rowCount() > 0:
            self.table_view.selectRow(0)
        else:
            self.vel_plot.plot_spectrum(None)

    def update_limit_boxes(self, v_min, v_max):
        self.vmin_in.blockSignals(True)
        self.vmax_in.blockSignals(True)
        self.vmin_in.setText(f"{v_min:.0f}")
        self.vmax_in.setText(f"{v_max:.0f}")
        self.vmin_in.blockSignals(False)
        self.vmax_in.blockSignals(False)

    def _apply_limits(self):
        try:
            vmin = float(self.vmin_in.text())
            vmax = float(self.vmax_in.text())
            if vmin < vmax: self.vel_plot.set_velocity_limits(vmin, vmax)
        except ValueError: pass

    def _on_selection_changed(self, sel, desel):
        indexes = self.table_view.selectionModel().selectedRows()
        if not indexes: return
        
        primary_row = indexes[0].row()
        primary_comp = self.table_model.get_component_at(primary_row)
        
        selected_comps = []
        for idx in indexes:
            c = self.table_model.get_component_at(idx.row())
            if c: selected_comps.append(c)

        if primary_comp and self.current_session:
            cur = self.trans_in.text()
            toks = [t.strip() for t in cur.split(',') if t.strip()]
            if primary_comp.series and primary_comp.series not in toks: toks.append(primary_comp.series)
            new_txt = ", ".join(toks)
            self.trans_in.setText(new_txt)
            self._parse(new_txt)
            self.vel_plot.plot_system(self.current_session, primary_comp, selected_comps)

    def _apply(self):
        self._parse(self.trans_in.text())
        self._apply_limits()
        self.vel_plot.plot_spectrum()

    def _parse(self, text):
        self.active_transitions = [t.strip() for t in text.split(',') if t.strip()]

    def _forward_recipe_request(self, recipe_name, params):
        if self.main_window: self.main_window._on_recipe_requested("absorbers", recipe_name, params, {})

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