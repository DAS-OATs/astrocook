import logging
import numpy as np
import astropy.units as au
from scipy.special import wofz 
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.ticker as ticker
from PySide6.QtCore import Qt, QAbstractTableModel, QModelIndex, Signal, QSortFilterProxyModel, QLocale
from PySide6.QtGui import QAction, QCursor, QDoubleValidator
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QTableView, 
    QPushButton, QHeaderView, QAbstractItemView, QMenu,
    QScrollArea, QSizePolicy, QMessageBox, QLabel, QLineEdit, QLayout, QFrame, QScrollBar, QCheckBox
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

# --- Proxy Model ---
class SystemSortFilterProxyModel(QSortFilterProxyModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.group_filtering_enabled = False
        self.target_group_id = None

    def set_group_filter(self, enabled: bool, group_id: int = None):
        self.group_filtering_enabled = enabled
        self.target_group_id = group_id
        self.invalidateFilter()

    def filterAcceptsRow(self, source_row, source_parent):
        if not self.group_filtering_enabled or self.target_group_id is None: return True
        model = self.sourceModel()
        comp = model.get_component_at(source_row)
        return True if comp and comp.id == self.target_group_id else False

    def lessThan(self, left, right):
        l_data = self.sourceModel().data(left, Qt.DisplayRole)
        r_data = self.sourceModel().data(right, Qt.DisplayRole)
        try: return float(l_data) < float(r_data)
        except (ValueError, TypeError): return str(l_data) < str(r_data)

# --- Helper for Voigt Profile ---
def calc_voigt_profile(wave_grid_ang, lambda_0, f_val, gamma, z, N, b_kms):
    if wave_grid_ang is None or len(wave_grid_ang) == 0: return np.array([])
    lambda_c = lambda_0 * (1.0 + z)
    c_kms = 2.99792458e5
    b_safe = max(b_kms, 0.1)
    dop_width_ang = (b_safe / c_kms) * lambda_c
    x = (wave_grid_ang - lambda_c) / dop_width_ang
    c_ang_s = 2.99792458e18
    a = (gamma * lambda_c**2) / (4.0 * np.pi * c_ang_s * dop_width_ang)
    H_ax = wofz(x + 1j * a).real
    tau = 1.4974e-15 * N * f_val * lambda_0 / b_safe * H_ax
    return np.exp(-tau)

# --- 2. The Paged Plot Widget ---

class VelocityPlotWidget(QWidget):
    def __init__(self, inspector_parent): 
        super().__init__()
        self.inspector = inspector_parent
        self.PAGE_SIZE = 3 
        self._xlim = (-300, 300)
        self.cursor_logN = 13.5
        self.cursor_b = 10.0
        
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
        self.axes = [] 
        
        self._current_session = None
        self._current_component = None
        self._selected_components = []
        self._plot_center_z = 0.0 # NEW: Controls the 0 km/s reference

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

    def update_cursor_params(self, logN, b):
        self.cursor_logN = logN
        self.cursor_b = b
        self.canvas.draw_idle()

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

    def set_center_redshift(self, z):
        self._plot_center_z = z
        # Reset velocity limits to default window around new center
        width = self._xlim[1] - self._xlim[0]
        self._xlim = (-width/2, width/2)
        # Update UI to reflect the reset limits (optional, but good for sync)
        self.inspector.update_limit_boxes(self._xlim[0], self._xlim[1])
        self._update_plot()

    def plot_system(self, session, component: ComponentDataV2, selected_components: List[ComponentDataV2] = None):
        is_new_component = (self._current_component != component)
        
        self._current_session = session
        self._current_component = component
        self._selected_components = selected_components if selected_components else [component]
        
        self._plot_center_z = component.z

        # If it's a new component, recenter the view on 0 km/s
        # while keeping the user's preferred zoom width.
        if is_new_component:
            width = self._xlim[1] - self._xlim[0]
            self._xlim = (-width/2, width/2)
            # Update the UI text boxes to match
            self.inspector.update_limit_boxes(self._xlim[0], self._xlim[1])

        if not session or not component:
            self._all_transitions = []
            self._update_plot()
            return

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

        if abs(xlim[1] - xlim[0]) < 10: return

        if not np.allclose(xlim, self._xlim, atol=0.1):
            self._xlim = xlim
            self.inspector.update_limit_boxes(xlim[0], xlim[1])
            # Back-calculate Z based on the center of the current view relative to _plot_center_z
            c_kms = 299792.458
            v_center = np.mean(xlim)
            z_view_center = (1 + self._plot_center_z) * (1 + v_center / c_kms) - 1
            self.inspector.update_z_box(z_view_center)

    def _update_plot(self):
        self.fig.clear()
        self.axes = []
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
        z_sys = self._plot_center_z
        
        has_cont = session.spec.cont is not None
        y = session.spec.y.value
        dy = session.spec.dy.value if session.spec.dy is not None else np.zeros_like(y)
        cont = session.spec.cont.value if has_cont else np.ones_like(y)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            y_norm = np.divide(y, cont, where=cont!=0)
            dy_norm = np.divide(dy, cont, where=cont!=0)
        
        show_resid = self.inspector.resid_cb.isChecked()
        y_resid = None
        if show_resid and session.spec.model is not None:
            mod = session.spec.model.value
            with np.errstate(divide='ignore', invalid='ignore'):
                mod_norm = np.divide(mod, cont, where=cont!=0) if has_cont else mod
                y_resid = y_norm - mod_norm
        
        colors = get_color_cycle(5, cmap='tab20')

        for i, trans_name in enumerate(visible_trans):
            ax = axs[i]
            self.axes.append(ax)
            
            cursor_line, = ax.plot([], [], color='purple', lw=1.5, alpha=0.8, zorder=10)
            self._panel_map[ax] = (trans_name, cursor_line)

            def make_fmt(z):
                def fmt(x, y):
                    val_z = (1 + z) * (1 + x / c_kms) - 1
                    return f"v={x:.1f}, y={y:.2f}, z={val_z:.5f}"
                return fmt
            ax.format_coord = make_fmt(z_sys)

            if trans_name not in xem_d:
                ax.text(0.5, 0.5, f"Unknown: {trans_name}", ha='center'); continue
            try: 
                lam_0 = xem_d[trans_name].to(session.spec.x.unit).value
            except: 
                ax.text(0.5, 0.5, f"Unit Error", ha='center'); continue

            lam_obs = lam_0 * (1.0 + z_sys)
            v = c_kms * (x_full - lam_obs) / lam_obs
            
            v_min, v_max = self._xlim
            mask = (v > v_min - 500) & (v < v_max + 500)
            v_p = v[mask]; y_p = y_norm[mask]; dy_p = dy_norm[mask]
            
            if len(v_p) == 0:
                ax.text(0.5, 0.5, "No Data", ha='center', transform=ax.transAxes)
                continue

            ax.step(v_p, y_p-dy_p, where='mid', color=PLOT_STYLE['error']['color'], lw=0.5, alpha=0.5)
            ax.step(v_p, y_p+dy_p, where='mid', color=PLOT_STYLE['error']['color'], lw=0.5, alpha=0.5)
            ax.step(v_p, y_p, where='mid', color=get_style_color('flux', colors), lw=0.8)

            if session.spec.model is not None:
                mod = session.spec.model.value
                with np.errstate(divide='ignore', invalid='ignore'):
                    mod_p = np.divide(mod[mask], cont[mask], where=cont[mask]!=0) if has_cont else mod[mask]
                ax.plot(v_p, mod_p, color=get_style_color('model', colors), lw=1.0)

            if show_resid and y_resid is not None:
                resid_p = y_resid[mask]
                offset = -0.2
                ax.axhline(offset, color='gray', ls=':', lw=0.8, alpha=0.5)
                ax.step(v_p, offset + resid_p - dy_p, where='mid', color='#cccccc', lw=0.3, alpha=0.4)
                ax.step(v_p, offset + resid_p + dy_p, where='mid', color='#cccccc', lw=0.3, alpha=0.4)
                ax.step(v_p, offset + resid_p, where='mid', color='black', lw=0.6, alpha=0.7)

            atom_info = ATOM_DATA.get(trans_name)
            if atom_info and self._selected_components:
                x_ang_p = session.spec.x.to(au.Angstrom).value[mask]
                for sel_c in self._selected_components:
                    is_match = (sel_c.series == trans_name) or \
                               (sel_c.series in STANDARD_MULTIPLETS and trans_name in STANDARD_MULTIPLETS[sel_c.series])
                    if not is_match: continue 
                    prof = calc_voigt_profile(
                        x_ang_p, atom_info['wave'], atom_info['f'], atom_info['gamma'],
                        sel_c.z, 10**sel_c.logN, sel_c.b
                    )
                    ax.plot(v_p, prof, color='orange', ls='--', lw=1.2, alpha=0.9)

            tick_ymin, tick_ymax = 0.02, 0.08 
            trans_axis = ax.get_xaxis_transform()
            for other_c in all_comps:
                v_shift = c_kms * (other_c.z - z_sys) / (1.0 + z_sys)
                if v_min <= v_shift <= v_max:
                    is_h = other_c.series.startswith('Ly') or other_c.series.startswith('H')
                    col = 'red' if is_h else 'gray'
                    alpha = 0.8 if is_h else 0.5
                    lw = 1.5 if is_h else 1.0
                    zorder = 5 if is_h else 4
                    ax.plot([v_shift, v_shift], [tick_ymin, tick_ymax], 
                            transform=trans_axis, color=col, lw=lw, alpha=alpha, zorder=zorder)

            ax.axvline(0, color='gray', ls='--', lw=0.8)
            ax.axhline(1.0, color='green', ls=':', alpha=0.5)
            ax.axhline(0.0, color='gray', lw=0.5)
            ax.text(0.98, 0.85, trans_name, transform=ax.transAxes, ha='right', fontweight='bold', 
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            
            ax.set_xlim(self._xlim)
            ax.set_ylim(-0.4 if show_resid else -0.2, 1.4)

        self.fig.supxlabel("Velocity (km/s)")
        self.fig.supylabel("Normalized Flux")
        self.canvas.draw()

    def on_mouse_move(self, event):
        if not event.inaxes or not self._current_component or not self._current_session: return
        
        v_mouse = event.xdata
        c_kms = 299792.458
        z_new = (1 + self._plot_center_z) * (1 + v_mouse / c_kms) - 1
        
        N = 10**self.cursor_logN
        b = self.cursor_b
        
        show_resid = self.inspector.resid_cb.isChecked()

        for ax in self.axes:
            if ax not in self._panel_map: continue
            trans_name, line_artist = self._panel_map[ax]
            # [FIX] Ensure the profile calculation also uses the correct center logic
            # The profile function takes z_new (calculated above), so it is correct.
            # However, we must ensure lam_obs_sys is defined relative to the PLOT CENTER
            # so the grid aligns.
            
            # Actually, atom_info['wave'] * (1+z) gives the obs wavelength of the line.
            # The grid 'v_grid' is relative to z_sys (_plot_center_z).
            # So lam_grid = lam_0 * (1+z_sys) * (1 + v/c).
            # This logic below needs to match _update_plot's grid definition:
            
            atom_info = ATOM_DATA.get(trans_name)
            if not atom_info: continue

            lam_0 = atom_info['wave']
            lam_obs_center = lam_0 * (1 + self._plot_center_z) # [FIX] Center on plot Z
            
            v_min, v_max = ax.get_xlim()
            v_grid = np.linspace(v_min, v_max, 300)
            
            # Grid in Angstroms corresponding to the plot's V pixels
            lam_grid = lam_obs_center * (1 + v_grid / c_kms)
            
            prof = calc_voigt_profile(
                lam_grid, lam_0, atom_info['f'], atom_info['gamma'],
                z_new, N, b
            )
            
            if show_resid:
                offset = -0.2
                resid_prof = offset - (1.0 - prof)
                line_artist.set_data(v_grid, resid_prof)
            else:
                line_artist.set_data(v_grid, prof)
            
        self.canvas.draw_idle()

    def on_press(self, event):
        if event.button == 3 and event.inaxes:
            if event.inaxes not in self._panel_map: return
            trans_name, _ = self._panel_map[event.inaxes]
            
            if self.inspector.main_window:
                v = event.xdata
                c_kms = 299792.458
                z_new = (1 + self._plot_center_z) * (1 + v / c_kms) - 1

                series = trans_name
                for g, l in STANDARD_MULTIPLETS.items():
                    if trans_name in l: series = g; break
                
                menu = QMenu(self)
                act = QAction(f"Add {series} at z={z_new:.5f}", menu)
                act.triggered.connect(lambda: self.inspector.main_window._on_recipe_requested(
                    "absorbers", "add_component", {
                        'series': series, 'z': z_new,
                        'logN': self.cursor_logN, 'b': self.cursor_b
                    }, {}))
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
        self._setup_ui()
        
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        
        # Left: Table
        self.table_view = QTableView()
        self.table_model = SystemTableModel()
        
        self.proxy_model = SystemSortFilterProxyModel()
        self.proxy_model.setSourceModel(self.table_model)
        
        self.table_view.setModel(self.proxy_model)
        self.table_view.setSortingEnabled(True)
        
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_view.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.table_view.setAlternatingRowColors(True)
        self.table_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table_view.selectionModel().selectionChanged.connect(self._on_selection_changed)
        self.table_view.customContextMenuRequested.connect(self._on_context_menu)
        self.table_model.data_changed_request.connect(self._forward)
        
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0,0,0,0)
        
        self.group_cb = QCheckBox("Group View (Same ID)")
        self.group_cb.toggled.connect(self._on_group_toggled)
        left_layout.addWidget(self.group_cb)
        left_layout.addWidget(self.table_view)
        
        splitter.addWidget(left_widget)
        
        # Right Side
        right = QWidget()
        r_layout = QVBoxLayout(right)
        r_layout.setContentsMargins(0,0,0,0)
        
        # --- Toolbar Controls ---
        c_layout = QHBoxLayout()
        c_layout.setContentsMargins(5,5,5,0)
        
        c_layout.addWidget(QLabel("Trans:"))
        self.trans_in = QLineEdit()
        self.trans_in.setPlaceholderText("CIV, SiIV")
        self.trans_in.returnPressed.connect(self._apply)
        c_layout.addWidget(self.trans_in, 1) 
        
        c_layout.addWidget(QLabel("z:"))
        self.z_in = QLineEdit()
        self.z_in.setFixedWidth(70)
        validator = QDoubleValidator()
        validator.setLocale(QLocale.C) # [FIX] Enforce dot decimal
        self.z_in.setValidator(validator)
        self.z_in.returnPressed.connect(self._apply_z)
        c_layout.addWidget(self.z_in)

        c_layout.addWidget(QLabel("V:"))
        self.vmin_in = QLineEdit("-300")
        self.vmin_in.setFixedWidth(50)
        self.vmin_in.setValidator(validator) # Use C locale
        self.vmin_in.returnPressed.connect(self._apply_limits)
        c_layout.addWidget(self.vmin_in)
        self.vmax_in = QLineEdit("300")
        self.vmax_in.setFixedWidth(50)
        self.vmax_in.setValidator(validator) # Use C locale
        self.vmax_in.returnPressed.connect(self._apply_limits)
        c_layout.addWidget(self.vmax_in)
        
        c_layout.addWidget(QLabel("N:"))
        self.logn_in = QLineEdit("13.5")
        self.logn_in.setFixedWidth(40)
        self.logn_in.setValidator(validator) # Use C locale
        self.logn_in.textChanged.connect(self._update_cursor)
        c_layout.addWidget(self.logn_in)
        
        c_layout.addWidget(QLabel("b:"))
        self.b_in = QLineEdit("10.0")
        self.b_in.setFixedWidth(40)
        self.b_in.setValidator(validator) # Use C locale
        self.b_in.textChanged.connect(self._update_cursor)
        c_layout.addWidget(self.b_in)
        
        self.resid_cb = QCheckBox("Resid")
        self.resid_cb.toggled.connect(self._apply) 
        c_layout.addWidget(self.resid_cb)
        
        pal = QApplication.palette()
        style = f"QLineEdit {{ padding: 3px; border-radius: 4px; background: {pal.color(pal.ColorRole.Base).name()}; color: {pal.color(pal.ColorRole.Text).name()}; }}"
        for w in [self.trans_in, self.z_in, self.vmin_in, self.vmax_in, self.logn_in, self.b_in]: w.setStyleSheet(style)
        
        r_layout.addLayout(c_layout)
        
        self.vel_plot = VelocityPlotWidget(self)
        r_layout.addWidget(self.vel_plot, 1) 
        
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 4); splitter.setStretchFactor(1, 6)
        layout.addWidget(splitter)

    def set_session(self, session):
        sel = self.table_view.selectionModel().selectedRows()
        prev_uuid = None
        if sel:
            idx_source = self.proxy_model.mapToSource(sel[0])
            prev_uuid = self.table_model.get_component_at(idx_source.row()).uuid

        self.current_session = session
        self.table_model.update_data(session.systs.components if session and session.systs else [])
        
        if prev_uuid:
            for r in range(self.table_model.rowCount()):
                if self.table_model.get_component_at(r).uuid == prev_uuid:
                    idx_proxy = self.proxy_model.mapFromSource(self.table_model.index(r, 0))
                    self.table_view.selectRow(idx_proxy.row())
                    return
        if self.proxy_model.rowCount() > 0:
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

    def update_z_box(self, z_val):
        self.z_in.blockSignals(True)
        self.z_in.setText(f"{z_val:.5f}")
        self.z_in.blockSignals(False)

    def _update_cursor(self):
        try:
            n = float(self.logn_in.text())
            b = float(self.b_in.text())
            self.vel_plot.update_cursor_params(n, b)
        except ValueError: pass

    def _apply_limits(self):
        try:
            vmin = float(self.vmin_in.text())
            vmax = float(self.vmax_in.text())
            if vmin < vmax: self.vel_plot.set_velocity_limits(vmin, vmax)
        except ValueError: pass

    def _apply_z(self):
        try:
            z_target = float(self.z_in.text())
            self.vel_plot.set_center_redshift(z_target)
        except ValueError as e:
            logging.error(f"Error applying z: {e}")

    def _on_group_toggled(self, checked):
        if checked:
            sel = self.table_view.selectionModel().selectedRows()
            if sel:
                idx_src = self.proxy_model.mapToSource(sel[0])
                comp = self.table_model.get_component_at(idx_src.row())
                if comp:
                    self.proxy_model.set_group_filter(True, comp.id)
        else:
            self.proxy_model.set_group_filter(False)

    def _on_selection_changed(self, sel, desel):
        indexes = self.table_view.selectionModel().selectedRows()
        if not indexes: return
        
        idx_proxy = indexes[0]
        idx_source = self.proxy_model.mapToSource(idx_proxy)
        
        primary_row = idx_source.row()
        primary_comp = self.table_model.get_component_at(primary_row)
        
        selected_comps = []
        for idx in indexes:
            src = self.proxy_model.mapToSource(idx)
            c = self.table_model.get_component_at(src.row())
            if c: selected_comps.append(c)

        if primary_comp and self.current_session:
            cur = self.trans_in.text()
            toks = [t.strip() for t in cur.split(',') if t.strip()]
            if primary_comp.series and primary_comp.series not in toks: toks.append(primary_comp.series)
            new_txt = ", ".join(toks)
            self.trans_in.setText(new_txt)
            self.z_in.setText(f"{primary_comp.z:.5f}")
            self._parse(new_txt)
            
            if self.group_cb.isChecked():
                self.proxy_model.set_group_filter(True, primary_comp.id)
            
            self.vel_plot.plot_system(self.current_session, primary_comp, selected_comps)

    def _apply(self):
        self._parse(self.trans_in.text())
        self._apply_limits()
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
        idx = self.table_view.currentIndex()
        src = self.proxy_model.mapToSource(idx)
        comp = self.table_model.get_component_at(src.row())
        if comp and self.main_window:
            if QMessageBox.question(self, "Delete", f"Delete {comp.series}?", QMessageBox.Yes|QMessageBox.No) == QMessageBox.Yes:
                self.main_window._on_recipe_requested("absorbers", "delete_component", {"uuid": comp.uuid}, {})

    def _refit(self):
        idx = self.table_view.currentIndex()
        src = self.proxy_model.mapToSource(idx)
        comp = self.table_model.get_component_at(src.row())
        if comp and self.main_window:
            self.main_window._on_recipe_requested("absorbers", "fit_component", {"uuid": comp.uuid}, {})