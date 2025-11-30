import logging
import numpy as np
import astropy.units as au
from scipy.special import wofz 
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from PySide6.QtCore import Qt, QAbstractTableModel, QModelIndex, Signal, QSortFilterProxyModel, QLocale, QPoint
from PySide6.QtGui import QAction, QCursor, QDoubleValidator, QBrush, QColor, QFont
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QTableView, 
    QPushButton, QHeaderView, QAbstractItemView, QMenu,
    QScrollArea, QSizePolicy, QMessageBox, QLabel, QLineEdit, QLayout, QFrame, QScrollBar, QCheckBox, QToolTip
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
        self._highlighted_uuids = set()
        # Store constraints: {uuid: {param: ConstraintObj}}
        self._constraints_map = {} 

    def update_data(self, components: List[ComponentDataV2], constraints_map: dict = None):
        """
        Updates the model with new components and constraint data.
        """
        self.beginResetModel()
        self._components = components
        self._constraints_map = constraints_map if constraints_map else {}
        self.endResetModel()

    def set_highlighted_uuids(self, uuids: set):
        if self._highlighted_uuids == uuids: return
        self._highlighted_uuids = uuids
        if self._components:
            top_left = self.index(0, 0)
            bottom_right = self.index(len(self._components) - 1, len(COLUMNS) - 1)
            self.dataChanged.emit(top_left, bottom_right, [Qt.BackgroundRole])

    def rowCount(self, parent=QModelIndex()): return len(self._components)
    def columnCount(self, parent=QModelIndex()): return len(COLUMNS)
    
    def flags(self, index):
        if not index.isValid(): return Qt.NoItemFlags
        return super().flags(index) | Qt.ItemIsEditable if COLUMNS[index.column()] in EDITABLE_COLS else super().flags(index)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid(): return None
        comp = self._components[index.row()]
        col_name = COLUMNS[index.column()]
        attr = COL_MAP[col_name]
        
        # 1. Background (Group Highlighting)
        if role == Qt.BackgroundRole:
            if comp.uuid in self._highlighted_uuids:
                #return QColor("#FFF8DC") # Cornsilk
                return QColor("deepskyblue")
            return None

        # --- CONSTRAINTS CHECK ---
        is_constrained_col = col_name in ['z', 'logN', 'b', 'btur']
        c_data = None
        is_free = True
        is_linked = False

        if is_constrained_col and comp.uuid in self._constraints_map:
            c_data = self._constraints_map[comp.uuid].get(attr)
            if c_data:
                is_free = c_data.is_free
                is_linked = (c_data.expression is not None) or (c_data.target_uuid is not None)

        # 2. Font Styles
        if role == Qt.FontRole and is_constrained_col:
            font = QFont()
            if is_linked:
                font.setBold(True); return font
            elif not is_free:
                font.setItalic(True); return font
            return None

        # 3. [FIX] Unified Tooltips
        if role == Qt.ToolTipRole:
            # Build Base Info (ID, Series, Z, Stats)
            c_chi2 = f"{comp.chi2:.2f}" if comp.chi2 is not None else "N/A"
            c_resol = f"{comp.resol:.0f}" if comp.resol is not None else "N/A"
            
            tooltip_html = (f"<b>ID:</b> {comp.id} | <b>{comp.series}</b><br>"
                            f"<b>z:</b> {comp.z:.6f}<br>"
                            f"<b>Chi2:</b> {c_chi2} | <b>Resol:</b> {c_resol}")

            # Append Parameter Info if applicable
            if is_constrained_col:
                tooltip_html += "<hr>" # Divider
                if is_linked:
                    target_name = "Unknown"; target_z = ""
                    if c_data.target_uuid:
                        for c_ref in self._components:
                            if c_ref.uuid == c_data.target_uuid:
                                target_name = c_ref.series; target_z = f" (z={c_ref.z:.5f})"; break
                    
                    import re
                    expr_display = c_data.expression if c_data.expression else "Direct Link"
                    def replacer(match):
                        uid = match.group(1)
                        for c_ref in self._components:
                            if c_ref.uuid == uid: return c_ref.series
                        return "Unknown"
                    expr_pretty = re.sub(r"p\['([^']+)'\]", replacer, expr_display)
                    expr_pretty = re.sub(r'p\["([^"]+)"\]', replacer, expr_pretty)
                    
                    tooltip_html += (f"<b>Status:</b> Linked<br>"
                                     f"<b>Target:</b> {target_name}{target_z}<br>"
                                     f"<b>Formula:</b> <code>{expr_pretty}</code>")
                elif not is_free:
                    tooltip_html += "<b>Status:</b> Frozen (Fixed)"
                else:
                    tooltip_html += "<b>Status:</b> Free"
            
            return tooltip_html

        # 4. Text Display
        if role in (Qt.DisplayRole, Qt.EditRole):
            val = getattr(comp, attr, None)
            if role == Qt.DisplayRole and isinstance(val, float):
                return f"{val:.6f}" if attr == 'z' else f"{val:.3f}" if attr in ['logN', 'dlogN'] else f"{val:.2f}"
            return str(val) if val is not None else ""
        if role == Qt.TextAlignmentRole: return Qt.AlignCenter
        return None

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

    def get_component_at(self, row: int): 
        return self._components[row] if 0 <= row < len(self._components) else None
    
# --- Proxy Model ---
class SystemSortFilterProxyModel(QSortFilterProxyModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.group_filtering_enabled = False
        self.allowed_uuids = set()

    def set_group_filter_enabled(self, enabled: bool):
        self.group_filtering_enabled = enabled
        self.invalidateFilter()
        
    def set_allowed_uuids(self, uuids: set):
        """Updates the list of visible UUIDs (Fluid Group)."""
        # Break recursion: If the group hasn't changed, do nothing.
        if self.allowed_uuids == uuids: 
            return
        
        self.allowed_uuids = uuids
        # Only re-filter if the checkbox is actually checked
        if self.group_filtering_enabled:
            self.invalidateFilter()

    def filterAcceptsRow(self, source_row, source_parent):
        # If filtering is off, show everything
        if not self.group_filtering_enabled: return True
        
        # If filtering is on, only show components in the Fluid Group
        model = self.sourceModel()
        comp = model.get_component_at(source_row)
        if not comp: return False
        
        return comp.uuid in self.allowed_uuids

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
        
        # Store the 'Home' view state
        self._home_xlim = (-300, 300)

        self.main_layout = QHBoxLayout(self)
        self.main_layout.setContentsMargins(0, 0, 0, 0)
        self.main_layout.setSpacing(2)

        # Left Side
        left_container = QWidget()
        self.left_layout = QVBoxLayout(left_container)
        self.left_layout.setContentsMargins(0, 0, 0, 0)
        self.left_layout.setSpacing(9)

        self.fig = Figure(figsize=(5, 6), dpi=100)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        self.toolbar = AstrocookToolbar(self.canvas, self)
        
        # --- Toolbar Customization ---
        # 1. Remove 'Select' and 'Subplots' if present
        # 2. Hook up the 'Home' button to our custom logic
        for action in self.toolbar.actions():
            txt = action.text()
            tip = action.toolTip()
            
            if txt == 'Select': 
                self.toolbar.removeAction(action)
            
            # Identify Home button by text or tooltip (standard mpl is 'Reset original view')
            elif txt == 'Home' or 'Reset original view' in tip:
                try:
                    # Disconnect the original matplotlib 'home' method
                    action.triggered.disconnect()
                except RuntimeError:
                    pass # Was not connected or already disconnected
                
                # Connect our custom handler
                action.triggered.connect(self.on_home)
        
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
        self.axes = [] 
        
        self._current_session = None
        self._current_component = None
        self._selected_components = []
        self._plot_center_z = 0.0

        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('draw_event', self.on_draw_update_limits)
        
        # Detect end of pan/zoom to refetch data
        self.canvas.mpl_connect('button_release_event', self.on_mouse_release)

    @property
    def main_window(self):
        return self.inspector.main_window

    def toggle_region_selector(self): pass

    # [FIX 2] Custom Home Logic
    def on_home(self):
        """Restores the view to the initial state for the current component."""
        if not self._current_component: return
        
        # 1. Reset Internal State
        self._xlim = self._home_xlim
        self._plot_center_z = self._current_component.z
        
        # 2. Update UI Boxes
        self.inspector.update_limit_boxes(self._xlim[0], self._xlim[1])
        self.inspector.update_z_box(self._plot_center_z)
        
        # 3. Force Re-plot
        self._update_plot()

    # [FIX 1] Smart Refetch on Pan Release
    def on_mouse_release(self, event):
        """Triggers a data refresh when the mouse is released during Panning."""
        # Check if the toolbar is in Pan or Zoom mode
        if self.toolbar.mode:
            # We delay slightly or just call update to ensure limits are finalized
            self._update_plot()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        row_height = 300 if self.inspector.resid_cb.isChecked() else 200
        self.PAGE_SIZE = max(1, event.size().height() // row_height)
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
        width = self._xlim[1] - self._xlim[0]
        self._xlim = (-width/2, width/2)
        self.inspector.update_limit_boxes(self._xlim[0], self._xlim[1])
        self._update_plot()

    def plot_system(self, session, component: ComponentDataV2, selected_components: List[ComponentDataV2] = None):
        is_new_component = (self._current_component != component)
        
        self._current_session = session
        self._current_component = component
        self._selected_components = selected_components if selected_components else [component]
        
        self._plot_center_z = component.z

        if is_new_component:
            # Center on 0 km/s, preserving current width preference
            width = self._xlim[1] - self._xlim[0]
            self._xlim = (-width/2, width/2)
            
            # [FIX 2] Capture this as the "Home" state
            self._home_xlim = self._xlim
            
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
        self._tick_map_visuals = [] # Store tick positions for tooltips

        if not self._all_transitions or not self._current_session:
            self.canvas.draw(); return

        start = self.scrollbar.value() if self.scrollbar.isEnabled() else 0
        end = min(start + self.PAGE_SIZE, len(self._all_transitions))
        visible_trans = self._all_transitions[start:end]
        
        num_trans = len(visible_trans)
        if num_trans == 0: self.canvas.draw(); return
        
        show_resid = self.inspector.resid_cb.isChecked()
        
        if show_resid:
            gs_ratios = []
            for _ in range(num_trans): gs_ratios.extend([3, 1])
            axs = self.fig.subplots(nrows=num_trans*2, ncols=1, sharex=True, 
                                    gridspec_kw={'height_ratios': gs_ratios, 'hspace': 0.05})
        else:
            axs = self.fig.subplots(nrows=num_trans, ncols=1, sharex=True)
            if num_trans == 1: axs = [axs]

        self.fig.set_layout_engine(layout='constrained')
        
        session = self._current_session
        x_full = session.spec.x.value
        c_kms = 299792.458
        z_sys = self._plot_center_z
        
        has_cont = session.spec.cont is not None
        y = session.spec.y.value
        dy = session.spec.dy.value if session.spec.dy is not None else np.ones_like(y)
        cont = session.spec.cont.value if has_cont else np.ones_like(y)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            y_norm = np.divide(y, cont, where=cont!=0)
            dy_norm = np.divide(dy, cont, where=cont!=0)
            
        y_resid = None
        if session.spec.model is not None:
            mod = session.spec.model.value
            with np.errstate(divide='ignore', invalid='ignore'):
                mod_norm = np.divide(mod, cont, where=cont!=0) if has_cont else mod
                y_resid = np.divide(y_norm - mod_norm, dy_norm, where=dy_norm!=0)

        colors = get_color_cycle(5, cmap='tab20')

        for i, trans_name in enumerate(visible_trans):
            if show_resid:
                ax_main = axs[i*2]
                ax_resid = axs[i*2+1]
                self.axes.extend([ax_main, ax_resid])
            else:
                ax_main = axs[i] if num_trans > 1 else axs[0]
                ax_resid = None
                self.axes.append(ax_main)

            cursor_line_main, = ax_main.plot([], [], color=get_style_color('model', colors), lw=2, alpha=0.3, zorder=10)
            cursor_line_resid = None
            if show_resid:
                cursor_line_resid, = ax_resid.plot([], [], color=get_style_color('model', colors), lw=2, alpha=0.3, zorder=10)
            
            self._panel_map[ax_main] = (trans_name, cursor_line_main, cursor_line_resid)

            def make_fmt(z):
                def fmt(x, y):
                    val_z = (1 + z) * (1 + x / c_kms) - 1
                    return f"Δv={x:.0f}, y={y:.2f}, z={val_z:.5f}"
                return fmt
            ax_main.format_coord = make_fmt(z_sys)

            if trans_name not in xem_d:
                ax_main.text(0.5, 0.5, f"Unknown: {trans_name}", ha='center'); continue
            
            # Panel Reference Physics
            lam_0 = xem_d[trans_name].to(session.spec.x.unit).value
            lam_obs = lam_0 * (1.0 + z_sys)
            v = c_kms * (x_full - lam_obs) / lam_obs
            
            v_min, v_max = self._xlim
            mask = (v > v_min - 500) & (v < v_max + 500)
            v_p = v[mask]; y_p = y_norm[mask]; dy_p = dy_norm[mask]
            
            if len(v_p) == 0: continue

            ax_main.step(v_p, y_p-dy_p, where='mid', color='#aaaaaa', lw=0.3)
            ax_main.step(v_p, y_p+dy_p, where='mid', color='#aaaaaa', lw=0.3)
            ax_main.step(v_p, y_p, where='mid', color=get_style_color('flux', colors), lw=0.8)
            
            if session.spec.model is not None:
                mod_p = mod_norm[mask]
                ax_main.plot(v_p, mod_p, color=get_style_color('model', colors), lw=1.0)

            # Selected Component Overlay
            if self._selected_components:
                atom_info = ATOM_DATA.get(trans_name)
                if atom_info:
                    x_ang_p = session.spec.x.to(au.Angstrom).value[mask]
                    for sel_c in self._selected_components:
                        is_match = (sel_c.series == trans_name) or \
                                   (sel_c.series in STANDARD_MULTIPLETS and trans_name in STANDARD_MULTIPLETS[sel_c.series])
                        if is_match:
                            prof = calc_voigt_profile(
                                x_ang_p, atom_info['wave'], atom_info['f'], atom_info['gamma'],
                                sel_c.z, 10**sel_c.logN, sel_c.b
                            )
                            ax_main.plot(v_p, prof, color='orange', ls='--', lw=1.2, alpha=0.9)

            # --- Ticks for ALL Components (Interlopers included) ---
            tick_ymin, tick_ymax = 0.02, 0.08 
            trans_axis = ax_main.get_xaxis_transform()
            
            # We need the Panel's observed center wavelength (in Angstroms for safety with ATOM_DATA)
            lam_rest_panel_ang = xem_d[trans_name].to(au.Angstrom).value
            lam_obs_panel_center_ang = lam_rest_panel_ang * (1 + z_sys)

            for other_c in session.systs.components:
                # 1. Determine which lines this component has
                comp_lines = []
                if other_c.series in STANDARD_MULTIPLETS:
                    comp_lines = STANDARD_MULTIPLETS[other_c.series]
                elif other_c.series in xem_d:
                    comp_lines = [other_c.series]
                
                # 2. Check each line to see if it falls in this panel
                for line_name in comp_lines:
                    if line_name not in xem_d: continue
                    
                    # Calculate velocity shift relative to THIS panel's center
                    if line_name == trans_name:
                        # Optimization: Same transition -> simple redshift diff
                        v_shift = c_kms * (other_c.z - z_sys) / (1.0 + z_sys)
                    else:
                        # Interloper: Calculate based on observed wavelength coincidence
                        lam_rest_c = xem_d[line_name].to(au.Angstrom).value
                        lam_obs_c = lam_rest_c * (1 + other_c.z)
                        v_shift = c_kms * (lam_obs_c - lam_obs_panel_center_ang) / lam_obs_panel_center_ang

                    # 3. Draw tick if visible
                    if v_min <= v_shift <= v_max:
                        self._tick_map_visuals.append((ax_main, v_shift, other_c))

                        is_h = other_c.series.startswith('Ly') or other_c.series.startswith('H')
                        col = 'red' if is_h else 'gray'
                        
                        # Optional: Use dashed line for interlopers (different transition)?
                        # For now, keep solid to indicate "Real Component Here"
                        ax_main.plot([v_shift, v_shift], [tick_ymin, tick_ymax], 
                                transform=trans_axis, color=col, lw=1.5 if is_h else 1.0, alpha=0.8)
            
            ax_main.axvline(0, color='gray', ls=':', lw=0.8)
            ax_main.axhline(1.0, color='gray', ls=':', lw=0.8)
            ax_main.axhline(0.0, color='gray', ls=':', lw=0.8)
            ax_main.text(0.98, 0.85, trans_name, transform=ax_main.transAxes, ha='right', fontweight='bold', 
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax_main.set_xlim(self._xlim)
            ax_main.set_ylim(-0.1, 1.4) 

            if show_resid and y_resid is not None:
                resid_p = y_resid[mask]
                ax_resid.fill_between(v_p, - 1.0, 1.0, step='mid', color='#aaaaaa', alpha=0.15)
                ax_resid.step(v_p, resid_p, where='mid', color='red', lw=0.8)
                ax_resid.axhline(0, color=get_style_color('model', colors), lw=1.0, alpha=0.5)
                ax_resid.set_ylim(-5, 5)
                ax_resid.set_ylabel(r"χ")
                ax_main.tick_params(labelbottom=False)

        self.fig.supxlabel("Relative velocity (Δv, km/s)")
        if not show_resid:
            self.fig.supylabel("Normalized Flux")
        self.canvas.draw()

    def on_mouse_move(self, event):
        if not event.inaxes or not self._current_component or not self._current_session: return
        
        # [NEW] Check for Tick Tooltips
        closest_dist = 10.0 # pixels tolerance
        tooltip_text = ""

        # Convert mouse data coords to display coords for distance check?
        # Simpler: check data coords width.
        # 10 pixels in data coords depends on zoom. 
        # Let's use a velocity threshold, e.g., 5 km/s or 1% of view.
        view_width = self._xlim[1] - self._xlim[0]
        tol_v = view_width * 0.015 
        
        for ax, v_pos, comp in self._tick_map_visuals:
            if event.inaxes == ax:
                dist = abs(event.xdata - v_pos)
                if dist < tol_v:
                    c_chi2 = f"{comp.chi2:.2f}" if comp.chi2 is not None else "N/A"
                    c_resol = f"{comp.resol:.0f}" if comp.resol is not None else "N/A"
                    tooltip_text = (f"ID: {comp.id} | {comp.series}\n"
                                    f"z: {comp.z:.5f}\n"
                                    f"Chi2: {c_chi2}\n"
                                    f"Res: {c_resol}")
                    break
        
        if tooltip_text:
            # Use global mouse position directly to avoid HiDPI scaling issues
            QToolTip.showText(QCursor.pos(), tooltip_text, self.canvas)
        else:
            QToolTip.hideText()

        main_ax = event.inaxes
        
        if main_ax not in self._panel_map:
            found = False
            for k, val in self._panel_map.items():
                resid_line = val[2]
                if resid_line is not None and resid_line.axes == event.inaxes:
                    main_ax = k
                    found = True
                    break
            if not found: return

        v_mouse = event.xdata
        c_kms = 299792.458
        z_new = (1 + self._plot_center_z) * (1 + v_mouse / c_kms) - 1
        
        N = 10**self.cursor_logN
        b = self.cursor_b
        
        session = self._current_session
        x_unit = session.spec.x.unit

        for ax in self._panel_map:
            trans_name, line_main, line_resid = self._panel_map[ax]
            
            atom_info = ATOM_DATA.get(trans_name)
            if not atom_info: continue

            lam_0_ang = atom_info['wave']
            lam_obs_center_ang = lam_0_ang * (1 + self._plot_center_z)
            
            v_min, v_max = ax.get_xlim()
            v_grid = np.linspace(v_min, v_max, 300)
            
            lam_grid_ang = lam_obs_center_ang * (1 + v_grid / c_kms)
            
            prof = calc_voigt_profile(
                lam_grid_ang, lam_0_ang, atom_info['f'], atom_info['gamma'],
                z_new, N, b
            )
            
            line_main.set_data(v_grid, prof)
            
            if line_resid:
                lam_0_spec = (lam_0_ang * au.Angstrom).to(x_unit).value
                lam_obs_center_spec = lam_0_spec * (1 + self._plot_center_z)

                x_full = session.spec.x.value
                v_full = c_kms * (x_full - lam_obs_center_spec) / lam_obs_center_spec
                
                dy = session.spec.dy.value if session.spec.dy is not None else np.ones_like(x_full)
                cont = session.spec.cont.value if session.spec.cont is not None else np.ones_like(x_full)
                
                with np.errstate(divide='ignore', invalid='ignore'):
                    dy_norm = np.divide(dy, cont, where=cont!=0)
                
                dy_interp = np.interp(v_grid, v_full, dy_norm, left=1.0, right=1.0)
                
                with np.errstate(divide='ignore', invalid='ignore'):
                    prof_sigma = np.divide(prof - 1.0, dy_interp, where=dy_interp>1e-9)
                    prof_sigma[~np.isfinite(prof_sigma)] = 0.0
                
                line_resid.set_data(v_grid, prof_sigma)
            
        self.canvas.draw_idle()

    def on_press(self, event):
        """
        Handles mouse clicks on the velocity plot.
        Right-Click (Button 3):
        1. If clicked on a Tick -> Focus on that component (Group View).
        2. If clicked on Background -> Add new component.
        """
        if event.button == 3 and event.inaxes: 
            menu = QMenu(self)
            has_actions = False
            
            # --- 1. Check for Component Tick Clicks (Focus) ---
            found_comp = None
            # Calculate view width for tolerance (e.g. 1.5%)
            view_width = self._xlim[1] - self._xlim[0]
            tol = view_width * 0.015

            # Iterate through stored ticks
            for ax, v_pos, comp in self._tick_map_visuals:
                if ax == event.inaxes:
                    dist = abs(event.xdata - v_pos)
                    if dist < tol:
                        found_comp = comp
                        break

            if found_comp:
                # "Focus" triggers the Group View logic in the Inspector
                act_focus = QAction(f"Focus on Component {found_comp.id}", menu)
                act_focus.triggered.connect(
                    lambda: self.inspector.focus_on_component(found_comp.uuid, force_group_view=False)
                )
                menu.addAction(act_focus)
                menu.addSeparator()
                has_actions = True

            # --- 2. Background Click (Add Component) ---
            ax_clicked = event.inaxes
            main_ax = None
            
            # Resolve the clicked panel
            if ax_clicked in self._panel_map:
                main_ax = ax_clicked
            else:
                # Reverse lookup (check if clicked on residual panel)
                for ax_m, val in self._panel_map.items():
                    resid_line = val[2]
                    if resid_line is not None and resid_line.axes == ax_clicked:
                        main_ax = ax_m
                        break
            
            if main_ax:
                trans_name, _, _ = self._panel_map[main_ax]
                
                # Calculate coordinates
                v = event.xdata
                c_kms = 299792.458
                z_new = (1 + self._plot_center_z) * (1 + v / c_kms) - 1

                # Determine series name (e.g. CIV_1548 -> CIV)
                series = trans_name
                if series not in STANDARD_MULTIPLETS:
                     for g, l in STANDARD_MULTIPLETS.items():
                        if series in l: series = g; break
                
                # Add Action
                act_add = QAction(f"Add {series} at z={z_new:.5f}", menu)
                act_add.triggered.connect(lambda: self.inspector.main_window._on_recipe_requested(
                    "absorbers", "add_component", {
                        'series': series, 'z': z_new,
                        'logN': self.cursor_logN, 'b': self.cursor_b
                    }, {}))
                
                menu.addAction(act_add)
                has_actions = True

            # --- Execute ---
            if has_actions:
                menu.exec(QCursor.pos())

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
        # [Task 3] Zero margins for the main window to match Main Window style
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(0)
        
        # [Task 1] Reduce splitter size
        splitter = QSplitter(Qt.Horizontal)
        splitter.setHandleWidth(9)
        # Optional: Light gray line to make the thin split visible
        #splitter.setStyleSheet("QSplitter::handle { background-color: #D3D3D3; }")
        
        # --- Left Side: Table ---
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
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(0)
        
        # [Task 2] Create a layout for the checkbox that matches the Right Toolbar's height/margins
        self.group_cb = QCheckBox("Group View (Selected Component Only)")
        self.group_cb.toggled.connect(self._on_group_toggled)
        
        l_toolbar = QHBoxLayout()
        l_toolbar.setContentsMargins(5, 5, 15, 15) # Exact match to Right Toolbar
        l_toolbar.addWidget(self.group_cb)
        l_toolbar.addStretch() # Push checkbox to the left
        
        left_layout.addLayout(l_toolbar)
        left_layout.addWidget(self.table_view)
        
        splitter.addWidget(left_widget)
        
        # --- Right Side: Plot & Controls ---
        right_widget = QWidget() # Renamed locally for clarity
        r_layout = QVBoxLayout(right_widget)
        r_layout.setContentsMargins(0, 0, 0, 0)
        r_layout.setSpacing(0)
        
        # Right Toolbar
        c_layout = QHBoxLayout()
        c_layout.setContentsMargins(5, 5, 5, 5) # Exact match to Left Toolbar
        c_layout.setSpacing(5)

        c_layout.addWidget(QLabel("Trans:"))
        self.trans_in = QLineEdit()
        self.trans_in.setPlaceholderText("CIV, SiIV")
        self.trans_in.returnPressed.connect(self._apply)
        c_layout.addWidget(self.trans_in, 1) 
        
        c_layout.addWidget(QLabel("z:"))
        self.z_in = QLineEdit()
        self.z_in.setFixedWidth(70)
        validator = QDoubleValidator()
        validator.setLocale(QLocale.C) 
        self.z_in.setValidator(validator)
        self.z_in.returnPressed.connect(self._apply_z)
        c_layout.addWidget(self.z_in)

        c_layout.addWidget(QLabel("Δv:"))
        self.vmin_in = QLineEdit("-300")
        self.vmin_in.setFixedWidth(50)
        self.vmin_in.setValidator(validator)
        self.vmin_in.returnPressed.connect(self._apply_limits)
        c_layout.addWidget(self.vmin_in)
        self.vmax_in = QLineEdit("300")
        self.vmax_in.setFixedWidth(50)
        self.vmax_in.setValidator(validator)
        self.vmax_in.returnPressed.connect(self._apply_limits)
        c_layout.addWidget(self.vmax_in)
        
        c_layout.addWidget(QLabel("logN:"))
        self.logn_in = QLineEdit("13.5")
        self.logn_in.setFixedWidth(40)
        self.logn_in.setValidator(validator)
        self.logn_in.textChanged.connect(self._update_cursor)
        c_layout.addWidget(self.logn_in)
        
        c_layout.addWidget(QLabel("b:"))
        self.b_in = QLineEdit("10.0")
        self.b_in.setFixedWidth(40)
        self.b_in.setValidator(validator)
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
        
        splitter.addWidget(right_widget)
        splitter.setStretchFactor(0, 6); splitter.setStretchFactor(1, 4)
        layout.addWidget(splitter)

    def set_session(self, session):
        """
        Updates the table with data from the session and handles row selection logic.
        Prioritizes the 'latest added' component if a new one is detected.
        """
        # 1. Capture previous state (Selection & Max ID)
        sel = self.table_view.selectionModel().selectedRows()
        prev_uuid = None
        if sel:
            idx_source = self.proxy_model.mapToSource(sel[0])
            comp = self.table_model.get_component_at(idx_source.row())
            if comp:
                prev_uuid = comp.uuid
        
        # Detect existing max ID to identify if a new component is added later
        old_comps = self.table_model._components
        max_id_old = max([c.id for c in old_comps], default=0) if old_comps else 0

        # 2. Update Model with new data
        self.current_session = session
        new_comps = session.systs.components if session and session.systs else []
        
        # Extract Constraints Map from the Session
        constraints = {}
        if session and session.systs and session.systs.constraint_model:
            constraints = session.systs.constraint_model.v2_constraints_by_uuid
            
        self.table_model.update_data(new_comps, constraints)
        
        # 3. Determine Selection Target
        max_id_new = max([c.id for c in new_comps], default=0) if new_comps else 0
        target_uuid = prev_uuid
        
        # Rule: If the max ID increased, a component was likely added. Select it.
        if max_id_new > max_id_old:
            for c in new_comps:
                if c.id == max_id_new:
                    target_uuid = c.uuid
                    break
        
        # 4. Apply Selection (Find row by UUID)
        selection_applied = False
        if target_uuid:
            for r in range(self.table_model.rowCount()):
                comp = self.table_model.get_component_at(r)
                if comp and comp.uuid == target_uuid:
                    idx_src = self.table_model.index(r, 0)
                    idx_proxy = self.proxy_model.mapFromSource(idx_src)
                    
                    if idx_proxy.isValid():
                        self.table_view.selectRow(idx_proxy.row())
                        self.table_view.scrollTo(idx_proxy) # Ensure it's visible
                        selection_applied = True
                    break

        # 5. Fallback (Default to top row if no target found or table empty)
        if not selection_applied:
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
        self.proxy_model.set_group_filter_enabled(checked)

    def _on_selection_changed(self, sel, desel):
        indexes = self.table_view.selectionModel().selectedRows()
        if not indexes: 
            # Clear group if selection is cleared
            self._update_group_definition([])
            return
        
        idx_proxy = indexes[0]
        idx_source = self.proxy_model.mapToSource(idx_proxy)
        
        primary_row = idx_source.row()
        primary_comp = self.table_model.get_component_at(primary_row)
        
        selected_comps = []
        for idx in indexes:
            src = self.proxy_model.mapToSource(idx)
            c = self.table_model.get_component_at(src.row())
            if c: selected_comps.append(c)

        # Update the Group Logic (Highlighter + Filter List)
        self._update_group_definition(selected_comps)

        if primary_comp and self.current_session:
            cur = self.trans_in.text()
            toks = [t.strip() for t in cur.split(',') if t.strip()]
            if primary_comp.series and primary_comp.series not in toks: toks.append(primary_comp.series)
            new_txt = ", ".join(toks)
            self.trans_in.setText(new_txt)
            self.z_in.setText(f"{primary_comp.z:.5f}")
            self._parse(new_txt)
            
            self.vel_plot.plot_system(self.current_session, primary_comp, selected_comps)
            
    def _apply(self):
        self._parse(self.trans_in.text())
        self._apply_limits()
        self.vel_plot.plot_spectrum()

    def _parse(self, text):
        self.active_transitions = [t.strip() for t in text.split(',') if t.strip()]

    def _forward(self, name, params):
        if self.main_window: self.main_window._on_recipe_requested("absorbers", name, params, {})

    def focus_on_component(self, uuid: str, force_group_view: bool = False): # <--- Added Arg
        """
        Public API: Selects a component and scrolls to it.
        Args:
            force_group_view: If True, checks 'Group View' at the end (Main Window behavior).
                              If False, restores previous state (Inspector Plot behavior).
        """
        # 1. Capture previous state
        was_checked = self.group_cb.isChecked()
        
        # 2. Temporarily disable grouping to ensure row is found/selectable
        if was_checked:
            self.group_cb.setChecked(False) 
        
        # 3. Find row in Source Model
        target_row = -1
        for r in range(self.table_model.rowCount()):
            c = self.table_model.get_component_at(r)
            if c.uuid == uuid:
                target_row = r
                break
        
        if target_row == -1: 
            if was_checked: self.group_cb.setChecked(True)
            return

        # 4. Map to Proxy Index
        source_idx = self.table_model.index(target_row, 0)
        proxy_idx = self.proxy_model.mapFromSource(source_idx)
        
        if not proxy_idx.isValid(): 
            if was_checked: self.group_cb.setChecked(True)
            return

        # 5. Select Row
        self.table_view.selectRow(proxy_idx.row())
        self.table_view.scrollTo(proxy_idx)
        
        # 6. [FIX] Logic for Group View State
        if force_group_view:
            # Case 1: Main Window -> Always Enforce
            self.group_cb.setChecked(True)
        else:
            # Case 2: Velocity Plot -> Restore User Preference
            if was_checked:
                self.group_cb.setChecked(True)
        
        # 7. Ensure window is active
        self.show()
        self.raise_()
        self.activateWindow()

    def _on_context_menu(self, pos):
        index = self.table_view.indexAt(pos)
        if not index.isValid(): return
        
        # 1. Identify Context (Clicked Cell)
        idx_source = self.proxy_model.mapToSource(index)
        clicked_comp = self.table_model.get_component_at(idx_source.row())
        col_name = COLUMNS[idx_source.column()]
        param_attr = COL_MAP.get(col_name)

        if not clicked_comp: return
        
        m = QMenu(self)
        
        # 2. Check Selection State (Is it a Link Operation?)
        # Get all selected rows (mapped to source model)
        sel_rows = self.table_view.selectionModel().selectedRows()
        source_rows = [self.proxy_model.mapToSource(idx).row() for idx in sel_rows]
        unique_rows = list(set(source_rows))
        
        # LINKING LOGIC: Exactly 2 rows selected, and we clicked on a parameter column
        if len(unique_rows) == 2 and col_name in ['z', 'logN', 'b', 'btur']:
            # Identify Dependent (Clicked) vs Independent (The other one)
            other_row = unique_rows[0] if unique_rows[0] != idx_source.row() else unique_rows[1]
            source_comp = self.table_model.get_component_at(other_row)
            
            if source_comp:
                m.addSection("Parameter Linking")
                
                # A. Link by Value (Identity)
                expr_val = f"p['{source_comp.uuid}'].{param_attr}"
                act_link = QAction(f"Link {col_name} to {source_comp.series} (Value)", m)
                act_link.triggered.connect(
                    lambda: self._toggle_constraint(
                        clicked_comp.uuid, param_attr, False, expr_val, source_comp.uuid # <--- Pass Source UUID
                    )
                )
                m.addAction(act_link)
                
                # B. Link by Temperature (Thermal) - Only for 'b'
                if col_name == 'b':
                    m_dep = self._get_mass(clicked_comp.series)
                    m_src = self._get_mass(source_comp.series)
                    
                    if m_dep and m_src:
                        mass_ratio = np.sqrt(m_src / m_dep)
                        expr_therm = f"p['{source_comp.uuid}'].b * {mass_ratio:.4f}"
                        
                        act_therm = QAction(f"Link {col_name} to {source_comp.series} (Thermal)", m)
                        act_therm.triggered.connect(
                            lambda: self._toggle_constraint(
                                clicked_comp.uuid, param_attr, False, expr_therm, source_comp.uuid # <--- Pass Source UUID
                            )
                        )
                        m.addAction(act_therm)
                
                m.addSeparator()

        # 3. Standard Constraint Actions (Freeze/Unfreeze)
        # (Only show if we are NOT in linking mode, or as secondary options)
        if col_name in ['z', 'logN', 'b', 'btur']:
            is_frozen = False
            is_linked = False
            
            if self.current_session and self.current_session.systs.constraint_model:
                c_map = self.current_session.systs.constraint_model.v2_constraints_by_uuid
                if clicked_comp.uuid in c_map and param_attr in c_map[clicked_comp.uuid]:
                    cons = c_map[clicked_comp.uuid][param_attr]
                    if not cons.is_free: is_frozen = True
                    if cons.expression: is_linked = True
            
            if is_linked:
                # Option to break link
                act = QAction(f"Unlink '{col_name}' (Set Free)", m)
                # Pass None for target_uuid to clear it
                act.triggered.connect(lambda: self._toggle_constraint(clicked_comp.uuid, param_attr, True, None, None))
                m.addAction(act)
            elif is_frozen:
                act = QAction(f"Unfreeze '{col_name}'", m)
                act.triggered.connect(lambda: self._toggle_constraint(clicked_comp.uuid, param_attr, True, None))
                m.addAction(act)
            else:
                act = QAction(f"Freeze '{col_name}'", m)
                act.triggered.connect(lambda: self._toggle_constraint(clicked_comp.uuid, param_attr, False, None))
                m.addAction(act)
            
            m.addSeparator()

        # Set Resolution Action
        act_resol = QAction("Set Resolution...", m)
        if self.main_window:
            act_resol.triggered.connect(lambda: self.main_window._launch_recipe_dialog("edit", "set_properties"))
        m.addAction(act_resol)

        m.addAction("Refit Selected", self._refit)
        m.addAction("Delete", self._delete)
        m.exec(self.table_view.mapToGlobal(pos))

    # Helper for Mass
    def _get_mass(self, series_name):
        # 1. Direct
        if series_name in ATOM_DATA: return ATOM_DATA[series_name]['mass']
        # 2. Multiplet
        if series_name in STANDARD_MULTIPLETS:
            prim = STANDARD_MULTIPLETS[series_name][0]
            if prim in ATOM_DATA: return ATOM_DATA[prim]['mass']
        return None

    def _toggle_constraint(self, uuid, param, set_free, expression=None, target_uuid=None): # <--- Added Arg
        if self.main_window:
            self.main_window._on_recipe_requested(
                "absorbers", "update_constraint", 
                {
                    "uuid": uuid, 
                    "param": param, 
                    "is_free": set_free,
                    "expression": expression,
                    "target_uuid": target_uuid # <--- Pass it
                }, 
                {}
            )
    
    def _delete(self):
        idx = self.table_view.currentIndex()
        src = self.proxy_model.mapToSource(idx)
        comp = self.table_model.get_component_at(src.row())
        if comp and self.main_window:
            confirm = self.main_window._show_custom_message(
                title="Delete Component",
                header="Delete this component?",
                text=f"Series: {comp.series}\nRedshift: {comp.z:.5f}",
                buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                default_btn=QMessageBox.StandardButton.No,
                parent=self
            )
            if confirm == QMessageBox.Yes:
                self.main_window._on_recipe_requested(
                    "absorbers", "delete_component", {"uuid": comp.uuid}, {})

    def _refit(self):
        idx = self.table_view.currentIndex()
        src = self.proxy_model.mapToSource(idx)
        comp = self.table_model.get_component_at(src.row())
        if comp and self.main_window:
            self.main_window._on_recipe_requested("absorbers", "fit_component", {"uuid": comp.uuid}, {})
    
    def _update_group_definition(self, selected_comps: List[ComponentDataV2]):
        """
        Calculates the 'Fluid Group' using the SSOT in SystemListV2.
        """
        if not selected_comps or not self.current_session:
            self.table_model.set_highlighted_uuids(set())
            self.proxy_model.set_allowed_uuids(set())
            return

        seed_uuids = [c.uuid for c in selected_comps]
        
        # Delegate to the API (Single Source of Truth)
        group_uuids = self.current_session.systs.get_connected_group(seed_uuids)

        # Apply to Models
        self.table_model.set_highlighted_uuids(group_uuids)
        self.proxy_model.set_allowed_uuids(group_uuids)