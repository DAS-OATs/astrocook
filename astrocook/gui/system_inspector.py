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
    QPushButton, QHeaderView, QAbstractItemView, QMenu, QInputDialog, QDialog, QDialogButtonBox,
    QScrollArea, QSizePolicy, QMessageBox, QLabel, QLineEdit, QLayout, QFrame, QScrollBar, QCheckBox, QToolTip
)
from typing import List, Optional

from ..core.structures import ComponentDataV2
from ..core.atomic_data import METAL_MULTIPLETS, STANDARD_MULTIPLETS, xem_d, ATOM_DATA
from .pyside_plot import AstrocookToolbar, get_color_cycle, PLOT_STYLE, get_style_color

# --- [CHANGE] Robust Convolution Import with Local Fallback ---
try:
    from ..fitting.voigt_fitter import convolve_flux
except ImportError:
    logging.warning("SystemInspector: Could not import convolve_flux. Using local fallback.")
    convolve_flux = None

# Local Fallback Function (Guaranteed to exist)
def _local_convolve_flux(flux, x, resol, resol_unit='R'):
    """
    Robust fallback convolution using Astropy Gaussian1DKernel.
    """
    try:
        from astropy.convolution import convolve, Gaussian1DKernel
        
        # 1. Determine effective Resolution (Scalar)
        r_eff = resol
        if np.ndim(resol) > 0:
            r_eff = np.nanmedian(resol)
        
        if r_eff is None or r_eff <= 0: 
            return flux

        # 2. Determine Sigma in Angstroms
        # FWHM = lambda / R   (if R)
        # FWHM = lambda * (km/s / c) (if km/s)
        lam_center = np.nanmedian(x)
        
        if resol_unit == 'R':
            fwhm_ang = lam_center / r_eff
        else: # 'km/s'
            fwhm_ang = (r_eff / 299792.458) * lam_center
        
        # Sigma = FWHM / 2.355
        sigma_ang = fwhm_ang / 2.35482 
        
        # 3. Determine Sigma in Pixels
        dx = np.nanmedian(np.diff(x))
        if dx == 0: return flux
        
        sigma_pix = sigma_ang / dx
        
        # 4. Convolve (if kernel is at least somewhat significant)
        if sigma_pix < 0.01: return flux # Relaxed threshold

        kernel = Gaussian1DKernel(stddev=sigma_pix)
        return convolve(flux, kernel, boundary='extend')
        
    except ImportError:
        logging.error("SystemInspector: Astropy not found for convolution.")
        return flux
    except Exception as e:
        logging.error(f"SystemInspector: Convolution failed: {e}")
        return flux

# --- 1. Table Configuration (UPDATED) ---
# We removed the separate error columns (dz, dlogN, db)
COL_MAP = {
    "ID": "id", 
    "Trans.": "series", 
    "z": "z", 
    "logN": "logN", 
    "b": "b", 
    "btur": "btur"
}
COLUMNS = list(COL_MAP.keys())
EDITABLE_COLS = {"z", "logN", "b", "btur"}

# Mapping to find the error attribute for a given value attribute
ERR_MAP = {
    "z": "dz",
    "logN": "dlogN",
    "b": "db",
    "btur": "dbtur"
}

class SystemTableModel(QAbstractTableModel):
    data_changed_request = Signal(str, dict)

    def __init__(self, components: List[ComponentDataV2] = None):
        super().__init__()
        self._components = components if components else []
        self._primary_uuids = set()
        self._secondary_uuids = set()

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

    def set_highlighted_uuids(self, primary_uuids: set, secondary_uuids: set):
        if self._primary_uuids == primary_uuids and self._secondary_uuids == secondary_uuids:
            return
            
        self._primary_uuids = primary_uuids
        self._secondary_uuids = secondary_uuids
        
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
            if comp.uuid in self._secondary_uuids:
                return QColor(245, 161, 0, 50) 
            return None

        # --- CONSTRAINTS CHECK ---
        # (We calculate this upfront so we can use it for Font, Tooltip, AND Display)
        is_constrained_col = col_name in ['z', 'logN', 'b', 'btur']
        c_data = None
        is_free = True
        is_linked = False

        if is_constrained_col and comp.uuid in self._constraints_map:
            c_data = self._constraints_map[comp.uuid].get(attr)
            if c_data:
                is_free = c_data.is_free
                is_linked = (c_data.expression is not None) or (c_data.target_uuid is not None)

        if role == Qt.FontRole and is_constrained_col:
            font = QFont()
            if is_linked:
                font.setBold(True); return font
            elif not is_free:
                font.setItalic(True); return font
            return None

        # 2. Tooltip (Keep existing logic)
        if role == Qt.ToolTipRole:
             c_chi2 = f"{comp.chi2:.2f}" if comp.chi2 is not None else "N/A"
             c_resol = f"{comp.resol:.0f}" if comp.resol is not None else "N/A"
             tooltip_html = (
                f"<div style='font-size: 13px'>"
                f"<b>ID:</b> {comp.id} | <b>{comp.series}</b><br>"
                f"<b>z:</b> {comp.z:.6f}<br>"
                f"<b>\u03c7\u00b2:</b> {c_chi2} | <b>R:</b> {c_resol}"
             )
             if is_constrained_col:
                tooltip_html += "<hr style='margin: 4px 0;'>"
                if is_linked:
                    target_name = "Unknown"; target_z = ""
                    if c_data.target_uuid:
                        for c_ref in self._components:
                            if c_ref.uuid == c_data.target_uuid:
                                target_name = c_ref.series; target_z = f" (z={c_ref.z:.5f})"; break
                    import re
                    expr_display = c_data.expression if c_data.expression else "Direct Link"
                    # Clean up expression for display
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
             tooltip_html += "</div>"
             return tooltip_html

        # 3. DISPLAY & EDIT LOGIC (Refined)
        if role == Qt.DisplayRole or role == Qt.EditRole:
            val = getattr(comp, attr, None)
            
            val_fmt = "{}"
            err_fmt = "{}"
            
            if attr == 'z':
                val_fmt = "{:.6f}"; err_fmt = "{:.2e}"
            elif attr == 'logN':
                val_fmt = "{:.3f}"; err_fmt = "{:.3f}"
            elif attr == 'b':
                val_fmt = "{:.2f}"; err_fmt = "{:.2f}"
            elif attr == 'btur':
                val_fmt = "{:.2f}"; err_fmt = "{:.2f}"

            # --- EDIT ROLE: Clean Value Only ---
            if role == Qt.EditRole:
                if isinstance(val, (int, float)):
                    return val_fmt.format(val)
                return str(val) if val is not None else ""
            
            # --- DISPLAY ROLE: Value ± Error ---
            if role == Qt.DisplayRole:
                if val is None: return ""
                
                if isinstance(val, (int, float)):
                    val_str = val_fmt.format(val)
                else:
                    val_str = str(val)

                # Append Error if it exists... AND if not linked!
                if attr in ERR_MAP and isinstance(val, (int, float)):
                    
                    # --- [FIX] Suppress error if linked ---
                    if not is_linked: 
                        err_attr = ERR_MAP[attr]
                        err_val = getattr(comp, err_attr, None)
                        
                        if err_val is not None and np.isfinite(err_val) and err_val > 0:
                            err_str = err_fmt.format(err_val)
                            return f"{val_str} ± {err_str}"
                    # --------------------------------------
                
                return val_str

        # 4. Sorting (Raw Float)
        if role == Qt.UserRole:
            val = getattr(comp, attr, 0)
            if val is None: return -np.inf
            try: return float(val)
            except: return str(val)

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
        if self.allowed_uuids == uuids: return
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
        l_data = self.sourceModel().data(left, Qt.UserRole)
        r_data = self.sourceModel().data(right, Qt.UserRole)
        
        # Handle Nones/Types safely
        if l_data is None: l_data = -np.inf
        if r_data is None: r_data = -np.inf
        
        try: 
            return float(l_data) < float(r_data)
        except (ValueError, TypeError): 
            return str(l_data) < str(r_data)
        
def calc_voigt_profile(wave_grid_ang, lambda_0, f_val, gamma, z, N, b_kms, resol=None, resol_unit='R', context="Unknown"):
    if wave_grid_ang is None or len(wave_grid_ang) == 0: return np.array([])
    
    # [FIX] Clean inputs to avoid Unit errors
    if hasattr(lambda_0, 'value'): lambda_0 = lambda_0.value
    if hasattr(wave_grid_ang, 'value'): wave_grid_ang = wave_grid_ang.value
    if hasattr(b_kms, 'value'): b_kms = b_kms.value

    # 1. Physics
    lambda_c = lambda_0 * (1.0 + z)
    c_kms = 2.99792458e5
    b_safe = max(b_kms, 0.1)
    dop_width_ang = (b_safe / c_kms) * lambda_c
    if dop_width_ang <= 0: dop_width_ang = 1e-5
        
    x = (wave_grid_ang - lambda_c) / dop_width_ang
    c_ang_s = 2.99792458e18
    a = (gamma * lambda_c**2) / (4.0 * np.pi * c_ang_s * dop_width_ang)
    
    H_ax = wofz(x + 1j * a).real
    tau = 1.4974e-15 * N * f_val * lambda_0 / b_safe * H_ax
    flux = np.exp(-tau)

    # 2. Convolution (SSOT)
    should_convolve = False
    if resol is not None:
        if np.ndim(resol) == 0: # Scalar
            if resol > 0: should_convolve = True
        else: # Array
            if np.any(np.nan_to_num(resol) > 0):
                should_convolve = True
    
    if should_convolve:
        try:
            # Try importing main function if defined
            if convolve_flux is not None:
                try:
                    return convolve_flux(flux, wave_grid_ang, resol, resol_unit=resol_unit)
                except TypeError:
                    return convolve_flux(flux, wave_grid_ang, resol)
            else:
                raise ImportError("Main convolve_flux not available")
                
        except Exception as e:
            # Force logging here to see if the main convolution is crashing
            logging.warning(f"DEBUG: Main convolution failed ({e}). Falling back to local.")
            return _local_convolve_flux(flux, wave_grid_ang, resol, resol_unit)
    else:
        return flux

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
        self._tick_map_visuals = []
        
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

    # Custom Home Logic
    def on_home(self):
        """Restores the view to the initial state for the current component."""
        if not self._current_component: return
        self._xlim = self._home_xlim
        self._plot_center_z = self._current_component.z
        self.inspector.update_limit_boxes(self._xlim[0], self._xlim[1])
        self.inspector.update_z_box(self._plot_center_z)
        self._update_plot()

    # Smart Refetch on Pan Release
    def on_mouse_release(self, event):
        """Triggers a data refresh when the mouse is released during Panning."""
        if self.toolbar.mode:
            # We must sync self._xlim with the axes' current state (set by the Zoom tool).
            # If we don't, _update_plot() will wipe the figure and rebuild it using 
            # the OLD limits, cancelling the zoom.
            if self.axes and len(self.axes) > 0:
                # Get the limits from the first axis (they are shared)
                self._xlim = self.axes[0].get_xlim()
                
                # Update the text boxes in the Inspector UI
                self.inspector.update_limit_boxes(self._xlim[0], self._xlim[1])

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
            width = self._xlim[1] - self._xlim[0]
            self._xlim = (-width/2, width/2)
            self._home_xlim = self._xlim
            self.inspector.update_limit_boxes(self._xlim[0], self._xlim[1])

        if not session or not component:
            self._all_transitions = []
            self._update_plot()
            return

        expanded = []
        for item in self.inspector.active_transitions:
            if item in xem_d: 
                expanded.append(item)
            elif item in STANDARD_MULTIPLETS: 
                expanded.extend(STANDARD_MULTIPLETS[item])
        
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
            c_kms = 299792.458
            v_center = np.mean(xlim)
            z_view_center = (1 + self._plot_center_z) * (1 + v_center / c_kms) - 1
            self.inspector.update_z_box(z_view_center)

    def _get_resolution_at(self, target_wavelength_nm: float) -> tuple[float, str]:
        """
        Robustly determines resolution for a specific wavelength.
        Returns: (Value, Unit) -> (float, 'R' or 'km/s')
        """
        sess = self._current_session
        if not sess or not sess.spec: return 0.0, 'R'

        # 1. Check for 'resol' Column (Variable Resolution)
        # This contains FWHM in km/s
        if sess.spec.has_aux_column('resol'):
            try:
                resol_col = sess.spec.get_column('resol').value
                x_arr = sess.spec.x.value # Assumed sorted
                
                # Find nearest pixel
                idx = np.searchsorted(x_arr, target_wavelength_nm)
                idx = np.clip(idx, 0, len(x_arr)-1)
                
                val = resol_col[idx]
                if np.isfinite(val) and val > 0:
                    return val, 'km/s'
            except Exception:
                pass

        # 2. Fallback to Global Metadata (Usually 'R')
        if hasattr(sess.spec, '_data'):
            r = getattr(sess.spec._data, 'resol', 0.0)
            if r > 0: return r, 'R'
            
        meta = sess.spec.meta
        if 'resol' in meta: return float(meta['resol']), 'R'
        if 'RESOL' in meta: return float(meta['RESOL']), 'R'
        
        return 0.0, 'R'

    def _update_plot(self):
        self.fig.clear()
        self.axes = []
        self._panel_map = {} 
        self._tick_map_visuals = [] 

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

        resol_col = None
        if session.spec.has_aux_column('resol'):
            resol_col = session.spec.get_column('resol').value

        has_norm = session.spec.norm is not None
        y = session.spec.y.value
        dy = session.spec.dy.value if session.spec.dy is not None else np.ones_like(y)
        norm = session.spec.norm.value if has_norm else np.ones_like(y)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            y_norm = np.divide(y, norm, where=norm!=0)
            dy_norm = np.divide(dy, norm, where=norm!=0)
            
        y_resid = None
        if session.spec.model is not None:
            mod = session.spec.model.value
            with np.errstate(divide='ignore', invalid='ignore'):
                mod_norm = np.divide(mod, norm, where=norm!=0) if has_norm else mod
                y_resid = np.divide(y_norm - mod_norm, dy_norm, where=dy_norm!=0)

        colors = get_color_cycle(5, cmap='tab20')
        #colors[0] = "#296bff"

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
            
            lam_0 = xem_d[trans_name].to(session.spec.x.unit).value
            lam_obs = lam_0 * (1.0 + z_sys)
            v = c_kms * (x_full - lam_obs) / lam_obs
            
            v_min, v_max = self._xlim
            mask = (v > v_min - 500) & (v < v_max + 500)
            v_p = v[mask]; y_p = y_norm[mask]; dy_p = dy_norm[mask]
            
            if len(v_p) == 0: continue
            
            # --- 1. Determine "Panel Resolution" (Column or Global) ---
            # This is the baseline resolution for the data in this panel.
            found_valid_col_data = False
            panel_res_arg = None
            panel_res_unit = 'R'
            
            if resol_col is not None:
                temp_res = resol_col[mask]
                if np.any(np.isfinite(temp_res)):
                    panel_res_arg = temp_res
                    found_valid_col_data = True
                    if np.nanmedian(panel_res_arg) > 500.0: panel_res_unit = 'R'
                    else: panel_res_unit = 'km/s'

            if not found_valid_col_data:
                # Fallback to Global (Metadata)
                panel_res_arg, panel_res_unit = self._get_resolution_at(lam_obs)

            ax_main.step(v_p, y_p-dy_p, where='mid', color='#aaaaaa', lw=0.3)
            ax_main.step(v_p, y_p+dy_p, where='mid', color='#aaaaaa', lw=0.3)
            ax_main.step(v_p, y_p, where='mid', color=get_style_color('flux', colors), lw=0.8)
            
            if session.spec.model is not None:
                mod_p = mod_norm[mask]
                ax_main.plot(v_p, mod_p, color=get_style_color('model', colors), lw=1.0)

            # --- Plot Components ---
            if self._selected_components:
                atom_info = ATOM_DATA.get(trans_name)
                if atom_info:
                    x_ang_p = session.spec.x.to(au.Angstrom).value[mask]
                    for sel_c in self._selected_components:
                        is_match = (sel_c.series == trans_name) or \
                                   (sel_c.series in STANDARD_MULTIPLETS and trans_name in STANDARD_MULTIPLETS[sel_c.series])
                        if is_match:
                            
                            # --- [FIX] Resolution Priority ---
                            # 1. Spectrum Column (if present)
                            # 2. Component Attribute (if present)
                            # 3. Global Attribute (Default/Panel)
                            
                            comp_res_arg = panel_res_arg
                            comp_res_unit = panel_res_unit
                            
                            # Only check for override if we are NOT using a spectrum column
                            if not found_valid_col_data:
                                # Check if component has specific resolution
                                c_res = getattr(sel_c, 'resol', None)
                                if c_res is not None and c_res > 0:
                                    comp_res_arg = c_res
                                    # Heuristic for unit
                                    comp_res_unit = 'R' if c_res > 500 else 'km/s'

                            b_eff = np.sqrt(sel_c.b**2 + sel_c.btur**2)

                            prof = calc_voigt_profile(
                                x_ang_p, atom_info['wave'], atom_info['f'], atom_info['gamma'],
                                sel_c.z, 10**sel_c.logN, b_eff, 
                                resol=comp_res_arg, resol_unit=comp_res_unit,
                                context="STATIC_PLOT" 
                            )
                            ax_main.plot(v_p, prof, color='#f5a100', ls='--', lw=1.2, alpha=0.9)

            # ... (Rest of plotting logic: Ticks, Labels, Residuals) ...
            tick_ymin, tick_ymax = 0.02, 0.08 
            trans_axis = ax_main.get_xaxis_transform()
            lam_rest_panel_ang = xem_d[trans_name].to(au.Angstrom).value
            lam_obs_panel_center_ang = lam_rest_panel_ang * (1 + z_sys)

            for other_c in session.systs.components:
                comp_lines = []
                if other_c.series in STANDARD_MULTIPLETS:
                    comp_lines = STANDARD_MULTIPLETS[other_c.series]
                elif other_c.series in xem_d:
                    comp_lines = [other_c.series]
                
                for line_name in comp_lines:
                    if line_name not in xem_d: continue
                    if line_name == trans_name:
                        v_shift = c_kms * (other_c.z - z_sys) / (1.0 + z_sys)
                    else:
                        lam_rest_c = xem_d[line_name].to(au.Angstrom).value
                        lam_obs_c = lam_rest_c * (1 + other_c.z)
                        v_shift = c_kms * (lam_obs_c - lam_obs_panel_center_ang) / lam_obs_panel_center_ang

                    if v_min <= v_shift <= v_max:
                        self._tick_map_visuals.append((ax_main, v_shift, other_c))
                        is_h = other_c.series.startswith('Ly') or other_c.series.startswith('H')
                        col = 'red' if is_h else 'gray'
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
        
        closest_dist = 10.0 
        tooltip_text = ""
        view_width = self._xlim[1] - self._xlim[0]
        tol_v = view_width * 0.015 
        
        for ax, v_pos, comp in self._tick_map_visuals:
            if event.inaxes == ax and abs(event.xdata - v_pos) < tol_v:
                c_chi2 = f"{comp.chi2:.2f}" if comp.chi2 is not None else "N/A"
                c_resol = f"{comp.resol:.0f}" if comp.resol is not None else "N/A"
                tooltip_text = (
                    f"<div style='font-size: 13px'>"
                    f"<b>ID:</b> {comp.id} | <b>{comp.series}</b><br>"
                    f"<b>z:</b> {comp.z:.5f}<br>"
                    f"<b>\u03c7\u00b2:</b> {c_chi2} | <b>R:</b> {c_resol}"
                    f"</div>"
                )
                break
        
        if tooltip_text: QToolTip.showText(QCursor.pos(), tooltip_text, self.canvas)
        else: QToolTip.hideText()

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
        x_full = session.spec.x.value
        resol_col = None
        if session.spec.has_aux_column('resol'):
            resol_col = session.spec.get_column('resol').value

        for ax in self._panel_map:
            trans_name, line_main, line_resid = self._panel_map[ax]
            
            siblings = [trans_name]
            for series_key, members in STANDARD_MULTIPLETS.items():
                if trans_name in members:
                    siblings = [m for m in members if m in self._all_transitions]
                    break
            
            atom_info = ATOM_DATA.get(trans_name)
            if not atom_info: continue
            
            lam_0_ang = atom_info['wave']
            lam_obs_center_ang = lam_0_ang * (1 + self._plot_center_z)
            v_min, v_max = ax.get_xlim()
            v_grid = np.linspace(v_min, v_max, 300)
            lam_grid_ang = lam_obs_center_ang * (1 + v_grid / c_kms)
            
            # --- [FIX] ROBUST RESOLUTION LOGIC FOR CURSOR ---
            res_arg = None
            res_unit = 'R'
            found_valid_res = False
            
            if resol_col is not None:
                # Interpolate resolution onto the cursor grid
                lam_grid_nm = lam_grid_ang / 10.0
                # Use np.interp but allow it to carry NaNs/junk if present
                res_temp = np.interp(lam_grid_nm, x_full, resol_col, left=np.nan, right=np.nan)
                
                # Check validity
                if np.any(np.isfinite(res_temp)):
                    res_arg = res_temp
                    found_valid_res = True
                    
                    # Smart Unit Detection
                    median_res = np.nanmedian(res_arg)
                    if median_res > 500.0:
                        res_unit = 'R'
                    else:
                        res_unit = 'km/s'

            if not found_valid_res:
                 # Fallback to scalar global resolution
                 lam_cursor_nm = (lam_0_ang * (1+z_new)) / 10.0
                 res_arg, res_unit = self._get_resolution_at(lam_cursor_nm)
            # ------------------------------------------------

            total_prof = np.ones_like(v_grid)
            for sibling in siblings:
                s_info = ATOM_DATA.get(sibling)
                if not s_info: continue
                
                prof = calc_voigt_profile(
                    lam_grid_ang, s_info['wave'], s_info['f'], s_info['gamma'],
                    z_new, N, b,
                    resol=res_arg, resol_unit=res_unit,
                    context="CURSOR" 
                )
                total_prof *= prof
            
            line_main.set_data(v_grid, total_prof)
            
            if line_resid:
                lam_0_spec = atom_info['wave'] / 10.0
                lam_obs_center_spec = lam_0_spec * (1 + self._plot_center_z)
                v_full = c_kms * (x_full - lam_obs_center_spec) / lam_obs_center_spec
                
                dy = session.spec.dy.value if session.spec.dy is not None else np.ones_like(x_full)
                norm = session.spec.norm.value if session.spec.norm is not None else np.ones_like(x_full)
                with np.errstate(divide='ignore', invalid='ignore'):
                    dy_norm = np.divide(dy, norm, where=norm!=0)
                
                dy_interp = np.interp(v_grid, v_full, dy_norm, left=1.0, right=1.0)
                with np.errstate(divide='ignore', invalid='ignore'):
                    prof_sigma = np.divide(total_prof - 1.0, dy_interp, where=dy_interp>1e-9)
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
        if not event.inaxes: return

        # --- 0. CHECK MODIFIERS ---
        # Use Qt directly to reliably detect 'Ctrl' key even if plot focus is fuzzy
        modifiers = QApplication.keyboardModifiers()
        is_ctrl_held = (modifiers & Qt.ControlModifier)

        # Check if a toolbar mode (Zoom/Pan) is currently active
        is_tool_active = bool(self.toolbar.mode)

        # Condition: Right Click (Button 3) AND (Ctrl held OR No tool active)
        if event.button == 3 and (is_ctrl_held or not is_tool_active):
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
            if ax_clicked in self._panel_map: main_ax = ax_clicked
            else:
                # Reverse lookup (check if clicked on residual panel)
                for ax_m, val in self._panel_map.items():
                    if val[2] is not None and val[2].axes == ax_clicked:
                        main_ax = ax_m
                        break
            
            if main_ax:
                trans_name, _, _ = self._panel_map[main_ax]
                
                # Calculate coordinates
                v = event.xdata
                c_kms = 299792.458
                z_new = (1 + self._plot_center_z) * (1 + v / c_kms) - 1

                # --- PREPARATION ---
                
                # 1. Get Active Transitions (Content of Trans. Box)
                active_trans = self._all_transitions

                # 2. Analyze Species Composition
                # Group active lines by their Species (e.g. {'FeII': [...], 'CIV': [...]})
                species_map = {}
                
                for t in active_trans:
                    # Identify Species Tag
                    g_found = None
                    for g, members in STANDARD_MULTIPLETS.items():
                        if t in members: 
                            g_found = g
                            break
                    # Fallback: extract prefix (e.g. 'FeII' from 'FeII_2374')
                    if not g_found:
                        g_found = t.split('_')[0]
                    
                    if g_found not in species_map:
                        species_map[g_found] = []
                    species_map[g_found].append(t)

                unique_species = list(species_map.keys())
                num_species = len(unique_species)

                # --- MENU CONSTRUCTION ---

                # 1. SINGLE LINE ADD (Always Available)
                # Adds the specific transition under cursor
                act_add_single = QAction(f"Add {trans_name} (Single Line) at z={z_new:.5f}", menu)
                act_add_single.triggered.connect(lambda: self.inspector.main_window._on_recipe_requested(
                    "absorbers", "add_component", {
                        'series': trans_name, 'z': z_new,
                        'logN': self.cursor_logN, 'b': self.cursor_b
                    }, {}))
                menu.addAction(act_add_single)

                # LOGIC BRANCHING BASED ON SCENARIOS

                # CASE A: SINGLE SPECIES (e.g. Just CIV or Just FeII)
                if num_species == 1:
                    sp_name = unique_species[0]
                    std_members = STANDARD_MULTIPLETS.get(sp_name, [])
                    
                    # Check definition size: Small (<=2, e.g. CIV) vs Large (>2, e.g. FeII)
                    is_large_multiplet = len(std_members) > 2
                    
                    if not is_large_multiplet:
                        # SCENARIO 1: Small Multiplet -> STANDARD ADD
                        # (Adds the full standard group, e.g. CIV)
                        act_std = QAction(f"Add {sp_name} (Standard) at z={z_new:.5f}", menu)
                        act_std.triggered.connect(lambda: self.inspector.main_window._on_recipe_requested(
                            "absorbers", "add_component", {
                                'series': sp_name, 'z': z_new,
                                'logN': self.cursor_logN, 'b': self.cursor_b
                            }, {}))
                        menu.addAction(act_std)
                    
                    else:
                        # SCENARIO 2: Large Multiplet -> ADD SUBSET
                        # (Adds exactly the active lines as a custom multiplet)
                        # Only show if subset > 1 (otherwise Single Line covers it)
                        if len(active_trans) > 1:
                            v1_series_key = ",".join(active_trans)
                            label_str = v1_series_key
                            if len(label_str) > 35: label_str = label_str[:32] + "..."

                            act_sub = QAction(f"Add {label_str} (Subset) at z={z_new:.5f}", menu)
                            act_sub.triggered.connect(lambda checked=False, s=v1_series_key, z=z_new: 
                                self.inspector.main_window._on_recipe_requested(
                                    "absorbers", "add_component", 
                                    {
                                        'series': s, 
                                        'z': z,
                                        'logN': self.cursor_logN, 
                                        'b': self.cursor_b
                                    }, {}
                                )
                            )
                            menu.addAction(act_sub)

                # CASE B: MIXED SPECIES (e.g. CIV + FeII)
                elif num_species > 1:
                    # SCENARIO 3: Mixed -> LINKED SYSTEM (Active Lines Only)
                    # Pass the raw list of active transitions to add_linked_system.
                    # This ensures we link specific visible lines, not whole groups.
                    
                    raw_list_str = ",".join(active_trans)
                    
                    # Generate label (e.g. "CIV_1548,FeII_2374...")
                    label_str = raw_list_str
                    if len(label_str) > 35: label_str = label_str[:32] + "..."
                    
                    act_link = QAction(f"Add {label_str} (Linked) at z={z_new:.5f}", menu)
                    act_link.triggered.connect(lambda checked=False, sl=raw_list_str, z=z_new: 
                        self.inspector.main_window._on_recipe_requested(
                            "absorbers", "add_linked_system",
                            {
                                'series_list': sl, 'z': z,
                                'logN': self.cursor_logN, 'b': self.cursor_b
                            }, {}
                        )
                    )
                    menu.addAction(act_link)

                has_actions = True

            # --- Execute ---
            if has_actions: menu.exec(QCursor.pos())

class SystemInspector(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.current_session = None
        self.active_transitions: List[str] = [] 
        self._manual_trans_edit = False 
        self._manual_base_text = ""  # Stores the persistent manual text
        self.setWindowTitle("System Inspector")
        self.resize(1200, 800) 

        # [CHANGE] Apply Rounded Tooltip Style
        self.setStyleSheet("""
            QTableView {
                selection-background-color: #f5a100;
                selection-color: black;
            }
            QToolTip {
                border: 1px solid gray;
                border-radius: 5px;
                background-color: #f2f8f8;
                color: black;
                padding: 4px;
                background-clip: border-box;
            }
        """)

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

        c_layout.addWidget(QLabel("Trans.:"))
        self.trans_in = QLineEdit()
        self.trans_in.setPlaceholderText("CIV, SiIV")
        self.trans_in.returnPressed.connect(self._apply)
        self.trans_in.textEdited.connect(self._on_manual_edit)
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
        
        self.resid_cb = QCheckBox("Resid.")
        self.resid_cb.toggled.connect(self._apply) 
        c_layout.addWidget(self.resid_cb)
        
        pal = QApplication.palette()
        style = f"""
            QLineEdit {{ 
                padding: 3px; 
                border-radius: 4px; 
                background: {pal.color(pal.ColorRole.Base).name()}; 
                color: {pal.color(pal.ColorRole.Text).name()}; }}
            QLineEdit:focus {{
                border: 1px solid #296bff; /* Highlight color */
            }}
        """
        for w in [self.trans_in, self.z_in, self.vmin_in, self.vmax_in, self.logn_in, self.b_in]: w.setStyleSheet(style)
        
        r_layout.addLayout(c_layout)
        
        self.vel_plot = VelocityPlotWidget(self)
        r_layout.addWidget(self.vel_plot, 1) 
        
        splitter.addWidget(right_widget)
        splitter.setStretchFactor(0, 6); splitter.setStretchFactor(1, 4)
        layout.addWidget(splitter)

        # Add helper to set style based on state
        self._update_trans_box_visuals()
    
    def _on_manual_edit(self, text):
        """
        Called ONLY when the user types in the box. 
        Updates the 'Base' text which serves as the anchor for fluid prepending.
        """
        if not text.strip():
            self._manual_trans_edit = False
            self._manual_base_text = ""
        else:
            self._manual_trans_edit = True
            self._manual_base_text = text

        # Add helper to set style based on state
        self._update_trans_box_visuals()

    def set_session(self, session):
        """
        Updates the table with data from the session and handles row selection logic.
        Prioritizes the 'latest added' component if a new one is detected.
        """
        if session:
            self.setWindowTitle(f"System Inspector: {session.name}")
        else:
            self.setWindowTitle("System Inspector")

        # 1. Capture previous state (Selection & Max ID)
        sel = self.table_view.selectionModel().selectedRows()
        prev_uuid = None
        if sel:
            idx_source = self.proxy_model.mapToSource(sel[0])
            comp = self.table_model.get_component_at(idx_source.row())
            if comp: prev_uuid = comp.uuid
        
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
            if self.proxy_model.rowCount() > 0: self.table_view.selectRow(0)
            else: self.vel_plot.plot_spectrum(None)
                
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
        except ValueError as e: logging.error(f"Error applying z: {e}")

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
            new_txt = ""
            
            # Use dict.fromkeys to preserve order while removing duplicates
            unique_series = list(dict.fromkeys([c.series for c in selected_comps if c.series]))
            series_str = ", ".join(unique_series)

            if not self._manual_trans_edit:
                # --- MODE A: Auto-Replace ---
                # The box reflects the FULL selection.
                new_txt = series_str
                # Update the BASE so if the user starts typing now, 
                # they append to the CURRENT selection.
                self._manual_base_text = new_txt 
                
            else:
                # --- MODE B: Fluid Prepend ---
                # We prepend the selection to the BASE text (what user typed),
                # NOT to the current box text (which might have old temp selections).
                
                # Check for duplicates in the base text to be clean
                base_tokens = [t.strip() for t in self._manual_base_text.split(',') if t.strip()]
                
                if primary_comp.series in base_tokens:
                    # Already in the manual list: just show the manual list
                    new_txt = self._manual_base_text
                else:
                    # Prepend: "Selection, Base"
                    new_txt = f"{primary_comp.series}, {self._manual_base_text}"
            
            # Update the widget (Programmatic setText does NOT fire textEdited)
            self.trans_in.setText(new_txt)
            self.z_in.setText(f"{primary_comp.z:.5f}")
            
            # Add helper to set style based on state
            self._update_trans_box_visuals()

            # Update Internal State & Plot
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
        
        force_group_view: If True, checks 'Group View' at the end (Main Window behavior). If False, restores previous state (Inspector Plot behavior).
        """
        # 1. Capture previous state
        was_checked = self.group_cb.isChecked()
        
        # 2. Temporarily disable grouping to ensure row is found/selectable
        if was_checked: self.group_cb.setChecked(False) 
        
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
            if was_checked: self.group_cb.setChecked(True)
        
        # 7. Ensure window is active
        self.show()
        self.raise_()
        self.activateWindow()

    def _get_selected_components(self) -> List[ComponentDataV2]:
        rows = self.table_view.selectionModel().selectedRows()
        comps = []
        for idx in rows:
            src = self.proxy_model.mapToSource(idx)
            c = self.table_model.get_component_at(src.row())
            if c: comps.append(c)
        return comps

    def _on_context_menu(self, pos):
        index = self.table_view.indexAt(pos)
        if not index.isValid(): return
        
        # 1. Identify Context (Clicked Cell)
        idx_source = self.proxy_model.mapToSource(index)
        clicked_comp = self.table_model.get_component_at(idx_source.row())
        col_name = COLUMNS[idx_source.column()]
        param_attr = COL_MAP.get(col_name)

        if not clicked_comp: return
        
        selected_comps = self._get_selected_components()
        
        # Guard: If right-clicking outside the selection, reset selection to clicked row
        # (Standard GUI behavior: Right-click selects the target if not already selected)
        if clicked_comp not in selected_comps:
            selected_comps = [clicked_comp]
            
        count = len(selected_comps)
        is_multi = count > 1
        suffix = f" ({count} items)" if is_multi else ""
        
        m = QMenu(self)
        
        # 1. Linking Logic (Strictly 2 items)
        if count == 2 and col_name in ['z', 'logN', 'b', 'btur']:
            # ... (Existing Link Logic - only works for exactly 2) ...
            # Identify Source vs Target (Clicked is Target)
            source_comp = selected_comps[0] if selected_comps[0] != clicked_comp else selected_comps[1]
            
            m.addSection("Parameter Linking")
            
            expr_val = f"p['{source_comp.uuid}'].{param_attr}"
            act_link = QAction(f"Link {col_name} to {source_comp.series} (Value)", m)
            act_link.triggered.connect(
                lambda: self._toggle_constraint(
                    clicked_comp.uuid, param_attr, False, expr_val, source_comp.uuid
                )
            )
            m.addAction(act_link)
            
            if col_name == 'b':
                m_dep = self._get_mass(clicked_comp.series)
                m_src = self._get_mass(source_comp.series)
                if m_dep and m_src:
                    mass_ratio = np.sqrt(m_src / m_dep)
                    expr_therm = f"p['{source_comp.uuid}'].b * {mass_ratio:.4f}"
                    act_therm = QAction(f"Link {col_name} to {source_comp.series} (Thermal)", m)
                    act_therm.triggered.connect(
                        lambda: self._toggle_constraint(
                            clicked_comp.uuid, param_attr, False, expr_therm, source_comp.uuid
                        )
                    )
                    m.addAction(act_therm)
            m.addSeparator()

        # 2. Constraint Actions (Multi-aware)
        if col_name in ['z', 'logN', 'b', 'btur']:
            
            if is_multi:
                # --- Batch Actions ---
                act_freeze_sel = QAction(f"Freeze '{col_name}'{suffix}", m)
                act_freeze_sel.triggered.connect(
                    lambda: self._toggle_constraint_list(selected_comps, param_attr, False)
                )
                m.addAction(act_freeze_sel)

                act_free_sel = QAction(f"Unfreeze '{col_name}'{suffix}", m)
                act_free_sel.triggered.connect(
                    lambda: self._toggle_constraint_list(selected_comps, param_attr, True)
                )
                m.addAction(act_free_sel)
                
            else:
                # --- Single Item Context (Smarter toggles) ---
                is_frozen = False
                is_linked = False
                if self.current_session and self.current_session.systs.constraint_model:
                    c_map = self.current_session.systs.constraint_model.v2_constraints_by_uuid
                    if clicked_comp.uuid in c_map and param_attr in c_map[clicked_comp.uuid]:
                        cons = c_map[clicked_comp.uuid][param_attr]
                        if not cons.is_free: is_frozen = True
                        if cons.expression: is_linked = True
                
                if is_linked:
                    act = QAction(f"Unlink '{col_name}' (Set Free)", m)
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

            # Global Column Actions (Always available)
            m.addSeparator()
            act_freeze_col = QAction(f"Freeze All '{col_name}' (Column)", m)
            act_freeze_col.triggered.connect(lambda: self._bulk_toggle_constraint(param_attr, False))
            m.addAction(act_freeze_col)
            
            act_free_col = QAction(f"Unfreeze All '{col_name}' (Column)", m)
            act_free_col.triggered.connect(lambda: self._bulk_toggle_constraint(param_attr, True))
            m.addAction(act_free_col)
            
            m.addSeparator()

        # Set Resolution (Multi-aware via internal logic)
        act_resol = QAction(f"Set R...{suffix}", m)
        act_resol.triggered.connect(self._set_component_resolution)
        m.addAction(act_resol)

        m.addSeparator()
        
        # 3. Operations (Multi-aware)
        act_refit = QAction(f"Refit Selected{suffix}", m)
        act_refit.triggered.connect(lambda: self._refit_list(selected_comps))
        m.addAction(act_refit)

        m.addAction("Refit All", self._refit_all)
        
        act_del = QAction(f"Delete{suffix}", m)
        act_del.triggered.connect(lambda: self._delete_list(selected_comps))
        m.addAction(act_del)

        m.exec(self.table_view.mapToGlobal(pos))

    # --- [NEW] Batch Helper Methods ---

    def _toggle_constraint_list(self, comps, param, is_free):
        """Iterates list and applies constraint update."""
        for c in comps:
            self._toggle_constraint(c.uuid, param, is_free, None, None)

    def _refit_list(self, comps):
        """Iterates list and triggers fit_component."""
        if self.main_window:
            for c in comps:
                # Note: This queues multiple recipes. 
                # Ideally we'd have a batch recipe, but this is safe and robust.
                self.main_window._on_recipe_requested(
                    "absorbers", "fit_component", {"uuid": c.uuid}, {}
                )

    def _delete_list(self, comps):
        """Batch delete with a single confirmation dialog."""
        if not comps or not self.main_window: return
        
        count = len(comps)
        if count == 1:
            c = comps[0]
            txt = f"Series: {c.series}\nRedshift: {c.z:.5f}"
        else:
            txt = f"You are about to delete {count} components."

        confirm = self.main_window._show_custom_message(
            title="Delete Components",
            header=f"Delete {count} component(s)?",
            text=txt,
            buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            default_btn=QMessageBox.StandardButton.No,
            parent=self
        )
        
        if confirm == QMessageBox.Yes:
            # --- [FIX] Send single batch request with ALL uuids ---
            uuid_list = [c.uuid for c in comps]
            
            self.main_window._on_recipe_requested(
                "absorbers", "delete_component", {"uuids": uuid_list}, {}
            )

    # Method to set resolution for specific components
    def _set_component_resolution(self):
        indexes = self.table_view.selectionModel().selectedRows()
        if not indexes: return

        # [NEW] Determine initial value to show
        initial_val = 0.0
        
        # 1. Try fetching from the first selected component
        first_src = self.proxy_model.mapToSource(indexes[0])
        first_comp = self.table_model.get_component_at(first_src.row())
        
        if first_comp and first_comp.resol is not None and first_comp.resol > 0:
            initial_val = first_comp.resol
        else:
            # 2. Fallback to Global Session Resolution
            if self.current_session:
                # Try scalar attribute
                if hasattr(self.current_session.spec._data, 'resol'):
                    r = self.current_session.spec._data.resol
                    if r > 0: initial_val = r
                
                # Try Metadata if still 0
                if initial_val == 0 and 'resol' in self.current_session.spec.meta:
                    try:
                        initial_val = float(self.current_session.spec.meta['resol'])
                    except: pass

        # 1. Create Custom Dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Set Component Resolving Power")
        dialog.setMinimumWidth(300)
        
        # Apply strict styling to match other forms
        pal = QApplication.palette()
        text_col = pal.color(pal.ColorRole.Text).name()
        base_col = pal.color(pal.ColorRole.Base).name()
        dialog.setStyleSheet(f"""
            QLineEdit {{
                padding: 4px;
                border-radius: 4px;
                background-color: {base_col};
                color: {text_col};
            }}
            QLineEdit:focus {{
                border: 1px solid #296bff; /* Highlight color */
            }}
        """)
        layout = QVBoxLayout(dialog)
        layout.setSpacing(10)

        # 2. Form Layout (Label Left, Box Right)
        form_layout = QHBoxLayout()
        label = QLabel("Resolving power R (e.g. 50000):")
        
        self.res_input = QLineEdit()
        if initial_val > 0:
            self.res_input.setText(f"{initial_val:.0f}")
        else:
            #self.res_input.setPlaceholderText("e.g. 50000")
            pass
        validator = QDoubleValidator()
        validator.setBottom(0.0)
        validator.setNotation(QDoubleValidator.StandardNotation)
        self.res_input.setValidator(validator)
        
        form_layout.addWidget(label)
        form_layout.addWidget(self.res_input)
        layout.addLayout(form_layout)
        
        # 3. Standard Buttons
        btn_box = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        btn_box.accepted.connect(dialog.accept)
        btn_box.rejected.connect(dialog.reject)
        layout.addWidget(btn_box)

        # 4. Execute
        if dialog.exec() == QDialog.Accepted:
            try:
                val = float(self.res_input.text())
                
                # Update each selected component
                for idx in indexes:
                    src = self.proxy_model.mapToSource(idx)
                    comp = self.table_model.get_component_at(src.row())
                    if comp:
                         if self.main_window:
                             self.main_window._on_recipe_requested(
                                 "absorbers", "update_component", 
                                 {'uuid': comp.uuid, 'resol': val}, {}
                             )
            except ValueError:
                pass

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
    
    def _bulk_toggle_constraint(self, param, is_free):
        if self.main_window:
            self.main_window._on_recipe_requested(
                "absorbers", "update_constraints_column",
                {"param": param, "is_free": is_free},
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
    
    def _refit_all(self):
        if self.main_window:
            # Trigger the newly created recipe
            self.main_window._on_recipe_requested(
                "absorbers", "refit_all", {}, {}
            )
    
    def _update_group_definition(self, selected_comps: List[ComponentDataV2]):
        if not selected_comps or not self.current_session:
            # Clear Highlights
            self.table_model.set_highlighted_uuids(set(), set())
            self.proxy_model.set_allowed_uuids(set())
            return

        # 1. Identify Primary UUIDs (The ones explicitly selected)
        primary_uuids = {c.uuid for c in selected_comps}
        
        # 2. Identify Full Group
        group_list = self.current_session.systs.get_connected_group(list(primary_uuids))
        all_group_uuids = set(group_list)
        
        # 3. Identify Secondary UUIDs (Group members NOT selected)
        secondary_uuids = all_group_uuids - primary_uuids

        # 4. Update Model Colors
        self.table_model.set_highlighted_uuids(primary_uuids, secondary_uuids)
        
        # 5. Update Proxy Filter (Show entire group)
        self.proxy_model.set_allowed_uuids(all_group_uuids)

    def _update_trans_box_visuals(self):
        """Updates the text color of the Trans box to indicate Auto vs Manual mode."""
        pal = QApplication.palette()
        base_col = pal.color(pal.ColorRole.Base).name()
        text_col = pal.color(pal.ColorRole.Text)
        
        # Calculate colors
        if self._manual_trans_edit:
            # Manual Mode: Standard Full Color
            col_str = text_col.name()
            font_weight = "bold"
        else:
            # Auto Mode: Dimmed (50% Alpha)
            # We can't use hex alpha easily in all Qt versions, so we use rgba
            col_str = f"rgba({text_col.red()},{text_col.green()},{text_col.blue()},0.6)"
            font_weight = "normal"

        self.trans_in.setStyleSheet(f"""
            QLineEdit {{
                padding: 3px;
                border-radius: 4px;
                background: {base_col};
                color: {col_str};
                font-weight: {font_weight};
            }}
            QLineEdit:focus {{
                border: 1px solid #296bff;
                color: {text_col.name()}; /* Always bright when focused for editing */
            }}
        """)