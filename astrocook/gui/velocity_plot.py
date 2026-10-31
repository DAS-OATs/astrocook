"""The velocity plot of the System Inspector: the panels of the transitions
of the selected components, their ticks, the profile under the cursor and
the menu of the right click; and ``calc_voigt_profile``, the profiles it
draws.

Moved out of ``system_inspector.py`` in phase 5 (split of the GUI).
"""
import numpy as np
import astropy.units as au
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.layout_engine import ConstrainedLayoutEngine
from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction, QCursor
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QMenu, QPushButton, QSizePolicy,
    QScrollBar, QStackedWidget, QToolTip,
)
from typing import List

from ..core.structures import ComponentDataV2
from ..core.atomic_data import STANDARD_MULTIPLETS, xem_d, ATOM_DATA
from ..core.resolution import resolution_at
from ..fitting.lsf import C_KMS as C_KMS_INSPECTOR, PixelLSF
from ..fitting.voigt_fitter import voigt_optical_depth
from .inspector_table import flag_cue
from .plot_style import ERROR_COLOR, get_color_cycle, get_style_color
from .theme import SELECTION
from .floating_cards import WrappedLabel
from .plot_toolbar import AstrocookToolbar


def calc_voigt_profile(wave_grid_ang, lambda_0, f_val, gamma, z, N, b_kms, resol=None,
                       resol_unit='R', context="Unknown", edges_ang=None):
    """Normalised profile of one line on ``wave_grid_ang`` (Angstrom), with
    the physics and the convolution of the fitter (#526): the optical depth
    on a fine grid, each point convolved with its own resolution (``resol``,
    scalar or one per point, R or FWHM in km/s by ``resol_unit``; None for
    none) and averaged over its pixel. ``edges_ang`` gives the edges of the
    pixels (xmin, xmax); by default they are halfway between the points."""
    if wave_grid_ang is None or len(wave_grid_ang) == 0:
        return np.array([])
    if hasattr(lambda_0, 'value'):
        lambda_0 = lambda_0.value
    if hasattr(wave_grid_ang, 'value'):
        wave_grid_ang = wave_grid_ang.value
    if hasattr(b_kms, 'value'):
        b_kms = b_kms.value
    x = np.asarray(wave_grid_ang, dtype=float)

    if edges_ang is not None:
        xmin, xmax = (np.asarray(e, dtype=float) for e in edges_ang)
    elif len(x) > 1:
        mid = 0.5 * (x[1:] + x[:-1])
        xmin = np.concatenate([[x[0] - (mid[0] - x[0])], mid])
        xmax = np.concatenate([mid, [x[-1] + (x[-1] - mid[-1])]])
    else:
        xmin, xmax = x * (1 - 1e-6), x * (1 + 1e-6)

    r = None
    if resol is not None:
        values = np.asarray(resol, dtype=float)
        with np.errstate(invalid='ignore', divide='ignore'):
            r = C_KMS_INSPECTOR / values if resol_unit != 'R' else values
        r = np.where(np.isfinite(r) & (r > 0), r, np.nan)
        if not np.any(np.isfinite(r)):
            r = None
        elif r.ndim == 0:
            r = float(r)
    r_fin = None if r is None else np.atleast_1d(r)[np.isfinite(np.atleast_1d(r))]
    lsf = PixelLSF(xmin, xmax,
                   r_min=None if r is None else float(np.min(r_fin)),
                   r_max=None if r is None else float(np.max(r_fin)),
                   b_min=float(b_kms), cache=False)
    tau = voigt_optical_depth(lsf.x_fine, np.array([lambda_0]), np.array([f_val]),
                              np.array([gamma]), np.array([z]), np.array([N]),
                              np.array([b_kms]))
    return lsf.apply(np.exp(-tau), r)

# --- 2. The Paged Plot Widget ---

#: The toolbar of the velocity plot shown or hidden (its tab), in the settings
INSPECTOR_TOOLBAR_KEY = "inspector_show_toolbar"

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
        
        # The plot, or a message saying why there is none: no continuum, or
        # no systems yet (visual review of phase 5-bis: the panels were
        # empty, or blank, without a word)
        self.message = QWidget()
        message_layout = QVBoxLayout(self.message)
        message_layout.addStretch()
        self.message_label = WrappedLabel()
        self.message_label.setWordWrap(True)
        self.message_label.setAlignment(Qt.AlignCenter)
        self.message_label.setFixedWidth(400)
        message_layout.addWidget(self.message_label, 0, Qt.AlignHCenter)
        self.continuum_button = QPushButton("Auto-estimate Continuum...")
        self.continuum_button.clicked.connect(self._estimate_continuum)
        message_layout.addWidget(self.continuum_button, 0, Qt.AlignHCenter)
        message_layout.addStretch()
        self.plot_stack = QStackedWidget()
        self.plot_stack.addWidget(self.canvas)
        self.plot_stack.addWidget(self.message)
        self.left_layout.addWidget(self.plot_stack)
        # The toolbar floats over the plot, as in the main window; shown or
        # hidden as the user left it (its own setting: the inspector has no
        # View menu)
        self.toolbar.float_over(self.plot_stack)
        self.toolbar.set_shown(QSettings("Astrocook", "AstrocookV2").value(
            INSPECTOR_TOOLBAR_KEY, True, type=bool))
        self.toolbar.on_toggled = lambda shown: QSettings("Astrocook", "AstrocookV2").setValue(
            INSPECTOR_TOOLBAR_KEY, shown)

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
        # Positions of the axes and of the common labels laid out by the
        # constrained layout, by (panels, residuals, size of the figure, x
        # limits): reused instead of laying out each draw again (half of a
        # draw on J0103-1305), see _update_plot
        self._layout_cache = {}
        self._layout_key = None
        
        self._current_session = None
        self._current_component = None
        self._selected_components = []
        self._plot_center_z = 0.0

        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('draw_event', self.on_draw_update_limits)
        self.canvas.mpl_connect('draw_event', self._keep_layout)
        
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
        session = (session_state or self._current_session
                   or getattr(self.inspector, 'current_session', None))
        comp = self._current_component
        if session and comp:
            self.plot_system(session, comp, self._selected_components)
        else:
            # No component: the transitions and the redshift typed above,
            # to look at the spectrum before any system is added
            self._current_session = session
            self._selected_components = []
            self._expand_transitions()
            self._update_scrollbar_range()
            self._update_plot()

    def _expand_transitions(self):
        """The transitions of the box above the plot, multiplets expanded."""
        expanded = []
        for item in self.inspector.active_transitions:
            if item in xem_d:
                expanded.append(item)
            elif item in STANDARD_MULTIPLETS:
                expanded.extend(STANDARD_MULTIPLETS[item])
        seen = set()
        self._all_transitions = [x for x in expanded if not (x in seen or seen.add(x))]

    #: The messages shown instead of the plot
    NO_CONTINUUM = ("No continuum. The velocity plot shows the flux divided by the "
                    "continuum: estimate it first.")
    NO_SYSTEMS = ("No absorption systems yet. Type transitions and a redshift above to "
                  "look at the spectrum, then add components from the menu of the plot "
                  "(right click), or run Absorbers > Auto-detect Doublets...")

    def _message(self):
        """The message to show instead of the plot, or None."""
        session = self._current_session or getattr(self.inspector, 'current_session', None)
        if session is None or session.spec is None:
            return None
        if session.spec.norm is None:
            return self.NO_CONTINUUM
        if not self._all_transitions and not (session.systs and session.systs.components):
            return self.NO_SYSTEMS
        return None

    def _show_message(self, text):
        if text is None:
            self.plot_stack.setCurrentWidget(self.canvas)
            return
        self.message_label.setText(text)
        self.continuum_button.setVisible(text == self.NO_CONTINUUM)
        self.plot_stack.setCurrentWidget(self.message)

    def _estimate_continuum(self):
        if self.main_window:
            self.main_window._launch_recipe_dialog("continuum", "estimate_auto")

    def set_velocity_limits(self, v_min, v_max):
        self._xlim = (v_min, v_max)
        self._update_plot()

    def set_center_redshift(self, z):
        self._plot_center_z = z
        width = self._xlim[1] - self._xlim[0]
        self._xlim = (-width/2, width/2)
        self.inspector.update_limit_boxes(self._xlim[0], self._xlim[1])
        if self._current_component is None:
            self.plot_spectrum()  # the transitions typed, without a component
        else:
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

        self._expand_transitions()
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
        Resolving power of the spectrum at a wavelength (nm), from its
        ``resol`` column (#526; before, the column was always read as km/s
        and the legacy sources as R). Returns (R, 'R'), or (0.0, 'R').
        """
        sess = self._current_session
        if not sess or not sess.spec: return 0.0, 'R'
        try:
            x_target = (target_wavelength_nm * au.nm).to_value(sess.spec.x.unit)
        except au.UnitConversionError:
            return 0.0, 'R'
        r = resolution_at(sess.spec, np.array([x_target]))[0]
        return (float(r), 'R') if np.isfinite(r) else (0.0, 'R')

    def _update_plot(self):
        self.fig.clear()
        self.axes = []
        self._panel_map = {} 
        self._tick_map_visuals = [] 

        # The canvas is drawn with draw_idle: opening the inspector or
        # changing its size updates the plot more than once in a row, and the
        # draws (~0.1 s each on J0103-1305) are coalesced into one
        self._show_message(self._message())
        if not self._all_transitions or not self._current_session:
            self.canvas.draw_idle(); return

        start = self.scrollbar.value() if self.scrollbar.isEnabled() else 0
        end = min(start + self.PAGE_SIZE, len(self._all_transitions))
        visible_trans = self._all_transitions[start:end]
        
        num_trans = len(visible_trans)
        if num_trans == 0: self.canvas.draw_idle(); return
        
        show_resid = self.inspector.resid_cb.isChecked()
        
        if show_resid:
            gs_ratios = []
            for _ in range(num_trans): gs_ratios.extend([3, 1])
            axs = self.fig.subplots(nrows=num_trans*2, ncols=1, sharex=True, 
                                    gridspec_kw={'height_ratios': gs_ratios, 'hspace': 0.05})
        else:
            axs = self.fig.subplots(nrows=num_trans, ncols=1, sharex=True)
            if num_trans == 1: axs = [axs]

        # The layout depends only on the number of panels, the residuals, the
        # size of the figure and the x limits (the tick labels): laid out
        # by matplotlib the first time, then reused (_apply_layout)
        self._layout_key = (num_trans, show_resid, tuple(self.fig.bbox.size),
                            self.fig.dpi, tuple(self._xlim))
        self.fig.set_layout_engine(
            layout='none' if self._layout_key in self._layout_cache else 'constrained')
        
        session = self._current_session
        x_full = session.spec.x.value
        c_kms = 299792.458
        z_sys = self._plot_center_z
        # Rest wavelengths in Angstrom, each converted once: the ticks of
        # every panel go through all the lines of all the components
        # (~1000 astropy conversions per panel on J0103-1305, half of a
        # rebuild, before phase 5)
        rest_ang = {}

        def rest_wavelength_ang(name):
            if name not in rest_ang:
                rest_ang[name] = xem_d[name].to(au.Angstrom).value
            return rest_ang[name]

        # Resolving power per pixel (#526), None if unknown
        resol_col = session.spec.resolution()
        if not np.isfinite(resol_col).any():
            resol_col = None

        has_norm = session.spec.norm is not None
        y = session.spec.y.value
        dy = session.spec.dy.value if session.spec.dy is not None else np.ones_like(y)
        norm = session.spec.norm.value if has_norm else np.ones_like(y)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            y_norm = np.divide(y, norm, out=np.full_like(y, np.nan), where=norm!=0)
            dy_norm = np.divide(dy, norm, out=np.full_like(dy, np.nan), where=norm!=0)
            
        y_resid = None
        if session.spec.model is not None:
            mod = session.spec.model.value
            with np.errstate(divide='ignore', invalid='ignore'):
                mod_norm = np.divide(mod, norm, out=np.full_like(mod, np.nan), where=norm!=0) if has_norm else mod
                y_resid = np.divide(y_norm - mod_norm, dy_norm, out=np.full_like(y_norm, np.nan), where=dy_norm!=0)

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
                    # As the coordinates of the main plot (status bar,
                    # phase 5-bis): ƒ is the normalized flux there too
                    return f"Δv = {x:.1f} km/s, ƒ = {y:.3f}, z = {val_z:.6f}"
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

            if not found_valid_col_data:
                # Fallback to Global (Metadata)
                panel_res_arg, panel_res_unit = self._get_resolution_at(lam_obs)

            # --- Error Plotting Logic ---
            # Check toggle state from Main Window
            show_err_as_line = False
            if self.inspector.main_window and hasattr(self.inspector.main_window, 'error_line_checkbox'):
                show_err_as_line = self.inspector.main_window.error_line_checkbox.isChecked()

            if show_err_as_line:
                # OPTION A: Plot Error as a separate line (e.g. Orange)
                ax_main.step(v_p, dy_p, where='mid', color=ERROR_COLOR, lw=1.0, alpha=0.8)
            else:
                # OPTION B: Traditional Band (Original)
                ax_main.step(v_p, y_p-dy_p, where='mid', color=ERROR_COLOR, lw=0.3)
                ax_main.step(v_p, y_p+dy_p, where='mid', color=ERROR_COLOR, lw=0.3)

            ax_main.step(v_p, y_p, where='mid', color=get_style_color('flux', colors), lw=0.8)
            
            if session.spec.model is not None:
                mod_p = mod_norm[mask]
                ax_main.plot(v_p, mod_p, color=get_style_color('model', colors), lw=1.6, alpha=0.8)

            # --- Plot Components ---
            if self._selected_components:
                atom_info = ATOM_DATA.get(trans_name)
                if atom_info:
                    x_ang_p = session.spec.x.to(au.Angstrom).value[mask]
                    edges_p = (session.spec.xmin.to(au.Angstrom).value[mask],
                               session.spec.xmax.to(au.Angstrom).value[mask])
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
                            
                            # Check if component has specific resolution FIRST
                            c_res = getattr(sel_c, 'resol', None)
                            if c_res is not None and c_res > 0:
                                comp_res_arg = c_res
                                comp_res_unit = 'R'  # always R (#526)

                            b_eff = np.sqrt(sel_c.b**2 + sel_c.btur**2)

                            prof = calc_voigt_profile(
                                x_ang_p, atom_info['wave'], atom_info['f'], atom_info['gamma'],
                                sel_c.z, 10**sel_c.logN, b_eff, 
                                resol=comp_res_arg, resol_unit=comp_res_unit,
                                context="STATIC_PLOT", edges_ang=edges_p
                            )
                            ax_main.plot(v_p, prof, color=SELECTION, ls='--', lw=1.2, alpha=0.9)

            # ... (Rest of plotting logic: Ticks, Labels, Residuals) ...
            tick_ymin, tick_ymax = 0.02, 0.08 
            trans_axis = ax_main.get_xaxis_transform()
            lam_rest_panel_ang = rest_wavelength_ang(trans_name)
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
                        lam_rest_c = rest_wavelength_ang(line_name)
                        lam_obs_c = lam_rest_c * (1 + other_c.z)
                        v_shift = c_kms * (lam_obs_c - lam_obs_panel_center_ang) / lam_obs_panel_center_ang

                    if v_min <= v_shift <= v_max:
                        self._tick_map_visuals.append((ax_main, v_shift, other_c))
                        is_h = other_c.series.startswith('Ly') or other_c.series.startswith('H')
                        col = 'red' if is_h else 'gray'
                        # Flags of the fit (#517): dashed, in the colour of the cue
                        flag_col, _ = flag_cue(other_c)
                        ax_main.plot([v_shift, v_shift], [tick_ymin, tick_ymax],
                                transform=trans_axis, color=flag_col or col,
                                ls='--' if flag_col else '-',
                                lw=1.5 if (is_h or flag_col) else 1.0, alpha=0.9 if flag_col else 0.8)
            
            ax_main.axvline(0, color='gray', ls=':', lw=0.8)
            ax_main.axhline(1.0, color='gray', ls=':', lw=0.8)
            ax_main.axhline(0.0, color='gray', ls=':', lw=0.8)

            ax_main.text(0.98, 0.85, trans_name, transform=ax_main.transAxes, ha='right', fontweight='bold', 
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            ax_main.set_xlim(self._xlim)
            ax_main.set_ylim(-0.1, 1.4) 

            if show_resid and y_resid is not None:
                resid_p = y_resid[mask]
                ax_resid.fill_between(v_p, - 1.0, 1.0, step='mid', color=ERROR_COLOR, alpha=0.15)
                ax_resid.step(v_p, resid_p, where='mid', color='red', lw=0.8)
                ax_resid.axhline(0, color=get_style_color('model', colors), lw=1.0, alpha=0.5)
                ax_resid.set_ylim(-5, 5)
                ax_resid.set_ylabel(r"χ")
                ax_main.tick_params(labelbottom=False)

        self.fig.supxlabel("Relative velocity (Δv, km/s)")
        if not show_resid:
            self.fig.supylabel("Normalized Flux")
        self._apply_layout()
        self.canvas.draw_idle()

    def _apply_layout(self):
        """Put the axes and the common labels where the constrained layout
        put them for the same key, if it was laid out already."""
        cached = self._layout_cache.get(self._layout_key)
        if cached is None:
            return
        positions, labels = cached
        for ax, pos in zip(self.fig.axes, positions):
            ax.set_position(pos)
        for name, pos in labels.items():
            label = getattr(self.fig, name, None)
            if label is not None:
                label.set_position(pos)

    def _keep_layout(self, event):
        """After a draw laid out by matplotlib, keep the layout: store it
        for its key and stop laying out the following draws (the cursor
        profile redraws the figure at each move of the mouse)."""
        if self._layout_key is None or self.fig.get_layout_engine() is None:
            return
        if not isinstance(self.fig.get_layout_engine(), ConstrainedLayoutEngine):
            return
        positions = [ax.get_position() for ax in self.fig.axes]
        labels = {name: getattr(self.fig, name).get_position()
                  for name in ('_supxlabel', '_supylabel')
                  if getattr(self.fig, name, None) is not None}
        self._layout_cache[self._layout_key] = (positions, labels)
        self.fig.set_layout_engine(layout='none')

    def on_mouse_move(self, event):
        if not event.inaxes or not self._current_session: return
        
        tooltip_text = ""
        view_width = self._xlim[1] - self._xlim[0]
        tol_v = view_width * 0.015 
        
        for ax, v_pos, comp in self._tick_map_visuals:
            if event.inaxes == ax and abs(event.xdata - v_pos) < tol_v:
                c_chi2 = f"{comp.chi2:.2f}" if comp.chi2 is not None else "N/A"
                c_resol = f"{comp.resol:.0f}" if comp.resol is not None else "spectrum"
                tooltip_text = (
                    f"<div style='font-size: 13px'>"
                    f"<b>ID:</b> {comp.id} | <b>{comp.series}</b><br>"
                    f"<b>z:</b> {comp.z:.5f}<br>"
                    f"<b>\u03c7\u00b2:</b> {c_chi2} | <b>R:</b> {c_resol}"
                )
                _, messages = flag_cue(comp)
                if messages:
                    tooltip_text += "<hr style='margin: 4px 0;'>" + "<br>".join(
                        f"\u26a0 {m}" for m in messages)
                tooltip_text += "</div>"
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
        # Resolving power per pixel (#526), None if unknown
        resol_col = session.spec.resolution()
        if not np.isfinite(resol_col).any():
            resol_col = None

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
                    dy_norm = np.divide(dy, norm, out=np.full_like(dy, np.nan), where=norm!=0)
                
                dy_interp = np.interp(v_grid, v_full, dy_norm, left=1.0, right=1.0)
                with np.errstate(divide='ignore', invalid='ignore'):
                    prof_sigma = np.divide(total_prof - 1.0, dy_interp, out=np.zeros_like(total_prof), where=dy_interp>1e-9)
                    prof_sigma[~np.isfinite(prof_sigma)] = 0.0
                line_resid.set_data(v_grid, prof_sigma)
            
        self.canvas.draw_idle()
    
    def on_press(self, event):
        """
        Handles mouse clicks on the velocity plot.
        Right-Click (Button 3):
        1. If clicked on a Tick -> Focus on that component (Only group).
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
                # "Focus" triggers the Only group logic in the Inspector
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
                        'logN': self.cursor_logN, 'b': self.cursor_b,
                        **self.inspector.fit_window_params()
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
                                'logN': self.cursor_logN, 'b': self.cursor_b,
                                **self.inspector.fit_window_params()
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
                                        'b': self.cursor_b,
                                        **self.inspector.fit_window_params()
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
                                'logN': self.cursor_logN, 'b': self.cursor_b,
                                **self.inspector.fit_window_params()
                            }, {}
                        )
                    )
                    menu.addAction(act_link)

                has_actions = True

            # The components are fitted when added: not without a continuum
            if not self.inspector.has_continuum():
                menu.setToolTipsVisible(True)
                for action in menu.actions():
                    if action.text().startswith("Add "):
                        action.setEnabled(False)
                        action.setToolTip(
                            "Needs a continuum (Continuum > Auto-estimate Continuum...)")

            # --- Execute ---
            if has_actions: menu.exec(QCursor.pos())
