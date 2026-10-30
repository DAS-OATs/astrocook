import astropy.constants as const
import astropy.units as au
import logging
import matplotlib
# No matplotlib.use('QtAgg') and no pyplot here: the canvases of the GUI are
# QtAgg whatever the backend of pyplot, which the launcher sets (it failed
# at import when pyplot had already picked a headless backend, F7)
import matplotlib.style
import matplotlib.ticker as plt_ticker
import numpy as np
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QVBoxLayout, QWidget
from dataclasses import dataclass
from typing import NamedTuple, Optional, TYPE_CHECKING

from astrocook.core.atomic_data import STANDARD_MULTIPLETS, trans_parse, xem_d
from astrocook.core.spectrum_operations import convert_axis_velocity

# Use TYPE_CHECKING to avoid circular import errors at runtime
if TYPE_CHECKING:
    from astrocook.core.session import SessionV2 

from .continuum_editor import ContinuumEditorMixin
from .plot_canvas import MatplotlibCanvas
from .plot_mouse import PlotMouseMixin
from .plot_style import (STRONG_EMISSION_LINES, PANEL_LEGEND_LABELS, decimate_x_for_min_max,
                         decimate_y_min_max, get_color_cycle)
from .plot_toolbar import AstrocookToolbar

#: When more pixels than this are in view, they are decimated (min-max) to
#: keep the plot fast (F76: at every plot, in the full view or in a wide
#: zoom). 1 million (Guido, phase 5-bis; it was 2 million, and the comment
#: of the first version said "> 20k points"): a redraw of 3 million pixels
#: took 3.9 s whole and 1.0 s decimated, of 300 000 pixels 0.77 and 0.57 s
#: (with the legend in a fixed corner, see F86)
DECIMATION_THRESHOLD = 1000000
DECIMATION_FACTOR = 10       # Plot 2 points for every 10

# Default limits of the normalized flux: linear, and with the logarithmic
# F scale (F78), which needs a positive lower limit (saturated lines are
# shown down to 1% of the continuum)
NORM_YLIM = (-0.3, 1.3)
NORM_YLIM_LOG = (0.01, 1.3)

# Symbols of the labels and of the coordinates under the mouse
SYM_LAMBDA = "λ"      # λ
SYM_LAMBDA_R = "λᵣ" # λᵣ (lambda sub r)
SYM_FLUX = "F"
SYM_FLUX_NORM = "ƒ"        # ƒ (hook f for flux/function)
SYM_ANGSTROM = "Å"    # Å
SYM_MICRON = "µm"     # µm


class _PlotArrays(NamedTuple):
    """The arrays of a plot of the spectrum, as shown (as they are,
    normalized, or divided by the error): wavelength, flux, pixel size,
    error, continuum and model (the last two None when missing)."""
    x: np.ndarray
    y: np.ndarray
    dx: np.ndarray
    dy: np.ndarray
    norm: Optional[np.ndarray]
    model: Optional[np.ndarray]

    def sliced(self, s: slice) -> '_PlotArrays':
        return _PlotArrays(self.x[s], self.y[s], self.dx[s], self.dy[s],
                           None if self.norm is None else self.norm[s],
                           None if self.model is None else self.model[s])

    def decimated(self, factor: int) -> '_PlotArrays':
        """Min-max decimation: 2 points every ``factor``."""
        return _PlotArrays(decimate_x_for_min_max(self.x, factor),
                           *(decimate_y_min_max(a, factor) for a in self[1:]))


@dataclass
class _PlotView:
    """What a plot of the spectrum is drawn with: the view settings of the
    main window and the limits of the previous plot."""
    is_norm_y: bool
    is_snr: bool
    is_log_x: bool
    is_log_y: bool
    x_unit: str
    snr_col: str  # the error column of the SNR view ("dy" for dF)
    was_zoomed: bool  # the limits of the previous plot are kept...
    force_autoscale: bool  # ... unless Home was pressed
    previous_xlim: tuple
    previous_ylim: tuple
    norm_state_changed: bool = False  # normalized or SNR switched since the last plot
    use_decimation: bool = False
    final_slice: slice = slice(None)  # the cut of the arrays to the zoom
    is_scatter: bool = False
    scatter_ms: float = 4.0

    def set_marker_style(self, spec, n_points: int):
        """Points instead of steps for the spectra that ask for them (the
        'plot_style' meta key), when they are not too many, larger when
        they are fewer."""
        plot_style = spec.meta.get('plot_style', 'step')
        if n_points < 20000:
            self.is_scatter = (plot_style == 'scatter')
            self.scatter_ms = 4.0
            if n_points < 600:
                self.scatter_ms = 10.0
            elif n_points < 2000:
                self.scatter_ms = 8.0
            elif n_points < 6000:
                self.scatter_ms = 6.0


def _divide(a, b):
    """a / b, NaN where b is 0."""
    return np.divide(a, b, out=np.full_like(a, np.nan), where=b!=0)


def _full_plot_arrays(spec) -> _PlotArrays:
    """The columns of the spectrum to plot, without units (pixel size and
    error 0 when missing)."""
    full_x_data = spec.x.value
    full_y_data = spec.y.value
    # We check that NEITHER of them is None before accessing .value
    full_dx_data = (
                    spec.xmax.value - spec.xmin.value
                    if spec.xmin is not None and spec.xmax is not None
                    else np.zeros_like(full_x_data)
                 )
    full_dy_data = spec.dy.value if spec.dy is not None else np.zeros_like(full_y_data)
    full_norm_data = spec.norm.value if hasattr(spec.norm, 'value') else spec.norm
    full_model_data = spec.model.value if hasattr(spec.model, 'value') else spec.model
    return _PlotArrays(full_x_data, full_y_data, full_dx_data, full_dy_data,
                       full_norm_data, full_model_data)


def _coordinate_formatter(spec, view: _PlotView, z_em, x_label_unit):
    """The text of the coordinates under the mouse (in the toolbar):
    wavelength, flux and rest-frame wavelength, in the units shown."""
    y_label_symbol = SYM_FLUX_NORM if view.is_norm_y else SYM_FLUX

    # The unit of y (none when normalized or SNR), with Unicode symbols
    y_label_unit = ""
    if view.is_norm_y:
        y_label_unit = ""
    elif view.is_snr:
        y_label_unit = ""
    elif spec.y.unit:
        unit_str = f"{spec.y.unit:s}"
        unit_str = unit_str.replace("Angstrom", SYM_ANGSTROM)
        unit_str = unit_str.replace("**", "^") # Es. cm^2
        y_label_unit = unit_str

    selected_x_unit = view.x_unit
    is_norm_y, is_snr = view.is_norm_y, view.is_snr

    def new_format_coord(x, y):
        # x in nm (the unit of the data), y as shown (flux, normalized or SNR)
        x_disp = x
        if selected_x_unit == 'Angstrom': x_disp = x * 10.0
        elif selected_x_unit == 'micron': x_disp = x / 1000.0

        if is_norm_y or is_snr:
            y_fmt = f"{y:.3f}"
        else:
            y_fmt = f"{y:.3e}"

        # Format: λ=1215.67 Å, ƒ=1.23e-17 ...
        coord_str = f"{SYM_LAMBDA} = {x_disp:.4f} {x_label_unit}, {y_label_symbol} = {y_fmt} {y_label_unit}"

        # Rest frame
        if z_em > 0.0 and spec.x.unit.is_equivalent(au.nm):
            x_rf_nm = x / (1 + z_em)

            if selected_x_unit == 'Angstrom': x_rf_disp = x_rf_nm * 10.0
            elif selected_x_unit == 'micron': x_rf_disp = x_rf_nm / 1000.0
            else: x_rf_disp = x_rf_nm

            coord_str += f", {SYM_LAMBDA_R} = {x_rf_disp:.4f} {x_label_unit}"

        return coord_str

    return new_format_coord


class SpectrumPlotWidget(PlotMouseMixin, ContinuumEditorMixin, QWidget):
    """
    Refactored widget that contains the plot canvas, toolbar, and toggles.
    This replaces the content that was previously inside SpectrumViewerPySide.
    """
    def __init__(self, initial_session_state: Optional['SessionV2'], main_window_ref):
        super().__init__()
        self.main_window = main_window_ref

        self._paused_tool_mode = None
        
        # This style applies to any QToolTip triggered from this widget
        self.setStyleSheet("""
            QToolTip {
                border: 1px solid gray;
                border-radius: 5px;     /* The rounded corners */
                background-color: #F2F8F8;
                color: black;
                padding: 4px;           /* internal padding */
                background-clip: border-box;
            }
        """)

        # Main Vertical Layout for the plot area
        self.main_layout = QVBoxLayout(self)

        # 1. Canvas
        self.canvas = MatplotlibCanvas(parent=self, plot_widget=self)
        self.main_layout.addWidget(self.canvas, 1)

        # 2. Toolbar (NO Checkboxes here anymore)
        self.toolbar = AstrocookToolbar(self.canvas, self, select=True) # Pass self as parent

        # ** Store the last-drawn normalization state **
        self._last_draw_norm_state = False # Default to non-normalized
        self._last_draw_snr_state = False

        # The arrays of the last plots ('raw', 'norm', 'snr'), the data and the
        # error column of the SNR view they were computed for (_plot_arrays)
        self._cached_data_id = None # Tracks the ID of the cached data
        self._plot_cache = {}
        self._cached_snr_error_col = "" # Tracks which error col was used

        # List to store active iso-mag curves
        self._show_isomag_grid = False
        self.isomag_artists = []

        # Selection Mode State (the span being dragged is on the canvas)
        self._is_selecting_region = False
        
        # Store the original Matplotlib coordinate formatter ONCE
        self._default_format_coord = self.canvas.axes.format_coord

        # Tooltip state
        self._tip_component_html = None
        self._tip_data_rows = None # List of (label, value, unit) tuples
        self._tip_region_html = None

        # Add ONLY toolbar to main layout
        self.main_layout.addWidget(self.toolbar, 0) # Stretch 0

        # --- Continuum Editor State ---
        self._edit_mode_active = False
        self._knots_x = None  # np.array
        self._knots_y = None  # np.array
        self._is_anchor = None
        self._active_knot_idx = None # Index of knot being dragged
        
        self.knot_artist = None  # Matplotlib Line2D (markers)
        self.spline_artist = None # Matplotlib Line2D (smooth line)

        # --- Highlight State ---
        self._highlighted_uuids = set()
        # The highlighted systems of the last plot, (marker line, guide
        # lines) or None, and the pixels of the axes without them and the
        # legend: a new selection is drawn over these (_blit_highlights)
        self._highlight_artists = None
        self._highlight_base = None
        self._drawing_highlight_base = False
        self._plotted_state = None
        # Any full draw may change what is under the highlights
        self.canvas.mpl_connect('draw_event', self._on_full_draw)

        # Draw initial plot (will read checkbox state from main window)
        self.plot_spectrum(session_state=initial_session_state, initial_draw=True)

    def plot_spectrum(self, session_state: Optional['SessionV2'], initial_draw=False, force_autoscale: bool = False):
        """
        Retrieves data from the immutable V2 Session and plots it.

        The axes are cleared and drawn again, layer by layer (other
        sessions selected, error, flux, continuum, model, auxiliary column,
        systems, redshift cursor), keeping the zoom of the previous plot
        unless ``force_autoscale`` (the Home button).
        """
        selected_snr_col = self.main_window.snr_col_combo.currentText()
        # [CHANGE] Map "dF" display name back to internal "dy" name
        if selected_snr_col == "dF":
            selected_snr_col = "dy"
        self._invalidate_plot_caches(session_state, force_autoscale, selected_snr_col)

        if session_state is None:
            logging.debug("plot_spectrum called with no session manager. Clearing plot.")
            self.toolbar.set_select_enabled(False)
            ax = self.canvas.axes
            ax.clear()
            ax.set_title("No Spectrum Data Loaded")
            # Ensure background capture logic is handled if needed after clear
            self.canvas.background = None # Reset background
            self._connect_background_capture()
            # Trigger draw
            if initial_draw: self.canvas.draw()
            else: self.canvas.draw_idle()
            return # Exit early

        spec = session_state.spec
        ax = self.canvas.axes
        # The Select tool needs λ (F73)
        self.toolbar.set_select_enabled(spec.x.unit.is_equivalent(au.nm))

        # ** Store current limits BEFORE clearing **
        previous_xlim = ax.get_xlim()
        previous_ylim = ax.get_ylim()
        was_zoomed = self._was_zoomed(previous_xlim, previous_ylim, force_autoscale)

        ax.clear() # Now clear axes
        self._reconnect_limit_callbacks()
        legend_loc = self._legend_location(initial_draw)

        # ** Reset background and clear artists for blitting **
        self.canvas.background = None
        self.canvas.cursor_artists = []

        self._tick_visuals = []
        self._highlight_artists = None
        self._plotted_state = session_state

        plot_occurred = False
        if spec and len(spec.x) > 0:
            plot_occurred = True
            view = self._read_view(spec, selected_snr_col, was_zoomed, force_autoscale,
                                   previous_xlim, previous_ylim)
            full = _full_plot_arrays(spec)
            arrays = self._plot_arrays(spec, full, view)
            source = id(arrays)  # the cached arrays the cut is made from
            arrays = self._slice_to_view(arrays, view)
            arrays = self._decimate(arrays, view, source)
            view.set_marker_style(spec, len(arrays.x))
            colors = get_color_cycle(5, cmap='tab20') # Get colors

            self._plot_other_sessions(ax, view)
            self._plot_error(ax, arrays, view)
            self._plot_flux(ax, arrays, view, colors)
            self._plot_continuum_and_model(ax, arrays, colors)
            self._plot_aux_column(ax, spec, arrays, view)
            self._plot_systems(ax, spec, session_state.systs)
            self._plot_cursor(session_state, colors)

            self._add_legend(ax, legend_loc)
            ax.set_title(f"Spectrum: {session_state.name}")
            self._format_axes(ax, spec, view, full.x)
            self._set_limits(ax, view)

        # Draw the iso-mag grid last. It will self-manage its impact on limits.
        self._draw_isomag_grid()

        if not plot_occurred:
            ax = self.canvas.axes # Ensure ax is defined
            ax.clear() # Clear axes even if no data
            ax.set_title("No Spectrum Data Loaded")

        self._connect_background_capture()

        # Redraw the canvas
        if initial_draw:
            self.canvas.draw()
        else:
            self.canvas.draw_idle()

        # Update the limit boxes in the main window *after* plot is drawn
        if self.main_window:
            # Schedule this to run just after the draw completes
            QTimer.singleShot(0, self.main_window._update_limit_boxes_from_plot)

        # --- Refresh Data Inspector if it's open ---
        if hasattr(self, '_data_inspector_dialog') and self._data_inspector_dialog is not None:
            if self._data_inspector_dialog.isVisible():
                self._data_inspector_dialog.refresh(session_state)

    # --- Parts of plot_spectrum: state of the plot ---

    def _invalidate_plot_caches(self, session_state, force_autoscale, selected_snr_col):
        """Empty the caches of the plot arrays when the data changed, on
        Home, or when the error column of the SNR view changed."""
        # Get the unique ID of the immutable data core
        current_data_id = None
        if session_state and session_state.spec:
            # The ._data object is the frozen SpectrumDataV2
            current_data_id = id(session_state.spec._data)

        # If data changed OR the SNR error col changed, invalidate caches
        if (current_data_id != self._cached_data_id or
            force_autoscale or
            (self.main_window.snr_checkbox.isChecked() and
             self._cached_snr_error_col != selected_snr_col)
           ):

            if force_autoscale:
                logging.debug("Home button pressed, invalidating plot caches.")
            elif current_data_id != self._cached_data_id:
                logging.debug("Data change detected, invalidating ALL plot caches.")
            else:
                logging.debug("SNR error column changed, invalidating plot caches.")

            self._cached_data_id = current_data_id
            self._plot_cache = {}
            self._cached_snr_error_col = "" # Clear the cached col name

    def _connect_background_capture(self):
        """Capture the background for blitting at the next draw."""
        # Disconnect any previous connection first
        if self.canvas.draw_event_cid is not None:
            try: self.canvas.mpl_disconnect(self.canvas.draw_event_cid)
            except Exception: pass # Ignore if already disconnected
        self.canvas.draw_event_cid = self.canvas.mpl_connect('draw_event', self.canvas._capture_background)

    @staticmethod
    def _was_zoomed(previous_xlim, previous_ylim, force_autoscale) -> bool:
        """Whether the limits of the previous plot are to be kept."""
        if force_autoscale:
            logging.debug("Home button pressed: Forcing autoscale.")
            return False
        was_zoomed = False
        try:
            if previous_xlim != (0.0, 1.0) or previous_ylim != (0.0, 1.0):
                if all(isinstance(v,(int,float)) for v in previous_xlim+previous_ylim) and \
                    previous_xlim[0] != previous_xlim[1] and previous_ylim[0] != previous_ylim[1]:
                    was_zoomed = True
        except Exception: pass # Ignore if axes not ready
        return was_zoomed

    def _reconnect_limit_callbacks(self):
        """Replot when the limits change (the callbacks are lost when the
        axes are cleared)."""
        # Disconnect previous first, just in case
        if self.canvas.xlim_cid is not None:
            try: self.canvas.axes.callbacks.disconnect(self.canvas.xlim_cid)
            except Exception: pass
        if self.canvas.ylim_cid is not None:
            try: self.canvas.axes.callbacks.disconnect(self.canvas.ylim_cid)
            except Exception: pass
        # Reconnect
        try:
            self.canvas.xlim_cid = self.canvas.axes.callbacks.connect('xlim_changed', self.canvas.on_lim_changed)
            self.canvas.ylim_cid = self.canvas.axes.callbacks.connect('ylim_changed', self.canvas.on_lim_changed)
            logging.debug("Reconnected limit change callbacks in plot_spectrum.")
        except Exception as e:
            logging.error(f"Failed to reconnect limit callbacks: {e}")

    def _legend_location(self, initial_draw) -> str:
        """The corner of the legend away from the open sidebars."""
        legend_loc = 'best' # Default
        if not initial_draw:
            try:
                # Check width as a proxy for visibility/expanded state
                left_sidebar_open = self.main_window.left_sidebar_widget.width() > 1
                right_sidebar_open = self.main_window.right_sidebar_widget.width() > 1

                if left_sidebar_open and not right_sidebar_open:
                    legend_loc = 'upper right'
                elif not left_sidebar_open and right_sidebar_open:
                    legend_loc = 'upper left'
                elif left_sidebar_open and right_sidebar_open:
                    # Maybe center top if both are open? Or stick to one side.
                    legend_loc = 'upper center'
                # else: use 'best' if both closed
                logging.debug(f"Legend location set to: {legend_loc} (L:{left_sidebar_open}, R:{right_sidebar_open})")
            except Exception as e:
                logging.warning(f"Could not determine sidebar state for legend: {e}")
                legend_loc = 'best' # Fallback
        return legend_loc

    def _read_view(self, spec, selected_snr_col, was_zoomed, force_autoscale,
                   previous_xlim, previous_ylim) -> '_PlotView':
        """The view settings of the main window, and whether the
        normalization changed since the last plot."""
        mw = self.main_window
        view = _PlotView(
            is_norm_y=mw.norm_y_checkbox.isChecked(),
            is_snr=mw.snr_checkbox.isChecked(),
            is_log_x=mw.log_x_checkbox.isChecked(),
            is_log_y=mw.log_y_checkbox.isChecked(),
            x_unit=mw.x_unit_combo.currentText(),
            snr_col=selected_snr_col,
            was_zoomed=was_zoomed,
            force_autoscale=force_autoscale,
            previous_xlim=previous_xlim,
            previous_ylim=previous_ylim)

        # ** Check if normalization state CHANGED since last draw **
        view.norm_state_changed = (view.is_norm_y != self._last_draw_norm_state) or \
                                  (view.is_snr != self._last_draw_snr_state)
        if view.norm_state_changed:
            logging.debug(f"View state changed (Norm: {view.is_norm_y}, SNR: {view.is_snr})")
        # ** Update the stored state for the next draw **
        self._last_draw_norm_state = view.is_norm_y
        self._last_draw_snr_state = view.is_snr

        return view

    # --- Parts of plot_spectrum: the arrays ---

    def _raw_plot_arrays(self, full: '_PlotArrays', view: '_PlotView') -> '_PlotArrays':
        """The arrays of the spectrum (cached until the data change). Not
        decimated: the decimation depends on the view (``_decimate``, F76)."""
        cached = self._plot_cache.get('raw')
        if cached is not None:
            logging.debug("Using cached RAW plot data.")
            return cached
        logging.debug("Computing RAW plot data...")
        self._plot_cache['raw'] = full
        return full

    def _plot_arrays(self, spec, full: '_PlotArrays', view: '_PlotView') -> '_PlotArrays':
        """The arrays to plot: as they are, normalized to the continuum, or
        divided by the error column of the SNR view (cached until the data
        or the error column change)."""
        if view.is_norm_y:
            cached = self._plot_cache.get('norm')
            if cached is not None:
                logging.debug("Using cached NORMALIZED plot data.")
                return cached
            logging.debug("Cache miss for normalized data. Computing...")
            raw = self._raw_plot_arrays(full, view)
            if raw.norm is None:
                logging.warning("Normalize checked, but continuum data is missing.")
                return raw
            norm = raw.norm
            arrays = _PlotArrays(
                raw.x, _divide(raw.y, norm), _divide(raw.dx, norm), _divide(raw.dy, norm),
                np.ones_like(norm), None if raw.model is None else _divide(raw.model, norm))
            self._plot_cache['norm'] = arrays
            return arrays

        if view.is_snr:
            cached = self._plot_cache.get('snr')
            if cached is not None:
                logging.debug("Using cached SNR plot data.")
                return cached
            logging.debug("Cache miss for SNR data. Computing...")
            raw = self._raw_plot_arrays(full, view)

            # The ERROR array to use
            # The name of the column ("dy" for the "dF" of the box, F75)
            error_col_name = view.snr_col
            try:
                if error_col_name == 'dy':
                    full_error_data = full.dy
                else:
                    # This path is for 'running_std', 'norm_err', etc.
                    full_error_data = spec.get_column(error_col_name).value
            except Exception as e:
                logging.error(f"Could not get SNR error column '{error_col_name}': {e}. Falling back to 'dy'.")
                full_error_data = full.dy # Fallback
                self.main_window.snr_col_combo.setCurrentIndex(0) # Reset combo
            error = full_error_data

            y = _divide(raw.y, error)
            arrays = _PlotArrays(
                raw.x, y, raw.dx,
                np.ones_like(y), # Error on SNR is (dy/dy) = 1.0
                None if raw.norm is None else _divide(raw.norm, error),
                None if raw.model is None else _divide(raw.model, error))
            self._plot_cache['snr'] = arrays
            self._cached_snr_error_col = error_col_name # Store the col used
            return arrays

        return self._raw_plot_arrays(full, view)

    @staticmethod
    def _slice_to_view(arrays: '_PlotArrays', view: '_PlotView') -> '_PlotArrays':
        """The arrays cut to the zoom kept from the previous plot, with a
        margin; ``view.final_slice`` records the cut."""
        if not (view.was_zoomed and not view.force_autoscale):
            return arrays
        try:
            # Find indices for the visible range *on the plot data*
            idx_start = np.searchsorted(arrays.x, view.previous_xlim[0], side='left')
            idx_end = np.searchsorted(arrays.x, view.previous_xlim[1], side='right')

            # Add a buffer (in *indices*, not data units)
            buffer = 50
            idx_start = max(0, idx_start - buffer)
            idx_end = min(len(arrays.x), idx_end + buffer)

            if idx_end > idx_start:
                logging.debug(f"Applying zoom-slice to plot data: {idx_end - idx_start} points")
                view.final_slice = slice(idx_start, idx_end)
                return arrays.sliced(view.final_slice)
            logging.debug("Zoom slice is empty, plotting full data.")
        except Exception as e:
            logging.warning(f"Failed to apply zoom slice: {e}")
        return arrays

    def _decimate(self, arrays: '_PlotArrays', view: '_PlotView', source: int) -> '_PlotArrays':
        """Min-max decimation of the arrays in view, when they are more than
        ``DECIMATION_THRESHOLD`` (F76): at every plot, in the full view or
        in a zoom, after the cut to the view. The whole arrays stay in the
        cache; the decimated ones are cached by view (the arrays they come
        from, ``source``, and the cut), which saves ~0.2 s at each plot of a
        spectrum of 3 million pixels."""
        if len(arrays.x) <= DECIMATION_THRESHOLD:
            return arrays
        view.use_decimation = True
        key = ('decimated', source, view.final_slice.start, view.final_slice.stop)
        cached = self._plot_cache.get(key)
        if cached is None:
            logging.debug(f"Using Min-Max decimation (factor {DECIMATION_FACTOR}): "
                          f"{len(arrays.x)} points in view.")
            # Only the last view: the cuts of a pan would pile up
            for old in [k for k in self._plot_cache if k[0] == 'decimated']:
                del self._plot_cache[old]
            cached = self._plot_cache[key] = arrays.decimated(DECIMATION_FACTOR)
        return cached

    # --- Parts of plot_spectrum: the layers ---

    def _plot_other_sessions(self, ax, view: '_PlotView'):
        """The other sessions selected in the list ("ghosts"), in the same
        view as the active one."""
        if not (hasattr(self.main_window, 'session_list_view') and self.main_window.session_list_view.selectionModel()):
            return
        indexes = self.main_window.session_list_view.selectionModel().selectedIndexes()
        ghost_states = []
        for idx in indexes:
            row = idx.row()
            if 0 <= row < len(self.main_window.session_histories):
                hist = self.main_window.session_histories[row]
                # Exclude the primary active history
                if hist != self.main_window.active_history:
                    ghost_states.append(hist.current_state)
        if not ghost_states:
            return

        cmap = matplotlib.colormaps.get_cmap('tab10') # Distinct, colorblind-friendly palette
        for i, ghost_state in enumerate(ghost_states):
            if not ghost_state.spec: continue

            ghost_x = ghost_state.spec.x.value
            ghost_y = ghost_state.spec.y.value
            ghost_dy = ghost_state.spec.dy.value if ghost_state.spec.dy is not None else None

            # 1. Apply view toggles (Norm/SNR) to ghost so it matches primary scale
            if view.is_norm_y:
                ghost_norm = ghost_state.spec.norm.value if ghost_state.spec.norm is not None else np.ones_like(ghost_y)
                ghost_y = _divide(ghost_y, ghost_norm)
                if ghost_dy is not None:
                    ghost_dy = _divide(ghost_dy, ghost_norm)
            elif view.is_snr:
                try:
                    err_col = view.snr_col
                    err_data = ghost_dy if err_col == 'dy' else ghost_state.spec.get_column(err_col).value
                    ghost_y = _divide(ghost_y, err_data)
                    ghost_dy = np.ones_like(ghost_y)
                except: pass

            # 2. Slice ghost data to current zoom to save rendering time
            if view.was_zoomed and not view.force_autoscale:
                g_start = np.searchsorted(ghost_x, view.previous_xlim[0], side='left')
                g_end = np.searchsorted(ghost_x, view.previous_xlim[1], side='right')
                g_start = max(0, g_start - 50) # Buffer
                g_end = min(len(ghost_x), g_end + 50)
                if g_end > g_start:
                    ghost_x = ghost_x[g_start:g_end]
                    ghost_y = ghost_y[g_start:g_end]
                    if ghost_dy is not None: ghost_dy = ghost_dy[g_start:g_end]
                else:
                    continue # Ghost is completely out of current view

            # 3. Plotting
            color = cmap((i + 1) % 10) # Grab a unique color from tab10
            label = f"Ghost: {ghost_state.name}"

            # Plot Ghost Flux (solid but slightly transparent, low z-order)
            ax.step(ghost_x, ghost_y, where='mid', color=color, alpha=0.7, lw=0.8, label=label, zorder=0.8, rasterized=True)

            # Plot Ghost Error (faint boundary lines of the same color)
            if ghost_dy is not None:
                ax.step(ghost_x, ghost_y - ghost_dy, where='mid', color=color, alpha=0.25, lw=0.4, zorder=0.7, rasterized=True)
                ax.step(ghost_x, ghost_y + ghost_dy, where='mid', color=color, alpha=0.25, lw=0.4, zorder=0.7, rasterized=True)

    def _plot_error(self, ax, arrays: '_PlotArrays', view: '_PlotView'):
        """The error: a band around the flux, or a line of its own."""
        if not self.main_window.error_checkbox.isChecked():
            return
        if arrays.dx is None or arrays.dy is None:
            logging.warning("Skipping error plot because dy_data is None.")
            return
        x_data, y_data, dx_data, dy_data = arrays.x, arrays.y, arrays.dx, arrays.dy
        label = '1-sigma error'

        show_err_as_line = False
        if hasattr(self.main_window, 'error_line_checkbox'):
            show_err_as_line = self.main_window.error_line_checkbox.isChecked()

        if show_err_as_line:
            # OPTION A: Plot Error as separate Line (dy)
            err_color = '#aaaaaa'
            if view.is_scatter:
                ax.plot(x_data, dy_data, label=label,
                        linestyle='None', marker='.', markersize=view.scatter_ms,
                        color=err_color, alpha=0.7, rasterized=True)
            elif view.use_decimation:
                ax.plot(x_data, dy_data, lw=0.8, color=err_color, rasterized=True, label=label)
            else:
                ax.step(x_data, dy_data, where='mid', lw=0.8, color=err_color, rasterized=True, label=label)
        else:
            # OPTION B: Plot Error as Band around Flux (y +/- dy)
            if view.is_scatter:
                # fmt='none' ensures no connecting lines
                ax.errorbar(x_data, y_data, xerr=dx_data,
                            fmt='none', ecolor='#aaaaaa', elinewidth=0.8,
                            alpha=0.5, rasterized=True, zorder=0.9, label='bin size')
                ax.errorbar(x_data, y_data, yerr=dy_data,
                            fmt='none', ecolor='#aaaaaa', elinewidth=0.8,
                            alpha=0.5, rasterized=True, zorder=0.9, label=label)
            elif view.use_decimation:
                # Use plot() for decimated data (fast)
                ax.plot(x_data, y_data - dy_data, lw=0.3, color='#aaaaaa', rasterized=True, label=label)
                ax.plot(x_data, y_data + dy_data, lw=0.3, color='#aaaaaa', rasterized=True)
            else:
                # Use step() for non-decimated data (correct)
                ax.step(x_data, y_data - dy_data, where='mid', lw=0.3, color='#aaaaaa', rasterized=True, label=label)
                ax.step(x_data, y_data + dy_data, where='mid', lw=0.3, color='#aaaaaa', rasterized=True)

    @staticmethod
    def _plot_flux(ax, arrays: '_PlotArrays', view: '_PlotView', colors):
        """The flux (uses Matplotlib style defaults for color)."""
        if view.is_scatter:
            # Use a fast marker plot for echelle data
            # 'ls="None"' ensures no connecting lines (sawtooth removal)
            # 'marker="."' is a pixel-like dot
            ax.plot(arrays.x, arrays.y, label="Spectrum",
                    linestyle='None', marker='.', markersize=view.scatter_ms, markeredgewidth=0,
                    color=colors[0], alpha=0.7, rasterized=True)
        elif view.use_decimation:
            ax.plot(arrays.x, arrays.y, label="Spectrum (decimated)", lw=0.5, color=colors[0], rasterized=True)
        else:
            ax.step(arrays.x, arrays.y, where='mid', label="Spectrum", lw=0.5, color=colors[0], rasterized=True)

    def _plot_continuum_and_model(self, ax, arrays: '_PlotArrays', colors):
        """The continuum (normalization) and the absorption model, when
        shown."""
        if self.main_window.continuum_checkbox.isChecked():
            if arrays.norm is not None:
                ax.plot(
                    arrays.x, arrays.norm, linestyle='--',
                    color='black', lw=0.8, label='Normalization', rasterized=True
                )
            else:
                logging.debug("Normalization requested but spec.norm is None.")

        if self.main_window.model_checkbox.isChecked():
            if arrays.model is not None:
                ax.plot(arrays.x, arrays.model, ls='-', color=colors[1], lw=1.6, alpha=0.8, label='Absorption model', rasterized=True)

    def _plot_aux_column(self, ax, spec, arrays: '_PlotArrays', view: '_PlotView'):
        """The auxiliary column chosen in the main window: the detected
        regions shaded, a mask as a band at the top, or a line."""
        aux_col_to_plot = self.main_window.aux_col_combo.currentText()
        if not aux_col_to_plot or aux_col_to_plot == "None":
            return
        x_data = arrays.x
        try:
            # Get the *full* original aux column data
            full_aux_col_data = spec.get_column(aux_col_to_plot)
            if full_aux_col_data is None:
                raise ValueError(f"Column '{aux_col_to_plot}' returned None.")

            full_aux_col_data = full_aux_col_data.value # Get numpy array

            # As the flux: cut to the view, then decimated
            aux_data = full_aux_col_data[view.final_slice]
            if view.use_decimation:
                if aux_data.dtype.kind == 'b':
                    # A mask: a decimated block is masked if any pixel is
                    aux_data = decimate_y_min_max(aux_data.astype(float),
                                                  DECIMATION_FACTOR) > 0
                else:
                    aux_data = decimate_y_min_max(aux_data, DECIMATION_FACTOR)

            if aux_col_to_plot in ['abs_mask', 'abs_ids']:
                # Create a boolean mask where regions exist
                region_mask = (aux_data > 0)

                # Only plot/label if there are actually regions to show
                if np.any(region_mask):
                    # Use a blended transform to shade the full vertical extent
                    # X in data coordinates, Y in axes coordinates (0 to 1)
                    trans = ax.get_xaxis_transform()
                    ax.fill_between(x_data, 0.05, 0.95, # 0 to 1 in axes coords covers full height
                                    where=region_mask,
                                    color='teal', alpha=0.05,
                                    transform=trans,
                                    label='Detected Regions',
                                    rasterized=True)

            # Plot based on data type
            elif aux_data.dtype.kind in 'fiu': # Float or Int
                if view.is_scatter:
                    ax.plot(x_data, aux_data, label=aux_col_to_plot,
                            linestyle='None', marker=view.scatter_ms, markersize=1,
                            color='purple', rasterized=True)
                elif view.use_decimation:
                    ax.plot(x_data, aux_data, label=aux_col_to_plot,
                            linestyle=':', lw=1.2, color='purple',
                            rasterized=True)
                else:
                    ax.step(x_data, aux_data, where='mid', label=aux_col_to_plot,
                            linestyle=':', lw=1.2, color='purple',
                            rasterized=True)

            elif aux_data.dtype.kind == 'b': # Boolean mask
                trans = ax.get_xaxis_transform()
                ax.fill_between(x_data, 0.95, 1.0,
                                where=aux_data.astype(bool),
                                color='cyan', alpha=0.5,
                                label=f'Mask: {aux_col_to_plot}',
                                transform=trans,
                                step='mid' if not view.use_decimation else None,
                                rasterized=True)

        except Exception as e:
            logging.warning(f"Could not plot auxiliary column '{aux_col_to_plot}': {e}")

    def _plot_systems(self, ax, spec, systs):
        """A tick at the bottom for each transition of each component in
        view (Lyman and metals), the selected ones highlighted; the ticks
        are kept in ``_tick_visuals`` for the tooltip and the menu."""
        if not (self.main_window.systems_checkbox.isChecked() and systs and systs.components):
            return
        # Get current plot limits to only draw visible lines
        xlim = ax.get_xlim()

        # Define Y position for markers (e.g., slightly below top, in axis coords)
        marker_y_axis_coord = 0.05 # 5% from the bottom

        trans = ax.get_xaxis_transform()

        # ** Create lists to batch coordinates **
        hi_lines_x = []
        metal_lines_x = []
        highlight_lines_x = []
        # The reference of a velocity axis: z_em, as in the conversion of
        # the axis by default and in the other conversions (F66: it read
        # '_zem', which the spectrum does not have, so it was always 0)
        zem = spec._data.z_em or 0.0

        # The observed wavelengths of all the transitions first, then
        # converted to the x unit at once: one conversion per
        # transition took most of the redraw (~950 on J0103-1305)
        ticks = []  # (comp, transition, is_highlighted)
        obs_nm = []
        rest_nm = {}
        for comp in systs.components:
            z = comp.z

            # V1/V2 Fallback Logic for Transitions
            transitions = []
            try:
                # 1. Try V2-native lookup (e.g., "CIV")
                transitions = STANDARD_MULTIPLETS[comp.series]
            except KeyError:
                # 2. Try the V1 series syntax (e.g., "CIV_1548,CIV_1550")
                try:
                    transitions = trans_parse(comp.series)
                except Exception:
                    pass # Will be caught by logging below

                # 3. Fallback for single lines (e.g., "Ly_a")
                if not transitions and comp.series in xem_d:
                    transitions = [comp.series]

            if not transitions:
                logging.warning(f"Could not parse series '{comp.series}' for system plot.")

            # Check if this component is highlighted
            is_highlighted = comp.uuid in self._highlighted_uuids

            for t in transitions:
                if t in xem_d:
                    if t not in rest_nm:
                        rest_nm[t] = xem_d[t].to_value(au.nm)
                    ticks.append((comp, t, is_highlighted))
                    obs_nm.append((1 + z) * rest_nm[t])
                else:
                    logging.warning(f"Transition '{t}' for series '{comp.series}' not found in xem_d.")

        try:
            x_ticks = convert_axis_velocity(np.array(obs_nm, dtype=float) * au.nm,
                                            zem, spec.x.unit).value
        except Exception:
            ticks, x_ticks = [], []

        for (comp, t, is_highlighted), x_plot in zip(ticks, x_ticks):
            # Check xlim *before* appending to keep lists small
            if xlim[0] <= x_plot <= xlim[1]:
                # If highlighted, add to special list (also drawn as a
                # standard tick)
                if is_highlighted:
                    highlight_lines_x.append(x_plot)

                if t.startswith('Ly_'):
                    hi_lines_x.append(x_plot)
                else:
                    metal_lines_x.append(x_plot)

                self._tick_visuals.append((ax, x_plot, comp))

        # ** Make only TWO plot calls, outside the loop **
        if hi_lines_x:
            ax.plot(hi_lines_x, [marker_y_axis_coord] * len(hi_lines_x),
                    marker='|', markersize=12, linestyle='None',
                    color='red', alpha=0.9, label="HI Systems",
                    transform=trans)
        if metal_lines_x:
            ax.plot(metal_lines_x, [marker_y_axis_coord] * len(metal_lines_x),
                    marker='|', markersize=12, linestyle='None',
                    color='darkgray', alpha=0.9, label="Metal Systems",
                    transform=trans)

        if highlight_lines_x:
            # 1. Bold Markers at top
            marker, = ax.plot(highlight_lines_x, [marker_y_axis_coord] * len(highlight_lines_x),
                    marker='|', markersize=20, markeredgewidth=2.0, linestyle='None',
                    color='#f5a100', alpha=1.0, label="Selected", zorder=5,
                    transform=trans)

            # 2. Faint Vertical Guide Lines (Gold)
            # These help trace the line down to the spectrum
            guides = [self._highlight_guide(ax, hx) for hx in highlight_lines_x]
            self._highlight_artists = (marker, guides)

    def _plot_cursor(self, session_state, colors):
        """The lines of the redshift cursor, animated (blitted when the
        cursor moves)."""
        if not self.main_window.cursor_show_checkbox.isChecked():
            return
        try:
            z_cursor = float(self.main_window.cursor_z_input.text())
            # Calculate positions using helper
            cursor_x_positions = self.get_cursor_line_positions_at_z(session_state, z_cursor)
            added_cursor_label = False

            for x_plot in cursor_x_positions:
                label = f"Cursor (z={z_cursor:.4f})" if not added_cursor_label else None
                line = self.canvas.axes.axvline(x_plot, ls='--', color=colors[3],
                                                 alpha=0.8, lw=1.0, label=label,
                                                 animated=True) # <<< ANIMATED
                self.canvas.cursor_artists.append(line)
                added_cursor_label = True # Only label first line potentially visible
            logging.debug(f"plot_spectrum: Added {len(self.canvas.cursor_artists)} cursor artists.")
        except (ValueError, AttributeError) as e:
            logging.warning(f"Could not draw cursor lines: {e}")
        except Exception as e:
            logging.error(f"Unexpected error drawing cursor: {e}", exc_info=True)

    # --- Parts of plot_spectrum: legend, axes and limits ---

    def _add_legend(self, ax, legend_loc):
        """The legend, without the entries shown as swatches by the plot
        controls when their section is open."""
        if not ax.has_data(): # Only add legend if something was plotted
            return
        handles, labels = ax.get_legend_handles_labels()
        panel_open = False
        try:
            self.main_window.update_plot_swatches(dict(zip(labels, handles)))
            panel_open = (self.main_window.right_sidebar_widget.width() > 1
                          and self.main_window.plot_elements_section.isExpanded())
        except AttributeError:
            pass  # No main window controls (e.g. in tests)
        if panel_open:
            # The plot controls show these entries as swatches
            # (when their Plot Elements section is open)
            kept = [(h, lab) for h, lab in zip(handles, labels)
                    if lab not in PANEL_LEGEND_LABELS]
            handles = [h for h, _ in kept]
            labels = [lab for _, lab in kept]
        if handles:
            ax.legend(handles, labels, loc=legend_loc, markerscale=0.2) # Smaller legend font

    def _format_axes(self, ax, spec, view: '_PlotView', full_x):
        """Scales, grid, units of the x axis, the text of the coordinates
        under the mouse, the rest-frame axis, the emission lines and the
        labels of the axes."""
        # 1. Scales (Log/Linear)
        ax.set_xscale('log' if view.is_log_x else 'linear')
        ax.set_yscale('log' if view.is_log_y else 'linear')

        ax.grid(True, linestyle=':')

        # 2. X-Axis Formatting (the data stay in nm)
        xaxis_formatter = None  # Store the formatter to re-use it
        if view.x_unit == 'Angstrom':
            xaxis_formatter = plt_ticker.FuncFormatter(lambda x, pos: f"{x * 10:g}")
            ax.xaxis.set_major_formatter(xaxis_formatter)
        elif view.x_unit == 'micron':
            xaxis_formatter = plt_ticker.FuncFormatter(lambda x, pos: f"{x / 1000:g}")
            ax.xaxis.set_major_formatter(xaxis_formatter)
        # else: no formatter needed for nm

        x_label_unit = "nm"
        if view.x_unit == 'Angstrom': x_label_unit = SYM_ANGSTROM
        elif view.x_unit == 'micron': x_label_unit = SYM_MICRON

        z_em = spec._data.z_em
        has_rest_frame = z_em > 0.0 and spec.x.unit.is_equivalent(au.nm)

        ax.format_coord = _coordinate_formatter(spec, view, z_em, x_label_unit)

        # Secondary Rest-Frame Axis (if applicable)
        if has_rest_frame:

            def obs_to_rf(x): return x / (1 + z_em)
            def rf_to_obs(x): return x * (1 + z_em)

            sec_ax = ax.secondary_xaxis('top', functions=(obs_to_rf, rf_to_obs))
            sec_ax.set_xlabel(f"Rest-frame ({SYM_LAMBDA_R}, {x_label_unit})")

            # Apply the unit formatter (or reset to default scalar if None)
            if xaxis_formatter:
                sec_ax.xaxis.set_major_formatter(xaxis_formatter)
            else:
                # Explicitly reset to default if we switched back to 'nm'
                sec_ax.xaxis.set_major_locator(plt_ticker.AutoLocator())

            self._plot_emission_lines(ax, view, z_em, full_x)

        # X/Y-Axis Labels
        ax.set_xlabel(f"Wavelength ({SYM_LAMBDA}, {x_label_unit})")

        if view.is_norm_y:
            y_label_text = f"Normalized Flux ({SYM_FLUX_NORM})"
        elif view.is_snr:
            y_label_text = f"Signal-to-Noise Ratio ({SYM_FLUX} / {view.snr_col})"
        else:
            # Check if we have a real physical unit
            if spec.y.unit is None or spec.y.unit == au.dimensionless_unscaled:
                y_label_text = f"Flux ({SYM_FLUX}, arbitrary units)"
            else:
                # Use Astropy's pretty formatter for the unit string
                unit_str = spec.y.unit.to_string('latex_inline')
                y_label_text = f"Flux ({SYM_FLUX}, {unit_str})"

        ax.set_ylabel(y_label_text)

    def _plot_emission_lines(self, ax, view: '_PlotView', z_em, full_x):
        """The emission lines at z_em in view: a label at the top and a
        faint vertical line (the strong ones darker, or only them)."""
        # The view: the full data range if autoscaling, the zoomed range
        # if kept
        if view.force_autoscale:
            view_xmin, view_xmax = full_x[0], full_x[-1]
        elif view.was_zoomed:
            view_xmin, view_xmax = view.previous_xlim
        else:
            view_xmin, view_xmax = full_x[0], full_x[-1]

        # Transformation: X is data coords, Y is axes coords (0=bottom, 1=top)
        trans = ax.get_xaxis_transform()

        strong_only = self.main_window.strong_lines_checkbox.isChecked()

        # Loop through standard emission lines
        for label, rest_wave_q in xem_d.items():
            if strong_only and label not in STRONG_EMISSION_LINES:
                continue

            if label not in STRONG_EMISSION_LINES:
                color = 'gray'
                alpha_text = 0.5
                alpha_line = 0.15
            else:
                color = 'black'
                alpha_text = 1.0
                alpha_line = 0.3

            # Calculate observed wavelength in native units (nm)
            try:
                rest_nm = rest_wave_q.to_value(au.nm)
                obs_nm = rest_nm * (1.0 + z_em)

                # Only plot if within the current view
                if view_xmin <= obs_nm <= view_xmax:
                    ax.text(obs_nm, 0.95, label,
                            transform=trans,
                            rotation=90,
                            verticalalignment='top',
                            horizontalalignment='center',
                            color=color, alpha=alpha_text,
                            fontsize=8, clip_on=True)
                    ax.axvline(obs_nm, color=color, alpha=alpha_line, lw=1.0, ls=':')

            except Exception:
                continue # Skip lines with incompatible units (rare)

    @staticmethod
    def _set_limits(ax, view: '_PlotView'):
        """The zoom of the previous plot kept, the default normalized range
        when the normalization was switched on, otherwise autoscale.

        With the logarithmic F scale a lower limit that is not positive
        (e.g. the zoom of a linear plot) is not kept: the default lower
        limit in the normalized view, autoscale otherwise (F78)."""
        # --- X Limits ---
        # Priority: Autoscale if requested (e.g., by unit change)
        if view.was_zoomed and not view.force_autoscale:
            logging.debug(f"Restoring previous X zoom: {view.previous_xlim}")
            ax.set_xlim(view.previous_xlim)
        # Else (initial draw or zoomed out), let Matplotlib autoscale
        else:
            logging.debug("Allowing default X autoscale.")
            # No command needed, autoscale is default after clear()

        # --- Y Limits ---
        if view.norm_state_changed:
            if view.is_norm_y:
                logging.debug("Normalization toggled ON: Setting default Y limits.")
                ax.set_ylim(NORM_YLIM_LOG if view.is_log_y else NORM_YLIM)
            else:
                logging.debug("Normalization toggled OFF: Autoscaling Y.")
                ax.autoscale(enable=True, axis='y')
        elif view.was_zoomed and not view.force_autoscale:
            ylim = view.previous_ylim
            if view.is_log_y and ylim[0] <= 0:
                if view.is_norm_y and ylim[1] > NORM_YLIM_LOG[0]:
                    ax.set_ylim(NORM_YLIM_LOG[0], ylim[1])
                else:
                    ax.autoscale(enable=True, axis='y')
            else:
                logging.debug(f"Restoring previous Y zoom: {ylim}")
                ax.set_ylim(ylim)
        # Else (initial draw or zoomed out), let Matplotlib autoscale
        else:
            logging.debug("Allowing default Y autoscale.")

    def _draw_isomag_grid(self):
        """
        Helper to draw the iso-magnitude grid.
        Crucially, it saves and restores axis limits so these lines
        do NOT affect autoscaling or the 'Home' view.
        """
        # 1. Clear old artists first
        for artist in self.isomag_artists:
            try: artist.remove()
            except: pass
        self.isomag_artists = []

        if not self._show_isomag_grid:
            return

        # 2. Get data and perform checks
        if not self.main_window.active_history: return
        spec = self.main_window.active_history.current_state.spec
        if not spec: return

        standard_flux = au.erg / (au.cm**2 * au.s * au.Angstrom)
        if spec.y.unit is None or not spec.y.unit.is_equivalent(standard_flux):
             return

        # 3. SAVE CURRENT LIMITS (The secret to acting like a passive grid)
        current_xlim = self.canvas.axes.get_xlim()
        current_ylim = self.canvas.axes.get_ylim()

        # 4. Calculate and Plot
        try:
            x_aa = spec.x.to(au.Angstrom).value
            x_plot = spec.x.value
            current_y_unit = spec.y.unit
            
            c_aa_s = const.c.to(au.Angstrom/au.s).value
            jy_to_cgs_fnu = (1.0 * au.Jy).to(au.erg / (au.cm**2 * au.s * au.Hz)).value
            
            grid_mags = [15, 16, 17, 18, 19, 20, 21]
            
            x_label_pos = min(current_xlim[1], x_plot[-1])

            for mag in grid_mags:
                # Calculate curve
                f_nu_val = 3631.0 * (10**(-0.4 * mag)) * jy_to_cgs_fnu
                f_lambda_val = f_nu_val * c_aa_s / (x_aa**2)
                y_plot_iso = (f_lambda_val * standard_flux).to(current_y_unit).value

                # Plot curve
                artist, = self.canvas.axes.plot(x_plot, y_plot_iso, 
                                                ls=':', lw=1.0, color='orange', alpha=0.4,
                                                zorder=0.1) 
                self.isomag_artists.append(artist)
                
                # --- NEW: Smart Labeling ---
                # Interpolate to find Y at our chosen X position
                # We can safely use the last data point as 'right' fill value because
                # we ensured x_label_pos <= x_plot[-1]
                y_label_pos = np.interp(x_label_pos, x_plot, y_plot_iso)
                
                # Only draw label if it's vertically within the current view
                if current_ylim[0] <= y_label_pos <= current_ylim[1]:
                    # Use annotate to place text exactly at (x_label_pos, y_label_pos)
                    # and then offset it slightly to the right (10 points) so it doesn't overlap.
                    lbl = self.canvas.axes.annotate(
                        f"m={mag}", 
                        xy=(x_label_pos, y_label_pos), 
                        xytext=(5, 0),                 # 5 points to the right
                        textcoords="offset points", 
                        color='orange', fontsize=9, alpha=0.8,
                        verticalalignment='center',
                        clip_on=False # Allow it to spill into the margin if needed
                    )
                    self.isomag_artists.append(lbl)
                # ----------------------------------------------------

        except Exception as e:
            logging.error(f"Error drawing iso-mag grid: {e}")

        # 5. RESTORE LIMITS
        self.canvas.axes.set_xlim(current_xlim)
        self.canvas.axes.set_ylim(current_ylim)   

    def toggle_isomag_grid(self, visible: bool):
        """ Toggles the state and triggers a redraw. """
        self._show_isomag_grid = visible
        # We call plot_spectrum to ensure everything (including limits) 
        # is in a consistent state.
        self.plot_spectrum(
            session_state=self.main_window.active_history.current_state,
            force_autoscale=False 
        )
    
    def update_plot(self, new_session_state: Optional['SessionV2'], force_autoscale: bool = False):
        """Called by MainWindowV2 to swap the immutable session object and redraw."""
        self.plot_spectrum(session_state=new_session_state, initial_draw=False, force_autoscale=force_autoscale)

    def set_highlights(self, uuids: set):
        """Updates the list of highlighted component UUIDs and refreshes the plot."""
        if self._highlighted_uuids != uuids:
            self._highlighted_uuids = uuids
            
            # Trigger a redraw using the current session state
            session_state = None
            if self.main_window and self.main_window.active_history:
                session_state = self.main_window.active_history.current_state

            # Only the highlights changed: drawn without rebuilding the plot
            if (session_state is not None and session_state is self._plotted_state
                    and self._blit_highlights()):
                return
            self.plot_spectrum(session_state, force_autoscale=False)

    @staticmethod
    def _highlight_guide(ax, x):
        """Faint vertical guide line of a highlighted system."""
        return ax.axvline(x, color='#f5a100', alpha=0.6, lw=1.0, linestyle=':', zorder=4)

    def _on_full_draw(self, event):
        if not self._drawing_highlight_base:
            self._highlight_base = None

    def _blit_highlights(self) -> bool:
        """Draw the highlights of ``self._highlighted_uuids`` over the plot,
        without rebuilding it and without a full draw (phase 5: a full draw
        of the main plot took ~0.15 s on a laptop screen, at each selection
        in the System Inspector).

        The pixels of the axes without the highlights and the legend are
        kept from a full draw (``_highlight_base``, made at the first
        selection after each full draw); a new selection restores them and
        draws the guide lines, the markers and the legend over them, in the
        order of a full draw (zorder 4, 5, then the legend): the same
        pixels as a full draw. Returns False when a full redraw is needed:
        the "Selected" entry of the legend appears or disappears, or the
        continuum is being edited.
        """
        if getattr(self, '_edit_mode_active', False):
            return False
        new_x = [x for _, x, comp in self._tick_visuals if comp.uuid in self._highlighted_uuids]
        if bool(new_x) != (self._highlight_artists is not None):
            return False
        if not new_x:
            return True  # nothing highlighted, before and after

        canvas = self.canvas
        ax = canvas.axes
        marker, guides = self._highlight_artists
        marker.set_data(new_x, [marker.get_ydata()[0]] * len(new_x))
        for guide in guides:
            guide.remove()
        guides = [self._highlight_guide(ax, x) for x in new_x]
        self._highlight_artists = (marker, guides)
        overlay = guides + [marker]
        if ax.get_legend() is not None:
            overlay.append(ax.get_legend())

        if self._highlight_base is None:
            for artist in overlay:
                artist.set_visible(False)
            self._drawing_highlight_base = True
            try:
                canvas.draw()
            finally:
                self._drawing_highlight_base = False
                for artist in overlay:
                    artist.set_visible(True)
            self._highlight_base = canvas.copy_from_bbox(ax.bbox)
        else:
            canvas.restore_region(self._highlight_base)
        for artist in overlay:
            ax.draw_artist(artist)
        # The background of the blitted cursor, which is drawn over it
        canvas.background = canvas.copy_from_bbox(ax.bbox)
        for artist in canvas.cursor_artists:
            ax.draw_artist(artist)
        canvas.blit(ax.bbox)
        return True
