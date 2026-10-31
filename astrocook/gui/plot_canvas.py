"""The canvas of the main plot: blitting of the cursor of z, the tooltip
after a pause of the mouse, the span of a selected region; the clicks are
passed to the plot widget.

Moved out of ``pyside_plot.py`` in phase 5 (split of the GUI).
"""
import logging
import matplotlib
import matplotlib.style
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import numpy as np
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QToolTip

from .plot_style import SELECTION_STYLE


def _snap_layout_to_pixels(fig):
    """Round the positions of the axes laid out by the constrained layout
    to whole pixels. The layout gives positions that differ by ~1e-13 px
    from one draw to the next: that changes how the long path of the
    spectrum is simplified, and the pixels kept for blitting (highlights,
    cursor) no longer matched a new full draw of the same plot."""
    engine = fig.get_layout_engine()
    execute = engine.execute

    def execute_and_snap(figure):
        result = execute(figure)
        w, h = figure.bbox.width, figure.bbox.height
        if w > 0 and h > 0:
            for ax in figure.axes:
                b = ax.get_position(original=False)
                x0, x1 = round(b.x0 * w) / w, round(b.x1 * w) / w
                y0, y1 = round(b.y0 * h) / h, round(b.y1 * h) / h
                ax._set_position([x0, y0, x1 - x0, y1 - y0], which='both')
        return result

    engine.execute = execute_and_snap


class MatplotlibCanvas(FigureCanvasQTAgg):
    """A custom Matplotlib canvas enabling blitting for cursor dragging."""
    def __init__(self, parent=None, width=5, height=4, dpi=100, plot_widget=None): # Added plot_widget ref
        try:
            #plt.style.use(['science', 'fast'])
            matplotlib.style.use('fast')
        except Exception as e:
            logging.warning(f"Could not apply 'scienceplots': {e}. Using default.")
            matplotlib.style.use('fast')

        fig = Figure(figsize=(width, height), dpi=dpi, constrained_layout=True)
        self.axes = fig.add_subplot(111)
        _snap_layout_to_pixels(fig)
        #fig.tight_layout()

        super(MatplotlibCanvas, self).__init__(fig)

        # We rely on the QtAgg backend's native speed for now.
        if parent:
            self.setParent(parent)

        self.plot_widget = plot_widget # Store reference

        # Blitting Attributes
        self.background = None # Store the clean background pixels
        self.cursor_artists = [] # Store references to the cursor axvline objects
        # The redshift of the cursor, written at the top of its first line
        # in view: the box in the card of the plot controls fades while the
        # mouse is on the plot (findings register F92)
        self.cursor_label = None
        self.draw_event_cid = None # Store connection ID for draw_event

        # Hover Timer for Tooltips
        self.hover_timer = QTimer(self)
        self.hover_timer.setSingleShot(True)
        self.hover_timer.setInterval(500) # 500ms delay before showing tooltip
        self.hover_timer.timeout.connect(self._on_hover_timeout)
        self._last_mouse_event = None

        # Region Selection State
        self.selection_artist = None # Will hold the axvspan
        self.selection_start_x = None

        # Connect Matplotlib events
        self.mpl_connect('motion_notify_event', self.on_motion)
        self.mpl_connect('button_press_event', self.on_press)
        self.mpl_connect('button_release_event', self.on_release)
        self.mpl_connect('resize_event', self._on_resize)
        self.xlim_cid = self.axes.callbacks.connect('xlim_changed', self.on_lim_changed)
        self.ylim_cid = self.axes.callbacks.connect('ylim_changed', self.on_lim_changed)


    def on_lim_changed(self, axes):
        """Callback: Invalidate background and schedule a full redraw."""
        if self.background is not None and self.plot_widget: # Check for plot_widget
            logging.debug(f"Limits changed ({axes.get_label()}): Invalidating background & scheduling plot_spectrum.")
            self.background = None # Invalidate background cache

            # Get the current state to pass to the scheduled plot call
            session_state = None
            if self.plot_widget.main_window and self.plot_widget.main_window.active_history:
                session_state = self.plot_widget.main_window.active_history.current_state

            # Schedule plot_spectrum to run soon via the event loop
            # ** Use lambda to pass the session_state argument **
            QTimer.singleShot(0, lambda: self.plot_widget.plot_spectrum(session_state=session_state))

    def blit_box(self):
        """The pixels kept and redrawn when the lines over the plot move
        (cursor, highlights): the axes and 2 pixels around them. The edges
        of the axes are at fractions of a pixel, and the vertical lines end
        there: a row of pixels just outside the axes was drawn over at each
        move and never restored."""
        return self.axes.bbox.padded(2)

    def update_cursor_label(self, z, x_positions):
        """Move the label of the cursor to its first line in view, and write
        its redshift; hidden when no line is in view. On the left of the
        line in the right half of the plot, where the card of the plot
        controls would cover it."""
        if self.cursor_label is None:
            return
        lo, hi = sorted(self.axes.get_xlim())
        in_view = [x for x in x_positions if lo <= x <= hi]
        self.cursor_label.set_visible(bool(in_view))
        if in_view:
            x = in_view[0]
            right_half = (x - lo) > 0.5 * (hi - lo)
            self.cursor_label.set_x(x)
            self.cursor_label.set_horizontalalignment('right' if right_half else 'left')
            self.cursor_label.set_text(f"z = {z:.5f} " if right_half else f" z = {z:.5f}")

    def _on_resize(self, event):
        """The pixels kept for blitting have the old size: captured again at
        the redraw that follows the resize. Before, the redshift cursor
        restored them over the resized plot, e.g. after View > Show Toolbar
        (reported by Guido in phase 5-bis)."""
        self.background = None
        if self.plot_widget is not None and hasattr(self.plot_widget, '_connect_background_capture'):
            self.plot_widget._connect_background_capture()

    def _capture_background(self, event):
        """Callback for draw_event to capture the background for blitting."""
        try:
            if self.background is None and self.figure.get_size_inches().any():
                if self.axes.bbox.width > 0 and self.axes.bbox.height > 0:
                    self.background = self.copy_from_bbox(self.blit_box())
        except Exception as e:
            logging.error(f"Failed during background copy: {e}")
            self.background = None
        finally:
            if self.draw_event_cid is not None:
                try: self.mpl_disconnect(self.draw_event_cid)
                except Exception: pass
                self.draw_event_cid = None
                 
    def draw_animated(self):
        if self.background is None:
            # This can happen in the small gap between on_lim_changed
            # and the full redraw completing. Calling draw_idle()
            # is a safe, lightweight way to request the redraw.
            self.draw_idle()
            return

        # Proceed with blitting
        self.restore_region(self.background)
        drawn_artists = []
        for artist in self.cursor_artists:
            try: self.axes.draw_artist(artist); drawn_artists.append(artist)
            except Exception as e: logging.error(f"Error drawing anim artist {artist}: {e}")
        if getattr(self.plot_widget, '_edit_mode_active', False):
            if hasattr(self.plot_widget, 'knot_artist') and self.plot_widget.knot_artist:
                self.axes.draw_artist(self.plot_widget.knot_artist)
            if hasattr(self.plot_widget, 'spline_artist') and self.plot_widget.spline_artist:
                self.axes.draw_artist(self.plot_widget.spline_artist)
        if self.cursor_label is not None and drawn_artists:
            self.axes.draw_artist(self.cursor_label)
        if drawn_artists:
            try: self.blit(self.blit_box())
            except Exception as e: logging.error(f"Blitting failed: {e}"); self.draw_idle()


    def on_motion(self, event):
        """Handles mouse motion: Updates cursor position if active."""
        # 1. Always restart the hover timer on any motion
        # Immediately hide any existing tooltip when the mouse moves.
        # It will reappear only when the user stops again (timer timeout).
        QToolTip.hideText()
        
        self._last_mouse_event = event
        self.hover_timer.start()

        # Safety checks
        if (not self.plot_widget or not hasattr(self.plot_widget, 'main_window')
            or not self.plot_widget.main_window):
            # logging.debug("on_motion: plot_widget or main_window invalid.")
            return

        main_window = self.plot_widget.main_window # type: MainWindowV2

        # Update Selection Span while dragging
        if (self.plot_widget and self.plot_widget._is_selecting_region and
            self.selection_start_x is not None and event.inaxes == self.axes):
            
            # Update the span's extent
            # axvspan is a Polygon, we convert it to a simple Rect for easy updating if needed,
            # but .set_xy is the standard way to update a Polygon.
            # Actually, for axvspan it's easier to just update the vertices.
            # A simpler way for axvspan: it works by setting xmin/xmax.
            # BUT axvspan doesn't have simple set_xlim methods.
            # Re-drawing might be easiest, or using a Rectangle instead.
            
            # Let's use a simpler approach: remove old, draw new (fast enough for this)
            if self.selection_artist:
                self.selection_artist.remove()
            
            self.selection_artist = self.axes.axvspan(
                self.selection_start_x, event.xdata, **SELECTION_STYLE
            )
            self.draw_idle()
            return # Don't do cursor updates while selecting
        
        # Check if inside axes and cursor checkbox is checked
        if (event.inaxes == self.axes and event.xdata is not None
            and main_window.cursor_show_checkbox.isChecked()):

            session_state = None
            if main_window.active_history:
                session_state = main_window.active_history.current_state

            if not session_state or not session_state.spec:
                logging.debug(
                    "on_motion: No active session state or spec. Skipping.")
                return  # Can't do conversions without a spec
            
            # --- Always update on motion if cursor is active ---
            spec = session_state.spec
            if not spec: return

            try:
                new_z = self.plot_widget.calculate_z_from_x(event.xdata)

                if new_z is not None:
                    # Update QLineEdit (no need for blockSignals if not dragging)
                    current_text_z = None
                    try: current_text_z = float(main_window.cursor_z_input.text())
                    except ValueError: pass
                    # Update text only if significantly different to avoid excessive signals
                    if current_text_z is None or not np.isclose(new_z, current_text_z, atol=1e-7):
                        main_window.cursor_z_input.setText(f"{new_z:.7f}")

                    # Update Artist Positions & Redraw using Blit
                    new_x_positions = self.plot_widget.get_cursor_line_positions_at_z(session_state, new_z)

                    # Ensure artists exist and match count
                    if not self.cursor_artists or len(new_x_positions) != len(self.cursor_artists):
                        # Force a full redraw if artists are missing/mismatched
                        logging.warning("Cursor artists missing or mismatched on motion. Forcing full redraw.")
                        self.plot_widget.plot_spectrum(session_state=session_state) # Recreate artists
                    else:
                        # Update existing artists
                        for artist, x_pos in zip(self.cursor_artists, new_x_positions):
                            artist.set_xdata([x_pos])
                        self.update_cursor_label(new_z, new_x_positions)
                        self.draw_animated() # Use blitting redraw

            except Exception as e:
                logging.error(f"Error during cursor hover update: {e}", exc_info=True)
        # else:
            # Optional: Clear coordinate display if outside axes or cursor off?
            # Toolbar should handle clearing coordinate display by default.
            # pass
        self.plot_widget.on_mouse_move(event)


    def _on_hover_timeout(self):
        """Called when the mouse has stopped moving for a while."""
        if self.plot_widget and self._last_mouse_event and self._last_mouse_event.inaxes == self.axes:
            self.plot_widget.show_data_tooltip(self._last_mouse_event)

    def on_press(self, event):
        """Delegate press to widget."""
        # The plot takes the focus: the keys of its toolbar act on it, and
        # the shortcuts of the window still work (they act on the whole
        # window). Before, the focus went to the window itself
        self.setFocus()

        if self.plot_widget:
            self.plot_widget.on_press(event)

    def on_release(self, event):
        """Delegate release to widget."""
        if self.plot_widget:
            self.plot_widget.on_release(event)

    def cleanup_selection(self):
        """Helper to remove selection artist and reset state."""
        self.selection_start_x = None
        if self.selection_artist:
            try:
                self.selection_artist.remove()
            except Exception: pass # Might already be gone
            self.selection_artist = None
            self.draw_idle()
