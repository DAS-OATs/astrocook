"""The toolbar of the plots (main plot, velocity plot, corner plot): the
matplotlib toolbar without the buttons Astrocook does not use, and a Home
that replots the current state.

Moved out of ``pyside_plot.py`` in phase 5 (split of the GUI). Redesigned in
phase 5-bis (Guido's request): icons of one family (Material Design 6, in
the colour of the text, also in dark mode), rounded buttons, the active tool
on the accent colour, keyboard shortcuts. The functions of the Matplotlib
toolbar are all kept.
"""
import logging
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from PySide6.QtCore import QEvent, QSize, Qt
from PySide6.QtGui import QAction, QKeySequence, QPalette, QShortcut
from PySide6.QtWidgets import QApplication, QSizePolicy
import qtawesome as qta

from .floating_cards import CardFader
from .theme import ACCENT

#: The icon of each tool (Material Design 6, outlined)
#: The distance (px) of the floating toolbar from the corner of its plot
FLOAT_MARGIN = 8

TOOL_ICONS = {
    'home': 'mdi6.home-outline',
    'back': 'mdi6.arrow-left',
    'forward': 'mdi6.arrow-right',
    'pan': 'mdi6.cursor-move',
    'zoom': 'mdi6.magnify-plus-outline',
    'select': 'mdi6.arrow-expand-horizontal',  # a range along the x axis
    'save_figure': 'mdi6.content-save-outline',
}
#: The keyboard shortcuts, active while the plot has the focus (a click on
#: it): Matplotlib's own keys work only in its windows, not in a plot
#: embedded in an application. S is "save" in Matplotlib, so Select is X;
#: Ctrl+S is left to the session
TOOL_SHORTCUTS = {
    'home': 'H', 'back': 'Left', 'forward': 'Right',
    'pan': 'P', 'zoom': 'Z', 'select': 'X',
}
TOOL_TIPS = {
    'home': 'Home: the whole spectrum',
    'back': 'Back to the previous view',
    'forward': 'Forward to the next view',
    'pan': 'Pan: drag to move, right-drag to zoom',
    'zoom': 'Zoom: draw a rectangle',
    'select': 'Select a region: drag along the wavelength axis, then choose what '
              'to do with it',
    'save_figure': 'Save the figure',
}
ICON_SIZE = QSize(24, 24)


class AstrocookToolbar(NavigationToolbar):
    """
    Custom Matplotlib toolbar that inherits from NavigationToolbar2QT 
    to remove unwanted buttons and simplify the interface.
    """
    # Define the order of toolbar items to keep (based on V1 functionality)
    toolitems = [
        # Label, Tooltip, Icon Filename, Tool ID
        ('Home', 'Reset original view', 'home', 'home'),
        ('Back', 'Back to previous view', 'back', 'back'),
        ('Forward', 'Forward to next view', 'forward', 'forward'),
        (None, None, None, None), # Separator placeholder
        
        # FIX: Use Matplotlib's internal symbolic names for icons (third element).
        # We rely on the base class finding 'hand' and 'magnifying_glass' equivalent icons.
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'), 
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
        
        (None, None, None, None), # Separator placeholder
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
    ]

    def __init__(self, canvas, parent, select: bool = False):
        """``select``: add the Select tool (selection of a region, main plot
        only, F73)."""
        # 1. Base Class Initialization: This CREATES and stores the C++ QAction objects.
        # This MUST run first.
        super().__init__(canvas, parent, coordinates=True) 
        self.viewer_parent = parent
        self.select_action = None

        # 2. Final cleanup (Remove unwanted actions created by the base class)
        actions_to_remove = ['Customize', 'Subplots', 'Edit axis, curve and image parameters']
        
        for action in self.actions():
            if action.text() in actions_to_remove or action.iconText() in actions_to_remove:
                self.removeAction(action)

        # 3. The Select tool, after Zoom, as a third mode exclusive with
        # Pan and Zoom
        if select:
            self.select_action = QAction('Select', self)
            self.select_action.setCheckable(True)
            self.select_action.triggered.connect(self.toggle_select_mode)
            actions = self.actions()
            after_zoom = actions[actions.index(self._actions['zoom']) + 1]
            self.insertAction(after_zoom, self.select_action)

        # 4. Look and shortcuts
        self.setIconSize(ICON_SIZE)
        self._tool_actions = dict(self._actions)
        if self.select_action is not None:
            self._tool_actions['select'] = self.select_action
        for key, action in self._tool_actions.items():
            tip = TOOL_TIPS[key]
            if key in TOOL_SHORTCUTS:
                tip += f" ({QKeySequence(TOOL_SHORTCUTS[key]).toString(QKeySequence.NativeText)})"
            action.setToolTip(tip)
        self.refresh_icons()
        self._add_shortcuts(parent)

    def set_message_sink(self, sink):
        """Send the text of the toolbar (the coordinates under the mouse) to
        ``sink`` (a status bar) instead of showing it: the toolbar keeps
        only its tools, and the coordinates stay visible when it is hidden
        (phase 5-bis)."""
        self._message_sink = sink
        # The label is in an action of the toolbar, which shows it again:
        # the action is hidden
        for action in self.actions():
            if self.widgetForAction(action) is getattr(self, 'locLabel', None):
                action.setVisible(False)

    def set_message(self, s):
        sink = getattr(self, '_message_sink', None)
        if sink is not None:
            sink(s)
        else:
            super().set_message(s)

    def _add_shortcuts(self, parent):
        """The keys of the tools, on the plot widget and its children: they
        act while the plot has the focus (a click on it), and do not take
        the keys of the fields or tables of the window. Esc switches the
        active tool off."""
        if parent is None:
            return
        # The canvas takes the focus when clicked
        self.canvas.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        slots = {'home': self.home, 'back': self.back, 'forward': self.forward,
                 'pan': self.pan, 'zoom': self.zoom}
        if self.select_action is not None:
            slots['select'] = self.toggle_select_mode
        self._shortcuts = []
        for key, slot in slots.items():
            shortcut = QShortcut(QKeySequence(TOOL_SHORTCUTS[key]), parent)
            shortcut.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            shortcut.activated.connect(slot)
            self._shortcuts.append(shortcut)
        escape = QShortcut(QKeySequence(Qt.Key.Key_Escape), parent)
        escape.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        escape.activated.connect(self._all_tools_off)
        self._shortcuts.append(escape)

    def _all_tools_off(self):
        self._leave_select()
        self.deactivate_tools()

    def refresh_icons(self):
        """The icons in the colour of the text (black, or light in dark
        mode), white on the accent of the active tool, grey when disabled
        (Back and Forward with no view to go to); again when the palette
        switches between light and dark (called by the window)."""
        # The palette of the application: when it changes, the window is
        # told before its children get the new one
        palette = QApplication.palette()
        text = palette.color(QPalette.ColorRole.WindowText)
        disabled = palette.color(QPalette.ColorGroup.Disabled, QPalette.ColorRole.WindowText)
        for key, action in self._tool_actions.items():
            action.setIcon(qta.icon(TOOL_ICONS[key], color=text, color_on='white',
                                    color_disabled=disabled, scale_factor=1.15))
        floating = getattr(self, '_anchor', None) is not None
        # Floating over the plot: a rounded card, as the tabs of the cards
        frame = (f"background: {palette.color(QPalette.ColorRole.Window).name()}; "
                 f"border: 1px solid {palette.color(QPalette.ColorRole.Mid).name()}; "
                 "border-radius: 10px;") if floating else ""
        self.setStyleSheet(f"""
            QToolBar {{ spacing: 2px; padding: 4px 6px; {frame} }}
            QToolButton {{ border: none; border-radius: 6px; padding: 5px; }}
            QToolButton:hover {{ background: rgba(128, 128, 128, 0.18); }}
            QToolButton:checked {{ background: {ACCENT}; }}
            QLabel {{ color: {palette.color(QPalette.ColorRole.PlaceholderText).name()}; }}
        """)

    def float_over(self, anchor):
        """Float the toolbar over the bottom left corner of ``anchor`` (the
        plot), out of the layout: the plot takes the room it had, and the
        toolbar fades while the mouse is away from it, as the cards of the
        main window (phase 5-bis, Guido's request). Its keys still work."""
        self._anchor = anchor
        self.setParent(anchor.parentWidget())
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._fader = CardFader(self)
        anchor.installEventFilter(self)
        self.refresh_icons()
        self.place()

    def place(self):
        """At the bottom left corner of the plot, above it."""
        anchor = getattr(self, '_anchor', None)
        if anchor is None:
            return
        self.adjustSize()
        geometry = anchor.geometry()
        self.move(geometry.left() + FLOAT_MARGIN,
                  geometry.bottom() - self.height() - FLOAT_MARGIN + 1)
        self.raise_()

    def eventFilter(self, obj, event):
        if obj is getattr(self, '_anchor', None) and event.type() in (
                QEvent.Type.Resize, QEvent.Type.Move, QEvent.Type.Show):
            self.place()
        return super().eventFilter(obj, event)

    def deactivate_tools(self):
        """Switch off the active tool (Pan or Zoom), if any.

        Matplotlib's ``pan`` and ``zoom`` toggle their tool and ignore their
        argument: ``pan(False)`` with no tool on switches pan on (findings
        register F81).
        """
        if self.mode == 'pan/zoom':
            self.pan()
        elif self.mode == 'zoom rect':
            self.zoom()

    def toggle_select_mode(self):
        """Toggles the region selection mode on/off (the Select button)."""
        if not self.viewer_parent: return
        self.viewer_parent.toggle_region_selector()

    def sync_select_action(self):
        """Check the Select button as the selection mode is."""
        if self.select_action is not None:
            self.select_action.setChecked(bool(self.viewer_parent._is_selecting_region))

    def set_select_enabled(self, enabled: bool):
        """The Select tool is offered only on a wavelength axis (the region
        becomes an expression on λ): switched off otherwise."""
        if self.select_action is None:
            return
        if not enabled and self.viewer_parent._is_selecting_region:
            self.viewer_parent.set_region_selection(False)
        self.select_action.setEnabled(enabled)

    def _leave_select(self):
        if self.select_action is not None and self.viewer_parent._is_selecting_region:
            self.viewer_parent.set_region_selection(False)

    def pan(self, *args):
        """Pan, switching the Select tool off (the three tools exclude each
        other)."""
        self._leave_select()
        super().pan(*args)

    def zoom(self, *args):
        """Zoom, switching the Select tool off."""
        self._leave_select()
        super().zoom(*args)

    def home(self):
        """ Overrides the default 'home' to reset to the *current*
            data state, respecting normalization.
        """
        logging.debug("Custom 'Home' button pressed.")
        if self.viewer_parent and self.viewer_parent.main_window:
            session_state = None
            if self.viewer_parent.main_window.active_history:
                session_state = self.viewer_parent.main_window.active_history.current_state
            
            # Call our plot function with the new flag
            self.viewer_parent.plot_spectrum(
                session_state=session_state,
                force_autoscale=True
            )
        else:
            # Fallback to original behavior
            super().home()
