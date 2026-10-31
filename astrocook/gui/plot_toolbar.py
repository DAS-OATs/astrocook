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
from PySide6.QtWidgets import QApplication, QPushButton, QSizePolicy

QWIDGETSIZE_MAX = (1 << 24) - 1  # the maximum size of a widget in Qt
import qtawesome as qta

from .floating_cards import CARD_WIDTH, CardFader
from .theme import ACCENT

#: The icon of each tool (Material Design 6, outlined)
#: The tab that shows or hides the floating toolbar, as those of the cards
TAB_SIZE = 28
TAB_ICON = 'mdi6.tools'
#: The icons of the floating toolbar, a little smaller than in a toolbar
#: row, to fit the width of the cards
FLOAT_ICON_SIZE = 22

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
        # Floating over the plot: the background of the side cards, no
        # border, rounded on the inner corner only (it is flush with the
        # bottom left corner of the plot area)
        card = palette.color(QPalette.ColorRole.Window)
        background = f"rgba({card.red()}, {card.green()}, {card.blue()}, 0.95)"
        frame = (f"background: {background}; border: none; border-radius: 0px; "
                 "border-top-right-radius: 8px;") if floating else ""
        if floating:
            self.tab.setIcon(qta.icon(TAB_ICON, color=text))
            self.tab.setIconSize(QSize(18, 18))
            self.tab.setStyleSheet(
                f"QPushButton#ToolbarTab {{ background: {background}; border: none; "
                "border-top-right-radius: 6px; border-bottom-right-radius: 0px; "
                "margin: 0px; padding: 0px; }")
        # Tighter when floating, to fit the width of the cards
        spacing, padding, button = ("1px", "4px 4px", "3px") if floating else ("2px", "4px 6px", "5px")
        self.setStyleSheet(f"""
            QToolBar {{ spacing: {spacing}; padding: {padding}; {frame} }}
            QToolButton {{ border: none; border-radius: 6px; padding: {button}; }}
            QToolButton:hover {{ background: rgba(128, 128, 128, 0.18); }}
            QToolButton:checked {{ background: {ACCENT}; }}
            QLabel {{ color: {palette.color(QPalette.ColorRole.PlaceholderText).name()}; }}
        """)

    def float_over(self, anchor):
        """Float the toolbar over the plot, as the side cards of the main
        window (phase 5-bis, Guido's request): out of the layout, in the
        bottom left corner of the area of ``anchor`` (the plot), flush with
        it, without a border, fading while the mouse is away; a tab at its
        right edge shows or hides it. Its keys work also when hidden."""
        self._anchor = anchor
        host = anchor.parentWidget()
        self.setParent(host)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.tab = QPushButton(host)
        self.tab.setObjectName("ToolbarTab")
        self.tab.setFixedSize(TAB_SIZE, TAB_SIZE)
        self.tab.setCursor(Qt.CursorShape.PointingHandCursor)
        self.tab.clicked.connect(lambda: self.set_shown(not self.isVisible(), notify=True))
        self.on_toggled = None  # called with the new state when the tab is clicked
        self._fader = CardFader(self)
        self._fader.add_handle(self.tab)
        anchor.installEventFilter(self)
        host.installEventFilter(self)
        self.refresh_icons()
        self.place()

    def _fit_to_width(self, width: int):
        size = FLOAT_ICON_SIZE
        self.ensurePolished()
        self.setMinimumWidth(0)
        self.setMaximumWidth(QWIDGETSIZE_MAX)
        while True:
            if self.iconSize().width() != size:
                self.setIconSize(QSize(size, size))
            if self.sizeHint().width() <= width or size <= 16:
                break
            size -= 2
        # The toolbar of the main plot as wide as the side cards next to it;
        # that of the System Inspector (no Select, no cards) as wide as its
        # tools, with the same icons, not to end with an empty space
        self.setFixedWidth(width if self.select_action is not None
                           else min(width, self.sizeHint().width()))
        self.adjustSize()

    def set_shown(self, shown: bool, notify: bool = False):
        """Show or hide the floating toolbar; its tab stays."""
        self.setVisible(shown)
        self.place()
        if notify and self.on_toggled is not None:
            self.on_toggled(shown)

    def place(self):
        """The toolbar in the bottom left corner of the plot area, its tab
        at its right edge, at the bottom (in the corner when hidden)."""
        anchor = getattr(self, '_anchor', None)
        if anchor is None:
            return
        # As wide as the side cards (Guido): tighter buttons (refresh_icons),
        # and smaller icons where they still do not fit (platform styles).
        # Again at each placement: the size of the buttons is known only
        # once the window is shown
        self._fit_to_width(CARD_WIDTH)
        host = anchor.parentWidget()
        bottom = host.rect().bottom() + 1
        self.move(0, bottom - self.height())
        right = self.width() if not self.isHidden() else 0
        self.tab.move(right, bottom - TAB_SIZE)
        self.tab.setToolTip("Hide the toolbar" if right else "Show the toolbar")
        self.raise_()
        self.tab.raise_()

    def eventFilter(self, obj, event):
        if event.type() in (QEvent.Type.Resize, QEvent.Type.Move, QEvent.Type.Show):
            try:
                anchor = getattr(self, '_anchor', None)
                if anchor is not None and obj in (anchor, anchor.parentWidget()):
                    self.place()
            except RuntimeError:  # the plot already deleted (window closing)
                pass
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
