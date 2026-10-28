"""The floating cards of the main window (the session list on the left, the
plot controls on the right): the fading of a card while the mouse is away,
its collapsible sections and their disclosure triangle, and a label that
wraps.

Moved out of ``main_window.py`` in phase 5 (split of the GUI).
"""
import qtawesome as qta
from PySide6.QtCore import QEvent, QObject, QPropertyAnimation, QSize, QTimer, Qt, Signal
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import (QApplication, QGraphicsOpacityEffect, QLabel, QToolButton,
                               QVBoxLayout, QWidget)

CARD_FADED = 0.35      # Opacity of a card while the mouse is away
CARD_FADE_DELAY = 500  # ms before fading out after the mouse leaves
CARD_FADE_TIME = 200   # ms of the fade animation


class CardFader(QObject):
    """Fades a floating card out while the mouse is away from it.

    The card comes back to full opacity when the mouse enters it. Its tab
    (``add_handle``) fades with it and counts as part of it. A closed card
    does not fade, so its tab stays opaque. It also calls ``on_layout`` when
    the card's contents change size (a section expanded or collapsed), so
    that the window can resize the card.
    """

    def __init__(self, card: QWidget, on_layout=None):
        super().__init__(card)
        self.card = card
        self.on_layout = on_layout
        self.effect = QGraphicsOpacityEffect(card)
        self.effect.setOpacity(1.0)
        card.setGraphicsEffect(self.effect)
        self.anim = QPropertyAnimation(self.effect, b"opacity", self)
        self.anim.setDuration(CARD_FADE_TIME)
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.setInterval(CARD_FADE_DELAY)
        self.timer.timeout.connect(self._fade_if_left)
        self.handles = []
        card.installEventFilter(self)

    def add_handle(self, widget: QWidget):
        """A widget (the card's tab) that fades with the card."""
        effect = QGraphicsOpacityEffect(widget)
        effect.setOpacity(self.effect.opacity())
        widget.setGraphicsEffect(effect)
        self.effect.opacityChanged.connect(effect.setOpacity)
        widget.installEventFilter(self)
        self.handles.append(widget)

    def eventFilter(self, obj, event):
        etype = event.type()
        if etype == QEvent.Type.Enter:
            self.timer.stop()
            self.set_faded(False)
        elif etype == QEvent.Type.Leave:
            self.timer.start()
        elif (etype == QEvent.Type.LayoutRequest and obj is self.card
              and self.on_layout is not None):
            QTimer.singleShot(0, self.on_layout)
        return super().eventFilter(obj, event)

    def _fade_if_left(self):
        # A combo box popup or a context menu takes the mouse away from the
        # card without the user leaving it: wait until it is closed
        if not self.card.isVisible():
            return
        if (QApplication.activePopupWidget() is not None or self.card.underMouse()
                or any(h.underMouse() for h in self.handles)):
            self.timer.start()
            return
        self.set_faded(True)

    def set_faded(self, faded: bool, animate: bool = True):
        target = CARD_FADED if faded else 1.0
        self.anim.stop()
        if animate:
            self.anim.setStartValue(self.effect.opacity())
            self.anim.setEndValue(target)
            self.anim.start()
        else:
            self.effect.setOpacity(target)


def disclosure_icon(direction: str) -> QIcon:
    """The triangle of the headers of the card sections.

    ``direction`` is 'right' or 'down'. A qtawesome glyph rather than the
    style's arrows, so that it looks the same on every platform.
    """
    color = QApplication.palette().color(QApplication.palette().ColorRole.ButtonText)
    # The glyph fills only the middle of its box: scaled to read at 16 px
    return qta.icon(f'mdi.menu-{direction}', color=color, scale_factor=1.8)


class WrappedLabel(QLabel):
    """A label that wraps and asks for the height its text needs.

    Layouts give a wrapped QLabel its height for some other width, so the
    last line can be cut; the card of the plot controls, as tall as its
    contents, would then cut it off.
    """

    def resizeEvent(self, event):
        super().resizeEvent(event)
        needed = self.heightForWidth(self.width())
        if needed > 0 and needed != self.minimumHeight():
            self.setMinimumHeight(needed)


class CollapsibleSection(QWidget):
    """A section of a floating card, shown or hidden by its header."""

    toggled = Signal(bool)

    def __init__(self, title: str, expanded: bool = True, parent=None):
        super().__init__(parent)
        self.header = QToolButton()
        self.header.setText(title)
        self.header.setCheckable(True)
        self.header.setChecked(expanded)
        self.header.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.header.setIcon(disclosure_icon('down' if expanded else 'right'))
        self.header.setIconSize(QSize(16, 16))
        self.header.setAutoRaise(True)
        self.header.setStyleSheet("QToolButton { border: none; font-weight: bold; padding: 0px; }")
        self.content = QWidget()
        self.content.setVisible(expanded)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        layout.addWidget(self.header)
        layout.addWidget(self.content)
        self.header.toggled.connect(self._on_toggled)

    def setContentLayout(self, layout):
        layout.setContentsMargins(4, 0, 0, 0)
        self.content.setLayout(layout)

    def isExpanded(self) -> bool:
        return self.header.isChecked()

    def setExpanded(self, expanded: bool):
        self.header.setChecked(expanded)

    def _on_toggled(self, expanded: bool):
        self.header.setIcon(disclosure_icon('down' if expanded else 'right'))
        self.content.setVisible(expanded)
        self.toggled.emit(expanded)
