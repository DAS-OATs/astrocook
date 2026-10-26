"""Tests of the floating cards of the main window (offscreen, pytest-qt).

The session list (left) and the plot controls (right) float over the plot,
so that opening them does not change the aspect of the spectrum. Each is a
card as tall as its contents, below the menu bar, with a tab just outside
its edge; it fades out while the mouse is away from it. The sections of the
plot controls are collapsible and remember their state. While the Plot
Elements section is open, its checkboxes show the plotted styles and those
entries are left out of the matplotlib legend.
"""
from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QEvent, QSettings
from PySide6.QtWidgets import QApplication, QMessageBox

DATA = Path(__file__).parent / "data"
SESSION_FILE = DATA / "J0103-1305_metallines_lya_new_3.acs2"
OTHER_SESSION_FILE = DATA / "J0103-1305_metallines_lya_new_2.acs2"
ANIMATION_TIMEOUT = 5000


@pytest.fixture
def window(qtbot, monkeypatch, tmp_path):
    # Keep the settings of the test run away from the user's ones
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    for kind in ("warning", "critical", "information", "question"):
        monkeypatch.setattr(QMessageBox, kind,
                            staticmethod(lambda *a, **k: QMessageBox.StandardButton.Ok))
    # Sections start from their defaults, whatever an earlier test stored
    QSettings("Astrocook", "AstrocookV2").remove("plot_controls")
    from astrocook.gui.main_window import MainWindowV2
    win = MainWindowV2(None, None)
    qtbot.addWidget(win)
    win.resize(1400, 850)
    win.show()
    win.open_session_from_path(str(SESSION_FILE))
    yield win
    QSettings("Astrocook", "AstrocookV2").remove("plot_controls")


def _open(qtbot, win, side):
    card = win.right_sidebar_widget if side == "right" else win.left_sidebar_widget
    toggle = win._toggle_right_sidebar if side == "right" else win._toggle_left_sidebar
    toggle()
    qtbot.waitUntil(lambda: card.width() > 0 and card.isVisible(), timeout=ANIMATION_TIMEOUT)
    animation = win.right_animation_group if side == "right" else win.left_animation_group
    qtbot.waitUntil(lambda: animation.state() == animation.State.Stopped, timeout=ANIMATION_TIMEOUT)
    return card


def _sections(win):
    from astrocook.gui.main_window import CollapsibleSection
    return {s.header.text(): s
            for s in win.right_sidebar_widget.findChildren(CollapsibleSection)}


def _legend_labels(win):
    legend = win.plot_viewer.canvas.axes.get_legend()
    return [] if legend is None else [t.get_text() for t in legend.get_texts()]


def _is_blank(pixmap):
    image = pixmap.toImage()
    return all(image.pixelColor(x, y).alpha() == 0
               for x in range(image.width()) for y in range(image.height()))


# --- Geometry --------------------------------------------------------------

def test_plot_controls_card_is_as_tall_as_its_contents(qtbot, window):
    from astrocook.gui.main_window import CARD_MARGIN
    card = _open(qtbot, window, "right")
    top = window.menuBar().height() + CARD_MARGIN
    assert card.y() == top
    assert card.height() < window.height() - top - CARD_MARGIN
    assert card.x() + card.width() == window.width()

    # Expanding a section makes the card taller, collapsing it shorter
    before = card.height()
    _sections(window)["Redshift Cursor"].setExpanded(True)
    qtbot.waitUntil(lambda: card.height() > before, timeout=ANIMATION_TIMEOUT)
    _sections(window)["Redshift Cursor"].setExpanded(False)
    qtbot.waitUntil(lambda: card.height() == before, timeout=ANIMATION_TIMEOUT)


def test_session_card_grows_with_the_sessions(qtbot, window):
    card = _open(qtbot, window, "left")
    one = card.height()
    window.open_session_from_path(str(OTHER_SESSION_FILE))
    qtbot.waitUntil(lambda: card.height() > one, timeout=ANIMATION_TIMEOUT)
    assert card.height() < window.height()


def test_tabs_are_just_outside_the_cards(qtbot, window):
    right = _open(qtbot, window, "right")
    left = _open(qtbot, window, "left")
    right_tab = window.plot_controls_collapse_button
    left_tab = window.session_collapse_button
    assert right_tab.x() + right_tab.width() == right.x()
    assert left_tab.x() == left.x() + left.width()
    assert right_tab.y() == right.y() and left_tab.y() == left.y()


def test_closed_tabs_stay_above_the_axes(qtbot, window):
    # A closed tab sits in the corner of the window and does not cover the
    # corner of the plot
    canvas = window.plot_viewer.canvas
    window._trigger_replot()
    qtbot.wait(200)
    axes_top = canvas.height() - canvas.axes.get_window_extent().y1 / canvas.devicePixelRatioF()
    for tab in (window.plot_controls_collapse_button, window.session_collapse_button):
        tab_bottom = canvas.mapFrom(window, tab.geometry().bottomLeft()).y()
        assert tab_bottom < axes_top


# --- Sections --------------------------------------------------------------

def test_sections_default_and_remembered_state(qtbot, window):
    sections = _sections(window)
    assert {name: s.isExpanded() for name, s in sections.items()} == {
        "Plot Elements": True, "View Controls": True, "Axis Limits": True,
        "Redshift Cursor": False, "Edit Continuum": False}
    sections["View Controls"].setExpanded(False)
    assert not sections["View Controls"].content.isVisible()

    from astrocook.gui.main_window import MainWindowV2
    other = MainWindowV2(None, None)
    qtbot.addWidget(other)
    assert not _sections(other)["View Controls"].isExpanded()
    assert _sections(other)["Axis Limits"].isExpanded()


# --- Legend and swatches ---------------------------------------------------

def test_panel_entries_leave_the_legend_while_shown(qtbot, window):
    assert "Spectrum" in _legend_labels(window)

    _open(qtbot, window, "right")
    qtbot.waitUntil(lambda: "Spectrum" not in _legend_labels(window), timeout=ANIMATION_TIMEOUT)

    # With Plot Elements closed, the legend is complete again
    _sections(window)["Plot Elements"].setExpanded(False)
    qtbot.waitUntil(lambda: "Spectrum" in _legend_labels(window), timeout=ANIMATION_TIMEOUT)


def test_swatches_show_what_is_drawn(qtbot, window):
    window._trigger_replot()
    # Spectrum, error and systems are drawn
    assert not _is_blank(window.spectrum_swatch.pixmap())
    assert not _is_blank(window.error_checkbox.icon().pixmap(window.error_checkbox.iconSize()))
    assert not _is_blank(window.systems_checkbox.icon().pixmap(window.systems_checkbox.iconSize()))
    # The session has no continuum: the checked box gets a blank swatch
    assert window.active_history.current_state.spec.cont is None
    assert window.continuum_checkbox.isChecked()
    assert _is_blank(window.continuum_checkbox.icon().pixmap(window.continuum_checkbox.iconSize()))


def test_swatch_takes_the_colour_of_the_line(window):
    from matplotlib.lines import Line2D
    pixmap = window._swatch([Line2D([0, 1], [0, 1], color="#ff0000", ls="--")])
    image = pixmap.toImage()
    colours = {image.pixelColor(x, y).name()
               for x in range(image.width()) for y in range(image.height())
               if image.pixelColor(x, y).alpha() == 255}
    assert colours == {"#ff0000"}


# --- Fading ----------------------------------------------------------------

def test_card_fades_out_and_back(qtbot, window):
    from astrocook.gui.main_window import CARD_FADED
    card = _open(qtbot, window, "right")
    fader = window._plot_controls_fader
    assert fader.effect.opacity() == 1.0

    QApplication.sendEvent(card, QEvent(QEvent.Type.Leave))
    qtbot.waitUntil(lambda: np.isclose(fader.effect.opacity(), CARD_FADED), timeout=ANIMATION_TIMEOUT)

    QApplication.sendEvent(card, QEvent(QEvent.Type.Enter))
    qtbot.waitUntil(lambda: np.isclose(fader.effect.opacity(), 1.0), timeout=ANIMATION_TIMEOUT)


def test_opening_a_card_shows_it_opaque(qtbot, window):
    window._set_sessions_faded(True, animate=False)
    _open(qtbot, window, "left")
    assert window._sessions_fader.effect.opacity() == 1.0


def test_tab_fades_with_its_card(qtbot, window):
    from astrocook.gui.main_window import CARD_FADED
    card = _open(qtbot, window, "right")
    tab = window.plot_controls_collapse_button
    QApplication.sendEvent(card, QEvent(QEvent.Type.Leave))
    qtbot.waitUntil(lambda: np.isclose(tab.graphicsEffect().opacity(), CARD_FADED),
                    timeout=ANIMATION_TIMEOUT)

    # Entering the tab counts as entering the card
    QApplication.sendEvent(tab, QEvent(QEvent.Type.Enter))
    qtbot.waitUntil(lambda: np.isclose(window._plot_controls_fader.effect.opacity(), 1.0),
                    timeout=ANIMATION_TIMEOUT)
    assert np.isclose(tab.graphicsEffect().opacity(), 1.0)


def test_closed_card_leaves_its_tab_opaque(qtbot, window):
    card = _open(qtbot, window, "right")
    window._set_plot_controls_faded(True, animate=False)
    window._toggle_right_sidebar()
    qtbot.waitUntil(lambda: not card.isVisible(), timeout=ANIMATION_TIMEOUT)
    tab = window.plot_controls_collapse_button
    assert tab.graphicsEffect().opacity() == 1.0
    # With the card closed, leaving the tab does not fade it
    QApplication.sendEvent(tab, QEvent(QEvent.Type.Leave))
    qtbot.wait(800)
    assert tab.graphicsEffect().opacity() == 1.0
