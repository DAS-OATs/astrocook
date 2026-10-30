"""The colours of the interface (Qt widgets), in one place.

Before phase 5-bis they were written by hand in a dozen stylesheets and
widgets. The colours of the plots are in ``plot_style.py``: the plots stay
white also in dark mode (Guido's decision), so they do not follow the
theme. The colours of the system palette (backgrounds, text, borders) are
read with ``palette_hex``.

Dark mode (phase 5-bis, theme PR B): the interface follows the palette of
the system. The constants below are the colours on a light palette; those
that would not read on a dark one have a dark variant (``DARK``), which
``themed`` returns when the palette is dark (``is_dark``). The widgets that
live long (the main window, the System Inspector, the toolbars) apply them
again when the system switches between light and dark
(``QEvent.ApplicationPaletteChange``).
"""
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication

#: The accent of Astrocook: the border of a field being edited, titles,
#: buttons that act on everything
ACCENT = "#296bff"
#: The selection: the highlighted group of the System Inspector, the
#: selected components and regions in the plots, the row of the Data
#: Inspector, the version in the About dialog
SELECTION = "#f5a100"
#: Text that warns (an oversampled grid, buttons that act on the systems)
WARNING = "#d35400"
#: Text of a choice that removes data
DANGER = "darkred"
#: A frozen parameter in the table of the System Inspector
FROZEN_TEXT = "#888888"
#: The hint of the empty main window ("Use 'File > Open Session...'")
HINT_TEXT = "#ff7f0e"
#: The colours of the links of the selected group in the System Inspector
#: (one per link, in order)
LINK_COLOURS = ["#296bff", "#d90077", "#008a00", "#7030a0", "#008080", "#800000"]
#: The flags of a component in the System Inspector (#517): violet for a
#: parameter at a bound of the fit, red for uncertainties that could not be
#: computed (Guido's choice; amber, proposed first, is the colour of the
#: highlighted group, SELECTION, and was confused with it in the screenshots)
FLAG_BOUND_COLOR = "#9467bd"
FLAG_NAN_COLOR = "#d62728"

#: Borders and buttons of the recipe dialogs
INPUT_BORDER = "#ccc"
BUTTON_BACKGROUND = "#e0e0e0"
BUTTON_BORDER = "#999"
BUTTON_PRESSED = "#d0d0d0"
#: The text on the selection, and the tooltips of the System Inspector
SELECTION_TEXT = "black"
TOOLTIP_BACKGROUND = "#f2f8f8"
TOOLTIP_TEXT = "black"
TOOLTIP_BORDER = "gray"
#: The console of the scripts
CONSOLE_BACKGROUND = "#222"
CONSOLE_TEXT = "#EEE"

#: The rule of a field being edited, shared by the stylesheets
FOCUS_BORDER = f"border: 1px solid {ACCENT};"


#: The variants on a dark palette of the colours that would not read on it
#: (light borders and buttons, dark red, the darker link colours). The
#: accent, the selection (with black text on it), the tooltips of the
#: System Inspector (light, with black text) and the console (dark) read
#: on both.
DARK = {
    DANGER: "#ff6b6b",
    WARNING: "#ff8c42",
    INPUT_BORDER: "#5a5a5a",
    BUTTON_BACKGROUND: "#4a4a4a",
    BUTTON_BORDER: "#777777",
    BUTTON_PRESSED: "#3a3a3a",
}
#: The link colours on a dark palette: the same hues, lighter
LINK_COLOURS_DARK = ["#6e9bff", "#ff5cc0", "#4cd964", "#c08cff", "#40d0d0", "#ff7070"]


def palette_hex(role: QPalette.ColorRole) -> str:
    """A colour of the system palette (light or dark), as ``#rrggbb``."""
    return QApplication.palette().color(role).name()


def is_dark() -> bool:
    """Whether the palette of the application is dark (the background of
    the windows darker than mid grey): on macOS in dark mode, or with a dark
    palette set by the user."""
    return QApplication.palette().color(QPalette.ColorRole.Window).lightness() < 128


def themed(colour: str) -> str:
    """``colour`` (one of the constants above), or its dark variant when
    the palette is dark."""
    return DARK.get(colour, colour) if is_dark() else colour


def link_colours() -> list:
    """The link colours for the current palette."""
    return LINK_COLOURS_DARK if is_dark() else LINK_COLOURS


# --- Appearance (View > Appearance) ------------------------------------------

#: The choices of View > Appearance, saved in the settings ("appearance")
APPEARANCES = ("system", "light", "dark")

#: The dark palette used where the platform does not switch by itself (Linux
#: without a dark theme, Qt offscreen); also used by the tests and the
#: screenshots
DARK_PALETTE = {
    QPalette.ColorRole.Window: "#2b2b2b", QPalette.ColorRole.WindowText: "#e6e6e6",
    QPalette.ColorRole.Base: "#1e1e1e", QPalette.ColorRole.AlternateBase: "#262626",
    QPalette.ColorRole.Text: "#e6e6e6", QPalette.ColorRole.Button: "#353535",
    QPalette.ColorRole.ButtonText: "#e6e6e6", QPalette.ColorRole.Highlight: "#2f65ca",
    QPalette.ColorRole.HighlightedText: "#ffffff", QPalette.ColorRole.Mid: "#555555",
    # Links (About) and the placeholders of the fields (script console):
    # left to the light palette they did not read on the dark background
    QPalette.ColorRole.Link: "#6e9bff", QPalette.ColorRole.LinkVisited: "#b28cff",
    QPalette.ColorRole.PlaceholderText: "#8a8a8a",
}

_palette_forced = False


def dark_palette() -> QPalette:
    """The palette of the application with the dark colours above."""
    palette = QPalette(QApplication.palette())
    for role, colour in DARK_PALETTE.items():
        palette.setColor(role, QColor(colour))
    return palette


def apply_appearance(mode: str):
    """Light, dark, or as the system (``APPEARANCES``).

    Qt is asked for the native appearance first (``QStyleHints.
    setColorScheme``, Qt >= 6.8: macOS and Windows switch their own
    controls); where the palette does not follow (Linux without a dark
    theme, offscreen), the palette is set here. "system" undoes a palette
    set here, giving back that of the system."""
    global _palette_forced
    app = QApplication.instance()
    scheme = {"system": Qt.ColorScheme.Unknown, "light": Qt.ColorScheme.Light,
              "dark": Qt.ColorScheme.Dark}[mode]
    hints = app.styleHints()
    if hasattr(hints, "setColorScheme"):
        hints.setColorScheme(scheme)
    if mode == "system":
        if _palette_forced:
            app.setPalette(QPalette())  # an empty palette: that of the system
            _palette_forced = False
    elif is_dark() != (mode == "dark"):
        if _palette_forced:
            app.setPalette(QPalette())
        if is_dark() != (mode == "dark"):
            app.setPalette(dark_palette() if mode == "dark" else app.style().standardPalette())
            _palette_forced = True
