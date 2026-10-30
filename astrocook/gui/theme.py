"""The colours of the interface (Qt widgets), in one place.

Before phase 5-bis they were written by hand in a dozen stylesheets and
widgets. The colours of the plots are in ``plot_style.py``: the plots stay
white also in dark mode (Guido's decision), so they do not follow the
theme. The colours of the system palette (backgrounds, text, borders) are
read with ``palette_hex``.
"""
from PySide6.QtGui import QPalette
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


def palette_hex(role: QPalette.ColorRole) -> str:
    """A colour of the system palette (light or dark), as ``#rrggbb``."""
    return QApplication.palette().color(role).name()
