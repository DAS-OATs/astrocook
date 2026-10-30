"""The theme of the interface (``gui/theme.py``) and the dark mode.

The colours of the interface are in one module (phase 5-bis, theme PR A);
the interface follows the palette of the system, light or dark, and applies
it again when the system switches (theme PR B). The plots stay white
(Guido's decision). A dark palette is set on the application to simulate
the dark mode of macOS.
"""
import pytest
from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication, QDialog, QMessageBox

from astrocook.gui import theme

from test_gui_smoke import SESSION_FILE

DARK_PALETTE = {
    QPalette.ColorRole.Window: "#2b2b2b", QPalette.ColorRole.WindowText: "#e6e6e6",
    QPalette.ColorRole.Base: "#1e1e1e", QPalette.ColorRole.AlternateBase: "#262626",
    QPalette.ColorRole.Text: "#e6e6e6", QPalette.ColorRole.Button: "#353535",
    QPalette.ColorRole.ButtonText: "#e6e6e6", QPalette.ColorRole.Highlight: "#2f65ca",
    QPalette.ColorRole.HighlightedText: "#ffffff", QPalette.ColorRole.Mid: "#555555",
}


def _dark_palette():
    palette = QPalette(QApplication.palette())
    for role, colour in DARK_PALETTE.items():
        palette.setColor(role, QColor(colour))
    return palette


@pytest.fixture
def dark(qapp):
    """The application with a dark palette, the light one restored after."""
    light = QPalette(qapp.palette())
    qapp.setPalette(_dark_palette())
    yield
    qapp.setPalette(light)


def _contrast(a: str, b: str) -> float:
    """WCAG contrast ratio of two colours."""
    def luminance(colour):
        c = QColor(colour)
        rgb = [v / 255 for v in (c.red(), c.green(), c.blue())]
        lin = [v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4 for v in rgb]
        return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    hi, lo = sorted((luminance(a), luminance(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


# --- The module ----------------------------------------------------------------

def test_light_by_default(qapp):
    assert not theme.is_dark()
    assert theme.themed(theme.DANGER) == theme.DANGER
    assert theme.link_colours() == theme.LINK_COLOURS


def test_dark_variants(dark):
    assert theme.is_dark()
    assert theme.themed(theme.DANGER) == theme.DARK[theme.DANGER]
    assert theme.themed(theme.ACCENT) == theme.ACCENT  # reads on both
    assert theme.link_colours() == theme.LINK_COLOURS_DARK
    assert len(theme.LINK_COLOURS_DARK) == len(theme.LINK_COLOURS)


@pytest.mark.parametrize("palette_dark", [False, True])
def test_text_colours_read_on_the_background(qapp, palette_dark):
    """The coloured texts (links, warnings, frozen parameters) against the
    background of the fields and tables: a contrast of at least 3 (WCAG
    for large or bold text; the link colours of the light palette on white
    are above 4.5)."""
    light = QPalette(qapp.palette())
    if palette_dark:
        qapp.setPalette(_dark_palette())
    try:
        base = theme.palette_hex(QPalette.ColorRole.Base)
        colours = theme.link_colours() + [theme.themed(theme.WARNING), theme.themed(theme.DANGER),
                                          theme.FROZEN_TEXT]
        for colour in colours:
            assert _contrast(colour, base) >= 3.0, (colour, base)
    finally:
        qapp.setPalette(light)


# --- The interface follows the palette ---------------------------------------

@pytest.fixture
def window(qtbot, monkeypatch, tmp_path):
    from astrocook.gui.main_window import MainWindowV2
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    for kind in ("warning", "critical", "information", "question"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(
            lambda *a, **k: QMessageBox.StandardButton.Ok))
    win = MainWindowV2(None, None)
    qtbot.addWidget(win)
    monkeypatch.setattr(win, "_show_custom_message", lambda *a, **k: None)
    win.show()
    win.open_session_from_path(str(SESSION_FILE))
    qtbot.wait(100)
    return win


def test_the_interface_switches_with_the_palette(qtbot, window, qapp):
    window._on_view_system_inspector()
    inspector = window.system_inspector
    light_base = theme.palette_hex(QPalette.ColorRole.Base)
    assert light_base in window.styleSheet() and light_base in inspector.z_in.styleSheet()
    icon_before = window.plot_viewer.toolbar.select_action.icon().pixmap(16, 16).toImage()
    light = QPalette(qapp.palette())
    qapp.setPalette(_dark_palette())
    try:
        qtbot.wait(50)
        dark_base = DARK_PALETTE[QPalette.ColorRole.Base]
        # The stylesheets of the main window and of the fields of the
        # System Inspector, made again with the dark palette
        assert dark_base in window.styleSheet()
        assert dark_base in inspector.z_in.styleSheet()
        assert dark_base in inspector.trans_in.styleSheet()
        # The icon of Select takes the colour of the text
        icon_after = window.plot_viewer.toolbar.select_action.icon().pixmap(16, 16).toImage()
        assert icon_after != icon_before
    finally:
        qapp.setPalette(light)
    qtbot.wait(50)
    assert light_base in window.styleSheet()


def test_link_colours_follow_the_palette(qapp):
    from PySide6.QtCore import Qt
    from astrocook.core.structures import ComponentDataV2, ParameterConstraintV2
    from astrocook.gui.system_inspector import COLUMNS, SystemTableModel
    comps = [ComponentDataV2(id=i, z=2.0 + 1e-4 * i, dz=None, logN=13.0, dlogN=None,
                             b=10.0, db=None, series="CIV") for i in (1, 2)]
    links = {comps[1].uuid: {"b": ParameterConstraintV2(
        is_free=False, target_uuid=comps[0].uuid, expression=f"p['{comps[0].uuid}'].b")}}
    model = SystemTableModel()
    model.update_data(comps, links)
    model.set_highlighted_uuids({comps[0].uuid}, {comps[1].uuid})

    def colour():
        brush = model.data(model.index(1, COLUMNS.index("b")), Qt.ForegroundRole)
        return brush.color().name()
    assert colour() == theme.LINK_COLOURS[0]
    light = QPalette(qapp.palette())
    qapp.setPalette(_dark_palette())
    try:
        # The same link, the same place in the list: its lighter colour
        assert colour() == theme.LINK_COLOURS_DARK[0]
    finally:
        qapp.setPalette(light)


def test_session_info_stylesheet(window, monkeypatch):
    # It was a plain string with the placeholders of an f-string, not valid
    # CSS (findings register F87)
    shown = []
    monkeypatch.setattr(QDialog, "exec", lambda self: shown.append(self.styleSheet()))
    window._on_session_info(window.active_history)
    assert shown and "{" + "palette" not in shown[0] and "{{" not in shown[0]
    assert theme.FOCUS_BORDER in shown[0]
    assert theme.palette_hex(QPalette.ColorRole.Base) in shown[0]
