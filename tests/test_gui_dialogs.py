"""The accessory widgets of the main window: the progress dialog of the
recipes and the log handlers that feed it and collect the warnings, the
About dialog, the collapsible sections of the cards and the label that
wraps.

Characterization tests written before these classes were moved out of
``main_window.py`` (phase 5, split of the GUI).
"""
import logging
import sys

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialogButtonBox, QLabel, QPushButton

from astrocook.gui import main_window
from astrocook.gui.about_dialog import AboutDialog
from astrocook.gui.floating_cards import CollapsibleSection, WrappedLabel
from astrocook.gui.recipe_progress import (LogSignalBridge, QtLogHandler, RecipeProgressDialog,
                                           WarningCaptureHandler)


def _module(cls):
    return sys.modules[cls.__module__]


# --- Progress of a recipe --------------------------------------------------------------

def test_progress_dialog(qtbot):
    dialog = RecipeProgressDialog("fit_systems")
    qtbot.addWidget(dialog)
    assert dialog.windowTitle() == "Running fit_systems" and dialog.isModal()
    assert not dialog.windowFlags() & Qt.WindowType.WindowCloseButtonHint
    assert dialog.progress_bar.value() == 0 and dialog.status_label.text() == "Initializing..."
    dialog.update_log(">> PROGRESS: 42")
    assert dialog.progress_bar.value() == 42
    assert dialog.status_label.text() == "Initializing..."  # not a status
    dialog.update_log(">> PROGRESS: many")  # not a number: ignored
    assert dialog.progress_bar.value() == 42
    dialog.update_log("  Fitting group 3 of 7\n")
    assert dialog.status_label.text() == "Fitting group 3 of 7"
    dialog.update_log("   ")  # blank: the last status stays
    assert dialog.status_label.text() == "Fitting group 3 of 7"


def test_progress_dialog_stop(qtbot):
    dialog = RecipeProgressDialog("fit_systems")
    qtbot.addWidget(dialog)
    with qtbot.waitSignal(dialog.stop_requested):
        dialog.stop_btn.click()
    assert not dialog.stop_btn.isEnabled() and dialog.stop_btn.text() == "Stopping..."
    assert dialog.status_label.text() == "Attempting to stop gracefully..."


def test_log_handler_feeds_the_dialog(qtbot):
    dialog = RecipeProgressDialog("fit_systems")
    qtbot.addWidget(dialog)
    bridge = LogSignalBridge()
    bridge.new_log_message.connect(dialog.update_log)
    handler = QtLogHandler(bridge)
    handler.setFormatter(logging.Formatter("%(levelname)s %(message)s"))
    logger = logging.getLogger("astrocook.test_gui_dialogs")
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    try:
        logger.info(">> PROGRESS: 7")
        logger.info("Group done")
    finally:
        logger.removeHandler(handler)
    assert dialog.progress_bar.value() == 7
    assert dialog.status_label.text() == "INFO Group done"


def test_warnings_are_captured():
    handler = WarningCaptureHandler()
    logger = logging.getLogger("astrocook.test_gui_dialogs.warnings")
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    try:
        logger.info("fine")
        logger.warning("b at the bound")
        logger.error("failed")
    finally:
        logger.removeHandler(handler)
    assert handler.warnings == ["b at the bound", "failed"]


# --- About ----------------------------------------------------------------------------------

def test_about_dialog(qtbot, monkeypatch):
    from astrocook import __version__
    opened = []
    services = _module(AboutDialog).QDesktopServices
    monkeypatch.setattr(services, "openUrl", staticmethod(lambda url: opened.append(url.toString())))
    dialog = AboutDialog()
    qtbot.addWidget(dialog)
    assert dialog.windowTitle() == "About Astrocook"
    texts = [label.text() for label in dialog.findChildren(QLabel)]
    assert f"Version: {__version__}" in texts
    assert "A thousand recipes to cook a spectrum" in texts
    assert any("Claudia Cupani" in t for t in texts)
    buttons = {b.text(): b for b in dialog.findChildren(QPushButton)}
    buttons["Report Issue"].click()
    buttons["Cite on ADS"].click()
    assert opened == ["https://github.com/DAS-OATs/astrocook/issues",
                      "https://ui.adsabs.harvard.edu/abs/2020SPIE11452E..1UC/exportcitation"]
    box = dialog.findChild(QDialogButtonBox)
    box.button(QDialogButtonBox.StandardButton.Ok).click()
    assert dialog.result() == dialog.DialogCode.Accepted


def test_about_from_the_window(qtbot, monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    shown = []
    monkeypatch.setattr(AboutDialog, "exec", lambda self: shown.append(self))
    window = main_window.MainWindowV2(None, None)
    qtbot.addWidget(window)
    window._on_about_requested()
    assert len(shown) == 1 and shown[0].parent() is window


# --- Cards ------------------------------------------------------------------------------------

def test_collapsible_section(qtbot):
    section = CollapsibleSection("Plot Elements", expanded=False)
    qtbot.addWidget(section)
    section.show()
    assert section.header.text() == "Plot Elements"
    assert not section.isExpanded() and not section.content.isVisible()
    with qtbot.waitSignal(section.toggled) as signal:
        section.header.click()
    assert signal.args == [True] and section.isExpanded() and section.content.isVisible()
    section.setExpanded(False)
    assert not section.content.isVisible()


@pytest.mark.parametrize("width", [80, 400])
def test_wrapped_label_asks_for_its_height(qtbot, width):
    label = WrappedLabel("Shift+click to add a knot, right click to remove it, drag to move it.")
    label.setWordWrap(True)
    qtbot.addWidget(label)
    label.resize(width, 10)
    label.show()
    qtbot.wait(10)
    assert label.minimumHeight() == label.heightForWidth(width)
