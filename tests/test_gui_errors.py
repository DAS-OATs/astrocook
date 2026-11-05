"""Errors of the GUI (phase 6, PR 6b).

The hook for the errors that nothing handles (``gui/errors.py``), the
entries of the session log written once and not dropped silently, the
boundaries of the GUI that remain.
"""
import logging
import sys

import pytest

from astrocook.recipes.errors import RecipeFailed
from PySide6.QtWidgets import QMessageBox

from astrocook.core.structures import HistoryLogV2
from astrocook.gui import errors
from astrocook.gui import main_recipes

from test_main_recipes import _messages, _save_as, ready, window  # noqa: F401


def _raise(error):
    try:
        raise error
    except type(error):
        return sys.exc_info()


@pytest.fixture
def boxes(monkeypatch):
    """The dialogs of the hook, recorded instead of shown."""
    shown = []

    def exec_(self):
        shown.append({"text": self.text(), "informative": self.informativeText(),
                      "details": self.detailedText(), "title": self.windowTitle()})
        return QMessageBox.StandardButton.Ok
    monkeypatch.setattr(QMessageBox, "exec", exec_)
    return shown


def test_unexpected_error_is_logged_and_shown(qapp, boxes, caplog):
    with caplog.at_level(logging.ERROR, logger="astrocook.gui.errors"):
        errors.report_unexpected_error(*_raise(KeyError("bug")))
    [record] = caplog.records
    assert record.getMessage() == "Unexpected error: KeyError: 'bug'"
    assert record.exc_info is not None          # with the traceback
    [box] = boxes
    assert box["title"] == "Unexpected Error"
    assert box["text"] == "An unexpected error occurred: KeyError: 'bug'"
    assert "Traceback" in box["details"] and "KeyError: 'bug'" in box["details"]


def test_no_dialog_while_one_is_open(qapp, boxes, monkeypatch):
    # An error in a redraw could otherwise open one dialog per frame
    monkeypatch.setattr(errors, "_showing", True)
    errors.report_unexpected_error(*_raise(ValueError("again")))
    assert boxes == []


def test_install_excepthook(qapp, boxes, monkeypatch):
    previous_calls = []
    monkeypatch.setattr(sys, "excepthook", lambda *a: previous_calls.append(a[0]))
    previous = errors.install_excepthook()
    try:
        sys.excepthook(*_raise(RuntimeError("slot")))
        assert [b["text"] for b in boxes] == ["An unexpected error occurred: RuntimeError: slot"]
        # Ctrl-C goes to the previous hook, without a dialog
        sys.excepthook(*_raise(KeyboardInterrupt()))
        assert previous_calls == [KeyboardInterrupt] and len(boxes) == 1
    finally:
        sys.excepthook = previous


def test_launcher_installs_the_hook():
    import inspect
    import astrocook._run_astrocook as launcher
    assert "install_excepthook()" in inspect.getsource(launcher.main)


# --- The session log ---

class _History:
    def __init__(self, log_manager):
        self.log_manager = log_manager


def test_log_entry_written():
    history = _History(HistoryLogV2())
    main_recipes._add_log_entry(history, "smooth", {"sigma_kms": "10.0"})
    assert [e.recipe_name for e in history.log_manager.entries] == ["smooth"]


def test_unknown_log_manager_is_an_error(caplog):
    # It was skipped silently
    with caplog.at_level(logging.ERROR):
        main_recipes._add_log_entry(_History(object()), "smooth", {})
    assert ("The log of the session has an unknown type (object): smooth is not "
            "recorded.") in caplog.text


def test_failed_log_entry_keeps_log_and_state_together(ready, monkeypatch):  # noqa: F811
    # The entry is written before the new state is applied: when it fails,
    # the result is not applied and the error is shown (it was dropped
    # silently, while the result was applied)
    win = ready
    history = win.active_history
    states, entries = len(history.states), len(history.log_manager.entries)

    def broken(*args):
        raise RuntimeError("log")
    monkeypatch.setattr(history.log_manager, "add_entry", broken)
    shown = []
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a: shown.append(a[2])))
    state = history.current_state
    win._process_recipe_result((state.with_new_spectrum(state.spec), "smooth",
                                   {"sigma_kms": "10.0"}))
    assert shown == ["Could not show the result of the recipe:\nlog"]
    assert len(history.states) == states
    assert len(history.log_manager.entries) == entries


# --- Boundaries of the GUI ---

def test_failed_export_is_shown(ready, monkeypatch, tmp_path):  # noqa: F811
    # export_ascii returned False and nothing was shown
    win = ready
    _messages(win, monkeypatch, "Yes")
    _save_as(win, monkeypatch, tmp_path)
    recipes = type(win.active_history.current_state.file)
    warned = []
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda *a: warned.append(a[2])))
    # Since phase 6-bis (PR 3) the recipe raises the reason, shown here

    def fail(self, **kw):
        raise RecipeFailed("No space left on device")
    monkeypatch.setattr(recipes, "export_ascii", fail)
    win._on_recipe_requested("file", "export_ascii", {"targets": "y, dy"}, {})
    assert warned == ["The data could not be exported:\nNo space left on device."]
    # False: the targets name nothing to export
    monkeypatch.setattr(recipes, "export_ascii", lambda self, **kw: False)
    win._on_recipe_requested("file", "export_ascii", {"targets": "y, dy"}, {})
    assert warned[-1].startswith("There is nothing to export")

