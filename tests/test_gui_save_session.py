"""GUI tests for saving a session (offscreen, pytest-qt)."""
import os

import numpy as np
import astropy.units as au
import pytest
from PySide6.QtWidgets import QFileDialog, QMessageBox

from astrocook.core.session import SessionV2
from astrocook.core.session_manager import SessionHistory
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, HistoryLogV2, SpectrumDataV2


@pytest.fixture
def window_with_v2_session(qtbot, monkeypatch, tmp_path):
    # Keep the settings of the test run away from the user's ones
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    # Import only once qtbot has created the QApplication: the GUI selects the
    # QtAgg matplotlib backend at import time, which fails without one.
    from astrocook.gui.main_window import MainWindowV2
    window = MainWindowV2(None, None)
    qtbot.addWidget(window)

    x = np.linspace(400.0, 410.0, 50)
    one = np.ones_like(x)
    data = SpectrumDataV2(
        x=DataColumnV2(x, au.nm),
        xmin=DataColumnV2(x - 0.1, au.nm),
        xmax=DataColumnV2(x + 0.1, au.nm),
        y=DataColumnV2(one, au.dimensionless_unscaled),
        dy=DataColumnV2(0.1 * one, au.dimensionless_unscaled),
    )
    session = SessionV2(name="gui_save_test", gui=None, spec=SpectrumV2(data))
    history = SessionHistory(session, HistoryLogV2())
    return window, history


def _choose_file(monkeypatch, path):
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        lambda *args, **kwargs: (str(path), ""))


def _answer_dialog(monkeypatch, window, answer):
    calls = []

    def fake_message(*args, **kwargs):
        calls.append(kwargs)
        return answer

    monkeypatch.setattr(window, "_show_custom_message", fake_message)
    return calls


def test_saving_v2_log_to_acs_asks_and_can_be_cancelled(
        window_with_v2_session, monkeypatch, tmp_path):
    window, history = window_with_v2_session
    path = tmp_path / "session.acs"
    _choose_file(monkeypatch, path)
    calls = _answer_dialog(monkeypatch, window, QMessageBox.StandardButton.Cancel)

    assert window._save_session(history) is False
    assert len(calls) == 1
    assert ".acs2" in calls[0]["text"]
    assert not os.path.exists(path)


def test_saving_v2_log_to_acs_proceeds_when_confirmed(
        window_with_v2_session, monkeypatch, tmp_path):
    window, history = window_with_v2_session
    path = tmp_path / "session.acs"
    _choose_file(monkeypatch, path)
    calls = _answer_dialog(monkeypatch, window, QMessageBox.StandardButton.Save)

    assert window._save_session(history) is True
    assert len(calls) == 1
    assert os.path.exists(path)


def test_saving_to_acs2_does_not_ask(window_with_v2_session, monkeypatch, tmp_path):
    window, history = window_with_v2_session
    path = tmp_path / "session.acs2"
    _choose_file(monkeypatch, path)
    calls = _answer_dialog(monkeypatch, window, QMessageBox.StandardButton.Cancel)

    assert window._save_session(history) is True
    assert calls == []
    assert os.path.exists(path)


def test_loaded_log_is_kept_by_the_gui(window_with_v2_session, tmp_path):
    # F26: the GUI used to give every loaded session a new, empty log
    window, history = window_with_v2_session
    history.log_manager.add_entry("smooth", {"sigma_kms": "20"})
    path = tmp_path / "session.acs2"
    history.current_state.save(str(path))
    loaded = SessionV2.open_new(str(path), "reloaded", None, "auto")

    window.add_session(loaded)
    log = window.active_history.log_manager
    assert [e.recipe_name for e in log.entries] == ["smooth"]

    # A duplicate does not share the log of its source history
    window.add_session(window.active_history.current_state)
    assert window.active_history.log_manager is not log
    assert window.active_history.log_manager.entries == []
