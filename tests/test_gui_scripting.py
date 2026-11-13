"""The script console: ``app.run_recipe`` (F101)."""
from pathlib import Path

import pytest

from astrocook.gui.scripting import AppController
from astrocook.session import Session
from astrocook.session.session import load_session_from_file

SESSION_FILE = Path(__file__).parent / "data" / "J0103-1305_metallines_lya_new_3.acs2"


@pytest.fixture
def controller(qtbot):
    return AppController(main_window=None)


@pytest.fixture
def session():
    return load_session_from_file(str(SESSION_FILE), "s", None, "auto")


def test_run_recipe_returns_the_new_session(controller, session):
    """The console's ``app.run_recipe`` used to fail on a keyword that the
    window's dispatcher does not take (F101)."""
    messages = []
    controller.log_message.connect(messages.append)
    new = controller.run_recipe(session, "edit", "set_properties", z_em="4.072")
    assert isinstance(new, Session)
    assert new is not session
    assert new.spec._data.z_em == pytest.approx(4.072)
    assert messages == ["Running edit.set_properties..."]


def test_run_recipe_unknown_recipe_raises(controller, session):
    with pytest.raises((KeyError, ValueError)):
        controller.run_recipe(session, "edit", "no_such_recipe")
