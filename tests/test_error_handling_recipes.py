"""Errors of the recipes reported once, with their context (phase 6, PR 5).

PR 5a: the recipes catch only the expected exceptions inside, their entry
points log an unexpected error once, with the traceback, and their
messages follow the rule of phase 6 (``ROADMAP.md``).
"""
import logging

import pytest

from astrocook.core.session import SessionV2

from test_recipes_absorbers import _systs
from test_recipes_edit_flux import _spec


def _two():
    session = SessionV2("s", None, spec=_spec(),
                        systs=_systs(("CIV", 2.0, 13.0, 10.0), ("SiIV", 2.0001, 12.5, 8.0)))
    return session, [c.uuid for c in session.systs.components]


def test_linked_value_uses_the_functions_of_the_fit(caplog):
    # F98: the update of the linked value offered 'sqrt' and not 'ln', the
    # fit the opposite; both use LINK_FUNCTIONS now
    session, (a, b) = _two()
    with caplog.at_level(logging.WARNING):
        new = session.absorbers.update_constraint(
            b, "logN", is_free=False, expression=f"ln(p['{a}'].logN)", target_uuid=a)
    assert new != 0
    assert new.systs.get_component_by_uuid(b).logN == pytest.approx(2.564949357)
    assert caplog.text == ""


def test_link_with_sqrt_is_refused(caplog):
    session, (a, b) = _two()
    with caplog.at_level(logging.ERROR):
        assert session.absorbers.update_constraint(
            b, "b", is_free=False, expression=f"sqrt(p['{a}'].b)", target_uuid=a) == 0
    assert "Could not set the constraint" in caplog.text
    assert "NameError" in caplog.text
