"""Issue #450: linked parameters must be visible without hovering.

The inspector table colours, within the selected group, each linked
parameter and its target with the same colour, and shows frozen
parameters in grey.
"""
import pytest
from PySide6.QtCore import Qt

from astrocook.core.structures import ComponentDataV2, ParameterConstraintV2

GREY = "#888888"


@pytest.fixture
def model(qapp):
    from astrocook.gui.system_inspector import SystemTableModel
    return SystemTableModel()


def _comps():
    return [ComponentDataV2(id=i, z=2.0 + 1e-4 * i, dz=None, logN=13.0,
                            dlogN=None, b=10.0, db=None, series="CIV")
            for i in (1, 2, 3)]


def _colour(model, row, col_name):
    from astrocook.gui.system_inspector import COLUMNS
    brush = model.data(model.index(row, COLUMNS.index(col_name)), Qt.ForegroundRole)
    return brush.color().name() if brush is not None else None


def _link_b(target, dependent):
    # dependent.b = target.b
    return {dependent.uuid: {"b": ParameterConstraintV2(
        is_free=False, target_uuid=target.uuid,
        expression=f"p['{target.uuid}'].b")}}


def test_linked_pair_shares_a_colour_in_the_selected_group(model):
    c1, c2, c3 = _comps()
    model.update_data([c1, c2, c3], _link_b(c1, c2))
    model.set_highlighted_uuids({c1.uuid}, {c2.uuid})

    target, dependent = _colour(model, 0, "b"), _colour(model, 1, "b")
    assert target is not None and target != GREY
    assert dependent == target
    # Unlinked parameters and components keep the default colour
    assert _colour(model, 0, "logN") is None
    assert _colour(model, 2, "b") is None


def test_links_are_not_coloured_outside_the_selected_group(model):
    c1, c2, c3 = _comps()
    model.update_data([c1, c2, c3], _link_b(c1, c2))
    model.set_highlighted_uuids({c3.uuid}, set())
    assert _colour(model, 0, "b") is None
    assert _colour(model, 1, "b") is None


def test_frozen_parameter_is_grey(model):
    c1, c2, c3 = _comps()
    frozen = {c3.uuid: {"logN": ParameterConstraintV2(is_free=False)}}
    model.update_data([c1, c2, c3], frozen)
    model.set_highlighted_uuids(set(), set())
    assert _colour(model, 2, "logN") == GREY


@pytest.mark.xfail(strict=True, reason=(
    "set_highlighted_uuids returns early when the selection is unchanged, so "
    "the link colours are not recomputed after update_data adds a link to the "
    "selected group. Findings register F13 (ROADMAP.md, section 3.2)."))
def test_link_created_with_group_selected_is_coloured(model):
    c1, c2, c3 = _comps()
    model.update_data([c1, c2, c3], {})
    model.set_highlighted_uuids({c1.uuid}, {c2.uuid})
    # The user links c2.b to c1.b while the group stays selected; the
    # inspector calls update_data and then set_highlighted_uuids again
    model.update_data([c1, c2, c3], _link_b(c1, c2))
    model.set_highlighted_uuids({c1.uuid}, {c2.uuid})
    assert _colour(model, 1, "b") is not None
