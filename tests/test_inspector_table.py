"""The table of the System Inspector: what the model gives Qt for each role,
and the columns fitted to their contents.

Characterization tests written before the table was made faster (phase 5:
the columns were fitted by Qt with ``ResizeToContents``, which queried every
row several times): the values of the model and the fit of the columns must
not change. The table moved to ``gui/inspector_table.py`` later in phase 5
(split of the GUI), with tests for the editing of a cell and the proxy.
"""
import pytest
from PySide6.QtCore import QModelIndex, Qt

from astrocook.core.structures import ComponentDataV2, ParameterConstraintV2


@pytest.fixture
def model(qapp):
    from astrocook.gui.inspector_table import SystemTableModel
    return SystemTableModel()


def _comps():
    c1 = ComponentDataV2(id=1, z=2.0, dz=1.5e-5, logN=13.0, dlogN=0.05,
                         b=10.0, db=1.25, series="CIV",
                         flags=("b_at_lower_bound",))
    c2 = ComponentDataV2(id=2, z=2.0001, dz=None, logN=12.5, dlogN=None,
                         b=10.0, db=None, series="CIV", flags=("nan_errors",))
    c3 = ComponentDataV2(id=3, z=2.5, dz=2e-6, logN=14.2, dlogN=0.1,
                         b=20.0, db=2.0, series="Ly_a")
    return c1, c2, c3


def _constraints(c1, c2, c3):
    return {
        # c2.b = c1.b
        c2.uuid: {"b": ParameterConstraintV2(
            is_free=False, target_uuid=c1.uuid,
            expression=f"p['{c1.uuid}'].b")},
        # c3.logN frozen
        c3.uuid: {"logN": ParameterConstraintV2(is_free=False)},
    }


def _data(model, row, col_name, role):
    from astrocook.gui.inspector_table import COLUMNS
    return model.data(model.index(row, COLUMNS.index(col_name)), role)


@pytest.fixture
def filled(model):
    c1, c2, c3 = _comps()
    model.update_data([c1, c2, c3], _constraints(c1, c2, c3))
    model.set_highlighted_uuids({c1.uuid}, {c2.uuid})
    return model


def test_display_and_edit(filled):
    assert _data(filled, 0, "ID", Qt.DisplayRole) == "⚠ 1"  # flagged
    assert _data(filled, 2, "ID", Qt.DisplayRole) == "3"
    assert _data(filled, 0, "Trans.", Qt.DisplayRole) == "CIV"
    # The error of z in velocity, c dz / (1 + z) (visual review of phase 5-bis)
    assert _data(filled, 0, "z", Qt.DisplayRole) == "2.000000 ± 1.50 km/s"
    assert "± 1.50e-05 (1.50 km/s)" in _data(filled, 0, "z", Qt.ToolTipRole)
    assert _data(filled, 0, "logN", Qt.DisplayRole) == "13.000 ± 0.050"
    assert _data(filled, 0, "b", Qt.DisplayRole) == "10.00 ± 1.25"
    # No error without one, for a linked or a frozen parameter
    assert _data(filled, 1, "z", Qt.DisplayRole) == "2.000100"
    assert _data(filled, 1, "b", Qt.DisplayRole) == "10.00"
    assert _data(filled, 2, "logN", Qt.DisplayRole) == "14.200"
    assert _data(filled, 0, "z", Qt.EditRole) == "2.000000"
    assert _data(filled, 0, "Trans.", Qt.EditRole) == "CIV"


def test_sorting_and_alignment(filled):
    assert _data(filled, 2, "z", Qt.UserRole) == 2.5
    assert _data(filled, 0, "ID", Qt.UserRole) == 1.0
    assert _data(filled, 0, "Trans.", Qt.UserRole) == "CIV"
    for col in ("ID", "Trans.", "z", "btur"):
        assert _data(filled, 0, col, Qt.TextAlignmentRole) == Qt.AlignCenter


def test_fonts_and_colours(filled):
    from astrocook.gui.inspector_table import FLAG_BOUND_COLOR, FLAG_NAN_COLOR
    assert _data(filled, 1, "b", Qt.FontRole).bold()          # linked
    assert _data(filled, 2, "logN", Qt.FontRole).italic()     # frozen
    assert _data(filled, 0, "b", Qt.FontRole) is None
    assert _data(filled, 0, "ID", Qt.FontRole) is None
    assert _data(filled, 2, "logN", Qt.ForegroundRole).color().name() == "#888888"
    target = _data(filled, 0, "b", Qt.ForegroundRole).color().name()
    assert _data(filled, 1, "b", Qt.ForegroundRole).color().name() == target
    assert _data(filled, 0, "ID", Qt.ForegroundRole) is None
    # Background: the flags, then the highlighted group
    assert _data(filled, 0, "b", Qt.BackgroundRole).name() == FLAG_BOUND_COLOR
    assert _data(filled, 1, "ID", Qt.BackgroundRole).name() == FLAG_NAN_COLOR
    assert _data(filled, 1, "z", Qt.BackgroundRole).getRgb() == (245, 161, 0, 50)
    assert _data(filled, 2, "z", Qt.BackgroundRole) is None


def test_tooltip(filled):
    tip = _data(filled, 1, "b", Qt.ToolTipRole)
    assert "<b>Status:</b> Linked" in tip and "<b>Target:</b> CIV (z=2.00000)" in tip
    assert "⚠" not in tip
    assert "Frozen" in _data(filled, 2, "logN", Qt.ToolTipRole)
    assert "Free" in _data(filled, 0, "z", Qt.ToolTipRole)
    assert "lower bound of the fit (1 km/s)" in _data(filled, 0, "b", Qt.ToolTipRole)
    assert "could not be computed" in _data(filled, 1, "ID", Qt.ToolTipRole)


@pytest.mark.parametrize("role", [Qt.DecorationRole, Qt.CheckStateRole,
                                  Qt.SizeHintRole, Qt.StatusTipRole,
                                  Qt.WhatsThisRole, Qt.AccessibleTextRole])
def test_other_roles_are_empty(filled, role):
    from astrocook.gui.inspector_table import COLUMNS
    for row in range(3):
        for col in COLUMNS:
            assert _data(filled, row, col, role) is None


def test_invalid_index(filled):
    from PySide6.QtCore import QModelIndex
    assert filled.data(QModelIndex(), Qt.DisplayRole) is None


# --- Editing a cell -------------------------------------------------------------------

def test_editing_a_cell_asks_for_the_update(filled):
    from astrocook.gui.inspector_table import COLUMNS
    requests = []
    filled.data_changed_request.connect(lambda name, params: requests.append((name, params)))
    comp = filled.get_component_at(0)
    index = filled.index(0, COLUMNS.index("b"))
    assert filled.setData(index, "12.5", Qt.EditRole)
    assert requests == [("update_component", {"b": 12.5, "uuid": comp.uuid})]
    # Not a number, or not the edit role: nothing asked
    assert not filled.setData(index, "abc", Qt.EditRole)
    assert not filled.setData(index, "13", Qt.DisplayRole)
    assert len(requests) == 1
    assert filled.setData(filled.index(0, COLUMNS.index("Trans.")), "SiIV", Qt.EditRole)
    assert requests[-1] == ("update_component", {"series": "SiIV", "uuid": comp.uuid})


def test_tooltip_of_a_link_by_expression(model):
    c1, c2, c3 = _comps()
    # c2.b = c1.b * 0.5, written with the uuid of c1
    cons = {c2.uuid: {"b": ParameterConstraintV2(
        is_free=False, target_uuid=None, expression=f"p['{c1.uuid}'].b * 0.5")}}
    model.update_data([c1, c2, c3], cons)
    tip = _data(model, 1, "b", Qt.ToolTipRole)
    assert "<b>Status:</b> Linked" in tip and "<b>Target:</b> Unknown" in tip
    assert "<code>CIV.b * 0.5</code>" in tip


# --- The proxy: filters and sorting -------------------------------------------------------

@pytest.fixture
def proxy(filled):
    from astrocook.gui.inspector_table import SystemSortFilterProxyModel
    p = SystemSortFilterProxyModel()
    p.setSourceModel(filled)
    return p


def _ids(proxy):
    from astrocook.gui.inspector_table import COLUMNS
    return [proxy.data(proxy.index(r, COLUMNS.index("ID")), Qt.UserRole)
            for r in range(proxy.rowCount())]


def test_proxy_filters(filled, proxy):
    c1, c2, c3 = (filled.get_component_at(r) for r in range(3))
    assert _ids(proxy) == [1.0, 2.0, 3.0]
    proxy.set_flagged_only(True)
    assert _ids(proxy) == [1.0, 2.0]  # the components with flags
    proxy.set_flagged_only(False)
    proxy.set_allowed_uuids({c3.uuid})
    assert _ids(proxy) == [1.0, 2.0, 3.0]  # the group filter is off
    proxy.set_group_filter_enabled(True)
    assert _ids(proxy) == [3.0]
    proxy.set_allowed_uuids({c1.uuid, c3.uuid})
    assert _ids(proxy) == [1.0, 3.0]


def test_proxy_sorts_numbers_and_text(filled, proxy):
    from astrocook.gui.inspector_table import COLUMNS
    proxy.sort(COLUMNS.index("z"), Qt.DescendingOrder)
    assert _ids(proxy) == [3.0, 2.0, 1.0]
    proxy.sort(COLUMNS.index("Trans."), Qt.DescendingOrder)  # by text
    assert _ids(proxy)[0] == 3.0  # Ly_a after CIV


def test_proxy_drops_a_row_without_component(filled, proxy):
    assert not proxy.filterAcceptsRow(99, QModelIndex())
