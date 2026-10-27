"""The table of the System Inspector: what the model gives Qt for each role,
and the columns fitted to their contents.

Characterization tests written before the table was made faster (phase 5:
the columns were fitted by Qt with ``ResizeToContents``, which queried every
row several times): the values of the model and the fit of the columns must
not change.
"""
import pytest
from PySide6.QtCore import Qt

from astrocook.core.structures import ComponentDataV2, ParameterConstraintV2


@pytest.fixture
def model(qapp):
    from astrocook.gui.system_inspector import SystemTableModel
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
    from astrocook.gui.system_inspector import COLUMNS
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
    assert _data(filled, 0, "z", Qt.DisplayRole) == "2.000000 ± 1.50e-05"
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
    from astrocook.gui.system_inspector import FLAG_BOUND_COLOR, FLAG_NAN_COLOR
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
    from astrocook.gui.system_inspector import COLUMNS
    for row in range(3):
        for col in COLUMNS:
            assert _data(filled, row, col, role) is None


def test_invalid_index(filled):
    from PySide6.QtCore import QModelIndex
    assert filled.data(QModelIndex(), Qt.DisplayRole) is None
