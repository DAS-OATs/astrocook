"""The Identifications window, a table of the identified regions (visual
review of phase 5-bis: it was a list of redshifts per transition), and the
formats of the Data Inspector."""
import json

import numpy as np
import pytest
from PySide6.QtCore import Qt

from astrocook.gui.data_inspector import format_value
from astrocook.gui.identification_viewer_dialog import (COLUMNS, IdentificationViewerDialog,
                                                        group_multiplets)

from test_acs2_round_trip import rich_session  # noqa: F401 (fixture)


@pytest.fixture
def dialog(qtbot, rich_session):  # noqa: F811
    meta = rich_session.spec._data.meta
    meta["region_identifications"] = json.dumps(
        {"1": [["CIV", 0.9], ["SiIV", 0.4]], "2": [["MgII", 0.7]], "3": []})
    zooms = []
    dlg = IdentificationViewerDialog(rich_session.spec, "rich",
                                     zoom_to=lambda a, b: zooms.append((a, b)))
    qtbot.addWidget(dlg)
    dlg.zooms = zooms
    return dlg


def _cell(dlg, row, name):
    return dlg.table.item(row, COLUMNS.index(name)).text()


def test_one_row_per_identified_region(dialog, rich_session):  # noqa: F811
    table = dialog.table
    assert table.isVisible() or not dialog.results_edit.isVisible()
    assert [table.horizontalHeaderItem(c).text() for c in range(table.columnCount())] == COLUMNS
    assert table.rowCount() == 2  # region 3 has no candidates
    assert "Found <b>3</b>" in dialog.summary_label.text()
    assert _cell(dialog, 0, "Regions") == "1" and _cell(dialog, 0, "Identification") == "CIV"
    assert _cell(dialog, 0, "Score") == "0.90"
    assert _cell(dialog, 0, "Other candidates") == "SiIV (0.40)"
    # The range of the region (abs_ids: 100 pixels each) and its redshift
    x = rich_session.spec.x.value
    assert float(_cell(dialog, 0, "λ from (nm)")) == pytest.approx(x[100], abs=1e-3)
    assert float(_cell(dialog, 0, "λ to (nm)")) == pytest.approx(x[199], abs=1e-3)
    assert _cell(dialog, 0, "z") != ""
    assert _cell(dialog, 1, "Identification") == "MgII"


def _row(region, series, z, score=0.9, xmin=500.0, xmax=500.1, others=""):
    return {"region": region, "series": series, "z": z, "score": score,
            "xmin": xmin, "xmax": xmax, "others": others}


def test_members_of_a_multiplet_are_one_row():
    # The members of a doublet at one redshift: one identification, seen
    # together (Guido); a line of its own stays alone
    z = 2.0
    dz = 5.0 * (1 + z) / 299792.458  # 5 km/s
    rows = group_multiplets([
        _row(6, "CIV_1548", z, 0.91, 464.4, 464.5, "SiIV_1393 (0.30)"),
        _row(8, "CIV_1550", z + dz, 0.89, 465.2, 465.3, "CIV_1548 (0.20)"),
        _row(9, "Ly_a", 2.5, 0.7)])
    assert len(rows) == 2
    civ = rows[0]
    assert civ["regions"] == [6, 8] and civ["series"] == "CIV"
    assert civ["transitions"] == "CIV_1548, CIV_1550"
    assert civ["score"] == pytest.approx(0.90)
    assert (civ["xmin"], civ["xmax"]) == (464.4, 465.3)  # both regions
    assert civ["others"] == "SiIV_1393 (0.30)"  # not its own members
    assert rows[1]["regions"] == [9] and rows[1]["transitions"] == "Ly_a"


def test_close_systems_are_not_merged():
    # Each transition once per row: two CIV systems 10 km/s apart are two
    # rows; members farther than MULTIPLET_DV_KMS are two rows too
    z = 2.0
    dv = lambda kms: kms * (1 + z) / 299792.458  # noqa: E731
    rows = group_multiplets([_row(1, "CIV_1548", z), _row(2, "CIV_1548", z + dv(10)),
                             _row(3, "CIV_1550", z + dv(2)), _row(4, "CIV_1550", z + dv(12))])
    assert sorted(r["regions"] for r in rows) == [[1, 3], [2, 4]]
    rows = group_multiplets([_row(1, "CIV_1548", z), _row(2, "CIV_1550", z + dv(100))])
    assert len(rows) == 2


def test_the_table_groups_the_members(qtbot, rich_session):  # noqa: F811
    x = rich_session.spec.x.value
    meta = rich_session.spec._data.meta
    meta["region_identifications"] = json.dumps(
        {"1": [["CIV_1548", 0.9]], "2": [["CIV_1550", 0.9]]})
    zooms = []
    dlg = IdentificationViewerDialog(rich_session.spec, "rich",
                                     zoom_to=lambda a, b: zooms.append((a, b)))
    qtbot.addWidget(dlg)
    # The regions are at their own redshifts in this synthetic spectrum:
    # grouped only if close, else one row each
    assert dlg.table.rowCount() in (1, 2)
    dlg.table.selectRow(0)
    assert zooms and zooms[0][0] < x[100]


def test_sorted_by_number(dialog):
    table = dialog.table
    table.sortItems(COLUMNS.index("Score"), Qt.SortOrder.DescendingOrder)
    assert [table.item(r, 0).text() for r in range(2)] == ["1", "2"]
    table.sortItems(COLUMNS.index("Score"), Qt.SortOrder.AscendingOrder)
    assert [table.item(r, 0).text() for r in range(2)] == ["2", "1"]


def test_a_row_shows_its_region(dialog, rich_session):  # noqa: F811
    dialog.table.selectRow(0)
    assert len(dialog.zooms) == 1
    xmin, xmax = dialog.zooms[0]
    x = rich_session.spec.x.value
    assert xmin < x[100] and xmax > x[199]


def test_no_identifications(qtbot, rich_session):  # noqa: F811
    rich_session.spec._data.meta.pop("region_identifications")
    dlg = IdentificationViewerDialog(rich_session.spec, "rich")
    qtbot.addWidget(dlg)
    assert "No identifications" in dlg.summary_label.text()
    assert "Identify Absorption Lines" in dlg.results_edit.toPlainText()


class _Column:
    def __init__(self, values):
        self.value = np.asarray(values)


@pytest.mark.parametrize("col, value, text", [
    ("x", 544.40512345, "544.40512"),          # wavelengths: 5 decimals
    ("y", 34521.34, "34521"),                  # 5 significant digits
    ("y", 0.872314, "0.87231"),
    ("y", -12.5, "-12.500"),
    ("y", 1.2345e-17, "1.2345e-17"),           # a flux in cgs: scientific
    ("resol", 135000.0, "135000"),
    ("y", 0.0, "0"),
    ("y", np.nan, "nan"),
])
def test_data_inspector_formats(col, value, text):
    assert format_value(col, value) == text


def test_data_inspector_integer_columns():
    assert format_value("abs_mask", np.True_, _Column([True])) == "1"
    assert format_value("abs_ids", np.int64(12), _Column([12])) == "12"
