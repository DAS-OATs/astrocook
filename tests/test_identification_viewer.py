"""The Identifications window, a table of the identified regions (visual
review of phase 5-bis: it was a list of redshifts per transition), and the
formats of the Data Inspector."""
import json

import numpy as np
import pytest
from PySide6.QtCore import Qt

from astrocook.gui.data_inspector import format_value
from astrocook.gui.identification_viewer_dialog import COLUMNS, IdentificationViewerDialog

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
    assert _cell(dialog, 0, "Region") == "1" and _cell(dialog, 0, "Identification") == "CIV"
    assert _cell(dialog, 0, "Score") == "0.90"
    assert _cell(dialog, 0, "Other candidates") == "SiIV (0.40)"
    # The range of the region (abs_ids: 100 pixels each) and its redshift
    x = rich_session.spec.x.value
    assert float(_cell(dialog, 0, "λ from (nm)")) == pytest.approx(x[100], abs=1e-3)
    assert float(_cell(dialog, 0, "λ to (nm)")) == pytest.approx(x[199], abs=1e-3)
    assert _cell(dialog, 0, "z") != ""
    assert _cell(dialog, 1, "Identification") == "MgII"


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
