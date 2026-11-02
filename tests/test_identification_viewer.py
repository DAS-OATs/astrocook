"""The Identifications window, a table of the identified regions (visual
review of phase 5-bis: it was a list of redshifts per transition), and the
formats of the Data Inspector."""
import json

import numpy as np
import pytest
from PySide6.QtCore import Qt

from astrocook.gui.data_inspector import format_value
from astrocook.gui.identification_viewer_dialog import (ASSUMED_LYA, COLUMNS,
                                                        IdentificationViewerDialog, doublet_rows)

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


class _Spec:
    """A spectrum with only an x axis in nm, for doublet_rows."""

    def __init__(self, x_nm):
        import astropy.units as au
        self.x = x_nm * au.nm


def _doublet_case(z=2.0, z2=None, score=0.9, score2=None):
    """Regions 1 (CIV_1548) and 2 (CIV_1550) at z (and z2), 3 a forest
    region with the default Ly_a."""
    from astrocook.core.atomic_data import xem_d
    import astropy.units as au
    x = np.linspace(360.0, 480.0, 120001)
    regions = np.zeros(len(x), dtype=int)
    for region, line, zz in ((1, "CIV_1548", z), (2, "CIV_1550", z2 or z), (3, "Ly_a", 2.1)):
        centre = xem_d[line].to_value(au.nm) * (1 + zz)
        regions[np.abs(x - centre) < 0.02] = region
    ident = {1: [["CIV_1548", score]], 2: [["CIV_1550", score2 or score]], 3: [["Ly_a", 0.0]]}
    return _Spec(x), ident, regions


def test_one_row_per_doublet():
    # As the recipe identifies: the doublet of a pair of regions, its
    # redshift range the envelope of both, its redshift the middle (Guido)
    spec, ident, regions = _doublet_case()
    rows = doublet_rows(spec, ident, regions)
    assert [r["series"] for r in rows] == ["CIV", ASSUMED_LYA]
    civ = rows[0]
    assert civ["regions"] == [1, 2] and civ["score"] == pytest.approx(0.9)
    assert civ["zmin"] < 2.0 < civ["zmax"]
    assert civ["z"] == pytest.approx(0.5 * (civ["zmin"] + civ["zmax"]))
    assert civ["xmin"] < 464.5 and civ["xmax"] > 465.2  # both regions, for the zoom
    # The forest region that no doublet identified: Ly_a assumed, no score
    assert rows[1]["regions"] == [3] and np.isnan(rows[1]["score"])


def test_pairs_need_overlap_and_the_same_score():
    # Regions at redshifts that do not overlap, or with different scores
    # (two different pairs), are not a doublet
    spec, ident, regions = _doublet_case(z2=2.01)
    assert all(len(r["regions"]) == 1 for r in doublet_rows(spec, ident, regions))
    spec, ident, regions = _doublet_case(score=0.9, score2=0.8)
    assert all(len(r["regions"]) == 1 for r in doublet_rows(spec, ident, regions))


def test_a_region_in_two_pairs():
    # As region 106 of the test session: two rows, one per pair
    spec, ident, regions = _doublet_case()
    ident[1].append(["CIV_1548", 0.7])
    ident[4] = [["CIV_1550", 0.7]]
    regions[np.where(regions == 2)[0] + 5] = 4  # a second region overlapping
    rows = [r for r in doublet_rows(spec, ident, regions) if r["series"] == "CIV"]
    assert sorted(r["regions"] for r in rows) == [[1, 2], [1, 4]]
    # Each is an alternative of the other (they share region 1)
    first = next(r for r in rows if r["regions"] == [1, 2])
    assert first["others"].startswith("CIV ") and first["others"].endswith("[1, 4] (0.70)")


def test_the_table(qtbot, rich_session):  # noqa: F811
    # The table on a session, its columns and the click on a row
    meta = rich_session.spec._data.meta
    meta["region_identifications"] = json.dumps({"1": [["CIV", 0.9]], "2": [["Ly_a", 0.0]]})
    zooms = []
    dlg = IdentificationViewerDialog(rich_session.spec, "rich",
                                     zoom_to=lambda a, b: zooms.append((a, b)))
    qtbot.addWidget(dlg)
    table = dlg.table
    assert [table.horizontalHeaderItem(c).text() for c in range(table.columnCount())] == COLUMNS
    assert table.rowCount() == 2
    assert _cell(dlg, 1, "Identification") == ASSUMED_LYA and _cell(dlg, 1, "Score") == ""
    assert "doublets found" in dlg.summary_label.text()
    # The last column takes the rest of the window: no empty column on the
    # right (Guido)
    from PySide6.QtWidgets import QHeaderView
    header = table.horizontalHeader()
    assert header.sectionResizeMode(len(COLUMNS) - 1) == QHeaderView.ResizeMode.Stretch
    dlg.resize(1200, 500)
    dlg.show()
    qtbot.wait(50)
    used = sum(header.sectionSize(c) for c in range(len(COLUMNS)))
    assert used >= table.viewport().width() - 2
    table.sortItems(COLUMNS.index("Regions"), Qt.SortOrder.DescendingOrder)
    assert _cell(dlg, 0, "Regions") == "2"
    table.selectRow(1)  # region 1 after sorting
    x = rich_session.spec.x.value
    assert zooms and zooms[-1][0] < x[100] and zooms[-1][1] > x[199]


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
