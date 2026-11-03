import astropy.units as au
import json
import logging
import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QDialog, QDialogButtonBox, QHeaderView, QLabel, QTableWidget,
    QTableWidgetItem,
    QTextEdit, QVBoxLayout
)
from typing import List, Optional, TYPE_CHECKING

from astrocook.core.atomic_data import STANDARD_MULTIPLETS, trans_parse, xem_d
if TYPE_CHECKING:
    from astrocook.core.spectrum import SpectrumV2

logger = logging.getLogger(__name__)


#: The columns of the table of the identifications
COLUMNS = ["Regions", "Identification", "z (mean)", "z min", "z max", "Score", "Alternatives"]

#: The label of a region of the forest that no doublet identified: the
#: recipe assigns it to Ly_a with score 0, it is not an identification
ASSUMED_LYA = "Ly_a (assumed)"

C_KMS = 299792.458
#: Tolerance (km/s) on the overlap of the redshift ranges of the two regions
#: of a doublet: the recipe finds the pairs on a shifted map of the regions,
#: interpolated to the nearest pixel
PAIR_TOLERANCE_KMS = 10.0


def _z_range(x_nm, region_map, region, transition):
    """The redshifts of the first and last pixel of a region, for a
    transition."""
    indices = np.where(region_map == region)[0]
    if transition in STANDARD_MULTIPLETS:  # a multiplet (older sessions): its first line
        transition = STANDARD_MULTIPLETS[transition][0]
    if not indices.size or transition not in xem_d:
        return None
    rest = xem_d[transition].to_value(au.nm)
    return x_nm[indices[0]] / rest - 1.0, x_nm[indices[-1]] / rest - 1.0


def doublet_rows(spec, ident_dict, region_map) -> List[dict]:
    """The identifications as the recipe makes them (Guido): one row per
    doublet found, i.e. per pair of regions whose redshift ranges overlap
    for the two lines of the doublet, with the same score (the score of
    the pair, given to both regions); the redshift range is the envelope
    of the two regions, the redshift its middle. A region can be in more
    than one pair. The regions of the forest that no doublet identified
    are Ly_a with score 0 in the recipe: shown as assumed, without a score.
    The pairs are found again from the candidates stored per region (the
    recipe does not store the pairs)."""
    x_nm = spec.x.to_value(au.nm)
    rows, paired = [], set()
    for name, members in STANDARD_MULTIPLETS.items():
        if len(members) != 2:
            continue
        first, second = members
        with_first = [(r, sc) for r, cands in ident_dict.items() for n, sc in cands if n == first]
        with_second = [(r, sc) for r, cands in ident_dict.items() for n, sc in cands if n == second]
        for r1, s1 in with_first:
            z1 = _z_range(x_nm, region_map, r1, first)
            for r2, s2 in with_second:
                if r2 == r1 or abs(float(s1) - float(s2)) > 1e-9:
                    continue
                z2 = _z_range(x_nm, region_map, r2, second)
                if z1 is None or z2 is None:
                    continue
                tol = PAIR_TOLERANCE_KMS / C_KMS * (1.0 + z1[0])
                if min(z1[1], z2[1]) < max(z1[0], z2[0]) - tol:
                    continue
                zmin, zmax = min(z1[0], z2[0]), max(z1[1], z2[1])
                indices = np.where((region_map == r1) | (region_map == r2))[0]
                rows.append({"region": min(r1, r2), "regions": sorted((r1, r2)),
                             "series": name, "score": float(s1),
                             "z": 0.5 * (zmin + zmax), "zmin": zmin, "zmax": zmax,
                             "xmin": float(x_nm[indices[0]]), "xmax": float(x_nm[indices[-1]])})
                paired.update(((r1, first), (r2, second)))
    for region, cands in ident_dict.items():
        for n, sc in cands:
            if (region, n) in paired:
                continue
            assumed = n == "Ly_a" and float(sc) == 0.0
            zr = _z_range(x_nm, region_map, region, n)
            indices = np.where(region_map == region)[0]
            if zr is None or not indices.size:
                continue
            rows.append({"region": region, "regions": [region],
                         "series": ASSUMED_LYA if assumed else n,
                         "score": np.nan if assumed else float(sc),
                         "z": 0.5 * (zr[0] + zr[1]), "zmin": zr[0], "zmax": zr[1],
                         "xmin": float(x_nm[indices[0]]), "xmax": float(x_nm[indices[-1]])})
    rows.sort(key=lambda r: (r["region"], r["series"]))
    # The alternatives of a row: the other rows with one of its regions (the
    # same region in another pair, of another doublet or of the same one at
    # another redshift), as the recipe kept them for that region
    for row in rows:
        others = [o for o in rows if o is not row and set(o["regions"]) & set(row["regions"])
                  and not np.isnan(o["score"])]
        row["others"] = ", ".join(
            f"{o['series']} {o['z']:.5f} [{', '.join(map(str, o['regions']))}] ({o['score']:.2f})"
                                  for o in sorted(others, key=lambda o: -o["score"]))
    return rows


class _NumberItem(QTableWidgetItem):
    """A cell sorted by its number, not by its text."""

    def __init__(self, text: str, value: float):
        super().__init__(text)
        self.value = value

    def __lt__(self, other):
        if isinstance(other, _NumberItem):
            return self.value < other.value
        return super().__lt__(other)


class IdentificationViewerDialog(QDialog):
    """
    A non-modal, read-only dialog with the results of the identification of
    the absorbed regions (identify_lines and the automatic pipelines): one
    row per region, sortable; a click on a row shows the region in the main
    plot (``zoom_to(xmin, xmax)``, in nm). It was a list of redshifts per
    transition (visual review of phase 5-bis).
    """
    def __init__(self, spec: Optional['SpectrumV2'], session_name: str, parent=None,
                 zoom_to=None):
        super().__init__(parent)
        self.setWindowTitle(f"Identifications: {session_name}")
        self.zoom_to = zoom_to

        self.layout = QVBoxLayout(self)

        # 1. Summary Label
        self.summary_label = QLabel("Parsing identification results...")
        self.summary_label.setWordWrap(True)
        self.layout.addWidget(self.summary_label)

        # 2. The table of the regions
        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setStretchLastSection(False)
        self.table.itemSelectionChanged.connect(self._on_row_selected)
        self.layout.addWidget(self.table)
        # Kept for the errors and the hint when there is nothing to show
        self.results_edit = QTextEdit()
        self.results_edit.setReadOnly(True)
        self.results_edit.setVisible(False)
        self.layout.addWidget(self.results_edit)

        # 3. Close Button
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        button_box.rejected.connect(self.reject)
        self.layout.addWidget(button_box)

        self._regions = {}  # row region ID -> (xmin, xmax) in nm
        self._populate_results(spec)
        self.resize(760, 520)

    def _show_text(self, text: str):
        self.table.setVisible(False)
        self.results_edit.setVisible(True)
        self.results_edit.setPlainText(text)

    def _populate_results(self, spec: Optional['SpectrumV2']):
        """Parses spec.meta and fills the table."""
        if spec is None:
            self.summary_label.setText("No spectrum loaded.")
            self._show_text("")
            return

        json_string = spec.meta.get('region_identifications')

        region_map = None
        if spec.has_aux_column('abs_ids'):
            try:
                region_map = spec.get_column('abs_ids').value
            except Exception as e:
                logger.warning(f"Could not read 'abs_ids': {e}")

        if not json_string:
            self.summary_label.setText("No identifications found in session metadata.")
            self._show_text("# Run 'Absorbers > Identify Absorption Lines' first.")
            return

        try:
            # 1. De-serialize: a JSON string as written by identify_lines, or
            # already a dict after an .acs2 round trip (findings register F27)
            if isinstance(json_string, str):
                ident_dict_raw = json.loads(json_string)
            else:
                ident_dict_raw = json_string
            ident_dict = {int(float(k)): v for k, v in ident_dict_raw.items()}

            rows = doublet_rows(spec, ident_dict, region_map) if region_map is not None else []
            pairs = sum(1 for r in rows if len(r["regions"]) == 2)
            self.summary_label.setText(
                f"<b>{pairs}</b> doublets found, in {len(ident_dict)} absorbed regions; "
                "the other regions of the forest are assumed Ly_a. The redshift is the "
                "middle of the range of the regions. Click a row to show its regions in "
                "the main plot.")
            self._fill_table(rows)

        except Exception as e:
            logger.error(f"Failed to parse identifications JSON: {e}")
            self.summary_label.setText("Error parsing metadata.")
            self._show_text(f"# Error:\n{e}\n\n# Raw Data:\n{json_string}")

    def _fill_table(self, rows: List[dict]):
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))

        def number(value, fmt):
            return _NumberItem("" if not np.isfinite(value) else fmt.format(value),
                               value if np.isfinite(value) else -np.inf)

        for r, row in enumerate(rows):
            regions = _NumberItem(", ".join(str(region) for region in row["regions"]),
                                  row["region"])
            items = [regions, QTableWidgetItem(row["series"]),
                     number(row["z"], "{:.5f}"), number(row["zmin"], "{:.5f}"),
                     number(row["zmax"], "{:.5f}"), number(row["score"], "{:.2f}"),
                     QTableWidgetItem(row.get("others", ""))]
            for c, item in enumerate(items):
                if c < len(items) - 1:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.table.setItem(r, c, item)
            self._regions[r] = (row["xmin"], row["xmax"])
            regions.setData(Qt.ItemDataRole.UserRole, r)
        self.table.setSortingEnabled(True)
        self.table.sortItems(0, Qt.SortOrder.AscendingOrder)
        # The columns as wide as their contents, the last one (the
        # alternatives) taking the rest of the window: no empty column on
        # the right (reported by Guido)
        header = self.table.horizontalHeader()
        for c in range(len(COLUMNS) - 1):
            header.setSectionResizeMode(c, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(len(COLUMNS) - 1, QHeaderView.ResizeMode.Stretch)

    def _on_row_selected(self):
        """The region of the selected row in the main plot, with a margin of
        half its width on each side."""
        rows = self.table.selectionModel().selectedRows()
        if not rows or self.zoom_to is None:
            return
        key = self.table.item(rows[0].row(), 0).data(Qt.ItemDataRole.UserRole)
        xmin, xmax = self._regions.get(key, (np.nan, np.nan))
        if np.isfinite(xmin) and np.isfinite(xmax):
            pad = max(0.5 * (xmax - xmin), 0.05)
            self.zoom_to(xmin - pad, xmax + pad)

    def _get_ref_transition(self, series_name: str) -> Optional[str]:
        """Helper to get the primary transition for a series."""
        if series_name in STANDARD_MULTIPLETS:
            return STANDARD_MULTIPLETS[series_name][0]
        if series_name in xem_d:
            return series_name
        try:
            return trans_parse(series_name)[0] # V1 series syntax
        except Exception:
            pass
        return None
