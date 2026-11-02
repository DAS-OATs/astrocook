import astropy.units as au
import json
import logging
import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QDialog, QDialogButtonBox, QLabel, QTableWidget, QTableWidgetItem,
    QTextEdit, QVBoxLayout
)
from typing import Dict, List, Optional, Tuple, TYPE_CHECKING

from astrocook.core.atomic_data import STANDARD_MULTIPLETS, trans_parse, xem_d
if TYPE_CHECKING:
    from astrocook.core.spectrum import SpectrumV2


#: The columns of the table of the identifications
COLUMNS = ["Regions", "Identification", "Transitions", "z", "Score", "\u03bb from (nm)",
           "\u03bb to (nm)", "Other candidates"]

#: The regions identified as members of the same multiplet (e.g. CIV_1548
#: and CIV_1550) are one row when their redshifts differ by less than this
#: (km/s): measured at the deepest pixel of each region, the members of one
#: system differ by a few km/s (pixels, blends). Each transition at most once
#: per row, so that two close systems are not merged
MULTIPLET_DV_KMS = 20.0
C_KMS = 299792.458


def multiplet_of(series: str) -> Optional[str]:
    """The multiplet of a transition (e.g. 'CIV' for 'CIV_1548'), the name
    itself for a multiplet, None for a line of its own (e.g. 'Ly_a')."""
    if series in STANDARD_MULTIPLETS:
        return series
    for name, members in STANDARD_MULTIPLETS.items():
        if series in members:
            return name
    return None


def group_multiplets(rows: List[dict]) -> List[dict]:
    """The rows of the members of one multiplet at one redshift merged into
    one (Guido: the members are one identification, and seen together they
    are easier to confirm): its regions, the transitions found, the mean
    redshift and score, the range of all the regions."""
    groups = []  # [multiplet, [rows]]
    singles = []
    for row in sorted(rows, key=lambda r: (np.nan_to_num(r["z"], nan=np.inf), r["region"])):
        name = multiplet_of(row["series"])
        if name is None or not np.isfinite(row["z"]):
            singles.append(row)
            continue
        for group_name, members in groups:
            if group_name != name or any(m["series"] == row["series"] for m in members):
                continue
            z0 = members[0]["z"]
            if abs(row["z"] - z0) * C_KMS / (1 + z0) <= MULTIPLET_DV_KMS:
                members.append(row)
                break
        else:
            groups.append([name, [row]])
    merged = [dict(row, regions=[row["region"]], transitions=row["series"]) for row in singles]
    for name, members in groups:
        if len(members) == 1:
            row = members[0]
            merged.append(dict(row, regions=[row["region"]], transitions=row["series"]))
            continue
        found = {m["series"] for m in members}
        others = []
        for m in members:
            others += [o for o in m["others"].split(", ") if o and o.split(" (")[0] not in found]
        merged.append({
            "region": min(m["region"] for m in members),
            "regions": sorted(m["region"] for m in members),
            "series": name,
            "transitions": ", ".join(sorted(found)),
            "z": float(np.mean([m["z"] for m in members])),
            "score": float(np.mean([m["score"] for m in members])),
            "xmin": float(np.nanmin([m["xmin"] for m in members])),
            "xmax": float(np.nanmax([m["xmax"] for m in members])),
            "others": ", ".join(dict.fromkeys(others)),
        })
    return sorted(merged, key=lambda r: r["region"])


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
        self.layout.addWidget(self.summary_label)

        # 2. The table of the regions
        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setStretchLastSection(True)
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
        self.resize(760, 480)

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
                logging.warning(f"Could not read 'abs_ids': {e}")

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

            rows = group_multiplets(self._rows(spec, ident_dict, region_map))
            self.summary_label.setText(
                f"Found <b>{len(ident_dict)}</b> identified absorption regions; "
                "the members of a multiplet are one row. "
                "Click a row to show its regions in the main plot.")
            self._fill_table(rows)

        except Exception as e:
            logging.error(f"Failed to parse identifications JSON: {e}")
            self.summary_label.setText("Error parsing metadata.")
            self._show_text(f"# Error:\n{e}\n\n# Raw Data:\n{json_string}")

    def _rows(self, spec, ident_dict, region_map) -> List[dict]:
        """One row per region: its best identification and its redshift
        (at the deepest pixel of the region), its score, its range and the
        other candidates."""
        rows = []
        x_nm = spec.x.to_value(au.nm)
        for region_id, candidates in sorted(ident_dict.items()):
            if not candidates:
                continue
            series_name, score = candidates[0]
            row = {"region": region_id, "series": series_name, "score": float(score),
                   "z": np.nan, "xmin": np.nan, "xmax": np.nan,
                   "others": ", ".join(f"{name} ({float(sc):.2f})" for name, sc in candidates[1:])}
            if region_map is not None:
                indices = np.where(region_map == region_id)[0]
                if indices.size:
                    row["xmin"], row["xmax"] = float(x_nm[indices[0]]), float(x_nm[indices[-1]])
                    ref_trans = self._get_ref_transition(series_name)
                    if ref_trans and ref_trans in xem_d:
                        center = indices[np.argmin(spec.y.value[indices])]
                        row["z"] = x_nm[center] / xem_d[ref_trans].to_value(au.nm) - 1.0
            rows.append(row)
        return rows

    def _fill_table(self, rows: List[dict]):
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))

        def number(value, fmt):
            return _NumberItem("" if not np.isfinite(value) else fmt.format(value),
                               value if np.isfinite(value) else -np.inf)

        for r, row in enumerate(rows):
            regions = _NumberItem(", ".join(str(r) for r in row["regions"]), row["region"])
            items = [regions, QTableWidgetItem(row["series"]),
                     QTableWidgetItem(row["transitions"]),
                     number(row["z"], "{:.5f}"), number(row["score"], "{:.2f}"),
                     number(row["xmin"], "{:.3f}"), number(row["xmax"], "{:.3f}"),
                     QTableWidgetItem(row["others"])]
            for c, item in enumerate(items):
                if c < len(items) - 1:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.table.setItem(r, c, item)
            self._regions[row["region"]] = (row["xmin"], row["xmax"])
        self.table.setSortingEnabled(True)
        self.table.sortItems(0, Qt.SortOrder.AscendingOrder)
        self.table.resizeColumnsToContents()

    def _on_row_selected(self):
        """The region of the selected row in the main plot, with a margin of
        half its width on each side."""
        rows = self.table.selectionModel().selectedRows()
        if not rows or self.zoom_to is None:
            return
        region = self.table.item(rows[0].row(), 0).value
        xmin, xmax = self._regions.get(region, (np.nan, np.nan))
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

    def _invert_identifications(self, 
                                spec: 'SpectrumV2', 
                                ident_dict: Dict[int, List[Tuple[str, float]]], 
                                region_map: np.ndarray) -> Dict[str, set]:
        """
        Inverts the region-based dict to a transition-based dict of redshifts.
        {1: [('CIV', 0.9)]}  =>  {'CIV': {2.12345}}
        """
        inverted_map = {}
        if region_map is None:
             logging.warning("Cannot invert identifications, region_map is None.")
             return {}
             
        x_nm = spec.x.to_value(au.nm)
        
        for region_id, id_list in ident_dict.items():
            if not id_list:
                continue
                
            # Use the *top* identification for this region
            series_name, score = id_list[0]
            
            # Find the center of this region to get its redshift
            try:
                region_indices = np.where(region_map == region_id)[0]
                if not region_indices.size:
                    continue
                
                # Get the reference transition wavelength
                ref_trans = self._get_ref_transition(series_name)
                if not ref_trans or ref_trans not in xem_d:
                    logging.warning(f"Could not find ref_trans for '{series_name}'")
                    continue
                
                xem_nm = xem_d[ref_trans].to_value(au.nm)
                
                # Find the pixel with max absorption (min flux) in this region
                region_flux = spec.y.value[region_indices]
                center_idx = region_indices[np.argmin(region_flux)]
                
                z_peak = (x_nm[center_idx] / xem_nm) - 1.0
                
                if series_name not in inverted_map:
                    inverted_map[series_name] = set()
                inverted_map[series_name].add(z_peak)

            except Exception as e:
                logging.warning(f"Could not invert region {region_id} for {series_name}: {e}")
                
        return inverted_map