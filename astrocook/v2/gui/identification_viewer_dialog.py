import astropy.units as au
import json
import logging
import numpy as np
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QTextEdit, QDialogButtonBox, QLabel
)
from typing import Dict, List, Optional, Tuple, TYPE_CHECKING

from ..atomic_data import STANDARD_MULTIPLETS, xem_d
try:
    from ...v1.functions import trans_parse
    V1_FUNCTIONS_AVAILABLE = True
except ImportError:
    V1_FUNCTIONS_AVAILABLE = False
if TYPE_CHECKING:
    from ..spectrum import SpectrumV2


class IdentificationViewerDialog(QDialog):
    """
    A non-modal, read-only dialog to display the results of the
    identify_lines recipe.
    """
    def __init__(self, spec: Optional['SpectrumV2'], session_name: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Identifications: {session_name}")

        self.layout = QVBoxLayout(self)

        # 1. Summary Label
        self.summary_label = QLabel("Parsing identification results...")
        self.layout.addWidget(self.summary_label)

        # 2. Results Text Edit
        self.results_edit = QTextEdit()
        self.results_edit.setReadOnly(True)
        self.layout.addWidget(self.results_edit)

        # 3. Close Button
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        button_box.rejected.connect(self.reject)
        self.layout.addWidget(button_box)

        self._populate_results(spec)
        self.resize(500, 400)

    def _populate_results(self, spec: Optional['SpectrumV2']):
        """
        Parses the spec.meta and populates the text edit.
        """
        if spec is None:
            self.summary_label.setText("No spectrum loaded.")
            self.results_edit.setPlainText("")
            return

        json_string = spec.meta.get('region_identifications')

        total_regions = 0
        region_map = None
        if spec.has_aux_column('abs_ids'): # <<< *** Use new column name 'abs_ids' ***
            try:
                region_map = spec.get_column('abs_ids').value
                # Find all unique non-zero IDs
                total_regions = len(np.unique(region_map[region_map > 0]))
            except Exception as e:
                logging.warning(f"Could not count regions from 'abs_ids': {e}")

        if not json_string:
            self.summary_label.setText("No identifications found in session metadata.")
            self.results_edit.setPlainText("# Run 'Absorbers > Identify Absorption Lines' first.")
            return

        try:
            # 1. De-serialize the JSON string
            ident_dict_raw = json.loads(json_string)
            ident_dict = {int(k): v for k, v in ident_dict_raw.items()}
            num_regions_identified = len(ident_dict)

            # 2. Update summary label
            self.summary_label.setText(
                f"Found <b>{num_regions_identified}</b> identified absorption regions."
            )
            
            # --- Invert dictionary for transition-based list *** ---
            inverted_summary = self._invert_identifications(spec, ident_dict, region_map)
            
            # 4. Format the final summary string
            summary_lines = []
            # Sort by series name
            for series_name, z_list in sorted(inverted_summary.items()):
                # Format redshifts to 5 decimal places
                z_str = ", ".join([f"{z:.5f}" for z in sorted(list(z_list))])
                summary_lines.append(f"<b>{series_name}</b> ({len(z_list)}):<br>{z_str}")
                
            self.results_edit.setHtml("<br><br>".join(summary_lines))

        except Exception as e:
            logging.error(f"Failed to parse identifications JSON: {e}")
            self.summary_label.setText("Error parsing metadata.")
            self.results_edit.setPlainText(f"# Error:\n{e}\n\n# Raw Data:\n{json_string}")
    
    def _get_ref_transition(self, series_name: str) -> Optional[str]:
        """Helper to get the primary transition for a series."""
        if series_name in STANDARD_MULTIPLETS:
            return STANDARD_MULTIPLETS[series_name][0]
        if series_name in xem_d:
            return series_name
        if V1_FUNCTIONS_AVAILABLE:
            try:
                return trans_parse(series_name)[0] # V1 fallback
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