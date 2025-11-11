import json
import logging
import numpy as np
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QTextEdit, QDialogButtonBox, QLabel
)
from typing import TYPE_CHECKING, Optional

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
        if not json_string:
            self.summary_label.setText("No identifications found in session metadata.")
            self.results_edit.setPlainText("# Run 'Absorbers > Identify Absorption Lines' first.")
            return

        try:
            # 1. De-serialize the JSON string
            ident_dict_raw = json.loads(json_string)
            # Convert string keys back to int
            ident_dict = {int(k): v for k, v in ident_dict_raw.items()}
            
            num_regions_identified = len(ident_dict)
            
            # 2. Get total number of regions
            total_regions = 0
            if spec.has_aux_column('region_id'):
                region_map = spec.get_column('region_id').value
                total_regions = int(np.max(region_map))
                
            # 3. Update summary
            self.summary_label.setText(
                f"Found identifications for <b>{num_regions_identified}</b> of "
                f"<b>{total_regions}</b> total absorption regions."
            )
            
            # 4. Pretty-print the dictionary
            pretty_json = json.dumps(ident_dict, indent=2, sort_keys=True)
            self.results_edit.setPlainText(pretty_json)

        except Exception as e:
            logging.error(f"Failed to parse identifications JSON: {e}")
            self.summary_label.setText("Error parsing metadata.")
            self.results_edit.setPlainText(f"# Error:\n{e}\n\n# Raw Data:\n{json_string}")