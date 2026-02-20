import logging
import os
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2


FILE_RECIPES_SCHEMAS = {
    "export_ascii": {
        "brief": "Export session data to ASCII.",
        "details": "Saves selected session components (Flux, Continuum, System List) as CSV files.",
        "params": [
            {"name": "targets", "type": str, "default": "y, cont, systems", "doc": "Comma-separated elements to export"},
            {"name": "path", "type": str, "default": "", "gui_hidden": True, "doc": "Destination folder"}
        ],
        "url": "file_cb.html#export"
    },
    "import_ascii": {
        "brief": "Import session data from ASCII.",
        "details": "Imports session components (Flux, Continuum, System List) from CSV files.",
        "params": [
            {"name": "path", "type": str, "default": "", "gui_hidden": True, "doc": "Source folder"}
        ],
        "url": "file_cb.html#export"
    }
}

class RecipeFileV2:
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cb'

    def export_ascii(self, targets: str, path: str, prefix: str) -> bool:
        target_list = [t.strip() for t in targets.split(',') if t.strip()]
        saved_any = False
        
        # Get authoritative column lists
        spec_colnames = self._session.spec.t.colnames if self._session.spec else []
        syst_colnames = self._session.systs.t.colnames if self._session.systs else []
        
        try:
            # 1. Spectrum Export
            # Any target that exists in spec table (and isn't the literal 'systems')
            spec_req = [t for t in target_list if t in spec_colnames and t != 'x']
            if spec_req:
                valid_cols = ['x'] + spec_req
                raw_data = {c: self._session.spec.t[c] for c in valid_cols}
                out_path = os.path.join(path, f"{prefix}_spec.csv")
                ascii.write(Table(raw_data), out_path, format='csv', overwrite=True)
                saved_any = True

            # 2. Systems Export (Refined)
            syst_req = [t for t in target_list if t in syst_colnames]
            if ("systems" in target_list or syst_req) and self._session.systs:
                syst_table = self._session.systs.t
                
                if "systems" not in target_list:
                    # [NEW LOGIC] Always include 'series' and 'z' for context.
                    # We use dict.fromkeys to maintain order and uniqueness.
                    # 'series' first, then 'z', then other requested params.
                    base_cols = ['series', 'z']
                    # Add any other columns requested that aren't already z/series
                    extra_cols = [t for t in syst_req if t not in base_cols]
                    
                    cols_to_keep = base_cols + extra_cols
                    # Verify they exist in the actual table before slicing
                    cols_to_keep = [c for c in cols_to_keep if c in syst_table.colnames]
                    syst_table = syst_table[cols_to_keep]
                
                out_path = os.path.join(path, f"{prefix}_systems.csv")
                ascii.write(syst_table, out_path, format='csv', overwrite=True)
                saved_any = True
                
            return saved_any

        except Exception as e:
            logging.error(f"Export failed: {e}")
            return False
        
    def import_ascii(self, file_path: str) -> 'SessionV2':

        try:
            # 1. Load the table
            table = ascii.read(file_path)
            current_len = len(self._session.spec.x)
            import_len = len(table)

            # 2. Safety Check: Pixel count must match
            if import_len != current_len:
                raise ValueError(f"Length mismatch: File has {import_len} pixels, but current session has {current_len}.")

            new_spec = self._session.spec
            imported_cols = []

            # 3. Iterate through columns in the ASCII file
            for col_name in table.colnames:
                # We skip 'x' to prevent grid misalignment, but use everything else
                if col_name.lower() == 'x':
                    continue
                
                values = np.array(table[col_name])
                # Update the column (this handles both core y/dy and aux columns)
                new_spec = new_spec.update_column(col_name, values)
                imported_cols.append(col_name)

            logging.info(f"Imported columns from ASCII: {', '.join(imported_cols)}")
            return self._session.with_new_spectrum(new_spec)

        except Exception as e:
            logging.error(f"Import from ASCII failed: {e}")
            return 0