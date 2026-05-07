import logging
import os
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2


FILE_RECIPES_SCHEMAS = {
    "export_ascii": {
        "brief": "Export session data to ASCII.",
        "details": "Saves selected session components (Flux, Continuum, System List) as CSV files.",
        "params": [
            {"name": "targets", "type": str, "default": "F, cont, systems", "doc": "Comma-separated elements to export"},
            {"name": "path", "type": str, "default": "", "gui_hidden": True, "doc": "Destination folder"}
        ],
        "url": "file_cb.html#export"
    },
    "import_ascii_spec": {
        "brief": "Import spectrum data from ASCII.",
        "details": "Imports spectrum columns from a CSV file. Overwrites existing columns with matching names.",
        "params": [
            {"name": "path", "type": str, "default": "", "gui_hidden": True, "doc": "Source file"}
        ],
        "url": "file_cb.html#import"
    },
    "import_ascii_systs": {
        "brief": "Import system list from ASCII.",
        "details": "Imports absorption systems from a CSV file containing 'series' and 'z' columns.",
        "params": [
            {"name": "path", "type": str, "default": "", "gui_hidden": True, "doc": "Source file"},
            {"name": "syst_mode", "type": str, "default": "append", "doc": "'append' to or 'replace' the current list."}
        ],
        "url": "file_cb.html#import"
    }
}

class RecipeFileV2:
    """
    Recipes for file input/output operations.

    This class manages exporting session data (spectra and system lists) to 
    external formats like ASCII/CSV, and importing modified data back into 
    the session. It operates immutably where applicable to support Undo/Redo.
    """
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'cb'

    def export_ascii(self, targets: str, path: str, prefix: str) -> bool:
        """
        Export selected components of the session to ASCII/CSV files.

        Generates distinct files for spectrum data (`_spec.csv`) and system 
        lists (`_systems.csv`) based on the requested targets. Spectrum 
        exports automatically include the wavelength array (`x`), while 
        system exports automatically include the transition name (`series`) 
        and redshift (`z`).

        Parameters
        ----------
        targets : str
            Comma-separated list of structures or columns to export 
            (e.g., ``'spec, systems'`` or ``'F, cont, logN'``).
        path : str
            The destination directory path.
        prefix : str
            The base filename prefix for the exported files.

        Returns
        -------
        bool
            ``True`` if at least one file was successfully saved, ``False`` otherwise.
        """
        mapping = {'λ': 'x', 'λmin': 'xmin', 'λmax': 'xmax', 'F': 'y', 'dF': 'dy'}
        processed_targets = targets
        for ui_name, internal_name in mapping.items():
            # Replace whole words only
            processed_targets = re.sub(r'\b' + re.escape(ui_name) + r'\b', internal_name, processed_targets)
            
        target_list = [t.strip() for t in processed_targets.split(',') if t.strip()]
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
        
    def import_ascii_spec(self, path: str) -> 'SessionV2':
        """Import columns from an ASCII/CSV file into the current spectrum."""
        try:
            table = ascii.read(path)
            current_len = len(self._session.spec.x)
            import_len = len(table)

            if import_len != current_len:
                raise ValueError(f"Length mismatch: File has {import_len} pixels, but current session has {current_len}.")

            new_spec = self._session.spec
            imported_cols = []

            for col_name in table.colnames:
                if col_name.lower() == 'x':
                    continue
                
                values = np.array(table[col_name])
                new_spec = new_spec.update_column(col_name, values)
                imported_cols.append(col_name)

            logging.info(f"Imported columns from ASCII: {', '.join(imported_cols)}")
            return self._session.with_new_spectrum(new_spec)

        except Exception as e:
            logging.error(f"Import spectrum from ASCII failed: {e}", exc_info=True)
            return 0

    def import_ascii_systs(self, path: str, syst_mode: str = "append") -> 'SessionV2':
        """Import systems from an ASCII/CSV file into the current session."""
        try:
            table = ascii.read(path)

            if 'series' not in table.colnames or 'z' not in table.colnames:
                raise ValueError("System import requires 'series' and 'z' columns in the CSV.")

            from astrocook.core.structures import ComponentDataV2, SystemListDataV2
            from astrocook.core.system_list import SystemListV2
            import uuid

            current_systs = self._session.systs
            components = []
            
            if syst_mode.lower() == 'append' and current_systs and current_systs.components:
                components = list(current_systs.components)
                next_id = max(c.id for c in components) + 1
            else:
                next_id = 1
                
            for row in table:
                def safe_float(col_name, default=None):
                    if col_name not in table.colnames: return default
                    val = row[col_name]
                    if np.ma.is_masked(val) or (isinstance(val, (float, np.floating)) and np.isnan(val)):
                        return default
                    try: return float(val)
                    except (ValueError, TypeError): return default

                comp = ComponentDataV2(
                    id=next_id,
                    z=safe_float('z', 0.0),
                    dz=safe_float('dz'),
                    logN=safe_float('logN', 13.5),
                    dlogN=safe_float('dlogN'),
                    b=safe_float('b', 10.0),
                    db=safe_float('db'),
                    btur=safe_float('btur', 0.0),
                    dbtur=safe_float('dbtur'),
                    func=str(row['func']) if 'func' in table.colnames else 'voigt',
                    series=str(row['series']),
                    chi2=safe_float('chi2'),
                    resol=safe_float('resol')
                )
                
                object.__setattr__(comp, 'uuid', str(uuid.uuid4()))
                components.append(comp)
                next_id += 1
            
            new_data = SystemListDataV2(components=components)
            new_systs = SystemListV2(data=new_data)

            # 1. Update the session with the new system list
            new_session = self._session.with_new_system_list(new_systs)
            
            # --- NEW: Automatically regenerate the model ---
            try:
                from astrocook.fitting.voigt_fitter import VoigtFitterV2
                
                # Initialize the fitter with the updated session data
                fitter = VoigtFitterV2(new_session.spec, new_session.systs)
                
                # Compute the theoretical flux (without fitting/changing parameters)
                _, model_flux = fitter.compute_model_flux()
                
                # Inject the new model array into the spectrum's columns
                new_spec = new_session.spec.update_column('model', model_flux)
                new_session = new_session.with_new_spectrum(new_spec)
                
                logging.info("Successfully generated 'model' column from imported systems.")
            except ImportError:
                logging.warning("VoigtFitterV2 not found. Skipping automatic model generation.")
            except Exception as e:
                logging.warning(f"Could not generate model flux during import: {e}")
            # -----------------------------------------------
            
            logging.info(f"Imported {len(table)} systems from ASCII ({syst_mode} mode).")
            return new_session

        except Exception as e:
            logging.error(f"Import systems from ASCII failed: {e}", exc_info=True)
            return 0