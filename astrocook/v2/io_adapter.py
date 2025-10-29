from astropy import units as au
from astropy.table import Table, Column
import dataclasses
import json
import logging
import numpy as np
import os
import tarfile
import tempfile
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

from .io_v1_stubs import load_v1_spec_object, load_v1_systs_object
if TYPE_CHECKING:
    from .session import SessionV2
from .structures import ( # <<< IMPORT ADDED
    SessionMetadataV2, DataColumnV2, SpectrumDataV2, 
    SystemListDataV2, ComponentDataV2
)
from .spectrum import SpectrumV2 
from .system_list import SystemListV2

def v1_table_to_data_v2(v1_spectrum_instance: Any) -> SpectrumDataV2:
    """
    Adapter: Converte un'istanza Spectrum V1 (basata su Frame/Table) in SpectrumDataV2 immutabile.
    """
    
    # Assumiamo che l'istanza V1 abbia l'attributo _t (astropy.Table) e _xunit, _yunit, _meta
    t_v1 = v1_spectrum_instance._t 
    meta_v1 = v1_spectrum_instance._meta
    x_unit = v1_spectrum_instance._xunit
    y_unit = v1_spectrum_instance._yunit

    # 1. Colonne Core
    x_col = DataColumnV2(t_v1['x'].value, x_unit, description="Channels")
    y_col = DataColumnV2(t_v1['y'].value, y_unit, description="Flux density")
    dy_col = DataColumnV2(t_v1['dy'].value, y_unit, description="Error on Flux")
    xmin_col = DataColumnV2(t_v1['xmin'].value, x_unit, description="Lower channel limit")
    xmax_col = DataColumnV2(t_v1['xmax'].value, x_unit, description="Upper channel limit")

    # 2. Colonne Ausiliarie (Cont, Resol, Model, etc.)
    aux_cols = {}

    if 'resol' not in t_v1.colnames:
        # Create a default 'resol' column with NaNs if missing
        nan_array = np.full_like(t_v1['x'].value, np.nan)
        # Use a safe default unit, e.g., dimensionless or km/s (as often used in V1)
        aux_cols['resol'] = DataColumnV2(nan_array, au.dimensionless_unscaled, description="Spectral Resolution")

    for colname in t_v1.colnames:
        # Check if column is a core column before processing aux_cols
        if colname not in ['x', 'xmin', 'xmax', 'y', 'dy']:
            # Use defensive coding to check if the column is accessible
            try:
                col_data = t_v1[colname]
                col_unit = col_data.unit if col_data.unit is not None else au.dimensionless_unscaled
                aux_cols[colname] = DataColumnV2(col_data.value, col_unit, description=f"Auxiliary column: {colname}")
            except Exception as e:
                logging.warning(f"Could not map auxiliary column {colname} from V1 table: {e}")
                # Continue mapping other columns
                continue

    return SpectrumDataV2(
        x=x_col, xmin=xmin_col, xmax=xmax_col, 
        y=y_col, dy=dy_col, 
        aux_cols=aux_cols, 
        meta=meta_v1,
        rf_z=getattr(v1_spectrum_instance, '_rfz', 0.0)
    )

def load_and_migrate_structure(
    archive_root: str, 
    structure_name: str, 
    gui_context: Any, 
    format_name: str, 
    spec_file_path: Optional[str] = None,
    v2_metadata: Optional[Dict] = None  # <<< NEW ARGUMENT
):
    """
    General purpose function to load a V1 or V2 structure from file 
    and return the V2 API-level object.
    """
    
    path_to_load = None
    # 1. Determine the exact path based on structure type
    if structure_name == 'spec':
        if not spec_file_path:
             raise ValueError("spec_file_path is required for loading 'spec'")
        path_to_load = spec_file_path
        
    elif structure_name == 'systs':
        # Construct the path for the associated structure using the archive root
        path_to_load = f"{archive_root}_systs.fits"
        
    else:
        logging.warning(f"Unsupported structure name requested: {structure_name}")
        return None
        
    # Check if the determined path actually exists, except for systs (optional)
    if structure_name != 'systs' and (not path_to_load or not os.path.exists(path_to_load)):
         raise FileNotFoundError(f"Required file for structure '{structure_name}' not found at {path_to_load}")
         
    if structure_name == 'spec':
        # --- Spectrum Loading/Migration ---
        logging.debug(f"Attempting to load spectrum from: {path_to_load}")
        v1_spec = load_v1_spec_object(path_to_load, format_name, gui_context)

        # CRITICAL CHANGE: Raise error if spec loading fails (returns None or 0)
        if not v1_spec: # Handles None, 0, or potentially other falsy values
             raise RuntimeError(f"Failed to load V1 spectrum object from {path_to_load}. load_v1_spec_object returned falsy value.")

        try:
            logging.debug("Migrating V1 spectrum object to V2 data core.")
            data_core_v2 = v1_table_to_data_v2(v1_spec)
            logging.debug("Wrapping V2 spectrum data core in API object.")
            return SpectrumV2(data=data_core_v2)
        except Exception as e:
             # Catch potential errors during migration (e.g., missing columns)
             raise RuntimeError(f"Failed to migrate V1 spectrum data to V2: {e}")

    elif structure_name == 'systs':
        # Check if the systs file exists first
        if not os.path.exists(path_to_load):
            logging.info(f"Optional system list file not found: {path_to_load}. Returning empty SystemListV2.")
            # Return an empty SystemListV2 if the file is missing
            return SystemListV2(data=SystemListDataV2()) # Return empty object


        # --- V2 NATIVE LOADING PATH ---
        if v2_metadata:
            try:
                logging.debug("Using V2 native system list loader.")
                return load_systs_data_v2_from_archive(path_to_load, v2_metadata)
            except Exception as e:
                logging.error(f"V2 native systs loading failed: {e}. Falling back to V1 migration path.")
                # Fallback handled below

        # --- V1 MIGRATION PATH (or fallback) ---
        from .system_list_migration import migrate_system_list_v1_to_v2
        # from .system_list import SystemListV2 # Already imported above

        # Add try-except around V1 systs loading
        try:
            logging.debug(f"Attempting to load V1 system list object from: {path_to_load}")
            v1_systs, syst_header = load_v1_systs_object(path_to_load)
        except Exception as e:
             logging.error(f"Failed during V1 systs object loading from {path_to_load}: {e}")
             # Return empty list on V1 load failure too
             return SystemListV2(data=SystemListDataV2())


        if v1_systs:
            try:
                logging.debug("Migrating V1 system list object to V2 data core.")
                data_core_v2 = migrate_system_list_v1_to_v2(v1_systs, syst_header)
                logging.debug("Wrapping V2 system list data core in API object.")
                return SystemListV2(data=data_core_v2)
            except Exception as e:
                 # Catch potential errors during migration
                 raise RuntimeError(f"Failed to migrate V1 system list data to V2: {e}")
        else:
             # If load_v1_systs_object returned None/False without error
             logging.warning(f"V1 system list loading from {path_to_load} returned no data.")
             return SystemListV2(data=SystemListDataV2()) # Return empty

    # This path should ideally not be reached if spec/systs logic is exhaustive
    logging.error(f"load_and_migrate_structure reached end unexpectedly for {structure_name}.")
    return None

# --- NEW V2 LOADING HELPERS ---

def _parse_v2_metadata_constraints(v2_metadata: Dict) -> Dict[Tuple[int, str], Dict]:
    """
    Converts the V2 JSON metadata constraints (string keys) back into the
    V1-style tuple-key map (parsed_constraints) that VoigtModelConstraintV2 expects.
    
    Input:  {"1__z": {"is_free": false, ...}}
    Output: {(1, 'z'): {"is_free": false, ...}}
    """
    # This data is stored in the 'constraints_by_uuid' field (which we temporarily used)
    v1_constraints_str_key = v2_metadata.get('constraints_by_uuid', {})
    parsed_constraints = {}
    
    for str_key, data in v1_constraints_str_key.items():
        parts = str_key.split('__')
        if len(parts) == 2:
            try:
                comp_id = int(parts[0])
                param_name = parts[1]
                # Re-create the V1-style tuple key
                parsed_constraints[(comp_id, param_name)] = data
            except ValueError:
                logging.warning(f"Could not parse V2 metadata constraint key: {str_key}")
                continue
                
    logging.info(f"Loaded {len(parsed_constraints)} constraints from V2 metadata.")
    return parsed_constraints

def _v2_table_to_component_list(systs_table: Table) -> (List[ComponentDataV2], Dict[int, str]):
    """
    Converts a V2-native Astropy Table (from _systs.fits) back into
    a list of ComponentDataV2 objects and the id-to-uuid map.
    """
    components = []
    id_to_uuid_map = {}
    
    for row in systs_table:
        # Helper to convert NaNs from FITS table back to None for Optional[float]
        def nan_to_none(val):
            return None if (val is None or np.isnan(val)) else float(val)

        comp = ComponentDataV2(
            id=int(row['id']),
            z=float(row['z']),
            dz=nan_to_none(row['dz']),
            logN=float(row['logN']),
            dlogN=nan_to_none(row['dlogN']),
            b=float(row['b']),
            db=nan_to_none(row['db']),
            btur=float(row['btur']),
            dbtur=nan_to_none(row['dbtur']),
            func=str(row['func']),
            series=str(row['series'])
        )
        
        # CRITICAL: Override the auto-generated UUID with the saved one
        # We must use object.__setattr__ because the dataclass is frozen
        object.__setattr__(comp, 'uuid', str(row['uuid']))
        
        components.append(comp)
        id_to_uuid_map[comp.id] = comp.uuid
        
    return components, id_to_uuid_map

def load_systs_data_v2_from_archive(systs_fits_path: str, v2_metadata: Dict) -> 'SystemListV2':
    """
    Loads a V2-native system list by combining the _systs.fits FITS table
    and the _meta.json metadata dictionary.
    """
    from .system_list import SystemListV2

    # 1. Load FITS table
    try:
        systs_table = Table.read(systs_fits_path)
    except FileNotFoundError:
        logging.warning(f"V2 load: _systs.fits file not found at {systs_fits_path}. Returning empty list.")
        return SystemListV2(data=SystemListDataV2())
    
    # 2. Get components and ID map from the FITS table
    components, id_to_uuid_map = _v2_table_to_component_list(systs_table)
    
    # 3. Get V1-style constraints from the JSON metadata
    parsed_constraints = _parse_v2_metadata_constraints(v2_metadata)
    
    # 4. Get other metadata
    meta = systs_table.meta
    
    # 5. Build the data core
    data_core = SystemListDataV2(
        components=components,
        parsed_constraints=parsed_constraints,
        v1_id_to_uuid_map=id_to_uuid_map,
        meta=meta
        # v1_header_constraints and v1_models_t are None, which is correct
    )
    return SystemListV2(data=data_core)


# --- METADATA SERIALIZATION (Updated) ---

def _serialize_v2_metadata(session: 'SessionV2') -> str:
    """Utility to extract and serialize V2 metadata (constraints, log) to JSON."""
    
    # 1. Extract V2 constraint definitions (UUID-keyed)
    # This map is Dict[str, Dict[str, ParameterConstraintV2]]
    constraint_map_obj = session.systs.v2_constraints_for_save
    
    # 2. Extract History
    history_log = session.log.str # Using the V1 logger's JSON string output

    # 3. Create the SessionMetadataV2 object
    metadata = SessionMetadataV2(
        constraints_by_uuid=constraint_map_obj,
        log_history_json=history_log,
        v1_reconstruction_data=None, # We do not save the V1 lmfit models
        component_uuids=[c.uuid for c in session.systs.components]
    )

    # 4. Serialize the dataclass to JSON
    # dataclasses.asdict will recursively convert the main object and all
    # nested ParameterConstraintV2 objects into dicts.
    metadata_dict = dataclasses.asdict(metadata)
    
    return json.dumps(metadata_dict, indent=2)

def _convert_spec_data_to_table(spec_data: SpectrumDataV2) -> Table:
    """Converts the V2 SpectrumDataV2 core into a FITS-compatible Astropy Table."""
    
    t = Table()
    t['x'] = Column(spec_data.x.values, unit=spec_data.x.unit)
    t['xmin'] = Column(spec_data.xmin.values, unit=spec_data.xmin.unit)
    t['xmax'] = Column(spec_data.xmax.values, unit=spec_data.xmax.unit)
    t['y'] = Column(spec_data.y.values, unit=spec_data.y.unit)
    t['dy'] = Column(spec_data.dy.values, unit=spec_data.dy.unit)
    
    # Add auxiliary columns
    for name, data_col in spec_data.aux_cols.items():
        t[name] = Column(data_col.values, unit=data_col.unit)
        
    # Add V1-compatible metadata (ORIGIN)
    t.meta.update(spec_data.meta)
    t.meta['ORIGIN'] = 'Astrocook V2'
    
    return t

def _convert_syst_list_to_table(systs_data: SystemListDataV2) -> Table:
    """Converts the V2 SystemListDataV2 core into a FITS-compatible Astropy Table."""
    
    # 1. Extract data from immutable ComponentDataV2 objects
    data_dict = {
        'uuid': [], 'id': [], 'z': [], 'dz': [], 'logN': [], 'dlogN': [], 
        'b': [], 'db': [], 'btur': [], 'dbtur': [], 'func': [], 'series': []
    }
    
    for c in systs_data.components:
        data_dict['uuid'].append(c.uuid)
        data_dict['id'].append(c.id)
        data_dict['z'].append(c.z)
        data_dict['dz'].append(c.dz if c.dz is not None else np.nan)
        data_dict['logN'].append(c.logN)
        data_dict['dlogN'].append(c.dlogN if c.dlogN is not None else np.nan)
        data_dict['b'].append(c.b)
        data_dict['db'].append(c.db if c.db is not None else np.nan)
        data_dict['btur'].append(c.btur)
        data_dict['dbtur'].append(c.dbtur if c.dbtur is not None else np.nan)
        data_dict['func'].append(c.func)
        data_dict['series'].append(c.series)

    # 2. Build the Astropy Table
    t = Table(data_dict)
    
    # 3. Add Units (V2 standard)
    t['z'].unit = au.dimensionless_unscaled
    t['dz'].unit = au.dimensionless_unscaled
    t['logN'].unit = au.dimensionless_unscaled
    t['dlogN'].unit = au.dimensionless_unscaled
    t['b'].unit = au.km/au.s
    t['db'].unit = au.km/au.s
    t['btur'].unit = au.km/au.s
    t['dbtur'].unit = au.km/au.s

    t.meta.update(systs_data.meta)
    return t

# --- The Main V2 Archive Writer ---

def save_archive_v2(session_v2: 'SessionV2', file_path: str):
    """
    Saves the SessionV2 to a V2-native .acs2 archive by serializing data cores directly.
    """
    temp_dir = None
    
    try:
        # 1. Setup Temporary Directory
        temp_dir = tempfile.mkdtemp()
        
        # Use the filename stem as the base name for files inside the archive
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        
        # --- 2. Serialize V2 Structures ---
        
        # A. Spectrum (as FITS)
        spec_fname = f"{base_name}_spec.fits"
        spec_path_temp = os.path.join(temp_dir, spec_fname)
        spec_table = _convert_spec_data_to_table(session_v2.spec._data)
        spec_table.write(spec_path_temp, format='fits', overwrite=True)
        logging.info(f"Wrote V2 spectrum FITS file: {spec_fname}")
        
        # B. System List (as FITS)
        systs_fname = None
        if session_v2.systs and session_v2.systs.components:
            systs_fname = f"{base_name}_systs.fits"
            systs_path_temp = os.path.join(temp_dir, systs_fname)
            systs_table = _convert_syst_list_to_table(session_v2.systs._data)
            systs_table.write(systs_path_temp, format='fits', overwrite=True)
            logging.info(f"Wrote V2 system list FITS file: {systs_fname}")

        # C. Metadata (Log, Constraints, etc. as JSON)
        meta_fname = f"{base_name}_meta.json"
        meta_path_temp = os.path.join(temp_dir, meta_fname)
        meta_json_str = _serialize_v2_metadata(session_v2)
        with open(meta_path_temp, 'w') as f:
            f.write(meta_json_str)
        logging.info(f"Wrote V2 metadata JSON: {meta_fname}")
            
        # 3. Create the .acs2 (tar.gz) Archive
        final_path = file_path if file_path.lower().endswith(('.acs2', '.tar.gz')) else file_path + '.acs2'
        
        with tarfile.open(final_path, 'w:gz') as tar:
            tar.add(spec_path_temp, arcname=spec_fname)
            if systs_fname:
                 tar.add(systs_path_temp, arcname=systs_fname)
            tar.add(meta_path_temp, arcname=meta_fname)
            
        logging.info(f"V2 archive created successfully at: {final_path}")

    except Exception as e:
        logging.error(f"FATAL: V2 archive creation failed during I/O: {e}")
        raise
        
    finally:
        # 4. Clean up temporary directory
        if temp_dir and os.path.exists(temp_dir):
            import shutil
            shutil.rmtree(temp_dir)