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

# --- V2 Imports ---
from .v1_stubs import load_v1_spec_object, load_v1_systs_object
if TYPE_CHECKING:
    from ..session import SessionV2
from ..structures import (
    SessionMetadataV2, DataColumnV2, SpectrumDataV2, 
    SystemListDataV2, ComponentDataV2,
    HistoryLogV2, V1LogArtifact, LogEntryV2
)
from ..spectrum import SpectrumV2 
from ..system_list import SystemListV2
# --- V1 Imports ---
from ...v1.gui_log import GUILog

def v1_table_to_data_v2(v1_spectrum_instance: Any) -> SpectrumDataV2:
    """
    Adapter: Converte un'istanza Spectrum V1 (basata su Frame/Table) in SpectrumDataV2 immutabile.
    """
    
    t_v1 = v1_spectrum_instance._t 
    meta_v1 = v1_spectrum_instance._meta
    x_unit = v1_spectrum_instance._xunit
    v1_y_unit = v1_spectrum_instance._yunit # <-- Get the *original* y_unit
    y_unit = au.dimensionless_unscaled

    x_col = DataColumnV2(t_v1['x'].value, x_unit, description="Channels")
    y_col = DataColumnV2(t_v1['y'].value, y_unit, description="Flux density")
    dy_col = DataColumnV2(t_v1['dy'].value, y_unit, description="Error on Flux")
    xmin_col = DataColumnV2(t_v1['xmin'].value, x_unit, description="Lower channel limit")
    xmax_col = DataColumnV2(t_v1['xmax'].value, x_unit, description="Upper channel limit")

    aux_cols = {}

    if 'resol' not in t_v1.colnames:
        nan_array = np.full_like(t_v1['x'].value, np.nan)
        aux_cols['resol'] = DataColumnV2(nan_array, au.dimensionless_unscaled, description="Spectral Resolution")

    for colname in t_v1.colnames:
        if colname not in ['x', 'xmin', 'xmax', 'y', 'dy']:
            try:
                col_data = t_v1[colname]
                col_unit = col_data.unit if col_data.unit is not None else au.dimensionless_unscaled
                if v1_y_unit is not None and col_unit.is_equivalent(v1_y_unit):
                    col_unit = au.dimensionless_unscaled
                    logging.debug(f"Forcing y-like column '{colname}' to dimensionless.")
                    
                aux_cols[colname] = DataColumnV2(col_data.value, col_unit, description=f"Auxiliary column: {colname}")
            except Exception as e:
                logging.warning(f"Could not map auxiliary column {colname} from V1 table: {e}")
                continue

    return SpectrumDataV2(
        x=x_col, xmin=xmin_col, xmax=xmax_col, 
        y=y_col, dy=dy_col, 
        aux_cols=aux_cols, 
        meta=meta_v1,
        z_rf=getattr(v1_spectrum_instance, '_rfz', 0.0), # Read old _rfz
        z_em=getattr(v1_spectrum_instance, '_zem', 0.0)  # Try to get _zem
    )

def load_and_migrate_structure(
    archive_root: str, 
    structure_name: str, 
    gui_context: Any, 
    format_name: str, 
    spec_file_path: Optional[str] = None,
    v2_metadata: Optional[Dict] = None
):
    """
    General purpose function to load a V1 or V2 structure from file 
    and return the V2 API-level object.
    """
    
    path_to_load = None
    if structure_name == 'spec':
        if not spec_file_path:
             raise ValueError("spec_file_path is required for loading 'spec'")
        path_to_load = spec_file_path
        
    elif structure_name == 'systs':
        path_to_load = f"{archive_root}_systs.fits"
    else:
        logging.warning(f"Unsupported structure name requested: {structure_name}")
        return None
        
    if structure_name != 'systs' and (not path_to_load or not os.path.exists(path_to_load)):
         raise FileNotFoundError(f"Required file for structure '{structure_name}' not found at {path_to_load}")
         
    if structure_name == 'spec':
        
        # --- *** START V2 NATIVE LOADING PATH *** ---
        if v2_metadata:
            try:
                logging.debug("Using V2 native spectrum loader.")
                return load_spec_data_v2_from_archive(path_to_load)
            except Exception as e:
                logging.error(f"V2 native spec loading failed: {e}. Falling back to V1 migration.")
        # --- *** END V2 NATIVE LOADING PATH *** ---

        # --- V1 MIGRATION PATH (or fallback) ---
        logging.debug(f"Attempting to load V1 spectrum object from: {path_to_load}")
        v1_spec = load_v1_spec_object(path_to_load, format_name, gui_context)

        if not v1_spec:
             raise RuntimeError(f"Failed to load V1 spectrum object from {path_to_load}. load_v1_spec_object returned falsy value.")

        try:
            logging.debug("Migrating V1 spectrum object to V2 data core.")
            data_core_v2 = v1_table_to_data_v2(v1_spec)
            logging.debug("Wrapping V2 spectrum data core in API object.")
            return SpectrumV2(data=data_core_v2)
        except Exception as e:
             raise RuntimeError(f"Failed to migrate V1 spectrum data to V2: {e}")

    elif structure_name == 'systs':
        if not os.path.exists(path_to_load):
            logging.info(f"Optional system list file not found: {path_to_load}. Returning empty SystemListV2.")
            return SystemListV2(data=SystemListDataV2()) 

        if v2_metadata:
            try:
                logging.debug("Using V2 native system list loader.")
                return load_systs_data_v2_from_archive(path_to_load, v2_metadata)
            except Exception as e:
                logging.error(f"V2 native systs loading failed: {e}. Falling back to V1 migration path.")

        from ..system_list_migration import migrate_system_list_v1_to_v2
        try:
            logging.debug(f"Attempting to load V1 system list object from: {path_to_load}")
            v1_systs, syst_header = load_v1_systs_object(path_to_load)
        except Exception as e:
            logging.error(f"Failed during V1 systs object loading from {path_to_load}: {e}")
            return SystemListV2(data=SystemListDataV2())
        if v1_systs:
            try:
                logging.debug("Migrating V1 system list object to V2 data core.")
                data_core_v2 = migrate_system_list_v1_to_v2(v1_systs, syst_header)
                logging.debug("Wrapping V2 system list data core in API object.")
                return SystemListV2(data=data_core_v2)
            except Exception as e:
                raise RuntimeError(f"Failed to migrate V1 system list data to V2: {e}")
        else:
             logging.warning(f"V1 system list loading from {path_to_load} returned no data.")
             return SystemListV2(data=SystemListDataV2())

    return None

# --- NEW V2 LOADING HELPERS ---

def _v2_table_to_spectrum_data(spec_table: Table) -> SpectrumDataV2:
    """
    Converts a V2-native Astropy Table (from _spec.fits) back into
    a SpectrumDataV2 object.
    """
    t = spec_table
    meta = t.meta
    
    x_col = DataColumnV2(t['x'].value, t['x'].unit)
    y_col = DataColumnV2(t['y'].value, t['y'].unit)
    dy_col = DataColumnV2(t['dy'].value, t['dy'].unit)
    xmin_col = DataColumnV2(t['xmin'].value, t['xmin'].unit)
    xmax_col = DataColumnV2(t['xmax'].value, t['xmax'].unit)
    
    aux_cols = {}
    core_cols = ['x', 'y', 'dy', 'xmin', 'xmax']
    for colname in t.colnames:
        if colname not in core_cols:
            col_data = t[colname]
            col_unit = col_data.unit if col_data.unit is not None else au.dimensionless_unscaled
            aux_cols[colname] = DataColumnV2(col_data.value, col_unit, description=f"Auxiliary column: {colname}")
            
    z_rf = float(meta.pop('Z_RF', 0.0))
    z_em = float(meta.pop('Z_EM', 0.0))
    
    # [FIX] Extract Resolution info and hydrate the scalar 'resol' field
    resol_val = 0.0
    
    # Check common keys for resolution
    if 'resol_fwhm' in meta:
        # If we have FWHM in km/s, calculate approx R
        fwhm = float(meta['resol_fwhm'])
        if fwhm > 0:
            c_kms = 299792.458
            resol_val = c_kms / fwhm
    elif 'RESOL' in meta: # Scalar R
        resol_val = float(meta['RESOL'])

    return SpectrumDataV2(
        x=x_col, xmin=xmin_col, xmax=xmax_col, 
        y=y_col, dy=dy_col, 
        aux_cols=aux_cols, 
        meta=meta,
        z_rf=z_rf,
        z_em=z_em,
        resol=resol_val # <--- Pass this to the constructor
    )

def load_spec_data_v2_from_archive(spec_fits_path: str) -> 'SpectrumV2':
    """
    Loads a V2-native spectrum from its _spec.fits FITS table.
    """
    try:
        spec_table = Table.read(spec_fits_path)
    except FileNotFoundError:
        logging.error(f"V2 load: _spec.fits file not found at {spec_fits_path}.")
        raise
    data_core = _v2_table_to_spectrum_data(spec_table)
    return SpectrumV2(data=data_core)


def _parse_v2_metadata_constraints(v2_metadata: Dict) -> Dict[str, Dict]:
    """
    Loads the V2-native constraint map from the metadata dictionary.
    This map is: {UUID-str: {param_name-str: {constraint_dict}}}
    """
    v2_constraints_map = v2_metadata.get('constraints_by_uuid', {})
    
    # This log line was redundant and confusing.
    # logging.info(f"Loaded {len(v2_constraints_map)} V2 constraint sets from metadata.")
    
    return v2_constraints_map

def _v2_table_to_component_list(systs_table: Table) -> (List[ComponentDataV2], Dict[int, str]):
    """
    Converts a V2-native Astropy Table (from _systs.fits) back into
    a list of ComponentDataV2 objects and the id-to-uuid map.
    """
    components = []
    id_to_uuid_map = {}
    
    for row in systs_table:
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
        object.__setattr__(comp, 'uuid', str(row['uuid']))
        components.append(comp)
        id_to_uuid_map[comp.id] = comp.uuid
        
    return components, id_to_uuid_map

def load_systs_data_v2_from_archive(systs_fits_path: str, v2_metadata: Dict) -> 'SystemListV2':
    """
    Loads a V2-native system list by combining the _systs.fits FITS table
    and the _meta.json metadata dictionary.
    """
    from ..system_list import SystemListV2

    try:
        systs_table = Table.read(systs_fits_path)
    except FileNotFoundError:
        logging.warning(f"V2 load: _systs.fits file not found at {systs_fits_path}. Returning empty list.")
        return SystemListV2(data=SystemListDataV2())
    
    components, id_to_uuid_map = _v2_table_to_component_list(systs_table)
    
    # --- *** FIX: We are loading V2 constraints, not parsing V1 *** ---
    # This now returns the V2-native map: Dict[UUID, Dict[Param, Dict]]
    v2_constraints_map = _parse_v2_metadata_constraints(v2_metadata)
    # --- *** END FIX *** ---
    
    meta = systs_table.meta
    
    data_core = SystemListDataV2(
        components=components,
        # We need a new field in SystemListDataV2 for V2 constraints
        # or we need to convert them.
        # Let's assume SystemListV2 constructor handles this.
        parsed_constraints={}, # V1 constraints are empty
        v2_constraints_map=v2_constraints_map, # Pass the V2 map
        v1_id_to_uuid_map=id_to_uuid_map,
        meta=meta
    )
    return SystemListV2(data=data_core)


# --- METADATA SERIALIZATION (Updated) ---

def _serialize_v2_metadata(session: 'SessionV2') -> str:
    """Utility to extract and serialize V2 metadata (constraints, log) to JSON."""
    
    # 1. Extract V2 constraint definitions (UUID-keyed)
    # This was fixed to `v2_constraints_by_uuid`
    constraint_map_obj = session.systs.constraint_model.v2_constraints_by_uuid

    log_manager_to_save = getattr(session, 'log_manager', None)
    
    # Check if this is a V1-loaded session (has GUILog)
    v1_log_stub = getattr(session, 'log', None)
    if isinstance(v1_log_stub, GUILog) and v1_log_stub.json['set_menu']:
        # This V1 stub has been updated. Use it as the source of truth.
        log_manager_to_save = v1_log_stub
        logging.debug("Using session.log (GUILog) as save source.")
    else:
        logging.debug(f"Using session.log_manager ({type(log_manager_to_save)}) as save source.")

    log_history_data = None
    log_type = "unknown"

    if isinstance(log_manager_to_save, HistoryLogV2):
        log_history_data = dataclasses.asdict(log_manager_to_save)
        log_type = "v2"
        logging.debug(f"Serializing HistoryLogV2 (type='v2') with {len(log_manager_to_save.entries)} entries.")
        
    elif isinstance(log_manager_to_save, V1LogArtifact):
        log_history_data = log_manager_to_save.v1_json
        log_type = "v1_artifact"
        logging.debug("Serializing V1LogArtifact (type='v1_artifact')")

    elif isinstance(log_manager_to_save, GUILog):
        try:
            log_history_data = json.loads(log_manager_to_save.str)
        except json.JSONDecodeError:
            log_history_data = {"set_menu": []} 
        log_type = "v1_legacy"
        logging.debug("Serializing V1 GUILog (type='v1_legacy')")
        
    else:
        logging.error(f"Unknown log type in session: {type(log_manager_to_save)}. Saving empty log.")
        log_history_data = {}
    
    # 3. Create the SessionMetadataV2 object
    metadata = SessionMetadataV2(
        constraints_by_uuid=constraint_map_obj,
        log_history_json=log_history_data,
        log_type=log_type,
        v1_reconstruction_data=None, 
        component_uuids=[c.uuid for c in session.systs.components]
    )

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
    
    for name, data_col in spec_data.aux_cols.items():
        t[name] = Column(data_col.values, unit=data_col.unit)
        
    t.meta.update(spec_data.meta)
    t.meta['ORIGIN'] = 'Astrocook V2'
    
    t.meta['Z_RF'] = spec_data.z_rf
    t.meta['Z_EM'] = spec_data.z_em

    return t

def _convert_syst_list_to_table(systs_data: SystemListDataV2) -> Table:
    """Converts the V2 SystemListDataV2 core into a FITS-compatible Astropy Table."""
    
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

    t = Table(data_dict)
    
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

def save_archive_v2(session_v2_final: 'SessionV2', 
                    session_v2_initial: Optional['SessionV2'], 
                    file_path: str):
    """
    Saves the SessionV2 to a V2-native .acs2 archive.
    
    Args:
        session_v2_final (SessionV2): The final, current state to save.
        session_v2_initial (Optional[SessionV2]): The initial, original state to save.
        file_path (str): The destination path.
    """
    temp_dir = None
    
    try:
        temp_dir = tempfile.mkdtemp()
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        
        # --- A. Final (Current) Data ---
        spec_fname = f"{base_name}_spec.fits"
        spec_path_temp = os.path.join(temp_dir, spec_fname)
        spec_table = _convert_spec_data_to_table(session_v2_final.spec._data)
        spec_table.write(spec_path_temp, format='fits', overwrite=True)
        logging.info(f"Wrote V2 final spectrum FITS file: {spec_fname}")
        
        systs_fname = None
        if session_v2_final.systs and session_v2_final.systs.components:
            systs_fname = f"{base_name}_systs.fits"
            systs_path_temp = os.path.join(temp_dir, systs_fname)
            systs_table = _convert_syst_list_to_table(session_v2_final.systs._data)
            systs_table.write(systs_path_temp, format='fits', overwrite=True)
            logging.info(f"Wrote V2 final system list FITS file: {systs_fname}")

        # --- B. Metadata (Log, Constraints, etc. from FINAL session) ---
        meta_fname = f"{base_name}_meta.json"
        meta_path_temp = os.path.join(temp_dir, meta_fname)
        # The log manager is attached to the final session state
        meta_json_str = _serialize_v2_metadata(session_v2_final) 
        with open(meta_path_temp, 'w') as f:
            f.write(meta_json_str)
        logging.info(f"Wrote V2 metadata JSON: {meta_fname}")
        
        # --- C. Original (Initial) Data (Point 2) ---
        spec_orig_fname = None
        systs_orig_fname = None
        if session_v2_initial:
            try:
                # Save original spectrum
                spec_orig_fname = f"{base_name}_spec_orig.fits"
                spec_orig_path_temp = os.path.join(temp_dir, spec_orig_fname)
                spec_orig_table = _convert_spec_data_to_table(session_v2_initial.spec._data)
                spec_orig_table.write(spec_orig_path_temp, format='fits', overwrite=True)
                logging.info(f"Wrote V2 original spectrum FITS file: {spec_orig_fname}")

                # Save original systems
                if session_v2_initial.systs and session_v2_initial.systs.components:
                    systs_orig_fname = f"{base_name}_systs_orig.fits"
                    systs_orig_path_temp = os.path.join(temp_dir, systs_orig_fname)
                    systs_orig_table = _convert_syst_list_to_table(session_v2_initial.systs._data)
                    systs_orig_table.write(systs_orig_path_temp, format='fits', overwrite=True)
                    logging.info(f"Wrote V2 original system list FITS file: {systs_orig_fname}")
            except Exception as e:
                logging.error(f"Failed to save original data, archive may be incomplete: {e}")
                spec_orig_fname = None
                systs_orig_fname = None
            
        # 3. Create the .acs2 (tar.gz) Archive
        final_path = file_path if file_path.lower().endswith(('.acs2', '.tar.gz')) else file_path + '.acs2'
        
        with tarfile.open(final_path, 'w:gz') as tar:
            # Add final data
            tar.add(spec_path_temp, arcname=spec_fname)
            if systs_fname:
                 tar.add(systs_path_temp, arcname=systs_fname)
            tar.add(meta_path_temp, arcname=meta_fname)
            
            # Add original data
            if spec_orig_fname:
                tar.add(spec_orig_path_temp, arcname=spec_orig_fname)
            if systs_orig_fname:
                tar.add(systs_orig_path_temp, arcname=systs_orig_fname)
            
        logging.info(f"V2 archive created successfully at: {final_path}")

    except Exception as e:
        logging.error(f"FATAL: V2 archive creation failed during I/O: {e}")
        raise
        
    finally:
        if temp_dir and os.path.exists(temp_dir):
            import shutil
            shutil.rmtree(temp_dir)