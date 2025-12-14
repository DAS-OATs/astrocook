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
    from astrocook.core.session import SessionV2
from astrocook.core.structures import (
    SessionMetadataV2, DataColumnV2, SpectrumDataV2, 
    SystemListDataV2, ComponentDataV2,
    HistoryLogV2, V1LogArtifact, LogEntryV2
)
from astrocook.core.spectrum import SpectrumV2 
from astrocook.core.system_list import SystemListV2
# --- V1 Imports ---
from astrocook.legacy.gui_log import GUILog

# ... (v1_table_to_data_v2 remains unchanged) ...
def v1_table_to_data_v2(v1_spectrum_instance: Any) -> SpectrumDataV2:
    t_v1 = v1_spectrum_instance._t 
    meta_v1 = v1_spectrum_instance._meta
    x_unit = v1_spectrum_instance._xunit
    v1_y_unit = v1_spectrum_instance._yunit 
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
                aux_cols[colname] = DataColumnV2(col_data.value, col_unit, description=f"Auxiliary column: {colname}")
            except Exception: continue

    return SpectrumDataV2(
        x=x_col, xmin=xmin_col, xmax=xmax_col, 
        y=y_col, dy=dy_col, 
        aux_cols=aux_cols, 
        meta=meta_v1,
        z_em=getattr(v1_spectrum_instance, '_zem', 0.0)
    )

def load_and_migrate_structure(archive_root, structure_name, gui_context, format_name, spec_file_path=None, v2_metadata=None):
    """
    Loads a V1 or V2 structure. 
    CRITICAL CHANGE: If v2_metadata is present, we DO NOT catch exceptions.
    We let them propagate so the user knows why V2 loading failed.
    """
    path_to_load = spec_file_path if structure_name == 'spec' else f"{archive_root}_systs.fits"
    
    if structure_name != 'systs' and (not path_to_load or not os.path.exists(path_to_load)):
         raise FileNotFoundError(f"File not found: {path_to_load}")
         
    if structure_name == 'spec':
        if v2_metadata:
            # [FIX] Removed try/except. Fail loudly if V2 fails.
            return load_spec_data_v2_from_archive(path_to_load)

        v1_spec = load_v1_spec_object(path_to_load, format_name, gui_context)
        if not v1_spec: raise RuntimeError("V1 load failed.")
        return SpectrumV2(data=v1_table_to_data_v2(v1_spec))

    elif structure_name == 'systs':
        if not os.path.exists(path_to_load): return SystemListV2(data=SystemListDataV2()) 

        if v2_metadata:
            # [FIX] Removed try/except. Fail loudly if V2 fails.
            return load_systs_data_v2_from_archive(path_to_load, v2_metadata)

        from astrocook.core.system_list_migration import migrate_system_list_v1_to_v2
        try:
            v1_systs, syst_header = load_v1_systs_object(path_to_load)
            return SystemListV2(data=migrate_system_list_v1_to_v2(v1_systs, syst_header))
        except Exception: return SystemListV2(data=SystemListDataV2())

    return None

# ... (_v2_table_to_spectrum_data, load_spec_data_v2_from_archive, _parse_v2_metadata_constraints remain same) ...
def _v2_table_to_spectrum_data(spec_table: Table) -> SpectrumDataV2:
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
            
    z_em = float(meta.pop('Z_EM', 0.0))
    
    resol_val = 0.0
    if 'resol_fwhm' in meta:
        fwhm = float(meta['resol_fwhm'])
        if fwhm > 0: resol_val = 299792.458 / fwhm
    elif 'RESOL' in meta:
        resol_val = float(meta['RESOL'])

    return SpectrumDataV2(
        x=x_col, xmin=xmin_col, xmax=xmax_col, 
        y=y_col, dy=dy_col, 
        aux_cols=aux_cols, meta=meta, z_em=z_em, resol=resol_val
    )

def load_spec_data_v2_from_archive(spec_fits_path: str) -> 'SpectrumV2':
    spec_table = Table.read(spec_fits_path)
    data_core = _v2_table_to_spectrum_data(spec_table)
    return SpectrumV2(data=data_core)

def _parse_v2_metadata_constraints(v2_metadata: Dict) -> Dict[str, Dict]:
    return v2_metadata.get('constraints_by_uuid', {})

def _v2_table_to_component_list(systs_table: Table) -> (List[ComponentDataV2], Dict[int, str]):
    """
    Converts a V2-native Astropy Table back into ComponentDataV2 objects.
    Robustly handles masked values and optional columns.
    """
    components = []
    id_to_uuid_map = {}
    
    has_chi2 = 'chi2' in systs_table.colnames
    has_resol = 'resol' in systs_table.colnames

    for row in systs_table:
        # [FIX] Robust value extractor for Astropy Rows
        def safe_float(col_name):
            val = row[col_name]
            # Handle MaskedConstant (Astropy) or NaN (Numpy)
            if np.ma.is_masked(val) or (isinstance(val, float) and np.isnan(val)):
                return None
            return float(val)

        # Mandatory fields (should not be None, but good to be safe)
        z_val = safe_float('z') or 0.0
        logN_val = safe_float('logN') or 0.0
        b_val = safe_float('b') or 0.0
        btur_val = safe_float('btur') or 0.0

        comp = ComponentDataV2(
            id=int(row['id']),
            z=z_val,
            dz=safe_float('dz'),
            logN=logN_val,
            dlogN=safe_float('dlogN'),
            b=b_val,
            db=safe_float('db'),
            btur=btur_val,
            dbtur=safe_float('dbtur'),
            
            # [FIX] Read metadata using safe extractor
            chi2=safe_float('chi2') if has_chi2 else None,
            resol=safe_float('resol') if has_resol else None,

            func=str(row['func']),
            series=str(row['series'])
        )
        object.__setattr__(comp, 'uuid', str(row['uuid']))
        components.append(comp)
        id_to_uuid_map[comp.id] = comp.uuid
        
    return components, id_to_uuid_map

# ... (load_systs_data_v2_from_archive, _serialize_v2_metadata remain same) ...
def load_systs_data_v2_from_archive(systs_fits_path: str, v2_metadata: Dict) -> 'SystemListV2':
    from astrocook.core.system_list import SystemListV2
    try:
        systs_table = Table.read(systs_fits_path)
    except FileNotFoundError:
        return SystemListV2(data=SystemListDataV2())
    
    components, id_to_uuid_map = _v2_table_to_component_list(systs_table)
    v2_constraints_map = _parse_v2_metadata_constraints(v2_metadata)
    meta = systs_table.meta
    
    data_core = SystemListDataV2(
        components=components,
        parsed_constraints={},
        v2_constraints_map=v2_constraints_map,
        v1_id_to_uuid_map=id_to_uuid_map,
        meta=meta
    )
    return SystemListV2(data=data_core)

def _serialize_v2_metadata(session: 'SessionV2') -> str:
    constraint_map_obj = session.systs.constraint_model.v2_constraints_by_uuid
    log_manager_to_save = getattr(session, 'log_manager', None)
    
    log_history_data = None
    log_type = "unknown"

    if isinstance(log_manager_to_save, HistoryLogV2):
        log_history_data = dataclasses.asdict(log_manager_to_save)
        log_type = "v2"
    elif isinstance(log_manager_to_save, V1LogArtifact):
        log_history_data = log_manager_to_save.v1_json
        log_type = "v1_artifact"
    elif isinstance(log_manager_to_save, GUILog):
        try: log_history_data = json.loads(log_manager_to_save.str)
        except json.JSONDecodeError: log_history_data = {"set_menu": []} 
        log_type = "v1_legacy"
    else:
        log_history_data = {}
    
    metadata = SessionMetadataV2(
        constraints_by_uuid=constraint_map_obj,
        log_history_json=log_history_data,
        log_type=log_type,
        v1_reconstruction_data=None, 
        component_uuids=[c.uuid for c in session.systs.components]
    )
    return json.dumps(dataclasses.asdict(metadata), indent=2)

def _convert_spec_data_to_table(spec_data: SpectrumDataV2) -> Table:
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
    t.meta['Z_EM'] = spec_data.z_em
    return t

def _convert_syst_list_to_table(systs_data: SystemListDataV2) -> Table:
    """
    Converts the V2 SystemListDataV2 core into a FITS-compatible Astropy Table.
    Now supports Chi2 and Resol columns.
    """
    data_dict = {
        'uuid': [], 'id': [], 'z': [], 'dz': [], 'logN': [], 'dlogN': [], 
        'b': [], 'db': [], 'btur': [], 'dbtur': [], 'func': [], 'series': [],
        'chi2': [], 'resol': []
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
        
        # [FIX] Write Metadata
        data_dict['chi2'].append(c.chi2 if c.chi2 is not None else np.nan)
        data_dict['resol'].append(c.resol if c.resol is not None else np.nan)

    t = Table(data_dict)
    
    t['z'].unit = au.dimensionless_unscaled
    t['dz'].unit = au.dimensionless_unscaled
    t['logN'].unit = au.dimensionless_unscaled
    t['dlogN'].unit = au.dimensionless_unscaled
    t['b'].unit = au.km/au.s
    t['db'].unit = au.km/au.s
    t['btur'].unit = au.km/au.s
    t['dbtur'].unit = au.km/au.s
    # [FIX] Units for metadata
    t['chi2'].unit = au.dimensionless_unscaled
    t['resol'].unit = au.dimensionless_unscaled

    t.meta.update(systs_data.meta)
    return t

# ... (save_archive_v2 remains same) ...
def save_archive_v2(session_v2_final: 'SessionV2', session_v2_initial: Optional['SessionV2'], file_path: str):
    temp_dir = None
    try:
        temp_dir = tempfile.mkdtemp()
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        
        spec_fname = f"{base_name}_spec.fits"
        spec_path_temp = os.path.join(temp_dir, spec_fname)
        spec_table = _convert_spec_data_to_table(session_v2_final.spec._data)
        spec_table.write(spec_path_temp, format='fits', overwrite=True)
        
        systs_fname = None
        if session_v2_final.systs and session_v2_final.systs.components:
            systs_fname = f"{base_name}_systs.fits"
            systs_path_temp = os.path.join(temp_dir, systs_fname)
            systs_table = _convert_syst_list_to_table(session_v2_final.systs._data)
            systs_table.write(systs_path_temp, format='fits', overwrite=True)

        meta_fname = f"{base_name}_meta.json"
        meta_path_temp = os.path.join(temp_dir, meta_fname)
        meta_json_str = _serialize_v2_metadata(session_v2_final) 
        with open(meta_path_temp, 'w') as f: f.write(meta_json_str)
        
        spec_orig_fname = None
        systs_orig_fname = None
        if session_v2_initial:
            try:
                spec_orig_fname = f"{base_name}_spec_orig.fits"
                spec_orig_path_temp = os.path.join(temp_dir, spec_orig_fname)
                spec_orig_table = _convert_spec_data_to_table(session_v2_initial.spec._data)
                spec_orig_table.write(spec_orig_path_temp, format='fits', overwrite=True)

                if session_v2_initial.systs and session_v2_initial.systs.components:
                    systs_orig_fname = f"{base_name}_systs_orig.fits"
                    systs_orig_path_temp = os.path.join(temp_dir, systs_orig_fname)
                    systs_orig_table = _convert_syst_list_to_table(session_v2_initial.systs._data)
                    systs_orig_table.write(systs_orig_path_temp, format='fits', overwrite=True)
            except Exception as e:
                logging.error(f"Failed to save original data: {e}")
            
        final_path = file_path if file_path.lower().endswith(('.acs2', '.tar.gz')) else file_path + '.acs2'
        
        with tarfile.open(final_path, 'w:gz') as tar:
            tar.add(spec_path_temp, arcname=spec_fname)
            if systs_fname: tar.add(systs_path_temp, arcname=systs_fname)
            tar.add(meta_path_temp, arcname=meta_fname)
            if spec_orig_fname: tar.add(spec_orig_path_temp, arcname=spec_orig_fname)
            if systs_orig_fname: tar.add(systs_orig_path_temp, arcname=systs_orig_fname)
            
        logging.info(f"V2 archive created successfully at: {final_path}")

    except Exception as e:
        logging.error(f"FATAL: V2 archive creation failed: {e}")
        raise
        
    finally:
        if temp_dir and os.path.exists(temp_dir):
            import shutil
            shutil.rmtree(temp_dir)