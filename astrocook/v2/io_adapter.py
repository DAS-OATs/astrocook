from astropy import units as au
from astropy.table import Table
import logging
import numpy as np
from typing import Any, Optional

from .io_v1_stubs import load_v1_spec_object, load_v1_systs_object
from .structures import DataColumnV2, SpectrumDataV2
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
                col_unit = col_data.unit if col_data.unit else au.dimensionless_unscaled
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

def load_and_migrate_structure(archive_root: str, structure_name: str, gui_context: Any, format_name: str, spec_file_path: Optional[str] = None):
    """
    General purpose function to load a V1 structure from file and migrate it to V2 immutable object.
    
    :param archive_root: The clean, extension-less root path (e.g., /tmp/xyz/session_root).
    :param spec_file_path: The specific, resolved path for the 'spec' structure (used only for 'spec').
    """
    
    # 1. Determine the exact path based on structure type
    if structure_name == 'spec':
        # Use the explicit file path passed from the orchestrator (required for archive resolution)
        path_to_load = spec_file_path
        
    elif structure_name == 'systs':
        # Construct the path for the associated structure using the archive root
        path_to_load = f"{archive_root}_systs.fits"
        
    else:
        # Placeholder for other structure names
        return None
        
    
    if structure_name == 'spec':
        # --- Spectrum Migration ---
        # Call V1 loading utility with the full, resolved path
        v1_spec = load_v1_spec_object(path_to_load, format_name, gui_context) 
        if v1_spec:
            data_core_v2 = v1_table_to_data_v2(v1_spec)
            return SpectrumV2(data=data_core_v2)
        
    elif structure_name == 'systs':
        # --- System List Migration ---
        from .system_list_migration import migrate_system_list_v1_to_v2 
        from .system_list import SystemListV2
        
        # Call V1 loading utility (load_v1_systs_object needs only the path)
        v1_systs, syst_header = load_v1_systs_object(path_to_load) 
        
        if v1_systs:
            data_core_v2 = migrate_system_list_v1_to_v2(v1_systs, syst_header)
            return SystemListV2(data=data_core_v2)
            
    # For all other structures or if loading fails
    return None