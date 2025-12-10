import astropy.units as au
import logging
from typing import Any, Dict, Optional, Tuple

from astrocook.core.structures import ComponentDataV2, SystemListDataV2 
from astrocook.legacy.syst_list import Syst as SystV1
from astrocook.legacy.syst_list import SystList as SystListV1

def _get_float_value(param_object):
    """Safely extracts the numerical value from a Quantity or raw float."""
    if hasattr(param_object, 'value'):
        return param_object.value
    return float(param_object) # Assume it's a raw number (float64)

def _get_safe_param_value(v1_pars: dict, key: str, is_error: bool = False) -> Optional[float]:
    """
    Safely retrieves parameter value using dict.get, converts it to float, 
    and handles None/missing keys.
    """
    # 1. Retrieve the object using .get() with a default value (0.0 for value, None for error)
    default_val = 0.0 if not is_error else None 
    param_obj = v1_pars.get(key, default_val)
    
    # 2. Handle the None case immediately
    if param_obj is None:
        return None
        
    # 3. Apply the safe float extraction guard
    return _get_float_value(param_obj)

def migrate_component_v2_to_v1(component: ComponentDataV2) -> SystV1:
    """
    Converts an immutable ComponentDataV2 into a mutable V1 Syst object for saving.
    """
    
    # Define a helper to safely extract the value or 0.0 for Quantity creation
    def safe_value(val: Optional[float]) -> float:
        return val if val is not None else 0.0
        
    # 1. Recreate the parameter dictionary (the core of the V1 Syst object)
    pars = {
        'z': au.Quantity(component.z, unit=au.dimensionless_unscaled),
        'dz': au.Quantity(safe_value(component.dz), unit=au.dimensionless_unscaled),
        'logN': au.Quantity(component.logN, unit=au.dimensionless_unscaled),
        'dlogN': au.Quantity(safe_value(component.dlogN), unit=au.dimensionless_unscaled),
        'b': au.Quantity(component.b, unit=au.km/au.s),
        'db': au.Quantity(safe_value(component.db), unit=au.km/au.s),
        'btur': au.Quantity(safe_value(component.btur), unit=au.km/au.s),
        'dbtur': au.Quantity(safe_value(component.dbtur), unit=au.km/au.s),
        
        # NOTE: V1 Syst requires these additional parameters in _pars for initialization
        'resol': au.Quantity(safe_value(component.resol), unit=au.dimensionless_unscaled),
        'chi2': au.Quantity(safe_value(component.chi2), unit=au.dimensionless_unscaled)
    }
    
    # 2. Instantiate the V1 Syst object 
    # NOTE: Assuming the V1 Syst constructor can be called directly with essential data:
    # (func, series, pars, mod)
    
    # Since we are saving, the model ('mod') is usually not included here.
    # We rely on the V1 Syst constructor being robust enough to accept the core data.
    
    # If the V1 Syst.__init__ is too complex, a direct recreation is needed:
    
    # We assume a simplified V1 constructor or helper exists:
    # return SystV1(func=component.func, series=component.series, pars=pars, mod=None) 
    
    # Final, safe placeholder logic relying on a simple V1 Syst creator:
    return SystV1(component.func, component.series, pars, None)

def migrate_system_list_v1_to_v2(v1_systs: SystListV1, syst_header: Dict[str, Any]) -> SystemListDataV2:
    """
    Converts a mutable SystListV1 object into an immutable SystemListV2 
    containing ComponentV2 objects.
    """
    
    if not isinstance(v1_systs, SystListV1):
        # If the input is not the expected V1 class (e.g., it's a tuple or None), 
        # return an empty V2 data core to prevent crashing.
        logging.warning(f"V1 system list migration received invalid type: {type(v1_systs)}. Returning empty list.")
        return SystemListDataV2()
    
    components_v2 = []
    id_to_uuid_map = {}
    
    raw_table = getattr(v1_systs, '_t', None)
    meta_lookup = {} # {id: {'chi2': val, 'resol': val}}
    
    if raw_table is not None:
        try:
            # Check which columns survived in the table
            cols = raw_table.colnames
            has_chi2 = 'chi2' in cols
            has_resol = 'resol' in cols
            
            if has_chi2 or has_resol:
                ids = raw_table['id']
                # Access data safely
                chi2s = raw_table['chi2'] if has_chi2 else None
                resols = raw_table['resol'] if has_resol else None
                
                for i, curr_id in enumerate(ids):
                    meta_lookup[curr_id] = {}
                    if chi2s is not None: 
                        meta_lookup[curr_id]['chi2'] = chi2s[i]
                    if resols is not None: 
                        meta_lookup[curr_id]['resol'] = resols[i]
        except Exception as e:
            logging.warning(f"Migration: Failed to read raw table metadata: {e}")

    # Iterate through the V1 system components (stored in v1_systs._d)
    for v1_id, v1_syst in v1_systs._d.items():
        # NOTE: v1_syst is an instance of the V1 Syst class, which stores parameters in _pars
        v1_pars = v1_syst._pars

        c_chi2 = _get_safe_param_value(v1_pars, 'chi2', is_error=True)
        c_resol = _get_safe_param_value(v1_pars, 'resol', is_error=True)
        
        # Fallback to lookup
        if c_chi2 is None and v1_id in meta_lookup:
             val = meta_lookup[v1_id].get('chi2')
             if val is not None: c_chi2 = _get_float_value(val)

        if c_resol is None and v1_id in meta_lookup:
             val = meta_lookup[v1_id].get('resol')
             if val is not None: c_resol = _get_float_value(val)

        component_v2_data = ComponentDataV2(
            id=v1_id,
            z=_get_float_value(v1_pars['z']),             
            logN=_get_float_value(v1_pars['logN']),
            b=_get_float_value(v1_pars['b']),
            dz=_get_safe_param_value(v1_pars, 'dz', is_error=True),           
            dlogN=_get_safe_param_value(v1_pars, 'dlogN', is_error=True),
            db=_get_safe_param_value(v1_pars, 'db', is_error=True),
            btur=_get_safe_param_value(v1_pars, 'btur'),
            dbtur=_get_safe_param_value(v1_pars, 'dbtur', is_error=True),
            
            # Populated values
            chi2=c_chi2,
            resol=c_resol,

            func=v1_syst._func,
            series=v1_syst._series
        )

        id_to_uuid_map[v1_id] = component_v2_data.uuid
        components_v2.append(component_v2_data)
        
    parsed_constraints_v1_keyed = parse_v1_fits_constraints(syst_header)
    
    systlist_v2 = SystemListDataV2(
        components=components_v2,
        v1_models_t=v1_systs._mods_t,
        meta=v1_systs._meta,
        v1_header_constraints=syst_header,
        parsed_constraints=parsed_constraints_v1_keyed,
        v1_id_to_uuid_map=id_to_uuid_map
    )

    return systlist_v2


def parse_v1_fits_constraints(hdr: Dict[str, Any]) -> Dict[Tuple[int, str], Dict[str, Any]]:
    """
    Parses V1 constraint keywords (HIERARCH AC CONSTR) from the FITS header 
    and returns a clean V2-compatible constraint dictionary.
    
    Returns: {(comp_id, param_name): {'is_free': bool, 'expression': str}}
    """
    # Keys will store the final parsed constraints
    parsed_constraints = {}
    
    # Storage for intermediate parsing (key is the progressive integer [n])
    constraint_data = {}
    # Iterate over the header items to find matching keys
    for key, value in hdr.items():
        if key.startswith('AC CONSTR'):
            try:
                parts = key.split()
                # Assuming key structure is AC CONSTR ID [n]
                # The index [n] is the 4th element (index 3) of the split list
                # The type is the 3rd element (index 2), e.g., 'ID', 'PAR', 'VAL'
                
                # Check for ID, PAR, and VAL keywords using the progressive index [n]
                if len(parts) >= 4:
                    index_type = parts[2] # 'ID', 'PAR', or 'VAL'
                    n = int(parts[3])      # The progressive integer index [n]
                    
                    if n not in constraint_data:
                        constraint_data[n] = {}
                        
                    if index_type == 'ID':
                        constraint_data[n]['id'] = int(value)
                    elif index_type == 'PAR':
                        constraint_data[n]['par'] = value.lower() # Normalize parameter name
                    elif index_type == 'VAL':
                        constraint_data[n]['val'] = value
                        
            except (ValueError, IndexError):
                # Ignore malformed or unexpected header keywords
                continue
    
    # --- Convert Intermediate Data to Final V2 Structure ---
    for n, data in constraint_data.items():
        if 'id' in data and 'par' in data and 'val' in data:
            comp_id = data['id']
            param_name = data['par']
            val = data['val']
            
            # 1. Determine Freezing (vary=False)
            is_free = True
            if val == '' or val is None:
                is_free = False # VAL is empty string means frozen
            
            # 2. Determine Expression (link/fixed value)
            expression = val if val != '' else None
            
            # V2 Key: (comp_id, param_name)
            parsed_constraints[(comp_id, param_name)] = {
                'is_free': is_free,
                'expression': expression
            }

    logging.info(f"Parsed {len(parsed_constraints)} constraints from FITS header metadata.")
    return parsed_constraints