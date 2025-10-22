from typing import Optional

from .structures import ComponentDataV2, SystemListDataV2 
from astrocook.v1.syst_list import SystList as SystListV1

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

def migrate_system_list_v1_to_v2(v1_systs: SystListV1) -> SystemListDataV2:
    """
    Converts a mutable SystListV1 object into an immutable SystemListV2 
    containing ComponentV2 objects.
    """
    
    components_v2 = []
    
    # Iterate through the V1 system components (stored in v1_systs._d)
    for v1_id, v1_syst in v1_systs._d.items():
        # NOTE: v1_syst is an instance of the V1 Syst class, which stores parameters in _pars
        v1_pars = v1_syst._pars

        # 1. Create the immutable ComponentV2 object
        component_v2_data = ComponentDataV2(
            # Apply the guard for every parameter extraction
            z=_get_float_value(v1_pars['z']),             
            logN=_get_float_value(v1_pars['logN']),
            b=_get_float_value(v1_pars['b']),
            
            # Safe extraction for ALL error terms (dz, dlogN, db) ---
            dz=_get_safe_param_value(v1_pars, 'dz', is_error=True),           
            dlogN=_get_safe_param_value(v1_pars, 'dlogN', is_error=True),
            db=_get_safe_param_value(v1_pars, 'db', is_error=True),

            # Safe extraction of btur/dbtur ---
            btur=_get_safe_param_value(v1_pars, 'btur'),
            dbtur=_get_safe_param_value(v1_pars, 'dbtur', is_error=True),

            func=v1_syst._func,
            series=v1_syst._series,
            # The 'id' is set automatically by ComponentV2.__post_init__ or will be assigned later
        )
        components_v2.append(component_v2_data)
        
    # 2. Create the V2 container
    # CRITICAL: We pass the mutable V1 lmfit models (v1_systs._mods_t) as a placeholder 
    # until the V2 ConstraintModelV2 is implemented.
    systlist_v2 = SystemListDataV2(
        components=components_v2,
        v1_models_t=v1_systs._mods_t, # Placeholder for V1 lmfit models
        meta=v1_systs._meta,
    )
    
    return systlist_v2