from .structures import ComponentDataV2, SystemListDataV2 
from astrocook.v1.syst_list import SystList as SystListV1

def migrate_system_list_v1_to_v2(v1_systs: SystListV1) -> SystemListDataV2:
    """
    Converts a mutable SystListV1 object into an immutable SystemListV2 
    containing ComponentV2 objects.
    """
    
    components_v2 = []
    
    # Iterate through the V1 system components (stored in v1_systs._d)
    for v1_id, v1_syst in v1_systs._d.items():
        # NOTE: v1_syst is an instance of the V1 Syst class, which stores parameters in _pars
        
        # 1. Create the immutable ComponentV2 object
        system_v2 = ComponentDataV2(
            z=v1_syst._pars['z'].value,             # Extract value from Astropy Quantity
            dz=v1_syst._pars['dz'].value,           # Extract value
            logN=v1_syst._pars['logN'].value,
            dlogN=v1_syst._pars['dlogN'].value,
            b=v1_syst._pars['b'].value,
            db=v1_syst._pars['db'].value,
            btur=v1_syst._pars['btur'].value,
            dbtur=v1_syst._pars['dbtur'].value,
            func=v1_syst._func,
            series=v1_syst._series,
            # The 'id' is set automatically by ComponentV2.__post_init__ or will be assigned later
        )
        components_v2.append(system_v2)
        
    # 2. Create the V2 container
    # CRITICAL: We pass the mutable V1 lmfit models (v1_systs._mods_t) as a placeholder 
    # until the V2 ConstraintModelV2 is implemented.
    systlist_v2 = SystemListDataV2(
        components=components_v2,
        v1_models_t=v1_systs._mods_t, # Placeholder for V1 lmfit models
        meta=v1_systs._meta,
    )
    
    return systlist_v2