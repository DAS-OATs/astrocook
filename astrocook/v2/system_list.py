from astropy.table import Table, Column
import astropy.units as au
import numpy as np
from typing import List

from .structures import ComponentDataV2, SystemListDataV2
from .system_list_migration import migrate_component_v2_to_v1
from ..v1.syst_list import SystList as SystListV1

# --- 3. System List API Layer (The orchestrator) ---
class SystemListV2:
    """API layer for managing component lists, constraints, and fitting."""
    
    def __init__(self, data: SystemListDataV2):
        self._data = data
        # NOTE: ConstraintModelV2 will be instantiated here later
        # self.constraints = ConstraintModelV2(data) 

    @property
    def components(self) -> List[ComponentDataV2]:
        return self._data.components

    @property
    def t(self) -> Table:
        """
        Adapter: Dynamically generates an Astropy Table (V1-compatible) 
        from the immutable list of ComponentDataV2 objects.
        """
        
        # 1. Initialize lists to hold data from all components
        z_list, dz_list, logN_list, dlogN_list, b_list, db_list = ([] for _ in range(6))
        btur_list, dbtur_list, func_list, series_list, id_list = ([] for _ in range(5))
        
        # 2. Extract data from immutable ComponentDataV2 objects
        for c in self._data.components: # Access components list from the Data Core
            z_list.append(c.z)
            dz_list.append(c.dz if c.dz is not None else np.nan)
            logN_list.append(c.logN)
            dlogN_list.append(c.dlogN if c.dlogN is not None else np.nan)
            b_list.append(c.b)
            db_list.append(c.db if c.db is not None else np.nan)
            btur_list.append(c.btur)
            dbtur_list.append(c.dbtur if c.dbtur is not None else np.nan)
            func_list.append(c.func)
            series_list.append(c.series)
            id_list.append(c.id)

        # 3. Build the Astropy Table (V1 structure)
        t = Table()
        t['z'] = Column(z_list, unit=au.dimensionless_unscaled)
        t['dz'] = Column(dz_list, unit=au.dimensionless_unscaled)
        t['logN'] = Column(logN_list, unit=au.dimensionless_unscaled)
        t['dlogN'] = Column(dlogN_list, unit=au.dimensionless_unscaled)
        t['b'] = Column(b_list, unit=au.km/au.s)
        t['db'] = Column(db_list, unit=au.km/au.s)
        t['btur'] = Column(btur_list, unit=au.km/au.s)
        t['dbtur'] = Column(dbtur_list, unit=au.km/au.s)
        t['func'] = Column(func_list)
        t['series'] = Column(series_list)
        t['id'] = Column(id_list)

        return t

    @property
    def z(self) -> List[float]:
        """Returns a flat list of component redshifts (required for plotting)."""
        return [c.z for c in self._data.components]

    @property
    def series(self) -> List[str]:
        """Returns a flat list of component series (required for plotting)."""
        return [c.series for c in self._data.components]
        
    @property
    def id(self) -> List[int]:
        """Returns a flat list of component IDs."""
        return [c.id for c in self._data.components]

    # TODO: Methods for constraints, fitting, and immutable updates will go here.

    def to_v1_systlist(self) -> SystListV1:
        """Converts the immutable SystemListV2 back to a mutable SystListV1 for saving."""
        
        # 1. Convert each V2 component into a V1 Syst object
        v1_syst_objects = []
        for component_data in self._data.components:
            # We assume migrate_component_v2_to_v1 returns a V1 Syst object
            v1_syst = migrate_component_v2_to_v1(component_data)
            v1_syst_objects.append(v1_syst)
            
        # 2. Extract necessary parallel arrays for V1 SystList constructor
        # NOTE: V1 SystList often requires constructing parallel arrays/columns from the V1 Syst objects.
        # This is complex and relies on the V1 SystList constructor signature.
        
        # SIMPLIFIED V2 APPROACH: Construct a list of dictionaries that the V1 SystList constructor can parse.
        # If the V1 SystList constructor is robust, it can take the raw data arrays.
        components = self._data.components
        z_list = [c.z for c in components]
        dz_list = [c.dz if c.dz is not None else np.nan for c in components]
        logN_list = [c.logN for c in components]
        dlogN_list = [c.dlogN if c.dlogN is not None else np.nan for c in components]
        b_list = [c.b for c in components]
        db_list = [c.db if c.db is not None else np.nan for c in components]
        btur_list = [c.btur for c in components]
        dbtur_list = [c.dbtur if c.dbtur is not None else np.nan for c in components]
        
        func_list = [c.func for c in components]
        series_list = [c.series for c in components]
        id_list = [c.id for c in components]
        
        # NOTE: V1 SystList constructor also expects arrays for state/results:
        # resol, chi2r, snr. We must provide placeholder arrays with the correct length.
        length = len(components)
        nan_array = [np.nan] * length
        
        # 2. Instantiate the V1 SystList with *all* parameters explicitly passed
        v1_systs = SystListV1(
            # Core Parameters
            z=z_list, dz=dz_list, 
            logN=logN_list, dlogN=dlogN_list, 
            b=b_list, db=db_list, 
            btur=btur_list, dbtur=dbtur_list,
            # State/Metadata Parameters (required by V1 __init__)
            func=func_list,
            series=series_list,
            id=id_list,
            resol=nan_array,
            chi2r=nan_array,
            snr=nan_array
            # Other potential parameters from the V1 constructor (e.g., id_start, meta) must also be checked
        )
        
        # 4. Attach the V1 lmfit models placeholder (critical for V1 saving logic)
        v1_systs._mods_t = self._data.v1_models_t
        
        return v1_systs