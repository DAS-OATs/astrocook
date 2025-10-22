from astropy.table import Table, Column
from astropy import units as au
import numpy as np
from typing import List, Any

from .structures import ComponentDataV2, SystemListDataV2

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