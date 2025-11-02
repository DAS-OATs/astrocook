from astropy.table import Table, Column
import astropy.units as au
import logging
import numpy as np
from typing import Any, Dict, List, Optional, Tuple

from .fitting.voigt_model import VoigtModelConstraintV2
from .structures import ComponentDataV2, SystemListDataV2
from .system_list_migration import migrate_component_v2_to_v1
from ..v1.syst_list import SystList as SystListV1

# --- 3. System List API Layer (The orchestrator) ---
class SystemListV2:
    """API layer for managing component lists, constraints, and fitting."""
    
    def __init__(self, data: SystemListDataV2):
        self._data = data
        
        # 1. Initialize the attribute to None first
        self.constraint_model = None 
        
        # 2. CRITICAL STEP: Call the initialization logic now
        self._initialize_constraint_model()

    def _initialize_constraint_model(self):
        """Initializes the constraint model once the SystemListV2 object is constructed."""
        
        # The VoigtModelConstraintV2 constructor
        # takes the *SystemListV2* instance (self) and reads the
        # data core (_data) from it to build the constraints.
        #
        # We must pass 'self' to the constructor.
        
        try:
            self.constraint_model = VoigtModelConstraintV2(self)
        except Exception as e:
            logging.error(f"Failed to initialize VoigtModelConstraintV2: {e}", exc_info=True)
            # Fallback to an empty model if initialization fails
            if self.constraint_model is None:
                 # Create an empty model if it failed during init
                 self.constraint_model = VoigtModelConstraintV2(SystemListV2(data=SystemListDataV2()))

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

    @property
    def v1_models_t(self) -> Any:
        """Adapter for temporary read-only access to the V1 lmfit models placeholder (_mods_t)."""
        # CRITICAL FIX: Return the attribute from the underlying data core
        return self._data.v1_models_t

    @property
    def constraint_status(self) -> Dict[str, np.ndarray]:
        """
        Returns the constraint status arrays (is_free, link_target_index) 
        for GUI display (e.g., coloring cells).
        """
        return self.constraint_model._param_map
    
    @property
    def parsed_constraints(self) -> Dict[Tuple[int, str], Dict[str, Any]]:
        """Provides read-only access to the V2-parsed constraint state."""
        return self._data.parsed_constraints
    
    @property
    def v2_constraints_for_save(self) -> Dict[str, Dict[str, Any]]:
        """
        Returns the V2-native constraint map (UUID-keyed) for serialization.
        """
        if self.constraint_model:
            return self.constraint_model.v2_constraints_by_uuid
        return {}

    def to_v1_systlist(self) -> Optional[Table]:
        """
        Converts the V2 component data back into a V1-style Astropy Table (SystList).
        
        NOTE: This is a placeholder. The full implementation needs to reconstruct
        the complex V1 Table format, including metadata and constraints in the header.
        """
        if not self.components:
            return None # Return None if empty

        # This is a simplified conversion, focusing on core data.
        # A full V1 save would require rebuilding the exact V1 Table structure.
        logging.warning("Simplified V2->V1 system list conversion. V1 constraints may be lost.")
        
        data_dict = {
            'id': [c.id for c in self.components],
            'z': [c.z for c in self.components],
            'dz': [c.dz if c.dz is not None else 0.0 for c in self.components],
            'logN': [c.logN for c in self.components],
            'dlogN': [c.dlogN if c.dlogN is not None else 0.0 for c in self.components],
            'b': [c.b for c in self.components],
            'db': [c.db if c.db is not None else 0.0 for c in self.components],
            'btur': [c.btur for c in self.components],
            'dbtur': [c.dbtur if c.dbtur is not None else 0.0 for c in self.components],
            'func': [c.func for c in self.components],
            'series': [c.series for c in self.components],
        }
        
        v1_table = Table(data_dict)
        # TODO: Add V1 constraint logic (e.g., self.constraint_model.to_v1_header())
        
        return v1_table