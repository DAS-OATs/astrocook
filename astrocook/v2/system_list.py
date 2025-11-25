from astropy.table import Table, Column
import astropy.units as au
import dataclasses
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
            
            # Create a "dummy" empty SystemListV2 to break the
            # recursion loop. We use __new__ to get an
            # un-initialized object, and then set its data manually.
            try:
                logging.warning("Creating fallback empty constraint model.")
                # 1. Create a dummy object *without* calling __init__
                dummy_systs = SystemListV2.__new__(SystemListV2)
                # 2. Set its data core manually
                dummy_systs._data = SystemListDataV2()
                dummy_systs.components = []
                # 3. Now, pass this safe object to the constructor
                self.constraint_model = VoigtModelConstraintV2(dummy_systs)
            except Exception as e_fallback:
                logging.error(f"FATAL: Could not create fallback constraint model: {e_fallback}")
                self.constraint_model = None # Give up
            

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
            # This must point to the *attribute* in the model,
            # not a property with the same name.
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
    
    def add_components_from_regions(self,
                                  spec: 'SpectrumV2',
                                  region_id_col: str,
                                  series_map: Dict[int, str],
                                  z_map: Dict[int, float]) -> 'SystemListV2':
        """
        Finds peaks in AOD within identified regions and adds them as
        new components to the system list.
        """
        from .spectrum_operations import find_signal_peaks # Local import for pureness

        if not spec.has_aux_column(region_id_col):
            logging.warning(f"add_components: '{region_id_col}' not in spectrum. No components added.")
            return self # Return self (no change)

        region_map = spec.get_column(region_id_col).value
        
        # 1. Get AOD profile for peak finding
        #cont = spec.cont.value if spec.cont is not None else np.ones_like(spec.y.value)
        #flux = spec.y.value
        #with np.errstate(divide='ignore', invalid='ignore'):
        #    aod = -np.log(flux / cont)
        #aod[~np.isfinite(aod)] = 0.0 # Clean NaNs/Infs

        new_components = list(self._data.components)
        
        # 2. Find next available component ID
        next_id = 1
        if new_components:
            next_id = max(c.id for c in new_components) + 1

        # 3. Replace peak finding with simple loop *** ---
        logging.info(f"Populating SystemList with {len(series_map)} components...")
        
        # 3. Loop through regions that have an ID
        for rid, series in series_map.items():
            if rid not in z_map:
                continue
            
            z_region = z_map[rid]
            
            # 4. Create ONE new ComponentDataV2 object per region
            new_comp = ComponentDataV2(
                id=next_id,
                z=z_region, # All components in region get same guess-z
                dz=0.001,
                logN=14.0,  # Default guess
                dlogN=0.5,
                b=10.0,     # Default guess
                db=5.0,
                btur=0.0,
                dbtur=None,
                func='voigt',
                series=series # Use the series ID (e.g., 'CIV_1548')
            )
            new_components.append(new_comp)
            next_id += 1

        # 6. Create new immutable SystemListDataV2
        new_data = dataclasses.replace(
            self._data, 
            components=new_components,
            # We'll need to handle constraints later
        )
        return SystemListV2(data=new_data)
    
    def add_component(self, 
                      series: str, 
                      z: float, 
                      logN: float = 13.5, 
                      b: float = 10.0, 
                      btur: float = 0.0) -> 'SystemListV2':
        """
        API: Adds a single component to the list.
        Returns a NEW SystemListV2 instance.
        """
        # 1. Get existing components
        current_components = self._data.components
        
        # 2. Determine new ID (V1 integer compatibility)
        next_id = 1
        if current_components:
            next_id = max(c.id for c in current_components) + 1
            
        # 3. Create the immutable Component object
        new_comp = ComponentDataV2(
            id=next_id,
            z=z,
            dz=None,
            logN=logN,
            dlogN=None,
            b=b,
            db=None,
            btur=btur,
            dbtur=None,
            func='voigt',
            series=series
        )
        
        # 4. Create the new list
        new_component_list = current_components + [new_comp]
        
        # 5. Create new Data Core (preserving meta/constraints)
        # We need to update the UUID map for the new component
        new_id_map = dict(self._data.v1_id_to_uuid_map)
        new_id_map[next_id] = new_comp.uuid
        
        new_data_core = dataclasses.replace(
            self._data,
            components=new_component_list,
            v1_id_to_uuid_map=new_id_map
        )
        
        # 6. Return new API Object
        return SystemListV2(new_data_core)
    
    def get_component_by_uuid(self, uuid: str) -> Optional[ComponentDataV2]:
        for c in self.components:
            if c.uuid == uuid:
                return c
        return None

    def update_component(self, uuid: str, **changes) -> 'SystemListV2':
        """
        API: Returns a new SystemListV2 with one component modified.
        """
        new_components = []
        found = False
        for c in self.components:
            if c.uuid == uuid:
                # Apply changes using dataclasses.replace
                # We need to filter 'changes' to match field names (e.g. remove 'id')
                valid_fields = {f.name for f in dataclasses.fields(c)}
                filtered_changes = {k: v for k, v in changes.items() if k in valid_fields}
                
                new_c = dataclasses.replace(c, **filtered_changes)
                new_components.append(new_c)
                found = True
            else:
                new_components.append(c)
        
        if not found:
            logging.warning(f"Component {uuid} not found for update.")
            return self

        # Return new container
        new_data = dataclasses.replace(self._data, components=new_components)
        return SystemListV2(new_data)

    def delete_component(self, uuid: str) -> 'SystemListV2':
        """
        API: Returns a new SystemListV2 with the component removed.
        """
        new_components = [c for c in self.components if c.uuid != uuid]
        
        # Return new container
        new_data = dataclasses.replace(self._data, components=new_components)
        return SystemListV2(new_data)