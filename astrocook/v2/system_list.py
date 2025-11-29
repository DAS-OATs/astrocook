import ast
from astropy.table import Table, Column
import astropy.units as au
import dataclasses
import logging
import numpy as np
from typing import Any, Dict, List, Optional, Set, Tuple

from .atomic_data import ATOM_DATA, STANDARD_MULTIPLETS
from .fitting.voigt_model import VoigtModelConstraintV2
from .structures import ComponentDataV2, SystemListDataV2
from .system_list_migration import migrate_component_v2_to_v1
from ..v1.syst_list import SystList as SystListV1

# --- 3. System List API Layer (The orchestrator) ---
class SystemListV2:
    """API layer for managing component lists, constraints, and fitting."""
    
    def __init__(self, data: SystemListDataV2):
        self._data = data
        
        # Sanitize data immediately (Fix JSON artifacts)
        # This ensures self._data always has proper Tuples/Ints as keys
        self._ensure_clean_data()

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
    
    def get_connected_group(self, seed_uuids: List[str]) -> Set[str]:
        """
        Returns a set of UUIDs that are 'Fluid Neighbors' of the seed components.
        
        Logic: 
        - Shallow Grouping (Seed + Immediate Neighbors only).
        - Dynamic Threshold: 4.0 * max(b_seed, b_candidate) (Min 30 km/s).
        - Multiplet-Aware: Checks ALL lines in a series for overlaps (e.g. 1548 vs 1550).
        """
        if not seed_uuids:
            return set()
        
        c_kms = 299792.458
        
        # 1. Identify Seed Objects
        seeds = []
        comp_map = {c.uuid: c for c in self.components}
        
        for uuid in seed_uuids:
            if uuid in comp_map:
                seeds.append(comp_map[uuid])
        
        group_uuids = set(seed_uuids)

        # Helper: Get ALL rest wavelengths for a series
        # Returns a list, e.g., [1548.2, 1550.7] for 'CIV'
        def get_all_w_rest(series_name: str) -> List[float]:
            # 1. Check if it's a specific line (e.g., 'CIV_1548') in ATOM_DATA
            if series_name in ATOM_DATA:
                return [ATOM_DATA[series_name]['wave']]
            
            # 2. Check if it's a Multiplet group (e.g., 'CIV') in STANDARD_MULTIPLETS
            if series_name in STANDARD_MULTIPLETS:
                waves = []
                for line_name in STANDARD_MULTIPLETS[series_name]:
                    if line_name in ATOM_DATA:
                        waves.append(ATOM_DATA[line_name]['wave'])
                if waves: return waves

            return []

        for seed in seeds:
            # Prepare Seed Wavelengths (List of all potential lines)
            seed_waves_rest = get_all_w_rest(seed.series)
            if not seed_waves_rest:
                logging.warning(f"[GROUP] No atomic data for {seed.series}")
            
            seed_z = seed.z
            seed_b = seed.b
            
            # Calculate observed wavelengths for ALL lines in the seed multiplet
            seed_waves_obs = [w * (1 + seed_z) for w in seed_waves_rest]

            for other in self.components:
                if other.uuid in group_uuids: continue 
                
                # --- THRESHOLD CALCULATION ---
                # 4.0 sigma + 30 km/s floor
                threshold_kms = max(4.0 * max(seed_b, other.b), 30.0)

                is_connected = False
                
                # --- Criterion 1: Kinematic Proximity (Same System) ---
                # Simple redshift check
                dz = abs(seed_z - other.z)
                dv_kin = c_kms * dz / (1 + seed_z)
                
                if dv_kin < threshold_kms:
                    is_connected = True
                
                # --- Criterion 2: Spectral Overlap (Blends) ---
                # Check EVERY line in Seed vs EVERY line in Other
                dv_blend_min = 99999.0
                
                if not is_connected and seed_waves_obs:
                    other_waves_rest = get_all_w_rest(other.series)
                    other_waves_obs = [w * (1 + other.z) for w in other_waves_rest]
                    
                    # Nested loop: Does ANY line of Seed overlap ANY line of Other?
                    for w1_obs in seed_waves_obs:
                        for w2_obs in other_waves_obs:
                            lam_avg = (w1_obs + w2_obs) / 2.0
                            dv = c_kms * abs(w1_obs - w2_obs) / lam_avg
                            
                            if dv < dv_blend_min:
                                dv_blend_min = dv
                            
                            if dv < threshold_kms:
                                is_connected = True
                                break # Found an overlap
                        if is_connected: break

                # --- DEBUG LOGGING ---
                # Log if it's a match OR if it's a close miss (e.g. within 200 km/s)
                debug_condition = is_connected or (dv_kin < 200.0) or (dv_blend_min < 200.0)
                
                if debug_condition:
                    status = "MATCH" if is_connected else "SKIP"
                    logging.debug(
                        f"[GROUP] {seed.series}({seed.z:.4f}) vs {other.series}({other.z:.4f}) | "
                        f"Limit: {threshold_kms:.1f} | dV_kin: {dv_kin:.1f} | dV_blend: {dv_blend_min:.1f} -> {status}"
                    )

                if is_connected:
                    group_uuids.add(other.uuid)
                    
        return group_uuids
    
    def _ensure_clean_data(self):
        """
        Robustly initializes data core.
        1. Sanitize artifacts in ID map.
        2. [NEW] Recover constraints from JSON-friendly 'v2_constraints_map' if present.
        3. [NEW] Backfill 'v2_constraints_map' if missing (so next save is clean).
        """
        dirty = False
        
        # A. Fix ID Map (Str "123" -> Int 123)
        id_map = self._data.v1_id_to_uuid_map
        clean_id_map = {}
        for k, v in id_map.items():
            if isinstance(k, str) and k.isdigit():
                clean_id_map[int(k)] = v
                dirty = True
            else:
                clean_id_map[k] = v
        
        # B. Robust Constraint Recovery
        # Strategy: Trust 'v2_constraints_map' (clean JSON) over 'parsed_constraints' (tuple keys)
        # if the former exists.
        
        clean_parsed = dict(self._data.parsed_constraints)
        clean_v2_map = dict(self._data.v2_constraints_map)
        constraints_dirty = False

        # 1. If we have the Clean Map (from a good save), rebuild Parsed Map from it.
        if clean_v2_map:
            # Rebuild parsed constraints from the clean nested source
            # This completely bypasses the "Stringified Tuple" issue.
            new_parsed = {}
            for u, p_map in clean_v2_map.items():
                for p, c in p_map.items():
                    new_parsed[(u, p)] = c
            
            # If they differ, update parsed (Logic: Clean Map is the Master)
            if len(new_parsed) != len(clean_parsed): # Simple check, could be more thorough
                clean_parsed = new_parsed
                constraints_dirty = True
                logging.info("SystemListV2: Restored constraints from clean JSON map.")
        
        # 2. If Clean Map is empty but Parsed Map exists (Legacy Load), 
        #    Sanitize Parsed Map and Populated Clean Map.
        elif clean_parsed:
            sanitized_parsed = {}
            new_v2_map = {}
            
            for k, v in clean_parsed.items():
                final_key = k
                
                # Fix Stringified Tuples "('u', 'z')"
                if isinstance(k, str) and k.startswith('('):
                    try:
                        val = ast.literal_eval(k)
                        if isinstance(val, tuple): final_key = val
                    except: pass
                
                # Fix Stringified Ints in Tuples ("123", 'z')
                if isinstance(final_key, tuple) and len(final_key) == 2:
                    p0, p1 = final_key
                    if isinstance(p0, str) and p0.isdigit():
                        final_key = (int(p0), p1)
                
                sanitized_parsed[final_key] = v
                
                # Populate the V2 Nested Map (uuid -> param -> constraint)
                # This ensures the NEXT save will be clean.
                uuid_key, param_key = final_key
                # Note: If key is (int_id, param), we can't put it in v2_map (needs uuid).
                # We skip legacy ID keys for the v2 map unless we can resolve the UUID.
                # (Skipping legacy backfill complexity for now to keep it safe)
                if isinstance(uuid_key, str): # valid UUID
                    if uuid_key not in new_v2_map: new_v2_map[uuid_key] = {}
                    new_v2_map[uuid_key][param_key] = v

            clean_parsed = sanitized_parsed
            clean_v2_map = new_v2_map
            constraints_dirty = True
            logging.info("SystemListV2: Sanitized legacy constraints and populated clean map.")

        if dirty or constraints_dirty:
            self._data = dataclasses.replace(
                self._data, 
                v1_id_to_uuid_map=clean_id_map,
                parsed_constraints=clean_parsed,
                v2_constraints_map=clean_v2_map
            )

    def update_constraint(self, uuid: str, param: str, 
                          is_free: Optional[bool] = None,
                          expression: Optional[str] = None) -> 'SystemListV2':
        """
        API: Updates the constraint for a specific component parameter.
        Updates BOTH the parsed_constraints (tuple key) and v2_constraints_map (clean structure).
        """
        current_parsed = dict(self._data.parsed_constraints)
        
        # 1. Determine Key (Use UUID for V2)
        key = (uuid, param)
        
        # 2. Retrieve Existing State (Check V2 key, then V1 ID key)
        existing = current_parsed.get(key)
        if not existing:
             comp = self.get_component_by_uuid(uuid)
             if comp:
                 legacy_key = (comp.id, param)
                 existing = current_parsed.get(legacy_key)
        
        from .structures import ParameterConstraintV2
        
        # 3. Handle Dict (Legacy) vs Object (V2)
        if existing:
            if isinstance(existing, dict):
                 old_free = existing.get('is_free', True)
                 old_expr = existing.get('expression')
                 base_constraint = ParameterConstraintV2(is_free=old_free, expression=old_expr)
            else:
                 base_constraint = existing

            new_free = is_free if is_free is not None else base_constraint.is_free
            new_expr = expression if expression is not None else base_constraint.expression
            new_constraint = dataclasses.replace(base_constraint, is_free=new_free, expression=new_expr)
        else:
            new_free = is_free if is_free is not None else True
            new_constraint = ParameterConstraintV2(is_free=new_free, expression=expression)

        # 4. Update Parsed Map (Tuple Key)
        current_parsed[key] = new_constraint
        
        # 5. [NEW] Update V2 Nested Map (Clean JSON Structure)
        # We assume the user saving logic will pick this up.
        current_v2_map = {k: v.copy() for k, v in self._data.v2_constraints_map.items()}
        if uuid not in current_v2_map:
            current_v2_map[uuid] = {}
        current_v2_map[uuid][param] = new_constraint
        
        new_data = dataclasses.replace(self._data, 
                                       parsed_constraints=current_parsed,
                                       v2_constraints_map=current_v2_map)
        
        return SystemListV2(new_data)