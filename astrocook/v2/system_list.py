from astropy.table import Table, Column
import astropy.units as au
import dataclasses
import logging
import numpy as np
import ast
from typing import Any, Dict, List, Optional, Set, Tuple

from .fitting.voigt_model import VoigtModelConstraintV2
from .structures import ComponentDataV2, SystemListDataV2, ParameterConstraintV2
from .system_list_migration import migrate_component_v2_to_v1
from ..legacy.syst_list import SystList as SystListV1
from .atomic_data import ATOM_DATA, STANDARD_MULTIPLETS

# --- 3. System List API Layer (The orchestrator) ---
class SystemListV2:
    """API layer for managing component lists, constraints, and fitting."""
    
    def __init__(self, data: SystemListDataV2):
        self._data = data
        
        # 1. Sanitize & Hydrate data immediately
        # (Fixes JSON artifacts and converts dicts -> objects)
        self._ensure_clean_data()
        
        self.constraint_model = None 
        self._initialize_constraint_model()

    def _ensure_clean_data(self):
        """
        Robustly initializes data core.
        1. Sanitize artifacts in ID map.
        2. HYDRATE 'v2_constraints_map' (Dict -> ParameterConstraintV2).
        3. Recover parsed constraints from the clean map if present.
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
        
        # B. Robust Constraint Recovery & Hydration
        clean_parsed = dict(self._data.parsed_constraints)
        raw_v2_map = dict(self._data.v2_constraints_map)
        hydrated_v2_map = {}
        constraints_dirty = False

        # 1. Hydrate V2 Map (Convert Dicts back to Objects)
        if raw_v2_map:
            for u, p_map in raw_v2_map.items():
                hydrated_v2_map[u] = {}
                for p, c in p_map.items():
                    if isinstance(c, dict):
                        # Convert Dict -> Dataclass
                        hydrated_v2_map[u][p] = ParameterConstraintV2(
                            is_free=c.get('is_free', True),
                            target_uuid=c.get('target_uuid'),
                            expression=c.get('expression'),
                            v1_target_id=c.get('v1_target_id')
                        )
                        constraints_dirty = True
                    else:
                        hydrated_v2_map[u][p] = c
        
        # 2. Rebuild Parsed Constraints from Clean Source (if available)
        if hydrated_v2_map:
            new_parsed = {}
            for u, p_map in hydrated_v2_map.items():
                for p, c in p_map.items():
                    new_parsed[(u, p)] = c
            
            # Use the clean rebuilt map
            clean_parsed = new_parsed
            constraints_dirty = True
            logging.info("SystemListV2: Hydrated constraints from clean JSON map.")
        
        # 3. Legacy Fallback: If V2 map empty, Sanitize Parsed Map & Populate V2
        elif clean_parsed:
            sanitized_parsed = {}
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
                
                # Ensure value is an object (defensive)
                val_obj = v
                if isinstance(v, dict):
                    val_obj = ParameterConstraintV2(
                        is_free=v.get('is_free', True),
                        target_uuid=v.get('target_uuid'),
                        expression=v.get('expression')
                    )
                    sanitized_parsed[final_key] = val_obj

                # Populate V2 Map for next save
                uuid_key, param_key = final_key
                if isinstance(uuid_key, str): 
                    if uuid_key not in hydrated_v2_map: hydrated_v2_map[uuid_key] = {}
                    hydrated_v2_map[uuid_key][param_key] = val_obj

            clean_parsed = sanitized_parsed
            constraints_dirty = True
            logging.info("SystemListV2: Sanitized legacy constraints and populated clean map.")

        if dirty or constraints_dirty:
            self._data = dataclasses.replace(
                self._data, 
                v1_id_to_uuid_map=clean_id_map,
                parsed_constraints=clean_parsed,
                v2_constraints_map=hydrated_v2_map
            )

    def _initialize_constraint_model(self):
        try:
            self.constraint_model = VoigtModelConstraintV2(self)
        except Exception as e:
            logging.error(f"Failed to initialize VoigtModelConstraintV2: {e}", exc_info=True)
            try:
                logging.warning("Creating fallback empty constraint model.")
                dummy_systs = SystemListV2.__new__(SystemListV2)
                dummy_systs._data = SystemListDataV2()
                dummy_systs.components = []
                self.constraint_model = VoigtModelConstraintV2(dummy_systs)
            except Exception:
                self.constraint_model = None

    @property
    def components(self) -> List[ComponentDataV2]:
        return self._data.components

    @property
    def t(self) -> Table:
        z_list, dz_list, logN_list, dlogN_list, b_list, db_list = ([] for _ in range(6))
        btur_list, dbtur_list, func_list, series_list, id_list = ([] for _ in range(5))
        
        for c in self._data.components: 
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
        return [c.z for c in self._data.components]

    @property
    def series(self) -> List[str]:
        return [c.series for c in self._data.components]
        
    @property
    def id(self) -> List[int]:
        return [c.id for c in self._data.components]

    @property
    def v1_models_t(self) -> Any:
        return self._data.v1_models_t

    @property
    def constraint_status(self) -> Dict[str, np.ndarray]:
        return self.constraint_model._param_map
    
    @property
    def parsed_constraints(self) -> Dict[Tuple[int, str], Dict[str, Any]]:
        return self._data.parsed_constraints
    
    @property
    def v2_constraints_for_save(self) -> Dict[str, Dict[str, Any]]:
        if self.constraint_model:
            return self.constraint_model.v2_constraints_by_uuid
        return {}

    def to_v1_systlist(self) -> Optional[Table]:
        if not self.components: return None
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
        return Table(data_dict)
    
    def add_components_from_regions(self, spec: 'SpectrumV2', region_id_col: str,
                                  series_map: Dict[int, str], z_map: Dict[int, float]) -> 'SystemListV2':
        if not spec.has_aux_column(region_id_col):
            logging.warning(f"add_components: '{region_id_col}' not in spectrum. No components added.")
            return self

        new_components = list(self._data.components)
        next_id = 1
        if new_components:
            next_id = max(c.id for c in new_components) + 1

        logging.info(f"Populating SystemList with {len(series_map)} components...")
        
        for rid, series in series_map.items():
            if rid not in z_map: continue
            z_region = z_map[rid]
            new_comp = ComponentDataV2(
                id=next_id,
                z=z_region, 
                dz=0.001,
                logN=14.0, 
                dlogN=0.5,
                b=10.0, 
                db=5.0,
                btur=0.0,
                dbtur=None,
                func='voigt',
                series=series
            )
            new_components.append(new_comp)
            next_id += 1

        new_data = dataclasses.replace(self._data, components=new_components)
        return SystemListV2(data=new_data)
    
    def add_component(self, series: str, z: float, logN: float = 13.5, 
                      b: float = 10.0, btur: float = 0.0) -> 'SystemListV2':
        current_components = self._data.components
        next_id = 1
        if current_components:
            next_id = max(c.id for c in current_components) + 1
            
        new_comp = ComponentDataV2(
            id=next_id,
            z=z, dz=None,
            logN=logN, dlogN=None,
            b=b, db=None,
            btur=btur, dbtur=None,
            func='voigt', series=series
        )
        
        new_component_list = current_components + [new_comp]
        new_id_map = dict(self._data.v1_id_to_uuid_map)
        new_id_map[next_id] = new_comp.uuid
        
        new_data_core = dataclasses.replace(
            self._data,
            components=new_component_list,
            v1_id_to_uuid_map=new_id_map
        )
        return SystemListV2(new_data_core)
    
    def get_component_by_uuid(self, uuid: str) -> Optional[ComponentDataV2]:
        for c in self.components:
            if c.uuid == uuid: return c
        return None

    def update_component(self, uuid: str, **changes) -> 'SystemListV2':
        new_components = []
        found = False
        for c in self.components:
            if c.uuid == uuid:
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

        new_data = dataclasses.replace(self._data, components=new_components)
        return SystemListV2(new_data)

    

    def update_constraint(self, uuid: str, param: str, 
                          is_free: Optional[bool] = None,
                          expression: Optional[str] = None,
                          target_uuid: Optional[str] = None) -> 'SystemListV2':
        """
        API: Updates the constraint for a specific component parameter.
        """
        current_parsed = dict(self._data.parsed_constraints)
        
        # 1. Resolve Existing
        key = (uuid, param)
        existing = current_parsed.get(key)
        if not existing:
             comp = self.get_component_by_uuid(uuid)
             if comp:
                 legacy_key = (comp.id, param)
                 existing = current_parsed.get(legacy_key)
        
        from .structures import ParameterConstraintV2
        
        # 2. Determine Base
        if existing:
            if isinstance(existing, dict):
                 base_constraint = ParameterConstraintV2(
                     is_free=existing.get('is_free', True), 
                     expression=existing.get('expression'),
                     target_uuid=existing.get('target_uuid')
                 )
            else:
                 base_constraint = existing
        else:
            base_constraint = ParameterConstraintV2(is_free=True)

        # 3. Apply Changes
        # [FIX] Logic for Unlinking: If we explicitly set is_free=True, 
        # we must wipe the expression and target, regardless of what was passed.
        if is_free is True:
            new_free = True
            new_expr = None
            new_tgt = None
        else:
            new_free = is_free if is_free is not None else base_constraint.is_free
            new_expr = expression if expression is not None else base_constraint.expression
            new_tgt = target_uuid if target_uuid is not None else base_constraint.target_uuid

        new_constraint = dataclasses.replace(base_constraint, is_free=new_free, expression=new_expr, target_uuid=new_tgt)

        # 4. Update Maps
        current_parsed[key] = new_constraint
        
        current_v2_map = {k: v.copy() for k, v in self._data.v2_constraints_map.items()}
        if uuid not in current_v2_map:
            current_v2_map[uuid] = {}
        current_v2_map[uuid][param] = new_constraint
        
        new_data = dataclasses.replace(self._data, 
                                       parsed_constraints=current_parsed,
                                       v2_constraints_map=current_v2_map)
        
        return SystemListV2(new_data)

    def delete_component(self, uuid: str) -> 'SystemListV2':
        # 1. Remove the component
        new_components = [c for c in self.components if c.uuid != uuid]
        
        # 2. [FIX] Clean up Hanging Links (Sanitize Constraints)
        # Scan all constraints. If any target the deleted UUID, reset them to Free.
        from .structures import ParameterConstraintV2
        
        current_v2_map = {k: v.copy() for k, v in self._data.v2_constraints_map.items()}
        current_parsed = dict(self._data.parsed_constraints)
        
        # We need to find which (u, p) pairs need resetting
        to_reset = []
        
        for u_src, params in current_v2_map.items():
            for p_name, constr in params.items():
                # Check for explicit target_uuid match
                if constr.target_uuid == uuid:
                    to_reset.append((u_src, p_name))
                # Check for expression reference (regex search)
                elif constr.expression and uuid in constr.expression:
                    to_reset.append((u_src, p_name))

        for u_src, p_name in to_reset:
            logging.info(f"Removing broken link from {u_src}.{p_name} (target deleted).")
            # Reset to Free
            clean_constr = ParameterConstraintV2(is_free=True, expression=None, target_uuid=None)
            
            # Update Nested Map
            current_v2_map[u_src][p_name] = clean_constr
            
            # Update Parsed Map
            current_parsed[(u_src, p_name)] = clean_constr

        # 3. Create new data
        new_data = dataclasses.replace(self._data, 
                                       components=new_components,
                                       parsed_constraints=current_parsed,
                                       v2_constraints_map=current_v2_map)
        return SystemListV2(new_data)

    def get_connected_group(self, seed_uuids: List[str]) -> Set[str]:
        """
        Returns a set of UUIDs that are 'Fluid Neighbors' of the seed components.
        Logic: Linkage -> Spectral Blends -> Kinematic Proximity (Same Species).
        """
        if not seed_uuids:
            return set()
        
        c_kms = 299792.458
        seeds = []
        comp_map = {c.uuid: c for c in self.components}
        
        for uuid in seed_uuids:
            if uuid in comp_map:
                seeds.append(comp_map[uuid])
        
        group_uuids = set(seed_uuids)

        def get_all_w_rest(series_name: str) -> List[float]:
            if series_name in ATOM_DATA:
                return [ATOM_DATA[series_name]['wave']]
            if series_name in STANDARD_MULTIPLETS:
                waves = []
                for line_name in STANDARD_MULTIPLETS[series_name]:
                    if line_name in ATOM_DATA:
                        waves.append(ATOM_DATA[line_name]['wave'])
                if waves: return waves
            return []

        def are_same_species(s1: str, s2: str) -> bool:
            if s1 == s2: return True
            parent1, parent2 = None, None
            for k, members in STANDARD_MULTIPLETS.items():
                if s1 in members: parent1 = k
                if s2 in members: parent2 = k
            if parent1 and parent2 and parent1 == parent2: return True
            if s1 in STANDARD_MULTIPLETS and s2 in STANDARD_MULTIPLETS[s1]: return True
            if s2 in STANDARD_MULTIPLETS and s1 in STANDARD_MULTIPLETS[s2]: return True
            return False

        constraints_map = self._data.v2_constraints_map

        for seed in seeds:
            seed_waves_rest = get_all_w_rest(seed.series)
            seed_z = seed.z
            seed_b = seed.b
            seed_waves_obs = [w * (1 + seed_z) for w in seed_waves_rest]

            for other in self.components:
                if other.uuid in group_uuids: continue 
                
                is_connected = False
                
                # 1. Parameter Links
                if seed.uuid in constraints_map:
                    for constr in constraints_map[seed.uuid].values():
                        if constr.target_uuid == other.uuid:
                            is_connected = True; break
                
                if not is_connected and other.uuid in constraints_map:
                    for constr in constraints_map[other.uuid].values():
                        if constr.target_uuid == seed.uuid:
                            is_connected = True; break
                
                if is_connected:
                    group_uuids.add(other.uuid); continue

                threshold_kms = max(4.0 * max(seed_b, other.b), 30.0)

                # 2. Spectral Overlap
                if seed_waves_obs:
                    other_waves_rest = get_all_w_rest(other.series)
                    other_waves_obs = [w * (1 + other.z) for w in other_waves_rest]
                    
                    for w1_obs in seed_waves_obs:
                        for w2_obs in other_waves_obs:
                            lam_avg = (w1_obs + w2_obs) / 2.0
                            dv = c_kms * abs(w1_obs - w2_obs) / lam_avg
                            if dv < threshold_kms:
                                is_connected = True; break
                        if is_connected: break
                
                if is_connected:
                    group_uuids.add(other.uuid); continue

                # 3. Kinematic Proximity (Same Species)
                if are_same_species(seed.series, other.series):
                    dz = abs(seed_z - other.z)
                    dv_kin = c_kms * dz / (1 + seed_z)
                    if dv_kin < threshold_kms:
                        is_connected = True
                        
                if is_connected:
                    group_uuids.add(other.uuid)
                    
        return group_uuids