from astropy.table import Table
import logging
import numpy as np
import re
from typing import List, Dict, Any, TYPE_CHECKING, Optional, Set, Tuple

from ..structures import ComponentDataV2, ParameterConstraintV2 

if TYPE_CHECKING:
    from ..system_list import SystemListV2

class VoigtModelConstraintV2:
    """
    Manages parameter constraints and converts component data into vectors.
    Supports "Selective Fitting" with Neighbor Activation.
    """
    
    PARAMETER_ORDER = ['z', 'logN', 'b', 'btur'] 
    V1_ID_REGEX = re.compile(r"p\[(\d+)\]")

    # Default Physical Limits (Min, Max)
    DEFAULT_BOUNDS = {
        'z':    (-0.5, 10.0),
        'logN': (10.0, 23.0),
        'b':    (1.0, 300.0),
        'btur': (0.0, 100.0)
    }

    def __init__(self, system_list: 'SystemListV2', constraints: Dict[str, Any] = None):
        self._system_list = system_list
        self.v2_constraints_by_uuid: Dict[str, Dict[str, ParameterConstraintV2]] = {}          

        self._active_uuids: Optional[Set[str]] = None
        self._component_groups = [list(self._system_list.components)]
        self._refresh_state()
        
        logging.info(f"Initialized VoigtModelConstraintV2 with {len(self.p_free_vector)} free parameters.")

    def set_active_components(self, target_uuids: List[str] = None):
        """
        Defines the "Fluid Group". 
        Uses the shared logic from SystemListV2 (Kinematic + Blends, Shallow).
        """
        if target_uuids is None:
            self._active_uuids = None
        else:
            # Use the shared logic from the SystemList
            self._active_uuids = self._system_list.get_connected_group(target_uuids)
            
            logging.info(f"Fluid Group: {len(target_uuids)} targets -> {len(self._active_uuids)} active components.")
        
        self._refresh_state()

    def _refresh_state(self):
        self._initial_p_vector = self._build_initial_vector()
        self._param_map = self._build_constraint_map()
        self._p_free_vector = self._initial_p_vector[self._param_map['is_free']]

    def _build_initial_vector(self) -> np.ndarray:
        initial_values = []
        for group in self._component_groups:
            for component in group:
                for param_name in self.PARAMETER_ORDER:
                    value = getattr(component, param_name)
                    initial_values.append(value if value is not None else 0.0) 
        return np.array(initial_values)
    
    # [NEW HELPER] Add this method to the class
    def _sanitize_data_structures(self):
        """
        Robustly fixes JSON serialization artifacts in the data core.
        1. Converts stringified ID keys "123" -> 123 in the ID map.
        2. Converts stringified Tuple keys "('uuid', 'param')" -> ('uuid', 'param') in constraints.
        """
        # A. Fix ID Map
        id_map = self._system_list._data.v1_id_to_uuid_map
        fixed_map = {}
        for k, v in id_map.items():
            if isinstance(k, str) and k.isdigit():
                fixed_map[int(k)] = v
            else:
                fixed_map[k] = v
        
        # B. Fix Parsed Constraints
        constraints = self._system_list._data.parsed_constraints
        fixed_constraints = {}
        for k, v in constraints.items():
            final_key = k
            
            # Handle Stringified Tuples: "('uuid', 'param')" or "('123', 'param')"
            if isinstance(k, str) and k.startswith('(') and k.endswith(')'):
                try:
                    # Strip parens and split
                    content = k[1:-1]
                    # Simple CSV parse respecting quotes
                    parts = [p.strip().strip("'").strip('"') for p in content.split(',')]
                    if len(parts) == 2:
                        p0, p1 = parts
                        # Try to restore Integer ID
                        if p0.isdigit(): p0 = int(p0)
                        final_key = (p0, p1)
                except:
                    pass # Keep original if parse fails

            # Handle Stringified Integers in Tuple: (123, 'z') vs ("123", 'z')
            if isinstance(final_key, tuple):
                 p0, p1 = final_key
                 if isinstance(p0, str) and p0.isdigit():
                     final_key = (int(p0), p1)

            fixed_constraints[final_key] = v

        # Apply fixes to the data object (in memory only, affects this session)
        # Note: We can't modify frozen dataclass, but we can update the local refs
        # used by this class. Ideally, we would update the system_list, but 
        # patching locally is sufficient for the fitter/view logic.
        return fixed_map, fixed_constraints

    def _build_constraint_map(self) -> Dict[str, np.ndarray]:
        # [FIX] Run sanitization first
        v1_id_map, v1_parsed = self._sanitize_data_structures()

        num_params = len(self._initial_p_vector)
        is_free = np.ones(num_params, dtype=bool)
        
        # Helper: Local import to ensure we can check types
        from ..structures import ParameterConstraintV2
        import dataclasses

        for i, (comp, param_name) in enumerate(self._get_component_param_iterator()):
            # Keys: Check both V1 (ID) and V2 (UUID) formats
            v1_key = (comp.id, param_name)
            v2_key = (comp.uuid, param_name)
            v2_uuid = comp.uuid
            
            # FLUID GROUP LOGIC
            force_freeze = False
            if self._active_uuids is not None:
                if v2_uuid not in self._active_uuids:
                    force_freeze = True

            # LOOKUP: Prioritize V2 key, fallback to V1 key
            raw_state = None
            if v2_key in v1_parsed:
                raw_state = v1_parsed[v2_key]
            elif v1_key in v1_parsed:
                raw_state = v1_parsed[v1_key]

            v2_constraint_obj: ParameterConstraintV2
            
            if raw_state:
                if isinstance(raw_state, dict):
                    c_free = raw_state.get('is_free', True)
                    v2_constraint_obj = self._convert_v1_to_v2(raw_state, v1_id_map, c_free)
                    c_expr = v2_constraint_obj.expression 
                else:
                    v2_constraint_obj = raw_state
                    c_free = v2_constraint_obj.is_free
                    c_expr = v2_constraint_obj.expression
                
                effective_is_free = c_free and (not force_freeze)
                if c_expr is not None:
                    effective_is_free = False 
                
                is_free[i] = effective_is_free
                
                if v2_constraint_obj.is_free != effective_is_free:
                     v2_constraint_obj = dataclasses.replace(v2_constraint_obj, is_free=effective_is_free)
            else:
                effective_is_free = True and (not force_freeze)
                is_free[i] = effective_is_free
                v2_constraint_obj = ParameterConstraintV2(is_free=effective_is_free)

            if v2_uuid not in self.v2_constraints_by_uuid:
                self.v2_constraints_by_uuid[v2_uuid] = {}
            self.v2_constraints_by_uuid[v2_uuid][param_name] = v2_constraint_obj

        return {'is_free': is_free}

    def _convert_v1_to_v2(self, v1_state, v1_id_map, effective_is_free: bool):
        expr = v1_state['expression']
        tgt_uuid = None
        if expr:
            match = self.V1_ID_REGEX.search(expr)
            if match:
                vid = int(match.group(1))
                tgt_uuid = v1_id_map.get(vid)
                if tgt_uuid:
                    expr = expr.replace(f"p[{vid}]", f"p['{tgt_uuid}']")
        return ParameterConstraintV2(effective_is_free, tgt_uuid, expr, None)

    def _get_component_param_iterator(self):
        for group in self._component_groups:
            for component in group:
                for param_name in self.PARAMETER_ORDER:
                    yield component, param_name

    @property
    def p_free_vector(self) -> np.ndarray:
        return self._p_free_vector
        
    def map_p_free_to_full(self, p_free_fitted: np.ndarray) -> np.ndarray:
        p_full = self._initial_p_vector.copy()
        p_full[self._param_map['is_free']] = p_free_fitted
        return p_full
    
    def get_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        lower = []
        upper = []
        for group in self._component_groups:
            for component in group:
                for param_name in self.PARAMETER_ORDER:
                    uuid = component.uuid
                    constraint = self.v2_constraints_by_uuid[uuid][param_name]
                    if constraint.is_free:
                        def_min, def_max = self.DEFAULT_BOUNDS[param_name]
                        lower.append(def_min)
                        upper.append(def_max)
        return np.array(lower), np.array(upper)

    @property
    def v2_constraints_for_save(self) -> Dict[str, Dict[str, ParameterConstraintV2]]:
        return self.v2_constraints_by_uuid