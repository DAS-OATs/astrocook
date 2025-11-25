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
    Supports "Selective Fitting" (Fluid Grouping).
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

        # If this set is populated, ONLY these UUIDs are free. Others are fixed.
        self._active_uuids: Optional[Set[str]] = None

        # 1. Group components
        self._component_groups = [list(self._system_list.components)]
        
        # 2. Build vectors
        self._refresh_state()
        
        logging.info(f"Initialized VoigtModelConstraintV2 with {len(self.p_free_vector)} free parameters.")

    def set_active_components(self, uuids: List[str] = None):
        """
        Defines the "Fluid Group". 
        If uuids is provided, only these components will be FREE.
        All others will be FIXED (but still contribute to the model).
        """
        if uuids is None:
            self._active_uuids = None
        else:
            self._active_uuids = set(uuids)
        
        self._refresh_state()

    def _refresh_state(self):
        """Rebuilds internal vectors and maps."""
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
    
    def _build_constraint_map(self) -> Dict[str, np.ndarray]:
        """
        Builds the mapping.
        """
        v1_parsed = self._system_list._data.parsed_constraints
        v1_id_map = self._system_list._data.v1_id_to_uuid_map

        num_params = len(self._initial_p_vector)
        is_free = np.ones(num_params, dtype=bool)
        
        for i, (comp, param_name) in enumerate(self._get_component_param_iterator()):
            
            v1_key = (comp.id, param_name)
            v2_uuid = comp.uuid
            
            # --- FLUID GROUP LOGIC ---
            force_freeze = False
            if self._active_uuids is not None:
                if v2_uuid not in self._active_uuids:
                    force_freeze = True

            v2_constraint_obj: ParameterConstraintV2
            
            if v1_key in v1_parsed:
                v1_state = v1_parsed[v1_key]
                
                # Calculate effective status
                effective_is_free = v1_state['is_free'] and (not force_freeze)
                if v1_state['expression'] is not None:
                    effective_is_free = False 
                
                is_free[i] = effective_is_free
                
                # FIX: Pass effective_is_free to helper
                v2_constraint_obj = self._convert_v1_to_v2(v1_state, v1_id_map, effective_is_free)
            else:
                is_free[i] = True and (not force_freeze)
                v2_constraint_obj = ParameterConstraintV2(is_free=is_free[i])

            if v2_uuid not in self.v2_constraints_by_uuid:
                self.v2_constraints_by_uuid[v2_uuid] = {}
            self.v2_constraints_by_uuid[v2_uuid][param_name] = v2_constraint_obj

        return {'is_free': is_free}

    def _convert_v1_to_v2(self, v1_state, v1_id_map, effective_is_free: bool):
        """Helper to convert V1 constraint dict to V2 object."""
        expr = v1_state['expression']
        tgt_uuid = None
        if expr:
            match = self.V1_ID_REGEX.search(expr)
            if match:
                vid = int(match.group(1))
                tgt_uuid = v1_id_map.get(vid)
                if tgt_uuid:
                    expr = expr.replace(f"p[{vid}]", f"p['{tgt_uuid}']")
        
        # FIX: Use the passed effective_is_free, NOT v1_state['is_free']
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
        """
        Returns the (lower_bounds, upper_bounds) vectors for the FREE parameters.
        """
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