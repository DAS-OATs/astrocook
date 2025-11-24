from astropy.table import Table
import logging
import numpy as np
import re
from typing import List, Dict, Any, TYPE_CHECKING, Tuple

# We assume SystemListV2 and ComponentDataV2 are importable from the V2 structure
from ..structures import ComponentDataV2, ParameterConstraintV2 
if TYPE_CHECKING:
    from ..system_list import SystemListV2

class VoigtModelConstraintV2:
    """
    Manages parameter constraints and converts component data into a flat P_free 
    vector required by SciPy optimization routines.
    """
    
    # Define the order of parameters in the P_free vector (V2 Standard)
    PARAMETER_ORDER = ['z', 'logN', 'b', 'btur'] 

    # Default Physical Limits (Min, Max)
    DEFAULT_BOUNDS = {
        'z':    (-0.5, 10.0),    # Redshift usually positive, but allow small neg for calibration
        'logN': (10.0, 23.0),    # Column density limits
        'b':    (1.0, 300.0),    # Doppler param (km/s)
        'btur': (0.0, 100.0)     # Turbulence
    }

    # V1 expression regex (e.g., p[12], p[4], etc.)
    V1_ID_REGEX = re.compile(r"p\[(\d+)\]")

    def __init__(self, system_list: 'SystemListV2', constraints: Dict[str, Any] = None):
        """
        Initializes the constraint model using data from SystemListV2.
        """
        self._system_list = system_list

        # This will hold the new, V2-native (UUID-keyed) constraint map
        self.v2_constraints_by_uuid: Dict[str, Dict[str, ParameterConstraintV2]] = {}          

        # If this set is populated, ONLY these UUIDs are free. Others are fixed.
        self._active_uuids: Optional[Set[str]] = None

        # 1. NEW: Group components
        self._component_groups = self._build_groups() 
        
        # 2. Build the initial parameter vector (P_full) based on groups
        self._initial_p_vector = self._build_initial_vector()
        
        # 3. Build the constraint map (CRITICAL: Must set self._param_map here)
        self._param_map = self._build_constraint_map()

        # 4. Cache the P_free vector (the result used by SciPy)
        self._p_free_vector = self._initial_p_vector[self._param_map['is_free']]
        
        logging.info(f"Initialized VoigtModelConstraintV2 with {len(self.p_free_vector)} free parameters.")

    def _build_groups(self, overlap_threshold=0.005) -> List[List[ComponentDataV2]]:
        """
        Groups components based on V1 criteria (constraints or overlapping profiles).
        Returns a list of lists, where each inner list is a group of components to be fitted together.
        """
        # --- PHASE 1: Simple Grouping (All components in one group) ---
        if not self._system_list.components:
            return []
            
        return [list(self._system_list.components)]
        
    def _build_initial_vector(self) -> np.ndarray:
        """
        Assembles a flat NumPy array based on the components WITHIN the groups.
        """
        initial_values = []
        
        for group in self._component_groups:
            for component in group:
                for param_name in self.PARAMETER_ORDER:
                    value = getattr(component, param_name)
                    initial_values.append(value if value is not None else 0.0) 
                    
        return np.array(initial_values)

    def _convert_v1_state_to_v2_constraint(
        self, 
        v1_state: Dict[str, Any], 
        v1_id_to_uuid_map: Dict[int, str]
    ) -> ParameterConstraintV2:
        """
        Translates a single V1 constraint dict into a V2-native ParameterConstraintV2 object.
        """
        v1_expression = v1_state['expression']
        target_uuid: str = None
        v1_target_id: int = None
        v2_expression: str = v1_expression

        if v1_expression:
            match = self.V1_ID_REGEX.search(v1_expression)
            if match:
                try:
                    # 1. Found a V1 ID link (e.g., p[12])
                    v1_target_id = int(match.group(1))
                    
                    # 2. Find the corresponding UUID
                    target_uuid = v1_id_to_uuid_map.get(v1_target_id)
                    
                    if target_uuid:
                        # 3. Create the new V2 expression string
                        v2_expression = v1_expression.replace(
                            f"p[{v1_target_id}]", f"p['{target_uuid}']"
                        )
                    else:
                        logging.warning(f"Constraint link error: V1 ID {v1_target_id} not found in UUID map.")
                        
                except Exception as e:
                    logging.error(f"Failed to parse V1 expression '{v1_expression}': {e}")
        
        return ParameterConstraintV2(
            is_free=v1_state['is_free'],
            target_uuid=target_uuid,
            expression=v2_expression,
            v1_target_id=v1_target_id
        )
    
    def _build_constraint_map(self) -> Dict[str, np.ndarray]:
        """
        Builds the mapping, incorporating loaded V1 constraints.
        Returns a dictionary containing boolean masks and index maps.
        """
        # Get the V1-style constraint map (keyed by V1 ID)
        v1_parsed_constraints = self._system_list._data.parsed_constraints
        
        # Get the V1->V2 ID-UUID map
        v1_id_to_uuid_map = self._system_list._data.v1_id_to_uuid_map

        num_params = len(self._initial_p_vector)
        is_free = np.ones(num_params, dtype=bool)
        link_target_index = np.arange(num_params) # Placeholder for future logic
        
        # Iterate over all parameters in the vector order
        for i, (comp, param_name) in enumerate(self._get_component_param_iterator()):
            
            # 1. Use the V1 ID to look up the constraint
            v1_key = (comp.id, param_name)
            v2_uuid = comp.uuid
            
            # Fluid Group Logic:
            # If we have an active set, and this component is NOT in it, FORCE FREEZE.
            force_freeze = False
            if self._active_uuids is not None:
                if v2_uuid not in self._active_uuids:
                    force_freeze = True

            v2_constraint_obj: ParameterConstraintV2
            
            if v1_key in v1_parsed_constraints:
                # 2a. A constraint EXISTS for this parameter
                v1_state = v1_parsed_constraints[v1_key]
                
                # 3a. Update the SciPy 'is_free' array
                is_free[i] = v1_state['is_free'] and (not force_freeze)
                if v1_state['expression'] is not None:
                    # Linked/fixed parameters are not "free" in the SciPy vector
                    is_free[i] = False 
                
                # 4a. Convert this V1 state to a V2 object
                v2_constraint_obj = self._convert_v1_state_to_v2_constraint(
                    v1_state, v1_id_to_uuid_map
                )
                
            else:
                # 2b. No constraint found. This is a default FREE parameter.
                is_free[i] = True and (not force_freeze)
                v2_constraint_obj = ParameterConstraintV2(is_free=True)

            # 5. Store the V2 object in the new UUID-keyed map
            if v2_uuid not in self.v2_constraints_by_uuid:
                self.v2_constraints_by_uuid[v2_uuid] = {}
                
            self.v2_constraints_by_uuid[v2_uuid][param_name] = v2_constraint_obj

        return {
            'is_free': is_free,
            'link_target_index': link_target_index
        }

    def _convert_v1_to_v2(self, v1_state, v1_id_map):
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
        return ParameterConstraintV2(v1_state['is_free'], tgt_uuid, expr, None)

    def _get_component_param_iterator(self):
        """Helper to yield (ComponentDataV2, param_name) tuple for P_full mapping."""
        for group in self._component_groups:
            for component in group:
                for param_name in self.PARAMETER_ORDER:
                    yield component, param_name

    @property
    def p_free_vector(self) -> np.ndarray:
        """
        Returns the final flat vector containing only the FREE parameters for SciPy.
        """
        return self._p_free_vector
        
    def map_p_free_to_full(self, p_free_fitted: np.ndarray) -> np.ndarray:
        """
        Maps the fitted P_free vector back to the original P_full vector size, 
        filling in fixed and linked parameters.
        """
        p_full = self._initial_p_vector.copy()
        
        # 1. Inject the fitted values back into the free slots
        p_full[self._param_map['is_free']] = p_free_fitted
        
        # 2. Apply linking logic (simplified, assuming direct indexing for now)
        # TODO: Implement expression evaluation for linked parameters here
        
        return p_full
    
    @property
    def v2_constraints_for_save(self) -> Dict[str, Dict[str, ParameterConstraintV2]]:
        """
        Exposes the V2-native, UUID-keyed constraint map for serialization.
        """
        return self.v2_constraints_by_uuid
    
    def set_active_components(self, uuids: List[str] = None):
        """
        Defines the "Fluid Group". 
        If uuids is provided, only these components will be FREE.
        All others will be FIXED (but still contribute to the model).
        Pass None to free all components.
        """
        if uuids is None:
            self._active_uuids = None
        else:
            self._active_uuids = set(uuids)
        
        # Rebuild the parameter map based on the new selection
        self._refresh_state()

    def _refresh_state(self):
        """Rebuilds internal vectors and maps."""
        self._initial_p_vector = self._build_initial_vector()
        self._param_map = self._build_constraint_map()
        self._p_free_vector = self._initial_p_vector[self._param_map['is_free']]

    def get_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Returns the (lower_bounds, upper_bounds) vectors for the FREE parameters.
        Used by scipy.optimize.least_squares.
        """
        lower = []
        upper = []
        
        # Iterate exactly as we do for building the p_vector
        for group in self._component_groups:
            for component in group:
                for param_name in self.PARAMETER_ORDER:
                    
                    # 1. Check if this parameter is FREE
                    # We need to check the map we built in __init__
                    # Efficient way: We iterate in the same order, so we can check the constraint obj
                    # in self.v2_constraints_by_uuid
                    
                    uuid = component.uuid
                    constraint = self.v2_constraints_by_uuid[uuid][param_name]
                    
                    if constraint.is_free:
                        # 2. Apply Bounds
                        # TODO: In the future, allow per-parameter custom bounds from V1 data
                        def_min, def_max = self.DEFAULT_BOUNDS[param_name]
                        lower.append(def_min)
                        upper.append(def_max)
        
        return np.array(lower), np.array(upper)