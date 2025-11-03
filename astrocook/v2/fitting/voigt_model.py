from astropy.table import Table
import logging
import numpy as np
import re
from typing import List, Dict, Any, TYPE_CHECKING

# We assume SystemListV2 and ComponentDataV2 are importable from the V2 structure
from ..structures import ComponentDataV2, ParameterConstraintV2 
if TYPE_CHECKING:
    from ..system_list import SystemListV2
class VoigtModelConstraintV2:
    """
    Manages parameter constraints and converts component data into a flat P_free 
    vector required by SciPy optimization routines.
    (Replaces lmfit's Parameters object functionality)
    """
    
    # Define the order of parameters in the P_free vector (V2 Standard)
    PARAMETER_ORDER = ['z', 'logN', 'b', 'btur'] 

    # V1 expression regex (e.g., p[12], p[4], etc.)
    V1_ID_REGEX = re.compile(r"p\[(\d+)\]")

    def __init__(self, system_list: 'SystemListV2', constraints: Dict[str, Any] = None):
        """
        Initializes the constraint model using data from SystemListV2.
        """
        self._system_list = system_list
        self._param_names = ['z', 'logN', 'b', 'btur']

        # This will hold the new, V2-native (UUID-keyed) constraint map
        self.v2_constraints_by_uuid: Dict[str, Dict[str, ParameterConstraintV2]] = {}          

        # 1. NEW: Group components
        self._component_groups = self._build_groups() 
        
        # 2. Build the initial parameter vector (P_full) based on groups
        self._initial_p_vector = self._build_initial_vector()
        
        # 3. Build the constraint map
        if self._system_list:
            self._build_constraint_map()
        
        #logging.info(f"Initialized VoigtModelConstraintV2 with {len(self.p_free_vector)} free parameters.")

    def _build_groups(self, overlap_threshold=0.005) -> List[List[ComponentDataV2]]:
        """
        Groups components based on V1 criteria (constraints or overlapping profiles).
        Returns a list of lists, where each inner list is a group of components to be fitted together.
        
        Note: This is simplified and currently assumes all components form a single group 
        until the full profile calculation logic is integrated.
        """
        
        # --- PHASE 1: Simple Grouping (All components in one group) ---
        # This will be replaced by the complex profile overlap check later.
        if not self._system_list.components:
            return []
            
        return [list(self._system_list.components)]
        # ---------------------------------------------------------------
        
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
        Translates a single V1 constraint dict (from parsed_constraints) into
        a V2-native ParameterConstraintV2 object, converting ID-based
        expressions to UUID-based expressions.
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
    
    def _build_constraint_map(self):
        """
        Initializes or loads the constraint map.
        - If V2 data is present, load it.
        - If V1 data is present, migrate it.
        - If neither, create a default (all free) map.
        """
        data = self._system_list._data # Get the immutable data core
        
        # 1. Initialize all parameters as free by default
        self._initialize_all_params()
        
        # 2. Check for and load V2 or V1 constraints
        if data.v2_constraints_map:
            logging.info(f"Loading {len(data.v2_constraints_map)} V2 constraint sets from metadata...")
            loaded_map = {}
            for comp_uuid, params_dict in data.v2_constraints_map.items():
                loaded_map[comp_uuid] = {}
                for param_name, constr_dict in params_dict.items():
                    try:
                        # Create the dataclass instance from the dictionary
                        loaded_map[comp_uuid][param_name] = ParameterConstraintV2(**constr_dict)
                    except TypeError as e:
                        logging.warning(f"Failed to load constraint for {comp_uuid}.{param_name}: {e}. Setting as 'free'.")
                        loaded_map[comp_uuid][param_name] = ParameterConstraintV2(is_free=True)
            
            # Now assign the map of *ParameterConstraintV2 objects*
            self.v2_constraints_by_uuid = loaded_map

        elif data.parsed_constraints:
            # This is a V1-migrated session.
            # Call the *existing* V1 migration helper.
            logging.info(f"Migrating {len(data.parsed_constraints)} V1 constraints...")
            self.load_v1_constraints_as_v2(data.parsed_constraints, data.v1_id_to_uuid_map)
            
        # 3. Count free parameters
        free_param_count = self._count_free_params()
        logging.info(f"Initialized VoigtModelConstraintV2 with {free_param_count} free parameters.")

    def _initialize_all_params(self):
        """Sets all parameters for all components to 'free' by default."""
        self.v2_constraints_by_uuid = {}
        default_constraint = ParameterConstraintV2(is_free=True)
        
        for comp in self._system_list.components:
            uuid = comp.uuid
            self.v2_constraints_by_uuid[uuid] = {}
            for param in self._param_names:
                self.v2_constraints_by_uuid[uuid][param] = default_constraint

    def _count_free_params(self) -> int:
        """Counts the number of parameters marked as 'is_free'."""
        count = 0
        for comp_uuid, params in self.v2_constraints_by_uuid.items():
            for param_name, constraint in params.items():
                if constraint.is_free:
                    count += 1
        return count

    def load_v1_constraints_as_v2(self, 
                                  v1_parsed_constraints: Dict[tuple, Dict], 
                                  v1_id_to_uuid_map: Dict[int, str]):
        """
        Migrates V1-style constraints (from .acs header) into the
        V2-native (UUID-based) constraint map.
        """
        self.v1_parsed_constraints = v1_parsed_constraints
        self.v1_id_to_uuid_map = v1_id_to_uuid_map

        for (v1_id, param_name), constr_data in v1_parsed_constraints.items():
            try:
                # Find the UUID for this V1 ID
                comp_uuid = self.v1_id_to_uuid_map.get(v1_id)
                if not comp_uuid:
                    logging.warning(f"V1 constraint for ID {v1_id} skipped: No matching UUID.")
                    continue
                
                # Find the UUID for the *target* V1 ID
                target_v1_id = constr_data.get('target_id') # V1 used 'target_id'
                target_uuid = None
                if target_v1_id is not None:
                    target_uuid = self.v1_id_to_uuid_map.get(target_v1_id)
                    if not target_uuid:
                        logging.warning(f"V1 constraint target ID {target_v1_id} skipped: No matching UUID.")
                
                # Create the V2 ParameterConstraintV2 object
                v2_constraint = ParameterConstraintV2(
                    is_free=bool(constr_data.get('is_free', True)),
                    target_uuid=target_uuid,
                    expression=constr_data.get('expr'),
                    v1_target_id=target_v1_id
                )
                
                # Overwrite the default "free" constraint in our map
                if comp_uuid in self.v2_constraints_by_uuid:
                    self.v2_constraints_by_uuid[comp_uuid][param_name] = v2_constraint
                else:
                    logging.warning(f"V1 constraint for UUID {comp_uuid} skipped: UUID not in component list.")

            except Exception as e:
                logging.error(f"Error migrating V1 constraint {(v1_id, param_name)}: {e}")

    def _get_component_param_iterator(self):
        """Helper to yield (ComponentDataV2, param_name) tuple for P_full mapping."""
        for component in self._system_list.components:
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
        # This step would apply linked values based on link_target_index.
        # For our default state (no links), this step is a simple identity function.
        
        return p_full
    
    @property
    def v2_constraints_for_save(self) -> Dict[str, Dict[str, ParameterConstraintV2]]:
        """
        Exposes the V2-native, UUID-keyed constraint map for serialization.
        """
        return self.v2_constraints_by_uuid