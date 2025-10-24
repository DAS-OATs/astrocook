from astropy.table import Table
from typing import List, Dict, Any, TYPE_CHECKING
import numpy as np
import logging

# We assume SystemListV2 and ComponentDataV2 are importable from the V2 structure
from ..structures import ComponentDataV2 # Now ComponentDataV2
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

    def __init__(self, system_list: 'SystemListV2', constraints: Dict[str, Any] = None):
        """
        Initializes the constraint model using data from SystemListV2.
        """
        self._system_list = system_list
        self._constraints = constraints if constraints is not None else {}
        
        # 1. NEW: Group components
        self._component_groups = self._build_groups() 
        
        # 2. Build the initial parameter vector (P_full) based on groups
        self._initial_p_vector = self._build_initial_vector()
        
        # 3. Build the constraint map
        self._param_map = self._build_constraint_map()

        # 4. Cache the P_free vector (the result used by SciPy)
        self._p_free_vector = self._initial_p_vector[self._param_map['is_free']]
        
        logging.info(f"Initialized VoigtModelConstraintV2 with {len(self.p_free_vector)} free parameters.")

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

    def _build_constraint_map(self) -> Dict[str, np.ndarray]:
        """
        Builds the mapping, incorporating loaded V1 constraints.
        """
        v2_constraints = self._system_list.parsed_constraints
        
        num_params = len(self._initial_p_vector)
        is_free = np.ones(num_params, dtype=bool)
        link_target_index = np.arange(num_params)
        
        # --- PHASE 2: CRITICAL LOGIC: Translate V1 constraints to V2 arrays ---
        
        # We need a reverse map from (Component ID, Parameter Name) to the P_full index.
        
        for i, (comp, param_name) in enumerate(self._get_component_param_iterator()):
            key = (comp.id, param_name)
            
            if key in v2_constraints:
                state = v2_constraints[key]
                
                # 1. Handle Freezing (is_free=False)
                if state['is_free'] is False:
                    is_free[i] = False
                    
                # 2. Handle Linking (expression != None)
                if state['expression'] is not None:
                    # Linking is a complex calculation; for P_free filtering, we flag it as not free.
                    is_free[i] = False
        
        return {
            'is_free': is_free,
            'link_target_index': link_target_index
        }

    def _load_v1_constraints(self, v1_models_t: Any) -> Dict[str, Any]:
        """
        Parses the V1 SystList._mods_t (lmfit models/IDs) to identify 
        initial constraints (freezes, links).
        """
        constraints = {}
        if not isinstance(v1_models_t, Table):
            return constraints
            
        # The V1 structure stores the lmfit model in the 'mod' column.
        for row in v1_models_t:
            model_obj = row['mod'] # The lmfit.Model object
            
            if not hasattr(model_obj, 'params'):
                continue

            params = model_obj.params # The lmfit.Parameters object
            
            # --- Iterate over lmfit Parameters to extract state ---
            for name, param in params.items():
                
                # Check for parameter components (e.g., lines_voigt_123_z)
                if any(param_name in name for param_name in self.PARAMETER_ORDER):
                    
                    # 1. Determine Component ID and Parameter Name (e.g., 'z')
                    # V1 lmfit prefix is typically lines_voigt_<ID>_<param>
                    parts = name.split('_')
                    
                    # Assuming the structure is PREF_ID_PARAM:
                    param_name = parts[-1] 
                    comp_id = int(parts[-2])
                    
                    # 2. Store Constraint State
                    
                    # Key for the constraints dictionary: (Component ID, Parameter Name)
                    key = (comp_id, param_name)
                    
                    # Value: True/False for vary, or the expression string for linking
                    constraints[key] = {
                        'is_free': param.vary,
                        'expression': param.expr if param.expr else None
                    }

        return constraints

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