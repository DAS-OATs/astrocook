from astropy.table import Table
import logging
import numpy as np
import re
from typing import List, Dict, Any, TYPE_CHECKING, Optional, Set, Tuple

from astrocook.core.structures import ComponentDataV2, ParameterConstraintV2 

if TYPE_CHECKING:
    from astrocook.core.system_list import SystemListV2

class VoigtModelConstraintV2:
    """
    Manages parameter constraints and vector mapping for Voigt fitting.

    This class acts as a bridge between the high-level, object-oriented representation
    of absorption systems (:class:`~astrocook.core.system_list.SystemListV2`) and the
    flat parameter arrays required by numerical optimizers.

    It handles:
    
    1.  **Selective Fitting**: Identifying which components are "active" (free to vary)
        and which are "frozen" based on user selection and grouping logic.
    2.  **Parameter Linking**: Evaluating mathematical expressions to link parameters
        (e.g., tying the redshift of two components together).
    3.  **Vectorization**: Converting component objects to a flat vector ``p`` and back.

    Parameters
    ----------
    system_list : SystemListV2
        The system list containing the components to model.
    constraints : dict, optional
        Initial constraints dictionary (legacy argument, usually handled internally).
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
        
        # Linking State
        self._compiled_links: List[Tuple[int, Any]] = [] # [(target_index, code_object)]
        self._uuid_to_idx: Dict[str, int] = {}
        
        self._refresh_state()
        
        logging.debug(f"Initialized VoigtModelConstraintV2 with {len(self.p_free_vector)} free parameters.")

    def set_active_components(self, target_uuids: List[str] = None, group_depth: int = 2):
        """
        Define the set of components to be fit (the "Fluid Group").

        When fitting a specific line, the fitter must also include any lines that are:
        1.  Physically linked (via parameter constraints).
        2.  Spectrally overlapping (blended).

        This method calculates that connected group and marks all other components
        as "frozen" (fixed parameters) for the upcoming fit.

        Parameters
        ----------
        target_uuids : list of str, optional
            The UUIDs of the components the user explicitly selected.
            If ``None``, ALL components are considered active.
        group_depth : int, optional
            How many degrees of separation to traverse when finding connected components.
            - ``1``: Direct neighbors/links only.
            - ``2``: Neighbors of neighbors (Friends of Friends).
            Defaults to ``2``.
        """
        if target_uuids is None:
            self._active_uuids = None
        else:
            # Use the shared logic from the SystemList
            self._active_uuids = self._system_list.get_connected_group(target_uuids)
            logging.debug(f"Fluid Group: {len(target_uuids)} targets -> {len(self._active_uuids)} active components.")
        
        self._refresh_state()

    def _refresh_state(self):
        """Rebuilds internal vectors and maps based on current constraints."""
        self._initial_p_vector = self._build_initial_vector()
        self._param_map = self._build_constraint_map()
        self._p_free_vector = self._initial_p_vector[self._param_map['is_free']]
        self._compile_links()

    def _build_initial_vector(self) -> np.ndarray:
        """
        Flattens all component parameters into a single 1D array.
        Also rebuilds the UUID -> Index map.
        """
        initial_values = []
        # Also build the UUID -> Index map here for linking
        self._uuid_to_idx = {}
        idx_counter = 0
        
        for group in self._component_groups:
            for component in group:
                self._uuid_to_idx[component.uuid] = idx_counter
                for param_name in self.PARAMETER_ORDER:
                    value = getattr(component, param_name)
                    initial_values.append(value if value is not None else 0.0) 
                    idx_counter += 1
        return np.array(initial_values)
    
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
            
            # Handle Stringified Tuples
            if isinstance(k, str) and k.startswith('(') and k.endswith(')'):
                try:
                    content = k[1:-1]
                    parts = [p.strip().strip("'").strip('"') for p in content.split(',')]
                    if len(parts) == 2:
                        p0, p1 = parts
                        if p0.isdigit(): p0 = int(p0)
                        final_key = (p0, p1)
                except: pass

            # Handle Stringified Integers in Tuple
            if isinstance(final_key, tuple):
                 p0, p1 = final_key
                 if isinstance(p0, str) and p0.isdigit():
                     final_key = (int(p0), p1)

            fixed_constraints[final_key] = v

        return fixed_map, fixed_constraints

    def _build_constraint_map(self) -> Dict[str, np.ndarray]:
        """
        Constructs the boolean mask identifying free parameters.
        
        Logic:
        1. Check explicit constraints (User said "Freeze").
        2. Check Fluid Group (If not active, force Freeze).
        3. Check Linking (If linked via expression, force Freeze).
        """
        v1_id_map, v1_parsed = self._sanitize_data_structures()

        num_params = len(self._initial_p_vector)
        is_free = np.ones(num_params, dtype=bool)
        
        from astrocook.core.structures import ParameterConstraintV2
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
                    effective_is_free = False # Linked parameters are not free
                
                is_free[i] = effective_is_free
                
                if v2_constraint_obj.is_free != effective_is_free:
                     v2_constraint_obj = dataclasses.replace(v2_constraint_obj, is_free=effective_is_free)
            else:
                effective_is_free = True and (not force_freeze)
                if param_name == 'btur':
                    effective_is_free = False
                is_free[i] = effective_is_free
                v2_constraint_obj = ParameterConstraintV2(is_free=effective_is_free)

            if v2_uuid not in self.v2_constraints_by_uuid:
                self.v2_constraints_by_uuid[v2_uuid] = {}
            self.v2_constraints_by_uuid[v2_uuid][param_name] = v2_constraint_obj

        return {'is_free': is_free}
    
    def _compile_links(self):
        """Pre-compiles expressions for linked parameters to speed up fitting."""
        self._compiled_links = []
        
        for uuid, p_map in self.v2_constraints_by_uuid.items():
            if uuid not in self._uuid_to_idx: continue
            comp_idx = self._uuid_to_idx[uuid]
            
            for param, constraint in p_map.items():
                if not constraint.is_free and constraint.expression:
                    # Determine target index in p_full
                    try:
                        param_offset = self.PARAMETER_ORDER.index(param)
                        target_idx = comp_idx + param_offset
                        
                        # Compile expression
                        # We use 'p' as the dictionary of components, e.g. p['uuid'].b
                        code = compile(constraint.expression, f"<link {uuid}.{param}>", 'eval')
                        self._compiled_links.append((target_idx, code))
                    except Exception as e:
                        logging.error(f"Failed to compile link for {uuid}.{param}: {e}")

    def _convert_v1_to_v2(self, v1_state, v1_id_map, effective_is_free: bool):
        expr = v1_state['expression']
        tgt_uuid = None
        if expr:
            # Convert "p[123]" -> "p['uuid']"
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
        """
        The 1D array of currently free (varying) parameters.

        This vector serves as the input ``x`` for the optimizer. It excludes
        fixed parameters and dependent (linked) parameters.
        """
        return self._p_free_vector
        
    def map_p_free_to_full(self, p_free_fitted: np.ndarray) -> np.ndarray:
        """
        Reconstruct the full parameter vector from the free parameters.

        This method:
        1. Places the free parameter values into their correct positions in the full vector.
        2. Leaves fixed parameters at their initial values.
        3. Evaluates and applies all linked parameter expressions (e.g., ``z2 = p['uuid1'].z``).

        Parameters
        ----------
        p_free_fitted : np.ndarray
            The vector of optimized free parameters (output from optimizer).

        Returns
        -------
        np.ndarray
            The full vector containing all parameters (z, logN, b, btur) for all components.
        """
        # 1. Fill Free Parameters
        p_full = self._initial_p_vector.copy()
        p_full[self._param_map['is_free']] = p_free_fitted
        
        # 2. Evaluate Links
        if self._compiled_links:
            # Create a lightweight context for 'p["uuid"].param' access
            # We map p['uuid'] to a proxy object that reads from p_full
            
            class ComponentProxy:
                __slots__ = ('_arr', '_idx')
                def __init__(self, arr, idx):
                    self._arr = arr
                    self._idx = idx
                @property
                def z(self): return self._arr[self._idx]
                @property
                def logN(self): return self._arr[self._idx+1]
                @property
                def b(self): return self._arr[self._idx+2]
                @property
                def btur(self): return self._arr[self._idx+3]

            # Create the 'p' dictionary
            # Note: uuid_to_idx is computed in _refresh_state
            p_dict = {
                u: ComponentProxy(p_full, i) 
                for u, i in self._uuid_to_idx.items()
            }
            
            # Common math functions for expressions
            eval_globals = {'p': p_dict, 'np': np, 'log': np.log10, 'ln': np.log}
            
            # Execute compiled expressions
            # Note: If chains exist (A->B->C), we might need multiple passes or topological sort.
            # For now, we assume standard depth-1 links or ordered calculation.
            for target_idx, code in self._compiled_links:
                try:
                    val = eval(code, eval_globals)
                    p_full[target_idx] = float(val)
                except Exception:
                    pass # Keep initial value if evaluation fails

        return p_full
    
    def get_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Return lower and upper bounds for the free parameters.

        Extracts the physical limits (e.g., b > 0) for every parameter currently
        in the free vector.

        Returns
        -------
        tuple
            A tuple ``(lower_bounds, upper_bounds)``, each being a numpy array
            matching the size of ``p_free_vector``.
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
    
    def get_active_param_labels(self) -> List[str]:
        """
        Return a list of formatting labels for all active free parameters.

        Iterates through the components and parameter constraints to build
        identifying strings (e.g., ``'z (1)'``, ``'logN (2)'``) for every
        parameter that is currently free to vary during the fit.
        Used primarily for visualization (e.g., Corner Plots).

        Returns
        -------
        list of str
            A list of parameter labels, matching the order of the ``p_free_vector``.
        """
        labels = []
        for i, (comp, param_name) in enumerate(self._get_component_param_iterator()):
            if self._param_map['is_free'][i]:
                labels.append(f"{param_name} ({comp.id})")
        return labels