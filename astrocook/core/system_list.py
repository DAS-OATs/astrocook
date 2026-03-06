from astropy.table import Table, Column
import astropy.units as au
import contextlib
import dataclasses
import logging
import numpy as np
import ast
from typing import Any, Dict, List, Optional, Set, Tuple

from astrocook.fitting.voigt_model import VoigtModelConstraintV2
from astrocook.core.structures import ComponentDataV2, SystemListDataV2, ParameterConstraintV2
from astrocook.core.system_list_migration import migrate_component_v2_to_v1
from astrocook.legacy.syst_list import SystList as SystListV1
from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS

def get_all_w_rest(series_name: str) -> List[float]:
    """
    Returns a list of rest wavelengths in Angstroms for a given series name.
    Handles single lines, standard multiplets, and custom comma-separated lists.
    """
    # 1. Check if it's a known single transition
    if series_name in ATOM_DATA:
        return [ATOM_DATA[series_name]['wave']]
    
    # 2. Check if it's a standard multiplet tag (e.g., 'CIV')
    if series_name in STANDARD_MULTIPLETS:
        waves = []
        for line_name in STANDARD_MULTIPLETS[series_name]:
            if line_name in ATOM_DATA:
                waves.append(ATOM_DATA[line_name]['wave'])
        return waves

    # 3. Handle custom comma-separated strings (e.g., 'CIV_1548,CIV_1550')
    if ',' in series_name:
        waves = []
        for part in series_name.split(','):
            name = part.strip()
            if name in ATOM_DATA:
                waves.append(ATOM_DATA[name]['wave'])
        return waves

    return []

def are_overlapping(c1: ComponentDataV2, c2: ComponentDataV2, constraints_map: dict) -> bool:
    """
    Checks if two components should be grouped together for fitting.
    """
    c_kms = 299792.458
    
    # 1. CONSTRAINT CHECK (Highest Priority)
    # If components are linked via expressions, they MUST be fitted together.
    for c_src, c_other in [(c1, c2), (c2, c1)]:
        if c_src.uuid in constraints_map:
            for constr in constraints_map[c_src.uuid].values():
                if hasattr(constr, 'target_uuid') and constr.target_uuid == c_other.uuid:
                    return True

    # 2. SPECTRAL OVERLAP (The "Blend" Logic)
    # Check if ANY transition from system 1 is near ANY transition from system 2.
    # We use a threshold based on the thermal width (b-parameter) or a 30 km/s floor.
    threshold_kms = max(4.0 * max(c1.b, c2.b), 30.0)
    
    obs1 = [w * (1 + c1.z) for w in get_all_w_rest(c1.series)]
    obs2 = [w * (1 + c2.z) for w in get_all_w_rest(c2.series)]
    
    for w1 in obs1:
        for w2 in obs2:
            dv = c_kms * abs(w1 - w2) / ((w1 + w2) / 2.0)
            if dv < threshold_kms:
                return True

    # 3. KINEMATIC PROXIMITY (Same Species)
    # If they are the same species and close in velocity, group them.
    if c1.series == c2.series:
        dz = abs(c1.z - c2.z)
        dv_kin = c_kms * dz / (1 + c1.z)
        if dv_kin < threshold_kms:
            return True
            
    return False

# --- 3. System List API Layer (The orchestrator) ---
class SystemListV2:
    """
    API layer for managing component lists, constraints, and fitting.

    This class wraps the immutable :class:`~astrocook.core.structures.SystemListDataV2`
    and provides methods to add, update, delete, and group components. It also
    initializes the constraint model used by the fitting engine.

    Parameters
    ----------
    data : SystemListDataV2
        The core data container holding the list of components and constraint maps.
    """
    
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
            logging.debug("SystemListV2: Hydrated constraints from clean JSON map.")
        
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
                    # If parsing 'btur', default to False. For everything else, True.
                    # We need to extract the param name from the key.
                    default_free = True
                    if isinstance(final_key, tuple) and len(final_key) == 2:
                        p_name = final_key[1]
                        if p_name == 'btur':
                            default_free = False

                    val_obj = ParameterConstraintV2(
                        is_free=v.get('is_free', default_free),
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
        """List of all component data objects."""
        return self._data.components

    @property
    def t(self) -> Table:
        """
        Returns an Astropy Table representation of the system list.
        Useful for quick inspection or V1 compatibility.
        """
        z_list, dz_list, logN_list, dlogN_list, b_list, db_list = ([] for _ in range(6))
        btur_list, dbtur_list, func_list, series_list, id_list = ([] for _ in range(5))
        chi2_list, resol_list = ([] for _ in range(2))
        
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
            chi2_list.append(c.chi2 if c.chi2 is not None else np.nan)
            resol_list.append(c.resol if c.resol is not None else np.nan)

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
        t['chi2'] = Column(chi2_list, unit=au.dimensionless_unscaled)
        t['resol'] = Column(resol_list, unit=au.dimensionless_unscaled)
        t['chi2'] = Column(chi2_list, unit=au.dimensionless_unscaled)
        t['resol'] = Column(resol_list, unit=au.dimensionless_unscaled)

        return t

    @property
    def z(self) -> List[float]:
        """List of redshifts for all components."""
        return [c.z for c in self._data.components]

    @property
    def series(self) -> List[str]:
        """List of series names for all components."""
        return [c.series for c in self._data.components]
        
    @property
    def id(self) -> List[int]:
        """List of integer IDs for all components."""
        return [c.id for c in self._data.components]

    @property
    def v1_models_t(self) -> Any:
        return self._data.v1_models_t

    @property
    def constraint_status(self) -> Dict[str, np.ndarray]:
        """Dictionary reflecting the free/fixed status of parameters."""
        return self.constraint_model._param_map
    
    @property
    def parsed_constraints(self) -> Dict[Tuple[int, str], Dict[str, Any]]:
        return self._data.parsed_constraints
    
    @property
    def v2_constraints_for_save(self) -> Dict[str, Dict[str, Any]]:
        if self.constraint_model:
            return self.constraint_model.v2_constraints_by_uuid
        return {}
    
    # --- Context Manager for Safe Fitting ---
    @contextlib.contextmanager
    def fitting_context(self, active_uuids: List[str] = None, group_depth: int = 2):
        """
        Context manager for safe fitting.

        Temporarily sets the active fitting mask on the constraint model and
        guarantees cleanup (resetting the mask) even if errors occur.

        Parameters
        ----------
        active_uuids : list of str, optional
            List of UUIDs to set as active.
        group_depth : int, optional
            Grouping depth for finding connected components. Defaults to ``2``.
        """
        if not self.constraint_model:
            yield
            return

        # Snapshot original state safely
        # Copy the set to avoid reference issues
        original_active = None
        if self.constraint_model._active_uuids is not None:
            original_active = set(self.constraint_model._active_uuids)

        try:
            self.constraint_model.set_active_components(active_uuids, group_depth=group_depth)
            yield
        finally:
            # Robust cleanup
            try:
                self.constraint_model.set_active_components(original_active)
            except Exception as e:
                logging.error(f"fitting_context cleanup failed: {e}. Forcing reset to None.")
                # Fail-safe: Unlock everything to prevent persistent 'Frozen' state
                try:
                    self.constraint_model.set_active_components(None)
                except Exception:
                    # Last resort: direct attribute clear
                    self.constraint_model._active_uuids = None

    def to_v1_systlist(self) -> Optional[Table]:
        """
        Convert the V2 system list to a V1-compatible Astropy Table.

        Used for saving archives in the legacy ``.acs`` format.

        Returns
        -------
        astropy.table.Table or None
            A table containing the system list in V1 format, or ``None`` if the list is empty.
        """
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
            'chi2': [c.chi2 if c.chi2 is not None else np.nan for c in self.components],
            'resol': [c.resol if c.resol is not None else np.nan for c in self.components]
        }
        return Table(data_dict)
    
    def add_components_from_regions(self, spec: 'SpectrumV2', region_id_col: str,
                                  series_map: Dict[int, str], z_map: Dict[int, float]) -> 'SystemListV2':
        """
        Populate the list with initial components based on identified regions.

        Parameters
        ----------
        spec : SpectrumV2
            The spectrum containing the region definitions.
        region_id_col : str
            Name of the column containing region IDs (e.g., ``'abs_ids'``).
        series_map : dict
            Mapping ``{region_id: series_name}``.
        z_map : dict
            Mapping ``{region_id: redshift}``.

        Returns
        -------
        SystemListV2
            A new instance populated with the identified components.
        """
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
        """
        Return a new SystemList with properties of a specific component updated.

        Parameters
        ----------
        uuid : str
            The UUID of the component to update.
        **changes : dict
            Keyword arguments for the fields to change (e.g. ``z=2.5``, ``logN=14.0``).

        Returns
        -------
        SystemListV2
            A new instance with the updated component.
        """
        from astrocook.core.structures import ComponentDataV2

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
            func='voigt', 
            series=series # Simply store 'CIV'
        )
        
        # Standard append logic (no expansion, no linking)
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
        """
        Retrieve a component object by its UUID.

        Parameters
        ----------
        uuid : str
            The unique identifier of the component.

        Returns
        -------
        ComponentDataV2 or None
            The component object, or None if not found.
        """
        for c in self.components:
            if c.uuid == uuid: return c
        return None

    def update_component(self, uuid: str, **changes) -> 'SystemListV2':
        """
        Return a new SystemList with properties of a specific component updated.

        Parameters
        ----------
        uuid : str
            The UUID of the component to update.
        **changes : dict
            Keyword arguments for the fields to change (e.g. ``z=2.5``, ``logN=14.0``).

        Returns
        -------
        SystemListV2
            A new instance with the updated component.
        """
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


    def update_components(self, uuids: List[str], **changes) -> 'SystemListV2':
        """
        Batch update parameters for multiple components in a single pass.

        This method is significantly faster than calling :meth:`update_component` 
        repeatedly because it creates a new immutable state only once after applying 
        all changes.

        Parameters
        ----------
        uuids : list of str
            A list of component UUIDs to update.
        **changes : dict
            Keyword arguments specifying the fields to update (e.g., ``z=2.5``).
            Changes are applied identically to all specified components.

        Returns
        -------
        SystemListV2
            A new instance with the components updated.
        """
        # Convert list to set for O(1) lookups
        target_uuids = set(uuids)
        new_components = []
        
        # Single pass through the list
        for c in self.components:
            if c.uuid in target_uuids:
                # Apply changes only to valid fields
                valid_fields = {f.name for f in dataclasses.fields(c)}
                filtered_changes = {k: v for k, v in changes.items() if k in valid_fields}
                new_c = dataclasses.replace(c, **filtered_changes)
                new_components.append(new_c)
            else:
                new_components.append(c)

        new_data = dataclasses.replace(self._data, components=new_components)
        return SystemListV2(new_data)
    

    def update_constraint(self, uuid: str, param: str, 
                          is_free: Optional[bool] = None,
                          expression: Optional[str] = None,
                          target_uuid: Optional[str] = None) -> 'SystemListV2':
        """
        Update the constraint for a specific component parameter.

        Parameters
        ----------
        uuid : str
            The component UUID.
        param : str
            The parameter name (e.g., ``'z'``, ``'logN'``).
        is_free : bool, optional
            If True, parameter varies freely. If False, it is fixed or linked.
        expression : str, optional
            Mathematical string for linking (e.g. ``"p['target_uuid'].z"``).
        target_uuid : str, optional
            UUID of the component this parameter is linked to.

        Returns
        -------
        SystemListV2
            A new instance with the updated constraints.
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
        
        from astrocook.core.structures import ParameterConstraintV2
        
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
            # Default btur to Frozen, others to Free
            default_free = False if param == 'btur' else True
            base_constraint = ParameterConstraintV2(is_free=default_free)

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

    def delete_component(self, uuid: str = None, uuids: list = None) -> 'SystemListV2':
        """
        Remove component(s) and clean up constraints linking to them.

        Parameters
        ----------
        uuid : str, optional
            Single UUID to delete.
        uuids : list, optional
            List of UUIDs to delete.

        Returns
        -------
        SystemListV2
            A new instance with components removed and constraints cleaned.
        """
        import dataclasses
        from astrocook.core.structures import ParameterConstraintV2

        # 1. Normalize input to a Set of targets
        targets = set()
        if uuid: targets.add(uuid)
        if uuids: targets.update(uuids)
        
        if not targets:
            return self

        # 2. Remove the components
        new_components = [c for c in self.components if c.uuid not in targets]
        
        # 3. Clean up Hanging Links (Sanitize Constraints)
        # Scan all constraints. If any target ANY deleted UUID, reset them to Free.
        
        current_v2_map = {k: v.copy() for k, v in self._data.v2_constraints_map.items()}
        current_parsed = dict(self._data.parsed_constraints)
        
        to_reset = []
        
        for u_src, params in current_v2_map.items():
            for p_name, constr in params.items():
                
                # Check for direct target match in the set of deleted UUIDs
                if constr.target_uuid in targets:
                    to_reset.append((u_src, p_name))
                
                # Check for expression reference (if expression contains any deleted UUID)
                elif constr.expression:
                    # Quick check: is any deleted UUID inside the expression string?
                    # (A regex would be safer, but uuid strings are usually unique enough)
                    if any(deleted_uuid in constr.expression for deleted_uuid in targets):
                        to_reset.append((u_src, p_name))

        for u_src, p_name in to_reset:
            logging.info(f"Removing broken link from {u_src}.{p_name} (target deleted).")
            # Reset to Free
            clean_constr = ParameterConstraintV2(is_free=True, expression=None, target_uuid=None)
            
            # Update Nested Map
            current_v2_map[u_src][p_name] = clean_constr
            
            # Update Parsed Map
            current_parsed[(u_src, p_name)] = clean_constr

        # 4. Create new data
        new_data = dataclasses.replace(self._data, 
                                       components=new_components,
                                       parsed_constraints=current_parsed,
                                       v2_constraints_map=current_v2_map)
        return SystemListV2(new_data)

    def get_connected_group(self, seed_uuids: List[str], max_depth: int = 2) -> Set[str]:
        """
        Find a set of components that must be fitted together.

        This method traverses relationships to find all "fluid neighbors" of
        the seed components (via links, spectral overlap, or kinematic proximity).

        Parameters
        ----------
        seed_uuids : list of str
            The UUIDs of the components initially selected for fitting.
        max_depth : int, optional
            How many degrees of separation to traverse. Defaults to ``2``.

        Returns
        -------
        set of str
            The full set of UUIDs in the connected group.
        """
        if not seed_uuids: return set()

        c_kms = 299792.458
        # 1. Map components and retrieve constraints from the data core
        comp_map = {c.uuid: c for c in self.components}
        constraints_map = self._data.v2_constraints_map
        group_uuids = set(seed_uuids)

        # --- STEP 1: Build a Spatial Wavelength Index ---
        # We map every individual line of every component to its observed wavelength
        flat_index = []
        for c in self.components:
            for w_rest in get_all_w_rest(c.series):
                w_obs = w_rest * (1 + c.z)
                flat_index.append((w_obs, c.uuid))

        # Sort for binary search efficiency
        flat_index.sort(key=lambda x: x[0])
        all_w_obs = np.array([x[0] for x in flat_index])
        all_uuids = [x[1] for x in flat_index]

        # --- STEP 2: Transitive Closure (Friends-of-Friends) ---
        for _ in range(max_depth):
            added_count = 0
            current_members = [comp_map[u] for u in group_uuids if u in comp_map]
            if not current_members: break

            # Identify candidate UUIDs nearby in observed wavelength space
            candidate_uuids = set()
            for m in current_members:
                for w_rest in get_all_w_rest(m.series):
                    w_obs = w_rest * (1 + m.z)
                    # Use a 500 km/s search window to find potential candidates
                    dv_ang = (500.0 / c_kms) * w_obs

                    idx_start = np.searchsorted(all_w_obs, w_obs - dv_ang)
                    idx_end = np.searchsorted(all_w_obs, w_obs + dv_ang)

                    for i in range(idx_start, idx_end):
                        candidate_uuids.add(all_uuids[i])

            # --- [FIX] Explicit Link Traversal ---
            # Components linked via constraints MUST be in the same group, 
            # even if they are far apart in wavelength (e.g. CIV and SiIV).
            for m in current_members:
                if m.uuid in constraints_map:
                    for constr in constraints_map[m.uuid].values():
                        if hasattr(constr, 'target_uuid') and constr.target_uuid:
                            candidate_uuids.add(constr.target_uuid)
            
            # Also check if anyone links TO our current members (Inverse search)
            for cand_uuid, c_map in constraints_map.items():
                for constr in c_map.values():
                    if hasattr(constr, 'target_uuid') and constr.target_uuid in group_uuids:
                        candidate_uuids.add(cand_uuid)

            # --- STEP 3: Detailed Physics and Constraint Validation ---
            for cand_uuid in candidate_uuids:
                if cand_uuid in group_uuids: continue

                other = comp_map[cand_uuid]
                # Check proximity and constraints using the updated signature
                is_connected = False
                for member in current_members:
                    if are_overlapping(member, other, constraints_map):
                        is_connected = True
                        break
                    
                if is_connected:
                    group_uuids.add(other.uuid)
                    added_count += 1

            if added_count == 0: break

        return group_uuids

    def _sanitize_after_fit(self, systs: 'SystemListV2') -> 'SystemListV2':
        """
        Removes the 'active_components' mask from a system list.
        This ensures that when we generate models from it, ALL components are included.
        """
        if systs.constraint_model:
            systs.constraint_model.set_active_components(None)
        return systs

    # --- HELPER 3: Compute Fixed AIC ---
    def _compute_fixed_aic(self, spec: 'SpectrumV2', systs: 'SystemListV2', 
                           static_mask: np.ndarray, k_free: int) -> float:
        """
        Calculates AIC on a FIXED set of pixels (static_mask).
        
        AIC = Chi2 + 2*k
        Chi2 = Sum( (Flux - Model)^2 / Error^2 ) masked
        """
        from astrocook.fitting.voigt_fitter import VoigtFitterV2
        
        # 1. Generate Full Model (No masking allowed here!)
        eval_fitter = VoigtFitterV2(spec, systs)
        _, model_flux = eval_fitter.compute_model_flux()
        
        # 2. Calculate Chi2 on the static mask
        y_win = spec.y.value[static_mask]
        mod_win = model_flux[static_mask]
        dy_win = spec.dy.value[static_mask]
        
        chi2_static = np.sum(((y_win - mod_win) / dy_win)**2)
        
        # 3. Calculate AIC
        aic = chi2_static + 2 * k_free
        return aic
    
    def merge(self, other: 'SystemListV2', copy_uuids: bool = False) -> 'SystemListV2':
        """
        Merge components from another system list into this one.

        Parameters
        ----------
        other : SystemListV2
            The source list to copy from.
        copy_uuids : bool, optional
            If True, keeps original UUIDs (risky if merging into self).
            If False, regenerates UUIDs to avoid conflicts. Defaults to ``False``.

        Returns
        -------
        SystemListV2
            A new merged system list.
        """
        from astrocook.core.structures import ComponentDataV2
        import uuid

        current_comps = list(self.components)
        
        # Determine next ID for V1 compatibility
        next_id = 1
        if current_comps:
            next_id = max(c.id for c in current_comps) + 1

        for c in other.components:
            # Prepare the arguments we want to change
            changes = {'id': next_id}
            
            if not copy_uuids:
                changes['uuid'] = str(uuid.uuid4())
            
            # [FIX] Apply changes using replace() instead of direct assignment
            new_c = dataclasses.replace(c, **changes)
            
            current_comps.append(new_c)
            next_id += 1

        # Create new data core
        new_data = dataclasses.replace(self._data, components=current_comps)
        
        return SystemListV2(new_data)
    
    def optimize_hierarchy(self, spec: 'SpectrumV2', uuid_seed: str, 
                           max_components: int = 5, threshold_sigma: float = 2.5,
                           aic_penalty: float = 0.0, z_window_kms: float = 100.0,
                           min_dv: float = 10.0,  # <--- NEW PARAMETER
                           group_depth: int = 2, patience: int = 2) -> 'SystemListV2':
        """
        Iteratively adds components to a system to minimize residuals using the
        Akaike Information Criterion (AIC).

        This method performs the following steps:
        
        1.  **Baseline Fit:** Fits the seed component (and its connected group) to establish a baseline AIC.
        2.  **Residual Scan:** Calculates the residuals (Flux - Model) over the specified velocity window.
        3.  **Candidate Detection:** Identifies the strongest un-modeled absorption feature (negative residual).
        4.  **Trial:** If the feature exceeds ``threshold_sigma``, adds a new component at that position.
        5.  **Evaluation:** Fits the new configuration. If AIC improves (decreases) by ``aic_penalty``, 
            the new component is kept.
        6.  **Termination:** Repeats until ``max_components`` is reached or ``patience`` is exhausted.

        Parameters
        ----------
        spec : SpectrumV2
            The spectrum data used for fitting and residual calculation.
        uuid_seed : str
            UUID of the seed component. This defines the series (e.g. CIV) and the center of the 
            velocity search window.
        max_components : int, optional
            Maximum number of additional components to attempt adding (default: 5).
        threshold_sigma : float, optional
            Significance level (in standard deviations) of the residual dip required to trigger 
            a new component trial (default: 2.5).
        aic_penalty : float, optional
            Minimum AIC improvement required to accept a new component (default: 0.0).
            Higher values favor simpler models (Occam's razor).
        z_window_kms : float, optional
            The velocity window size (km/s) around the seed component to search for residuals (default: 100.0).
        min_dv : float, optional
            Minimum velocity separation (km/s) required between a candidate and 
            existing components. Lower this to de-blend close doublets (default: 10.0).
        group_depth : int, optional
            Depth for finding connected neighbors during the fit (default: 2).
        patience : int, optional
            Number of consecutive trials allowed without AIC improvement before stopping (default: 2).

        Returns
        -------
        SystemListV2
            A new system list containing the optimized components.
        """
        from astrocook.fitting.voigt_fitter import VoigtFitterV2
        from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS

        # 1. Setup Static Window
        current_systs = self
        target_comp = current_systs.get_component_by_uuid(uuid_seed)
        if not target_comp: return self

        series_name = target_comp.series
        lam_target_ang = 1215.67
        if series_name in ATOM_DATA: lam_target_ang = ATOM_DATA[series_name]['wave']
        elif series_name in STANDARD_MULTIPLETS: lam_target_ang = ATOM_DATA[STANDARD_MULTIPLETS[series_name][0]]['wave']
        
        lam_target_nm = lam_target_ang / 10.0
        c_kms = 299792.458

        # Define STATIC MASK for AIC comparison
        lam_obs_center = lam_target_nm * (1 + target_comp.z)
        dv_arr = c_kms * (spec.x.value - lam_obs_center) / lam_obs_center
        static_mask = (np.abs(dv_arr) < z_window_kms) & np.isfinite(spec.y.value) & (spec.dy.value > 0)
        
        if np.sum(static_mask) == 0:
            logging.warning("optimize_hierarchy: Window contains no valid data.")
            return self

        # [FIX] Initialize Taboo Mask (Rejected Pixels)
        rejected_mask = np.zeros_like(spec.x.value, dtype=bool)

        # 2. Local Fit Orchestrator
        def fit_and_evaluate(systs_obj, target_uuids):
            with systs_obj.fitting_context(target_uuids, group_depth=group_depth):
                fitter = VoigtFitterV2(spec, systs_obj)
                new_systs, _, res = fitter.fit(max_nfev=1000, z_window_kms=z_window_kms, verbose=0)
                
                # Cleanup
                if new_systs.constraint_model:
                    new_systs.constraint_model.set_active_components(None)
                
                if not res or not res.success:
                    return float('inf'), systs_obj 

                # Fixed AIC
                eval_fitter = VoigtFitterV2(spec, new_systs)
                _, model_flux = eval_fitter.compute_model_flux()
                
                y_win = spec.y.value[static_mask]
                mod_win = model_flux[static_mask]
                dy_win = spec.dy.value[static_mask]
                
                chi2_static = np.sum(((y_win - mod_win) / dy_win)**2)
                k = len(res.x)
                aic = chi2_static + 2 * k
                
                return aic, new_systs

        # 3. Baseline Fit
        logging.info(f"Optimizing {series_name} hierarchy...")
        best_aic, best_systs = fit_and_evaluate(current_systs, [uuid_seed])
        logging.info(f"Baseline AIC: {best_aic:.1f}")
        
        current_systs = best_systs
        trials_without_imp = 0
        
        # 4. Iteration Loop
        for i in range(max_components):
            # A. Calculate Residuals
            eval_fitter = VoigtFitterV2(spec, best_systs)
            _, model_phys = eval_fitter.compute_model_flux()
            
            with np.errstate(divide='ignore', invalid='ignore'):
                sigma_resid = (spec.y.value - model_phys) / spec.dy.value
            
            # Apply Masks (Window + Taboo)
            sigma_resid[~static_mask] = 0.0
            sigma_resid[rejected_mask] = 0.0 # Ignore previously failed spots!
            
            # B. Find Dip
            min_idx = np.argmin(sigma_resid)
            min_val = sigma_resid[min_idx]
            
            if min_val > -threshold_sigma:
                logging.info(f"Converged. Max resid {min_val:.1f}s > {threshold_sigma}")
                break
                
            # C. Create Candidate
            z_new = spec.x.value[min_idx] / lam_target_nm - 1.0
            
            # D. Proximity Check
            too_close = False
            for c in best_systs.components:
                if c.series == series_name:
                    dz = abs(c.z - z_new)
                    dv_sep = c_kms * dz / (1+z_new)
                    if dv_sep < min_dv: 
                        too_close = True; break
            
            if too_close:
                logging.warning(f"Iter {i+1}: Candidate too close ({dv_sep:.1f} < {min_dv} km/s). Masking and retrying.")
                # [FIX] Mask this spot and CONTINUE (Do not Break!)
                sl = slice(max(0, min_idx-5), min(len(spec.x.value), min_idx+6))
                rejected_mask[sl] = True
                continue 

            logging.info(f"Iter {i+1}: Candidate at {z_new:.5f} (Sigma={min_val:.1f})")
            
            # E. Add & Fit
            cand_systs = best_systs.add_component(series=series_name, z=z_new, logN=13.0, b=15.0)
            new_uuid = cand_systs.components[-1].uuid
            
            trial_aic, trial_systs = fit_and_evaluate(cand_systs, [uuid_seed, new_uuid])
            
            # F. Compare
            delta = trial_aic - best_aic
            if delta < -aic_penalty:
                logging.info(f"  -> Accepted. AIC {best_aic:.1f} -> {trial_aic:.1f} (Delta: {delta:.1f})")
                best_aic = trial_aic
                best_systs = trial_systs
                trials_without_imp = 0
            else:
                logging.info(f"  -> Rejected. AIC Delta {delta:.1f}")
                trials_without_imp += 1
                
                # [FIX] Mask this spot so we don't loop forever
                sl = slice(max(0, min_idx-10), min(len(spec.x.value), min_idx+11))
                rejected_mask[sl] = True
            
            if trials_without_imp >= patience:
                logging.info("Patience exhausted.")
                break
        
        return best_systs
    
    def remove_in_region(self, obs_min: au.Quantity, obs_max: au.Quantity) -> 'SystemListV2':
        """
        Remove components whose center falls within the observed wavelength range.

        Parameters
        ----------
        obs_min, obs_max : Quantity
            The wavelength bounds.

        Returns
        -------
        SystemListV2
            A new instance with components removed.
        """
        surviving_comps = []
        deleted_count = 0
        
        for c in self.components:
            # Resolve Rest Wavelength
            lam_rest = None
            if c.series in ATOM_DATA:
                lam_rest = ATOM_DATA[c.series]['wave']
            elif c.series in STANDARD_MULTIPLETS:
                 # Use primary line for location check
                 prim = STANDARD_MULTIPLETS[c.series][0]
                 if prim in ATOM_DATA:
                     lam_rest = ATOM_DATA[prim]['wave']
            
            should_delete = False
            if lam_rest:
                # Convert to Angstrom Quantity for comparison
                lam_obs = lam_rest * (1 + c.z) * au.Angstrom
                if obs_min <= lam_obs <= obs_max:
                    should_delete = True
            
            if should_delete:
                deleted_count += 1
            else:
                surviving_comps.append(c)
        
        if deleted_count == 0:
            return self

        logging.info(f"Smart clean: Removed {deleted_count} components in region {obs_min:.1f}-{obs_max:.1f}")
        new_data = dataclasses.replace(self._data, components=surviving_comps)
        return SystemListV2(new_data)

    def filter_by_criteria(self, condition_func) -> 'SystemListV2':
        """
        Remove components whose center falls within the observed wavelength range.

        Parameters
        ----------
        obs_min, obs_max : Quantity
            The wavelength bounds.

        Returns
        -------
        SystemListV2
            A new instance with components removed.
        """
        valid_comps = []
        deleted_count = 0
        
        for c in self.components:
            if condition_func(c):
                valid_comps.append(c)
            else:
                deleted_count += 1
                
        if deleted_count == 0:
            return self
            
        logging.info(f"Filtered {deleted_count} components.")
        new_data = dataclasses.replace(self._data, components=valid_comps)
        return SystemListV2(new_data)
    
    def filter_by_range(self, xmin_nm: float, xmax_nm: float) -> 'SystemListV2':
        """
        Return a new SystemListV2 containing only components within the range.

        Parameters
        ----------
        xmin_nm, xmax_nm : float
            Wavelength range in nanometers.

        Returns
        -------
        SystemListV2
            A new instance with filtered components.
        """
        import numpy as np
        import astropy.units as au
        import dataclasses
        from astrocook.core.atomic_data import xem_d, STANDARD_MULTIPLETS
        
        valid_components = []
        
        # Ensure limits are valid
        if not np.isfinite(xmin_nm) or not np.isfinite(xmax_nm):
            return self

        for c in self.components:
            is_in_range = False
            
            # 1. Identify which lines this component represents
            lines_to_check = []
            
            # Case A: It's a specific line (e.g., "CIV 1548")
            if c.series in xem_d:
                lines_to_check.append(c.series)
            
            # Case B: It's a Multiplet Tag (e.g., "CIV")
            elif c.series in STANDARD_MULTIPLETS:
                lines_to_check.extend(STANDARD_MULTIPLETS[c.series])
                
            # Case C: Try normalized name (e.g., "CIV 1548" vs "CIV_1548")
            else:
                norm = c.series.replace(' ', '_')
                if norm in xem_d:
                    lines_to_check.append(norm)

            # 2. Check if ANY line falls in range
            if not lines_to_check:
                # If we can't identify the line physics, we default to KEEPING it
                # (Safe fallback for custom markers)
                is_in_range = True
            else:
                for line in lines_to_check:
                    if line in xem_d:
                        lam_rest_nm = xem_d[line].to(au.nm).value
                        lam_obs = lam_rest_nm * (1.0 + c.z)
                        
                        # Check bounds (1.0 nm padding)
                        if (xmin_nm - 1.0) <= lam_obs <= (xmax_nm + 1.0):
                            is_in_range = True
                            break
            
            if is_in_range:
                valid_components.append(c)

        # Return new instance with filtered data
        new_data = dataclasses.replace(self._data, components=valid_components)
        return SystemListV2(new_data)