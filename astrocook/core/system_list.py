from astrocook._renamed import renamed
from astropy.table import Table, Column
import astropy.units as au
import contextlib
import dataclasses
import logging
import numpy as np
import ast
from typing import TYPE_CHECKING, Any, Collection, cast, Dict, List, Optional, Set, Tuple

from astrocook.fitting.voigt_model import LINK_ERRORS, LINK_FUNCTIONS, VoigtModelConstraint
from astrocook.core.structures import ComponentData, SystemListData, ParameterConstraint
from astrocook.core.atomic_data import ATOM_DATA, is_multiplet, series_lines

if TYPE_CHECKING:
    from astrocook.core.spectrum import Spectrum

logger = logging.getLogger(__name__)

def get_all_w_rest(series_name: str) -> List[float]:
    """
    Returns a list of rest wavelengths in Angstroms for a given series name.
    Handles single lines, standard multiplets, and custom comma-separated lists.
    """
    # 1. Check if it's a known single transition
    if series_name in ATOM_DATA:
        return [ATOM_DATA[series_name]['wave']]
    
    # 2. Check if it's a standard multiplet tag (e.g., 'CIV')
    if is_multiplet(series_name):
        waves = []
        for line_name in series_lines(series_name):
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

#: Blend threshold in units of the instrumental FWHM (findings register F14):
#: two lines are grouped when their centres are closer than
#: ``FWHM_BLEND_FACTOR`` times the FWHM, i.e. when their instrumental
#: profiles overlap above half maximum (Guido's choice, October 2026). The
#: thresholds ``4 b`` and 30 km/s are kept as floors (they were the whole
#: rule before F14, and the only one when the resolution is unknown).
FWHM_BLEND_FACTOR = 1.0


def instrumental_fwhm(spec: Optional['Spectrum']):
    """Instrumental FWHM (km/s) of a spectrum as a function of the observed
    wavelength in Angstrom, from its ``resol`` column (#526).

    Returns ``None`` if there is no spectrum or it has no resolution; the
    returned function gives NaN outside the range of the spectrum.
    """
    if spec is None:
        return None
    from astrocook.core.resolution import C_KMS

    r = spec.resolution()
    ok = np.isfinite(r)
    if not ok.any():
        return None
    try:
        x_ang = (spec.x.value[ok] * spec.x.unit).to_value(au.Angstrom)
    except au.UnitConversionError:
        return None
    order = np.argsort(x_ang)
    x_ang, fwhm = x_ang[order], C_KMS / r[ok][order]

    def fwhm_at(w_ang):
        return np.interp(w_ang, x_ang, fwhm, left=np.nan, right=np.nan)
    return fwhm_at


def _component_fwhm(comp: ComponentData, w_obs_ang: float, fwhm_at) -> float:
    """FWHM (km/s) for a line of a component: its own resolution if set
    (an override, always R; the same priority as the fitter), otherwise the
    spectrum's at the line; NaN if unknown."""
    if comp.resol is not None and comp.resol > 0:
        return 299792.458 / comp.resol
    if fwhm_at is None:
        return np.nan
    return float(fwhm_at(w_obs_ang))


def are_overlapping(c1: ComponentData, c2: ComponentData, constraints_map: dict,
                    fwhm_at=None) -> bool:
    """
    Checks if two components should be grouped together for fitting.

    ``fwhm_at`` (from :func:`instrumental_fwhm`) gives the instrumental FWHM
    at a wavelength: two lines are blended when their centres are closer
    than ``max(4 b, 30 km/s, FWHM_BLEND_FACTOR * FWHM)`` (F14), with the
    larger FWHM of the two lines. Without it (or where the resolution is
    unknown) the threshold is ``max(4 b, 30 km/s)``.
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
    # We use a threshold based on the thermal width (b-parameter) or a 30 km/s
    # floor, raised to the instrumental FWHM where it is larger (F14).
    base_kms = max(4.0 * max(c1.b, c2.b), 30.0)

    def threshold(w1, w2):
        fwhm = np.nanmax([_component_fwhm(c1, w1, fwhm_at),
                          _component_fwhm(c2, w2, fwhm_at), -np.inf])
        return max(base_kms, FWHM_BLEND_FACTOR * fwhm) if np.isfinite(fwhm) else base_kms

    obs1 = [w * (1 + c1.z) for w in get_all_w_rest(c1.series)]
    obs2 = [w * (1 + c2.z) for w in get_all_w_rest(c2.series)]
    
    for w1 in obs1:
        for w2 in obs2:
            dv = c_kms * abs(w1 - w2) / ((w1 + w2) / 2.0)
            if dv < threshold(w1, w2):
                return True

    # 3. KINEMATIC PROXIMITY (Same Species)
    # If they are the same species and close in velocity, group them
    # (threshold at their first transition).
    if c1.series == c2.series:
        dz = abs(c1.z - c2.z)
        dv_kin = c_kms * dz / (1 + c1.z)
        threshold_kms = (threshold(obs1[0], obs2[0]) if obs1 and obs2 else base_kms)
        if dv_kin < threshold_kms:
            return True
            
    return False

# --- 3. System List API Layer (The orchestrator) ---
#: Parameters whose change clears the flags of the last fit (#517)
_FITTED_PARAMETERS = {'z', 'logN', 'b', 'btur'}


def _valid_changes(comp, changes: dict) -> dict:
    """The changes that apply to a component; editing a fitted parameter
    clears the flags left by the last fit (#517), unless they are set too."""
    valid_fields = {f.name for f in dataclasses.fields(comp)}
    out = {k: v for k, v in changes.items() if k in valid_fields}
    if _FITTED_PARAMETERS & out.keys() and 'flags' not in out:
        out['flags'] = ()
    return out


class SystemList:
    """
    API layer for managing component lists, constraints, and fitting.

    This class wraps the immutable :class:`~astrocook.core.structures.SystemListData`
    and provides methods to add, update, delete, and group components. It also
    initializes the constraint model used by the fitting engine.

    Parameters
    ----------
    data : SystemListData
        The core data container holding the list of components and constraint maps.
    """
    
    def __init__(self, data: SystemListData):
        self._data = data
        
        # 1. Sanitize & Hydrate data immediately
        # (Fixes JSON artifacts and converts dicts -> objects)
        self._ensure_clean_data()
        
        # Created by the call below, so that it is never None when used
        self.constraint_model: VoigtModelConstraint = None  # type: ignore[assignment]
        self._initialize_constraint_model()

    def _ensure_clean_data(self):
        """
        Robustly initializes data core.
        1. Sanitize artifacts in ID map.
        2. HYDRATE 'v2_constraints_map' (Dict -> ParameterConstraint).
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
                        hydrated_v2_map[u][p] = ParameterConstraint(
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
            logger.debug("SystemList: Hydrated constraints from clean JSON map.")
        
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
                    except (ValueError, SyntaxError, TypeError, MemoryError, RecursionError):
                        pass  # not a literal: the key is kept as it is
                
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

                    val_obj = ParameterConstraint(
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
            logger.debug("Old constraints converted to the current format.")

        if dirty or constraints_dirty:
            self._data = dataclasses.replace(
                self._data, 
                v1_id_to_uuid_map=clean_id_map,
                parsed_constraints=clean_parsed,
                v2_constraints_map=hydrated_v2_map
            )

    def _initialize_constraint_model(self):
        # A failure reaches the caller: an empty fallback model lost the
        # links and the frozen parameters at the next save (the archive
        # reads them from here), and None made the save crash (phase 6)
        self.constraint_model = VoigtModelConstraint(self)

    @property
    def components(self) -> List[ComponentData]:
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
    def parsed_constraints(self) -> Dict[Tuple[str, str], ParameterConstraint]:
        return self._data.parsed_constraints
    
    @property
    def v2_constraints_for_save(self) -> Dict[str, Dict[str, Any]]:
        if self.constraint_model:
            return self.constraint_model.v2_constraints_by_uuid
        return {}
    
    # --- Context Manager for Safe Fitting ---
    @contextlib.contextmanager
    def fitting_context(self, active_uuids: Optional[List[str]] = None, group_depth: int = 2,
                        spec: Optional['Spectrum'] = None):
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
        spec : Spectrum, optional
            The spectrum to be fitted: its resolution sets the blend
            threshold of the grouping (see :meth:`get_connected_group`).
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
            self.constraint_model.set_active_components(active_uuids, group_depth=group_depth,
                                                        spec=spec)
            yield
        finally:
            # Robust cleanup
            try:
                # Restore exactly the previous set (no further expansion)
                self.constraint_model.set_active_components(original_active, group_depth=0)
            except (KeyError, ValueError, TypeError, AttributeError, IndexError) as e:
                # A failure here would leave the components outside the
                # group frozen for the next fits
                logger.error(f"Could not restore the active components after the fit "
                             f"({e}): all the components are made active.")
                try:
                    self.constraint_model.set_active_components(None)
                except (KeyError, ValueError, TypeError, AttributeError, IndexError):
                    # Last resort: direct attribute clear
                    self.constraint_model._active_uuids = None

    def add_components_from_regions(self, spec: 'Spectrum', region_id_col: str,
                                  series_map: Dict[int, str], z_map: Dict[int, float]) -> 'SystemList':
        """
        Populate the list with initial components based on identified regions.

        Parameters
        ----------
        spec : Spectrum
            The spectrum containing the region definitions.
        region_id_col : str
            Name of the column containing region IDs (e.g., ``'abs_ids'``).
        series_map : dict
            Mapping ``{region_id: series_name}``.
        z_map : dict
            Mapping ``{region_id: redshift}``.

        Returns
        -------
        SystemList
            A new instance populated with the identified components.
        """
        if not spec.has_aux_column(region_id_col):
            logger.warning(f"There is no column '{region_id_col}': no components are added.")
            return self

        new_components = list(self._data.components)
        next_id = 1
        if new_components:
            next_id = max(c.id for c in new_components) + 1

        logger.info(f"Adding {len(series_map)} components from the identified regions.")
        
        for rid, series in series_map.items():
            if rid not in z_map: continue
            z_region = z_map[rid]
            new_comp = ComponentData(
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
        return SystemList(data=new_data)
    
    def add_component(self, series: str, z: float, logN: float = 13.5, 
                      b: float = 10.0, btur: float = 0.0) -> 'SystemList':
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
        SystemList
            A new instance with the updated component.
        """
        from astrocook.core.structures import ComponentData

        current_components = self._data.components
        next_id = 1
        if current_components:
            next_id = max(c.id for c in current_components) + 1
            
        new_comp = ComponentData(
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
        return SystemList(new_data_core)
    
    def get_component_by_uuid(self, uuid: str) -> Optional[ComponentData]:
        """
        Retrieve a component object by its UUID.

        Parameters
        ----------
        uuid : str
            The unique identifier of the component.

        Returns
        -------
        ComponentData or None
            The component object, or None if not found.
        """
        for c in self.components:
            if c.uuid == uuid: return c
        return None

    def update_component(self, uuid: str, **changes) -> 'SystemList':
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
        SystemList
            A new instance with the updated component.
        """
        new_components = []
        found = False
        for c in self.components:
            if c.uuid == uuid:
                new_c = dataclasses.replace(c, **_valid_changes(c, changes))
                new_components.append(new_c)
                found = True
            else:
                new_components.append(c)
        
        if not found:
            logger.warning(f"There is no component {uuid}: nothing is updated.")
            return self

        new_data = dataclasses.replace(self._data, components=new_components)
        return SystemList(new_data)


    def update_components(self, uuids: List[str], **changes) -> 'SystemList':
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
        SystemList
            A new instance with the components updated.
        """
        # Convert list to set for O(1) lookups
        target_uuids = set(uuids)
        new_components = []
        
        # Single pass through the list
        for c in self.components:
            if c.uuid in target_uuids:
                # Apply changes only to valid fields
                new_c = dataclasses.replace(c, **_valid_changes(c, changes))
                new_components.append(new_c)
            else:
                new_components.append(c)

        new_data = dataclasses.replace(self._data, components=new_components)
        return SystemList(new_data)
    

    def _check_link(self, uuid: str, param: str, expression: str) -> None:
        """Raise ``ValueError`` if the expression of a link cannot be evaluated.

        It is evaluated with the current values of the components, as the
        fit evaluates it (``VoigtModelConstraint.map_p_free_to_full``).
        """
        class _Params:
            __slots__ = ('z', 'logN', 'b', 'btur')

            def __init__(self, c):
                self.z, self.logN, self.b, self.btur = c.z, c.logN, c.b, c.btur

        p = {c.uuid: _Params(c) for c in self.components}
        try:
            code = compile(expression, f"<link {uuid}.{param}>", 'eval')
            float(eval(code, {'p': p, **LINK_FUNCTIONS}))
        except (SyntaxError,) + LINK_ERRORS as e:
            raise ValueError(f"The link of {param} ('{expression}') cannot be "
                             f"evaluated ({type(e).__name__}: {e}).") from e

    def update_constraint(self, uuid: str, param: str, 
                          is_free: Optional[bool] = None,
                          expression: Optional[str] = None,
                          target_uuid: Optional[str] = None) -> 'SystemList':
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
        SystemList
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
                 existing = current_parsed.get(cast(Any, legacy_key))  # legacy keys use the int id
        
        from astrocook.core.structures import ParameterConstraint
        
        # 2. Determine Base
        if existing:
            if isinstance(existing, dict):
                 base_constraint = ParameterConstraint(
                     is_free=existing.get('is_free', True), 
                     expression=existing.get('expression'),
                     target_uuid=existing.get('target_uuid')
                 )
            else:
                 base_constraint = existing
        else:
            # Default btur to Frozen, others to Free
            default_free = False if param == 'btur' else True
            base_constraint = ParameterConstraint(is_free=default_free)

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

        if expression and is_free is not True:
            # A new link is checked here, so that a wrong one is refused
            # with a message instead of freezing the parameter (phase 6).
            # Also with ``is_free`` left out: the fit treats an expression as
            # a link either way, and the example of the description of the
            # recipe, ``p["z"]``, was stored unchecked and kept the parameter
            # frozen (F117, trial 3 of the MCP server)
            self._check_link(uuid, param, expression)

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
        
        return SystemList(new_data)

    def delete_component(self, uuid: Optional[str] = None, uuids: Optional[list] = None) -> 'SystemList':
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
        SystemList
            A new instance with components removed and constraints cleaned.
        """
        import dataclasses
        from astrocook.core.structures import ParameterConstraint

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
            logger.info(f"The link of {p_name} of component {u_src} pointed to a "
                        "deleted component: the parameter is now free.")
            # Reset to Free
            clean_constr = ParameterConstraint(is_free=True, expression=None, target_uuid=None)
            
            # Update Nested Map
            current_v2_map[u_src][p_name] = clean_constr
            
            # Update Parsed Map
            current_parsed[(u_src, p_name)] = clean_constr

        # 4. Create new data
        new_data = dataclasses.replace(self._data, 
                                       components=new_components,
                                       parsed_constraints=current_parsed,
                                       v2_constraints_map=current_v2_map)
        return SystemList(new_data)

    def get_connected_group(self, seed_uuids: Collection[str], max_depth: int = 2,
                            spec: Optional['Spectrum'] = None) -> Set[str]:
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
        spec : Spectrum, optional
            The spectrum the components are fitted to. Lines closer than its
            instrumental FWHM are blended (findings register F14, see
            :func:`are_overlapping`); without it, only ``max(4 b, 30 km/s)``.

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
        fwhm_at = instrumental_fwhm(spec)

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

                    # --- Explicit Link Traversal ---
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
                    if are_overlapping(member, other, constraints_map, fwhm_at):
                        is_connected = True
                        break
                    
                if is_connected:
                    group_uuids.add(other.uuid)
                    added_count += 1

            if added_count == 0: break

        return group_uuids

    def fit_groups(self, seed_uuids, max_depth: int, spec: Optional['Spectrum'] = None,
                   max_size: Optional[int] = None) -> List[List[str]]:
        """
        The groups of components that are fitted together, as *Refit All*
        fits them (moved here from ``refit_all`` in phase 6-ter, PR 3, so
        that the MCP server shows the same groups).

        The seeds are taken by increasing redshift (then series): the group
        of a seed is :meth:`get_connected_group` (it can take in components
        that are not seeds), and its members are not seeds again. A group
        larger than ``max_size`` is split in chunks of ``max_size``
        components, by redshift.

        Parameters
        ----------
        seed_uuids : iterable of str
            The components to group (e.g. those in the region of the fit).
        max_depth : int
            Depth of the friends-of-friends expansion
            (:meth:`get_connected_group`).
        spec : Spectrum, optional
            The spectrum the components are fitted to (its resolution sets
            which lines are blended).
        max_size : int, optional
            The largest group fitted at once; None: no limit.

        Returns
        -------
        list of list of str
            The groups, each sorted by redshift (then series), in the order
            of their seeds. A component pulled in by two seeds is in both
            groups, as in ``refit_all``.
        """
        comp_map = {c.uuid: c for c in self.components}
        todo = {u for u in seed_uuids if u in comp_map}
        # Seeds in a deterministic order (increasing z): the depth-limited
        # group depends on the seed, and the iteration order of a set of
        # random UUIDs changes from run to run
        seed_order = sorted(todo, key=lambda u: (comp_map[u].z, comp_map[u].series))
        groups = []
        for seed in seed_order:
            if seed not in todo:
                continue
            group = self.get_connected_group([seed], max_depth=max_depth, spec=spec)
            if not group:
                todo.discard(seed)
                continue
            group_list = sorted(group, key=lambda u: (comp_map[u].z, comp_map[u].series))
            if max_size is None or len(group_list) <= max_size:
                groups.append(group_list)
            else:
                by_z = [c.uuid for c in sorted((comp_map[u] for u in group_list),
                                               key=lambda c: c.z)]
                for k in range(0, len(by_z), max_size):
                    groups.append(by_z[k:k + max_size])
                logger.debug(f"  -> Split massive group ({len(group_list)} comps) into "
                             f"{len(by_z)//max_size + 1} chunks.")
            todo.difference_update(group)
        return groups

    def _sanitize_after_fit(self, systs: 'SystemList') -> 'SystemList':
        """
        Removes the 'active_components' mask from a system list.
        This ensures that when we generate models from it, ALL components are included.
        """
        if systs.constraint_model:
            systs.constraint_model.set_active_components(None)
        return systs

    # --- HELPER 3: Compute Fixed AIC ---
    def _compute_fixed_aic(self, spec: 'Spectrum', systs: 'SystemList', 
                           static_mask: np.ndarray, k_free: int) -> float:
        """
        Calculates AIC on a FIXED set of pixels (static_mask).
        
        AIC = Chi2 + 2*k
        Chi2 = Sum( (Flux - Model)^2 / Error^2 ) masked
        """
        from astrocook.fitting.voigt_fitter import VoigtFitter
        
        # 1. Generate Full Model (No masking allowed here!)
        eval_fitter = VoigtFitter(spec, systs)
        _, model_flux = eval_fitter.compute_model_flux()
        
        # 2. Calculate Chi2 on the static mask
        y_win = spec.y.value[static_mask]
        mod_win = model_flux[static_mask]
        dy_win = spec.dy.value[static_mask]
        
        chi2_static = np.sum(((y_win - mod_win) / dy_win)**2)
        
        # 3. Calculate AIC
        aic = chi2_static + 2 * k_free
        return aic
    
    def merge(self, other: 'SystemList', copy_uuids: bool = False) -> 'SystemList':
        """
        Merge components from another system list into this one.

        Parameters
        ----------
        other : SystemList
            The source list to copy from.
        copy_uuids : bool, optional
            If True, keeps original UUIDs (risky if merging into self).
            If False, regenerates UUIDs to avoid conflicts. Defaults to ``False``.

        Returns
        -------
        SystemList
            A new merged system list.
        """
        import uuid

        current_comps = list(self.components)
        
        # Determine next ID for V1 compatibility
        next_id = 1
        if current_comps:
            next_id = max(c.id for c in current_comps) + 1

        for c in other.components:
            # Prepare the arguments we want to change
            changes: Dict[str, Any] = {'id': next_id}
            
            if not copy_uuids:
                changes['uuid'] = str(uuid.uuid4())
            
            # [FIX] Apply changes using replace() instead of direct assignment
            new_c = dataclasses.replace(c, **changes)
            
            current_comps.append(new_c)
            next_id += 1

        # Create new data core
        new_data = dataclasses.replace(self._data, components=current_comps)
        
        return SystemList(new_data)
    
    def optimize_hierarchy(self, spec: 'Spectrum', uuid_seed: str, 
                           max_components: int = 5, threshold_sigma: float = 2.5,
                           aic_penalty: float = 0.0, z_window_kms: float = 100.0,
                           min_dv: float = 10.0,  # <--- NEW PARAMETER
                           group_depth: int = 2, patience: int = 2,
                           trace: Optional[dict] = None,
                           max_nfev: int = 1000) -> 'SystemList':
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
        spec : Spectrum
            The spectrum data used for fitting and residual calculation.
        uuid_seed : str
            UUID of the seed component. This defines the series (e.g. CIV) and the center of the 
            velocity search window.
        max_components : int, optional
            Maximum number of candidates tried, each counting also when rejected or
            too close to a component (default: 5).
        threshold_sigma : float, optional
            Significance level (in standard deviations) of the residual dip required to trigger 
            a new component trial (default: 2.5).
        aic_penalty : float, optional
            Minimum AIC improvement required to accept a new component (default: 0.0).
            Higher values favor simpler models (Occam's razor).
        z_window_kms : float, optional
            Half-width (km/s) of the window around the first line of the seed component,
            where the residuals are searched and the AIC is computed; also passed to the
            fits, where it bounds how far the redshifts can move (default: 100.0).
        min_dv : float, optional
            Minimum velocity separation (km/s) required between a candidate and 
            existing components. Lower this to de-blend close doublets (default: 10.0).
        group_depth : int, optional
            Depth for finding connected neighbors during the fit (default: 2).
        patience : int, optional
            Number of consecutive trials allowed without AIC improvement before stopping (default: 2).
        max_nfev : int, optional
            Maximum number of evaluations of the model in each fit (default:
            1000). A fit that does not converge within it counts as failed
            (AIC infinite): its candidate is rejected (F149).
        trace : dict, optional
            Filled with what the search did (F118, trial 3 of the MCP server):
            ``start_aic``, ``tried`` (one dict per candidate: ``z``,
            ``residual_sigma``, ``outcome``, and ``aic_before``/``aic_after``
            when it was fitted) and ``stopped`` (the reason).

        Returns
        -------
        SystemList
            A new system list containing the optimized components.
        """
        from astrocook.fitting.voigt_fitter import VoigtFitter
        from astrocook.core.atomic_data import ATOM_DATA

        # 1. Setup Static Window
        current_systs = self
        target_comp = current_systs.get_component_by_uuid(uuid_seed)
        if not target_comp: return self

        series_name = target_comp.series
        lam_target_ang = 1215.67
        if series_name in ATOM_DATA: lam_target_ang = ATOM_DATA[series_name]['wave']
        elif is_multiplet(series_name): lam_target_ang = ATOM_DATA[series_lines(series_name)[0]]['wave']
        
        lam_target_nm = lam_target_ang / 10.0
        c_kms = 299792.458

        # Define STATIC MASK for AIC comparison
        lam_obs_center = lam_target_nm * (1 + target_comp.z)
        dv_arr = c_kms * (spec.x.value - lam_obs_center) / lam_obs_center
        static_mask = (np.abs(dv_arr) < z_window_kms) & np.isfinite(spec.y.value) & (spec.dy.value > 0)
        
        if np.sum(static_mask) == 0:
            logger.warning("No valid pixels around the component: the system is "
                           "not optimized.")
            return self

        # [FIX] Initialize Taboo Mask (Rejected Pixels)
        rejected_mask = np.zeros_like(spec.x.value, dtype=bool)

        # 2. Local Fit Orchestrator
        def fit_and_evaluate(systs_obj, target_uuids):
            with systs_obj.fitting_context(target_uuids, group_depth=group_depth, spec=spec):
                fitter = VoigtFitter(spec, systs_obj)
                new_systs, _, res = fitter.fit(max_nfev=max_nfev, z_window_kms=z_window_kms, verbose=0)
                
                # Cleanup
                if new_systs.constraint_model:
                    new_systs.constraint_model.set_active_components(None)
                
                if not res or not res.success:
                    return float('inf'), systs_obj 

                # Fixed AIC
                eval_fitter = VoigtFitter(spec, new_systs)
                _, model_flux = eval_fitter.compute_model_flux()
                
                y_win = spec.y.value[static_mask]
                mod_win = model_flux[static_mask]
                dy_win = spec.dy.value[static_mask]
                
                chi2_static = np.sum(((y_win - mod_win) / dy_win)**2)
                k = len(res.x)
                aic = chi2_static + 2 * k
                
                return aic, new_systs

        # 3. Baseline Fit
        logger.debug(f"Optimizing the components of {series_name}.")
        best_aic, best_systs = fit_and_evaluate(current_systs, [uuid_seed])
        logger.info(f"AIC of the starting model: {best_aic:.1f}.")
        tried: List[Any] = []
        stopped = f"the {max_components} components allowed were tried"
        if trace is not None:
            trace.update(start_aic=float(best_aic), tried=tried)
        
        current_systs = best_systs
        trials_without_imp = 0
        
        # 4. Iteration Loop
        for i in range(max_components):
            # A. Calculate Residuals
            eval_fitter = VoigtFitter(spec, best_systs)
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
                logger.info(f"No residual deeper than {threshold_sigma} sigma "
                            f"(the deepest is {-min_val:.1f} sigma): done.")
                stopped = (f"no residual deeper than threshold_sigma = "
                           f"{threshold_sigma} (the deepest is {-min_val:.1f} sigma)")
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
                logger.debug(f"Iteration {i+1}: candidate {dv_sep:.1f} km/s from a "
                             f"component (less than {min_dv} km/s): masked.")
                # [FIX] Mask this spot and CONTINUE (Do not Break!)
                sl = slice(max(0, min_idx-5), min(len(spec.x.value), min_idx+6))
                rejected_mask[sl] = True
                tried.append({'z': float(z_new), 'residual_sigma': float(-min_val),
                              'outcome': f"{dv_sep:.1f} km/s from a component "
                                         f"(less than min_dv = {min_dv}): not fitted"})
                continue 

            logger.debug(f"Iteration {i+1}: candidate at z = {z_new:.5f} "
                         f"(residual {min_val:.1f} sigma).")
            
            # E. Add & Fit
            cand_systs = best_systs.add_component(series=series_name, z=z_new, logN=13.0, b=15.0)
            new_uuid = cand_systs.components[-1].uuid
            
            trial_aic, trial_systs = fit_and_evaluate(cand_systs, [uuid_seed, new_uuid])
            
            # F. Compare
            delta = trial_aic - best_aic
            if delta < -aic_penalty:
                logger.info(f"Component added at z = {z_new:.5f}: AIC {best_aic:.1f} "
                            f"to {trial_aic:.1f}.")
                tried.append({'z': float(z_new), 'residual_sigma': float(-min_val),
                              'outcome': 'added', 'aic_before': float(best_aic),
                              'aic_after': float(trial_aic)})
                best_aic = trial_aic
                best_systs = trial_systs
                trials_without_imp = 0
            else:
                logger.debug(f"Candidate at z = {z_new:.5f} rejected: AIC changes by "
                             f"{delta:.1f}.")
                tried.append({'z': float(z_new), 'residual_sigma': float(-min_val),
                              'outcome': f"rejected: AIC changes by {delta:.1f}, "
                                         f"less than aic_penalty = {aic_penalty}",
                              'aic_before': float(best_aic),
                              'aic_after': float(trial_aic)})
                trials_without_imp += 1
                
                # [FIX] Mask this spot so we don't loop forever
                sl = slice(max(0, min_idx-10), min(len(spec.x.value), min_idx+11))
                rejected_mask[sl] = True
            
            if trials_without_imp >= patience:
                logger.info(f"{patience} candidates in a row rejected: done.")
                stopped = f"{patience} candidates in a row rejected (patience)"
                break
        
        if trace is not None:
            trace['stopped'] = stopped
        return best_systs
    
    def remove_in_region(self, obs_min: au.Quantity, obs_max: au.Quantity) -> 'SystemList':
        """
        Remove components whose center falls within the observed wavelength range.

        Parameters
        ----------
        obs_min, obs_max : Quantity
            The wavelength bounds.

        Returns
        -------
        SystemList
            A new instance with components removed.
        """
        surviving_comps = []
        deleted_count = 0
        
        for c in self.components:
            # Resolve Rest Wavelength
            lam_rest = None
            if c.series in ATOM_DATA:
                lam_rest = ATOM_DATA[c.series]['wave']
            elif is_multiplet(c.series):
                 # Use primary line for location check
                 prim = series_lines(c.series)[0]
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

        logger.info(f"{deleted_count} components removed between {obs_min:.1f} "
                    f"and {obs_max:.1f}.")
        new_data = dataclasses.replace(self._data, components=surviving_comps)
        return SystemList(new_data)

    def filter_by_criteria(self, condition_func) -> 'SystemList':
        """
        Remove components whose center falls within the observed wavelength range.

        Parameters
        ----------
        obs_min, obs_max : Quantity
            The wavelength bounds.

        Returns
        -------
        SystemList
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
            
        logger.info(f"{deleted_count} components removed.")
        new_data = dataclasses.replace(self._data, components=valid_comps)
        return SystemList(new_data)
    
    def filter_by_range(self, xmin_nm: float, xmax_nm: float) -> 'SystemList':
        """
        Return a new SystemList containing only components within the range.

        Parameters
        ----------
        xmin_nm, xmax_nm : float
            Wavelength range in nanometers.

        Returns
        -------
        SystemList
            A new instance with filtered components.
        """
        import numpy as np
        import astropy.units as au
        import dataclasses
        from astrocook.core.atomic_data import xem_d
        
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
            elif is_multiplet(c.series):
                lines_to_check.extend(series_lines(c.series))
                
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
        return SystemList(new_data)

# Old names (phase 6-quater, PR 6): see astrocook/_renamed.py
__getattr__ = renamed(__name__, {"SystemListV2": "SystemList"})
