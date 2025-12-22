from dataclasses import dataclass, field
from astropy import units as au
from astropy.table import Table
import logging
import numpy as np
from typing import Any, Dict, List, Optional, Tuple
import uuid


@dataclass(frozen=True)
class SessionMetadataV2:
    """
    Core immutable structure for holding constraints, configuration, and history.

    This structure is designed to be serialized to a JSON or YAML file within
    the .acs archive. It separates the lightweight metadata from the heavy
    spectral arrays.

    Parameters
    ----------
    constraints_by_uuid : dict
        A mapping of Component UUIDs to constraint dictionaries.
    v1_reconstruction_data : Any
        Legacy data structure preserved for technical debt logging or V1 reconstruction.
    log_history_json : Any
        The serialized log object (either a V2 dict or V1 JSON format).
    log_type : str
        Type of log stored: ``'v2'``, ``'v1_artifact'``, or ``'v1_legacy'``.
        Defaults to ``'v1_legacy'``.
    component_uuids : list of str
        An ordered list of all component UUIDs in the session.
    """
    
    # Stores constraints linked by Component UUID
    constraints_by_uuid: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    
    # Stores the original system reconstruction data (for technical debt logging)
    v1_reconstruction_data: Any = None 

    # Improved Session History/Log
    # Stores the serialized log object (either a V2 dict or V1 JSON)
    log_history_json: Any = field(default_factory=dict) # Changed default
    # Stores the type ('v2', 'v1_artifact', 'v1_legacy')
    log_type: str = "v1_legacy" # Default to old V1 log
    
    # Stores a list of all component UUIDs in the order they appear
    component_uuids: List[str] = field(default_factory=list)

@dataclass(frozen=True)
class DataColumnV2:
    """
    Immutable container for a single column of data with physical units.

    Parameters
    ----------
    values : np.ndarray
        The numerical data array.
    unit : astropy.units.Unit
        The physical unit associated with the values.
    description : str, optional
        A short description of the column's content.
    """
    values: np.ndarray
    unit: au.Unit
    description: str = ""
    
    @property
    def quantity(self) -> au.Quantity:
        """
        Return the column data as an Astropy Quantity.

        Combines the numerical values and the physical unit.
        """
        return au.Quantity(self.values, self.unit)

    # --- [FIX] Compatibility Aliases ---
    @property
    def value(self) -> np.ndarray:
        """
        Return the raw numerical array.

        Alias for ``.values`` to support legacy/Astropy-style access.
        """
        return self.values

    def __len__(self):
        """Return the number of elements in the column."""
        return len(self.values)
    
    def __array__(self):
        """Return the raw array representation (for NumPy compatibility)."""
        return self.values

@dataclass(frozen=True)
class SpectrumDataV2:
    r"""
    Immutable core container for all spectral data.

    This class holds the "heavy" numerical arrays (flux, wavelength, error)
    and auxiliary columns. It is the data payload for :class:`~astrocook.core.spectrum.SpectrumV2`.

    Parameters
    ----------
    x : DataColumnV2
        The independent variable (wavelength or velocity).
    xmin : DataColumnV2
        Lower bound of the x-bins.
    xmax : DataColumnV2
        Upper bound of the x-bins.
    y : DataColumnV2
        The dependent variable (flux density).
    dy : DataColumnV2
        The 1-sigma uncertainty on ``y``.
    aux_cols : dict of DataColumnV2
        Additional columns (e.g. ``'cont'``, ``'model'``, ``'abs_mask'``).
    meta : dict
        Header keywords and other metadata.
    z_em : float
        Emission redshift of the target object.
    resol : float
        Median resolving power :math:`R = \lambda / \Delta\lambda`.
    """
    
    # Dati Spaziali (Assi X e Binning)
    x: DataColumnV2    # Lunghezza d'onda / Velocità
    xmin: DataColumnV2 # Limite inferiore del bin
    xmax: DataColumnV2 # Limite superiore del bin
    
    # Dati Verticali (Flusso e Errore)
    y: DataColumnV2    # Flusso / Densità di flusso
    dy: DataColumnV2   # Errore sul Flusso

    # Colonne Aggiuntive (Cont, Resol, Model, ecc.)
    # Vengono gestite come un dizionario di DataColumnV2
    aux_cols: Dict[str, DataColumnV2] = field(default_factory=dict)
    
    # Metadata di Alto Livello
    meta: Dict[str, Any] = field(default_factory=dict)
    
    # Session properties
    z_em: float = 0.0 # Emission redshift
    resol: float = 0.0 # Median resolution

    @classmethod
    def from_flat_data(cls, 
                       x_values: np.ndarray, x_unit: au.Unit,
                       y_values: np.ndarray, y_unit: au.Unit,
                       dy_values: np.ndarray, 
                       xmin_values: Optional[np.ndarray] = None,
                       xmax_values: Optional[np.ndarray] = None,
                       aux_data: Optional[Dict[str, tuple]] = None, # {'col_name': (values, unit, desc)}
                       meta: Optional[Dict[str, Any]] = None,
                       z_em: float = 0.0) -> 'SpectrumDataV2':
        """
        Factory method to create a SpectrumDataV2 from raw NumPy arrays.

        Parameters
        ----------
        x_values : np.ndarray
            Wavelength/Velocity data.
        x_unit : astropy.units.Unit
            Unit for the x-axis.
        y_values : np.ndarray
            Flux data.
        y_unit : astropy.units.Unit
            Unit for the y-axis (and dy).
        dy_values : np.ndarray
            Error data.
        xmin_values : np.ndarray, optional
            Lower bin edges. Defaults to ``x_values`` if None.
        xmax_values : np.ndarray, optional
            Upper bin edges. Defaults to ``x_values`` if None.
        aux_data : dict, optional
            Dictionary of auxiliary columns: ``{'name': (values, unit, description)}``.
        meta : dict, optional
            Metadata dictionary.
        z_em : float, optional
            Emission redshift. Defaults to ``0.0``.

        Returns
        -------
        SpectrumDataV2
            The initialized data object.
        """

        # Assicurare che xmin/xmax abbiano valori (es. copiando x se assenti)
        if xmin_values is None: xmin_values = x_values
        if xmax_values is None: xmax_values = x_values

        # Costruire le colonne principali
        x_col = DataColumnV2(x_values, x_unit)
        y_col = DataColumnV2(y_values, y_unit)
        dy_col = DataColumnV2(dy_values, y_unit) # dy usa la stessa unità di y
        xmin_col = DataColumnV2(xmin_values, x_unit)
        xmax_col = DataColumnV2(xmax_values, x_unit)
        
        # Costruire le colonne ausiliarie
        aux_cols = {}
        if aux_data:
            for name, (vals, unit, desc) in aux_data.items():
                 aux_cols[name] = DataColumnV2(vals, unit, desc)
        
        return cls(x=x_col, xmin=xmin_col, xmax=xmax_col, 
                   y=y_col, dy=dy_col, 
                   aux_cols=aux_cols, 
                   meta=meta if meta is not None else {},
                   z_em=z_em)
    
@dataclass(frozen=True)
class ComponentDataV2:
    """
    Immutable data structure for a single absorption component (Voigt profile).

    Attributes
    ----------
    id : int
        Internal integer ID (from V1 legacy logic).
    z : float
        Redshift of the component center.
    dz : float, optional
        Uncertainty on redshift.
    logN : float
        Logarithmic column density (cm^-2).
    dlogN : float, optional
        Uncertainty on logN.
    b : float
        Doppler broadening parameter (km/s).
    db : float, optional
        Uncertainty on b.
    btur : float
        Turbulent broadening contribution (km/s).
    func : str
        Profile function name (default: ``'voigt'``).
    series : str
        Atomic series name (e.g., ``'Ly_a'``, ``'CIV'``).
    chi2 : float, optional
        Reduced Chi-Squared of the fit group this component belongs to.
    resol : float, optional
        Resolution used for the fit.
    uuid : str
        A globally unique, stable identifier for V2 constraints and linking.
    """
    
    # ID is required
    id: int

    # Core Physical Parameters (z, logN, b)
    z: float
    dz: Optional[float] 
    logN: float
    dlogN: Optional[float]
    b: float
    db: Optional[float]
    
    # Turbulance Parameter (from V1)
    btur: float = 0.0
    dbtur: Optional[float] = None
    
    # Type and identification
    func: str = 'voigt'
    series: str = 'Ly_a'    

    # Fit Metadata
    chi2: Optional[float] = None # Reduced Chi-Squared of the fit group
    resol: Optional[float] = None # Resolution used for the fit

    # V2 IDENTIFIER: Stable, Global, Immutable
    uuid: str = field(default_factory=lambda: str(uuid.uuid4()))

@dataclass(frozen=True)
class ParameterConstraintV2:
    """
    Immutable structure defining a constraint for a single parameter.

    Used by the fitting engine to determine if a parameter is free, fixed,
    or linked to another component.

    Attributes
    ----------
    is_free : bool
        If True, the parameter varies during minimization.
    target_uuid : str, optional
        If set, this parameter is linked to the component with this UUID.
    expression : str, optional
        A mathematical string expression for linking (e.g. ``"p['{uuid}'].z"``).
    v1_target_id : int, optional
        Legacy ID of the target component (for debugging/reconstruction only).
    """
    # True if the parameter should be varied during minimization.
    is_free: bool
    
    # If not None, indicates a link. Must use UUID of the target component.
    target_uuid: Optional[str] = None
    
    # Stores the mathematical expression linking this parameter to the target UUID.
    expression: Optional[str] = None
    
    # The original V1 integer ID of the target (for debugging/reconstruction only)
    v1_target_id: Optional[int] = None

@dataclass(frozen=True)
class SystemListDataV2:
    """
    Immutable container for system components and constraints.

    This class replaces the mutable tables used in V1. It holds the list
    of components and the maps defining their relationships.

    Attributes
    ----------
    components : list of ComponentDataV2
        The list of absorption components.
    v1_header_constraints : dict
        Legacy constraints from V1 file headers.
    parsed_constraints : dict
        Map of (Component UUID, Parameter Name) -> ConstraintDataV2.
    v2_constraints_map : dict
        Nested map {UUID: {ParamName: Constraint}}.
    v1_id_to_uuid_map : dict
        Map of integer IDs to V2 UUIDs.
    meta : dict
        Metadata for the system list.
    """
    
    components: List[ComponentDataV2] = field(default_factory=list)
    v1_header_constraints: Dict[str, Any] = field(default_factory=dict)
    
    # Stores a dictionary mapping (Component UUID, Parameter Name) -> ConstraintDataV2
    parsed_constraints: Dict[Tuple[str, str], ParameterConstraintV2] = field(default_factory=dict)

    # V2-style constraints: {UUID: {ParamName: {...}}}
    v2_constraints_map: Dict[str, Dict[str, ParameterConstraintV2]] = field(default_factory=dict)
    
    # Add the missing field
    v1_id_to_uuid_map: Dict[int, str] = field(default_factory=dict)


    # Placeholder for complex V1 structures (mutable state reference, if needed)
    v1_models_t: Any = None 
    meta: Dict[str, Any] = field(default_factory=dict)

@dataclass
class LogEntryV2:
    """
    A single entry in the V2 analysis log.
    
    Corresponds 1-to-1 with a state change in the Session history.
    """
    recipe_name: str
    params: Dict[str, Any]
    # We can add more data later, like execution time, success, etc.

@dataclass
class HistoryLogV2:
    """
    Manages the linear history of operations (Log V2).

    This class maintains a list of :class:`LogEntryV2` objects matching the
    stack of immutable states in the session manager. It supports truncation
    to handle branching history (undo/redo logic).

    Attributes
    ----------
    entries : list of LogEntryV2
        The list of operations performed.
    current_index : int
        Pointer to the currently active log entry (state).
    """
    def __init__(self):
        self.entries: List[LogEntryV2] = []
        # current_index points to the *last valid entry*
        self.current_index: int = -1 

    @property
    def is_v2_log(self) -> bool:
        """Return True, indicating this is a V2-native log structure."""
        return True

    def add_entry(self, recipe_name: str, params: Dict[str, Any]):
        """
        Add a new log entry, truncating any 'future' history (from Undo actions).

        Parameters
        ----------
        recipe_name : str
            The name of the recipe executed.
        params : dict
            The parameters used for the recipe.
        """
        # 1. Truncate stale future history
        if self.current_index < len(self.entries) - 1:
            logging.debug(f"Truncating V2 log from index {self.current_index + 1}")
            del self.entries[self.current_index + 1:]
        
        # 2. Add the new entry
        entry = LogEntryV2(recipe_name=recipe_name, params=params)
        self.entries.append(entry)
        
        # 3. Update the index to point to the new entry
        self.current_index = len(self.entries) - 1
        logging.debug(f"Added log entry: {recipe_name}. Index is now {self.current_index}.")

    def undo(self) -> bool:
        """
        Add a new log entry, truncating any 'future' history (from Undo actions).

        Parameters
        ----------
        recipe_name : str
            The name of the recipe executed.
        params : dict
            The parameters used for the recipe.
        """
        if self.current_index > -1:
            self.current_index -= 1
            logging.debug(f"Log index undone to {self.current_index}")
            return True
        return False

    def redo(self) -> bool:
        """
        Move the history pointer forward one step.

        Returns
        -------
        bool
            True if redo was successful, False if already at the latest entry.
        """
        if self.current_index < len(self.entries) - 1:
            self.current_index += 1
            logging.debug(f"Log index redone to {self.current_index}")
            return True
        return False

# --- V1 Compatibility Wrapper ---

class V1LogArtifact:
    """
    A read-only wrapper for legacy V1 logs.

    Used to display the history of a session loaded from a V1 archive
    without converting it to the full V2 structure.

    Parameters
    ----------
    v1_log_json : dict
        The raw JSON dictionary from the V1 archive.
    """
    def __init__(self, v1_log_json: dict):
        self.v1_json = v1_log_json
    
    @property
    def is_v2_log(self) -> bool:
        """Return False, indicating this is a legacy V1 log wrapper."""
        return False