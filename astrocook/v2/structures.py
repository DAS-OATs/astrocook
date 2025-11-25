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
    This will be serialized to a JSON/YAML file within the .acs archive.
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
    """Contenitore per una singola colonna di dati con unità."""
    values: np.ndarray
    unit: au.Unit
    description: str = ""
    
    @property
    def quantity(self) -> au.Quantity:
        return au.Quantity(self.values, self.unit)

@dataclass(frozen=True)
class SpectrumDataV2:
    """Contenitore immutabile per tutti i dati di uno Spettro (Core Logico)."""
    
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
    z_rf: float = 0.0 # Redshift del rest frame (was rf_z)
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
                       z_rf: float = 0.0,
                       z_em: float = 0.0) -> 'SpectrumDataV2':

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
                   z_rf=z_rf,
                   z_em=z_em)
    
@dataclass(frozen=True)
class ComponentDataV2:
    """Immutable data structure for a single absorption component (Voigt profile)."""
    
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

    # V2 IDENTIFIER: Stable, Global, Immutable
    uuid: str = field(init=False, default_factory=lambda: str(uuid.uuid4()))

    def __post_init__(self):
        # NOTE: The dataclass decorator handles the default_factory for frozen objects,
        # but we use object.__setattr__ here for safety/explicit control if other fields 
        # were using post_init (and to ensure the UUID is set if not via default_factory)
        
        if not hasattr(self, 'uuid'):
             object.__setattr__(self, 'uuid', str(uuid.uuid4()))

@dataclass(frozen=True)
class ParameterConstraintV2:
    """
    Immutable data structure defining the constraint for a single parameter.
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
    """Immutable container for system components (list of ComponentDataV2)."""
    
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
    """A V2-native log entry. 1-to-1 with a state change."""
    recipe_name: str
    params: Dict[str, Any]
    # We can add more data later, like execution time, success, etc.

@dataclass
class HistoryLogV2:
    """
    Manages the V2-native log. This is a 1-to-1 list of states
    and supports truncation for Option 2.
    """
    def __init__(self):
        self.entries: List[LogEntryV2] = []
        # current_index points to the *last valid entry*
        self.current_index: int = -1 

    @property
    def is_v2_log(self) -> bool:
        return True

    def add_entry(self, recipe_name: str, params: Dict[str, Any]):
        """
        Adds a new entry, truncating any "future" (undone) logs.
        This implements Option 2's truncation.
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
        """Moves the index back one step. Does not delete."""
        if self.current_index > -1:
            self.current_index -= 1
            logging.debug(f"Log index undone to {self.current_index}")
            return True
        return False

    def redo(self) -> bool:
        """Moves the index forward one step."""
        if self.current_index < len(self.entries) - 1:
            self.current_index += 1
            logging.debug(f"Log index redone to {self.current_index}")
            return True
        return False

# --- V1 Compatibility Wrapper ---

class V1LogArtifact:
    """
    A read-only wrapper for a V1 log string.
    It mimics just enough of the HistoryLogV2 interface
    for the LogViewerDialog.
    """
    def __init__(self, v1_log_json: dict):
        self.v1_json = v1_log_json
    
    @property
    def is_v2_log(self) -> bool:
        return False