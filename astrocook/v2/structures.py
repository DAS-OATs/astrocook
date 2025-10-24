from dataclasses import dataclass, field
from astropy import units as au
from astropy.table import Table
import numpy as np
from typing import Any, Dict, List, Optional, Tuple
import uuid

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
    rf_z: float = 0.0 # Redshift del rest frame

    @classmethod
    def from_flat_data(cls, 
                       x_values: np.ndarray, x_unit: au.Unit,
                       y_values: np.ndarray, y_unit: au.Unit,
                       dy_values: np.ndarray, 
                       xmin_values: Optional[np.ndarray] = None,
                       xmax_values: Optional[np.ndarray] = None,
                       aux_data: Optional[Dict[str, tuple]] = None, # {'col_name': (values, unit, desc)}
                       meta: Optional[Dict[str, Any]] = None,
                       rf_z: float = 0.0) -> 'SpectrumDataV2':

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
                   rf_z=rf_z)
    
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
class SystemListDataV2:
    """Immutable container for system components (list of ComponentDataV2)."""
    
    components: List[ComponentDataV2] = field(default_factory=list)
    v1_header_constraints: Dict[str, Any] = field(default_factory=dict)
    parsed_constraints: Dict[Tuple[int, str], Dict[str, Any]] = field(default_factory=dict)

    # Placeholder for complex V1 structures (mutable state reference, if needed)
    v1_models_t: Any = None 
    meta: Dict[str, Any] = field(default_factory=dict)