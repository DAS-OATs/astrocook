from .spectrum_operations import convert_x_axis, convert_y_axis, rebin_spectrum
from .structures import SpectrumDataV2, DataColumnV2

import astropy.units as au
from copy import deepcopy
import numpy as np
from typing import Self, Dict, Any, Optional, Union # Usiamo Self per i metodi che restituiscono una nuova istanza

class TableAdapterV2(object):
    """Adapter to wrap the V2 .t dictionary and expose Astropy Table-like properties."""
    def __init__(self, data_dict):
        self._data_dict = data_dict
        self.colnames = list(data_dict.keys())
        
    def __getitem__(self, key):
        return self._data_dict[key]
        
    # Crucially, add a method to simulate row access (needed for row iteration)
    def __len__(self):
        # Assumes 'x' is present and holds the length
        return len(self._data_dict['x']) 

    # Implement a simple iterator for row access (the "r" in "for r in t")
    def __iter__(self):
        self._current_row = 0
        return self

    def __next__(self):
        if self._current_row < len(self):
            row_data = {col: self._data_dict[col][self._current_row] for col in self.colnames}
            self._current_row += 1
            return row_data
        raise StopIteration
    
    def __getitem__(self, key: Union[str, slice]) -> Union[au.Quantity, 'TableAdapterV2']:
        """
        Handles access via column name (string) OR slicing (slice object).
        """
        if isinstance(key, str):
            # Case 1: Column name lookup (original working logic)
            return self._data_dict[key]
        
        elif isinstance(key, slice):
            # Case 2: Row slicing (NEW LOGIC)
            
            # Create a new dictionary to hold the sliced data
            new_data_dict = {}
            
            # Iterate through all columns and apply the slice to their arrays
            for col_name, quantity in self._data_dict.items():
                
                # Check if the object is a Quantity (most V2 columns are)
                if hasattr(quantity, 'value'):
                    # Apply slice to the underlying NumPy array for efficiency
                    new_values = quantity.value[key]
                    new_unit = quantity.unit
                    # Recreate the Quantity for consistency (or just use the sliced array if safe)
                    new_data_dict[col_name] = au.Quantity(new_values, new_unit)
                else:
                    # Handle non-Quantity items (like metadata lists if any)
                    new_data_dict[col_name] = quantity[key]
                    
            # Return a NEW TableAdapterV2 instance containing the sliced data
            return TableAdapterV2(new_data_dict)
            
        else:
            raise KeyError(f"Invalid key type: {type(key)}. Expected string or slice.")

class SpectrumV2:
    """
    API pubblica e contenitore composito per l'analisi spettrale.
    Tutti i metodi operativi (che trasformano i dati) restituiscono una NUOVA istanza SpectrumV2.
    """
    
    def __init__(self, data: SpectrumDataV2, history: list = None):
        self._data = data
        self.history = history if history is not None else []

    # Accessori (API pubblica: NESSUN setter per i dati)
    @property
    def x(self) -> au.Quantity:
        return self._data.x.quantity

    @property
    def xmin(self) -> au.Quantity:
        return self._data.xmin.quantity

    @property
    def xmax(self) -> au.Quantity:
        return self._data.xmax.quantity

    @property
    def y(self) -> au.Quantity:
        return self._data.y.quantity

    @property
    def dy(self) -> au.Quantity:
        return self._data.dy.quantity

    @property
    def meta(self) -> Dict[str, Any]:
        return deepcopy(self._data.meta) # Restituisce una copia per prevenire modifiche esterne

    @property
    def t(self) -> TableAdapterV2:
        """
        Adapter: Replaces the raw table/dictionary access.
        Returns a TableAdapterV2 instance to simulate an Astropy Table (V1)
        for compatibility with GUI/V1 logic (like gui_table.py).
        """
        # 1. Build the raw data dictionary (Integration of existing logic)
        output = {
            'x': self.x,
            'xmin': self.xmin,
            'xmax': self.xmax,
            'y': self.y,
            'dy': self.dy,
        }
        
        # 2. Add auxiliary columns
        for name, col_data in self._data.aux_cols.items():
            output[name] = col_data.quantity
            
        # 3. Return the Adapter Instance (Replacement)
        return TableAdapterV2(output)

    @property
    def _rfz(self) -> float:
        """
        Adapter GETTER per l'attributo V1 '_rfz'.
        Restituisce l'attributo immutabile rf_z dal contenitore dati.
        """
        return self._data.rf_z

    @_rfz.setter
    def _rfz(self, val: float):
        """
        Adapter SETTER per l'attributo V1 '_rfz'.
        Ignora silenziosamente la modifica, poiché SpectrumV2 è immutabile.
        """
        # Se il codice V1 chiama self.spec._rfz = 1.5, l'oggetto V2 non viene modificato.
        # L'unica vera modifica dovrebbe avvenire tramite un metodo che restituisce
        # una NUOVA istanza SpectrumV2 (come shift_rest_frame che avevamo abbozzato).
        pass # Ignora la modifica in-place

    @property
    def _xunit(self) -> au.Unit:
        """
        Adapter GETTER per l'attributo V1 '_xunit'.
        Restituisce l'unità X corrente del contenitore dati.
        """
        return self._data.x.unit
    
    @_xunit.setter
    def _xunit(self, val: au.Unit):
        """
        Adapter SETTER per l'attributo V1 '_xunit'.
        Ignora silenziosamente la modifica (immutabilità V2).
        """
        pass # L'unità X dovrebbe essere aggiornata solo creando un nuovo SpectrumV2
    
    @property
    def _xunit_old(self) -> au.Unit:
        """
        Adapter GETTER per l'attributo V1 '_xunit_old'.
        Questo attributo non esiste nel dato V2, quindi lo emuliamo.
        """
        # Usiamo un fallback nel meta o semplicemente l'unità X stessa.
        # Se il codice V1 ha bisogno di una versione "vecchia" persistente,
        # dovrai modificarne la logica, ma per ora lo emuliamo:
        
        # Tentiamo di recuperare l'unità iniziale dai meta se disponibile, altrimenti l'unità corrente.
        # Se il codice V1 tenta di assegnarlo, il setter vuoto preverrà il crash.
        return self._data.meta.get('initial_x_unit', self._data.x.unit)
    
    @_xunit_old.setter
    def _xunit_old(self, val: au.Unit):
        """
        Adapter SETTER per l'attributo V1 '_xunit_old'.
        Ignora silenziosamente la modifica.
        """
        pass

    @property
    def _yunit(self) -> au.Unit:
        """
        Adapter GETTER per l'attributo V1 '_yunit'.
        Restituisce l'unità Y corrente dal contenitore dati.
        """
        return self._data.y.unit # Restituisce l'unità della colonna y
    
    @_yunit.setter
    def _yunit(self, val: au.Unit):
        """
        Adapter SETTER per l'attributo V1 '_yunit'.
        Ignora silenziosamente la modifica (immutabilità V2).
        """
        pass 

    def _safe(self, quantity: au.Quantity) -> au.Quantity:
        """
        Adapter per il metodo V1 _safe:
        Restituisce una nuova Quantity contenente solo i valori non-NaN.
        (Non modifica lo stato interno come faceva la V1)
        """
        if quantity is None:
            return None
            
        where_safe = ~np.isnan(quantity.value)
        return quantity[where_safe]

    def has_aux_column(self, name: str) -> bool:
        """Checks if an auxiliary column exists."""
        return name in self._data.aux_cols
    
    def get_column(self, name: str) -> Optional[au.Quantity]:
        """Accede a una colonna ausiliaria per la GUI (es. 'cont')."""
        col = self._data.aux_cols.get(name)
        return col.quantity if col else None
    
    # Methods implementing the API

    def rebin(self, xstart: Optional[au.Quantity], xend: Optional[au.Quantity], 
              dx: au.Quantity, kappa: Optional[float], 
              filling: float) -> 'SpectrumV2':
        """
        API: Rebins the spectrum and returns a NEW SpectrumV2 instance.
        """
        
        # 1. Execute pure rebinning logic
        x_new, xmin_new, xmax_new, y_new, dy_new = rebin_spectrum(
            self._data, xstart, xend, dx, kappa, filling
        )
        
        # 2. Build new SpectrumDataV2
        new_data = SpectrumDataV2(
            x=DataColumnV2(x_new.value, x_new.unit),
            xmin=DataColumnV2(xmin_new.value, x_new.unit),
            xmax=DataColumnV2(xmax_new.value, x_new.unit),
            y=DataColumnV2(y_new.value, y_new.unit),
            dy=DataColumnV2(dy_new.value, dy_new.unit),
            aux_cols=self._data.aux_cols, # Auxiliary columns are NOT handled in core rebin
            meta=self._data.meta, rf_z=self._data.rf_z 
        )
        
        # 3. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Rebinned spectrum (dx={dx})"]
        return SpectrumV2(data=new_data, history=new_history)

    def x_convert(self, zem: float, xunit: au.Unit) -> 'SpectrumV2':
        """
        API: Converts the X-axis units and returns a NEW SpectrumV2 instance.
        (Replaces V1's mutable _x_convert)
        """
        
        # 1. Execute pure conversion logic
        x_new, xmin_new, xmax_new = convert_x_axis(self._data, zem, xunit)
        
        # 2. Build new SpectrumDataV2, preserving the core vertical data
        new_data = SpectrumDataV2(
            x=DataColumnV2(x_new.value, xunit),
            xmin=DataColumnV2(xmin_new.value, xunit),
            xmax=DataColumnV2(xmax_new.value, xunit),
            y=self._data.y, 
            dy=self._data.dy,
            aux_cols=self._data.aux_cols,
            meta=self._data.meta,
            rf_z=zem # Update rest frame Z in the new object
        )
        
        # 3. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"X-axis converted to {xunit} at zem={zem}"]
        return SpectrumV2(data=new_data, history=new_history)

    def y_convert(self, yunit: au.Unit) -> 'SpectrumV2':
        """
        API: Converts the Y-axis units and returns a NEW SpectrumV2 instance.
        (Replaces V1's mutable _y_convert)
        """
        
        # 1. Execute pure conversion logic
        new_y, new_dy, new_aux_cols = convert_y_axis(self._data, yunit)
        
        # 2. Build new SpectrumDataV2, preserving the X-axis and metadata
        new_data = SpectrumDataV2(
            x=self._data.x, 
            xmin=self._data.xmin,
            xmax=self._data.xmax,
            y=new_y, # Updated
            dy=new_dy, # Updated
            aux_cols=new_aux_cols, # Auxiliary columns should also be converted if unit-dependent
            meta=self._data.meta,
            rf_z=self._data.rf_z 
        )
        
        # 3. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Y-axis converted to {yunit}"]
        return SpectrumV2(data=new_data, history=new_history)