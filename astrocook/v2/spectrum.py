
import astropy.units as au
from copy import deepcopy
import numpy as np
from typing import Dict, Any, Optional, Union

from .spectrum_operations import convert_axis_velocity, convert_x_axis, convert_y_axis, rebin_spectrum
from .structures import SpectrumDataV2, DataColumnV2
from ..v1.frame import Frame as FrameV1
from ..v1.spectrum import Spectrum as SpectrumV1

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
    def cont(self) -> Optional[au.Quantity]:
        """
        Returns the 'cont' auxiliary column as an Astropy Quantity,
        or None if the column does not exist.
        """
        cont_col = self._data.aux_cols.get('cont')
        if cont_col is not None:
            # Combine values and unit into a Quantity
            return au.Quantity(cont_col.values, cont_col.unit)
        return None

    @property
    def model(self) -> Optional[au.Quantity]:
        """
        Returns the 'model' auxiliary column as an Astropy Quantity,
        or None if the column does not exist.
        """
        model_col = self._data.aux_cols.get('model')
        if model_col is not None:
            # Combine values and unit into a Quantity
            return au.Quantity(model_col.values, model_col.unit)
        return None

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

    def to_v1_spectrum(self) -> SpectrumV1:
        """
        Converts the immutable SpectrumV2 (API + Data Core) back into a mutable 
        SpectrumV1 object for saving.
        """
        # NOTE: This is complex because we need to extract the data and auxiliary columns
        # and feed them to the V1 SpectrumV1 constructor (which calls FrameV1.__init__).
        
        data_core = self._data
        
        # 1. Prepare core data lists/arrays
        x_list = data_core.x.values
        xmin_list = data_core.xmin.values
        xmax_list = data_core.xmax.values
        y_list = data_core.y.values
        dy_list = data_core.dy.values
        
        # 2. Instantiate the V1 Spectrum object
        v1_spec = SpectrumV1(
            x_list, xmin_list, xmax_list, 
            y_list, dy_list, 
            xunit=data_core.x.unit, 
            yunit=data_core.y.unit, 
            meta=data_core.meta
        )

        # NOTE: V1 Spectrum objects often store meta in an internal dict/Astropy Header.
        # We assume the V1 constructor places 'meta' into an accessible attribute like v1_spec._meta
        
        # Ensure the V1 object's metadata is explicitly set to include 'ORIGIN'
        if not hasattr(v1_spec, 'meta'):
             # If V1 object doesn't have 'meta', we cannot easily inject it.
             # Assume V1 objects have a standard attribute for metadata.
             # We rely on the V1 object having a setter or an accessible dictionary.
             
             # If using a direct attribute:
             # v1_spec._meta['ORIGIN'] = 'Astrocook V2' 
             
             # If using the provided dictionary directly (best approach):
             v1_spec.meta['ORIGIN'] = 'Astrocook V2' 

        else:
             # If V1 object exposes metadata via a property (safer):
             v1_spec.meta['ORIGIN'] = 'Astrocook V2'

             
        # 3. Add auxiliary columns (e.g., cont, model, etc.)
        for name, data_col in data_core.aux_cols.items():
            # Directly add the column to the V1 object's internal table (_t)
            v1_spec._t[name] = data_col.quantity # Quantity ensures unit is preserved
            
        return v1_spec
    
    
    # --- Methods implementing the API ---

    def arithmetics(self, col_target: str, col_left: str, op: str, col_right: str) -> 'SpectrumV2':
        """
        API: Performs column arithmetics and returns a NEW SpectrumV2 instance.
        """
        
        # 1. Get all data columns via the adapter
        all_cols = self.t._data_dict # Access the raw dict from the adapter
        
        # 2. Get data arrays
        if col_left not in all_cols: raise ValueError(f"Column '{col_left}' not found.")
        if col_right not in all_cols: raise ValueError(f"Column '{col_right}' not found.")
        # Note: col_target does not need to exist yet
            
        left_q = all_cols[col_left]
        right_q = all_cols[col_right]

        # 3. Perform operation
        new_quantity = None
        if op == '+': new_quantity = left_q + right_q
        elif op == '-': new_quantity = left_q - right_q
        elif op == '*': new_quantity = left_q * right_q
        elif op == '/': 
            new_values = np.divide(left_q.value, right_q.value, 
                                   out=np.full_like(left_q.value, np.nan), 
                                   where=right_q.value!=0)
            new_quantity = au.Quantity(new_values, left_q.unit / right_q.unit)
        else:
            raise ValueError(f"Invalid operator '{op}'")
            
        # 4. Create new SpectrumDataV2 using the immutable pattern
        # Start with all original data
        core_cols = {
            'x': self._data.x,
            'xmin': self._data.xmin,
            'xmax': self._data.xmax,
            'y': self._data.y,
            'dy': self._data.dy
        }
        aux_cols = deepcopy(self._data.aux_cols) # Copy the aux_cols dict

        # 5. Create the new DataColumnV2
        new_data_col = DataColumnV2(new_quantity.value, new_quantity.unit)

        # 6. Overwrite the target column in the *correct dictionary*
        if col_target in core_cols:
            core_cols[col_target] = new_data_col
        elif col_target in aux_cols:
            # This allows overwriting an existing aux col
            aux_cols[col_target] = new_data_col
        else:
            # This was a new column, add it to aux_cols
            logging.debug(f"Arithmetics: creating new aux_col '{col_target}'")
            aux_cols[col_target] = new_data_col

        # 7. Create the new data core
        new_data_core = SpectrumDataV2(
            x=core_cols['x'],
            xmin=core_cols['xmin'],
            xmax=core_cols['xmax'],
            y=core_cols['y'],
            dy=core_cols['dy'],
            aux_cols=aux_cols,
            meta=deepcopy(self._data.meta),
            rf_z=self._data.rf_z
        )

        # 8. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Arithmetics: {col_target} = {col_left} {op} {col_right}"]
        return SpectrumV2(data=new_data_core, history=new_history)

    def mask(self, col_target: str, col_cond: str, op: str, value: float) -> 'SpectrumV2':
        """
        API: Masks a column based on a condition and returns a NEW SpectrumV2 instance.
        """
        
        # 1. Get all data columns
        all_cols = self.t._data_dict
        
        # 2. Get data arrays
        if col_target not in all_cols: raise ValueError(f"Target column '{col_target}' not found.")
        if col_cond not in all_cols: raise ValueError(f"Condition column '{col_cond}' not found.")
            
        target_q = all_cols[col_target]
        cond_q = all_cols[col_cond]
        
        # 3. Find indices to mask
        mask_indices = None
        cond_vals = cond_q.value # Get raw numpy array
        
        if op == '<=': mask_indices = np.where(cond_vals <= value)
        elif op == '<': mask_indices = np.where(cond_vals < value)
        elif op == '==': mask_indices = np.where(cond_vals == value)
        elif op == '>': mask_indices = np.where(cond_vals > value)
        elif op == '>=': mask_indices = np.where(cond_vals >= value)
        elif op == '!=': mask_indices = np.where(cond_vals != value)
        else:
            raise ValueError(f"Invalid operator '{op}'")
            
        # 4. Create new data array
        new_values = target_q.value.copy()
        new_values[mask_indices] = np.nan
        new_quantity = au.Quantity(new_values, target_q.unit)
        new_data_col = DataColumnV2(new_quantity.value, new_quantity.unit)

        # 5. Create new SpectrumDataV2 using the immutable pattern
        # Start with all original data
        core_cols = {
            'x': self._data.x,
            'xmin': self._data.xmin,
            'xmax': self._data.xmax,
            'y': self._data.y,
            'dy': self._data.dy
        }
        aux_cols = deepcopy(self._data.aux_cols) # Copy the aux_cols dict

        # 6. Overwrite the target column in the *correct dictionary*
        if col_target in core_cols:
            core_cols[col_target] = new_data_col
        elif col_target in aux_cols:
            aux_cols[col_target] = new_data_col
        else:
            raise ValueError(f"Target column '{col_target}' is not a core or aux column.")

        # 7. Create the new data core
        new_data_core = SpectrumDataV2(
            x=core_cols['x'],
            xmin=core_cols['xmin'],
            xmax=core_cols['xmax'],
            y=core_cols['y'],
            dy=core_cols['dy'],
            aux_cols=aux_cols,
            meta=deepcopy(self._data.meta),
            rf_z=self._data.rf_z
        )

        # 8. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Mask: {col_target} where {col_cond} {op} {value}"]
        return SpectrumV2(data=new_data_core, history=new_history)

    def split(self, col_check: str, min_val: float, max_val: float) -> 'SpectrumV2':
        """
        API: Extracts a slice of the spectrum and returns a NEW SpectrumV2 instance.
        """
        
        # 1. Get all data columns
        all_cols = self.t._data_dict
        if col_check not in all_cols: raise ValueError(f"Check column '{col_check}' not found.")
        
        check_q = all_cols[col_check]
        
        # 2. Find indices for the slice
        # This assumes the check column is sorted (like 'x')
        idx_start = np.searchsorted(check_q.value, min_val, side='left')
        idx_end = np.searchsorted(check_q.value, max_val, side='right')
        
        if idx_start >= idx_end:
            raise ValueError("Split range resulted in an empty spectrum.")
            
        data_slice = slice(idx_start, idx_end)

        # 3. Create a new data core by slicing all columns
        # We can use the TableAdapterV2's slicing capability
        
        # Create a new dictionary to hold the sliced data
        new_data_cols_dict = {}
        
        # Iterate through all columns and apply the slice
        for col_name, quantity in all_cols.items():
            new_values = quantity.value[data_slice]
            new_unit = quantity.unit
            new_data_cols_dict[col_name] = DataColumnV2(new_values, new_unit)

        # 4. Build new SpectrumDataV2
        new_data_core = SpectrumDataV2(
            x=new_data_cols_dict.pop('x'),
            xmin=new_data_cols_dict.pop('xmin'),
            xmax=new_data_cols_dict.pop('xmax'),
            y=new_data_cols_dict.pop('y'),
            dy=new_data_cols_dict.pop('dy'),
            aux_cols=new_data_cols_dict, # Remaining items are aux_cols
            meta=deepcopy(self._data.meta), 
            rf_z=self._data.rf_z 
        )
        
        # 5. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Split: {col_check} from {min_val} to {max_val}"]
        return SpectrumV2(data=new_data_core, history=new_history)
    

    def rebin(self, xstart: Optional[au.Quantity], xend: Optional[au.Quantity], 
              dx: au.Quantity, kappa: Optional[float], 
              filling: float) -> 'SpectrumV2':
        """
        API: Rebins the spectrum and returns a NEW SpectrumV2 instance.
        """
        
        original_x_unit = self._data.x.unit 
        zem = self._data.rf_z # Use the current rest-frame Z as the conversion zero point

        # Determine target unit for rebinning (km/s if the step is in km/s, otherwise original unit)
        if dx.unit.is_equivalent(au.km/au.s):
            target_calc_unit = au.km/au.s
        else:
            target_calc_unit = original_x_unit

        # --- Unit Conversion to Calculation Space ---
        x_calc = convert_axis_velocity(self._data.x.quantity, zem, target_calc_unit)
        xmin_calc = convert_axis_velocity(self._data.xmin.quantity, zem, target_calc_unit)
        xmax_calc = convert_axis_velocity(self._data.xmax.quantity, zem, target_calc_unit)

        # 1. Execute pure rebinning logic
        x_new, xmin_new, xmax_new, y_new, dy_new = rebin_spectrum(
            x_calc, xmin_calc, xmax_calc, # Pass the converted X-arrays
            self._data, 
            xstart, xend, dx, target_calc_unit,
            kappa, filling
        )
        
        x_final = convert_axis_velocity(x_new, zem, original_x_unit)
        xmin_final = convert_axis_velocity(xmin_new, zem, original_x_unit)
        xmax_final = convert_axis_velocity(xmax_new, zem, original_x_unit)

        # 2. Build new SpectrumDataV2
        new_data = SpectrumDataV2(
            x=DataColumnV2(x_final.value, x_final.unit),
            xmin=DataColumnV2(xmin_final.value, xmin_final.unit),
            xmax=DataColumnV2(xmax_final.value, xmax_final.unit),
            y=DataColumnV2(y_new.value, y_new.unit),
            dy=DataColumnV2(dy_new.value, dy_new.unit),
            #aux_cols=self._data.aux_cols, # Auxiliary columns are NOT handled in core rebin
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