
import astropy.units as au
from copy import deepcopy
import dataclasses
import json
import logging
import numexpr as ne
import numpy as np
from typing import Any, Dict, List, Optional, Union

from .atomic_data import get_multiplet_velocity_lags, is_hydrogen_line, STANDARD_MULTIPLETS, xem_d
from .photometry import generate_calibration_curve
from .spectrum_operations import (
    compute_identification_signal, convert_axis_velocity, 
    convert_x_axis, convert_y_axis, detect_regions, find_absorbed_regions, 
    find_signal_peaks, fit_continuum_interp, fit_powerlaw_to_regions, smooth_spectrum, rebin_spectrum, running_std,
)
from .structures import SpectrumDataV2, DataColumnV2
from ..v1.frame import Frame as FrameV1
from ..v1.spectrum import Spectrum as SpectrumV1

_NE_GLOBAL_DICT = {
    'min': np.nanmin,
    'max': np.nanmax,
    'median': np.nanmedian,
    'mean': np.nanmean,
    'std': np.nanstd,
    # numexpr *does* have: sum, prod, mean, std, min, max, log, log10, etc.
}
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
    def _z_rf(self) -> float:
        """
        Adapter GETTER per l'attributo V1 '_rfz'.
        Restituisce l'attributo immutabile z_rf.
        """
        return self._data.z_rf

    @_z_rf.setter
    def _z_rf(self, val: float):
        """
        Adapter SETTER per l'attributo V1 '_rfz'.
        Ignora silenziosamente la modifica.
        """
        pass 
        
    @property
    def _z_em(self) -> float:
        """
        Adapter GETTER per l'attributo V1 '_zem'.
        Restituisce l'attributo immutabile z_em.
        """
        return self._data.z_em

    @_z_em.setter
    def _z_em(self, val: float):
        """
        Adapter SETTER per l'attributo V1 '_zem'.
        Ignora silenziosamente la modifica.
        """
        pass

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
    
    def with_properties(self, z_em: float) -> 'SpectrumV2':
        """
        API: Returns a NEW SpectrumV2 instance with updated properties.
        """
        new_data = dataclasses.replace(self._data, z_em=z_em)
        new_history = self.history + [f"Set properties (z_em={z_em})"]
        return SpectrumV2(data=new_data, history=new_history)

    def _get_ne_contexts(self) -> (Dict[str, au.Quantity], Dict[str, np.ndarray]):
        """
        Helper to create the contexts for numexpr.
        Returns:
            all_cols (dict): Map of name -> Astropy Quantity (for unit logic)
            local_dict (dict): Map of name -> numpy array (for numexpr variables)
        """
        all_cols = self.t._data_dict
        # 1. Pass only the raw numpy arrays to local_dict
        local_dict = {k: v.value for k, v in all_cols.items()}
        
        # 2. The functions are now passed via _NE_GLOBAL_DICT
        
        return all_cols, local_dict
    
    def _renormalize_model(self, 
                           old_cont_vals: np.ndarray, 
                           new_cont_vals: np.ndarray, 
                           old_model_col: DataColumnV2) -> Optional[DataColumnV2]:
        """
        Private helper to re-normalize an old model to a new continuum
        by applying the *old* normalized profile to the *new* continuum.
        
        Used by fit_continuum and smooth_column.
        """
        if old_cont_vals is None or new_cont_vals is None or old_model_col is None:
            return None
            
        logging.info("Re-normalizing 'model' column to new 'cont'...")
        try:
            # 1. Get old normalized profile
            normalized_model = np.divide(
                old_model_col.values,
                old_cont_vals,
                out=np.ones_like(old_model_col.values),
                where=old_cont_vals != 0
            )
            # 2. Recalculate model against *new* continuum
            new_model_values = new_cont_vals * normalized_model
            
            return DataColumnV2(
                values=new_model_values,
                unit=old_model_col.unit, # Use old model's unit
                description=old_model_col.description
            )
        except Exception as e:
            logging.warning(f"Failed to renormalize model column: {e}")
            return None
    
    # --- Methods implementing the API ---

    def apply_expression(self, target_col: str, expression: str, 
                         extra_vars: Dict[str, np.ndarray] = None) -> 'SpectrumV2': # <<< MODIFIED
        """
        API: Applies a numerical expression and returns a NEW SpectrumV2 instance.
        """
        
        # 1. Get column data and the numexpr variable dict
        all_cols, local_dict = self._get_ne_contexts()
        
        if extra_vars:
            local_dict.update(extra_vars)

        try:
            # 2. Evaluate the expression
            new_values = ne.evaluate(expression, 
                                     global_dict=_NE_GLOBAL_DICT,  # <<< *** PASS GLOBALS ***
                                     local_dict=local_dict)       # <<< PASS LOCALS
        except Exception as e:
            logging.error(f"Numexpr evaluation failed for expression: '{expression}'")
            raise e # Re-raise for the recipe to catch

        # 3. Determine the unit for the new column
        target_unit = au.dimensionless_unscaled
        
        # If the result is a scalar (e.g., from median(y)),
        # we must broadcast it to the shape of the spectrum
        if not hasattr(new_values, "shape") or new_values.shape == ():
            logging.debug("Expression result is a scalar, broadcasting to array.")
            # Use shape of 'x' as the reference
            new_values = np.full(self._data.x.values.shape, new_values)

        if target_col in all_cols:
            target_unit = all_cols[target_col].unit
            logging.debug(f"Expression result will inherit unit '{target_unit}' from target '{target_col}'")
        else:
            logging.debug(f"Expression result for new column '{target_col}' will be dimensionless.")

        # 4. Create new SpectrumDataV2 using the immutable pattern
        core_cols = {
            'x': self._data.x, 'xmin': self._data.xmin, 'xmax': self._data.xmax,
            'y': self._data.y, 'dy': self._data.dy
        }
        aux_cols = deepcopy(self._data.aux_cols)

        # 5. Create the new DataColumnV2
        new_data_col = DataColumnV2(new_values, target_unit)

        # 6. Overwrite or add the target column
        if target_col in core_cols:
            core_cols[target_col] = new_data_col
        else:
            aux_cols[target_col] = new_data_col

        # 7. Create the new data core
        new_data_core = SpectrumDataV2(
            x=core_cols['x'], xmin=core_cols['xmin'], xmax=core_cols['xmax'],
            y=core_cols['y'], dy=core_cols['dy'],
            aux_cols=aux_cols,
            meta=deepcopy(self._data.meta),
            z_rf=self._data.z_rf,  # Propagate
            z_em=self._data.z_em   # Propagate
        )

        # 8. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Apply: {target_col} = {expression}"]
        return SpectrumV2(data=new_data_core, history=new_history)

    def mask_expression(self, target_col: str, expression: str, 
                        extra_vars: Dict[str, np.ndarray] = None) -> 'SpectrumV2': # <<< MODIFIED
        """
        API: Masks a column based on a boolean expression and returns a NEW SpectrumV2 instance.
        """
        
        # 1. Get column data and the numexpr variable dict
        all_cols, local_dict = self._get_ne_contexts()
        
        if extra_vars:
            local_dict.update(extra_vars)

        # 2. Check that target column exists
        if target_col not in all_cols: 
            raise ValueError(f"Target column '{target_col}' not found.")
            
        target_q = all_cols[target_col]

        try:
            # 3. Evaluate the boolean expression to get the mask
            mask_indices = ne.evaluate(expression, 
                                       global_dict=_NE_GLOBAL_DICT,  # <<< *** PASS GLOBALS ***
                                       local_dict=local_dict)       # <<< PASS LOCALS
            if mask_indices.dtype != bool:
                raise TypeError("Mask expression did not return a boolean array.")
        except Exception as e:
            logging.error(f"Numexpr evaluation failed for mask expression: '{expression}'")
            raise e # Re-raise for the recipe to catch
            
        # ... (rest of the method is unchanged) ...
        # 4. Create new data array
        new_values = target_q.value.copy()
        new_values[mask_indices] = np.nan
        new_data_col = DataColumnV2(new_values, target_q.unit)

        # 5. Create new SpectrumDataV2 using the immutable pattern
        core_cols = {
            'x': self._data.x, 'xmin': self._data.xmin, 'xmax': self._data.xmax,
            'y': self._data.y, 'dy': self._data.dy
        }
        aux_cols = deepcopy(self._data.aux_cols) 

        # 6. Overwrite the target column
        if target_col in core_cols:
            core_cols[target_col] = new_data_col
        elif target_col in aux_cols:
            aux_cols[target_col] = new_data_col

        # 7. Create the new data core
        new_data_core = SpectrumDataV2(
            x=core_cols['x'], xmin=core_cols['xmin'], xmax=core_cols['xmax'],
            y=core_cols['y'], dy=core_cols['dy'],
            aux_cols=aux_cols,
            meta=deepcopy(self._data.meta),
            z_rf=self._data.z_rf,  # Propagate
            z_em=self._data.z_em   # Propagate
        )

        # 8. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Mask: {target_col} where {expression}"]
        return SpectrumV2(data=new_data_core, history=new_history)

    def split(self, expression: str, extra_vars: Dict[str, np.ndarray] = None) -> 'SpectrumV2':
        """
        API: Splits the spectrum based on a boolean expression.
        Returns a NEW SpectrumV2 instance containing only the rows where
        the expression evaluates to True.
        """
        # 1. Get contexts for numexpr
        all_cols, local_dict = self._get_ne_contexts()
        if extra_vars:
            local_dict.update(extra_vars)

        try:
            # 2. Evaluate expression to get a boolean mask
            mask = ne.evaluate(expression, 
                               global_dict=_NE_GLOBAL_DICT,
                               local_dict=local_dict)
            if mask.dtype.kind != 'b':
                 raise TypeError("Split expression must return a boolean mask.")
        except Exception as e:
            logging.error(f"Numexpr evaluation failed for split: '{expression}'")
            raise e

        # 3. Check if the resulting spectrum would be empty
        if np.sum(mask) == 0:
            raise ValueError("Split expression resulted in an empty spectrum. Aborted.")

        # 4. Create new data core by boolean-indexing all columns
        new_data_cols = {}
        for col_name, quantity in all_cols.items():
            # Apply the mask to get the subset of data
            new_values = quantity.value[mask]
            new_data_cols[col_name] = DataColumnV2(new_values, quantity.unit)

        # 5. Build new SpectrumDataV2
        new_data_core = SpectrumDataV2(
            x=new_data_cols.pop('x'),
            xmin=new_data_cols.pop('xmin'),
            xmax=new_data_cols.pop('xmax'),
            y=new_data_cols.pop('y'),
            dy=new_data_cols.pop('dy'),
            aux_cols=new_data_cols, # Remaining items are aux_cols
            meta=deepcopy(self._data.meta),
            z_rf=self._data.z_rf,
            z_em=self._data.z_em
        )
        
        # 6. Return NEW SpectrumV2
        new_history = self.history + [f"Split with expression: {expression}"]
        return SpectrumV2(data=new_data_core, history=new_history)
    
    def calculate_running_std(self, 
                              input_col: str, 
                              output_col: str, 
                              window_pix: int) -> 'SpectrumV2':
        """
        API: Calculates a running RMS on a column and saves it as a new column.
        Returns a NEW SpectrumV2 instance.
        """
        
        all_cols = self.t._data_dict
        
        # 1. Find the target column
        if input_col not in all_cols:
            raise ValueError(f"Input column '{input_col}' not found.")
            
        input_col_data = all_cols[input_col]
        
        # 2. Check if it's numerical
        if input_col_data.dtype.kind not in 'fiu':
            raise TypeError(f"Cannot calculate RMS on non-numerical column '{input_col}'.")

        # 3. Call the pure function (h is half-window)
        h_window = int(window_pix / 2)
        rms_values = running_std(
            y=input_col_data.value,
            h=h_window
        )
        
        new_std_col = DataColumnV2(
            values=rms_values,
            unit=input_col_data.unit, # RMS has same unit as the input
            description=f"Running StdDev (window={window_pix}pix) of {input_col}"
        )

        # 4. Create new SpectrumDataV2 by adding the new column
        new_aux_cols = deepcopy(self._data.aux_cols)
        new_aux_cols[output_col] = new_std_col
        
        new_data = dataclasses.replace(self._data, aux_cols=new_aux_cols)
        
        # 5. Return new SpectrumV2
        new_history = self.history + [f"Calculated running RMS (window={window_pix}pix) on {input_col}"]
        return SpectrumV2(data=new_data, history=new_history)

    def smooth(self, sigma_kms: float) -> 'SpectrumV2':
        """
        API: Applies Gaussian smoothing to all flux-like columns (y, dy, cont, model).
        Returns a NEW SpectrumV2 instance.
        """
        
        # 1. Smooth the core y and dy columns
        new_y_values = smooth_spectrum(
            x=self._data.x.quantity,
            y=self._data.y.values,
            sigma_kms=sigma_kms,
            z_rf=self._data.z_rf
        )
        new_y_col = DataColumnV2(new_y_values, self._data.y.unit, "Smoothed Flux")

        new_dy_values = smooth_spectrum(
            x=self._data.x.quantity,
            y=self._data.dy.values,
            sigma_kms=sigma_kms,
            z_rf=self._data.z_rf
        )
        new_dy_col = DataColumnV2(new_dy_values, self._data.dy.unit, "Smoothed Error")

        # 2. Smooth flux-like auxiliary columns
        new_aux_cols = deepcopy(self._data.aux_cols)
        # --- RE-ADD 'model' TO THIS LIST ---
        cols_to_smooth = ['cont', 'model', 'telluric_model']
        
        for col_name in cols_to_smooth:
            if col_name in new_aux_cols:
                try:
                    old_col = new_aux_cols[col_name]
                    if old_col.values.dtype.kind in 'fiu':
                        logging.debug(f"Smoothing auxiliary column: {col_name}")
                        new_vals = smooth_spectrum(
                            self._data.x.quantity,
                            old_col.values,
                            sigma_kms,
                            self._data.z_rf
                        )
                        
                        new_aux_cols[col_name] = DataColumnV2(
                            new_vals, old_col.unit, old_col.description
                        )
                except Exception as e:
                    logging.warning(f"Failed to smooth aux col '{col_name}': {e}")
        
        # --- RE-NORMALIZATION BLOCK HAS BEEN REMOVED ---

        # 3. Create new SpectrumDataV2
        new_data = dataclasses.replace(
            self._data, 
            y=new_y_col, 
            dy=new_dy_col, 
            aux_cols=new_aux_cols
        )
        
        # 4. Return new SpectrumV2
        new_history = self.history + [f"Smoothed spectrum (sigma={sigma_kms} km/s)"]
        return SpectrumV2(data=new_data, history=new_history)

    def smooth_column(self, target_col: str, sigma_kms: float, renorm_model: bool = False) -> 'SpectrumV2':
        """
        API: Applies Gaussian smoothing to a single target column.
        If target_col is 'cont' and renorm_model is True, also updates 'model'.
        Returns a NEW SpectrumV2 instance.
        """
        all_cols = self.t._data_dict
        
        # 1. Find the target column
        if target_col not in all_cols:
            raise ValueError(f"Target column '{target_col}' not found.")
            
        old_col_data = all_cols[target_col]
        
        # 2. Check if it's numerical
        if old_col_data.dtype.kind not in 'fiu':
            raise TypeError(f"Cannot smooth non-numerical column '{target_col}'.")

        # 3. Call the pure function
        new_values = smooth_spectrum(
            x=self._data.x.quantity,
            y=old_col_data.value,
            sigma_kms=sigma_kms,
            z_rf=self._data.z_rf
        )
        
        new_data_col = DataColumnV2(
            values=new_values,
            unit=old_col_data.unit,
            description=f"Smoothed {target_col}"
        )
        
        # 4. Create new SpectrumDataV2 using the immutable pattern
        core_cols = {
            'x': self._data.x, 'xmin': self._data.xmin, 'xmax': self._data.xmax,
            'y': self._data.y, 'dy': self._data.dy
        }
        aux_cols = deepcopy(self._data.aux_cols) 

        # 5. Overwrite the target column
        if target_col in core_cols:
            core_cols[target_col] = new_data_col
        elif target_col in aux_cols:
            aux_cols[target_col] = new_data_col

        # --- ADD THIS BLOCK (Copied from fit_continuum) ---
        # 6. If we just smoothed 'cont', check if we should re-normalize 'model'
        if target_col == 'cont' and renorm_model:
            old_model_col = self._data.aux_cols.get('model')
            old_cont_vals = old_col_data.value 
            
            new_model_col = self._renormalize_model(
                old_cont_vals,
                new_values, # The new smoothed continuum
                old_model_col
            )
            if new_model_col:
                aux_cols['model'] = new_model_col
        # --- END ADD ---

        # 7. Create the new data core
        new_data_core = SpectrumDataV2(
            x=core_cols['x'], xmin=core_cols['xmin'], xmax=core_cols['xmax'],
            y=core_cols['y'], dy=core_cols['dy'],
            aux_cols=aux_cols,
            meta=deepcopy(self._data.meta),
            z_rf=self._data.z_rf,
            z_em=self._data.z_em
        )

        # 8. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Smoothed column: {target_col} (sigma={sigma_kms} km/s)"]
        return SpectrumV2(data=new_data_core, history=new_history)

    def rebin(self, xstart: Optional[au.Quantity], xend: Optional[au.Quantity], 
              dx: au.Quantity, kappa: Optional[float], 
              filling: float) -> 'SpectrumV2':
        """
        API: Rebins the spectrum to a new *uniform* grid.
        This now correctly interpolates all auxiliary columns.
        """
        
        original_x_unit = self._data.x.unit 
        z_rf = self._data.z_rf 

        if dx.unit.is_equivalent(au.km/au.s):
            target_calc_unit = au.km/au.s
        else:
            target_calc_unit = original_x_unit

        x_calc = convert_axis_velocity(self._data.x.quantity, z_rf, target_calc_unit)
        xmin_calc = convert_axis_velocity(self._data.xmin.quantity, z_rf, target_calc_unit)
        xmax_calc = convert_axis_velocity(self._data.xmax.quantity, z_rf, target_calc_unit)

        # 1. Rebin the flux-like columns
        x_new, xmin_new, xmax_new, y_new, dy_new = rebin_spectrum(
            x_calc, xmin_calc, xmax_calc, 
            self._data, 
            xstart, xend, dx, target_calc_unit,
            kappa, filling
        )
        
        # 2. Convert new grid back to original units (e.g., nm)
        x_final = convert_axis_velocity(x_new, z_rf, original_x_unit)
        xmin_final = convert_axis_velocity(xmin_new, z_rf, original_x_unit)
        xmax_final = convert_axis_velocity(xmax_new, z_rf, original_x_unit)
        
        # --- *** 3. THIS IS THE FIX: Interpolate Aux Columns *** ---
        new_aux_cols = {}
        x_final_nm = x_final.to_value(au.nm) # Get the new grid in nm for interpolation
        
        for name, col_data in self._data.aux_cols.items():
            
            # --- *** FIX 1: Check if 'values' is None from a *previous* rebin *** ---
            if col_data.values is None:
                logging.debug(f"Skipping aux col '{name}' as its data is None.")
                continue
            # --- *** END FIX 1 *** ---

            if col_data.values.dtype.kind in 'fiub': # Check for float, int, unsigned, or bool
                try:
                    logging.debug(f"Rebin: Interpolating aux column '{name}'...")
                    interp_values = self.get_resampled_column(
                        name, 
                        x_final_nm 
                    )
                    
                    # --- *** FIX 2: Check if interpolation returned None (e.g., for string) *** ---
                    if interp_values is not None:
                        new_aux_cols[name] = DataColumnV2(interp_values, col_data.unit)
                    # (If it is None, get_resampled_column already warned us, so we just drop it)
                    # --- *** END FIX 2 *** ---

                except Exception as e:
                    logging.warning(f"Failed to interpolate aux col '{name}' during rebin: {e}. Dropping.")
            else:
                logging.warning(f"Cannot rebin non-numerical aux col '{name}'. Dropping.")
        # --- *** END FIX *** ---

        # 4. Build new SpectrumDataV2
        new_data = SpectrumDataV2(
            x=DataColumnV2(x_final.value, x_final.unit),
            xmin=DataColumnV2(xmin_final.value, xmin_final.unit),
            xmax=DataColumnV2(xmax_final.value, xmax_final.unit),
            y=DataColumnV2(y_new.value, y_new.unit),
            dy=DataColumnV2(dy_new.value, dy_new.unit),
            aux_cols=new_aux_cols, 
            meta=self._data.meta, 
            z_rf=self._data.z_rf,  # Propagate
            z_em=self._data.z_em   # Propagate
        )
        
        # 5. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Rebinned spectrum (dx={dx})"]
        return SpectrumV2(data=new_data, history=new_history)

    def resample_on_grid(self, target_grid_spec: 'SpectrumV2', 
                         fill_value: float = np.nan) -> 'SpectrumV2':
        # ... (this method is unchanged, it correctly resamples *everything*) ...
        x_old = self.x.to_value(au.nm)
        x_new = target_grid_spec.x.to_value(au.nm)
        all_cols_old = self.t._data_dict
        x_col = DataColumnV2(x_new, target_grid_spec.x.unit)
        xmin_col = DataColumnV2(target_grid_spec.xmin.value, target_grid_spec.xmin.unit)
        xmax_col = DataColumnV2(target_grid_spec.xmax.value, target_grid_spec.xmax.unit)
        y_col_new = None
        dy_col_new = None
        aux_cols_new = {}
        for name, col_q in all_cols_old.items():
            if name in ('x', 'xmin', 'xmax'):
                continue
            if col_q.dtype.kind in 'fiu': 
                logging.debug(f"Resampling column '{name}'...")
                y_interp = np.interp(x_new, x_old, col_q.value,
                                     left=fill_value, right=fill_value)
                new_data_col = DataColumnV2(y_interp, col_q.unit)
                if name == 'y':
                    y_col_new = new_data_col
                elif name == 'dy':
                    dy_col_new = new_data_col
                else:
                    aux_cols_new[name] = new_data_col
            else:
                logging.warning(f"Cannot resample non-numerical column '{name}', it will be dropped.")
        new_data_core = SpectrumDataV2(
            x=x_col, xmin=xmin_col, xmax=xmax_col, 
            y=y_col_new, 
            dy=dy_col_new, 
            aux_cols=aux_cols_new, 
            meta=deepcopy(self.meta),
            z_rf=self._data.z_rf,  # Propagate
            z_em=self._data.z_em   # Propagate
        )
        new_history = self.history + [f"Resampled on grid of '{target_grid_spec.meta.get('name', 'unnamed')}'"]
        return SpectrumV2(data=new_data_core, history=new_history)

    def get_resampled_column(self, col_name: str, 
                             target_x_grid: np.ndarray, 
                             fill_value: float = np.nan) -> np.ndarray:
        """
        Lightweight API to interpolate a single column onto a new x-grid.
        Assumes target_x_grid is already in 'nm'.
        """
        all_cols = self.t._data_dict
        
        # 1. Get column to interpolate
        if col_name not in all_cols:
            raise ValueError(f"Column '{col_name}' not found in spectrum '{self.name}'.")
        
        col_to_interp = all_cols[col_name]
        
        # 2. Check that it's numerical
        if col_to_interp.dtype.kind not in 'fiu':
            raise TypeError(f"Cannot resample non-numerical column '{col_name}'.")
            
        # 3. Get old grid (must be in nm)
        x_old = self.x.to_value(au.nm)
        y_old = col_to_interp.value
        
        # 4. Perform interpolation (target_x_grid is already in nm)
        y_new = np.interp(target_x_grid, x_old, y_old,
                          left=fill_value, right=fill_value)
                          
        return y_new
    
    def calibrate(self, magnitudes: str) -> 'SpectrumV2':
        """
        API: Calibrates (and warps) the flux to match target magnitudes.
        """
        # 1. Get the pure numerical correction curve
        try:
            correction_curve = generate_calibration_curve(
                self.x, self.y, magnitudes
            )
        except Exception as e:
            logging.error(f"Calibration failed: {e}")
            raise e

        # 2. Define standard physical unit
        new_unit = au.erg / (au.cm**2 * au.s * au.Angstrom)

        # 3. Helper to apply curve and unit to a column
        def apply_corr(col: DataColumnV2) -> DataColumnV2:
            return DataColumnV2(col.values * correction_curve, new_unit, col.description)

        # 4. Apply to all relevant columns
        new_y = apply_corr(self._data.y)
        new_dy = apply_corr(self._data.dy)
        
        new_aux_cols = deepcopy(self._data.aux_cols)
        cols_to_scale = ['cont', 'model', 'cont_pl', 'telluric_model']
        for col_name in cols_to_scale:
            if col_name in new_aux_cols:
                 new_aux_cols[col_name] = apply_corr(new_aux_cols[col_name])

        # 5. Return new state
        new_data = dataclasses.replace(self._data, y=new_y, dy=new_dy, aux_cols=new_aux_cols)
        new_history = self.history + [f"Flux calibrated with: {magnitudes}"]
        return SpectrumV2(data=new_data, history=new_history)
    
    def find_absorbed(self, smooth_len_lya: au.Quantity, smooth_len_out: au.Quantity, 
                      kappa: float, template: bool) -> 'SpectrumV2':
        """
        API: Finds unabsorbed regions and adds/updates the 'mask_unabs' column.
        """
        # 1. Call the pure function with the new signature
        mask = find_absorbed_regions(
            x=self._data.x.quantity,
            y=self._data.y.values,
            dy=self._data.dy.values,
            z_em=self._data.z_em, # Get z_em from self
            smooth_len_lya=smooth_len_lya,
            smooth_len_out=smooth_len_out,
            kappa=kappa,
            template=template
        )
        
        # 2. Create new data core
        new_aux_cols = deepcopy(self._data.aux_cols)
        new_aux_cols['abs_mask'] = DataColumnV2(
            values=mask,
            unit=au.dimensionless_unscaled,
            description="Mask of absorbed regions (True=absorbed)"
        )
        
        new_data = dataclasses.replace(self._data, aux_cols=new_aux_cols)
        
        # 3. Return new SpectrumV2
        new_history = self.history + [f"Found absorbed regions (kappa={kappa})"]
        return SpectrumV2(data=new_data, history=new_history)
        
    def fit_continuum(self, fudge: float, smooth_std_kms: float, 
                      mask_col: str = 'abs_mask', renorm_model: bool = False) -> 'SpectrumV2':
        """
        API: Fits a continuum to a mask using V1 'interp-and-smooth' logic.
        """
        # 1. Get the mask
        if mask_col not in self._data.aux_cols:
            raise ValueError(f"Mask column '{mask_col}' not found. Run 'find_absorbed' first.")
        
        mask_abs = self._data.aux_cols[mask_col].values
        if mask_abs.dtype.kind != 'b':
            raise TypeError(f"Mask column '{mask_col}' is not a boolean mask.")

        # 2. Call the pure "interp-and-smooth" function
        cont_values = fit_continuum_interp(
            x=self._data.x.values,
            y=self._data.y.values,
            mask_abs=mask_abs,
            fudge=fudge,
            smooth_std_kms=smooth_std_kms,
            z_rf=self._data.z_rf # Pass z_rf for smoothing
        )
        
        # 3. Create new data core
        new_aux_cols = deepcopy(self._data.aux_cols)
        new_aux_cols['cont'] = DataColumnV2(
            values=cont_values,
            unit=self._data.y.unit, 
            description=f"Continuum fit (interp+smooth) to '{mask_col}'"
        )
        
        old_cont_col = self._data.aux_cols.get('cont')
        old_model_col = self._data.aux_cols.get('model')

        if renorm_model and old_cont_col is not None and old_model_col is not None:
            new_model_col = self._renormalize_model(
                old_cont_col.values,
                cont_values, # The new continuum
                old_model_col
            )
            if new_model_col:
                new_aux_cols['model'] = new_model_col

        new_data = dataclasses.replace(self._data, aux_cols=new_aux_cols)
        
        # 4. Return new SpectrumV2
        new_history = self.history + [f"Fitted continuum (smooth={smooth_std_kms} km/s)"]
        return SpectrumV2(data=new_data, history=new_history)

    def fit_powerlaw(self, regions_str: str, kappa: float = 3.0) -> 'SpectrumV2':
        """
        API: Fits a power-law to specified rest-frame regions and adds/updates
        the 'cont_pl' column.
        """
        
        # 1. Call the pure function
        pl_values = fit_powerlaw_to_regions(
            x=self._data.x.values,
            y=self._data.y.values,
            dy=self._data.dy.values,
            z_em=self._data.z_em, # Get z_em from self
            regions_str=regions_str,
            kappa=kappa
        )
        
        # 2. Create new data core
        new_aux_cols = deepcopy(self._data.aux_cols)
        new_aux_cols['cont_pl'] = DataColumnV2(
            values=pl_values,
            unit=self._data.y.unit, # Power-law has same unit as flux
            description=f"Power-law fit to regions: {regions_str}"
        )
        
        new_data = dataclasses.replace(self._data, aux_cols=new_aux_cols)
        
        # 3. Return new SpectrumV2
        new_history = self.history + [f"Fitted power-law (regions={regions_str})"]
        return SpectrumV2(data=new_data, history=new_history)

    def detect_absorptions(self, mask_col: str = 'abs_mask', min_pix: int = 3) -> 'SpectrumV2':
        """
        API: Detects absorption regions based on a continuum mask.
        Creates/updates the 'region_id' auxiliary column.
        """
        # 1. Get the mask
        if mask_col not in self._data.aux_cols:
            raise ValueError(f"Mask column '{mask_col}' not found. Run 'find_absorbed' first.")
        
        mask_data = self._data.aux_cols[mask_col].values
        
        if mask_data.dtype.kind == 'b':
            # It's already boolean, perfect.
            abs_mask = mask_data
        elif mask_data.dtype.kind in 'iu': # Integer or Unsigned Integer
            # It's likely 0/1 from a file load. Cast it.
            logging.debug(f"Casting integer column '{mask_col}' to boolean mask.")
            abs_mask = mask_data.astype(bool)
        else:
            # It's float or string, which is probably wrong for a mask.
            raise TypeError(f"Column '{mask_col}' has dtype '{mask_data.dtype}', expected boolean or integer mask.")
             
        # 2. Call pure function
        region_map, num_regions = detect_regions(abs_mask, min_pix=min_pix)
        
        logging.info(f"Detected {num_regions} absorption regions (min_pix={min_pix}).")

        # 3. Create new DataColumnV2 for the IDs
        # We use a dimensionless unit for IDs
        region_col = DataColumnV2(
            values=region_map,
            unit=au.dimensionless_unscaled,
            description=f"Detected regions from '{mask_col}' (ID 0 = continuum)"
        )
        
        # 4. Create new state
        new_aux_cols = deepcopy(self._data.aux_cols)
        new_aux_cols['region_id'] = region_col
        
        new_data = dataclasses.replace(self._data, aux_cols=new_aux_cols)
        new_hist = self.history + [f"Detected {num_regions} absorbers from {mask_col}"]
        
        return SpectrumV2(data=new_data, history=new_hist)
    
    def identify_lines(self, 
                       multiplet_list: List[str],
                       z_grid_dz: float,
                       modul: float,
                       sigma: float,
                       distance_pix: int,
                       prominence: Optional[float],
                       mask_col: str = 'abs_mask', 
                       min_pix_region: int = 3) -> 'SpectrumV2':
        """
        API: Orchestrator method to identify absorption lines.
        1. Detects absorption regions (creates 'region_id').
        2. Computes identification signals for all multiplets.
        3. Finds significant peaks in the signals.
        4. Cross-matches peaks to regions.
        5. Saves results to spec.meta['region_identifications'].
        """
        
        # This is the 'self' spectrum object
        spec = self
        z_em = spec._data.z_em # <<< *** GET z_em ***
        if z_em == 0.0:
            # This should be caught by the GUI, but as a safeguard:
            raise ValueError("z_em must be set before running identify_lines.")
        
        logging.debug(f"identify_lines: Starting identification...")
        logging.debug(f"  Params: sigma={sigma}, distance={distance_pix}, prominence={prominence}")
        logging.debug(f"  Params: mask_col={mask_col}, min_pix_region={min_pix_region}")

        print(mask_col)

        # 1. Silently create abs_mask if it doesn't exist
        if not spec.has_aux_column(mask_col):
            logging.info(f"Mask '{mask_col}' not found. Running 'find_absorbed' with defaults.")
            try:
                # Use reasonable defaults from the continuum recipe
                spec = spec.find_absorbed(
                    smooth_len_lya=5000.0 * au.km/au.s,
                    smooth_len_out=400.0 * au.km/au.s,
                    kappa=2.0,
                    template=False
                )
            except Exception as e:
                logging.error(f"Failed to auto-generate absorbed mask: {e}", exc_info=True)
                raise e

        # 2. Detect Absorption Regions (using the now-guaranteed mask)
        try:
            # This returns a *new* spec with the 'region_id' column
            # Call detect_regions *internally*, don't create the full 'region_id' column
            mask_data = spec.get_column(mask_col).value
            if mask_data.dtype.kind != 'b':
                 mask_data = mask_data.astype(bool)
            
            region_id_map_full, num_regions = detect_regions(mask_data, min_pix_region)
            x_nm = spec.x.to_value(au.nm)
            logging.info(f"identify_lines: Found {num_regions} regions to analyze.")
        except Exception as e:
            logging.error(f"Failed to detect absorption regions: {e}", exc_info=True)
            raise e # Re-raise for the recipe to catch
            
        # 3. Compute Identification Signals
        # Define a common, high-resolution z_grid
        z_min_spec = np.min((x_nm / 600.0) - 1.0) # Approx z_min at 600nm
        z_max_spec = np.max((x_nm / 121.6) - 1.0) # Approx z_max at Lya
        z_grid = np.arange(max(0.0, z_min_spec), z_max_spec, z_grid_dz)
        
        signals_dict = {}
        
        logging.info(f"Computing identification signals for: {multiplet_list}")
        for i, series_name in enumerate(multiplet_list):
            if series_name not in STANDARD_MULTIPLETS:
                logging.warning(f"Skipping '{series_name}': Not in STANDARD_MULTIPLETS.")
                continue
            
            signal = compute_identification_signal(
                spec, series_name, z_grid, modul
            )
            signals_dict[series_name] = signal

            if i == 0: # Log stats for the *first* signal only
                logging.debug(f" Stats for first signal ('{series_name}'):")
                try:
                    valid_signal = signal[np.isfinite(signal)]
                    logging.debug(f"   Signal (min/mean/max): {np.min(valid_signal):.2f} / {np.mean(valid_signal):.2f} / {np.max(valid_signal):.2f}")
                except Exception as e:
                    logging.debug(f"   Could not compute signal stats: {e}")

        # 4. Find Significant Peaks
        logging.info("Finding peaks in signals...")
        peaks_list = find_signal_peaks(
            signals_dict, z_grid, sigma, distance_pix, prominence
        )
        logging.info(f"Found {len(peaks_list)} candidate peaks.")
        
        # 4b - Apply Physical Constraints
        logging.info(f"Found {len(peaks_list)} raw peaks. Applying physical constraints...")
        
        # Get boundary wavelengths in OBSERVED frame (nm)
        lya_limit_obs_nm = (1.0 + z_em) * xem_d['Ly_lim'].to_value(au.nm)
        lya_obs_nm = (1.0 + z_em) * xem_d['Ly_a'].to_value(au.nm)
        
        filtered_peaks = []
        for (series_name, z_peak, score) in peaks_list:
            try:
                ref_trans = STANDARD_MULTIPLETS[series_name][0]
                x_peak_nm = (1.0 + z_peak) * xem_d[ref_trans].to_value(au.nm)

                # Rule 1: Reject *anything* below the Lyman Limit
                if x_peak_nm < lya_limit_obs_nm:
                    continue
                    
                # Rule 2: Penalize metals in the Ly-alpha forest
                is_h = is_hydrogen_line(ref_trans)
                is_in_forest = x_peak_nm < lya_obs_nm
                
                # Rule 2: Penalize metals in the Ly-alpha forest
                if not is_h and is_in_forest:
                    # Require a much higher score (e.g., 3x sigma)
                    if score < (sigma * 10.0): 
                        continue # Reject low-confidence metal in forest
                
                # Rule 3: Reject H lines redward of Lya emission
                if is_h and not is_in_forest:
                    continue # Reject H lines found redward of Lya
                
                filtered_peaks.append((series_name, z_peak, score))
                
            except Exception:
                continue # Skip peaks that fail lookup

        logging.info(f"Found {len(filtered_peaks)} physically valid peaks.")

        # 5. Cross-match Peaks with Regions
        region_identifications: Dict[int, List[Tuple[str, float]]] = {}
        regions_found: set[int] = set() # <<< *** ADDED for Point 3 ***
            
        for (series_name, z_peak, score) in filtered_peaks:
            try:
                # --- *** REFACTORED: Point 3 *** ---
                # Iterate over ALL transitions in the multiplet
                transitions_in_multiplet = STANDARD_MULTIPLETS[series_name]
                
                for ref_trans in transitions_in_multiplet:
                    x_peak_nm = (1.0 + z_peak) * xem_d[ref_trans].to_value(au.nm)
                    
                    # Find the spectrum index corresponding to this peak
                    idx = np.abs(x_nm - x_peak_nm).argmin()
                    
                    # Get the region ID at that index
                    rid = int(region_id_map_full[idx])
                    
                    # Assign this peak to the region
                    if rid > 0: # (Ignore background region 0)
                        regions_found.add(rid) # Add to set of all identified regions
                        
                        if rid not in region_identifications:
                            region_identifications[rid] = []
                        
                        # Add the main series name, not the individual transition
                        region_identifications[rid].append((series_name, float(score)))

            except Exception as e:
                logging.warning(f"Failed to cross-match peak {series_name} at z={z_peak}: {e}")

        # 5b. De-duplicate entries per region
        for rid in region_identifications:
            # Sort identifications by score (descending)
            region_identifications[rid].sort(key=lambda x: x[1], reverse=True)
            # Keep only unique series names (highest score version)
            unique_ids = []
            seen_series = set()
            for (series, score) in region_identifications[rid]:
                if series not in seen_series:
                    unique_ids.append((series, score))
                    seen_series.add(series)
            region_identifications[rid] = unique_ids

        # 6. Create Visualization Column
        new_meta = spec.meta
        new_aux_cols = deepcopy(spec._data.aux_cols)

        if regions_found:
            # Create a copy of the region map
            vis_map = region_id_map_full.copy()
            
            # Find RIDs to set to 0 (unidentified)
            all_rids = np.unique(vis_map)
            rids_to_zero = [rid for rid in all_rids if rid > 0 and rid not in regions_found]

            if rids_to_zero:
                mask_to_zero = np.isin(vis_map, rids_to_zero)
                vis_map[mask_to_zero] = 0
            
            # Add this new map as the 'abs_ids' column
            vis_col = DataColumnV2(
                values=vis_map,
                unit=au.dimensionless_unscaled,
                description="Absorption regions with identifications"
            )
            new_aux_cols['abs_ids'] = vis_col # <<< *** RENAMED ***

        # 7. Save Metadata (as JSON)
        if not region_identifications:
            # If the dict is empty, save None (which FITS headers can handle)
            new_meta['region_identifications'] = None
        else:
            try:
                # If dict is not empty, serialize to a JSON string
                new_meta['region_identifications'] = json.dumps(region_identifications)
            except Exception as e:
                logging.error(f"Failed to serialize region identifications to JSON: {e}")
                new_meta['region_identifications'] = None # Save None on error
        
        final_data_core = dataclasses.replace(
            spec._data, 
            meta=new_meta,
            aux_cols=new_aux_cols
        )
        
        # 8. Create Final Data Core
        final_data_core = dataclasses.replace(
            spec._data, # <<< *** Use original spec's data ***
            meta=new_meta,
            aux_cols=new_aux_cols 
        )
        
        # 9. Return new SpectrumV2
        new_history = self.history + [f"Identified {len(filtered_peaks)} candidates"]
        return SpectrumV2(data=final_data_core, history=new_history)

    def x_convert(self, z_rf: float, xunit: au.Unit) -> 'SpectrumV2':
        """
        API: Converts the X-axis units and returns a NEW SpectrumV2 instance.
        (Replaces V1's mutable _x_convert)
        """
        
        # 1. Execute pure conversion logic
        x_new, xmin_new, xmax_new = convert_x_axis(self._data, z_rf, xunit)
        
        # 2. Build new SpectrumDataV2, preserving the core vertical data
        new_data = SpectrumDataV2(
            x=DataColumnV2(x_new.value, xunit),
            xmin=DataColumnV2(xmin_new.value, xunit),
            xmax=DataColumnV2(xmax_new.value, xunit),
            y=self._data.y, 
            dy=self._data.dy,
            aux_cols=self._data.aux_cols,
            meta=self._data.meta,
            z_rf=self._data.z_rf,  # Propagate
            z_em=self._data.z_em   # Propagate
        )
        
        # 3. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"X-axis converted to {xunit} at z_rf={z_rf}"]
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
            z_rf=self._data.z_rf,  # Propagate
            z_em=self._data.z_em   # Propagate
        )
        
        # 3. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Y-axis converted to {yunit}"]
        return SpectrumV2(data=new_data, history=new_history)