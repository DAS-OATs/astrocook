import astropy.constants as const
import astropy.units as au
from copy import deepcopy
import dataclasses
import json
import logging
import numexpr as ne
import numpy as np
from scipy.ndimage import gaussian_filter1d
from typing import Any, Dict, List, Optional, Union

from astrocook.core.atomic_data import STANDARD_MULTIPLETS, xem_d
from astrocook.core.photometry import generate_calibration_curve
from astrocook.core.spectrum_operations import (
    convert_axis_velocity, 
    convert_x_axis, convert_y_axis, detect_regions, find_absorbed_regions, _find_kinematic_doublet_candidates, fit_powerlaw_to_regions, 
    merge_regions_by_velocity, rate_doublet_candidate, smooth_spectrum, rebin_spectrum, running_std,
)
from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.structures import SpectrumDataV2, DataColumnV2
from astrocook.legacy.frame import Frame as FrameV1
from astrocook.legacy.spectrum import Spectrum as SpectrumV1

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
    Composite container for spectral data and operations.
    
    This class wraps the raw data (flux, wavelength, error) and provides
    immutable operations that return new SpectrumV2 instances. It serves as
    the primary data structure passed between Recipes and the GUI.

    Parameters
    ----------
    data : SpectrumDataV2
        The core data container holding the Astropy Quantities.
    history : list, optional
        A list of strings describing the modification history of this spectrum.
    """
    
    def __init__(self, data: SpectrumDataV2, history: list = None):
        self._data = data
        self.history = history if history is not None else []

    # Accessors
    @property
    def x(self) -> au.Quantity:
        """The wavelength (or velocity) array."""
        return self._data.x.quantity

    @property
    def xmin(self) -> au.Quantity:
        """The lower bound of the wavelength bins."""
        return self._data.xmin.quantity

    @property
    def xmax(self) -> au.Quantity:
        """The upper bound of the wavelength bins."""
        return self._data.xmax.quantity

    @property
    def y(self) -> au.Quantity:
        """The flux density array."""
        return self._data.y.quantity

    @property
    def dy(self) -> au.Quantity:
        """The flux error (1-sigma) array."""
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
    def norm(self) -> Optional[au.Quantity]:
        """
        Returns the normalization vector (Continuum * Telluric).
        
        Logic:
        1. If cont exists: 
           - Returns cont * telluric (if telluric exists).
           - Returns cont (otherwise).
        2. If cont is missing:
           - Returns 1.0 (if flux median is ~1, assuming pre-normalized).
           - Returns None (if flux is physical, signaling need for continuum fitting).
        """
        cont = self.get_column('cont')
        telluric = self.get_column('telluric_model')
        
        # Case 1: Continuum Exists
        if cont is not None:
            if telluric is not None:
                # Multiply Continuum (Quantity) by Telluric (Array/Quantity)
                # Assuming telluric is dimensionless 0..1 transmission
                return cont * telluric.value
            return cont
            
        # Case 2: No Continuum (Check for implicit normalization)
        # Heuristic: If median flux is order of unity (0.1 - 10), assume normalized.
        # If median flux is 1e-17, assume physical.
        y_median = np.nanmedian(self.y.value)
        if 0.1 < y_median < 10.0:
             return np.ones_like(self.y.value) * self.y.unit
             
        # Case 3: Physical Flux with no Continuum -> Needs Fitting
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
        """A dictionary of metadata (header keywords, redshifts, etc.)."""
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
    
    def update_column(self, name: str, values: np.ndarray, 
                      unit: Optional[au.Unit] = None,
                      description: str = "") -> 'SpectrumV2':
        """
        Returns a NEW SpectrumV2 instance with the specified column updated or added.
        Handles both core columns (y, dy) and auxiliary columns (cont, model, etc.).
        """
        # 1. Determine Unit (inherit if not provided and col exists)
        target_unit = unit
        if target_unit is None:
            if name in self._data.aux_cols:
                target_unit = self._data.aux_cols[name].unit
            elif name in ['y', 'dy']:
                target_unit = getattr(self._data, name).unit
            else:
                # Default fallback: match Flux unit
                target_unit = self._data.y.unit
        
        # 2. Create the DataColumnV2 object
        # (DataColumnV2 is already imported in your spectrum.py)
        new_col = DataColumnV2(values, target_unit, description)
        
        # 3. Prepare new data structure (Immutable Pattern)
        new_data = None
        
        # Check Core Columns
        if name == 'y':
            new_data = dataclasses.replace(self._data, y=new_col)
        elif name == 'dy':
            new_data = dataclasses.replace(self._data, dy=new_col)
        # We purposely exclude x/xmin/xmax here to prevent grid mismatches
        else:
            # Auxiliary Column
            new_aux = deepcopy(self._data.aux_cols)
            new_aux[name] = new_col
            new_data = dataclasses.replace(self._data, aux_cols=new_aux)
            
        # 4. Return new SpectrumV2
        new_history = self.history + [f"Updated column '{name}'"]
        return SpectrumV2(data=new_data, history=new_history)

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
    
    def with_properties(self, z_em: float, resol: float = 0.0, 
                        object_name: Optional[str] = None) -> 'SpectrumV2':
        """
        API: Returns a NEW SpectrumV2 instance with updated properties.
        """
        new_meta = deepcopy(self._data.meta)
        
        history_entries = [f"z_em={z_em}"]

        # Update Resolution
        if resol > 0:
            new_meta['resol'] = resol
            c_kms = 299792.458
            new_meta['resol_fwhm'] = c_kms / resol
            history_entries.append(f"R={resol}")

        # Update Object Name
        if object_name is not None:
            new_meta['OBJECT'] = object_name
            # Also update standard ESO keyword if present, for consistency
            if 'HIERARCH ESO OBS TARG NAME' in new_meta:
                new_meta['HIERARCH ESO OBS TARG NAME'] = object_name
            history_entries.append(f"Object='{object_name}'")

        # Update Core Data
        new_data = dataclasses.replace(self._data, z_em=z_em, resol=resol, meta=new_meta)
        
        new_history = self.history + [f"Set properties ({', '.join(history_entries)})"]
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
    
    def sanitize_legacy_tellurics(self) -> 'SpectrumV2':
        """
        Patches legacy spectra where 'cont' included tellurics and 'cont_no_telluric' 
        stored the intrinsic continuum.
        
        Logic:
        1. If 'cont_no_telluric' exists:
           - Overwrite 'cont' with 'cont_no_telluric' (Intrinsic).
           - Remove 'cont_no_telluric'.
        2. Returns a new SpectrumV2 instance (or self if no changes needed).
        """
        # Check if the legacy column exists
        if 'cont_no_telluric' not in self._data.aux_cols:
            return self
            
        logging.info("Sanitizing legacy spectrum: Swapping 'cont_no_telluric' into 'cont' to fix double-counting.")
        
        new_aux_cols = deepcopy(self._data.aux_cols)
        
        # 1. Promote 'cont_no_telluric' to be the canonical 'cont'
        # This ensures 'cont' is purely intrinsic.
        new_aux_cols['cont'] = new_aux_cols['cont_no_telluric']
        #new_aux_cols['cont'].description = "Intrinsic Continuum (Sanitized)"
        
        # 2. Remove the redundant legacy column
        del new_aux_cols['cont_no_telluric']
        
        # 3. Create new state
        new_data = dataclasses.replace(self._data, aux_cols=new_aux_cols)
        
        # Note: We add a history entry so the user knows this happened
        return SpectrumV2(data=new_data, history=self.history + ["Sanitized legacy telluric columns"])
    
    # --- Methods implementing the API ---

    def apply_expression(self, target_col: str, expression: str, 
                         extra_vars: Dict[str, np.ndarray] = None) -> 'SpectrumV2': # <<< MODIFIED
        """
        Applies a numerical expression to columns and returns a new spectrum.

        Parameters
        ----------
        target_col : str
            The name of the column to create or overwrite.
        expression : str
            A NumExpr-compatible string (e.g. ``"y / cont"``).
        extra_vars : dict, optional
            Additional variables to pass to the evaluation context.

        Returns
        -------
        SpectrumV2
            A new instance with the computed column.
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
            z_em=self._data.z_em   # Propagate
        )

        # 8. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Apply: {target_col} = {expression}"]
        return SpectrumV2(data=new_data_core, history=new_history)

    def mask_expression(self, target_col: str, expression: str, 
                        extra_vars: Dict[str, np.ndarray] = None) -> 'SpectrumV2': # <<< MODIFIED
        """
        Masks a column based on a boolean expression (sets values to NaN).

        Parameters
        ----------
        target_col : str
            The name of the column to mask.
        expression : str
            A boolean expression (e.g. ``"x > 500"``). True values will be masked.
        extra_vars : dict, optional
            Additional context variables.

        Returns
        -------
        SpectrumV2
            A new instance with the masked column.
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
            z_em=self._data.z_em   # Propagate
        )

        # 8. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Mask: {target_col} where {expression}"]
        return SpectrumV2(data=new_data_core, history=new_history)
    
    def calc_intersection_expression(self, others_specs: List['SpectrumV2'], 
                                     z_target: float, trans_self: str, 
                                     trans_others: List[str], window_kms: float) -> Optional[str]:
        """
        Calculates the wavelength range corresponding to the velocity intersection 
        of this spectrum and a list of others.

        Returns a string expression suitable for `edit.split`, e.g., "(x > 400) & (x < 500)".
        """
        
        # Helper logic (previously in batch_driver)
        def get_bounds(spec, trans, z):
            if trans not in ATOM_DATA: return -np.inf, np.inf
            x_vals = spec.x.value
            if len(x_vals) == 0: return 0.0, 0.0
            
            lam_rest = ATOM_DATA[trans]['wave'] / 10.0 # nm
            lam_obs = lam_rest * (1 + z)
            c_kms = 299792.458
            
            v = c_kms * (x_vals - lam_obs) / lam_obs
            return np.min(v), np.max(v)

        v_mins, v_maxs = [], []
        
        # 1. Self
        v1, v2 = get_bounds(self, trans_self, z_target)
        v_mins.append(v1); v_maxs.append(v2)
        
        # 2. Others
        for spec, trans in zip(others_specs, trans_others):
            v1, v2 = get_bounds(spec, trans, z_target)
            v_mins.append(v1); v_maxs.append(v2)
            
        # 3. Intersection
        common_min = max(v_mins)
        common_max = min(v_maxs)
        
        # 4. Windowing
        final_min = max(common_min, -window_kms)
        final_max = min(common_max, window_kms)
        
        if final_max <= final_min: return None
        
        # 5. Convert back to Wavelength (Self)
        lam_rest_self = ATOM_DATA[trans_self]['wave'] / 10.0
        lam_obs_self = lam_rest_self * (1 + z_target)
        c_kms = 299792.458
        
        l1 = lam_obs_self * (1 + final_min/c_kms)
        l2 = lam_obs_self * (1 + final_max/c_kms)
        
        return f"(x > {l1:.5f}) & (x < {l2:.5f})"

    def stitch(self, other_specs: List['SpectrumV2'], sort: bool = True) -> 'SpectrumV2':
        """
        Merges this spectrum with a list of other spectra.
        Concatenates x, y, dy, and all auxiliary columns.
        """
        from copy import deepcopy
        
        all_specs = [self] + other_specs
        
        xs, ys, dys = [], [], []
        aux_map = {k: [] for k in self._data.aux_cols.keys()}
        
        # Target Units
        x_u = self.x.unit
        y_u = self.y.unit

        for s in all_specs:
            # Core Units
            xs.append(s.x.to(x_u).value)
            ys.append(s.y.to(y_u).value)
            dys.append(s.dy.to(y_u).value)
            
            # Aux Columns
            for col_name in aux_map:
                if col_name in s._data.aux_cols:
                    t_unit = self._data.aux_cols[col_name].unit
                    val = s.get_column(col_name).to(t_unit).value
                    aux_map[col_name].append(val)
                else:
                    # Pad missing columns with NaN
                    aux_map[col_name].append(np.full(len(s.x), np.nan))

        # Concatenate
        x_final = np.concatenate(xs)
        y_final = np.concatenate(ys)
        dy_final = np.concatenate(dys)
        
        if sort:
            idx = np.argsort(x_final)
            x_final = x_final[idx]
            y_final = y_final[idx]
            dy_final = dy_final[idx]
        else:
            idx = np.arange(len(x_final))

        # Build Aux Columns
        from astrocook.core.structures import DataColumnV2
        new_aux = {}
        for col_name, arrays in aux_map.items():
            arr = np.concatenate(arrays)
            if sort: arr = arr[idx]
            # Copy metadata from self
            orig_col = self._data.aux_cols[col_name]
            new_aux[col_name] = DataColumnV2(arr, orig_col.unit, orig_col.description)

        # Build New Object
        from ..io.loaders import _auto_limits
        xmin, xmax = _auto_limits(x_final, x_u)
        new_meta = deepcopy(self.meta)
        new_meta['stitched'] = True
        
        new_data = SpectrumDataV2(
            x=DataColumnV2(x_final, x_u),
            xmin=xmin, xmax=xmax,
            y=DataColumnV2(y_final, y_u),
            dy=DataColumnV2(dy_final, y_u),
            aux_cols=new_aux,
            meta=new_meta,
            z_em=self._data.z_em,
            resol=self._data.resol
        )
        
        return SpectrumV2(new_data)

    def split(self, expression: str, extra_vars: Dict[str, np.ndarray] = None) -> 'SpectrumV2':
        """
        Splits the spectrum based on a boolean expression.

        Parameters
        ----------
        expression : str
            A boolean expression defining the rows to keep (e.g. ``"x > 121.6"``).
        extra_vars : dict, optional
            Additional context variables.

        Returns
        -------
        SpectrumV2
            A new instance containing only the rows where the expression evaluates to True.
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
            z_em=self._data.z_em
        )
        
        # 6. Return NEW SpectrumV2
        new_history = self.history + [f"Split with expression: {expression}"]
        return SpectrumV2(data=new_data_core, history=new_history)
    
    def scale_flux(self, scale_factor: Union[float, np.ndarray]) -> 'SpectrumV2':
        """
        Scales the flux (y, dy) and relevant auxiliary columns (cont, model) 
        by a multiplicative factor. Returns a new SpectrumV2 instance.
        """
        # 1. Scale core columns
        # Handle cases where scale_factor is scalar or array
        new_y = DataColumnV2(self.y.value * scale_factor, self.y.unit)
        new_dy = DataColumnV2(self.dy.value * scale_factor, self.dy.unit)

        # 2. Scale flux-like auxiliary columns
        new_aux_cols = deepcopy(self._data.aux_cols)
        cols_to_scale = ['cont', 'model']
        
        for name in cols_to_scale:
            if name in new_aux_cols:
                col = new_aux_cols[name]
                # Ensure we propagate description and unit
                new_aux_cols[name] = DataColumnV2(
                    col.values * scale_factor, 
                    col.unit, 
                    col.description
                )
        
        # 3. Create new immutable data core
        new_data = dataclasses.replace(
            self._data,
            y=new_y,
            dy=new_dy,
            aux_cols=new_aux_cols
        )
        
        # 4. Return new SpectrumV2
        # (Optional: summarize scalar or mean of array for history)
        if np.size(scale_factor) == 1:
            info = f"{scale_factor:.3f}"
        else:
            info = f"array(mean={np.mean(scale_factor):.3f})"
            
        new_history = self.history + [f"Scaled flux by {info}"]
        return SpectrumV2(data=new_data, history=new_history)

    def calculate_running_std(self, 
                              input_col: str, 
                              output_col: str, 
                              window_pix: int) -> 'SpectrumV2':
        """
        Calculates a running RMS on a column and saves it as a new column.

        Parameters
        ----------
        input_col : str
            Name of the column to calculate statistics on (e.g. 'y').
        output_col : str
            Name of the output column (e.g. 'running_std').
        window_pix : int
            Width of the window in pixels.

        Returns
        -------
        SpectrumV2
            A new instance with the added RMS column.
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
        Applies Gaussian smoothing to all flux-like columns (y, dy, cont, model).

        Parameters
        ----------
        sigma_kms : float
            Standard deviation of the Gaussian kernel in km/s.

        Returns
        -------
        SpectrumV2
            A new instance with smoothed columns.
        """
        
        # 1. Smooth the core y and dy columns
        new_y_values = smooth_spectrum(
            x=self._data.x.quantity,
            y=self._data.y.values,
            sigma_kms=sigma_kms,
            z_ref=self._data.z_em
        )
        new_y_col = DataColumnV2(new_y_values, self._data.y.unit, "Smoothed Flux")

        new_dy_values = smooth_spectrum(
            x=self._data.x.quantity,
            y=self._data.dy.values,
            sigma_kms=sigma_kms,
            z_ref=self._data.z_em
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
                            self._data.z_em
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
        Applies Gaussian smoothing to a single target column.

        Parameters
        ----------
        target_col : str
            Name of the column to smooth.
        sigma_kms : float
            Standard deviation of the Gaussian kernel in km/s.
        renorm_model : bool, optional
            If True and ``target_col`` is 'cont', the 'model' column will be re-normalized.

        Returns
        -------
        SpectrumV2
            A new instance with the smoothed column.
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
            z_ref=self._data.z_em
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
            z_em=self._data.z_em
        )

        # 8. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Smoothed column: {target_col} (sigma={sigma_kms} km/s)"]
        return SpectrumV2(data=new_data_core, history=new_history)

    def rebin(self, xstart: Optional[au.Quantity], xend: Optional[au.Quantity], 
              dx: au.Quantity, kappa: Optional[float], 
              filling: float) -> 'SpectrumV2':
        """
        Rebins the spectrum to a new uniform grid.

        Parameters
        ----------
        xstart : Quantity, optional
            Start wavelength. If None, uses data minimum.
        xend : Quantity, optional
            End wavelength. If None, uses data maximum.
        dx : Quantity
            Step size (bin width). Can be wavelength or velocity.
        kappa : float, optional
            Sigma clipping threshold for rejecting outliers in bins.
        filling : float
            Value to fill gaps (default: NaN).

        Returns
        -------
        SpectrumV2
            A new instance with the rebinned data.
        """
        
        original_x_unit = self._data.x.unit 
        z_ref = self._data.z_em if self._data.z_em is not None else 0.0 

        if dx.unit.is_equivalent(au.km/au.s):
            target_calc_unit = au.km/au.s
        else:
            target_calc_unit = original_x_unit

        x_calc = convert_axis_velocity(self._data.x.quantity, z_ref, target_calc_unit)
        xmin_calc = convert_axis_velocity(self._data.xmin.quantity, z_ref, target_calc_unit)
        xmax_calc = convert_axis_velocity(self._data.xmax.quantity, z_ref, target_calc_unit)

        # 1. Rebin the flux-like columns
        x_new, xmin_new, xmax_new, y_new, dy_new = rebin_spectrum(
            x_calc, xmin_calc, xmax_calc, 
            self._data, 
            xstart, xend, dx, target_calc_unit,
            kappa, filling
        )
        
        # 2. Convert new grid back to original units (e.g., nm)
        x_final = convert_axis_velocity(x_new, z_ref, original_x_unit)
        xmin_final = convert_axis_velocity(xmin_new, z_ref, original_x_unit)
        xmax_final = convert_axis_velocity(xmax_new, z_ref, original_x_unit)
        
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

        # Create a copy of the metadata to avoid modifying the original session
        new_meta = self._data.meta.copy()

        # Remove the 'plot_style' tag so the GUI reverts to the default (Line)
        if 'plot_style' in new_meta:
            del new_meta['plot_style']

        # 4. Build new SpectrumDataV2
        new_data = SpectrumDataV2(
            x=DataColumnV2(x_final.value, x_final.unit),
            xmin=DataColumnV2(xmin_final.value, xmin_final.unit),
            xmax=DataColumnV2(xmax_final.value, xmax_final.unit),
            y=DataColumnV2(y_new.value, y_new.unit),
            dy=DataColumnV2(dy_new.value, dy_new.unit),
            aux_cols=new_aux_cols, 
            meta=new_meta, 
            z_em=self._data.z_em   # Propagate
        )
        
        # 5. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Rebinned spectrum (dx={dx})"]
        return SpectrumV2(data=new_data, history=new_history)
    
    def coadd(self, others: List['SpectrumV2'], 
              xstart: Optional[au.Quantity], xend: Optional[au.Quantity], 
              dx: au.Quantity, kappa: Optional[float] = 5.0, 
              equalize_order: int = 0) -> 'SpectrumV2':
        """
        Co-adds this spectrum with a list of other spectra onto a new common grid.
        
        Handles scaling, unit conversion, "drizzle" rebinning, and 
        converts the result back to the original wavelength unit.
        """
        from astrocook.core.spectrum_operations import (
            compute_flux_scaling, 
            rebin_spectrum, 
            convert_axis_velocity
        )
        from astrocook.core.structures import SpectrumDataV2, DataColumnV2

        # 1. Setup Reference and List
        all_specs = [self] + others
        
        # Use self as the reference for units, metadata, and scaling baseline
        ref_spec = self 
        ref_x = ref_spec.x.value
        ref_y = ref_spec.y.value
        original_x_unit = ref_spec.x.unit
        
        # Calculation Unit (Target) comes from dx (usually km/s)
        calc_unit = dx.unit
        z_em = ref_spec.meta.get('z_em', 0.0)

        # 2. Accumulators
        big_x_list = []
        big_xmin_list = []
        big_xmax_list = []
        big_y_list = []
        big_dy_list = []

        for i, spec in enumerate(all_specs):
            # A. Flux Scaling (Skip for self/index 0)
            if i > 0 and equalize_order >= 0:
                logging.info(f"  Scaling spectrum {i} to reference...")
                # Compute scaling in original frame
                scale_model = compute_flux_scaling(
                    ref_x, ref_y, 
                    spec.x.value, spec.y.value, 
                    order=equalize_order
                )
                y_val = spec.y.value * scale_model
                dy_val = spec.dy.value * scale_model
            else:
                y_val = spec.y.value
                dy_val = spec.dy.value

            # B. Pre-Convert to Calculation Unit (e.g. km/s)
            # This handles the "negative pixel" issue for overlapping orders
            # by converting safely before stitching.
            x_conv = convert_axis_velocity(spec.x, z_em, calc_unit).value
            xmin_conv = convert_axis_velocity(spec.xmin, z_em, calc_unit).value
            xmax_conv = convert_axis_velocity(spec.xmax, z_em, calc_unit).value
            
            # C. Append
            big_x_list.append(x_conv)
            big_xmin_list.append(xmin_conv)
            big_xmax_list.append(xmax_conv)
            big_y_list.append(y_val)
            big_dy_list.append(dy_val)

        # 3. Concatenate (The "Manual Stitch")
        flat_x = np.concatenate(big_x_list)
        flat_xmin = np.concatenate(big_xmin_list)
        flat_xmax = np.concatenate(big_xmax_list)
        flat_y = np.concatenate(big_y_list)
        flat_dy = np.concatenate(big_dy_list)

        # 4. Convert Bounds to Calculation Unit
        if xstart is not None:
            xstart_conv = convert_axis_velocity(xstart, z_em, calc_unit)
        else:
            xstart_conv = None

        if xend is not None:
            xend_conv = convert_axis_velocity(xend, z_em, calc_unit)
        else:
            xend_conv = None

        # 5. Create Temp Data Container (in Calculation Units)
        temp_data = SpectrumDataV2(
            x=DataColumnV2(flat_x, calc_unit),
            xmin=DataColumnV2(flat_xmin, calc_unit),
            xmax=DataColumnV2(flat_xmax, calc_unit),
            y=DataColumnV2(flat_y, self.y.unit),
            dy=DataColumnV2(flat_dy, self.dy.unit),
            meta=self.meta
        )

        # 6. Rebin (Vectorized Drizzle)
        # Returns new arrays in `calc_unit` (km/s)
        reb_x, reb_xmin, reb_xmax, reb_y, reb_dy = rebin_spectrum(
            temp_data.x.quantity, temp_data.xmin.quantity, temp_data.xmax.quantity,
            temp_data,
            xstart_conv, xend_conv,
            dx, calc_unit,
            kappa, filling=np.nan
        )

        # 7. Convert Back to Original Unit (e.g. nm)
        # We transform the result back to the domain of the input data
        if not reb_x.unit.is_equivalent(original_x_unit):
            logging.info(f"  Converting result back to {original_x_unit}...")
            final_x = convert_axis_velocity(reb_x, z_em, original_x_unit)
            final_xmin = convert_axis_velocity(reb_xmin, z_em, original_x_unit)
            final_xmax = convert_axis_velocity(reb_xmax, z_em, original_x_unit)
        else:
            final_x = reb_x
            final_xmin = reb_xmin
            final_xmax = reb_xmax

        # Create a copy of the metadata to avoid modifying the original session
        new_meta = self.meta.copy()
        
        # Remove the 'plot_style' tag so the GUI reverts to the default (Line)
        if 'plot_style' in new_meta:
            del new_meta['plot_style']
            
        # Optional: Add a note about the co-add origin
        new_meta['ORIGIN'] = 'coadd'

        # 8. Wrap in Logic Class
        final_spec_data = SpectrumDataV2(
            x=DataColumnV2(final_x.value, final_x.unit), 
            xmin=DataColumnV2(final_xmin.value, final_xmin.unit), 
            xmax=DataColumnV2(final_xmax.value, final_xmax.unit),
            y=DataColumnV2(reb_y.value, reb_y.unit), 
            dy=DataColumnV2(reb_dy.value, reb_dy.unit),
            meta=new_meta 
        )
        
        return self.__class__(final_spec_data)

    def resample_on_grid(self, target_grid_spec: 'SpectrumV2', 
                         fill_value: float = np.nan) -> 'SpectrumV2':
        """
        Resamples this spectrum onto the exact wavelength grid of another spectrum.

        Parameters
        ----------
        target_grid_spec : SpectrumV2
            The spectrum whose grid defines the target wavelengths.
        fill_value : float
            Value to fill outside the interpolation range (default: NaN).

        Returns
        -------
        SpectrumV2
            A new instance interpolated to the target grid.
        """
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
            z_em=self._data.z_em   # Propagate
        )
        new_history = self.history + [f"Resampled on grid of '{target_grid_spec.meta.get('name', 'unnamed')}'"]
        return SpectrumV2(data=new_data_core, history=new_history)

    def get_resampled_column(self, col_name: str, 
                             target_x_grid: np.ndarray, 
                             fill_value: float = np.nan) -> np.ndarray:
        """
        Lightweight API to interpolate a single column onto a new x-grid.
        
        Parameters
        ----------
        col_name : str
            Name of the column to interpolate.
        target_x_grid : ndarray
            Target grid in nanometers.
        
        Returns
        -------
        ndarray
            Interpolated values.
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
        Calibrates (and warps) the flux to match target photometric magnitudes.

        Parameters
        ----------
        magnitudes : str
            Comma-separated pairs of Filter=Mag (e.g. ``"SDSS_g=17.5, SDSS_r=17.0"``).

        Returns
        -------
        SpectrumV2
            A new instance with calibrated flux.
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
        API Method: Identifies absorbed regions.
        Calls pure function: find_absorbed_regions
        """
        # Call Pure Function
        mask = find_absorbed_regions(
            x=self._data.x.quantity,
            y=self._data.y.values,
            dy=self._data.dy.values,
            z_em=self._data.z_em,
            smooth_len_lya=smooth_len_lya,
            smooth_len_out=smooth_len_out,
            kappa=kappa,
            template=template
        )
        
        # State Update
        new_aux_cols = deepcopy(self._data.aux_cols)
        new_aux_cols['abs_mask'] = DataColumnV2(
            values=mask,
            unit=au.dimensionless_unscaled,
            description="Mask of absorbed regions (True=absorbed)"
        )
        new_data = dataclasses.replace(self._data, aux_cols=new_aux_cols)
        return SpectrumV2(data=new_data, history=self.history + [f"Found absorbed (kappa={kappa})"])
        
    def fit_continuum(self, fudge: Union[float, str], smooth_std_kms: float, 
                      mask_col: str = 'abs_mask', renorm_model: bool = False) -> 'SpectrumV2':
        """
        Fits a continuum to the unabsorbed regions.
        Orchestrates: Interpolate -> Auto-Fudge -> Smooth -> Renormalize Model.
        """
        from astrocook.core.spectrum_operations import (
            interpolate_continuum_mask, 
            compute_auto_fudge,
            smooth_spectrum
        )

        # 1. Validation & Setup
        if mask_col not in self._data.aux_cols:
            raise ValueError(f"Mask column '{mask_col}' not found. Run 'find_absorbed' first.")
        
        mask_abs = self._data.aux_cols[mask_col].values.astype(bool)
        
        # Check for Telluric Model
        telluric_col = self._data.aux_cols.get('telluric_model')
        telluric_vals = telluric_col.values if telluric_col is not None else None

        # 2. Step A: Interpolate (Raw Guess)
        # Note: Ideally we would fit (y / telluric), but standard practice 
        # often fits y directly (assuming tellurics are masked).
        try:
            cont_raw = interpolate_continuum_mask(
                x=self._data.x.values,
                y=self._data.y.values,
                mask_abs=mask_abs
            )
        except ValueError as e:
            logging.error(e)
            return self

        # 3. Step B: Auto-Fudge
        applied_fudge = 1.0
        if isinstance(fudge, str) and fudge.lower() == 'auto':
            # Calculate fudge using the raw guess
            applied_fudge = compute_auto_fudge(
                y=self._data.y.values,
                mask_abs=mask_abs,
                cont_model=cont_raw
            )
            logging.info(f"Auto-fudge calculated: {applied_fudge:.4f}")
        else:
            try:
                applied_fudge = float(fudge)
            except (ValueError, TypeError):
                logging.warning(f"Invalid fudge '{fudge}', defaulting to 1.0")

        cont_fudged = cont_raw * applied_fudge
        
        # 4. Step C: Smooth
        if smooth_std_kms > 0:
            cont_intrinsic = smooth_spectrum(
                x=self._data.x.quantity,
                y=cont_fudged,
                sigma_kms=smooth_std_kms,
                z_ref=self._data.z_em
            )
        else:
            cont_intrinsic = cont_fudged
            
        # 5. Build New Auxiliary Columns
        new_aux_cols = deepcopy(self._data.aux_cols)
        
        # [CRITICAL] Store the Intrinsic Continuum
        # Since we added the .norm property, 'cont' should remain clean of tellurics
        # so that .norm = cont * telluric works correctly.
        new_aux_cols['cont'] = DataColumnV2(
            values=cont_intrinsic,
            unit=self._data.y.unit, 
            description=f"Continuum (smooth={smooth_std_kms} km/s, fudge={applied_fudge:.3f})"
        )
        
        # 6. Step D: Renormalize Model (Corrected Logic)
        old_cont_col = self._data.aux_cols.get('cont')
        old_model_col = self._data.aux_cols.get('model')

        if renorm_model and old_cont_col is not None and old_model_col is not None:
            # Construct the FULL Normalization vectors (Cont * Telluric)
            # This ensures we preserve the optical depth (e^-tau) correctly.
            
            old_norm_vals = old_cont_col.values
            new_norm_vals = cont_intrinsic
            
            if telluric_vals is not None:
                # Apply telluric to both old and new baselines
                old_norm_vals = old_norm_vals * telluric_vals
                new_norm_vals = new_norm_vals * telluric_vals

            # Pass the *Effective Norms* to the helper
            new_model_col = self._renormalize_model(
                old_cont_vals=old_norm_vals,
                new_cont_vals=new_norm_vals,
                old_model_col=old_model_col
            )
            if new_model_col:
                new_aux_cols['model'] = new_model_col
                logging.info("Model column re-normalized to new continuum (incorporating tellurics).")

        # 7. Create New Data Core
        new_data = dataclasses.replace(self._data, aux_cols=new_aux_cols)
        
        # 8. Return New SpectrumV2
        hist_entry = f"Fitted continuum (fudge={applied_fudge:.3f}, smooth={smooth_std_kms} km/s)"
        return SpectrumV2(data=new_data, history=self.history + [hist_entry])
    

    def fit_powerlaw(self, regions_str: str, kappa: float = 3.0) -> 'SpectrumV2':
        """
        Fits a power-law continuum to specified rest-frame regions.

        Parameters
        ----------
        regions_str : str
            Comma-separated list of regions (e.g. ``'128.0-129.0, 131.5-132.5'``).
        kappa : float
            Sigma-clipping threshold for removing outliers within regions.

        Returns
        -------
        SpectrumV2
            A new instance containing the 'cont_pl' column.
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
        Detects absorption regions based on a boolean mask.

        Groups adjacent pixels marked as 'absorbed' into discrete regions.

        Parameters
        ----------
        mask_col : str
            Name of the mask column (True=Absorbed).
        min_pix : int
            Minimum number of pixels to form a valid region.

        Returns
        -------
        SpectrumV2
            A new instance with the 'region_id' column updated.
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
                       mask_col: str = 'abs_mask', 
                       min_pix_region: int = 3,
                       merge_dv: float = 10.0,
                       score_threshold: float = 0.5,
                       bypass_scoring: bool = False,
                       debug_rating: bool = False) -> 'SpectrumV2':
        """
        Identifies absorption lines by correlating regions with multiplet templates.

        This method employs a multi-step algorithm:
        1.  Detects and merges absorption regions in velocity space.
        2.  Scans for kinematic doublets (e.g. CIV, SiIV) redward of Ly-alpha.
        3.  Scans for Ly-alpha/Ly-beta pairs in the forest.
        4.  Scores candidates using an R^2 correlation metric.
        5.  Populates metadata with candidate identifications.

        Parameters
        ----------
        multiplet_list : list of str
            List of multiplets to search for (e.g. ``['CIV', 'SiIV']``).
        mask_col : str
            Name of the absorption mask column.
        min_pix_region : int
            Minimum pixels to consider a region.
        merge_dv : float
            Maximum velocity separation (km/s) to merge adjacent regions.
        score_threshold : float
            Minimum R^2 score (0-1) to accept a candidate.
        bypass_scoring : bool
            If True, accept all kinematic matches regardless of score.
        debug_rating : bool
            If True, produce debug plots for the rating process.

        Returns
        -------
        SpectrumV2
            A new instance with 'abs_ids' column and identification metadata.
        """
        
        # This is the 'self' spectrum object
        spec = self
        z_em = spec._data.z_em # <<< *** GET z_em ***
        if z_em == 0.0:
            raise ValueError("z_em must be set before running identify_lines.")
        
        logging.info(f"identify_lines (New Algorithm): Starting identification...")
        logging.debug(f"  Params: merge_dv={merge_dv}, score_threshold={score_threshold}")

        # 1. Silently create abs_mask if it doesn't exist
        if not spec.has_aux_column(mask_col):
            logging.info(f"Mask '{mask_col}' not found. Running 'find_absorbed' with defaults.")
            try:
                spec = spec.find_absorbed(
                    smooth_len_lya=5000.0 * au.km/au.s,
                    smooth_len_out=400.0 * au.km/au.s,
                    kappa=2.0,
                    template=False
                )
            except Exception as e:
                logging.error(f"Failed to auto-generate absorbed mask: {e}", exc_info=True)
                raise e
        
        # --- *** START: NEW ALGORITHM (Steps 1-5) *** ---
        
        # 1. Get Merged Regions
        try:
            mask_data = spec.get_column(mask_col).value
            if mask_data.dtype.kind != 'b':
                 mask_data = mask_data.astype(bool)
            
            region_id_map_raw, num_raw = detect_regions(mask_data, min_pix_region)
            logging.info(f"identify_lines: Found {num_raw} raw regions.")

            region_id_map_merged, num_merged = merge_regions_by_velocity(
                region_id_map_raw,
                spec.x,
                spec._data.z_em,
                merge_dv
            )
            logging.info(f"identify_lines: {num_merged} regions remain after merging.")
        except Exception as e:
            logging.error(f"Failed to detect or merge absorption regions: {e}", exc_info=True)
            raise e

        # Define AOD profile (for Step 5)
        cont = spec.cont.value if spec.cont is not None else np.ones_like(spec.y.value)
        flux = spec.y.value
        with np.errstate(divide='ignore', invalid='ignore'):
            aod = -np.log(flux / cont)
        aod[~np.isfinite(aod)] = 0.0 # Clean NaNs/Infs

        # Define boundaries
        x_nm = spec.x.to_value(au.nm)
        lya_limit_obs_nm = (1.0 + z_em) * xem_d['Ly_lim'].to_value(au.nm)
        lya_obs_nm = (1.0 + z_em) * xem_d['Ly_a'].to_value(au.nm)
        
        # Split region map into Forest and Redward
        forest_mask = (x_nm < lya_obs_nm) & (x_nm > lya_limit_obs_nm)
        redward_mask = (x_nm >= lya_obs_nm)
        
        forest_map = region_id_map_merged.copy()
        forest_map[~forest_mask] = 0 # Zero out non-forest regions
        
        redward_map = region_id_map_merged.copy()
        redward_map[~redward_mask] = 0 # Zero out non-redward regions

        # ---
        # 2. Kinematic check for Metal Doublets (Red Side)
        # ---
        all_candidates = [] # List of (rid, z_test, series_name)
        metal_doublets = [
            m for m in multiplet_list 
            if m in STANDARD_MULTIPLETS and len(STANDARD_MULTIPLETS[m]) == 2 and m != 'Ly_ab'
        ]
        logging.info(f"Checking for metal doublets: {metal_doublets}")
        
        for series_name in metal_doublets:
            # --- *** MODIFIED: Only check redward_map *** ---
            candidates_red = _find_kinematic_doublet_candidates(
                redward_map, spec.x, spec._data.z_em, series_name, z_em
            )
            for (rid_1, rid_2, z_test) in candidates_red:
                # --- *** Append the threshold *** ---
                all_candidates.append((rid_1, rid_2, z_test, series_name, score_threshold))
        
        # ---
        # 3. Kinematic check for Ly_ab (Forest)
        # ---
        if 'Ly_ab' in multiplet_list:
            logging.info("Checking for Ly_ab doublet...")
            lyab_candidates = _find_kinematic_doublet_candidates(
                forest_map, spec.x, spec._data.z_ref, 'Ly_ab', z_em
            )
            for (rid_1, rid_2, z_test) in lyab_candidates:
                # --- *** Append the threshold *** ---
                all_candidates.append((rid_1, rid_2, z_test, 'Ly_ab', score_threshold))

        logging.info(f"Found {len(all_candidates)} total kinematic candidates.")

        # ---
        # 4. Rate Candidates with R^2 Score
        # ---
        reliable_ids = {} # {rid: (series_name, score, z_test)}
        logging.info(f"Rating {len(all_candidates)} kinematic candidates with R^2 score...")

        for (rid_1, rid_2, z_test, series_name, threshold) in all_candidates:
            
            # --- *** Get the masks for the PAIR *** ---
            if series_name == 'Ly_ab':
                threshold = score_threshold # Use normal threshold for Ly_ab
            else:
                # It's a metal line. Check if it's in the forest.
                region_mask_for_check = (region_id_map_merged == rid_1)
                x_region_mean = np.mean(x_nm[region_mask_for_check])
                if x_region_mean < lya_obs_nm:
                    threshold = 1.1 # Set threshold to 1.1 (impossible to pass)
                else:
                    threshold = score_threshold # Use normal threshold
                
            mask_1 = (region_id_map_merged == rid_1)
            mask_2 = (region_id_map_merged == rid_2)
            
            if not np.any(mask_1) or not np.any(mask_2):
                continue # Skip if a region is empty (shouldn't happen)
            
            # --- *** Call the R^2 rater on the PAIR *** ---
            r2_score = rate_doublet_candidate(
                spec, mask_1, mask_2, series_name, z_test,
                debug_rating=debug_rating
            )
            
            if r2_score > threshold or bypass_scoring:
                log_msg = f"PASSED > {threshold}"
                if bypass_scoring and r2_score <= threshold:
                    log_msg = f"BYPASSED (Score: {r2_score:.3f})"
                
                logging.debug(f"  > Pair ({rid_1}, {rid_2}): {series_name} at z={z_test:.4f} scored R^2 = {r2_score:.3f} ({log_msg})")

                # Get the component names
                lines = STANDARD_MULTIPLETS[series_name]
                comp_1, comp_2 = lines[0], lines[1]
                rid_1_int = int(rid_1) # rid_1 is the red region
                rid_2_int = int(rid_2) # rid_2 is the blue region
                
                # Add Red Line (comp_2) to Red Region (rid_1)
                if rid_1_int not in reliable_ids: reliable_ids[rid_1_int] = []
                reliable_ids[rid_1_int].append((comp_2, r2_score, z_test))

                # Add Blue Line (comp_1) to Blue Region (rid_2)
                if rid_2_int not in reliable_ids: reliable_ids[rid_2_int] = []
                reliable_ids[rid_2_int].append((comp_1, r2_score, z_test))
            
            # (If score < threshold, we just drop the candidate)

        # ---
        # 5. Forest Cleanup (Label un-ID'd regions as Ly_a)
        # ---
        forest_rids = np.unique(forest_map)
        forest_rids = forest_rids[forest_rids > 0]
        
        for rid in forest_rids:
            rid_int = int(rid)
            if rid_int not in reliable_ids:
                # This region is in the forest but was not ID'd as Ly_ab
                # Label it as Ly_a
                # We need a z_test. We'll use the AOD peak.
                region_mask = (region_id_map_merged == rid_int)
                aod_region = aod.copy()
                aod_region[~region_mask] = -np.inf
                
                aod_peak_idx = np.argmax(aod_region)
                if aod_peak_idx > 0:
                    x_peak_nm = x_nm[aod_peak_idx]
                    z_test = (x_peak_nm / xem_d['Ly_a'].to_value(au.nm)) - 1.0
                    reliable_ids[rid_int] = [('Ly_a', 0.0, z_test)] # Give it 0 score
                    logging.debug(f"  > Region {rid_int}: Labeled as Ly_a (cleanup) at z={z_test:.4f}")

        logging.info(f"Finalized {len(reliable_ids)} identified regions.")
        
        # ---
        # 6. Save results for populating
        # ---
        final_vis_map = np.zeros_like(region_id_map_merged)
        series_map_for_populating = {} # {rid: series_name}
        z_map_for_populating = {}      # {rid: z_test}
        
        identifications_for_meta = {} # {rid: [(series, score)]}

        comp_id_counter = 1
        for rid, candidates in reliable_ids.items():
            final_vis_map[region_id_map_merged == rid] = rid
            
            # --- *** MODIFIED: Flatten list for component creation *** ---
            meta_cands = []
            for (component_name, score, z_test) in candidates:
                # Use a *new* unique ID for each component, not the region ID
                series_map_for_populating[comp_id_counter] = component_name
                z_map_for_populating[comp_id_counter] = z_test
                meta_cands.append((component_name, score))
                comp_id_counter += 1
            
            identifications_for_meta[rid] = meta_cands

        # ---
        # 7. Create Final Data Core (saving maps for Step 6)
        # ---
        new_meta = spec.meta
        new_aux_cols = deepcopy(spec._data.aux_cols)

        new_meta['num_regions_raw'] = int(num_raw)
        new_meta['num_regions_merged'] = int(num_merged)
        new_meta['num_regions_identified'] = len(reliable_ids)

        # Add the 'abs_ids' column
        vis_col = DataColumnV2(
            values=final_vis_map,
            unit=au.dimensionless_unscaled,
            description="Absorption regions with identifications"
        )
        new_aux_cols['abs_ids'] = vis_col

        # Save metadata
        try:
            new_meta['region_identifications'] = json.dumps(identifications_for_meta)
        except Exception as e:
            logging.error(f"Failed to serialize region identifications: {e}")
            new_meta['region_identifications'] = None

        # --- *** Store maps for Step 6 (Populating) *** ---
        try:
            new_meta['series_map_json'] = json.dumps(series_map_for_populating)
            new_meta['z_map_json'] = json.dumps(z_map_for_populating)
        except Exception as e:
             logging.error(f"Failed to serialize component maps: {e}")
        # --- *** END *** ---
        
        final_data_core = dataclasses.replace(
            spec._data, 
            meta=new_meta,
            aux_cols=new_aux_cols
        )
        
        # 9. Return new SpectrumV2
        new_history = self.history + [f"Identified {len(reliable_ids)} regions"]
        return SpectrumV2(data=final_data_core, history=new_history)
    
    def get_velocity_bounds(self, z_target: float, trans: str) -> tuple:
        """
        Calculates the min/max velocity coverage of this spectrum relative to a target.
        """
        if trans not in ATOM_DATA: return -np.inf, np.inf
        
        x = self.x.value # Assuming nm or Angstrom based on loader
        if len(x) == 0: return 0.0, 0.0
        
        # Determine rest wavelength in same unit as x (assuming nm for internal)
        lam_rest = ATOM_DATA[trans]['wave'] / 10.0 
        lam_obs = lam_rest * (1 + z_target)
        c_kms = 299792.458
        
        v = c_kms * (x - lam_obs) / lam_obs
        return np.min(v), np.max(v)

    def detect_doublet_z(self, other_spec: 'SpectrumV2', 
                         trans_self: str, trans_other: str) -> float:
        """
        Finds the redshift of maximum coincidence (product of absorption depths)
        between this spectrum and another.
        """
        x1, y1 = self.x.value, self.y.value
        x2, y2 = other_spec.x.value, other_spec.y.value
        
        if len(x1) == 0 or len(x2) == 0: return 0.0

        lam1 = ATOM_DATA[trans_self]['wave'] / 10.0
        lam2 = ATOM_DATA[trans_other]['wave'] / 10.0

        z1 = x1 / lam1 - 1.0
        z2 = x2 / lam2 - 1.0
        
        z_min, z_max = max(z1.min(), z2.min()), min(z1.max(), z2.max())
        if z_max <= z_min: return 0.0
        
        # Grid based on self resolution
        dz = np.nanmedian(np.diff(z1))
        if dz <= 0: dz = 1e-5
        z_grid = np.arange(z_min, z_max, dz)
        
        y1_int = np.interp(z_grid, z1, y1)
        y2_int = np.interp(z_grid, z2, y2)
        
        # Product of Depths
        depth1 = np.maximum(0.0, 1.0 - gaussian_filter1d(y1_int, 2.0))
        depth2 = np.maximum(0.0, 1.0 - gaussian_filter1d(y2_int, 2.0))
        
        score = depth1 * depth2
        best_idx = np.argmax(score)
        
        return float(z_grid[best_idx])
    
    def x_convert(self, z_ref: float, xunit: au.Unit) -> 'SpectrumV2':
        """
        API: Converts the X-axis units and returns a NEW SpectrumV2 instance.
        (Replaces V1's mutable _x_convert)
        """
        
        # 1. Execute pure conversion logic
        x_new, xmin_new, xmax_new = convert_x_axis(self._data, z_ref, xunit)
        
        # 2. Build new SpectrumDataV2, preserving the core vertical data
        new_data = SpectrumDataV2(
            x=DataColumnV2(x_new.value, xunit),
            xmin=DataColumnV2(xmin_new.value, xunit),
            xmax=DataColumnV2(xmax_new.value, xunit),
            y=self._data.y, 
            dy=self._data.dy,
            aux_cols=self._data.aux_cols,
            meta=self._data.meta,
            z_em=self._data.z_em   # Propagate
        )
        
        # 3. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"X-axis converted to {xunit} at z_ref={z_ref}"]
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
            z_em=self._data.z_em   # Propagate
        )
        
        # 3. Return a NEW SpectrumV2 instance
        new_history = self.history + [f"Y-axis converted to {yunit}"]
        return SpectrumV2(data=new_data, history=new_history)