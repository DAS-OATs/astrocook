
from astrocook._renamed import renamed
import astropy.units as au
from copy import deepcopy
import dataclasses
import logging
import numpy as np
from typing import Any, Dict, List, Optional, Union, overload

from astrocook.core.spectrum_operations import convert_x_axis, convert_y_axis
from astrocook.core.structures import SpectrumData, DataColumn
from astrocook.core.resolution import column_to_resolving_power, migrate_resolution
from astrocook.core.spectrum_combine import CombineMixin
from astrocook.core.spectrum_continuum import ContinuumMixin
from astrocook.core.spectrum_expressions import ExpressionsMixin
from astrocook.core.spectrum_flux import FluxMixin
from astrocook.core.spectrum_lines import LinesMixin


#: The names of the core columns in the GUI and in the expressions of its
#: dialogs, and the ones in the table (the MCP reads accept both, F120)
GUI_COLUMN_NAMES = {'λ': 'x', 'λmin': 'xmin', 'λmax': 'xmax', 'F': 'y', 'dF': 'dy'}


def _block_means(values: np.ndarray, n: int, error: bool = False) -> np.ndarray:
    """Means over blocks of ``n`` values (the last block may be shorter);
    an error column (``error``) as the error of the mean. NaN are skipped."""
    if n == 1:
        return values
    starts = np.arange(0, len(values), n)
    out = np.empty(len(starts))
    for k, start in enumerate(starts):
        block = values[start:start + n]
        good = block[np.isfinite(block)]
        if not good.size:
            out[k] = np.nan
        elif error:
            out[k] = np.sqrt(np.sum(good ** 2)) / good.size
        else:
            out[k] = np.mean(good)
    return out

logger = logging.getLogger(__name__)


class TableAdapter(object):
    """Adapter to wrap the V2 .t dictionary and expose Astropy Table-like properties."""
    def __init__(self, data_dict):
        self._data_dict = data_dict
        self.colnames = list(data_dict.keys())
        
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
    
    @overload
    def __getitem__(self, key: str) -> au.Quantity: ...
    @overload
    def __getitem__(self, key: slice) -> 'TableAdapter': ...
    def __getitem__(self, key: Union[str, slice]) -> Union[au.Quantity, 'TableAdapter']:
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
                    
            # Return a NEW TableAdapter instance containing the sliced data
            return TableAdapter(new_data_dict)
        
        elif isinstance(key, list):
            # Selection of multiple columns: return a new TableAdapter
            # with only the requested subset of the dictionary.
            new_data_dict = {col: self._data_dict[col] for col in key if col in self._data_dict}
            return TableAdapter(new_data_dict)
        
        else:
            raise KeyError(f"Invalid key type: {type(key)}. Expected string or slice.")

def _edges_with_unit(data: SpectrumData) -> SpectrumData:
    """Pixel edges (``xmin``, ``xmax``) without a unit, from sessions saved
    by earlier versions, take the unit of ``x`` (#623: the fit converts the
    edges to Angstrom since #526, and failed with "'' (dimensionless) and
    'Angstrom' (length) are not convertible"). If their values do not
    bracket ``x`` (at least 90% of the pixels), they are not edges in the
    unit of ``x`` and are recomputed halfway between the pixels."""
    x_unit = au.Unit(data.x.unit) if data.x.unit is not None else au.dimensionless_unscaled
    if x_unit == au.dimensionless_unscaled:
        return data

    def no_unit(col):
        return col.unit is None or au.Unit(col.unit) == au.dimensionless_unscaled

    if not (no_unit(data.xmin) or no_unit(data.xmax)):
        return data
    x = np.asarray(data.x.values, dtype=float)
    xmin = np.asarray(data.xmin.values, dtype=float)
    xmax = np.asarray(data.xmax.values, dtype=float)
    with np.errstate(invalid='ignore'):
        inside = (xmin <= x) & (x <= xmax)
    finite = np.isfinite(x) & np.isfinite(xmin) & np.isfinite(xmax)
    if finite.any() and np.mean(inside[finite]) >= 0.9:
        logger.info("Pixel edges without a unit: read in the unit of x (%s).", x_unit)
    else:
        logger.warning("Pixel edges without a unit and not around x: recomputed "
                        "halfway between the pixels.")
        if len(x) > 1:
            mid = 0.5 * (x[1:] + x[:-1])
            xmin = np.concatenate([[x[0] - (mid[0] - x[0])], mid])
            xmax = np.concatenate([mid, [x[-1] + (x[-1] - mid[-1])]])
        else:
            xmin, xmax = x.copy(), x.copy()
    return dataclasses.replace(
        data,
        xmin=DataColumn(xmin, x_unit, getattr(data.xmin, 'description', '')),
        xmax=DataColumn(xmax, x_unit, getattr(data.xmax, 'description', '')))


class Spectrum(CombineMixin, ContinuumMixin, ExpressionsMixin, FluxMixin, LinesMixin):
    """
    Composite container for spectral data and operations.

    This class wraps the raw data (flux, wavelength, error) and provides
    immutable operations that return new ``Spectrum`` instances. It serves as
    the primary data structure passed between Recipes and the GUI.

    Parameters
    ----------
    data : SpectrumData
        The core data container holding the Astropy Quantities.
    history : list, optional
        A list of strings describing the modification history of this spectrum.
    """
    
    def __init__(self, data: SpectrumData, history: Optional[list] = None):
        # The resolution lives in the 'resol' column only: legacy sources
        # (attribute, meta keys) are migrated into it here (#526)
        self._data = migrate_resolution(_edges_with_unit(data))
        self.history = history if history is not None else []

    def resolution(self) -> np.ndarray:
        """
        Resolving power R of each pixel, from the ``resol`` column.

        The single source of the resolution (#526). The unit of the column
        decides how it is read: dimensionless for R, a velocity for the FWHM.

        Returns
        -------
        np.ndarray
            R for each pixel; NaN where it is not known (all NaN if the
            spectrum has no resolution).
        """
        col = self._data.aux_cols.get('resol')
        if col is None:
            return np.full(len(self._data.x.values), np.nan)
        return column_to_resolving_power(col.values, col.unit)

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
    def t(self) -> TableAdapter:
        """
        Adapter: Replaces the raw table/dictionary access.
        Returns a TableAdapter instance to simulate an Astropy Table (V1)
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
        return TableAdapter(output)

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
        pass # L'unità X dovrebbe essere aggiornata solo creando un nuovo Spectrum
    
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
        """
        Check if an auxiliary column exists.

        Parameters
        ----------
        name : str
            The name of the column to check (e.g., ``'cont'``).

        Returns
        -------
        bool
            ``True`` if the column is present in the auxiliary data, ``False`` otherwise.
        """
        return name in self._data.aux_cols
    
    def columns_in_range(self, names: List[str], x_min: Optional[float] = None,
                         x_max: Optional[float] = None,
                         max_points: Optional[int] = None) -> Dict[str, Any]:
        """
        The values of some columns in a range of x, to read them (phase
        6-ter, PR 3: the MCP server, or a script).

        A range of more than ``max_points`` pixels is averaged in blocks of
        ``ceil(n / max_points)`` pixels: the columns whose name starts with
        ``d`` (errors) as the error of the mean, the others as the mean.
        This is for reading only: the spectrum is not rebinned (that is
        *Rebin*).

        Parameters
        ----------
        names : list of str
            The columns (``x`` is always given). The names of the GUI are
            accepted for the core ones (``λ``, ``λmin``, ``λmax``, ``F``,
            ``dF``: ``x``, ``xmin``, ``xmax``, ``y``, ``dy``) and the values
            are given under the name asked (F120).
        x_min, x_max : float, optional
            The range, in the unit of x; None: the start (end) of the
            spectrum.
        max_points : int, optional
            The largest number of points; None: every pixel.

        Returns
        -------
        dict
            ``columns`` (name -> ndarray), ``units`` (name -> str),
            ``n_pixels`` (pixels in the range), ``binned_by``.

        Raises
        ------
        KeyError
            If a column does not exist (the message lists those that do).
        ValueError
            If the range is empty or outside the spectrum.
        """
        table = self.t
        available = table.colnames
        internal = {n: GUI_COLUMN_NAMES.get(n, n) for n in names}
        missing = [n for n in names if internal[n] not in available]
        if missing:
            raise KeyError(f"No column {', '.join(missing)}. The columns: "
                           f"{', '.join(available)}.")
        x = self.x.value
        lo = np.nanmin(x) if x_min is None else x_min
        hi = np.nanmax(x) if x_max is None else x_max
        if lo >= hi:
            raise ValueError(f"Empty range: the start ({lo}) must be less than the "
                             f"end ({hi}).")
        sel = (x >= lo) & (x <= hi)
        n = int(sel.sum())
        if not n:
            raise ValueError(f"No pixel between {lo} and {hi} {self.x.unit}: the "
                             f"spectrum covers {np.nanmin(x):.4g}-{np.nanmax(x):.4g} "
                             f"{self.x.unit}.")
        binned_by = 1 if not max_points else max(1, int(np.ceil(n / max_points)))
        columns, units = {}, {}
        for name in ["x", *[c for c in names if c != "x"]]:
            col = internal.get(name, name)
            values = np.asarray(table[col].value, dtype=float)[sel]
            columns[name] = _block_means(values, binned_by, error=col.startswith("d"))
            units[name] = str(table[col].unit or "")
        return {"columns": columns, "units": units, "n_pixels": n,
                "binned_by": binned_by}

    def get_column(self, name: str) -> Optional[au.Quantity]:
        """
        Get an auxiliary column as an Astropy Quantity.

        Parameters
        ----------
        name : str
            The name of the auxiliary column to retrieve (e.g., ``'cont'``).

        Returns
        -------
        astropy.units.Quantity or None
            The column data with units, or ``None`` if the column does not exist.
        """
        col = self._data.aux_cols.get(name)
        return col.quantity if col else None

    def with_properties(self, z_em: float, resol: float = 0.0, 
                        object_name: Optional[str] = None,
                        overwrite_resol_col: bool = False) -> 'Spectrum':
        """
        Return a NEW Spectrum instance with updated properties.

        Parameters
        ----------
        z_em : float
            The emission redshift.
        resol : float, optional
            Resolving power R; ``0.0`` (the default) leaves the resolution
            unchanged.
        object_name : str, optional
            The object name (header keyword).
        overwrite_resol_col : bool, optional
            If ``True``, a ``resol`` column that varies with wavelength is
            replaced by a flat one at ``resol``. If ``False`` (the default,
            #526) it is kept, and ``resol`` only replaces a flat column or
            creates a missing one.

        Returns
        -------
        Spectrum
            A new instance with updated metadata.
        """
        new_meta = deepcopy(self._data.meta)
        
        history_entries = [f"z_em={z_em}"]

        # The resolution goes into the column only (below; #526)
        if resol > 0:
            history_entries.append(f"R={resol}")

        # Update Object Name
        if object_name is not None:
            new_meta['OBJECT'] = object_name
            # Also update standard ESO keyword if present, for consistency
            if 'HIERARCH ESO OBS TARG NAME' in new_meta:
                new_meta['HIERARCH ESO OBS TARG NAME'] = object_name
            history_entries.append(f"Object='{object_name}'")

        # Update Core Data
        new_data = dataclasses.replace(self._data, z_em=z_em, meta=new_meta)
        
        new_history = self.history + [f"Set properties ({', '.join(history_entries)})"]
        spec = Spectrum(data=new_data, history=new_history)

        if resol > 0:
            current = self.resolution()
            finite = current[np.isfinite(current)]
            variable = finite.size > 0 and not np.allclose(finite, finite[0], rtol=1e-6)
            if variable and not overwrite_resol_col:
                logger.warning("The spectrum has a resolution that varies with "
                                "wavelength: kept (R=%g not applied; set "
                                "overwrite_resol_col to replace it).", resol)
            else:
                spec = spec.update_column(
                    'resol', np.full(len(spec.x), float(resol)),
                    unit=au.dimensionless_unscaled,
                    description="Resolving power (set by the user)")

        return spec

    def _renormalize_model(self, 
                           old_cont_vals: np.ndarray, 
                           new_cont_vals: np.ndarray, 
                           old_model_col: DataColumn) -> Optional[DataColumn]:
        """
        Private helper to re-normalize an old model to a new continuum
        by applying the *old* normalized profile to the *new* continuum.
        
        Used by fit_continuum and smooth_column.
        """
        if old_cont_vals is None or new_cont_vals is None or old_model_col is None:
            return None
            
        logger.debug("Renormalizing the model to the new continuum.")
        # 1. Get old normalized profile
        normalized_model = np.divide(
            old_model_col.values,
            old_cont_vals,
            out=np.ones_like(old_model_col.values),
            where=old_cont_vals != 0
        )
        # 2. Recalculate model against *new* continuum
        new_model_values = new_cont_vals * normalized_model
        logger.info("Model renormalized to the new continuum.")
        return DataColumn(
            values=new_model_values,
            unit=old_model_col.unit, # Use old model's unit
            description=old_model_col.description
        )
    
    def sanitize_legacy_tellurics(self) -> 'Spectrum':
        """
        Patches legacy spectra where 'cont' included tellurics and 'cont_no_telluric' 
        stored the intrinsic continuum.
        
        Logic:
        
        1. If 'cont_no_telluric' exists:
           - Overwrite 'cont' with 'cont_no_telluric' (Intrinsic).
           - Remove 'cont_no_telluric'.
        2. Returns a new Spectrum instance (or self if no changes needed).
        """
        # Check if the legacy column exists
        if 'cont_no_telluric' not in self._data.aux_cols:
            return self
            
        logger.info("Old spectrum with 'cont_no_telluric': it becomes 'cont', so "
                    "that the telluric transmission is not counted twice.")
        
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
        return Spectrum(data=new_data, history=self.history + ["Sanitized legacy telluric columns"])
    
    def update_column(self, name: str, values: np.ndarray, 
                      unit: Optional[au.Unit] = None,
                      description: str = "") -> 'Spectrum':
        """
        Return a NEW Spectrum instance with the specified column updated or added.

        Handles both core columns (y, dy) and auxiliary columns (cont, model, etc.).

        Parameters
        ----------
        name : str
            Name of the column to update (e.g., ``'cont'``, ``'model'``).
        values : np.ndarray
            The numerical data for the column.
        unit : astropy.units.Unit, optional
            The physical unit. If None, inherits from existing column or defaults to flux unit.
        description : str, optional
            Description of the column content.

        Returns
        -------
        Spectrum
            A new instance with the updated column data.
        """
        name = name.strip() # --- Enforce clean column names ---

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
        
        # 2. Create the DataColumn object
        # (DataColumn is already imported in your spectrum.py)
        new_col = DataColumn(values, target_unit, description)
        
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
            
        # 4. Return new Spectrum
        new_history = self.history + [f"Updated column '{name}'"]
        return Spectrum(data=new_data, history=new_history)

    def remove_meta_keys(self, keys: List[str]) -> 'Spectrum':
        """
        Remove specified keys from the metadata dictionary.

        Returns a new Spectrum instance with the specified metadata keys removed.
        This is commonly used to clean up temporary large objects (like JSON
        mapping dictionaries) before saving the session to disk.

        Parameters
        ----------
        keys : list of str
            A list of metadata keys to remove.

        Returns
        -------
        Spectrum
            A new instance with the cleaned metadata.
        """
        new_meta = deepcopy(self._data.meta)
        for k in keys:
            if k in new_meta:
                del new_meta[k]
        
        new_data = dataclasses.replace(self._data, meta=new_meta)
        return Spectrum(data=new_data, history=self.history)

    def remove_columns(self, col_names: List[str]) -> 'Spectrum':
        """
        Return a NEW Spectrum with specified auxiliary columns removed.

        Parameters
        ----------
        col_names : list of str
            Names of the columns to delete.

        Returns
        -------
        Spectrum
            A new instance with the columns removed.
        """
        new_aux_cols = deepcopy(self._data.aux_cols)
        removed = []
        
        # Create a mapping of stripped names to actual internal keys.
        # This allows users to delete columns that were imported from ASCII 
        # or manually created with accidental trailing whitespace.
        key_map = {k.strip(): k for k in new_aux_cols.keys()}

        for col in col_names:
            clean_col = col.strip()
            if clean_col in key_map:
                actual_key = key_map[clean_col]
                del new_aux_cols[actual_key]
                removed.append(actual_key) # Log the exact key removed
            elif clean_col in ['x', 'y', 'dy', 'xmin', 'xmax']:
                logger.warning(f"The column '{clean_col}' is a core column and "
                               "cannot be removed.")
            else:
                logger.warning(f"There is no column '{clean_col}' to remove.")
        
        if not removed:
            return self

        new_data = dataclasses.replace(self._data, aux_cols=new_aux_cols)
        new_history = self.history + [f"Deleted columns: {removed}"]
        return Spectrum(data=new_data, history=new_history)

    def merge_column(self, col_name: str, other_values: np.ndarray) -> 'Spectrum':
        """
        Merge new values into an existing column.

        Strategy: Use new value if it's non-zero, otherwise keep old value.

        Parameters
        ----------
        col_name : str
            Name of the column.
        other_values : np.ndarray
            Array of new values.

        Returns
        -------
        Spectrum
            A new instance with merged data.
        """
        old_col = self.get_column(col_name)
        if old_col is None:
            # If column doesn't exist, just create it with the new values
            return self.update_column(col_name, other_values, unit=au.dimensionless_unscaled)
            
        old_values = old_col.value
        
        # Merge Logic: new > 0 ? new : old
        merged = np.where(other_values != 0, other_values, old_values)
        
        return self.update_column(col_name, merged)
    
    def update_model(self, model_flux: np.ndarray) -> 'Spectrum':
        """
        Update the 'model' auxiliary column with new flux values.

        Parameters
        ----------
        model_flux : np.ndarray
            The new model flux array.

        Returns
        -------
        Spectrum
            A new instance with the updated model.
        """
        # Create DataColumn (inheriting unit from flux)
        new_col = DataColumn(
            values=model_flux,
            unit=self._data.y.unit,
            description="Voigt Fit Model"
        )
        
        # Immutable Update
        new_aux = deepcopy(self._data.aux_cols)
        new_aux['model'] = new_col
        new_data = dataclasses.replace(self._data, aux_cols=new_aux)
        
        # We don't necessarily add a history entry for every model update 
        # as it happens very frequently during interactive fitting.
        return Spectrum(data=new_data, history=self.history)
    
    def x_convert(self, z_ref: float, xunit: au.Unit) -> 'Spectrum':
        """
        Convert the X-axis units and return a NEW Spectrum instance.

        Parameters
        ----------
        z_ref : float
            Reference redshift for velocity conversion.
        xunit : astropy.units.Unit
            Target unit (e.g., ``nm``, ``km/s``).

        Returns
        -------
        Spectrum
            A new instance with the converted X-axis.
        """
        
        # 1. Execute pure conversion logic
        x_new, xmin_new, xmax_new = convert_x_axis(self._data, z_ref, xunit)
        
        # 2. Build new SpectrumData, preserving the core vertical data
        new_data = SpectrumData(
            x=DataColumn(x_new.value, xunit),
            xmin=DataColumn(xmin_new.value, xunit),
            xmax=DataColumn(xmax_new.value, xunit),
            y=self._data.y, 
            dy=self._data.dy,
            aux_cols=self._data.aux_cols,
            meta=self._data.meta,
            z_em=self._data.z_em   # Propagate
        )
        
        # 3. Return a NEW Spectrum instance
        new_history = self.history + [f"X-axis converted to {xunit} at z_ref={z_ref}"]
        return Spectrum(data=new_data, history=new_history)

    def y_convert(self, yunit: au.Unit) -> 'Spectrum':
        """
        Convert the Y-axis units and return a NEW Spectrum instance.

        Parameters
        ----------
        yunit : astropy.units.Unit
            Target flux unit.

        Returns
        -------
        Spectrum
            A new instance with converted Y-axis and error.
        """
        
        # 1. Execute pure conversion logic
        new_y, new_dy, new_aux_cols = convert_y_axis(self._data, yunit)
        
        # 2. Build new SpectrumData, preserving the X-axis and metadata
        new_data = SpectrumData(
            x=self._data.x, 
            xmin=self._data.xmin,
            xmax=self._data.xmax,
            y=new_y, # Updated
            dy=new_dy, # Updated
            aux_cols=new_aux_cols, # Auxiliary columns should also be converted if unit-dependent
            meta=self._data.meta,
            z_em=self._data.z_em   # Propagate
        )
        
        # 3. Return a NEW Spectrum instance
        new_history = self.history + [f"Y-axis converted to {yunit}"]
        return Spectrum(data=new_data, history=new_history)

# Old names (phase 6-quater, PR 6): see astrocook/_renamed.py
__getattr__ = renamed(__name__, {"TableAdapterV2": "TableAdapter", "SpectrumV2": "Spectrum"})
