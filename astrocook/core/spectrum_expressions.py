"""Expressions, masks and regions of :class:`~astrocook.core.spectrum.SpectrumV2`:
columns computed or masked with NumExpr expressions (also across spectra),
splitting by an expression, masks limited to a range, preset regions.

Part of the split of ``spectrum.py`` into thematic modules (phase 5): the
methods are defined here, in a mixin class that ``SpectrumV2`` inherits, and
are used as before (``spec.apply_expression(...)``).
"""
import astropy.units as au
from copy import deepcopy
import logging
import numexpr as ne
import numpy as np
import re
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.structures import DataColumnV2, SpectrumDataV2

if TYPE_CHECKING:
    from astrocook.core.spectrum import SpectrumV2

_NE_GLOBAL_DICT = {
    'min': np.nanmin,
    'max': np.nanmax,
    'median': np.nanmedian,
    'mean': np.nanmean,
    'std': np.nanstd,
    # numexpr *does* have: sum, prod, mean, std, min, max, log, log10, etc.
}

REGION_PRESETS = {
    'lya_forest': (102.57 * au.nm, 121.567 * au.nm),   # Ly-beta to Ly-alpha
    'all_ly_forest': (91.18 * au.nm, 121.567 * au.nm), # Ly-limit to Ly-alpha
    'red_side': (121.567 * au.nm, 10000.0 * au.nm)     # Redward of Ly-alpha
}


class ExpressionsMixin:
    """Expressions on the columns, masks and regions.

    Mixed into ``SpectrumV2``, whose data and methods it uses: ``_data``,
    ``history``, ``x``, ``t``, ``get_column``, ``update_column`` and
    ``get_resampled_column`` (of ``CombineMixin``). The results are new
    instances of the class (``self.__class__``).
    """

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

    def get_region_bounds(self, region: str) -> Tuple[au.Quantity, au.Quantity]:
        """
        Get observed wavelength bounds for a standard spectral region.

        Calculates the observed wavelength range for a preset region (e.g.,
        Ly-alpha forest) based on the object's emission redshift.

        Parameters
        ----------
        region : str
            The region identifier. Supported presets:
            - ``'lya_forest'``: Ly-beta to Ly-alpha.
            - ``'all_ly_forest'``: Ly-limit to Ly-alpha.
            - ``'red_side'``: Redward of Ly-alpha.

        Returns
        -------
        tuple
            A tuple ``(obs_min, obs_max)`` as Astropy Quantities.

        Raises
        ------
        ValueError
            If ``z_em`` is 0.0 or the region preset is unknown.
        """
        z_em = self._data.z_em
        if z_em == 0.0:
            raise ValueError("z_em is 0.0. Cannot calculate rest-frame regions.")
            
        if region not in REGION_PRESETS:
            raise ValueError(f"Unknown region '{region}'. Options: {list(REGION_PRESETS.keys())}")
            
        rf_min, rf_max = REGION_PRESETS[region]
        
        obs_min = rf_min * (1.0 + z_em)
        obs_max = rf_max * (1.0 + z_em)
        
        return obs_min, obs_max

    # --- Methods implementing the API ---

    def prepare_expression_context(self, expression: str, 
                                   other_specs: Dict[str, 'SpectrumV2']) -> Tuple[str, Dict[str, np.ndarray]]:
        """
        Prepare context for cross-spectrum arithmetic expressions.

        Parses an expression for aliases (e.g., ``"s1.F"``), resolves them against
        provided spectrum objects, and resamples their data onto the current
        spectrum's grid.

        Parameters
        ----------
        expression : str
            The raw expression string (e.g., ``"F - s1.F"``).
        other_specs : dict
            A dictionary mapping alias names to :class:`SpectrumV2` objects.

        Returns
        -------
        tuple
            A tuple containing:
            1. **final_expr** (*str*): The processed expression with safe variable names.
            2. **extra_vars** (*dict*): A dictionary mapping variable names to resampled numpy arrays.
        """
        # Regex for "alias.column"
        regex = re.compile(r'([a-zA-Z0-9_]+)\s*\.\s*([a-zA-Z0-9_]+)')
        extra_vars = {}
        
        # We need the current grid in nm for resampling
        current_grid_nm = self.x.to_value(au.nm)
        
        def replacer(match):
            alias = match.group(1)
            col_name = match.group(2)
            
            if alias not in other_specs:
                return match.group(0) # Not a known alias, ignore
                
            safe_var_name = f"{alias}_{col_name}"
            if safe_var_name in extra_vars:
                return safe_var_name
            
            target_spec = other_specs[alias]
            
            # Resample target column onto current grid
            try:
                resampled = target_spec.get_resampled_column(col_name, current_grid_nm)
                if resampled is None:
                    raise ValueError(f"Column '{col_name}' in session '{alias}' cannot be used (non-numerical).")
                
                extra_vars[safe_var_name] = resampled
                return safe_var_name
            except Exception as e:
                # We raise here to stop the process if a variable is invalid
                raise ValueError(f"Failed to resample {alias}.{col_name}: {e}")

        final_expr = regex.sub(replacer, expression)
        return final_expr, extra_vars

    def apply_expression(self, target_col: str, expression: str, 
                         extra_vars: Dict[str, np.ndarray] = None) -> 'SpectrumV2': # <<< MODIFIED
        """
        Apply a numerical expression to columns and return a new spectrum.

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
        target_col = target_col.strip() # --- Enforce clean column names ---
        
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
        return self.__class__(data=new_data_core, history=new_history)

    def mask_expression(self, target_col: str, expression: str, 
                        extra_vars: Dict[str, np.ndarray] = None) -> 'SpectrumV2': # <<< MODIFIED
        """
        Mask a column based on a boolean expression (sets values to NaN).

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
        target_col = target_col.strip() # --- Enforce clean column names ---
        
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
        return self.__class__(data=new_data_core, history=new_history)
    
    def calc_intersection_expression(self, others_specs: List['SpectrumV2'], 
                                     z_target: float, trans_self: str, 
                                     trans_others: List[str], window_kms: float) -> Optional[str]:
        """
        Calculate the wavelength expression for the velocity intersection of spectra.

        Determines the wavelength range where this spectrum and a list of others
        overlap in velocity space, relative to specific transitions.

        Parameters
        ----------
        others_specs : list of SpectrumV2
            List of other spectra to intersect with.
        z_target : float
            Target redshift for the velocity conversion.
        trans_self : str
            Transition name for the current spectrum (e.g., ``'Ly_a'``).
        trans_others : list of str
            List of transition names for the other spectra, corresponding to ``others_specs``.
        window_kms : float
            Symmetric velocity window width in km/s.

        Returns
        -------
        str or None
            A boolean expression string (e.g., ``"(x > 400.0) & (x < 500.0)"``) defining the
            intersection, or ``None`` if no overlap exists.
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

    def split(self, expression: str, extra_vars: Dict[str, np.ndarray] = None) -> 'SpectrumV2':
        """
        Split the spectrum based on a boolean expression.

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
        return self.__class__(data=new_data_core, history=new_history)
    
    def mask_region(self, mask_col: str, 
                    obs_min: au.Quantity, obs_max: au.Quantity) -> 'SpectrumV2':
        """
        Update a mask column, setting it to False outside the specified bounds.

        Parameters
        ----------
        mask_col : str
            Name of the mask column.
        obs_min : Quantity
            Minimum observed wavelength.
        obs_max : Quantity
            Maximum observed wavelength.

        Returns
        -------
        SpectrumV2
            A new instance with the updated mask.
        """
        col = self.get_column(mask_col)
        if col is None:
            logging.warning(f"Mask column '{mask_col}' not found.")
            return self

        # Work in Angstroms for consistency
        x_ang = self.x.to(au.Angstrom).value
        min_ang = obs_min.to(au.Angstrom).value
        max_ang = obs_max.to(au.Angstrom).value
        
        # Copy and Modify
        new_mask = col.value.astype(bool).copy()
        new_mask[x_ang < min_ang] = False
        new_mask[x_ang > max_ang] = False
        
        return self.update_column(mask_col, new_mask)
