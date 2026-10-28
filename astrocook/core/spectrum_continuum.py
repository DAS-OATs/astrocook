"""Continuum methods of :class:`~astrocook.core.spectrum.SpectrumV2`.

Part of the split of ``spectrum.py`` into thematic modules (phase 5): the
methods are defined here, in a mixin class that ``SpectrumV2`` inherits, and
are used as before (``spec.fit_continuum(...)``).
"""
import astropy.units as au
from copy import deepcopy
import dataclasses
import logging
import numpy as np
from scipy.interpolate import PchipInterpolator
from typing import TYPE_CHECKING, List, Union

from astrocook.core.spectrum_operations import find_absorbed_regions, fit_powerlaw_to_regions
from astrocook.core.structures import DataColumnV2

if TYPE_CHECKING:
    from astrocook.core.spectrum import SpectrumV2


class ContinuumMixin:
    """Estimate of the continuum: absorbed regions, interpolated or
    power-law continuum, continuum from manual knots.

    Mixed into ``SpectrumV2``, whose data and methods it uses: ``_data``,
    ``history``, ``x``, ``get_column``, ``has_aux_column``,
    ``update_column`` and ``_renormalize_model``. The results are new
    instances of the class (``self.__class__``).
    """

    def find_absorbed(self, smooth_len_lya: au.Quantity, smooth_len_out: au.Quantity, 
                      kappa: float, template: bool) -> 'SpectrumV2':
        """
        Identify absorbed regions using V1 logic.

        Parameters
        ----------
        smooth_len_lya : Quantity
            Smoothing length for the Ly-alpha forest.
        smooth_len_out : Quantity
            Smoothing length outside the forest.
        kappa : float
            Sigma clipping threshold.
        template : bool
            Use PCA template (Not Implemented).

        Returns
        -------
        SpectrumV2
            A new instance with the ``abs_mask`` column added/updated.
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
        return self.__class__(data=new_data, history=self.history + [f"Found absorbed (kappa={kappa})"])
        
    def fit_continuum(self, fudge: Union[float, str], smooth_std_kms: float, 
                      mask_col: str = 'abs_mask', renorm_model: bool = False,
                      template: bool = False) -> 'SpectrumV2':
        """
        Fit a continuum to the unabsorbed regions.

        Orchestrates: Interpolate -> Auto-Fudge -> Smooth -> Renormalize Model.

        Parameters
        ----------
        fudge : float or str
            Fudge factor. Can be a float or ``'auto'``. With ``'auto'`` the
            factor is computed from the residuals of the unabsorbed pixels
            around the interpolated continuum, which passes through those
            pixels: it is 1 by construction. In V1 the automatic fudge
            corrected the running mean from which V1 builds the continuum;
            V2 interpolates instead, so there is nothing to correct
            (findings register F3, accepted). The bias of the final
            smoothing is not corrected.
        smooth_std_kms : float
            Smoothing width for the continuum spline in km/s.
        mask_col : str, optional
            Name of the mask column to use. Defaults to ``'abs_mask'``.
        renorm_model : bool, optional
            If True, re-normalizes the 'model' column to the new continuum.
        template : bool, optional
            If True, uses a QSO template to trace Lyman-alpha emission.

        Returns
        -------
        SpectrumV2
            A new instance with the ``cont`` column added/updated.
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

        y_fit = self._data.y.values.copy()
        y_interp = None
        
        # --- NEW: Apply QSO Template Normalization ---
        if template:
            from astrocook.core.templates import qso_composite
            from scipy.interpolate import interp1d

            # The composite is in Angstroms. Convert to nm and shift to observed frame.
            x_template, y_template = qso_composite()
            x_template_nm = (x_template / 10.0) * (1.0 + self._data.z_em)

            # Interpolate template onto the current spectrum's grid (in nm)
            interp_fn = interp1d(x_template_nm, y_template, bounds_error=False, fill_value=1.0)
            y_interp = interp_fn(self.x.to_value(au.nm))

            # Normalize the working flux to flatten the emission peak
            valid = y_interp > 0
            y_fit[valid] = y_fit[valid] / y_interp[valid]
            logging.info("Applied QSO template for continuum fitting.")
        # ---------------------------------------------

        # 2. Step A: Interpolate (Raw Guess)
        try:
            cont_raw = interpolate_continuum_mask(
                x=self._data.x.values,
                y=y_fit, # Use flattened flux
                mask_abs=mask_abs
            )
        except ValueError as e:
            logging.error(e)
            return self

        # 3. Step B: Auto-Fudge
        applied_fudge = 1.0
        if isinstance(fudge, str) and fudge.lower() == 'auto':
            applied_fudge = compute_auto_fudge(
                y=y_fit, # Use flattened flux
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
        # We smooth the FLATTENED continuum to avoid smearing the template's Ly-a peak
        if smooth_std_kms > 0:
            cont_intrinsic = smooth_spectrum(
                x=self._data.x.quantity,
                y=cont_fudged,
                sigma_kms=smooth_std_kms,
                z_ref=self._data.z_em
            )
        else:
            cont_intrinsic = cont_fudged
            
        # --- NEW: Restore the Emission Peak ---
        if template and y_interp is not None:
            cont_intrinsic = cont_intrinsic * y_interp
        # --------------------------------------
        
        # 5. Build New Auxiliary Columns
        new_aux_cols = deepcopy(self._data.aux_cols)
        
        # [CRITICAL] Store the Intrinsic Continuum
        new_aux_cols['cont'] = DataColumnV2(
            values=cont_intrinsic,
            unit=self._data.y.unit, 
            description=f"Continuum (smooth={smooth_std_kms} km/s, fudge={applied_fudge:.3f}, template={template})"
        )
        
        # 6. Step D: Renormalize Model
        old_cont_col = self._data.aux_cols.get('cont')
        old_model_col = self._data.aux_cols.get('model')

        if renorm_model and old_cont_col is not None and old_model_col is not None:
            old_norm_vals = old_cont_col.values
            new_norm_vals = cont_intrinsic
            
            if telluric_vals is not None:
                old_norm_vals = old_norm_vals * telluric_vals
                new_norm_vals = new_norm_vals * telluric_vals

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
        hist_entry = f"Fitted continuum (fudge={applied_fudge:.3f}, smooth={smooth_std_kms} km/s, template={template})"
        return self.__class__(data=new_data, history=self.history + [hist_entry])
    

    def fit_powerlaw(self, regions_str: str, kappa: float = 3.0) -> 'SpectrumV2':
        """
        Fit a power-law continuum to specified rest-frame regions.

        Parameters
        ----------
        regions_str : str
            Comma-separated list of regions (e.g. ``'128.0-129.0, 131.5-132.5'``).
        kappa : float, optional
            Sigma-clipping threshold. Defaults to ``3.0``.

        Returns
        -------
        SpectrumV2
            A new instance containing the ``cont_pl`` column.
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
        return self.__class__(data=new_data, history=new_history)
    
    def update_continuum_from_knots(self, knots_x: List[float], knots_y: List[float], 
                                    renorm_model: bool = True) -> 'SpectrumV2':
        """
        Update the continuum using a Cubic Spline through provided knots.

        Parameters
        ----------
        knots_x : list of float
            X-coordinates of the knots.
        knots_y : list of float
            Y-coordinates of the knots.
        renorm_model : bool, optional
            If True, re-normalizes the 'model' column.

        Returns
        -------
        SpectrumV2
            A new instance with the updated continuum.
        """
        # 1. Validation (Safe for numpy arrays)
        # Check for None first, then length. Avoid "if not knots_x" on arrays.
        if knots_x is None or len(knots_x) == 0:
            logging.warning("Invalid knots provided (empty).")
            return self
            
        if len(knots_x) != len(knots_y):
            logging.warning(f"Knots length mismatch: x({len(knots_x)}) != y({len(knots_y)})")
            return self

        # 2. Create Spline
        try:
            # Sort knots to prevent spline errors
            knots_x_arr = np.array(knots_x)
            knots_y_arr = np.array(knots_y)
            sorter = np.argsort(knots_x_arr)
            
            cs = PchipInterpolator(knots_x_arr[sorter], knots_y_arr[sorter])
            
            # Evaluate on current grid
            x_axis = self.x.value
            new_cont = cs(x_axis)
        except Exception as e:
            logging.error(f"Failed to generate spline from knots: {e}")
            return self

        # 3. Update 'cont' column
        new_spec = self.update_column('cont', new_cont, description="Spline from manual knots")
        
        # 4. Handle Model Renormalization
        if renorm_model and new_spec.has_aux_column('model'):
            old_cont_col = self.get_column('cont')
            
            # We can use the helper we used in smooth_column/fit_continuum
            # Note: We access _data directly inside the class for the helper
            new_model_data_col = self._renormalize_model(
                old_cont_vals=old_cont_col.value,
                new_cont_vals=new_cont,
                old_model_col=self._data.aux_cols['model']
            )
            
            if new_model_data_col:
                # We use internal update logic to avoid creating two SpectrumV2 objects
                new_aux = deepcopy(new_spec._data.aux_cols)
                new_aux['model'] = new_model_data_col
                new_data = dataclasses.replace(new_spec._data, aux_cols=new_aux)
                new_spec = self.__class__(new_data, new_spec.history)
                logging.info("Model renormalized to match new continuum.")

        # 5. History
        new_spec.history.append("Updated continuum from knots")
        return new_spec
