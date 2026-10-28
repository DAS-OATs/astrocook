"""Methods of :class:`~astrocook.core.spectrum.SpectrumV2` that combine or
regrid spectra: stitch, rebin, equalization, co-addition, resampling.

Part of the split of ``spectrum.py`` into thematic modules (phase 5): the
methods are defined here, in a mixin class that ``SpectrumV2`` inherits, and
are used as before (``spec.rebin(...)``).
"""
import astropy.units as au
from copy import deepcopy
import logging
import numpy as np
from typing import TYPE_CHECKING, List, Optional

from astrocook.core.spectrum_operations import convert_axis_velocity, rebin_spectrum
from astrocook.core.structures import DataColumnV2, SpectrumDataV2

if TYPE_CHECKING:
    from astrocook.core.spectrum import SpectrumV2


class CombineMixin:
    """Stitching, rebinning, equalization, co-addition and resampling of
    spectra.

    Mixed into ``SpectrumV2``, whose data and methods it uses: ``_data``,
    ``history``, ``meta``, ``x``, ``y``, ``dy``, ``t``, ``resolution``,
    ``get_column`` and ``scale_flux`` (of ``FluxMixin``). The results are
    new instances of the class (``self.__class__``).
    """

    def stitch(self, other_specs: List['SpectrumV2'], sort: bool = True) -> 'SpectrumV2':
        """
        Merge this spectrum with a list of other spectra.

        Concatenates x, y, dy, and all auxiliary columns.

        Parameters
        ----------
        other_specs : list of SpectrumV2
            List of other spectrum objects to stitch.
        sort : bool, optional
            If True, sorts the final arrays by wavelength. Defaults to ``True``.

        Returns
        -------
        SpectrumV2
            A new merged spectrum.
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
        # Overlapping inputs give interleaved pixels, which a step plot draws
        # as a zigzag: ask the GUI for a scatter plot (issue #442). Rebin and
        # coadd remove the tag again.
        ranges = sorted((np.nanmin(x), np.nanmax(x)) for x in xs if len(x))
        if any(lo < prev_hi for (_, prev_hi), (lo, _) in zip(ranges, ranges[1:])):
            new_meta['plot_style'] = 'scatter'
        
        new_data = SpectrumDataV2(
            x=DataColumnV2(x_final, x_u),
            xmin=xmin, xmax=xmax,
            y=DataColumnV2(y_final, y_u),
            dy=DataColumnV2(dy_final, y_u),
            aux_cols=new_aux,
            meta=new_meta,
            z_em=self._data.z_em
        )
        
        return self.__class__(new_data)

    def rebin(self, xstart: Optional[au.Quantity], xend: Optional[au.Quantity], 
              dx: au.Quantity, kappa: Optional[float], 
              filling: float) -> 'SpectrumV2':
        """
        Rebin the spectrum to a new uniform grid.

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
        x_new, xmin_new, xmax_new, y_new, dy_new, resol_new = rebin_spectrum(
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
            # --- Bypass interpolation if it's the carefully drizzled resolution ---
            if name == 'resol' and resol_new is not None:
                # rebin_spectrum returns R (#526)
                new_aux_cols['resol'] = DataColumnV2(resol_new, au.dimensionless_unscaled,
                                                     description="Resolving power (rebinned)")
                continue

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
        return self.__class__(data=new_data, history=new_history)
    
    def equalize_to_reference(self, reference_spec: 'SpectrumV2', order: int = 0) -> 'SpectrumV2':
        """
        Scale this spectrum's flux to match a reference spectrum.

        Computes a scaling polynomial by comparing the flux of this spectrum
        to a reference spectrum, and applies this scaling to all flux-like
        columns (y, dy, cont, model).

        Parameters
        ----------
        reference_spec : SpectrumV2
            The reference spectrum object to match.
        order : int, optional
            Polynomial order for the scaling model. ``0`` applies a constant
            scalar factor, ``-1`` applies a spline interpolation, ``>0`` applies a
            polynomial of the given order. Defaults to ``0``.

        Returns
        -------
        SpectrumV2
            A new instance with scaled flux.
        """
        from astrocook.core.spectrum_operations import compute_flux_scaling
        
        # Compute the scaling model (array) based on the reference
        scale_model = compute_flux_scaling(
            reference_spec.x.value, reference_spec.y.value, 
            self.x.value, self.y.value, 
            order=order
        )
        
        # Use the existing scale_flux method to apply it
        return self.scale_flux(scale_model)
    
    def coadd(self, others: List['SpectrumV2'], 
              xstart: Optional[au.Quantity], xend: Optional[au.Quantity], 
              dx: au.Quantity, kappa: Optional[float] = 5.0, 
              equalize_order: int = 0) -> 'SpectrumV2':
        """
        Co-add this spectrum with a list of other spectra onto a new common grid.

        Handles scaling, unit conversion, "drizzle" rebinning, and
        converts the result back to the original wavelength unit.

        Parameters
        ----------
        others : list of SpectrumV2
            List of spectra to add.
        xstart, xend : Quantity, optional
            Bounds of the new grid.
        dx : Quantity
            Step size of the new grid.
        kappa : float, optional
            Sigma clipping threshold. Defaults to ``5.0``.
        equalize_order : int, optional
            Polynomial order for flux scaling (-1=Off, 0=Scalar). Defaults to ``0``.

        Returns
        -------
        SpectrumV2
            The resulting co-added spectrum.
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
        big_resol_list = []
        resol_sources_count = 0

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
            
            # --- Capture Resolution (as R, from the column only, #526; an
            # all-NaN column no longer counts as a source) ---
            res_arr = spec.resolution()
            has_res = bool(np.isfinite(res_arr).any())
            res_arr = np.where(np.isfinite(res_arr), res_arr, 0.0)
                
            if has_res:
                resol_sources_count += 1
                
            big_resol_list.append(res_arr)

            # C. Append
            big_x_list.append(x_conv)
            big_xmin_list.append(xmin_conv)
            big_xmax_list.append(xmax_conv)
            big_y_list.append(y_val)
            big_dy_list.append(dy_val)

        # --- All-or-Nothing Resolution Decision ---
        track_resolution = (resol_sources_count == len(all_specs))

        if resol_sources_count > 0 and not track_resolution:
            logging.warning("Coadd: Resolution mismatch! Some input spectra are missing resolution data. "
                            "Resolution tracking has been disabled for this co-add.")

        # 3. Concatenate (The "Manual Stitch")
        flat_x = np.concatenate(big_x_list)
        flat_xmin = np.concatenate(big_xmin_list)
        flat_xmax = np.concatenate(big_xmax_list)
        flat_y = np.concatenate(big_y_list)
        flat_dy = np.concatenate(big_dy_list)
        flat_resol = np.concatenate(big_resol_list)

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
        # --- Pack the resolution column for the drizzler ---
        temp_aux = {}
        if track_resolution:
            temp_aux['resol'] = DataColumnV2(flat_resol, au.dimensionless_unscaled)

        temp_data = SpectrumDataV2(
            x=DataColumnV2(flat_x, calc_unit),
            xmin=DataColumnV2(flat_xmin, calc_unit),
            xmax=DataColumnV2(flat_xmax, calc_unit),
            y=DataColumnV2(flat_y, self.y.unit),
            dy=DataColumnV2(flat_dy, self.dy.unit),
            aux_cols=temp_aux,
            meta=self.meta
        )

        # 6. Rebin (Vectorized Drizzle)
        # Returns new arrays in `calc_unit` (km/s)
        reb_x, reb_xmin, reb_xmax, reb_y, reb_dy, reb_resol = rebin_spectrum(
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

        # --- Restore the co-added resolution array ---
        new_aux_cols = {}
        
        if reb_resol is not None:
            new_aux_cols['resol'] = DataColumnV2(reb_resol, au.dimensionless_unscaled, description="Co-added IVW resolution")

        # 8. Wrap in Logic Class
        final_spec_data = SpectrumDataV2(
            x=DataColumnV2(final_x.value, final_x.unit), 
            xmin=DataColumnV2(final_xmin.value, final_xmin.unit), 
            xmax=DataColumnV2(final_xmax.value, final_xmax.unit),
            y=DataColumnV2(reb_y.value, reb_y.unit), 
            dy=DataColumnV2(reb_dy.value, reb_dy.unit),
            aux_cols=new_aux_cols,
            meta=new_meta, 
            z_em=self._data.z_em
        )
        
        return self.__class__(final_spec_data)

    def resample_on_grid(self, target_grid_spec: 'SpectrumV2', 
                         fill_value: float = np.nan) -> 'SpectrumV2':
        """
        Resample this spectrum onto the exact wavelength grid of another spectrum.

        Parameters
        ----------
        target_grid_spec : SpectrumV2
            The spectrum whose grid defines the target wavelengths.
        fill_value : float, optional
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
        return self.__class__(data=new_data_core, history=new_history)

    def get_resampled_column(self, col_name: str, 
                             target_x_grid: np.ndarray, 
                             fill_value: float = np.nan) -> np.ndarray:
        """
        Lightweight API to interpolate a single column onto a new x-grid.

        Parameters
        ----------
        col_name : str
            Name of the column to interpolate.
        target_x_grid : np.ndarray
            Target grid in nanometers.
        fill_value : float, optional
            Value to fill outside the interpolation range. Defaults to ``NaN``.

        Returns
        -------
        np.ndarray
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
