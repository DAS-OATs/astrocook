"""Methods of :class:`~astrocook.core.spectrum.Spectrum` that combine or
regrid spectra: stitch, rebin, equalization, co-addition, resampling.

Part of the split of ``spectrum.py`` into thematic modules (phase 5): the
methods are defined here, in a mixin class that ``Spectrum`` inherits, and
are used as before (``spec.rebin(...)``).
"""
import astropy.units as au
from copy import deepcopy
from itertools import pairwise
import logging
import numpy as np
from typing import TYPE_CHECKING, Self, cast, Dict, List, Optional

from astrocook.core.spectrum_operations import convert_axis_velocity, rebin_spectrum
from astrocook.core.structures import DataColumn, SpectrumData
from astrocook.core._spectrum_host import HostBase

if TYPE_CHECKING:
    from astrocook.core.spectrum import Spectrum

logger = logging.getLogger(__name__)


def _flux_values(flux: au.Quantity, unit) -> np.ndarray:
    """The values of ``flux`` in ``unit`` when it can be converted, else as
    they are (as they are also when the unit is the same)."""
    if flux.unit == unit or not au.Unit(flux.unit).is_equivalent(unit):
        return flux.value
    return flux.to_value(unit)


def _in_units_of_this(spec, x_unit, y_unit):
    """The values of the axis and of the flux of ``spec`` in ``x_unit`` and
    ``y_unit`` (phase 6-quater, PR 9, F133, F134, F138): the axis converted
    between wavelengths (ValueError for an axis of another kind, e.g. a
    velocity against a wavelength, which cannot be compared), the flux when
    it can be (see :func:`_flux_values`)."""
    x = spec.x
    if x.unit != x_unit:
        if not (au.Unit(x.unit).is_equivalent(au.nm) and au.Unit(x_unit).is_equivalent(au.nm)):
            raise ValueError(f"an axis in {x.unit} cannot be compared with one in "
                             f"{x_unit}: convert them to wavelengths with x_convert")
        x_values = x.to_value(x_unit)
    else:
        x_values = x.value
    return x_values, _flux_values(spec.y, y_unit)


class CombineMixin(HostBase):
    """Stitching, rebinning, equalization, co-addition and resampling of
    spectra.

    Mixed into ``Spectrum``, whose data and methods it uses: ``_data``,
    ``history``, ``meta``, ``x``, ``y``, ``dy``, ``t``, ``resolution``,
    ``get_column`` and ``scale_flux`` (of ``FluxMixin``). The results are
    new instances of the class (``self.__class__``).
    """

    def stitch(self, other_specs: List['Spectrum'], sort: bool = True) -> Self:
        """
        Merge this spectrum with a list of other spectra.

        Concatenates x, y, dy, and all auxiliary columns.

        Parameters
        ----------
        other_specs : list of Spectrum
            List of other spectrum objects to stitch.
        sort : bool, optional
            If True, sorts the final arrays by wavelength. Defaults to ``True``.

        Returns
        -------
        Spectrum
            A new merged spectrum.
        """
        from copy import deepcopy
        
        all_specs = [self] + other_specs
        
        xs, ys, dys = [], [], []
        aux_map: Dict[str, list] = {k: [] for k in self._data.aux_cols.keys()}
        
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
                    val = cast(au.Quantity, s.get_column(col_name)).to(t_unit).value
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
        from astrocook.core.structures import DataColumn
        new_aux = {}
        for col_name, arrays in aux_map.items():
            arr = np.concatenate(arrays)
            if sort: arr = arr[idx]
            # Copy metadata from self
            orig_col = self._data.aux_cols[col_name]
            new_aux[col_name] = DataColumn(arr, orig_col.unit, orig_col.description)

        # Build New Object
        from ..io.loaders import _auto_limits
        xmin, xmax = _auto_limits(x_final, x_u)
        new_meta = deepcopy(self.meta)
        new_meta['stitched'] = True
        # Overlapping inputs give interleaved pixels, which a step plot draws
        # as a zigzag: ask the GUI for a scatter plot (issue #442). Rebin and
        # coadd remove the tag again.
        ranges = sorted((np.nanmin(x), np.nanmax(x)) for x in xs if len(x))
        if any(lo < prev_hi for (_, prev_hi), (lo, _) in pairwise(ranges)):
            new_meta['plot_style'] = 'scatter'
        
        new_data = SpectrumData(
            x=DataColumn(x_final, x_u),
            xmin=xmin, xmax=xmax,
            y=DataColumn(y_final, y_u),
            dy=DataColumn(dy_final, y_u),
            aux_cols=new_aux,
            meta=new_meta,
            z_em=self._data.z_em
        )
        
        return self.__class__(new_data)

    def rebin(self, xstart: Optional[au.Quantity], xend: Optional[au.Quantity], 
              dx: au.Quantity, kappa: Optional[float], 
              filling: float) -> Self:
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
        Spectrum
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
                new_aux_cols['resol'] = DataColumn(resol_new, au.dimensionless_unscaled,
                                                     description="Resolving power (rebinned)")
                continue

            # --- *** FIX 1: Check if 'values' is None from a *previous* rebin *** ---
            if col_data.values is None:
                logger.debug(f"Skipping aux col '{name}' as its data is None.")
                continue
            # --- *** END FIX 1 *** ---

            if col_data.values.dtype.kind in 'fiub': # Check for float, int, unsigned, or bool
                try:
                    logger.debug(f"Rebin: Interpolating aux column '{name}'...")
                    interp_values = self.get_resampled_column(
                        name, 
                        x_final_nm 
                    )
                    
                    # --- *** FIX 2: Check if interpolation returned None (e.g., for string) *** ---
                    if interp_values is not None:
                        new_aux_cols[name] = DataColumn(interp_values, col_data.unit)
                    # (If it is None, get_resampled_column already warned us, so we just drop it)
                    # --- *** END FIX 2 *** ---

                except (ValueError, TypeError, KeyError) as e:
                    logger.warning(f"The column '{name}' could not be rebinned "
                                   f"({e}): it is dropped.")
            else:
                logger.warning(f"The column '{name}' is not numerical and cannot "
                               "be rebinned: it is dropped.")
        # --- *** END FIX *** ---

        # Create a copy of the metadata to avoid modifying the original session
        new_meta = self._data.meta.copy()

        # Remove the 'plot_style' tag so the GUI reverts to the default (Line)
        if 'plot_style' in new_meta:
            del new_meta['plot_style']

        # 4. Build new SpectrumData
        new_data = SpectrumData(
            x=DataColumn(x_final.value, x_final.unit),
            xmin=DataColumn(xmin_final.value, xmin_final.unit),
            xmax=DataColumn(xmax_final.value, xmax_final.unit),
            y=DataColumn(y_new.value, y_new.unit),
            dy=DataColumn(dy_new.value, dy_new.unit),
            aux_cols=new_aux_cols, 
            meta=new_meta, 
            z_em=self._data.z_em   # Propagate
        )
        
        # 5. Return a NEW Spectrum instance
        new_history = self.history + [f"Rebinned spectrum (dx={dx})"]
        return self.__class__(data=new_data, history=new_history)
    
    def equalize_to_reference(self, reference_spec: 'Spectrum', order: int = 0) -> Self:
        """
        Scale this spectrum's flux to match a reference spectrum.

        Computes a scaling polynomial by comparing the flux of this spectrum
        to a reference spectrum, and applies this scaling to all flux-like
        columns (y, dy, cont, model).

        Parameters
        ----------
        reference_spec : Spectrum
            The reference spectrum object to match.
        order : int, optional
            Polynomial order for the scaling model. ``0`` applies a constant
            scalar factor, ``-1`` applies a spline interpolation, ``>0`` applies a
            polynomial of the given order. Defaults to ``0``.

        Returns
        -------
        Spectrum
            A new instance with scaled flux, in the unit of this spectrum.

        Notes
        -----
        The two spectra are compared in the units of this one (F138, phase
        6-quater, PR 9): the axis of the reference converted to the unit of
        this axis, its flux to the unit of this flux when they are of one
        kind. Before, the values were compared as they were, so a reference
        with the axis in Å against this one in nm had no overlap (no
        scaling), and one per Å against this one per nm gave a flux off by
        their ratio. A flux of another kind (arbitrary, against a physical
        one) is compared by its values: the scaling calibrates it, and the
        recipe declares the unit. Nothing is converted when the units are
        the same.
        """
        from astrocook.core.spectrum_operations import compute_flux_scaling

        ref_x, ref_y = _in_units_of_this(reference_spec, self.x.unit, self.y.unit)
        # Compute the scaling model (array) based on the reference
        scale_model = compute_flux_scaling(
            ref_x, ref_y,
            self.x.value, self.y.value, 
            order=order
        )
        
        # Use the existing scale_flux method to apply it
        return self.scale_flux(scale_model)
    
    def coadd(self, others: List['Spectrum'], 
              xstart: Optional[au.Quantity], xend: Optional[au.Quantity], 
              dx: au.Quantity, kappa: Optional[float] = 5.0, 
              equalize_order: int = 0) -> Self:
        """
        Co-add this spectrum with a list of other spectra onto a new common grid.

        Handles scaling, unit conversion, "drizzle" rebinning, and
        converts the result back to the original wavelength unit.

        Parameters
        ----------
        others : list of Spectrum
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
        Spectrum
            The resulting co-added spectrum.
        """
        from astrocook.core.spectrum_operations import (
            compute_flux_scaling, 
            rebin_spectrum, 
            convert_axis_velocity
        )
        from astrocook.core.structures import SpectrumData, DataColumn

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
            # The axis and the flux in the units of the reference (F133,
            # F134, phase 6-quater, PR 9: the values were taken as they
            # were); a flux of another kind by its values, which the
            # scaling calibrates. As they are when the units are the same
            spec_y = _flux_values(spec.y, self.y.unit)
            spec_dy = _flux_values(spec.dy, self.y.unit)
            # A. Flux Scaling (Skip for self/index 0)
            if i > 0 and equalize_order >= 0:
                logger.debug(f"Scaling spectrum {i} to the reference.")
                spec_x, _ = _in_units_of_this(spec, original_x_unit, self.y.unit)
                # Compute scaling in original frame
                scale_model = compute_flux_scaling(
                    ref_x, ref_y, 
                    spec_x, spec_y, 
                    order=equalize_order
                )
                y_val = spec_y * scale_model
                dy_val = spec_dy * scale_model
            else:
                y_val = spec_y
                dy_val = spec_dy

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
            logger.warning("Some of the spectra to combine have no resolution: "
                           "the combined spectrum has no resolution column.")

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
            temp_aux['resol'] = DataColumn(flat_resol, au.dimensionless_unscaled)

        temp_data = SpectrumData(
            x=DataColumn(flat_x, calc_unit),
            xmin=DataColumn(flat_xmin, calc_unit),
            xmax=DataColumn(flat_xmax, calc_unit),
            y=DataColumn(flat_y, self.y.unit),
            dy=DataColumn(flat_dy, self.dy.unit),
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
            logger.debug(f"Converting the result back to {original_x_unit}.")
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
            new_aux_cols['resol'] = DataColumn(reb_resol, au.dimensionless_unscaled, description="Co-added IVW resolution")

        # 8. Wrap in Logic Class
        final_spec_data = SpectrumData(
            x=DataColumn(final_x.value, final_x.unit), 
            xmin=DataColumn(final_xmin.value, final_xmin.unit), 
            xmax=DataColumn(final_xmax.value, final_xmax.unit),
            y=DataColumn(reb_y.value, reb_y.unit), 
            dy=DataColumn(reb_dy.value, reb_dy.unit),
            aux_cols=new_aux_cols,
            meta=new_meta, 
            z_em=self._data.z_em
        )
        
        return self.__class__(final_spec_data)

    def resample_on_grid(self, target_grid_spec: 'Spectrum', 
                         fill_value: float = np.nan) -> Self:
        """
        Resample this spectrum onto the exact wavelength grid of another spectrum.

        Parameters
        ----------
        target_grid_spec : Spectrum
            The spectrum whose grid defines the target wavelengths.
        fill_value : float, optional
            Value to fill outside the interpolation range (default: NaN).

        Returns
        -------
        Spectrum
            A new instance interpolated to the target grid.
        """
        x_old = self.x.to_value(au.nm)
        x_new = target_grid_spec.x.to_value(au.nm)
        all_cols_old = self.t._data_dict
        x_col = DataColumn(x_new, target_grid_spec.x.unit)
        xmin_col = DataColumn(target_grid_spec.xmin.value, target_grid_spec.xmin.unit)
        xmax_col = DataColumn(target_grid_spec.xmax.value, target_grid_spec.xmax.unit)
        y_col_new = None
        dy_col_new = None
        aux_cols_new = {}
        for name, col_q in all_cols_old.items():
            if name in ('x', 'xmin', 'xmax'):
                continue
            if col_q.dtype.kind in 'fiu': 
                logger.debug(f"Resampling column '{name}'...")
                y_interp = np.interp(x_new, x_old, col_q.value,
                                     left=fill_value, right=fill_value)
                new_data_col = DataColumn(y_interp, col_q.unit)
                if name == 'y':
                    y_col_new = new_data_col
                elif name == 'dy':
                    dy_col_new = new_data_col
                else:
                    aux_cols_new[name] = new_data_col
            else:
                logger.warning(f"The column '{name}' is not numerical and cannot "
                               "be resampled: it is dropped.")
        assert y_col_new is not None and dy_col_new is not None
        new_data_core = SpectrumData(
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
            raise ValueError(f"Column '{col_name}' not found in the spectrum.")
        
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
