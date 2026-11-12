"""Flux methods of :class:`~astrocook.core.spectrum.Spectrum`: scaling,
running standard deviation, smoothing, calibration.

Part of the split of ``spectrum.py`` into thematic modules (phase 5): the
methods are defined here, in a mixin class that ``Spectrum`` inherits, and
are used as before (``spec.smooth(...)``).
"""
import astropy.units as au
from copy import deepcopy
import dataclasses
import logging
import numpy as np
from typing import Self, Union

from astrocook.core.photometry import generate_calibration_curve
from astrocook.core.spectrum_operations import running_std, smooth_spectrum
from astrocook.core.structures import DataColumn, SpectrumData
from astrocook.core._spectrum_host import HostBase
from astrocook.core.units import FLAM, FLUX_AUX_COLUMNS, has_no_unit, is_arbitrary_flux

logger = logging.getLogger(__name__)

# The auxiliary columns that the flux methods change along with y and dy
# (findings register F69). Chosen by name, not by unit: many loaders give
# the flux no unit, and then the unit cannot tell a continuum from a
# transmission, a resolution or a quality flag.
# Flux columns: scaled and calibrated as the flux, and smoothed.
# (FLUX_AUX_COLUMNS, defined in core/units.py since PR 8-bis, where
# y_convert uses it too)
# Transmissions (dimensionless, multiplied by the continuum in
# Spectrum.norm): smoothed, as the flux they apply to, but never scaled.
TRANSMISSION_AUX_COLUMNS = ('telluric_model',)


class FluxMixin(HostBase):
    """Operations on the flux: scaling, running standard deviation,
    smoothing (of the spectrum or of one column), calibration on magnitudes.

    Mixed into ``Spectrum``, whose data and methods it uses: ``_data``,
    ``history``, ``x``, ``y``, ``dy``, ``t`` and ``_renormalize_model``. The
    results are new instances of the class (``self.__class__``).
    """

    def scale_flux(self, scale_factor: Union[float, np.ndarray]) -> Self:
        """
        Scale the flux (y, dy) and the flux auxiliary columns (cont, model,
        cont_pl).

        Parameters
        ----------
        scale_factor : float or np.ndarray
            Multiplicative factor.

        Returns
        -------
        Spectrum
            A new instance with scaled flux.
        """
        # 1. Scale core columns
        # Handle cases where scale_factor is scalar or array
        new_y = DataColumn(self.y.value * scale_factor, self.y.unit)
        new_dy = DataColumn(self.dy.value * scale_factor, self.dy.unit)

        # 2. Scale flux-like auxiliary columns
        new_aux_cols = deepcopy(self._data.aux_cols)
        for name in FLUX_AUX_COLUMNS:
            if name in new_aux_cols:
                col = new_aux_cols[name]
                # Ensure we propagate description and unit
                new_aux_cols[name] = DataColumn(
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
        
        # 4. Return new Spectrum
        # (Optional: summarize scalar or mean of array for history)
        if np.size(scale_factor) == 1:
            info = f"{scale_factor:.3f}"
        else:
            info = f"array(mean={np.mean(scale_factor):.3f})"
            
        new_history = self.history + [f"Scaled flux by {info}"]
        return self.__class__(data=new_data, history=new_history)

    def calculate_running_std(self, 
                              input_col: str, 
                              output_col: str, 
                              window_pix: int) -> Self:
        """
        Calculate a running RMS on a column and save it as a new column.

        Parameters
        ----------
        input_col : str
            Name of the column to calculate statistics on (e.g., ``'y'``).
        output_col : str
            Name of the output column (e.g., ``'running_std'``).
        window_pix : int
            Width of the window in pixels.

        Returns
        -------
        Spectrum
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
        
        new_std_col = DataColumn(
            values=rms_values,
            unit=input_col_data.unit, # RMS has same unit as the input
            description=f"Running StdDev (window={window_pix}pix) of {input_col}"
        )

        # 4. Create new SpectrumData by adding the new column
        new_aux_cols = deepcopy(self._data.aux_cols)
        new_aux_cols[output_col] = new_std_col
        
        new_data = dataclasses.replace(self._data, aux_cols=new_aux_cols)
        
        # 5. Return new Spectrum
        new_history = self.history + [f"Calculated running RMS (window={window_pix}pix) on {input_col}"]
        return self.__class__(data=new_data, history=new_history)

    def smooth(self, sigma_kms: float) -> Self:
        """
        Apply Gaussian smoothing to the flux (y, dy), the flux auxiliary
        columns (cont, model, cont_pl) and the transmission (telluric_model).

        Parameters
        ----------
        sigma_kms : float
            Standard deviation of the Gaussian kernel in km/s.

        Returns
        -------
        Spectrum
            A new instance with smoothed columns.
        """
        
        # 1. Smooth the core y and dy columns
        new_y_values = smooth_spectrum(
            x=self._data.x.quantity,
            y=self._data.y.values,
            sigma_kms=sigma_kms,
            z_ref=self._data.z_em
        )
        new_y_col = DataColumn(new_y_values, self._data.y.unit, "Smoothed Flux")

        new_dy_values = smooth_spectrum(
            x=self._data.x.quantity,
            y=self._data.dy.values,
            sigma_kms=sigma_kms,
            z_ref=self._data.z_em
        )
        new_dy_col = DataColumn(new_dy_values, self._data.dy.unit, "Smoothed Error")

        # 2. Smooth flux-like auxiliary columns
        new_aux_cols = deepcopy(self._data.aux_cols)
        for col_name in FLUX_AUX_COLUMNS + TRANSMISSION_AUX_COLUMNS:
            if col_name in new_aux_cols:
                try:
                    old_col = new_aux_cols[col_name]
                    if old_col.values.dtype.kind in 'fiu':
                        logger.debug(f"Smoothing auxiliary column: {col_name}")
                        new_vals = smooth_spectrum(
                            self._data.x.quantity,
                            old_col.values,
                            sigma_kms,
                            self._data.z_em
                        )
                        
                        new_aux_cols[col_name] = DataColumn(
                            new_vals, old_col.unit, old_col.description
                        )
                except (ValueError, TypeError) as e:
                    logger.warning(f"The column '{col_name}' could not be smoothed "
                                   f"({e}): it is kept as it was.")
        
        # --- RE-NORMALIZATION BLOCK HAS BEEN REMOVED ---

        # 3. Create new SpectrumData
        new_data = dataclasses.replace(
            self._data, 
            y=new_y_col, 
            dy=new_dy_col, 
            aux_cols=new_aux_cols
        )
        
        # 4. Return new Spectrum
        new_history = self.history + [f"Smoothed spectrum (sigma={sigma_kms} km/s)"]
        return self.__class__(data=new_data, history=new_history)

    def smooth_column(self, target_col: str, sigma_kms: float, renorm_model: bool = False) -> Self:
        """
        Apply Gaussian smoothing to a single target column.

        Parameters
        ----------
        target_col : str
            Name of the column to smooth.
        sigma_kms : float
            Standard deviation of the Gaussian kernel in km/s.
        renorm_model : bool, optional
            If True and ``target_col`` is ``'cont'``, the ``'model'`` column will be re-normalized.

        Returns
        -------
        Spectrum
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
        
        new_data_col = DataColumn(
            values=new_values,
            unit=old_col_data.unit,
            description=f"Smoothed {target_col}"
        )
        
        # 4. Create new SpectrumData using the immutable pattern
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
        new_data_core = SpectrumData(
            x=core_cols['x'], xmin=core_cols['xmin'], xmax=core_cols['xmax'],
            y=core_cols['y'], dy=core_cols['dy'],
            aux_cols=aux_cols,
            meta=deepcopy(self._data.meta),
            z_em=self._data.z_em
        )

        # 8. Return a NEW Spectrum instance
        new_history = self.history + [f"Smoothed column: {target_col} (sigma={sigma_kms} km/s)"]
        return self.__class__(data=new_data_core, history=new_history)

    def with_flux_unit(self, unit: au.UnitBase) -> Self:
        """
        Declare the unit of the flux, without converting the values (phase
        6-quater, PR 8-bis; ``y_convert`` converts).

        The unit is given to ``y``, ``dy``, the flux columns (``cont``,
        ``model``, ``cont_pl``) in the unit of the flux or without one, and
        the other columns in the unit of the flux when this is not
        dimensionless (a dimensionless quality flag is not a flux).

        Parameters
        ----------
        unit : astropy.units.Unit
            The unit of the flux, e.g. ``erg/(Angstrom s cm2)`` or
            ``arb_flux``.

        Returns
        -------
        Spectrum
            A new instance with the same values and the unit declared.
        """
        old = self._data.y.unit

        def relabel(col: DataColumn) -> DataColumn:
            return DataColumn(col.values, unit, col.description)

        new_aux_cols = dict(self._data.aux_cols)
        for name, col in new_aux_cols.items():
            same = col.unit is not None and old is not None and au.Unit(col.unit) == au.Unit(old)
            if name in FLUX_AUX_COLUMNS and (same or has_no_unit(col.unit)):
                new_aux_cols[name] = relabel(col)
            elif same and not has_no_unit(old):
                new_aux_cols[name] = relabel(col)
        new_data = dataclasses.replace(self._data, y=relabel(self._data.y),
                                       dy=relabel(self._data.dy), aux_cols=new_aux_cols)
        return self.__class__(data=new_data,
                              history=self.history + [f"Flux unit declared: {unit}"])

    def calibrate(self, magnitudes: str) -> Self:
        """
        Calibrate (and warp) the flux to match target photometric magnitudes.

        Parameters
        ----------
        magnitudes : str
            Comma-separated pairs of ``Filter=Mag`` (e.g., ``"SDSS_g=17.5, SDSS_r=17.0"``).

        Returns
        -------
        Spectrum
            A new instance with calibrated flux.
        """
        # 1. Get the pure numerical correction curve
        correction_curve = generate_calibration_curve(
            self.x, self.y, magnitudes
        )

        # 2. The unit of the result: the curve scales the values in their
        #    own unit, so a physical unit is kept; a flux without one (or in
        #    arb_flux), read as erg/(Angstrom s cm2) by the photometry, takes
        #    it (F131: the result was always labelled erg/(Angstrom s cm2),
        #    so a flux in erg/(nm s cm2) came out 10 times brighter than it said)
        old_unit = self._data.y.unit
        if has_no_unit(old_unit) or is_arbitrary_flux(old_unit):
            new_unit = FLAM
        else:
            new_unit = old_unit

        # 3. Helper to apply curve and unit to a column
        def apply_corr(col: DataColumn) -> DataColumn:
            return DataColumn(col.values * correction_curve, new_unit, col.description)

        # 4. Apply to all relevant columns
        new_y = apply_corr(self._data.y)
        new_dy = apply_corr(self._data.dy)
        
        new_aux_cols = deepcopy(self._data.aux_cols)
        # Not the transmission: it stays dimensionless (F69)
        for col_name in FLUX_AUX_COLUMNS:
            if col_name in new_aux_cols:
                 new_aux_cols[col_name] = apply_corr(new_aux_cols[col_name])

        # 5. Return new state
        new_data = dataclasses.replace(self._data, y=new_y, dy=new_dy, aux_cols=new_aux_cols)
        new_history = self.history + [f"Flux calibrated with: {magnitudes}"]
        return self.__class__(data=new_data, history=new_history)
