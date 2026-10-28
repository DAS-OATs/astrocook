"""Line methods of :class:`~astrocook.core.spectrum.SpectrumV2`: absorbed
regions, identification of the lines, velocity bounds, redshift of a
doublet across two spectra.

Part of the split of ``spectrum.py`` into thematic modules (phase 5): the
methods are defined here, in a mixin class that ``SpectrumV2`` inherits, and
are used as before (``spec.identify_lines(...)``).
"""
import astropy.units as au
from copy import deepcopy
import dataclasses
import json
import logging
import numpy as np
from scipy.ndimage import gaussian_filter1d
from typing import TYPE_CHECKING, List

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS, xem_d
from astrocook.core.spectrum_operations import (
    _find_kinematic_doublet_candidates, detect_regions, merge_regions_by_velocity,
    rate_doublet_candidate,
)
from astrocook.core.structures import DataColumnV2

if TYPE_CHECKING:
    from astrocook.core.spectrum import SpectrumV2

#: Proximity margin of identify_lines (km/s): no candidate system beyond
#: z_em by more than this. Absorbers associated with the quasar (outflows,
#: infall) and the uncertainty of z_em from broad emission lines (up to
#: ~1000 km/s) allow absorption slightly beyond z_em; anything farther is
#: a misidentification, e.g. a SiIV at 16000 km/s in a telluric band (F55)
Z_EM_MARGIN_KMS = 3000.0


class LinesMixin:
    """Absorbed regions and identification of the lines.

    Mixed into ``SpectrumV2``, whose data and methods it uses: ``_data``,
    ``history``, ``meta``, ``x``, ``y``, ``dy``, ``cont``, ``get_column``,
    ``has_aux_column``, ``with_properties`` and ``find_absorbed`` (of
    ``ContinuumMixin``). The results are new instances of the class
    (``self.__class__``).
    """

    def detect_absorptions(self, mask_col: str = 'abs_mask', min_pix: int = 3) -> 'SpectrumV2':
        """
        Detect absorption regions based on a boolean mask.

        Groups adjacent pixels marked as 'absorbed' into discrete regions.

        Parameters
        ----------
        mask_col : str, optional
            Name of the mask column. Defaults to ``'abs_mask'``.
        min_pix : int, optional
            Minimum pixels to form a region. Defaults to ``3``.

        Returns
        -------
        SpectrumV2
            A new instance with the ``region_id`` column updated.
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
        
        return self.__class__(data=new_data, history=new_hist)
    
    def identify_lines(self, 
                       multiplet_list: List[str],
                       mask_col: str = 'abs_mask', 
                       min_pix_region: int = 3,
                       merge_dv: float = 10.0,
                       score_threshold: float = 0.5,
                       bypass_scoring: bool = False,
                       debug_rating: bool = False) -> 'SpectrumV2':
        """
        Identify absorption lines by correlating regions with multiplet templates.

        Uses a de-blending strategy to separate merged features. A doublet
        gives one candidate system per absorption component, from its
        stronger line, and no candidate lies beyond ``z_em`` by more than
        ``Z_EM_MARGIN_KMS`` (3000 km/s).

        Parameters
        ----------
        multiplet_list : list of str
            List of multiplets to search for.
        mask_col : str, optional
            Name of the mask column.
        min_pix_region : int, optional
            Minimum pixels for a region.
        merge_dv : float, optional
            Velocity threshold for merging regions (km/s).
        score_threshold : float, optional
            R^2 score threshold for identification.
        bypass_scoring : bool, optional
            Accept all candidates if True.
        debug_rating : bool, optional
            Show debug plots if True.

        Returns
        -------
        SpectrumV2
            A new instance with the ``abs_ids`` column and identification metadata.
        """
        # This is the 'self' spectrum object
        spec = self
        z_em = spec._data.z_em
        if z_em == 0.0:
            raise ValueError("z_em must be set before running identify_lines.")
        
        # Import SciPy signal processing tools
        from scipy.signal import find_peaks, peak_widths
        
        logging.info("identify_lines (De-blending): Starting identification...")

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
        
        # 2. Get Merged Regions
        try:
            mask_data = spec.get_column(mask_col).value
            if mask_data.dtype.kind != 'b':
                 mask_data = mask_data.astype(bool)
            
            region_id_map_raw, num_raw = detect_regions(mask_data, min_pix_region)
            region_id_map_merged, num_merged = merge_regions_by_velocity(
                region_id_map_raw, spec.x, spec._data.z_em, merge_dv
            )
            logging.info(f"identify_lines: {num_merged} regions remain after merging.")
        except Exception as e:
            logging.error(f"Failed to detect or merge absorption regions: {e}", exc_info=True)
            raise e

        # Define AOD profile
        cont = spec.cont.value if spec.cont is not None else np.ones_like(spec.y.value)
        flux = spec.y.value
        with np.errstate(divide='ignore', invalid='ignore'):
            aod = -np.log(flux / cont)
        aod[~np.isfinite(aod)] = 0.0 

        # Define boundaries
        x_nm = spec.x.to_value(au.nm)
        lya_limit_obs_nm = (1.0 + z_em) * xem_d['Ly_lim'].to_value(au.nm)
        lya_obs_nm = (1.0 + z_em) * xem_d['Ly_a'].to_value(au.nm)
        
        forest_mask = (x_nm < lya_obs_nm) & (x_nm > lya_limit_obs_nm)
        redward_mask = (x_nm >= lya_obs_nm)
        
        forest_map = region_id_map_merged.copy()
        forest_map[~forest_mask] = 0 
        
        redward_map = region_id_map_merged.copy()
        redward_map[~redward_mask] = 0

        # --- Kinematic Checks (Steps 2-4) ---
        all_candidates = [] 
        metal_doublets = [
            m for m in multiplet_list 
            if m in STANDARD_MULTIPLETS and len(STANDARD_MULTIPLETS[m]) == 2 and m != 'Ly_ab'
        ]
        
        # Candidate systems at most Z_EM_MARGIN_KMS beyond z_em (F55)
        z_max = z_em + Z_EM_MARGIN_KMS / 299792.458 * (1.0 + z_em)
        for series_name in metal_doublets:
            candidates_red = _find_kinematic_doublet_candidates(
                redward_map, spec.x, spec._data.z_em, series_name, z_em
            )
            for (rid_1, rid_2, z_test) in candidates_red:
                if z_test > z_max:
                    continue
                all_candidates.append((rid_1, rid_2, z_test, series_name, score_threshold))
        
        if 'Ly_ab' in multiplet_list:
            lyab_candidates = _find_kinematic_doublet_candidates(
                forest_map, spec.x, spec._data.z_em, 'Ly_ab', z_em
            )
            for (rid_1, rid_2, z_test) in lyab_candidates:
                all_candidates.append((rid_1, rid_2, z_test, 'Ly_ab', score_threshold))

        # --- Rate Candidates ---
        reliable_ids = {} 
        for (rid_1, rid_2, z_test, series_name, threshold) in all_candidates:
            if series_name == 'Ly_ab':
                threshold = score_threshold 
            else:
                region_mask_for_check = (region_id_map_merged == rid_1)
                x_region_mean = np.mean(x_nm[region_mask_for_check])
                if x_region_mean < lya_obs_nm:
                    threshold = 1.1 # Impossible threshold for metals in forest
                else:
                    threshold = score_threshold 
                
            mask_1 = (region_id_map_merged == rid_1)
            mask_2 = (region_id_map_merged == rid_2)
            
            if not np.any(mask_1) or not np.any(mask_2): continue
            
            r2_score = rate_doublet_candidate(
                spec, mask_1, mask_2, series_name, z_test, debug_rating=debug_rating
            )
            
            if r2_score > threshold or bypass_scoring:
                lines = STANDARD_MULTIPLETS[series_name]
                comp_1, comp_2 = lines[0], lines[1]
                rid_1_int = int(rid_1) 
                rid_2_int = int(rid_2) 
                
                # Both regions are labelled, but only the region of the
                # first (stronger) line gives the candidate systems: the
                # other one would give the same systems again, ~1 km/s
                # apart from the pixel sampling (F49)
                if rid_1_int not in reliable_ids: reliable_ids[rid_1_int] = []
                reliable_ids[rid_1_int].append((comp_2, r2_score, z_test, False))

                if rid_2_int not in reliable_ids: reliable_ids[rid_2_int] = []
                reliable_ids[rid_2_int].append((comp_1, r2_score, z_test, True))

        # --- Forest Cleanup ---
        forest_rids = np.unique(forest_map)
        forest_rids = forest_rids[forest_rids > 0]
        
        for rid in forest_rids:
            rid_int = int(rid)
            if rid_int not in reliable_ids:
                region_mask = (region_id_map_merged == rid_int)
                aod_region = aod.copy()
                aod_region[~region_mask] = -np.inf
                
                aod_peak_idx = np.argmax(aod_region)
                if aod_peak_idx > 0:
                    x_peak_nm = x_nm[aod_peak_idx]
                    z_test = (x_peak_nm / xem_d['Ly_a'].to_value(au.nm)) - 1.0
                    reliable_ids[rid_int] = [('Ly_a', 0.0, z_test, True)]

        logging.info(f"Finalized {len(reliable_ids)} identified regions (pre-deblending).")

        # --- 6. Save results with SMART DE-BLENDING ---
        
        final_vis_map = np.zeros_like(region_id_map_merged)
        system_series_map = {}
        system_z_map = {}
        system_logN_map = {}
        system_b_map = {}
        
        # [NEW] Rich metadata for UI Tooltips
        region_info = {}

        sys_id_counter = 1
        processed_systems = set()
        c_kms = 299792.458
        
        # Pre-calculations
        with np.errstate(divide='ignore', invalid='ignore'):
            dv_pix = c_kms * np.abs(np.gradient(spec.x.value)) / spec.x.value
            aod_clean = np.nan_to_num(aod, nan=0.0, posinf=0.0, neginf=0.0)
            aod_clean = np.maximum(aod_clean, 0.0)

        for rid, candidates in reliable_ids.items():
            rid_int = int(rid)
            region_mask = (region_id_map_merged == rid_int)
            final_vis_map[region_mask] = rid_int
            
            # Populate region_info for tooltips: {id: [[series, score], ...]}
            # We convert score to float to ensure JSON serializability
            region_info[rid_int] = [[c[0], float(c[1])] for c in candidates]

            # --- DE-BLENDING ---
            flux_seg = spec.y.value[region_mask]
            err_seg = spec.dy.value[region_mask]
            x_seg = spec.x.value[region_mask]
            
            if len(flux_seg) < 3:
                peaks = [np.argmin(flux_seg)]
                widths_pix = [len(flux_seg)]
            else:
                inv_flux = 1.0 - flux_seg
                inv_flux_smooth = gaussian_filter1d(inv_flux, sigma=1.0)
                local_sigma = np.nanmedian(err_seg) if np.all(err_seg > 0) else 0.01
                dynamic_prominence = max(0.02, 3.0 * local_sigma)
                peaks, _ = find_peaks(inv_flux_smooth, prominence=dynamic_prominence, width=1.0)
                
                if len(peaks) == 0:
                    peaks = [np.argmax(inv_flux_smooth)]
                    widths_pix = [len(flux_seg)]
                else:
                    widths_pix = peak_widths(inv_flux_smooth, peaks, rel_height=0.5)[0]

            for i, peak_idx in enumerate(peaks):
                x_peak_obs = x_seg[peak_idx]
                local_dv = dv_pix[region_mask][peak_idx]
                fwhm_kms = widths_pix[i] * local_dv
                b_guess = max(5.0, min(100.0, fwhm_kms / 1.665))

                for (component_name, score, z_test_orig, gives_systems) in candidates:
                    if not gives_systems:
                        continue
                    system_name = component_name
                    if component_name not in STANDARD_MULTIPLETS:
                        for mult, members in STANDARD_MULTIPLETS.items():
                            if component_name in members:
                                system_name = mult
                                break
                    
                    if component_name in ATOM_DATA:
                        lam_rest = ATOM_DATA[component_name]['wave'] / 10.0 
                        z_peak = (x_peak_obs / lam_rest) - 1.0
                    else:
                        z_peak = z_test_orig
                        
                    p_start = max(0, int(peak_idx - widths_pix[i]))
                    p_end = min(len(aod_clean[region_mask]), int(peak_idx + widths_pix[i] + 1))
                    local_tau = np.sum(aod_clean[region_mask][p_start:p_end] * dv_pix[region_mask][p_start:p_end])
                    
                    if component_name in ATOM_DATA:
                        atom = ATOM_DATA[component_name]
                        if atom['f'] > 0:
                            N_val = 3.768e14 * local_tau / (atom['f'] * atom['wave'])
                            logN_guess = max(11.0, min(19.0, np.log10(N_val) if N_val > 0 else 13.0))
                        else: logN_guess = 13.5
                    else: logN_guess = 13.5

                    if z_peak > z_max:
                        continue  # a peak of a region that straddles z_max (F55)
                    z_round = round(z_peak, 5)
                    uniq_key = (system_name, z_round)
                    
                    if uniq_key not in processed_systems:
                        system_series_map[sys_id_counter] = system_name
                        system_z_map[sys_id_counter] = z_peak
                        system_logN_map[sys_id_counter] = logN_guess
                        system_b_map[sys_id_counter] = b_guess
                        sys_id_counter += 1
                        processed_systems.add(uniq_key)

        # --- 7. Create Final Data Core ---
        new_meta = self.meta.copy()
        new_aux_cols = deepcopy(self._data.aux_cols)

        # Storage for Identification
        new_meta['series_map_json'] = json.dumps(system_series_map)
        new_meta['z_map_json'] = json.dumps(system_z_map)
        new_meta['logN_map_json'] = json.dumps(system_logN_map)
        new_meta['b_map_json'] = json.dumps(system_b_map)
        
        # rich UI Metadata [THE FIX]
        new_meta['region_identifications'] = json.dumps(region_info)

        vis_col = DataColumnV2(
            values=final_vis_map,
            unit=au.dimensionless_unscaled,
            description="Absorption regions with identifications"
        )
        new_aux_cols['abs_ids'] = vis_col

        final_data_core = dataclasses.replace(self._data, meta=new_meta, aux_cols=new_aux_cols)
        new_history = self.history + [f"Identified {len(processed_systems)} systems (de-blended)"]
        return self.__class__(data=final_data_core, history=new_history)
    
    def get_velocity_bounds(self, z_target: float, trans: str) -> tuple:
        """
        Calculate the min/max velocity coverage of this spectrum relative to a target.

        Parameters
        ----------
        z_target : float
            The target redshift.
        trans : str
            Transition name (e.g., ``'Ly_a'``).

        Returns
        -------
        tuple
            (min_velocity, max_velocity) in km/s. Returns ``(-inf, inf)`` if invalid.
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
        Find the redshift of maximum coincidence between this spectrum and another.

        Calculates the cross-correlation (product of absorption depths) to align
        two spectra covering different lines of a doublet.

        Parameters
        ----------
        other_spec : SpectrumV2
            The other spectrum object.
        trans_self : str
            Transition in this spectrum (e.g., ``'OVI_1031'``).
        trans_other : str
            Transition in the other spectrum (e.g., ``'OVI_1037'``).

        Returns
        -------
        float
            The detected redshift, or 0.0 on failure.
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
