# v2/photometry.py

import numpy as np
import astropy.units as au
import astropy.constants as const
import logging
from typing import Dict, Optional

# --- Approximate Standard Filter Definitions ---
# Format: 'Name': (Central Wavelength in nm, FWHM in nm)
STANDARD_FILTERS = {
    'Johnson_U': (365.0, 66.0),
    'Johnson_B': (445.0, 94.0),
    'Johnson_V': (551.0, 88.0),
    'Johnson_R': (658.0, 138.0),
    'Johnson_I': (806.0, 149.0),
    'SDSS_u': (354.3, 57.0),
    'SDSS_g': (477.0, 140.0),
    'SDSS_r': (623.1, 140.0),
    'SDSS_i': (762.5, 150.0),
    'SDSS_z': (913.4, 95.0),
}

def get_filter_transmission(filter_name: str, target_x_nm: np.ndarray) -> np.ndarray:
    """
    Generates a synthetic Gaussian transmission curve for a standard filter
    on the given wavelength grid.
    """
    if filter_name not in STANDARD_FILTERS:
        raise ValueError(f"Unknown filter '{filter_name}'. Available: {list(STANDARD_FILTERS.keys())}")
        
    cen_nm, fwhm_nm = STANDARD_FILTERS[filter_name]
    sigma_nm = fwhm_nm / 2.3548 # FWHM to Sigma
    
    # Simple Gaussian profile (un-normalized, peak at 1.0)
    transmission = np.exp(-0.5 * ((target_x_nm - cen_nm) / sigma_nm)**2)
    
    # Clip very small values to avoid numerical noise in integration
    transmission[transmission < 1e-5] = 0.0
    
    return transmission

def calculate_synthetic_ab_mag(x: au.Quantity, y: au.Quantity, filter_name: str) -> float:
    """
    Calculates the synthetic AB magnitude of a spectrum through a filter.
    Uses full unit conversion for accuracy.
    """
    # 1. Ensure we are working with standard units internally
    try:
        # Convert X to Angstroms for standard photometry math
        wave_aa = x.to(au.Angstrom)
        # Ensure Y is a flux density
        flux_flam = y.to(au.erg / (au.cm**2 * au.s * au.Angstrom), 
                         equivalencies=au.spectral_density(wave_aa))
    except au.UnitConversionError as e:
        # If units are missing/arbitrary, we can't do rigorous photometry
        if y.unit == au.dimensionless_unscaled:
             logging.warning("Spectrum has arbitrary units. Assuming erg/cm2/s/A for calculation.")
             flux_flam = y.value * (au.erg / (au.cm**2 * au.s * au.Angstrom))
        else:
             logging.error(f"Photometry failed: incompatible units.\n{e}")
             raise
        
    # 2. Get filter profile on the same grid
    x_nm_vals = x.to_value(au.nm)
    T = get_filter_transmission(filter_name, x_nm_vals)
    
    if np.sum(T) == 0:
        raise ValueError(f"Filter '{filter_name}' does not overlap with the spectrum.")

    # 3. Integration (Weighted average of F_nu)
    # Convert f_lambda to f_nu
    flux_fnu = flux_flam.to(au.Jy, equivalencies=au.spectral_density(wave_aa))
    nu = wave_aa.to(au.Hz, equivalencies=au.spectral())
    
    # Integrate: \int f_nu * T / nu dnu  /  \int T / nu dnu
    numerator = np.trapz(flux_fnu.value * T / nu.value, x=nu.value)
    denominator = np.trapz(T / nu.value, x=nu.value)
    
    avg_fnu_jy = np.abs(numerator / denominator)
    
    if avg_fnu_jy <= 0:
         return np.nan

    # 4. Convert Jy to AB mag (0 AB = 3631 Jy)
    ab_mag = -2.5 * np.log10(avg_fnu_jy / 3631.0)
    
    return ab_mag