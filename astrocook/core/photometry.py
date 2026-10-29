# v2/photometry.py

import numpy as np
import astropy.units as au
import logging
from scipy.integrate import trapezoid

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

FLAM = au.erg / (au.cm**2 * au.s * au.Angstrom)
AB_ZERO_POINT_JY = 3631.0


def _to_flam(x: au.Quantity, y: au.Quantity):
    """The wavelengths in Angstrom and the flux in erg/cm2/s/A."""
    # 1. Ensure we are working with standard units internally
    try:
        # Convert X to Angstroms for standard photometry math
        wave_aa = x.to(au.Angstrom)
        # Ensure Y is a flux density
        flux_flam = y.to(FLAM, equivalencies=au.spectral_density(wave_aa))
    except au.UnitConversionError as e:
        # If units are missing/arbitrary, we can't do rigorous photometry
        if y.unit == au.dimensionless_unscaled:
             logging.warning("Spectrum has arbitrary units. Assuming erg/cm2/s/A for calculation.")
             flux_flam = y.value * FLAM
        else:
             logging.error(f"Photometry failed: incompatible units.\n{e}")
             raise
    return wave_aa, flux_flam


def _mean_fnu_jy(wave_aa: au.Quantity, flux_flam: au.Quantity,
                 filter_name: str) -> float:
    """The mean f_nu (Jy) of a spectrum through a filter, with its sign:
    linear in the flux."""
    # 2. Get filter profile on the same grid
    T = get_filter_transmission(filter_name, wave_aa.to_value(au.nm))
    
    if np.sum(T) == 0:
        raise ValueError(f"Filter '{filter_name}' does not overlap with the spectrum.")

    # 3. Integration (Weighted average of F_nu)
    # Convert f_lambda to f_nu
    flux_fnu = flux_flam.to(au.Jy, equivalencies=au.spectral_density(wave_aa))
    nu = wave_aa.to(au.Hz, equivalencies=au.spectral())
    
    # Integrate: \int f_nu * T / nu dnu  /  \int T / nu dnu
    # (scipy's trapezoid: np.trapz is gone from recent NumPy, and
    # np.trapezoid does not exist in NumPy 1.26, which Astrocook supports)
    numerator = trapezoid(flux_fnu.value * T / nu.value, x=nu.value)
    denominator = trapezoid(T / nu.value, x=nu.value)
    return numerator / denominator


def _ab_mag(avg_fnu_jy: float) -> float:
    if avg_fnu_jy <= 0:
         return np.nan
    # 4. Convert Jy to AB mag (0 AB = 3631 Jy)
    return -2.5 * np.log10(avg_fnu_jy / AB_ZERO_POINT_JY)


def calculate_synthetic_ab_mag(x: au.Quantity, y: au.Quantity, filter_name: str) -> float:
    """
    Calculates the synthetic AB magnitude of a spectrum through a filter.
    Uses full unit conversion for accuracy.
    """
    wave_aa, flux_flam = _to_flam(x, y)
    return _ab_mag(np.abs(_mean_fnu_jy(wave_aa, flux_flam, filter_name)))


def generate_calibration_curve(x: au.Quantity, y: au.Quantity, 
                               magnitudes_str: str) -> np.ndarray:
    """
    Parses a 'Filter=Mag,Filter=Mag' string and returns a correction curve
    (NumPy array over the x-axis) that brings the synthetic magnitudes of
    the spectrum to the targets.

    With one band the curve is a constant. With two or more it is linear
    between the effective wavelengths of the filters and constant beyond
    them; its values there are solved together, because each filter also
    sees the curve near the other filters (findings register F70: the
    factor of each band alone missed the targets by ~0.1 mag).
    """
    # 1. Parse Input
    targets = []
    try:
        for pair in magnitudes_str.split(','):
            filt, mag = pair.split('=')
            targets.append((filt.strip(), float(mag.strip())))
    except ValueError:
        raise ValueError("Invalid format. Use 'Filter=Mag, Filter=Mag'.")

    if not targets:
         raise ValueError("No valid magnitudes provided.")

    wave_aa, flux_flam = _to_flam(x, y)

    # 2. Calculate Scale Factors
    bands = []  # (effective wavelength, filter, target, factor of the band alone)

    for filter_name, target_mag in targets:
        # _mean_fnu_jy handles validation of filter_name
        current_mag = _ab_mag(np.abs(_mean_fnu_jy(wave_aa, flux_flam, filter_name)))
        
        if np.isnan(current_mag):
             logging.warning(f"Could not calculate mag for {filter_name}, skipping.")
             continue
             
        # Factor = 10^(-0.4 * (target - current))
        factor = 10**(-0.4 * (target_mag - current_mag))
        
        lambda_eff = STANDARD_FILTERS[filter_name][0]
        bands.append((lambda_eff, filter_name, target_mag, factor))
        logging.info(f"Calibration point: {filter_name} (@{lambda_eff:.0f}nm) factor={factor:.3e}")

    if not bands:
        raise RuntimeError("Could not determine any valid scale factors.")

    # 3. Generate Curve
    x_nm = x.to_value(au.nm)
    
    if len(bands) == 1:
        # Constant scaling
        return np.full_like(x_nm, bands[0][3])

    # Warping: the curve is sum_j f_j * hat_j(x), hat_j the linear
    # interpolation of a unit value at the j-th effective wavelength
    # (constant beyond the first and the last); the mean f_nu in each band
    # is linear in the f_j, so A f = b gives the f_j that reach the targets
    bands.sort()
    waves = np.array([b[0] for b in bands])
    single = np.array([b[3] for b in bands])
    n = len(bands)
    A = np.empty((n, n))
    for j in range(n):
        hat = np.interp(x_nm, waves, np.eye(n)[j])
        for i, b in enumerate(bands):
            A[i, j] = _mean_fnu_jy(wave_aa, flux_flam * hat, b[1])
    b_vec = AB_ZERO_POINT_JY * 10**(-0.4 * np.array([b[2] for b in bands]))
    try:
        factors = np.linalg.solve(A, b_vec)
    except np.linalg.LinAlgError:
        factors = None
    if factors is None or not np.all(np.isfinite(factors)) or np.any(factors <= 0):
        # E.g. the same filter twice, or targets that only a curve with
        # negative values would reach: the factors of the bands alone
        logging.warning("The magnitudes cannot be matched together by a positive "
                        "curve: the factor of each band is used alone, and the "
                        "magnitudes are only approximately reached.")
        factors = single
    for b, f in zip(bands, factors):
        logging.info(f"Calibration curve at {b[1]} (@{b[0]:.0f}nm): {f:.3e}")
    return np.interp(x_nm, waves, factors)
