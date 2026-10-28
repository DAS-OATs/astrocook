"""Synthetic photometry and flux calibration (``core/photometry.py``,
``SpectrumV2.calibrate``, *Flux > Calibrate*).

No test reached this code until phase 5, when ``np.trapz`` was found to be
gone from recent NumPy: the calibration failed at every call.
"""
import astropy.units as au
import numpy as np
import pytest

from astrocook.core.photometry import calculate_synthetic_ab_mag, generate_calibration_curve
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, SpectrumDataV2

FLAM = au.erg / (au.cm**2 * au.s * au.Angstrom)
FLUX_NM = au.erg / au.cm**2 / au.s / au.nm


def _flat_fnu(jy=1e-3):
    """A spectrum with a flat f_nu, in erg/cm2/s/nm, 400-700 nm."""
    x = np.linspace(400.0, 700.0, 3001) * au.nm
    fnu = np.full(x.size, jy) * au.Jy
    return x, fnu.to(FLUX_NM, equivalencies=au.spectral_density(x))


@pytest.mark.parametrize("band", ["SDSS_g", "SDSS_r", "Johnson_V"])
def test_ab_magnitude_of_a_flat_spectrum(band):
    x, y = _flat_fnu(1e-3)
    # 0 AB = 3631 Jy: 1 mJy is 16.40
    assert calculate_synthetic_ab_mag(x, y, band) == pytest.approx(
        -2.5 * np.log10(1e-3 / 3631.0), abs=1e-9)


def test_filter_outside_the_spectrum():
    x, y = _flat_fnu()
    with pytest.raises(ValueError, match="does not overlap"):
        calculate_synthetic_ab_mag(x, y, "SDSS_z")


def test_one_band_calibration_is_exact():
    x, y = _flat_fnu()
    curve = generate_calibration_curve(x, y, "SDSS_g=17.5")
    np.testing.assert_allclose(curve, curve[0])  # a constant
    assert calculate_synthetic_ab_mag(x, y * curve, "SDSS_g") == pytest.approx(17.5, abs=1e-9)


def test_two_band_calibration():
    # The factors of the two bands are interpolated between their effective
    # wavelengths: the calibrated spectrum is ~0.1 mag from the targets
    # (findings register F70)
    x, y = _flat_fnu()
    curve = generate_calibration_curve(x, y, "SDSS_g=17.5, SDSS_r=17.0")
    assert calculate_synthetic_ab_mag(x, y * curve, "SDSS_g") == pytest.approx(17.404, abs=1e-3)
    assert calculate_synthetic_ab_mag(x, y * curve, "SDSS_r") == pytest.approx(17.083, abs=1e-3)


@pytest.mark.parametrize("magnitudes, error", [("SDSS_g 17.5", "Invalid format"),
                                               ("Nope=17.5", "Unknown filter")])
def test_invalid_magnitudes(magnitudes, error):
    x, y = _flat_fnu()
    with pytest.raises(ValueError, match=error):
        generate_calibration_curve(x, y, magnitudes)


def test_spectrum_calibrate():
    x, y = _flat_fnu()
    data = SpectrumDataV2(
        x=DataColumnV2(x.value, au.nm), xmin=DataColumnV2(x.value - 0.05, au.nm),
        xmax=DataColumnV2(x.value + 0.05, au.nm), y=DataColumnV2(y.value, FLUX_NM),
        dy=DataColumnV2(0.01 * y.value, FLUX_NM),
        aux_cols={"cont": DataColumnV2(y.value, FLUX_NM)}, z_em=3.0)
    spec = SpectrumV2(data)
    new = spec.calibrate("SDSS_g=17.5")
    assert new.y.unit == FLAM and new.get_column("cont").unit == FLAM
    # The values are multiplied by the curve; the unit is relabelled
    curve = generate_calibration_curve(spec.x, spec.y, "SDSS_g=17.5")
    np.testing.assert_array_equal(new.y.value, spec.y.value * curve)
    assert new.history[-1] == "Flux calibrated with: SDSS_g=17.5"
