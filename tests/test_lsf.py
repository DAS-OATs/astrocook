"""Exact convolution with the line spread function (#526, F9, F10):
``PixelLSF`` against the analytic result, and the model of the fitter
against the exact flux of mocks, with a flat and a varying resolution."""
import astropy.units as au
import numpy as np
import pytest
from scipy.special import erf

from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import (
    ComponentDataV2, DataColumnV2, SpectrumDataV2, SystemListDataV2,
)
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.lsf import C_KMS, FWHM_TO_SIGMA, PixelLSF
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.mocks import MockComponent
from astrocook.mocks.physics import convolve_velocity, optical_depth


def _pixels(width_kms, n=200, lam0=6000.0):
    edges = lam0 * np.exp(np.arange(n + 1) * width_kms / C_KMS)
    return edges[:-1], edges[1:]


def _gaussian_line(lsf, resol, b, depth=0.9):
    """Model of a Gaussian line (Doppler parameter b) through ``lsf`` and
    its exact value: the pixel mean of the line convolved with the LSF."""
    centre = 0.5 * (lsf.v_lo[100] + lsf.v_hi[100]) + 0.37 * (lsf.v_hi[100] - lsf.v_lo[100])
    s_line = b / np.sqrt(2.0)
    flux = 1.0 - depth * np.exp(-0.5 * ((lsf.v_fine - centre) / s_line) ** 2)
    model = lsf.apply(flux, resol)
    sigma = 0.0 if resol is None else C_KMS / np.asarray(resol) * FWHM_TO_SIGMA
    s_tot = np.sqrt(s_line**2 + sigma**2)
    amp = depth * s_line / s_tot
    lo = (lsf.v_lo - centre) / (np.sqrt(2.0) * s_tot)
    hi = (lsf.v_hi - centre) / (np.sqrt(2.0) * s_tot)
    exact = 1.0 - amp * s_tot * np.sqrt(np.pi / 2.0) * (erf(hi) - erf(lo)) / (lsf.v_hi - lsf.v_lo)
    return model, exact


@pytest.mark.parametrize("width, resol", [(2.66, 45000.0), (13.5, 8900.0)])
@pytest.mark.parametrize("b", [1.0, 12.0])
def test_resolved_kernel_is_exact(width, resol, b):
    lsf = PixelLSF(*_pixels(width), r_min=resol, r_max=resol, b_min=b)
    model, exact = _gaussian_line(lsf, resol, b)
    np.testing.assert_allclose(model, exact, atol=1e-6)


def test_resolution_per_pixel():
    # Each pixel takes its own kernel: R from 7000 to 11000 in the window
    resol = np.linspace(7000.0, 11000.0, 200)
    lsf = PixelLSF(*_pixels(13.5), r_min=resol.min(), r_max=resol.max(), b_min=12.0)
    model, exact = _gaussian_line(lsf, resol, 12.0)
    np.testing.assert_allclose(model, exact, atol=1e-6)


@pytest.mark.parametrize("width", [1.0, 2.66])
def test_no_resolution_is_the_pixel_mean(width):
    # Without a resolution the kernel is the pixel: the samples at the cell
    # centres would shift each pixel by up to half a step (3e-3 here)
    lsf = PixelLSF(*_pixels(width), b_min=12.0)
    model, exact = _gaussian_line(lsf, None, 12.0)
    np.testing.assert_allclose(model, exact, atol=5e-4)


# --- The fitter against mocks ---------------------------------------------------

COMPONENTS = [MockComponent("CIV", 2.97, 13.6, 12.0), MockComponent("CIV", 2.97015, 13.1, 3.0),
              MockComponent("Ly_a", 4.059, 14.5, 25.0), MockComponent("Ly_a", 4.12, 13.0, 6.0)]


def _mock(resol_at, dv_kms, k=20, lam_min=600.0, lam_max=640.0):
    """Exact normalised flux on pixels uniform in velocity: the model on k
    sub-pixels, convolved with the resolution of each pixel
    (``resol_at(x_nm)``), averaged over the pixel. Written with the
    physics of the mocks, independent of the fitter."""
    step = np.log1p(dv_kms / C_KMS)
    n = int(np.log(lam_max / lam_min) / step)
    xmin = lam_min * np.exp(step * np.arange(n))
    xmax = xmin * np.exp(step)
    x = np.sqrt(xmin * xmax)
    sub = np.repeat(xmin, k) * np.exp(step * (np.tile(np.arange(k), n) + 0.5) / k)
    tau = sum(optical_depth(sub * 10.0, c.series, c.z, c.logN, c.b) for c in COMPONENTS)
    resol = resol_at(x)
    flux = convolve_velocity(np.exp(-tau), dv_kms / k, np.repeat(resol, k))
    flux = flux.reshape(n, k).mean(axis=1)
    one = np.ones(n)
    unit = au.dimensionless_unscaled
    spec = SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(xmin, au.nm),
        xmax=DataColumnV2(xmax, au.nm), y=DataColumnV2(flux, unit),
        dy=DataColumnV2(0.01 * one, unit),
        aux_cols={"cont": DataColumnV2(one, unit), "resol": DataColumnV2(resol, unit)},
        z_em=4.2))
    return spec, flux


def _model(spec):
    systs = SystemListV2(SystemListDataV2(components=[
        ComponentDataV2(id=i + 1, z=c.z, dz=None, logN=c.logN, dlogN=None, b=c.b,
                        db=None, series=c.series)
        for i, c in enumerate(COMPONENTS)]))
    return VoigtFitterV2(spec, systs).compute_model_flux()[1]


# Before #526 the model of these mocks was off by up to 5.9e-3 (R 45000)
# and 5.7e-2 (R 8900, flat or varying); with #526 but before F61 (damping
# parameter (1 + z) too large), by 1.3e-3 to 1.8e-3. Now within 3e-5
@pytest.mark.parametrize("resol_at", [
    lambda x: np.full(len(x), 45000.0),
    lambda x: np.full(len(x), 8900.0),
    lambda x: 8900.0 * (x / 620.0) ** 4,       # 7000 to 10900 in the window
], ids=["R45000", "R8900", "R8900_varying"])
def test_fitter_model_matches_the_exact_flux(resol_at):
    r_mid = resol_at(np.array([620.0]))[0]
    spec, truth = _mock(resol_at, C_KMS / r_mid / 2.5)
    np.testing.assert_allclose(_model(spec), truth, atol=1e-4)
