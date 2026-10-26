"""Generation of mock spectra."""
from typing import Tuple

import numpy as np
import astropy.units as au

from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, SpectrumDataV2
from astrocook.mocks.physics import (
    C_KMS, FWHM_TO_SIGMA, convolve_velocity, is_hydrogen,
    lyman_limit_optical_depth, optical_depth,
)
from astrocook.mocks.truth import MockTruth


def make_grid(lam_min_nm: float, lam_max_nm: float,
              dv_kms: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pixel centres and edges (nm), uniform in velocity.

    Returns
    -------
    x, xmin, xmax : np.ndarray
        Centres and edges; the edges of adjacent pixels coincide.
    """
    step = np.log1p(dv_kms / C_KMS)
    n = int(np.floor(np.log(lam_max_nm / lam_min_nm) / step))
    edges = lam_min_nm * np.exp(step * np.arange(n + 1))
    x = np.sqrt(edges[:-1] * edges[1:])
    return x, edges[:-1], edges[1:]


def normalized_flux(truth: MockTruth) -> Tuple[np.ndarray, np.ndarray]:
    """Noiseless normalised flux, convolved and averaged over each pixel.

    The model is computed on ``oversample`` sub-pixels per pixel, convolved
    with the instrumental profile on that grid, then averaged over the
    sub-pixels of each pixel.

    Returns
    -------
    x_nm, flux : np.ndarray
    """
    x, xmin, _ = make_grid(truth.lam_min_nm, truth.lam_max_nm, truth.dv_kms)
    k = int(truth.oversample)
    step = np.log1p(truth.dv_kms / C_KMS)
    sub = np.repeat(xmin, k) * np.exp(step * (np.tile(np.arange(k), len(x)) + 0.5) / k)

    tau = np.zeros_like(sub)
    for comp in truth.components:
        tau += optical_depth(sub * 10.0, comp.series, comp.z, comp.logN,
                             comp.b, comp.btur)
        if truth.lyman_limit and is_hydrogen(comp.series):
            tau += lyman_limit_optical_depth(sub * 10.0, comp.z, comp.logN)
    flux_sub = convolve_velocity(np.exp(-tau), truth.dv_kms / k, truth.resol)
    return x, flux_sub.reshape(len(x), k).mean(axis=1)


def continuum(truth: MockTruth, x_nm: np.ndarray) -> np.ndarray:
    """Power law plus the Gaussian emission lines of ``truth``."""
    power_law = truth.cont_amplitude * (x_nm / truth.cont_ref_nm) ** truth.cont_slope
    bumps = np.zeros_like(x_nm)
    for line in truth.emission_lines:
        centre = line.rest_nm * (1.0 + truth.z_em)
        sigma = line.fwhm_kms / C_KMS * centre * FWHM_TO_SIGMA
        bumps += line.amplitude * np.exp(-0.5 * ((x_nm - centre) / sigma) ** 2)
    return power_law * (1.0 + bumps)


def noise_sigma(truth: MockTruth, cont: np.ndarray,
                true_flux: np.ndarray) -> np.ndarray:
    """Error per pixel of ``truth`` (see :class:`MockTruth`, ``snr`` and
    ``noise_floor``), from the continuum and the noiseless flux."""
    sigma_cont = cont / truth.snr
    if truth.noise_floor is None:
        return sigma_cont
    a = float(truth.noise_floor)
    if not 0.0 <= a <= 1.0:
        raise ValueError(f"noise_floor must be between 0 and 1, not {a}")
    ratio = np.clip(true_flux / cont, 0.0, None)
    return sigma_cont * np.sqrt(a + (1.0 - a) * ratio)


def spike_indices(truth: MockTruth, n_pixels: int) -> np.ndarray:
    """Pixel indices of the spikes of ``truth`` (reproducible from the seed)."""
    rng = np.random.default_rng([truth.seed, 1])
    return np.sort(rng.choice(n_pixels, size=truth.n_spikes, replace=False))


def generate(truth: MockTruth, truth_columns: bool = True) -> SpectrumV2:
    """Mock spectrum for ``truth``.

    The spectrum has the true continuum (with emission lines) in the
    ``cont`` column and the noiseless flux in ``mock_flux``; ``z_em`` and the resolving power are
    set (the latter also in the metadata, so that it survives a save).

    With ``truth_columns=False`` the two columns are left out, as for an
    observed spectrum: the recipes would otherwise use the true continuum
    (in the MCP pilot an agent fitted on it). The flux is the same.
    """
    x, xmin, xmax = make_grid(truth.lam_min_nm, truth.lam_max_nm, truth.dv_kms)
    _, flux = normalized_flux(truth)
    cont = continuum(truth, x)
    true_flux = cont * flux

    if truth.snr is None:
        y = true_flux.copy()
        dy = np.zeros_like(x)
    else:
        rng = np.random.default_rng(truth.seed)
        dy = noise_sigma(truth, cont, true_flux)
        y = true_flux + rng.standard_normal(len(x)) * dy
    if truth.n_spikes:
        spikes = spike_indices(truth, len(x))
        y[spikes] += truth.spike_amplitude * (dy[spikes] if truth.snr else cont[spikes])

    unit = au.dimensionless_unscaled
    aux = {
        "cont": DataColumnV2(cont, unit, description="True continuum"),
        "mock_flux": DataColumnV2(true_flux, unit, description="Noiseless flux"),
    } if truth_columns else {}
    if truth.resol is not None:
        # The resolution goes into the resol column only (#526)
        aux["resol"] = DataColumnV2(np.full(len(x), float(truth.resol)), unit,
                                    description="Resolving power")
    data = SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(xmin, au.nm),
        xmax=DataColumnV2(xmax, au.nm),
        y=DataColumnV2(y, unit), dy=DataColumnV2(dy, unit),
        aux_cols=aux, meta={"MOCK": True}, z_em=truth.z_em)
    return SpectrumV2(data)
