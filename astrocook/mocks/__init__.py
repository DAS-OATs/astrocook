"""Synthetic quasar spectra with known parameters (mocks).

The generator is independent of the V2 fitting code: it has its own Voigt
profile, convolution and pixel integration, so that the fitting code can be
tested against it. It shares only the atomic data table
(:data:`astrocook.core.atomic_data.ATOM_DATA`).

Typical use::

    from astrocook.mocks import MockComponent, MockTruth, generate

    truth = MockTruth(
        components=[MockComponent("CIV", z=2.97, logN=13.5, b=10.0)],
        lam_min_nm=612.0, lam_max_nm=620.0, dv_kms=2.5,
        resol=45000.0, snr=50.0, z_em=3.1, seed=1)
    spectrum = generate(truth)          # SpectrumV2
    truth.to_json("mock_truth.json")    # ground truth, for later comparison
"""
from astrocook.mocks.truth import (
    EmissionLine, MockComponent, MockTruth, QSO_EMISSION_LINES,
)
from astrocook.mocks.generate import (
    continuum, generate, make_grid, normalized_flux, spike_indices,
)
from astrocook.mocks.physics import optical_depth, convolve_velocity
from astrocook.mocks.populations import (
    expected_forest_lines, lya_forest, mean_transmission_fg08, metal_systems,
)

__all__ = [
    "EmissionLine", "MockComponent", "MockTruth", "QSO_EMISSION_LINES",
    "continuum", "generate", "make_grid", "normalized_flux", "spike_indices",
    "optical_depth", "convolve_velocity",
    "expected_forest_lines", "lya_forest", "mean_transmission_fg08",
    "metal_systems",
]
