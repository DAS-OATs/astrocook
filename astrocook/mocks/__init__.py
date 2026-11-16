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
    spectrum = generate(truth)          # Spectrum
    truth.to_json("mock_truth.json")    # ground truth, for later comparison

The cases of the robustness grid of phase 7 (PR 7a.1): BAL troughs
(:class:`BALTrough`), a richer set of emission lines with shifts and
asymmetries (:data:`QSO_EMISSION_LINES_FULL`), tellurics
(:class:`TelluricModel`, :mod:`astrocook.mocks.telluric`), a proximity zone
of the forest and associated absorbers (:func:`lya_forest`,
:func:`associated_absorbers`), and the truth map, what absorbs each pixel
(:func:`truth_map`).
"""
from astrocook.mocks.truth import (
    BALTrough, EmissionLine, MockComponent, MockTruth, QSO_EMISSION_LINES,
    QSO_EMISSION_LINES_FULL, TelluricModel,
)
from astrocook.mocks.generate import (
    continuum, generate, make_grid, noise_sigma, normalized_flux,
    spike_indices,
)
from astrocook.mocks.physics import (
    LYMAN_SERIES, bal_optical_depth, convolve_velocity,
    lyman_limit_optical_depth, optical_depth,
)
from astrocook.mocks.populations import (
    DLA_LOGN, LLS_LOGN, SUBDLA_LOGN, associated_absorbers,
    expected_forest_lines, lya_forest, mean_transmission_fg08, metal_systems,
    proximity_omega, strong_hi_systems,
)
from astrocook.mocks.telluric import O2_BANDS_NM, telluric_optical_depth
from astrocook.mocks.truthmap import (
    CATEGORIES, TruthMap, component_category, truth_map,
)

__all__ = [
    "BALTrough", "EmissionLine", "MockComponent", "MockTruth",
    "QSO_EMISSION_LINES", "QSO_EMISSION_LINES_FULL", "TelluricModel",
    "continuum", "generate", "make_grid", "noise_sigma", "normalized_flux",
    "spike_indices",
    "LYMAN_SERIES", "optical_depth", "convolve_velocity",
    "lyman_limit_optical_depth", "bal_optical_depth",
    "DLA_LOGN", "LLS_LOGN", "SUBDLA_LOGN", "expected_forest_lines",
    "lya_forest", "mean_transmission_fg08", "metal_systems",
    "strong_hi_systems", "associated_absorbers", "proximity_omega",
    "O2_BANDS_NM", "telluric_optical_depth",
    "CATEGORIES", "TruthMap", "component_category", "truth_map",
]
