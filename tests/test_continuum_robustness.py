"""Robustness of the automatic continuum (#465) on mocks with a known
continuum: upward spikes, low S/N, and the masked running mean of
``find_absorbed_regions``."""
import dataclasses
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "benchmarks"))

import validate_mocks as vm  # noqa: E402

from astrocook.core.session import SessionV2  # noqa: E402
from astrocook.core.spectrum_operations import _masked_window_mean  # noqa: E402
from astrocook.mocks import generate  # noqa: E402


def _continuum_ratio(snr, resol, n_spikes_fraction=0.0, amplitude=30.0, seed=2):
    """Median ratio of the estimated to the true continuum on the red side
    and in the forest, for a mock of the benchmark (z_em = 3.5)."""
    truth, _, _ = vm.make_truth(snr, resol, 3.5, seed)
    n = len(generate(truth, truth_columns=False).x)
    truth = dataclasses.replace(truth, n_spikes=int(n_spikes_fraction * n),
                                spike_amplitude=amplitude)
    spec = generate(truth, truth_columns=False)
    cont = SessionV2("m", None, spec=spec).continuum.estimate_auto(
        **vm.schema_params("continuum", "estimate_auto"))
    metrics = vm.continuum_metrics(truth, cont.spec)
    return metrics["red"]["median"], metrics["forest"]["median"]


def test_spikes_do_not_pull_the_continuum_up():
    # 1% of the pixels with spikes of 30 sigma, low resolution and S/N as
    # in NIR spectra: before #465 the continuum came out ~7 times too high
    red, forest = _continuum_ratio(10.0, 8900.0, n_spikes_fraction=0.01)
    assert red == pytest.approx(1.0, abs=0.02)
    assert 0.8 < forest < 1.0


def test_low_snr_red_side():
    # At S/N 3 the test of each pixel (y < 3 dy) dropped the lower half of
    # the noise: continuum 26% too high (the mean of a half Gaussian)
    red, _ = _continuum_ratio(3.0, 8900.0)
    assert red == pytest.approx(1.0, abs=0.05)


def test_masked_window_mean_is_exact_at_the_edges():
    values = np.arange(10.0)
    valid = np.ones(10, dtype=bool)
    mean = _masked_window_mean(values, valid, 2)
    assert mean[0] == pytest.approx(1.0)       # (0 + 1 + 2) / 3
    assert mean[5] == pytest.approx(5.0)
    assert mean[-1] == pytest.approx(8.0)      # (7 + 8 + 9) / 3


def test_masked_window_mean_skips_invalid_pixels():
    values = np.array([1.0, 100.0, 3.0, np.nan, 5.0])
    valid = np.array([True, False, True, False, True])
    mean = _masked_window_mean(values, valid, 1)
    assert mean[1] == pytest.approx(2.0)       # (1 + 3) / 2
    assert mean[3] == pytest.approx(4.0)       # (3 + 5) / 2
    assert np.isnan(_masked_window_mean(values, np.zeros(5, dtype=bool), 1)).all()
