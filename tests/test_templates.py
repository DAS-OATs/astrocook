"""Templates moved out of ``astrocook.legacy`` (phase 3).

``qso_composite`` must give the same arrays as V1 ``vars.qso_composite``;
the continuum fitted with the template is covered by the snapshot
``continuum_j0103_fit_template`` (``test_snapshots_continuum.py``).
"""
import numpy as np
import pytest

from astrocook.core.templates import qso_composite
from astrocook.legacy import vars as vars_v1


def test_qso_composite_matches_v1():
    x, y = qso_composite()
    np.testing.assert_array_equal(x, np.asarray(vars_v1.qso_composite["x"]))
    np.testing.assert_array_equal(y, np.asarray(vars_v1.qso_composite["y"]))


def test_qso_composite_is_read_only():
    x, y = qso_composite()
    assert x is qso_composite()[0]  # read once
    with pytest.raises(ValueError):
        x[0] = 0.0
    with pytest.raises(ValueError):
        y[0] = 0.0
