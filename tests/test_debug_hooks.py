"""The core does not depend on the GUI (findings register F6).

``core/spectrum_operations.py`` imported ``gui.debug_utils`` for the debug
plots of the doublet rating, so importing the core imported PySide6 too.
The core now hands the plot data to a handler set by the GUI
(``core/debug_hooks.py``).
"""
import os
import subprocess
import sys

import numpy as np
import pytest

from astrocook.core import debug_hooks
from astrocook.core.atomic_data import xem_d
from astrocook.core.spectrum_operations import rate_doublet_candidate
from astrocook.mocks import MockComponent, MockTruth, generate


def test_import_does_not_load_the_gui():
    code = ("import sys, astrocook, astrocook.core.session, "
            "astrocook.core.spectrum_operations, astrocook.recipes.absorbers, "
            "astrocook.fitting.voigt_fitter, astrocook.io.loaders, astrocook.mocks\n"
            "print(sorted(m for m in sys.modules "
            "if m.split('.')[0] == 'PySide6' or m.startswith('astrocook.gui')))")
    env = dict(os.environ, MPLBACKEND="Agg")
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=120, env=env)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "[]"


@pytest.fixture
def civ():
    """A CIV doublet on a noiseless mock, with the masks of its two lines
    (the redder first, as the identification passes them)."""
    z = 2.95
    truth = MockTruth(components=[MockComponent("CIV", z, 13.5, 10.0)],
                      lam_min_nm=608.0, lam_max_nm=616.0, dv_kms=2.5, snr=None)
    spec = generate(truth)
    x = spec.x.value
    masks = [np.abs(x - xem_d[t].to_value("nm") * (1 + z)) < 0.05
             for t in ("CIV_1550", "CIV_1548")]
    return spec, masks, z


@pytest.fixture
def no_handler(monkeypatch):
    monkeypatch.setattr(debug_hooks, "_debug_plot_handler", None)


def test_debug_rating_reaches_the_handler(civ, no_handler):
    spec, (mask_1, mask_2), z = civ
    plots = []
    debug_hooks.set_debug_plot_handler(plots.append)
    score = rate_doublet_candidate(spec, mask_1, mask_2, "CIV", z, debug_rating=True)
    assert score > 0.9
    (plot,) = plots
    assert plot["title"].startswith("CIV z=2.9500")
    assert set(plot) == {"title", "v_compare", "Y_data", "Y_model"}
    assert len(plot["v_compare"]) == len(plot["Y_data"]) == len(plot["Y_model"])
    # The debug plot does not change the rating
    assert score == rate_doublet_candidate(spec, mask_1, mask_2, "CIV", z)


def test_debug_rating_without_handler(civ, no_handler):
    spec, (mask_1, mask_2), z = civ
    assert (rate_doublet_candidate(spec, mask_1, mask_2, "CIV", z, debug_rating=True)
            == rate_doublet_candidate(spec, mask_1, mask_2, "CIV", z))
