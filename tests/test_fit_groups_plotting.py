"""What the tools of inspection of the MCP server compute, in the package
for scripts too (phase 6-ter, PR 3): the fit groups of Refit All
(``SystemList.fit_groups``), the reading of a range
(``Spectrum.columns_in_range``), the plots without the GUI
(``astrocook.plotting``)."""
import subprocess
import sys

import numpy as np
import pytest

from astrocook.session.session import Session
from astrocook.core.structures import ComponentData, SystemListData
from astrocook.core.system_list import SystemList
from astrocook.mocks import MockComponent, MockTruth, generate

C_KMS = 299792.458


def _comp(i, z, series="CIV"):
    return ComponentData(id=i, z=z, dz=None, logN=13.0, dlogN=None, b=10.0, db=None,
                           series=series, uuid=f"u{i}")


def _systs(*zs):
    return SystemList(SystemListData(
        components=[_comp(i, z) for i, z in enumerate(zs)]))


# --- fit_groups --------------------------------------------------------------

def test_fit_groups_blends_and_isolated():
    # Two CIV 10 km/s apart are one group; one far away its own
    dz = 10 / C_KMS * 2.9
    systs = _systs(1.95, 1.9, 1.9 + dz)
    assert systs.fit_groups(["u0", "u1", "u2"], max_depth=2) == [["u1", "u2"], ["u0"]]


def test_fit_groups_split_by_max_size_and_seeds():
    dz = 10 / C_KMS * 2.9
    systs = _systs(*(1.9 + k * dz for k in range(5)))
    groups = systs.fit_groups([f"u{k}" for k in range(5)], max_depth=10, max_size=2)
    assert groups == [["u0", "u1"], ["u2", "u3"], ["u4"]]
    # Only the seeds start a group, but the group takes in its neighbours
    assert systs.fit_groups(["u4"], max_depth=10) == [[f"u{k}" for k in range(5)]]
    # Unknown seeds are ignored
    assert systs.fit_groups(["nope"], max_depth=2) == []


def test_refit_all_fits_the_groups_of_fit_groups(monkeypatch):
    # refit_all builds its groups with fit_groups (moved there from the
    # recipe; the numbers are covered by the snapshots of the full chain)
    truth = MockTruth(components=[MockComponent("CIV", 1.9, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5,
                      resol=45000.0, snr=50.0, z_em=2.0, seed=1)
    session = Session("mock", None, spec=generate(truth))
    session = session.absorbers.add_component(series="CIV", z=1.9)
    seen = []
    original = SystemList.fit_groups

    def spy(self, *args, **kwargs):
        out = original(self, *args, **kwargs)
        seen.append(out)
        return out

    monkeypatch.setattr(SystemList, "fit_groups", spy)
    session.absorbers.refit_all()
    assert seen and seen[0] == [[session.systs.components[0].uuid]]


# --- columns_in_range ----------------------------------------------------------

@pytest.fixture(scope="module")
def spec():
    truth = MockTruth(components=[MockComponent("CIV", 1.9, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5,
                      resol=45000.0, snr=50.0, z_em=2.0, seed=1)
    return generate(truth)


def test_columns_in_range(spec):
    x, y, dy = spec.x.value, spec.y.value, spec.dy.value
    part = spec.columns_in_range(["y", "dy"], 449, 450)
    sel = (x >= 449) & (x <= 450)
    assert part["binned_by"] == 1 and part["n_pixels"] == sel.sum()
    assert part["columns"]["y"] == pytest.approx(y[sel])
    # The axis as λ when it is not asked (phase 6-quater, PR 8)
    assert part["units"]["λ"] == "nm" and part["axis"] == "λ"
    whole = spec.columns_in_range(["y", "dy"], max_points=100)
    n = whole["binned_by"]
    assert n == int(np.ceil(len(x) / 100)) and len(whole["columns"]["λ"]) <= 100
    assert whole["columns"]["y"][0] == pytest.approx(np.mean(y[:n]))
    # The error of the mean
    assert whole["columns"]["dy"][0] == pytest.approx(np.sqrt(np.sum(dy[:n] ** 2)) / n)


def test_columns_in_range_errors(spec):
    with pytest.raises(KeyError, match="No column nope. The columns: λ"):
        spec.columns_in_range(["nope"])
    with pytest.raises(ValueError, match="Empty range"):
        spec.columns_in_range(["y"], 450, 449)
    with pytest.raises(ValueError, match="No pixel between 900"):
        spec.columns_in_range(["y"], 900, 901)


# --- plotting ----------------------------------------------------------------

def test_plot_spectrum(spec):
    from astrocook.plotting import figure_png, plot_spectrum
    session = Session("mock", None, spec=spec).absorbers.add_component(series="CIV", z=1.9)
    fig = plot_spectrum(session, 445, 452, normalized=True, palette="colorblind")
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert labels == ["flux", "error", "continuum", "model"]
    assert fig.axes[0].get_ylabel() == "normalized flux"
    # A tick per line of the doublet, labelled with the series
    assert [t.get_text() for t in fig.axes[0].texts] == ["CIV", "CIV"]
    assert figure_png(fig)[:8] == b"\x89PNG\r\n\x1a\n"
    with pytest.raises(ValueError, match="No pixel"):
        plot_spectrum(session, 900, 901)


def test_plotting_needs_no_gui():
    # The headless image of the server has no PySide6: plotting must not
    # import Qt (nor pyplot, which would choose a backend)
    code = ("import sys, astrocook.plotting; "
            "assert 'PySide6' not in sys.modules, 'PySide6'; "
            "assert 'matplotlib.pyplot' not in sys.modules, 'pyplot'")
    subprocess.run([sys.executable, "-c", code], check=True)


def test_the_gui_uses_the_same_palettes():
    from astrocook import plotting
    from astrocook.gui import plot_style
    assert plot_style.PLOT_PALETTES is plotting.PLOT_PALETTES
    assert plot_style.ERROR_COLOR == plotting.ERROR_COLOR
