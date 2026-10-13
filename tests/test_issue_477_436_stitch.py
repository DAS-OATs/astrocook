"""Issues #477 and #436: combine two arms (e.g. X-Shooter VIS and NIR).

#436 asked for a way to combine the arms, scaling one to the other; #477
reported that with two arms the recipe asked for two cut-off wavelengths.
With two arms, ``equalize_and_stitch`` must accept exactly one cut-off,
scale the redder arm to the bluer one and join them without overlap.
"""
from types import SimpleNamespace

import numpy as np
import astropy.units as au
import pytest

from astrocook.core.session import SessionV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, SpectrumDataV2
from astrocook.recipes.edit import RecipeEditV2

CUT_NM = 1000.0
NIR_MISCALIBRATION = 0.7  # the NIR arm is 30% too faint
FLUX_UNIT = au.erg / (au.cm**2 * au.s * au.nm)


def _true_flux(x):
    return (x / 1000.0) ** -1.5


def _arm(x, scale, seed):
    rng = np.random.default_rng(seed)
    dx = x[1] - x[0]
    y = scale * _true_flux(x) * (1 + 0.005 * rng.standard_normal(len(x)))
    return SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - dx / 2, au.nm),
        xmax=DataColumnV2(x + dx / 2, au.nm), y=DataColumnV2(y, FLUX_UNIT),
        dy=DataColumnV2(0.005 * np.abs(y), FLUX_UNIT)))


@pytest.fixture
def arms():
    vis = _arm(np.arange(900.0, 1020.0, 0.02), 1.0, seed=1)
    nir = _arm(np.arange(990.0, 1100.0, 0.06), NIR_MISCALIBRATION, seed=2)
    # Minimal stand-in for the GUI: the recipe looks up the other sessions
    # by display name among the open session histories.
    gui = SimpleNamespace(session_histories=[
        SimpleNamespace(display_name="NIR",
                        current_state=SimpleNamespace(spec=nir))])
    session = SessionV2(name="VIS", gui=gui, spec=vis)
    return session, vis, nir


@pytest.mark.parametrize("kwargs", [
    dict(equalize_ranges="990-1020"),
    dict(equalize_ranges="auto"),
    dict(manual_factors=f"{1 / NIR_MISCALIBRATION}"),
])
def test_two_arms_one_cut(arms, kwargs):
    session, vis, nir = arms
    result = RecipeEditV2(session).equalize_and_stitch(
        other_sessions="NIR", stitch_wavelengths=f"{CUT_NM}", **kwargs)
    assert result != 0

    x = result.spec.x.value
    y = result.spec.y.value
    # Joined without overlap, covering both arms
    assert np.all(np.diff(x) > 0)
    assert x.min() == pytest.approx(vis.x.value.min())
    assert x.max() == pytest.approx(nir.x.value.max())

    blue = x <= CUT_NM
    # Blue side: the VIS data, untouched
    np.testing.assert_array_equal(y[blue], vis.y.value[vis.x.value <= CUT_NM])
    # Red side: the NIR data, rescaled to the VIS flux level
    ratio = np.median(y[~blue] / _true_flux(x[~blue]))
    assert ratio == pytest.approx(1.0, abs=0.01)


def test_two_arms_two_cuts_is_rejected(arms):
    session, _, _ = arms
    result = RecipeEditV2(session).equalize_and_stitch(
        other_sessions="NIR", stitch_wavelengths="1000, 1010")
    assert result == 0


def test_gui_dialog_asks_for_one_cut_with_two_arms(qtbot, monkeypatch, tmp_path):
    # GUI side of #477: with two arms the dialog must ask for one value and
    # pass the other session itself (the field is hidden from the user)
    import copy
    from astrocook.core.session_manager import SessionHistory
    from astrocook.core.structures import HistoryLogV2
    import astrocook.recipes.edit as edit

    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    # The handler edits the recipe schema in place: restore it afterwards
    monkeypatch.setattr(edit, "EDIT_RECIPES_SCHEMAS",
                        copy.deepcopy(edit.EDIT_RECIPES_SCHEMAS))
    from astrocook.gui.main_window import MainWindowV2  # needs the QApplication
    window = MainWindowV2(None, None)
    qtbot.addWidget(window)

    vis = _arm(np.arange(900.0, 1020.0, 0.02), 1.0, seed=1)
    history = SessionHistory(SessionV2(name="VIS", gui=None, spec=vis), HistoryLogV2())
    monkeypatch.setattr(window, "_update_view_for_session", lambda *a, **k: None)
    launched = []
    monkeypatch.setattr(window, "_launch_recipe_dialog",
                        lambda *a, **k: launched.append((a, k)))

    window._on_equalized_stitch_dialog_requested(history, ["NIR"])

    assert launched and launched[0][0] == ("edit", "equalize_and_stitch")
    assert window._hidden_stitch_sessions == "NIR"
    docs = {p["name"]: p["doc"]
            for p in edit.EDIT_RECIPES_SCHEMAS["equalize_and_stitch"]["params"]}
    assert "requires 1 value)" in docs["stitch_wavelengths"]
