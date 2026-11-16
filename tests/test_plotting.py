"""The plots of a session without the GUI (``astrocook.plotting``): what
``plot_spectrum`` draws (characterized before F158, phase 7-bis PR B11, added
the residuals), the residuals panel, and the plot of transitions in velocity."""
from types import SimpleNamespace

import numpy as np
import pytest

from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.core.atomic_data import ATOM_DATA
from astrocook.plotting import (PLOT_PALETTES, figure_png, line_velocity, plot_spectrum,
                                plot_velocity, residuals_sigma, velocity_transitions)
from astrocook.session.session import Session

C_KMS = 299792.458


@pytest.fixture(scope="module")
def spec():
    truth = MockTruth(components=[MockComponent("CIV", 1.9, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5,
                      resol=45000.0, snr=50.0, z_em=2.0, seed=1)
    return generate(truth)


@pytest.fixture(scope="module")
def session(spec):
    return Session("mock", None, spec=spec).absorbers.add_component(series="CIV", z=1.9)


# --- plot_spectrum, as it was before the residuals --------------------------------

def test_plot_spectrum_is_one_panel_with_its_layers(session, spec):
    fig = plot_spectrum(session, 445, 452)
    assert len(fig.axes) == 1
    ax = fig.axes[0]
    assert tuple(fig.get_size_inches()) == (10, 3.6)
    # step flux, step error, continuum, model: four lines, plus one tick per line
    # of the doublet (axvline)
    assert len(ax.lines) == 4 + 2
    assert [t.get_text() for t in ax.get_legend().get_texts()] == [
        "flux", "error", "continuum", "model"]
    assert ax.get_ylabel().startswith("flux (")
    assert ax.get_xlabel() == f"x ({spec.x.unit})"
    x = spec.x.value
    sel = (x >= 445) & (x <= 452)
    assert ax.get_xlim() == pytest.approx((x[sel].min(), x[sel].max()))
    # The flux drawn is the flux of the data (colour of the palette)
    flux = ax.lines[0]
    assert flux.get_color() == PLOT_PALETTES["standard"]["flux"]
    assert np.allclose(flux.get_ydata(), spec.y.value[sel])


def test_plot_spectrum_normalized_divides_by_the_continuum(session, spec):
    fig = plot_spectrum(session, 445, 452, normalized=True)
    ax = fig.axes[0]
    x = spec.x.value
    sel = (x >= 445) & (x <= 452)
    assert np.allclose(ax.lines[0].get_ydata(), spec.y.value[sel] / spec.cont.value[sel])
    assert ax.get_ylabel() == "normalized flux"


def test_plot_spectrum_without_a_model_has_no_model_line(spec):
    bare = Session("bare", None, spec=spec)
    if bare.spec.has_aux_column("model"):
        pytest.skip("the mock carries a model")
    fig = plot_spectrum(bare, 445, 452)
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert "model" not in labels


def test_plot_spectrum_palette_and_unknown_palette(session):
    colours = PLOT_PALETTES["colorblind"]
    fig = plot_spectrum(session, 445, 452, palette="colorblind")
    assert fig.axes[0].lines[0].get_color() == colours["flux"]
    # An unknown palette is the standard one
    fig = plot_spectrum(session, 445, 452, palette="nope")
    assert fig.axes[0].lines[0].get_color() == PLOT_PALETTES["standard"]["flux"]


# --- the residuals (F158) ----------------------------------------------------------

def test_residuals_panel(session, spec):
    fig = plot_spectrum(session, 445, 452, residuals=True)
    assert len(fig.axes) == 2
    main, res = fig.axes
    assert res.get_ylim() == (-5, 5)
    assert main.get_xlim() == res.get_xlim()  # x shared
    x = spec.x.value
    sel = (x >= 445) & (x <= 452)
    model = session.spec.get_column("model").value[sel]
    expected = (spec.y.value[sel] - model) / spec.dy.value[sel]
    step = [ln for ln in res.lines if len(ln.get_ydata()) == sel.sum()][0]
    assert np.allclose(step.get_ydata(), expected, equal_nan=True)
    # in sigma: the same normalized
    norm = plot_spectrum(session, 445, 452, normalized=True, residuals=True).axes[1]
    step_n = [ln for ln in norm.lines if len(ln.get_ydata()) == sel.sum()][0]
    assert np.allclose(step_n.get_ydata(), expected, equal_nan=True)
    assert figure_png(fig)[:8] == b"\x89PNG\r\n\x1a\n"


def test_residuals_need_a_model(spec):
    bare = Session("bare", None, spec=spec)
    if bare.spec.has_aux_column("model"):
        pytest.skip("the mock carries a model")
    with pytest.raises(ValueError, match="no model"):
        plot_spectrum(bare, 445, 452, residuals=True)
    with pytest.raises(ValueError, match="no model"):
        plot_velocity(bare, 1.9, ["CIV_1548"], residuals=True)


def test_residuals_sigma_with_a_zero_error():
    r = residuals_sigma([1.0, 2.0], [0.0, 0.5], [0.0, 1.0])
    assert np.isnan(r[0]) and r[1] == 2.0


# --- transitions in velocity ---------------------------------------------------

def test_line_velocity_same_transition_and_doublet():
    w = 1548.204
    assert line_velocity(w, 1.9 + 1e-4, w, 1.9) == pytest.approx(C_KMS * 1e-4 / 2.9)
    # the other line of the doublet is where the wavelengths say, at v = 0 of its own
    w2 = 1550.781
    assert line_velocity(w2, 1.9, w, 1.9) == pytest.approx(C_KMS * (w2 - w) / w)


def test_the_inspector_ticks_agree_with_the_shared_formula():
    # the formula of VelocityPlotWidget._update_plot, for two components
    # of a multiplet and of the same line
    z_sys, w_p = 1.9, ATOM_DATA["CIV_1548"]["wave"]
    lam_p = w_p * (1 + z_sys)
    for name, z_c in (("CIV_1548", 1.9004), ("CIV_1550", 1.9)):
        if name == "CIV_1548":
            inspector = C_KMS * (z_c - z_sys) / (1.0 + z_sys)
        else:
            lam_c = ATOM_DATA[name]["wave"] * (1 + z_c)
            inspector = C_KMS * (lam_c - lam_p) / lam_p
        assert line_velocity(ATOM_DATA[name]["wave"], z_c, w_p, z_sys) == pytest.approx(inspector)


def test_velocity_transitions_default_explicit_and_errors(session):
    assert velocity_transitions(session, 1.9)[0] == ["CIV_1548", "CIV_1550"]
    plotted, skipped = velocity_transitions(session, 1.9, ["CIV", "Ly_a", "CIV_1548"])
    assert plotted == ["CIV_1548", "CIV_1550"]  # multiplet expanded, no duplicates
    assert skipped == {"Ly_a": "outside the spectrum"}
    with pytest.raises(ValueError, match="Unknown transitions: XX_1"):
        velocity_transitions(session, 1.9, ["XX_1"])
    with pytest.raises(ValueError, match="No component has a line within 300 km/s of z = 1.5"):
        velocity_transitions(session, 1.5)
    with pytest.raises(ValueError, match="No transition has pixels"):
        velocity_transitions(session, 1.9, ["Ly_a"])
    with pytest.raises(ValueError, match="2 transitions: at most 1 per plot"):
        velocity_transitions(session, 1.9, ["CIV"], max_panels=1)


def test_plot_velocity_panels(session, spec):
    fig = plot_velocity(session, 1.9, ["CIV_1548", "CIV_1550"], half_width_kms=200)
    assert len(fig.axes) == 2
    for ax, name in zip(fig.axes, ("CIV_1548", "CIV_1550"), strict=True):
        assert ax.get_xlim() == (-200, 200)
        assert [t.get_text() for t in ax.texts] == [name]
    # flux in km/s: the minimum of the line at v ~ 0, normalized
    flux = fig.axes[0].lines[0]
    v = flux.get_xdata()
    assert abs(v[np.argmin(flux.get_ydata())]) < 10
    assert flux.get_ydata().max() < 1.3
    # with residuals: two axes per transition
    assert len(plot_velocity(session, 1.9, ["CIV"], residuals=True).axes) == 4
    # not normalized: raw flux
    raw = plot_velocity(session, 1.9, ["CIV_1548"], normalized=False)
    assert "Flux" in raw.texts[1].get_text()
    assert figure_png(fig)[:8] == b"\x89PNG\r\n\x1a\n"


def test_plot_velocity_ticks_the_components(session):
    # the component is at v = 0 in the panel of its line; the other line of
    # the doublet (~500 km/s away) shows only in a window that reaches it
    def ticks(half_width):
        ax = plot_velocity(session, 1.9, ["CIV_1548"], half_width_kms=half_width).axes[0]
        return sorted(float(ln.get_xdata()[0]) for ln in ax.lines
                      if len(ln.get_xdata()) == 2 and ln.get_ydata()[0] == 0.02)
    assert ticks(300) == pytest.approx([0.0], abs=0.1)
    wide = ticks(700)
    assert wide[1] == pytest.approx(line_velocity(
        ATOM_DATA["CIV_1550"]["wave"], 1.9, ATOM_DATA["CIV_1548"]["wave"], 1.9), abs=0.2)


def test_flag_colours_agree_with_the_inspector(session):
    pytest.importorskip("PySide6")
    from types import SimpleNamespace

    from astrocook.gui.inspector_table import flag_cue
    from astrocook.plotting import BOUND_FLAGS, flag_colour
    # flag_cue takes its colours from the palette in use (the standard one)
    colours = PLOT_PALETTES["standard"]
    for flags in [(), *[(f,) for f in (*BOUND_FLAGS, "nan_errors", "not_fitted")],
                  ("nan_errors", "b_at_upper_bound")]:
        comp = SimpleNamespace(flags=flags)
        assert flag_colour(comp, colours) == flag_cue(comp)[0], flags


def test_plot_velocity_flagged_ticks_are_dashed(session):
    from dataclasses import replace
    flagged = replace(session.systs.components[0], flags=("b_at_upper_bound",))
    stub = SimpleNamespace(spec=session.spec, systs=SimpleNamespace(components=[flagged]))
    ax = plot_velocity(stub, 1.9, ["CIV_1548"]).axes[0]
    ticks = [ln for ln in ax.lines if len(ln.get_xdata()) == 2 and ln.get_ydata()[0] == 0.02]
    assert [(t.get_linestyle(), t.get_color()) for t in ticks] == [
        ("--", PLOT_PALETTES["standard"]["flag_bound"])]
    clean = plot_velocity(session, 1.9, ["CIV_1548"]).axes[0]
    assert [ln.get_linestyle() for ln in clean.lines
            if len(ln.get_xdata()) == 2 and ln.get_ydata()[0] == 0.02] == ["-"]
