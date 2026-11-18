"""The evaluation of alternative models in a window of the spectrum (phase 7,
PR 7c.4a): the statistics and the short fit of a hypothesis on the window,
with the fixed components that reach it, against those on the whole spectrum
with the whole list. Data from the mock generator."""
import time

import numpy as np
import pytest

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.structures import SystemListData
from astrocook.core.system_list import SystemList
from astrocook.fitting.evaluation_window import (WINDOW_B, WINDOW_FWHM, EvaluationWindow,
                                                 window_mask)
from astrocook.fitting.lsf import C_KMS
from astrocook.fitting.model_evaluation import evaluate_model
from astrocook.fitting.voigt_fitter import VoigtFitter
from astrocook.mocks import MockComponent, MockTruth, generate

Z = 2.0
LAM_1548 = ATOM_DATA["CIV_1548"]["wave"]
LAM_1550 = ATOM_DATA["CIV_1550"]["wave"]
LAM_LYA = ATOM_DATA["Ly_a"]["wave"]
#: A short fit, as in the experiment of 7c.2 and the pipelines
SHORT_FIT = dict(max_nfev=300, b_max=50.0, logN_max=19.0)
#: The tolerance of the statistics of the window against the whole spectrum
#: (relative, on chi2 and the criteria): the fine grid of the model follows
#: the narrowest line in its pixels, which on the whole spectrum can be
#: another one (measured: about 1e-7)
RTOL = 1e-5


def _forest(resol, n_lines=400, seed=1):
    """A CIV doublet at z = 2 in a forest of Ly-a lines (360-470 nm), some
    of them blended with it, and the forest as the list of the fixed
    components."""
    rng = np.random.default_rng(seed)
    comps = [MockComponent("CIV", z=Z, logN=13.5, b=15.0)]
    for lam in rng.uniform(3600.0, 4700.0, n_lines):
        comps.append(MockComponent("Ly_a", z=lam / LAM_LYA - 1.0,
                                   logN=rng.uniform(12.5, 14.5), b=rng.uniform(15.0, 35.0)))
    truth = MockTruth(components=comps, lam_min_nm=360.0, lam_max_nm=470.0,
                      dv_kms=2.5 if resol > 20000 else 10.0, resol=resol, snr=50.0,
                      z_em=2.9, seed=3)
    fixed = SystemList(SystemListData(components=[]))
    for c in comps[1:]:
        fixed = fixed.add_component(c.series, c.z, logN=c.logN, b=c.b)
    return generate(truth), fixed


def _empty():
    return SystemList(SystemListData(components=[]))


def _hypotheses():
    """A CIV system (true) and two Ly-a on the two lines of the doublet."""
    civ = _empty().add_component("CIV", Z, logN=13.3, b=12.0)
    lya = (_empty()
           .add_component("Ly_a", LAM_1548 * (1 + Z) / LAM_LYA - 1, logN=13.4, b=12.0)
           .add_component("Ly_a", LAM_1550 * (1 + Z) / LAM_LYA - 1, logN=13.1, b=12.0))
    return civ, lya


def _with(fixed, systs):
    return SystemList(SystemListData(components=list(fixed.components) + list(systs.components)))


@pytest.fixture(scope="module", params=[45000.0, 8900.0], ids=["R45000", "R8900"])
def forest(request):
    return _forest(request.param)


# --- The window ---------------------------------------------------------------

def test_the_window_is_3_fwhm_plus_5_b_around_every_line():
    spec, _ = _forest(45000.0, n_lines=0)
    civ, lya = _hypotheses()
    mask = window_mask(spec, [civ, lya])
    x = spec.x.to_value("Angstrom")
    half = WINDOW_FWHM * C_KMS / 45000.0 + WINDOW_B * 12.0
    expected = np.zeros(len(x), dtype=bool)
    for lam in (LAM_1548 * (1 + Z), LAM_1550 * (1 + Z)):
        expected |= np.abs(C_KMS * np.log(x / lam)) <= half
    # The two hypotheses put their lines at the same wavelengths
    assert np.sum(mask ^ expected) <= 2
    assert mask.sum() == pytest.approx(expected.sum(), abs=2)


def test_only_the_fixed_components_that_reach_the_window_are_kept(forest):
    spec, fixed = forest
    win = EvaluationWindow.around(spec, _hypotheses(), fixed)
    assert 0 < len(win.fixed) < len(fixed.components) / 10
    # The others leave the model of the window as it is
    model_all = VoigtFitter(win.spectrum, fixed).compute_model_flux()[1]
    model_kept = VoigtFitter(win.spectrum, SystemList(
        SystemListData(components=win.fixed))).compute_model_flux()[1]
    assert np.max(np.abs(model_all - model_kept) / win.spectrum.norm.value) < 1e-5


def test_the_window_is_a_spectrum_of_its_pixels(forest):
    spec, fixed = forest
    win = EvaluationWindow.around(spec, _hypotheses(), fixed)
    assert len(win.spectrum.x) == win.mask.sum() == len(win.index)
    np.testing.assert_array_equal(win.spectrum.y.value, spec.y.value[win.mask])
    np.testing.assert_array_equal(win.spectrum.norm.value, spec.norm.value[win.mask])
    np.testing.assert_array_equal(win.spectrum.resolution(), spec.resolution()[win.mask])
    assert np.isnan(win.to_parent(np.zeros(len(win.index)))[~win.mask]).all()


def test_no_pixel_is_an_error(forest):
    spec, _ = forest
    with pytest.raises(ValueError, match="no pixel"):
        EvaluationWindow.from_mask(spec, np.zeros(len(spec.x), dtype=bool))


# --- The statistics, against the whole spectrum ---------------------------------

def test_the_evaluation_equals_that_on_the_whole_spectrum(forest):
    spec, fixed = forest
    hyps = _hypotheses()
    win = EvaluationWindow.around(spec, hyps, fixed)
    for systs in hyps:
        uuids = [c.uuid for c in systs.components]
        here = win.evaluate(systs, min_influence=0.0)
        whole = evaluate_model(spec, _with(fixed, systs), win.mask, uuids=uuids,
                               min_influence=0.0)
        assert here.n_pixels == whole.n_pixels
        assert here.k == whole.k == 3 * len(uuids)
        assert sorted(here.components) == sorted(whole.components)
        for name in ("chi2", "aic", "aicc", "bic", "spic", "spic_a", "spic_b"):
            assert getattr(here, name) == pytest.approx(getattr(whole, name), rel=RTOL), name
        np.testing.assert_allclose(here.residuals, whole.residuals, atol=1e-4)


def test_the_fixed_components_are_in_the_model_but_not_counted(forest):
    spec, fixed = forest
    civ, _ = _hypotheses()
    win = EvaluationWindow.around(spec, [civ], fixed)
    alone = EvaluationWindow.around(spec, [civ])
    with_fixed = win.evaluate(civ, min_influence=0.0)
    without = alone.evaluate(civ, min_influence=0.0)
    assert with_fixed.k == without.k == 3
    assert with_fixed.components == without.components == [civ.components[0].uuid]
    # The forest blended with the doublet is part of the model
    assert with_fixed.chi2 != pytest.approx(without.chi2, rel=1e-3)


def test_the_selection_and_the_influence_as_in_evaluate_model(forest):
    spec, fixed = forest
    _, lya = _hypotheses()
    win = EvaluationWindow.around(spec, [lya], fixed)
    first = [lya.components[0].uuid]
    here = win.evaluate(lya, uuids=first)
    whole = evaluate_model(spec, _with(fixed, lya), win.mask, uuids=first)
    assert here.components == whole.components == first
    assert here.spic == pytest.approx(whole.spic, rel=RTOL)


def test_the_short_fit_equals_that_on_the_whole_spectrum(forest):
    spec, fixed = forest
    hyps = _hypotheses()
    win = EvaluationWindow.around(spec, hyps, fixed)
    crit = []
    for systs in hyps:
        uuids = [c.uuid for c in systs.components]
        here = win.fit(systs, **SHORT_FIT)
        assert [c.uuid for c in here.components] == uuids
        full = _with(fixed, systs)
        with full.fitting_context(uuids, group_depth=0):
            whole, _, _ = VoigtFitter(spec, full).fit(**SHORT_FIT)
        by_uuid = {c.uuid: c for c in whole.components}
        for c in here.components:
            w = by_uuid[c.uuid]
            assert c.z == pytest.approx(w.z, abs=1e-8)
            assert c.logN == pytest.approx(w.logN, abs=1e-4)
            assert c.b == pytest.approx(w.b, abs=1e-3)
        # The fixed components did not move
        fixed_after = {c.uuid: c for c in whole.components if c.uuid not in uuids}
        assert all(fixed_after[c.uuid].logN == c.logN for c in fixed.components)
        here_ev = win.evaluate(here, min_influence=0.0)
        whole_ev = evaluate_model(spec, whole, win.mask, uuids=uuids, min_influence=0.0)
        assert here_ev.spic == pytest.approx(whole_ev.spic, rel=RTOL)
        crit.append(here_ev.spic)
    # The true hypothesis wins
    assert crit[0] < crit[1]


def test_a_fit_keeps_the_constraints_of_the_hypothesis(forest):
    spec, fixed = forest
    _, lya = _hypotheses()
    a, b = (c.uuid for c in lya.components)
    linked = lya.update_constraint(b, "b", is_free=False, expression=f"p['{a}'].b",
                                   target_uuid=a)
    win = EvaluationWindow.around(spec, [linked], fixed)
    fitted = win.fit(linked, **SHORT_FIT)
    fa, fb = fitted.components
    assert fb.b == pytest.approx(fa.b)
    assert win.evaluate(fitted, min_influence=0.0).k == 5


def test_a_hypothesis_without_components(forest):
    spec, fixed = forest
    civ, _ = _hypotheses()
    win = EvaluationWindow.around(spec, [civ], fixed)
    none = _empty()
    assert win.fit(none) is none
    ev = win.evaluate(none)
    assert ev.k == 0 and ev.n_pixels == len(win.index)
    whole = evaluate_model(spec, fixed, win.mask, uuids=[])
    assert ev.chi2 == pytest.approx(whole.chi2, rel=RTOL)


# --- The time -------------------------------------------------------------------

def test_the_window_is_faster_than_the_whole_spectrum():
    spec, fixed = _forest(45000.0)
    civ, _ = _hypotheses()
    uuids = [civ.components[0].uuid]
    t0 = time.perf_counter()
    win = EvaluationWindow.around(spec, [civ], fixed)
    win.evaluate(win.fit(civ, **SHORT_FIT), min_influence=0.0)
    t_window = time.perf_counter() - t0
    full = _with(fixed, civ)
    t0 = time.perf_counter()
    with full.fitting_context(uuids, group_depth=0):
        fitted, _, _ = VoigtFitter(spec, full).fit(**SHORT_FIT)
    evaluate_model(spec, fitted, win.mask, uuids=uuids, min_influence=0.0)
    t_whole = time.perf_counter() - t0
    # Measured: about 15 ms against 250 ms on 32 thousand pixels
    assert t_window < t_whole / 3
