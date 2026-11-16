"""The evaluation of a model on fixed pixels (F154, phase 7, PR 7a.5):
chi2, k, AIC, AICc, BIC and the SpIC of Webb et al. (2021), and the recipe
``absorbers.evaluate_model``. Data from the mock generator."""
import logging

import numpy as np
import pytest

from astrocook.core.structures import ComponentData, SystemListData
from astrocook.core.system_list import SystemList
from astrocook.fitting.model_evaluation import (VelocityComponent, evaluate_model,
                                                information_criteria, spic)
from astrocook.fitting.voigt_fitter import VoigtFitter
from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.recipes.errors import InvalidParameter
from astrocook.recipes.result import run
from astrocook.session.session import Session

Z = 2.0


@pytest.fixture(scope="module")
def spec():
    truth = MockTruth(
        components=[MockComponent("CIV", z=Z, logN=13.5, b=15.0),
                    MockComponent("SiIV", z=Z, logN=13.0, b=12.0)],
        lam_min_nm=415.0, lam_max_nm=470.0, dv_kms=2.5, resol=45000.0,
        snr=50.0, z_em=2.2, seed=3)
    return generate(truth)


def _systs(link=False):
    comps = [ComponentData(id=1, z=Z, dz=None, logN=13.5, dlogN=None, b=15.0, db=None,
                           series="CIV"),
             ComponentData(id=2, z=Z, dz=None, logN=13.0, dlogN=None, b=12.0, db=None,
                           series="SiIV")]
    systs = SystemList(SystemListData(components=comps))
    if link:
        civ, siiv = (c.uuid for c in systs.components)
        systs = systs.update_constraint(siiv, "z", is_free=False,
                                        expression=f"p['{civ}'].z", target_uuid=civ)
    return systs


def _pixels(spec, *ranges):
    x = spec.x.to_value("nm")
    return np.any([(x >= lo) & (x <= hi) for lo, hi in ranges], axis=0)


#: The SiIV doublet and the CIV doublet at z = 2
RANGES = ((417.5, 420.0), (463.8, 465.2))


# --- The criteria -------------------------------------------------------------

def test_information_criteria():
    aic, aicc, bic = information_criteria(100.0, 50, 4)
    assert aic == 108.0
    assert aicc == pytest.approx(108.0 + 2 * 4 * 5 / 45)
    assert bic == pytest.approx(100.0 + 4 * np.log(50))
    # Too few pixels for the AICc
    assert information_criteria(1.0, 5, 4)[1] is None


def test_spic_is_equation_3_of_webb_et_al_2021():
    terms = [VelocityComponent(["a"], 3, 40.0), VelocityComponent(["b", "c"], 6, 500.0)]

    def by_hand(f):
        return 100.0 + sum(2 * f * k * r / (r - k - 1) + (1 - f) * k * np.log(r)
                           for k, r in ((3, 40.0), (6, 500.0)))
    for f in (0.0, 0.5, 1.0):
        assert spic(100.0, terms, f) == pytest.approx(by_hand(f), rel=1e-14)
    # f = 1 tends to AICc, f = 0 to BIC, with R_a in place of N
    assert spic(0.0, [VelocityComponent(["a"], 3, 1e9)], 1.0) == pytest.approx(6.0, rel=1e-6)
    assert spic(0.0, [VelocityComponent(["a"], 3, 40.0)], 0.0) == pytest.approx(3 * np.log(40))


def test_spic_bounds_r_from_below_and_skips_components_with_no_free_parameter():
    # R_a >= k_a + 2 (section 2.4): the first term would diverge at k_a + 1
    weak = spic(0.0, [VelocityComponent(["a"], 3, 0.1)], 0.5)
    assert weak == pytest.approx(spic(0.0, [VelocityComponent(["a"], 3, 5.0)], 0.5))
    assert spic(7.0, [VelocityComponent(["a"], 0, 0.0)], 0.5) == 7.0


# --- The evaluation -----------------------------------------------------------

def test_chi2_and_residuals_on_the_pixels(spec):
    systs = _systs()
    pixels = _pixels(spec, *RANGES)
    ev = evaluate_model(spec, systs, pixels)
    _, model = VoigtFitter(spec, systs).compute_model_flux()
    y, dy = spec.y.value[pixels], spec.dy.value[pixels]
    resid = (y - model[pixels]) / dy
    assert ev.n_pixels == int(pixels.sum())
    assert ev.chi2 == pytest.approx(np.sum(resid ** 2), rel=1e-10)
    assert ev.k == 6 and ev.dof == ev.n_pixels - 6
    assert ev.chi2_reduced == pytest.approx(ev.chi2 / ev.dof)
    assert ev.resid_mean == pytest.approx(np.mean(resid), rel=1e-8, abs=1e-12)
    assert ev.frac_beyond_3sigma == pytest.approx(np.mean(np.abs(resid) > 3))
    assert np.allclose(ev.residuals[pixels], resid, rtol=1e-10)
    assert np.all(np.isnan(ev.residuals[~pixels]))
    # The true model of the mock fits: reduced chi2 about 1
    assert 0.85 < ev.chi2_reduced < 1.15


def test_r_of_a_component_alone_is_its_depth_over_the_error(spec):
    # With one component the model is its profile: R_a = sum (1 - M/C) / (dy/C)
    full = _systs()
    systs = full.delete_component(uuid=full.components[1].uuid)
    pixels = _pixels(spec, *RANGES)
    ev = evaluate_model(spec, systs, pixels)
    _, model = VoigtFitter(spec, systs).compute_model_flux()
    cont = spec.norm.value[pixels]
    expected = np.sum((1 - model[pixels] / cont) / (spec.dy.value[pixels] / cont))
    (vc,) = ev.velocity_components
    assert vc.k == 3 and vc.r == pytest.approx(expected, rel=1e-6)


def test_linked_and_frozen_parameters_are_not_counted(spec):
    systs = _systs(link=True)
    civ, siiv = (c.uuid for c in systs.components)
    systs = systs.update_constraint(siiv, "b", is_free=False)
    ev = evaluate_model(spec, systs, _pixels(spec, *RANGES))
    # CIV z, logN, b; SiIV logN only (z linked, b frozen; btur frozen by default)
    assert ev.k == 4
    # One velocity component: the redshifts are linked
    (vc,) = ev.velocity_components
    assert vc.uuids == [civ, siiv] and vc.k == 4
    alone = evaluate_model(spec, _systs(), _pixels(spec, *RANGES))
    assert vc.r == pytest.approx(sum(v.r for v in alone.velocity_components), rel=1e-12)
    assert ev.spic != pytest.approx(alone.spic)


def test_the_selection_counts_components_but_the_model_is_whole(spec):
    systs = _systs(link=True)
    civ, siiv = (c.uuid for c in systs.components)
    pixels = _pixels(spec, *RANGES)
    whole = evaluate_model(spec, systs, pixels)
    only_siiv = evaluate_model(spec, systs, pixels, [siiv])
    assert only_siiv.chi2 == whole.chi2
    assert only_siiv.k == 2
    (vc,) = only_siiv.velocity_components
    assert vc.uuids == [siiv]
    none = evaluate_model(spec, systs, pixels, [])
    assert none.k == 0 and none.velocity_components == []
    assert none.aic == none.bic == none.spic == none.chi2


def test_a_better_model_has_lower_criteria(spec):
    pixels = _pixels(spec, *RANGES)
    full = _systs()
    with_siiv = evaluate_model(spec, full, pixels)
    without = evaluate_model(spec, full.delete_component(uuid=full.components[1].uuid),
                             pixels)
    for name in ("aic", "aicc", "bic", "spic", "spic_a", "spic_b"):
        assert getattr(with_siiv, name) < getattr(without, name), name


def test_no_valid_pixel_is_an_error(spec):
    with pytest.raises(ValueError, match="no valid pixel"):
        evaluate_model(spec, _systs(), np.zeros(len(spec.x), dtype=bool))


# --- The recipe ---------------------------------------------------------------

@pytest.fixture(scope="module")
def session(spec):
    return Session("eval", None, spec=spec, systs=_systs(link=True))


def test_recipe_reports_the_evaluation(session, caplog):
    caplog.set_level(logging.INFO, logger="astrocook")
    result = run(session, "absorbers.evaluate_model",
                 x_ranges="417.5-420.0, 463.8-465.2")
    ev = result.extra["evaluation"]
    expected = evaluate_model(session.spec, session.systs, _pixels(session.spec, *RANGES))
    assert ev["chi2"] == expected.chi2 and ev["k"] == 5 and ev["spic"] == expected.spic
    assert ev["velocity_components"][0]["k"] == 5
    assert result.session is session
    assert any(r.getMessage().startswith(f"Model evaluated on {expected.n_pixels} pixels, "
                                         "2 components counted: chi2 = ")
               for r in caplog.records)


def test_recipe_selection_and_residuals_column(session):
    result = run(session, "absorbers.evaluate_model", x_ranges="463.8-465.2",
                 series="SiIV", z_min=1.9, z_max=2.1, residuals_col="resid")
    assert result.extra["evaluation"]["k"] == 2
    assert "resid" in result.changes.columns_added
    col = result.session.spec.get_column("resid").value
    pixels = _pixels(session.spec, (463.8, 465.2))
    assert np.all(np.isfinite(col[pixels])) and np.all(np.isnan(col[~pixels]))


def test_recipe_mask_column(session):
    pixels = _pixels(session.spec, *RANGES)
    masked = session.with_new_spectrum(session.spec.update_column("sel", pixels.astype(float)))
    by_mask = run(masked, "absorbers.evaluate_model", mask_col="sel").extra["evaluation"]
    by_ranges = run(session, "absorbers.evaluate_model",
                    x_ranges="417.5-420.0, 463.8-465.2").extra["evaluation"]
    assert by_mask["chi2"] == by_ranges["chi2"] and by_mask["n_pixels"] == by_ranges["n_pixels"]


def test_recipe_warns_when_the_selection_matches_nothing(session):
    result = run(session, "absorbers.evaluate_model", x_ranges="463.8-465.2", series="MgII")
    assert result.extra["evaluation"]["k"] == 0
    assert "No component matches the selection: none is counted (k = 0)." in result.warnings


@pytest.mark.parametrize("params, message", [
    ({"x_ranges": "417-416"}, "min < max"),
    ({"x_ranges": "abc"}, "is not a range"),
    ({"x_ranges": "300-310"}, "no valid pixel"),
    ({"mask_col": "nope"}, "no column 'nope'"),
    ({"residuals_col": "y"}, "column of the spectrum itself"),
    ({"uuids": ["nope"]}, "no component 'nope'"),
    ({"z_min": 2.1, "z_max": 2.0}, "is above"),
])
def test_recipe_refuses_bad_inputs(session, params, message):
    with pytest.raises(InvalidParameter, match=message):
        session.absorbers.evaluate_model(**params)
