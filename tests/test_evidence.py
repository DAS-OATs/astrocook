"""The Laplace evidence of a hypothesis (phase 7, PR 7c.4d-quater):
:mod:`astrocook.fitting.evidence` and the option ``evidence`` of the
identification by hypotheses."""
import json
import math

import numpy as np
import pytest

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.identification import META_KEY
from astrocook.core.structures import SystemListData
from astrocook.core.system_list import SystemList
from astrocook.fitting import evidence as E
from astrocook.fitting import hypotheses as H
from astrocook.fitting.evaluation_window import EvaluationWindow
from astrocook.fitting.voigt_fitter import VoigtFitter
from astrocook.mocks import MockComponent, MockTruth, generate

Z = 500.0 / (ATOM_DATA["CIV_1548"]["wave"] / 10.0) - 1.0


def _civ(snr=30.0, logn=13.5, b=8.0):
    return generate(MockTruth(components=[MockComponent("CIV", Z, logn, b)], lam_min_nm=498.0,
                              lam_max_nm=503.5, dv_kms=2.5, resol=45000.0, snr=snr,
                              z_em=3.5, seed=3))


def _fit(win, series_z):
    systs = SystemList(SystemListData(components=[]))
    for series, z in series_z:
        systs = systs.add_component(series, z, logN=13.0, b=10.0)
    return win.fit(systs, max_nfev=300, b_max=100.0, logN_max=19.0)


# --- The pieces -------------------------------------------------------------

def test_laplace_is_exact_for_a_linear_gaussian_model():
    """A linear model with a Gaussian prior: the evidence is the Gaussian
    integral, which the formula gives exactly (the step from the maximum of
    the likelihood to that of the posterior included)."""
    rng = np.random.default_rng(1)
    a = rng.normal(size=(40, 2))
    y = rng.normal(size=40)
    theta_ml = np.linalg.lstsq(a, y, rcond=None)[0]
    resid = y - a @ theta_ml
    mu, s = np.array([0.3, -0.2]), np.array([0.5, 2.0])
    lp = float(np.sum(-0.5 * ((theta_ml - mu) / s) ** 2 - np.log(s * np.sqrt(2 * np.pi))))
    grad = -(theta_ml - mu) / s ** 2
    ev = E.laplace(float(resid @ resid), -a, resid, lp, grad, 1.0 / s ** 2)
    # The exact integral of exp(-|y - A t|^2 / 2) N(t; mu, s)
    cov = a @ np.diag(s ** 2) @ a.T + np.eye(40)
    r = y - a @ mu
    # (the likelihood without its normalization, as in the module)
    exact = -0.5 * r @ np.linalg.solve(cov, r) - 0.5 * np.linalg.slogdet(cov)[1]
    assert ev.ln_z == pytest.approx(exact, abs=1e-9)
    assert ev.criterion == pytest.approx(-2 * ev.ln_z)


@pytest.mark.parametrize("prior", [E.B_FOREST, E.B_METAL, E.B_ANY])
@pytest.mark.parametrize("b", [3.0, 9.0, 30.0, 70.0])
def test_the_priors_of_b_their_derivatives(prior, b):
    lp, d1, c = prior.terms(b)
    h = 1e-4
    num = (prior.terms(b + h)[0] - prior.terms(b - h)[0]) / (2 * h)
    num2 = (prior.terms(b + h)[0] - 2 * lp + prior.terms(b - h)[0]) / h ** 2
    assert d1 == pytest.approx(num, rel=1e-5, abs=1e-8)
    assert c == pytest.approx(max(0.0, -num2), rel=1e-3, abs=1e-6)


@pytest.mark.parametrize("prior", [E.B_FOREST, E.B_METAL, E.B_ANY])
def test_the_priors_of_b_are_normalized(prior):
    b = np.linspace(1.0, 100.0, 200001)
    dens = np.exp([prior.terms(x)[0] for x in b])
    total = np.trapezoid(dens, b)
    # The Gaussian of the forest and the lognormal have tails below 1 km/s
    assert total == pytest.approx(1.0, abs=0.01)


def test_the_column_prior_is_the_differential_incidence():
    """rate * phi(log N) does not depend on the reference column."""
    pop = E.POPULATIONS["HI"]
    lp, slope, curv = pop.logn_terms(13.0)
    shifted = E.Population("HI", pop.rate * 10 ** ((pop.beta - 1) * 1.0), pop.z_ref, pop.gamma,
                           pop.logn_ref - 1.0, pop.beta, pop.b, pop.source)
    assert pop.log_rate(3.0) + lp == pytest.approx(shifted.log_rate(3.0)
                                                    + shifted.logn_terms(13.0)[0])
    assert slope == pytest.approx(-(pop.beta - 1) * math.log(10)) and curv == 0.0


def test_the_species_of_a_series():
    assert [E.species_of(s) for s in ("Ly_a", "Ly_ab", "CIV", "OI_1039,OI_1302", "CIII_977",
                                      "SiII")] == ["HI", "HI", "CIV", "OI", "CIII", "SiII"]


def test_a_known_system_is_uniform_over_its_range():
    p = E.at_system("CII", (2.0, 2.001))
    assert p.log_rate == pytest.approx(math.log(E.P_AT_SYSTEM / 0.001))
    assert E.at_system("HI", (2.0, 2.001)).log_rate == pytest.approx(math.log(1000.0))


# --- On a spectrum ----------------------------------------------------------

def test_the_fit_context_on_given_pixels():
    """``prepare_fit_context(pixels=...)`` takes the residuals on those
    pixels; without them, those of the data (unchanged)."""
    spec = _civ()
    systs = SystemList(SystemListData(components=[])).add_component("CIV", Z, logN=13.5, b=8.0)
    f1 = VoigtFitter(spec, systs)
    f1.prepare_fit_context()
    default = f1._fit_mask.copy()
    f2 = VoigtFitter(spec, systs)
    f2.prepare_fit_context(pixels=np.ones(len(spec.x), dtype=bool))
    assert f2._fit_mask.sum() > default.sum()
    np.testing.assert_array_equal(f2._fit_mask, f2._valid_mask)


def test_the_evidence_on_a_window():
    """The chi2 of the evidence is that of ``evaluate``; a narrow CIV doublet
    is more probable as a CIV than as two Ly-a lines with the same chi2."""
    spec = _civ(logn=13.2, b=7.0)
    lam = [ATOM_DATA[t]["wave"] * (1 + Z) for t in ("CIV_1548", "CIV_1550")]
    x = spec.x.to_value("nm") * 10
    mask = np.zeros(len(x), dtype=bool)
    for la in lam:
        mask |= np.abs(np.log(x / la)) * 299792.458 < 200
    win = EvaluationWindow.from_mask(spec, mask)
    civ = _fit(win, [("CIV", Z)])
    lya = _fit(win, [("Ly_a", la / ATOM_DATA["Ly_a"]["wave"] - 1) for la in lam])
    ev = {}
    for name, systs in (("civ", civ), ("lya", lya)):
        priors = {c.uuid: E.ComponentPrior(0.0, E.free_line(E.species_of(c.series), c.z)
                                           .population) for c in systs.components}
        ev[name] = E.window_evidence(win, systs, priors)
        assert ev[name].chi2 == pytest.approx(win.evaluate(systs, min_influence=0).chi2,
                                              rel=1e-9)
        assert ev[name].k == 3 * len(systs.components)
    assert ev["lya"].chi2 <= ev["civ"].chi2 + 1.0
    # The b of the forest makes two narrow Ly-a improbable
    assert ev["civ"].criterion < ev["lya"].criterion - 10


def test_a_component_without_prior_is_refused():
    spec = _civ()
    win = EvaluationWindow.from_mask(spec, np.ones(len(spec.x), dtype=bool))
    civ = _fit(win, [("CIV", Z)])
    with pytest.raises(KeyError):
        E.window_evidence(win, civ, {})


def test_no_component_is_the_chi2():
    spec = _civ()
    win = EvaluationWindow.from_mask(spec, np.ones(len(spec.x), dtype=bool))
    empty = SystemList(SystemListData(components=[]))
    ev = E.window_evidence(win, empty, {})
    assert ev.ln_z == pytest.approx(-0.5 * win.evaluate(empty).chi2)


# --- The option of the identification --------------------------------------

def _red_civ():
    """A CIV doublet and an unrelated line on the red side."""
    comps = [MockComponent("CIV", Z, 13.4, 8.0), MockComponent("CIV", Z + 0.004, 13.0, 6.0)]
    return generate(MockTruth(components=comps, lam_min_nm=490.0, lam_max_nm=510.0,
                              dv_kms=2.5, resol=45000.0, snr=50.0, z_em=2.4, seed=2))


def test_the_evidence_is_recorded_without_changing_the_decisions():
    spec = _red_civ()
    base = H.identify_hypotheses(spec, ["CIV"], sweep="none").spectrum
    rec = H.identify_hypotheses(spec, ["CIV"], sweep="none", evidence="record").spectrum
    for key in ("series_map_json", "z_map_json"):
        assert json.loads(base.meta[key]) == json.loads(rec.meta[key])
    structure = rec.meta[META_KEY]
    hyps = [h for g in structure.groups for h in g.hypotheses if h.criteria.get("spic")]
    assert hyps and all("evidence" in h.criteria for h in hyps)
    assert structure.params["evidence"] == "record"


def test_the_evidence_can_decide():
    spec = _red_civ()
    out = H.identify_hypotheses(spec, ["CIV"], sweep="none", evidence="decide").spectrum
    structure = out.meta[META_KEY]
    assert structure.params["criterion"] == "evidence"
    accepted = [h for g in structure.groups for h in g.hypotheses
                if str(h.status) == "accepted" and "evidence" in h.reason]
    assert accepted


def test_an_unknown_evidence_is_refused():
    with pytest.raises(ValueError, match="evidence"):
        H.identify_hypotheses(_civ(), ["CIV"], evidence="always")


def test_the_evidence_is_off_by_default():
    assert H.EVIDENCE is None
