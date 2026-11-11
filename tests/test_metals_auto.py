"""metals_auto and the fit of components with no pixel (F111, phase 6-ter,
PR 7-quater).

On the mock of trial 1 of PR 7 ``metals_auto`` kept no ion and missed the
SiIV at z = 1.90 (logN 13.3, b 12). The trial fit found it; the recipe then
refitted the session it started from, so what it found was lost. The ions
whose lines fall outside the spectrum had no pixel to fit, and the fitter
fitted the whole spectrum instead (hours on J0103-1305).
"""
import logging

import pytest

from astrocook.session.session import SessionV2
from astrocook.fitting.voigt_fitter import NoPixelsToFit, VoigtFitterV2
from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.recipes.errors import RecipeCancelled, RecipeFailed
from astrocook.recipes.progress import cancel_requested, cancellable, parse_progress
from astrocook.recipes.result import run

Z = 1.90


@pytest.fixture(scope="module")
def civ_and_siiv():
    """CIV and SiIV at z = 1.90 (the values of the full-chain mock), 400-455
    nm: the SiIV doublet at 404-407 nm, the CIV doublet at 449-450 nm. The
    session has the CIV, fitted, as doublets_auto leaves it."""
    truth = MockTruth(components=[MockComponent("CIV", Z, 14.2, 15.0),
                                  MockComponent("SiIV", Z, 13.3, 12.0)],
                      lam_min_nm=400.0, lam_max_nm=455.0, dv_kms=2.5,
                      resol=45000.0, snr=30.0, z_em=2.0, seed=7)
    session = SessionV2("civ_siiv", None, spec=generate(truth, truth_columns=False))
    session = session.continuum.estimate_auto()
    return session.absorbers.add_component("CIV", str(Z), "14.0", "12.0")


def _found(session, series):
    return [c for c in session.systs.components if c.series == series]


def test_the_ions_found_are_kept(civ_and_siiv, caplog):
    # It was: "Found SiIV at z = 1.9 (logN = 13.30)" in the log, then the
    # session it started from refitted, with the CIV only; the ions out of
    # the spectrum fitted on the whole spectrum and "found" at their first
    # values, logN 13.00 (F111)
    caplog.set_level(logging.INFO, logger="astrocook")
    result = run(civ_and_siiv, "absorbers.metals_auto")
    (siiv,) = _found(result.session, "SiIV")
    assert siiv.z == pytest.approx(Z, abs=2e-5)
    assert siiv.logN == pytest.approx(13.3, abs=0.1) and siiv.b == pytest.approx(12.0, abs=1.5)
    assert [c.series for c in result.session.systs.components] == ["CIV", "SiIV"]
    messages = [r.getMessage() for r in caplog.records]
    assert not any("whole spectrum" in m for m in messages)
    assert ("4 trial ions not tried, their lines outside the spectrum: NV at z = 1.9, "
            "OVI at z = 1.9, CII at z = 1.9, MgII at z = 1.9.") in result.warnings
    # The ions tried, with their AIC, and the progress of the trials
    (tried,) = result.extra["tried"]
    assert tried["ion"] == "SiIV" and tried["kept"] and tried["delta_aic"] < -100
    assert tried["reason"] == ""
    assert any(m.startswith("Tried SiIV at z = 1.9: logN = 13.") and m.endswith(": kept.")
               for m in messages)
    trials = [parse_progress(m) for m in messages[:messages.index("Added 1 metal components.")]
              if parse_progress(m) is not None]
    assert trials == [0, 20, 40, 60, 80, 100]


@pytest.mark.parametrize("params, reason", [
    # b of the SiIV is 12: a bound of 10 holds it there
    ({"max_b": "10"}, "b at its bound (10 km/s)"),
    # Its AIC change is some -7000 on its 42 pixels
    ({"aic_penalty": "1e6"}, "AIC change above -1e+06"),
])
def test_an_ion_at_a_bound_or_with_too_little_aic_is_not_kept(civ_and_siiv, params, reason):
    # F121 (Guido's choice, phase 6-ter PR 7-quater): the values at the
    # bounds were kept (67 of 192 ions on J0103-1305, b = 50), and
    # aic_penalty was not used
    result = run(civ_and_siiv, "absorbers.metals_auto", **params)
    assert not _found(result.session, "SiIV")
    (tried,) = result.extra["tried"]
    assert not tried["kept"] and tried["reason"] == reason


def test_a_component_with_no_pixel_is_not_fitted(civ_and_siiv):
    # NV at z = 1.9 is at 359 nm, outside the spectrum (400-455 nm): the
    # fitter fitted the whole spectrum, with a WARNING (F111)
    systs = civ_and_siiv.systs.add_component("NV", Z, logN=13.0, b=10.0)
    new = systs.components[-1].uuid
    with systs.fitting_context([new], spec=civ_and_siiv.spec):
        fitter = VoigtFitterV2(civ_and_siiv.spec, systs)
        with pytest.raises(NoPixelsToFit, match=r"\(NV at z = 1\.9000\): their lines "
                                                 r"fall outside the spectrum \(400\.0-455\.0 nm\)"):
            fitter.fit(max_nfev=10)
    # A recipe says so
    with pytest.raises(RecipeFailed, match="No valid pixels around the components"):
        civ_and_siiv.absorbers.add_component("NV", str(Z), "13.0", "10.0")


def test_a_cancelled_recipe_stops_at_its_next_trial(civ_and_siiv):
    # Cancelled during the first trial: it stops before the next one, with
    # an error, and makes no new state (phase 6-ter, PR 7-quater)
    lines = []

    class Cancel(logging.Handler):
        def emit(self, record):
            lines.append(record.getMessage())
            if record.getMessage().startswith("Tried SiIV"):
                event.set()

    handler = Cancel(level=logging.INFO)
    logger = logging.getLogger("astrocook")
    level = logger.level
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    try:
        with cancellable() as event:
            with pytest.raises(RecipeCancelled, match="cancelled while it ran"):
                run(civ_and_siiv, "absorbers.metals_auto")
    finally:
        logger.removeHandler(handler)
        logger.setLevel(level)
    assert lines[-1] == "Found SiIV at z = 1.9 (logN = 13.30)."
    assert not cancel_requested()  # the event belongs to the block


def test_a_cancelled_fit_stops(civ_and_siiv):
    # The fitter checks at each evaluation of the residuals
    with cancellable() as event:
        event.set()
        with pytest.raises(RecipeCancelled):
            civ_and_siiv.absorbers.refit_all()
