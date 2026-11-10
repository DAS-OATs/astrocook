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

from astrocook.core.session import SessionV2
from astrocook.mocks import MockComponent, MockTruth, generate

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


def test_characterization_found_then_lost(civ_and_siiv, caplog):
    caplog.set_level(logging.INFO, logger="astrocook")
    out = civ_and_siiv.absorbers.metals_auto()
    messages = [r.getMessage() for r in caplog.records]
    assert "Found SiIV at z = 1.9 (logN = 13.30)." in messages
    # The ions out of the spectrum, fitted on the whole spectrum, kept at
    # their first values
    assert "Found MgII at z = 1.9 (logN = 13.00)." in messages
    assert sum("the whole spectrum is used" in m for m in messages) == 4
    assert [c.series for c in out.systs.components] == ["CIV"]
