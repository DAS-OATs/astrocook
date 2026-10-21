"""Physics of the Voigt fitter not pinned by the snapshots.

Written after mutation testing of ``fitting/`` (``mutmut``, phase 1, item 6
of ``ROADMAP.md``): the turbulent broadening could be changed (e.g.
sqrt(b^2 - btur^2)) without any test failing, and no test fitted a line
broader than 50 km/s. Data from the mock generator (independent of the
fitting code). F32 (findings register): a penalty on b > 50 km/s biased the
fits of broad lines; it is replaced by an explicit upper bound on b
(``b_max``), used only by the automatic pipelines.
"""
import numpy as np
import pytest

from astrocook.core.structures import ComponentDataV2, SystemListDataV2
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.mocks import MockComponent, MockTruth, generate


def _systs(**params):
    base = dict(id=1, z=4.0, dz=None, logN=13.8, dlogN=None, b=20.0, db=None,
                series="Ly_a")
    base.update(params)
    return SystemListV2(SystemListDataV2(components=[ComponentDataV2(**base)]))


def _mock(snr=None, seed=1, **params):
    comp = MockComponent("Ly_a", params.get("z", 4.0), params.get("logN", 13.8),
                         params.get("b", 20.0), params.get("btur", 0.0))
    return generate(MockTruth(components=[comp], lam_min_nm=605.0,
                              lam_max_nm=611.0, dv_kms=2.5, resol=45000.0,
                              snr=snr, z_em=4.1, seed=seed))


def _model(spec, systs):
    return VoigtFitterV2(spec, systs).compute_model_flux()[1]


def test_turbulent_b_adds_in_quadrature():
    spec = _mock()
    with_btur = _model(spec, _systs(b=12.0, btur=16.0))
    thermal_only = _model(spec, _systs(b=20.0, btur=0.0))  # sqrt(12^2 + 16^2)
    np.testing.assert_allclose(with_btur, thermal_only, rtol=1e-12)
    assert not np.allclose(with_btur, _model(spec, _systs(b=12.0)))


def test_fit_with_turbulent_b():
    # btur frozen (default) at the true value: the thermal b is recovered,
    # through the vectorized model used during the fit
    spec = _mock(snr=200.0, b=12.0, btur=16.0)
    fitted, _, _ = VoigtFitterV2(spec, _systs(b=15.0, btur=16.0, logN=13.5)).fit(
        max_nfev=2000, z_window_kms=20.0)
    comp = fitted.components[0]
    assert comp.btur == 16.0
    assert comp.b == pytest.approx(12.0, abs=0.3)
    assert comp.logN == pytest.approx(13.8, abs=0.01)


def _free_btur(systs):
    # b frozen at the true thermal value, btur free: for a single line only
    # sqrt(b^2 + btur^2) is constrained, so both cannot be free
    uuid = systs.components[0].uuid
    return systs.update_constraint(uuid, "b", is_free=False).update_constraint(
        uuid, "btur", is_free=True)


def test_fit_stores_the_btur_error():
    # Regression test for F36 (#559): dbtur was computed and not stored, so
    # it kept its value from before the fit (here 0.7)
    spec = _mock(snr=50.0, b=12.0, btur=16.0)
    systs = _free_btur(_systs(b=12.0, btur=10.0, dbtur=0.7, logN=13.6))
    fitted, _, _ = VoigtFitterV2(spec, systs).fit(max_nfev=2000, z_window_kms=20.0)
    comp = fitted.components[0]
    assert comp.b == 12.0 and comp.db == 0.0
    assert 0.0 < comp.dbtur < 0.7
    assert abs(comp.btur - 16.0) < 3 * comp.dbtur


def test_frozen_btur_has_no_error():
    # btur is frozen unless the user frees it: its error is 0 after a fit,
    # as for the other frozen parameters
    spec = _mock(snr=50.0, b=12.0, btur=16.0)
    systs = _systs(b=15.0, btur=16.0, dbtur=0.7, logN=13.6)
    fitted, _, _ = VoigtFitterV2(spec, systs).fit(max_nfev=2000, z_window_kms=20.0)
    assert fitted.components[0].dbtur == 0.0


def test_spectator_keeps_its_btur_error():
    spec = _mock(snr=50.0)
    comps = [ComponentDataV2(id=i, z=z, dz=0.1, logN=13.6, dlogN=0.1, b=18.0,
                             db=0.1, dbtur=0.7, series="Ly_a")
             for i, z in ((1, 4.0), (2, 4.02))]
    systs = SystemListV2(SystemListDataV2(components=comps))
    with systs.fitting_context([comps[0].uuid], group_depth=0):
        fitted, _, _ = VoigtFitterV2(spec, systs).fit(max_nfev=2000, z_window_kms=20.0)
    fitted_comp, spectator = fitted.components
    assert fitted_comp.dbtur == 0.0
    assert (spectator.dz, spectator.dbtur) == (0.1, 0.7)


def test_evaluate_uncertainties_updates_the_btur_error():
    # The errors recomputed at the end of the pipelines (#534) include dbtur
    import dataclasses
    spec = _mock(snr=50.0, b=12.0, btur=16.0)
    systs = _free_btur(_systs(b=12.0, btur=10.0, logN=13.6))
    fitted, _, _ = VoigtFitterV2(spec, systs).fit(max_nfev=2000, z_window_kms=20.0)
    dbtur = fitted.components[0].dbtur
    stale = SystemListV2(dataclasses.replace(fitted._data, components=[
        dataclasses.replace(fitted.components[0], dbtur=99.0)]))
    evaluated = VoigtFitterV2(spec, stale).evaluate_uncertainties(20.0)
    assert evaluated.components[0].dbtur == pytest.approx(dbtur, rel=0.1)


def test_fit_narrow_line():
    spec = _mock(snr=20.0, b=40.0)
    fitted, _, _ = VoigtFitterV2(spec, _systs(b=32.0, logN=13.5)).fit(
        max_nfev=3000, z_window_kms=20.0)
    assert fitted.components[0].b == pytest.approx(40.0, rel=0.1)


@pytest.mark.parametrize("b", [60.0, 80.0, 120.0])
def test_fit_broad_line(b):
    # Regression test for F32: with the penalty on b > 50 km/s these lines
    # came out with b ~ 50 km/s and logN low by 0.1-0.3 dex
    spec = _mock(snr=100.0, b=b)
    fitted, _, _ = VoigtFitterV2(spec, _systs(b=0.8 * b, logN=13.5)).fit(
        max_nfev=3000, z_window_kms=20.0)
    comp = fitted.components[0]
    assert abs(comp.b - b) < 3 * comp.db
    assert abs(comp.logN - 13.8) < 3 * comp.dlogN


def test_b_max_bounds_the_fit():
    spec = _mock(snr=100.0, b=80.0)
    fitted, _, _ = VoigtFitterV2(spec, _systs(b=64.0, logN=13.5)).fit(
        max_nfev=3000, z_window_kms=20.0, b_max=50.0)
    assert fitted.components[0].b == pytest.approx(50.0, abs=1e-6)
    # A starting value above the bound is brought inside it
    fitted, _, _ = VoigtFitterV2(spec, _systs(b=70.0, logN=13.5)).fit(
        max_nfev=3000, z_window_kms=20.0, b_max=50.0)
    assert fitted.components[0].b <= 50.0


def test_nan_uncertainties_are_reported(monkeypatch, caplog):
    # F20: a negative diagonal element of the covariance (ill-conditioned
    # fit) gave a NaN error with only a numpy RuntimeWarning. The error stays
    # NaN, but a warning names the component and numpy is silenced.
    import logging
    import warnings
    import astrocook.fitting.voigt_fitter as vf

    real_pinv = vf.pinv

    def pinv_negative(a):
        cov = real_pinv(a)
        cov[0, 0] = -abs(cov[0, 0])
        return cov

    monkeypatch.setattr(vf, "pinv", pinv_negative)
    spec = _mock(snr=100.0, b=20.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with caplog.at_level(logging.WARNING):
            fitted, _, _ = VoigtFitterV2(spec, _systs(b=18.0, logN=13.5)).fit(
                max_nfev=200, z_window_kms=20.0)
    comp = fitted.components[0]
    assert np.isnan(comp.dz)
    assert np.isfinite(comp.dlogN) and np.isfinite(comp.db)
    assert "Uncertainties could not be computed" in caplog.text
    assert "Ly_a at z = 4.0" in caplog.text


def test_evaluate_uncertainties_matches_the_fit():
    # Issue #534: errors recomputed at fixed parameters (used at the end of
    # the pipelines) are those of the fit when nothing has changed
    spec = _mock(snr=50.0, b=20.0)
    fitted, _, _ = VoigtFitterV2(spec, _systs(b=18.0, logN=13.5)).fit(
        max_nfev=2000, z_window_kms=20.0)
    evaluated = VoigtFitterV2(spec, fitted).evaluate_uncertainties(20.0)
    a, b = fitted.components[0], evaluated.components[0]
    for key in ("dz", "dlogN", "db", "chi2"):
        assert getattr(b, key) == pytest.approx(getattr(a, key), rel=1e-3), key
    assert (b.z, b.logN, b.b) == (a.z, a.logN, a.b)
