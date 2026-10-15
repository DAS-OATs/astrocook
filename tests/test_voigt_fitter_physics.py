"""Physics of the Voigt fitter not pinned by the snapshots.

Written after mutation testing of ``fitting/`` (``mutmut``, phase 1, item 6
of ``ROADMAP.md``): the turbulent broadening could be changed (e.g.
sqrt(b^2 - btur^2)) without any test failing, and no test fitted a line
broader than 50 km/s. Data from the mock generator (independent of the
fitting code). Records F32 (findings register): a penalty on b > 50 km/s
biases the fits of broad lines.
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


def test_fit_narrow_line():
    spec = _mock(snr=20.0, b=40.0)
    fitted, _, _ = VoigtFitterV2(spec, _systs(b=32.0, logN=13.5)).fit(
        max_nfev=3000, z_window_kms=20.0)
    assert fitted.components[0].b == pytest.approx(40.0, rel=0.1)


@pytest.mark.xfail(strict=True, reason=(
    "_residual_function appends a penalty 5 (b - 50) for b > 50 km/s to the "
    "residuals: at S/N 20 a line with b = 80 km/s is fitted with b ~ 50 km/s "
    "and logN 0.1 dex low (b = 120: logN 0.35 dex low). Findings register F32 "
    "(ROADMAP.md, section 3.2)."))
def test_fit_broad_line():
    spec = _mock(snr=20.0, b=80.0)
    fitted, _, _ = VoigtFitterV2(spec, _systs(b=64.0, logN=13.5)).fit(
        max_nfev=3000, z_window_kms=20.0)
    assert fitted.components[0].b == pytest.approx(80.0, rel=0.1)
