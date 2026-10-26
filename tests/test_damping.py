"""Damping wings of the Voigt profile (F61): the damping parameter is
a = Gamma lambda_0 / (4 pi b), with the rest wavelength. Before, the fitter
used the observed one, (1 + z) times larger, and the column of a damped
Ly-alpha came out low by log(1 + z). The benchmark on mocks cannot measure
this (no damped systems in its mocks), so this test does."""
import pytest

from astrocook.core.structures import ComponentDataV2, SystemListDataV2
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.mocks import MockComponent, MockTruth, generate


@pytest.mark.parametrize("z", [2.5, 3.5])
def test_fit_recovers_the_column_of_a_damped_lya(z):
    # Before F61: logN 19.75 at z = 2.5 and 19.64 at z = 3.5
    lam = 121.567 * (1 + z)
    spec = generate(MockTruth(
        components=[MockComponent("Ly_a", z, 20.3, 20.0)], lam_min_nm=lam - 6.0,
        lam_max_nm=lam + 6.0, dv_kms=5.0, resol=45000.0, snr=30.0, z_em=z + 0.5,
        seed=3))
    systs = SystemListV2(SystemListDataV2(components=[ComponentDataV2(
        id=1, z=z + 1e-5, dz=None, logN=19.8, dlogN=None, b=30.0, db=None,
        series="Ly_a")]))
    fitted, _, _ = VoigtFitterV2(spec, systs).fit(max_nfev=500, z_window_kms=50.0)
    assert fitted.components[0].logN == pytest.approx(20.3, abs=0.02)
