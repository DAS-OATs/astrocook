"""Flags left by a fit on the components (#517, F21): a parameter at a
bound of the fit, or uncertainties that could not be computed."""
import dataclasses

import pytest

from astrocook.core.session import SessionV2
from astrocook.core.structures import ComponentDataV2, SystemListDataV2
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import COMPONENT_FLAGS, VoigtFitterV2
from astrocook.mocks import MockComponent, MockTruth, generate


def _spec(series, z, logN, b, half_width_nm=1.5, lam_rest_nm=121.567):
    lam = lam_rest_nm * (1 + z)
    return generate(MockTruth(
        components=[MockComponent(series, z, logN, b)], lam_min_nm=lam - half_width_nm,
        lam_max_nm=lam + half_width_nm, dv_kms=2.5, resol=45000.0, snr=50.0,
        z_em=z + 0.5, seed=4))


def _systs(series, z, logN, b):
    return SystemListV2(SystemListDataV2(components=[ComponentDataV2(
        id=1, z=z, dz=None, logN=logN, dlogN=None, b=b, db=None, series=series)]))


def test_broad_line_at_the_upper_bound_of_b():
    spec = _spec("Ly_a", 3.0, 13.5, 70.0)
    capped, _, _ = VoigtFitterV2(spec, _systs("Ly_a", 3.0, 13.3, 30.0)).fit(b_max=50.0)
    assert capped.components[0].flags == ("b_at_upper_bound",)
    free, _, _ = VoigtFitterV2(spec, _systs("Ly_a", 3.0, 13.3, 30.0)).fit()
    assert free.components[0].flags == ()


def test_damped_line_at_the_upper_bound_of_logN():
    spec = _spec("Ly_a", 3.0, 20.3, 20.0, half_width_nm=6.0)
    fitted, _, _ = VoigtFitterV2(spec, _systs("Ly_a", 3.0 + 1e-5, 19.8, 30.0)).fit(
        logN_max=19.0, max_nfev=500, z_window_kms=50.0)
    assert "logN_at_upper_bound" in fitted.components[0].flags


def test_flags_are_known():
    assert set(COMPONENT_FLAGS) >= {"b_at_upper_bound", "logN_at_upper_bound", "nan_errors"}


def test_editing_a_parameter_clears_the_flags():
    systs = _systs("CIV", 2.97, 13.0, 10.0)
    comp = dataclasses.replace(systs.components[0], flags=("b_at_upper_bound",))
    systs = SystemListV2(dataclasses.replace(systs._data, components=[comp]))
    assert systs.update_component(comp.uuid, b=20.0).components[0].flags == ()
    assert systs.update_component(comp.uuid, chi2=1.0).components[0].flags == (
        "b_at_upper_bound",)


def test_flags_survive_an_acs2_round_trip(tmp_path):
    systs = _systs("CIV", 2.97, 13.0, 10.0)
    comp = dataclasses.replace(systs.components[0], flags=("b_at_upper_bound", "nan_errors"))
    systs = SystemListV2(dataclasses.replace(systs._data, components=[comp]))
    session = SessionV2("s", None, spec=_spec("CIV", 2.97, 13.0, 10.0, lam_rest_nm=154.8),
                        systs=systs)
    path = tmp_path / "s.acs2"
    session.save(str(path))
    reloaded = SessionV2.open_new(str(path), "s", None, "auto")
    assert reloaded.systs.components[0].flags == ("b_at_upper_bound", "nan_errors")


@pytest.mark.parametrize("flags", [(), ("z_at_window",)])
def test_spectators_keep_their_flags(flags):
    spec = _spec("CIV", 2.97, 13.5, 10.0, lam_rest_nm=154.8)
    systs = SystemListV2(SystemListDataV2(components=[
        ComponentDataV2(id=1, z=2.97, dz=None, logN=13.3, dlogN=None, b=12.0, db=None,
                        series="CIV"),
        ComponentDataV2(id=2, z=2.975, dz=None, logN=12.5, dlogN=None, b=8.0, db=None,
                        series="CIV", flags=flags)]))
    active, spectator = systs.components
    with systs.fitting_context([active.uuid], group_depth=0, spec=spec):
        fitted, _, _ = VoigtFitterV2(spec, systs).fit()
    assert fitted.components[1].flags == flags
