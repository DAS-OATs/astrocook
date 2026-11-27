"""The flags of the values that are not measured (phase 7, PR 7d.2; F52,
F116): ``unresolved``, ``saturated``, ``large_errors``,
``telluric_overlap``, set by the fitter with those of #517."""
import dataclasses

import numpy as np
import pytest

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.structures import ComponentData, DataColumn, SystemListData
from astrocook.core.system_list import SystemList
from astrocook.fitting.component_flags import (
    MEASUREMENT_FLAGS, FitPixels, b_instrumental, central_tau, flag_quantities,
    large_error_params, measurement_flags, with_measurement_flags)
from astrocook.fitting.voigt_fitter import COMPONENT_FLAGS, VoigtFitter
from astrocook.methods_text import FLAG_NAMES
from astrocook.mocks import MockComponent, MockTruth, generate

C_KMS = 299792.458
CIV = 154.8204


def _comp(series="CIV", z=2.97, logN=13.0, b=10.0, dz=1e-6, dlogN=0.05, db=0.5, **kw):
    return ComponentData(id=1, z=z, dz=dz, logN=logN, dlogN=dlogN, b=b, db=db,
                         series=series, **kw)


def _pixels(z=2.97, resol=45000.0, telluric=None, rest_nm=CIV, fit_mask=None):
    """Pixels of 2.5 km/s around a line, with a flat R."""
    centre = rest_nm * 10 * (1 + z)
    wave = centre * np.exp(np.arange(-2000, 2000) * 2.5 / C_KMS)
    return FitPixels(wave=wave, resol=np.full(wave.size, resol), fit_mask=fit_mask,
                     telluric=telluric(wave) if telluric else None)


# --- The pure function --------------------------------------------------------

def test_b_instrumental():
    assert b_instrumental(45000.0) == pytest.approx(4.0, abs=0.01)
    assert b_instrumental(float("nan")) == 0.0 and b_instrumental(0.0) == 0.0


def test_central_tau_is_that_of_the_fitter():
    """tau_0 of the intrinsic profile, as the fitter computes it."""
    from astrocook.fitting.voigt_fitter import voigt_optical_depth
    at = ATOM_DATA["CIV_1550"]
    centre = at["wave"] * 3.97
    tau = voigt_optical_depth(np.array([centre]), np.array([at["wave"]]),
                              np.array([at["f"]]), np.array([at["gamma"]]),
                              np.array([2.97]), np.array([10 ** 14.0]), np.array([8.0]))
    assert central_tau(at, 14.0, 8.0) == pytest.approx(tau[0], rel=1e-6)


def test_unresolved():
    pix = _pixels()  # b_inst = 4.0 km/s
    assert "unresolved" in measurement_flags(_comp(b=1.9), pix)
    assert "unresolved" not in measurement_flags(_comp(b=2.1), pix)
    # A resolution of the component's own overrides the spectrum's
    assert "unresolved" not in measurement_flags(_comp(b=1.9, resol=100000.0), pix)
    # No resolution: no flag
    nan = dataclasses.replace(pix, resol=np.full(pix.wave.size, np.nan))
    assert "unresolved" not in measurement_flags(_comp(b=1.0), nan)


def test_saturated_on_the_weakest_line():
    """The weaker line of CIV decides: a doublet whose weaker line is not
    saturated is not flagged; the resolution does not matter (the intrinsic
    profile)."""
    pix = _pixels()
    weak = ATOM_DATA["CIV_1550"]
    logN_sat = next(n for n in np.arange(12.0, 16.0, 0.01) if central_tau(weak, n, 8.0) > 3)
    assert "saturated" in measurement_flags(_comp(logN=logN_sat + 0.02, b=8.0), pix)
    assert "saturated" not in measurement_flags(_comp(logN=logN_sat - 0.02, b=8.0), pix)
    low_r = _pixels(resol=8900.0)
    assert "saturated" in measurement_flags(_comp(logN=logN_sat + 0.02, b=8.0), low_r)
    q, lines = flag_quantities(_comp(logN=14.0, b=8.0), pix)
    assert q.tau0_weakest == pytest.approx(central_tau(weak, 14.0, 8.0))
    assert len(lines) == 2


def test_saturated_counts_only_the_lines_in_the_fit():
    """With only CIV 1548 among the pixels of the fit, it is the weakest."""
    pix = _pixels()
    strong = ATOM_DATA["CIV_1548"]
    mask = np.abs(pix.wave - strong["wave"] * 3.97) < 2.0
    only_strong = dataclasses.replace(pix, fit_mask=mask)
    logN = 13.9  # CIV 1548 saturated, CIV 1550 not
    assert central_tau(strong, logN, 8.0) > 3 > central_tau(ATOM_DATA["CIV_1550"], logN, 8.0)
    assert "saturated" not in measurement_flags(_comp(logN=logN, b=8.0), pix)
    assert "saturated" in measurement_flags(_comp(logN=logN, b=8.0), only_strong)


def test_damped_HI_is_not_saturated():
    pix = _pixels(rest_nm=121.567, z=3.0)
    assert "saturated" not in measurement_flags(_comp("Ly_a", z=3.0, logN=20.3, b=20.0), pix)
    assert "saturated" in measurement_flags(_comp("Ly_a", z=3.0, logN=18.0, b=20.0), pix)


def test_large_errors():
    assert large_error_params(_comp(dlogN=1.5)) == ("logN",)
    assert large_error_params(_comp(b=5.0, db=6.0)) == ("b",)
    # The error of z in velocity above max(b, b_inst)
    dz = 12.0 / C_KMS * 3.97
    assert large_error_params(_comp(b=10.0, dz=dz)) == ("z",)
    assert large_error_params(_comp(b=10.0, dz=dz), b_inst=15.0) == ()
    # Not finite: nan_errors, not large_errors
    assert large_error_params(_comp(dlogN=5.0, dz=float("nan"))) == ()
    assert large_error_params(_comp(dlogN=5.0, dz=None)) == ()
    assert measurement_flags(_comp(dlogN=1.5), _pixels()) == ("large_errors",)


def test_telluric_overlap():
    centre = CIV * 10 * 3.97

    def band(depth, at):
        return lambda w: 1.0 - depth * np.exp(-0.5 * ((w - at) / 0.05) ** 2)

    assert "telluric_overlap" in measurement_flags(_comp(), _pixels(telluric=band(0.3, centre)))
    assert "telluric_overlap" not in measurement_flags(
        _comp(), _pixels(telluric=band(0.03, centre)))
    # Far from the lines (more than 2 max(b, b_inst) from their centres)
    far = centre * (1 + 300 / C_KMS)
    assert "telluric_overlap" not in measurement_flags(_comp(), _pixels(telluric=band(0.5, far)))
    # Without the column, no flag
    assert "telluric_overlap" not in measurement_flags(_comp(), _pixels())


def test_unknown_series_and_lines_outside():
    pix = _pixels()
    assert measurement_flags(_comp("nothing"), pix) == ()
    assert measurement_flags(_comp(z=1.0, b=1.0), pix) == ()  # outside the pixels
    assert measurement_flags(_comp("nothing", dlogN=3.0), pix) == ("large_errors",)


def test_with_measurement_flags():
    assert with_measurement_flags(("b_at_lower_bound", "saturated"), ("unresolved",)) == (
        "b_at_lower_bound", "unresolved")


def test_every_flag_has_a_name():
    """The methods text, the MCP glossary and the filter of components_list
    know every flag."""
    from astrocook.mcp_server.outputs import FLAG_GLOSSARY
    for flag in (*COMPONENT_FLAGS, *MEASUREMENT_FLAGS, "not_fitted"):
        assert flag in FLAG_NAMES
        assert f"'{flag}'" in FLAG_GLOSSARY


# --- Set by the fitter ------------------------------------------------------------

def _spec(logN, b, resol=45000.0, snr=50.0, z=2.97, telluric=None):
    lam = CIV * (1 + z)
    spec = generate(MockTruth(components=[MockComponent("CIV", z, logN, b)],
                              lam_min_nm=lam - 1.0, lam_max_nm=lam + 2.0, dv_kms=2.5,
                              resol=resol, snr=snr, z_em=z + 0.5, seed=4))
    if telluric is not None:
        x = spec.x.to("Angstrom").value
        data = spec._data
        cols = dict(data.aux_cols, telluric=DataColumn(telluric(x), None))
        spec = spec.__class__(dataclasses.replace(data, aux_cols=cols))
    return spec


def _fit(spec, logN, b, z=2.97):
    systs = SystemList(SystemListData(components=[ComponentData(
        id=1, z=z, dz=None, logN=logN, dlogN=None, b=b, db=None, series="CIV")]))
    fitted, _, _ = VoigtFitter(spec, systs).fit()
    return fitted


@pytest.mark.parametrize("logN, b, guess, expected", [
    (13.0, 1.5, (13.0, 5.0), ("unresolved",)),
    (15.0, 8.0, (14.0, 10.0), ("saturated",)),
    (13.5, 10.0, (13.0, 10.0), ()),
])
def test_the_fit_sets_the_flags(logN, b, guess, expected):
    fitted = _fit(_spec(logN, b), *guess)
    assert fitted.components[0].flags == expected


def test_the_fit_changes_no_number():
    """The flags are annotations: the values and errors are those the fit
    gave before 7d.2 (the same with the flags cleared)."""
    import astrocook.fitting.voigt_fitter as vf
    spec = _spec(15.0, 8.0)
    flagged = _fit(spec, 14.0, 10.0).components[0]
    original = vf.measurement_flags
    try:
        vf.measurement_flags = lambda comp, pix: ()
        plain = _fit(spec, 14.0, 10.0).components[0]
    finally:
        vf.measurement_flags = original
    assert flagged.flags == ("saturated",) and plain.flags == ()
    assert dataclasses.replace(flagged, flags=()) == dataclasses.replace(plain, uuid=flagged.uuid)


def test_telluric_overlap_in_the_fit():
    centre = CIV * 10 * 3.97

    def tell(x):
        return 1.0 - 0.2 * np.exp(-0.5 * ((x - centre) / 0.05) ** 2)

    fitted = _fit(_spec(13.5, 10.0, telluric=tell), 13.0, 10.0)
    assert fitted.components[0].flags == ("telluric_overlap",)


def test_evaluate_uncertainties_updates_the_flags():
    """The errors computed again at the end of a pipeline (#534) decide
    large_errors; the flags of #517 stay."""
    spec = _spec(13.5, 10.0)
    fitted = _fit(spec, 13.0, 10.0)
    comp = dataclasses.replace(fitted.components[0], dlogN=3.0,
                               flags=("z_at_window", "large_errors", "saturated"))
    systs = SystemList(dataclasses.replace(fitted._data, components=[comp]))
    updated = VoigtFitter(spec, systs).evaluate_uncertainties()
    assert updated.components[0].flags == ("z_at_window",)
    assert updated.components[0].dlogN < 1
