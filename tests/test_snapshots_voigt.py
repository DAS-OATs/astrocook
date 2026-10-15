"""Numerical snapshots of the Voigt model and fit.

Two inputs:

- a synthetic normalised spectrum with known components, to cover the
  resolution sources of the fitter (global value, per-pixel column,
  per-component value, none);
- the J0103-1305 window with the continuum and system list saved in the
  session as frozen inputs, and R = 45000 (see ``conftest.py``).
"""
import numpy as np
import astropy.units as au
import pytest

from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import (
    ComponentDataV2, DataColumnV2, SpectrumDataV2, SystemListDataV2,
)
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2

R_SYNTH = 45000.0
# CIV doublet (two blended components) and an unrelated Ly-alpha line
TRUE_COMPONENTS = [
    dict(series="CIV", z=2.9700, logN=13.6, b=12.0),
    dict(series="CIV", z=2.97015, logN=13.1, b=6.0),
    dict(series="Ly_a", z=4.0590, logN=13.8, b=25.0),
]


def _components(params, resol=None):
    return [ComponentDataV2(id=i + 1, z=p["z"], dz=None, logN=p["logN"],
                            dlogN=None, b=p["b"], db=None, series=p["series"],
                            resol=resol[i] if resol else None)
            for i, p in enumerate(params)]


def _system_list(params, resol=None):
    return SystemListV2(SystemListDataV2(components=_components(params, resol)))


def _synthetic_spectrum(y=None, resol=0.0, resol_col=None, seed=None):
    """Normalised spectrum, 612-620 nm, about 2.4 km/s per pixel."""
    x = np.arange(612.0, 620.0, 0.005)
    n = len(x)
    if y is None:
        y = np.ones(n)
    aux = {"cont": DataColumnV2(np.ones(n), au.dimensionless_unscaled)}
    if resol_col is not None:
        aux["resol"] = DataColumnV2(resol_col, au.dimensionless_unscaled)
    return SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - 0.0025, au.nm),
        xmax=DataColumnV2(x + 0.0025, au.nm),
        y=DataColumnV2(y, au.dimensionless_unscaled),
        dy=DataColumnV2(np.full(n, 0.02), au.dimensionless_unscaled),
        aux_cols=aux, z_em=4.1, resol=resol))


def _model(spec, systs):
    _, flux = VoigtFitterV2(spec, systs).compute_model_flux()
    return flux


def _n(spec):
    return len(spec.x.value)


RESOLUTION_CASES = {
    "global": lambda n: dict(resol=R_SYNTH),
    "column_flat": lambda n: dict(resol_col=np.full(n, R_SYNTH)),
    "none": lambda n: dict(),
}


@pytest.mark.parametrize("case", RESOLUTION_CASES)
def test_model_resolution_sources(snapshot, case):
    n = _n(_synthetic_spectrum())
    spec = _synthetic_spectrum(**RESOLUTION_CASES[case](n))
    snapshot.check(f"voigt_model_synth_{case}",
                   {"model": _model(spec, _system_list(TRUE_COMPONENTS))})


def test_model_global_and_flat_column_agree():
    # Same resolving power from the attribute or from a flat column
    n = _n(_synthetic_spectrum())
    systs = _system_list(TRUE_COMPONENTS)
    np.testing.assert_allclose(
        _model(_synthetic_spectrum(resol=R_SYNTH), systs),
        _model(_synthetic_spectrum(resol_col=np.full(n, R_SYNTH)), systs),
        rtol=1e-10)


def _nan_middle(n):
    col = np.full(n, R_SYNTH)
    col[n // 3: n // 2] = np.nan
    return col


@pytest.mark.parametrize("label, resol_col", [
    ("flat_not_round", lambda n: np.full(n, 45050.0)),
    ("varying", lambda n: np.linspace(40000.0, 50000.0, n)),
    ("with_nan", _nan_middle),
])
def test_model_column_resolution_is_not_zero(label, resol_col):
    # Regression test for F10: pixels were grouped by the rounded resolution
    # but selected by the raw value, so those with a non-round resolution
    # (and those with NaN) were never filled and the model was 0 there
    n = _n(_synthetic_spectrum())
    spec = _synthetic_spectrum(resol_col=resol_col(n))
    model = _model(spec, _system_list(TRUE_COMPONENTS))
    # Far from the lines the normalised model must be close to 1
    assert np.all(model > 0.1)


def test_model_column_resolution_is_rounded():
    # The kernel uses the resolution rounded to 100: a flat column at
    # R = 45050 gives the same model as one at 45000
    n = _n(_synthetic_spectrum())
    systs = _system_list(TRUE_COMPONENTS)
    np.testing.assert_allclose(
        _model(_synthetic_spectrum(resol_col=np.full(n, 45050.0)), systs),
        _model(_synthetic_spectrum(resol_col=np.full(n, 45000.0)), systs),
        rtol=1e-12)


def test_model_component_resolution(snapshot):
    # Per-component resolution, with a different global value: records which
    # source wins in the current code
    spec = _synthetic_spectrum(resol=R_SYNTH)
    systs = _system_list(TRUE_COMPONENTS, resol=[20000.0, 20000.0, None])
    snapshot.check("voigt_model_synth_component_resol",
                   {"model": _model(spec, systs)})


def _fit_arrays(systs, result, uuids=None):
    comps = [c for c in systs.components if uuids is None or c.uuid in uuids]
    out = {key: np.array([getattr(c, key) for c in comps], dtype=float)
           for key in ("z", "dz", "logN", "dlogN", "b", "db")}
    out["cost"] = np.float64(result.cost)
    return out


def test_fit_synthetic(snapshot):
    clean = _synthetic_spectrum(resol=R_SYNTH)
    true_model = _model(clean, _system_list(TRUE_COMPONENTS))
    rng = np.random.default_rng(42)
    data = _synthetic_spectrum(
        y=true_model + rng.normal(0.0, 0.02, len(true_model)), resol=R_SYNTH)

    # Start from a guess inside the default z window (20 km/s)
    guess = [dict(p, z=p["z"] + 5e-5, logN=p["logN"] - 0.3, b=p["b"] * 1.3)
             for p in TRUE_COMPONENTS]
    systs, _, result = VoigtFitterV2(data, _system_list(guess)).fit(
        max_nfev=2000, z_window_kms=20.0)

    assert result.success
    snapshot.check("voigt_fit_synth", _fit_arrays(systs, result), rtol=1e-6)


@pytest.fixture(scope="module")
def j0103_systs(j0103_session):
    return j0103_session.systs


def test_model_j0103(snapshot, j0103_window_for_fit, j0103_systs):
    spec = j0103_window_for_fit
    model_norm = _model(spec, j0103_systs) / spec.get_column("cont").value
    # float32 keeps the file small; the tolerance is set accordingly
    snapshot.check("voigt_model_j0103",
                   {"model_norm": model_norm.astype(np.float32)}, rtol=1e-6,
                   atol=1e-7)


J0103_CIV_Z = 3.22549  # isolated CIV system on the red side


def test_fit_j0103_civ(snapshot, j0103_window_for_fit, j0103_systs):
    target = min(j0103_systs.components,
                 key=lambda c: abs(c.z - J0103_CIV_Z) if c.series == "CIV" else np.inf)
    assert abs(target.z - J0103_CIV_Z) < 1e-4
    with j0103_systs.fitting_context([target.uuid], group_depth=1):
        fitter = VoigtFitterV2(j0103_window_for_fit, j0103_systs)
        active = set(fitter._constraints._active_uuids)
        systs, _, result = fitter.fit(max_nfev=2000, z_window_kms=20.0)
    assert result.success
    snapshot.check("voigt_fit_j0103_civ", _fit_arrays(systs, result, active),
                   rtol=1e-6)
