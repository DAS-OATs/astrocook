"""The resolution as a single source, the ``resol`` column (#526): units,
migration of the legacy sources, ``with_properties``, the component
overrides, and the operations that must keep it."""
import dataclasses

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.resolution import (
    C_KMS, column_to_resolving_power, migrate_component_resolutions,
)
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import (
    ComponentDataV2, DataColumnV2, SpectrumDataV2, SystemListDataV2,
)
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2


def _spectrum(resol=0.0, meta=None, resol_col=None, resol_unit=au.dimensionless_unscaled,
              y=None):
    """Normalised spectrum, 612-620 nm, ~2.4 km/s per pixel."""
    x = np.arange(612.0, 620.0, 0.005)
    n = len(x)
    unit = au.dimensionless_unscaled
    aux = {"cont": DataColumnV2(np.ones(n), unit)}
    if resol_col is not None:
        aux["resol"] = DataColumnV2(np.asarray(resol_col, dtype=float), resol_unit)
    return SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - 0.0025, au.nm),
        xmax=DataColumnV2(x + 0.0025, au.nm),
        y=DataColumnV2(np.ones(n) if y is None else y, unit),
        dy=DataColumnV2(np.full(n, 0.02), unit), aux_cols=aux,
        meta=dict(meta or {}), z_em=4.1, resol=resol))


N = len(np.arange(612.0, 620.0, 0.005))


# --- Units -------------------------------------------------------------------

def test_column_units():
    assert column_to_resolving_power([45000.0], au.dimensionless_unscaled)[0] == 45000.0
    assert column_to_resolving_power([C_KMS / 45000], au.km / au.s)[0] == pytest.approx(45000.0)
    assert column_to_resolving_power([6.66], au.m / au.s)[0] == pytest.approx(C_KMS / 6.66e-3)
    assert np.isnan(column_to_resolving_power([0.0, -1.0, np.nan], None)).all()


def test_no_resolution():
    assert np.isnan(_spectrum().resolution()).all()


# --- Migration of the legacy sources -----------------------------------------

@pytest.mark.parametrize("kwargs, expected", [
    (dict(resol=45000.0), 45000.0),                         # attribute, R
    (dict(resol=6.0), C_KMS / 6.0),                         # attribute below 500: km/s
    (dict(meta={"resol": 40000.0}), 40000.0),               # meta 'resol'
    (dict(meta={"RESOL": 30000.0}), 30000.0),               # meta 'RESOL' (V1)
    (dict(meta={"resol_fwhm": 7.5}), C_KMS / 7.5),          # meta 'resol_fwhm' (km/s)
    (dict(resol=45000.0, meta={"resol": 1000.0}), 45000.0),  # attribute first
    (dict(resol=45000.0, resol_col=np.full(N, np.nan)), 45000.0),  # all-NaN column replaced
])
def test_legacy_sources_migrate_into_the_column(kwargs, expected):
    spec = _spectrum(**kwargs)
    col = spec._data.aux_cols["resol"]
    assert spec.resolution() == pytest.approx(np.full(N, expected))
    assert col.unit == au.dimensionless_unscaled


def test_the_attribute_is_cleared_by_the_migration():
    # PR 4: once in the column, the attribute is cleared, so that deleting
    # the column removes the resolution for good
    spec = _spectrum(resol=45000.0)
    assert spec._data.resol == 0.0
    assert np.isnan(spec.remove_columns(["resol"]).resolution()).all()


def test_a_valid_column_wins_over_the_legacy_sources():
    variable = np.linspace(40000.0, 50000.0, N)
    spec = _spectrum(resol=10000.0, meta={"resol": 20000.0}, resol_col=variable)
    assert spec.resolution() == pytest.approx(variable)


def test_velocity_column_and_legacy_dimensionless_fwhm():
    fwhm = np.full(N, 6.0)
    assert _spectrum(resol_col=fwhm, resol_unit=au.km / au.s).resolution() == pytest.approx(C_KMS / fwhm)
    # Saved before #526 as a dimensionless column of FWHM values
    spec = _spectrum(resol_col=fwhm)
    assert spec._data.aux_cols["resol"].unit == au.km / au.s
    assert spec.resolution() == pytest.approx(C_KMS / fwhm)


# --- with_properties -------------------------------------------------------------

def test_with_properties_sets_a_flat_column():
    spec = _spectrum(resol=45000.0).with_properties(z_em=4.1, resol=30000.0)
    assert spec.resolution() == pytest.approx(np.full(N, 30000.0))


def test_with_properties_keeps_a_variable_column_by_default():
    variable = np.linspace(40000.0, 50000.0, N)
    spec = _spectrum(resol_col=variable)
    assert spec.with_properties(z_em=4.1, resol=30000.0).resolution() == pytest.approx(variable)
    flat = spec.with_properties(z_em=4.1, resol=30000.0, overwrite_resol_col=True)
    assert flat.resolution() == pytest.approx(np.full(N, 30000.0))


def test_with_properties_without_resolution_leaves_it():
    variable = np.linspace(40000.0, 50000.0, N)
    spec = _spectrum(resol_col=variable).with_properties(z_em=3.0)
    assert spec.resolution() == pytest.approx(variable)
    assert spec._data.z_em == 3.0


def test_set_properties_current_keeps_a_variable_column():
    # Before #526 '_current_' wrote the attribute into a flat column
    from astrocook.core.session import SessionV2
    variable = np.linspace(40000.0, 50000.0, N)
    session = SessionV2("s", None, spec=_spectrum(resol=45000.0, resol_col=variable))
    out = session.edit.set_properties(z_em="3.0")
    assert out.spec.resolution() == pytest.approx(variable)


# --- Component overrides -----------------------------------------------------------

def _systs(*resols):
    comps = [ComponentDataV2(id=i + 1, z=2.97 + 0.001 * i, dz=None, logN=13.5,
                             dlogN=None, b=10.0, db=None, series="CIV", resol=r)
             for i, r in enumerate(resols)]
    return SystemListV2(SystemListDataV2(components=comps))


def test_stamped_component_resolutions_are_cleared():
    spec = _spectrum(resol=45000.0)
    systs = migrate_component_resolutions(
        _systs(45000.0, 45300.0, 30000.0, C_KMS / 45000.0, None), spec)
    resols = [c.resol for c in systs.components]
    # within 1% of the spectrum (stamped), a deliberate override, a legacy
    # FWHM equal to the spectrum (stamped too), none
    assert resols[0] is None and resols[1] is None
    assert resols[2] == 30000.0
    assert resols[3] is None
    assert resols[4] is None


def test_legacy_component_fwhm_is_converted():
    systs = migrate_component_resolutions(_systs(10.0), _spectrum(resol=45000.0))
    assert systs.components[0].resol == pytest.approx(C_KMS / 10.0)


def test_the_fitter_does_not_stamp_the_components():
    spec = _spectrum(resol=45000.0)
    x = spec.x.value
    y = np.ones(N)
    for lam in (154.8204 * 3.97, 155.0781 * 3.97):
        y -= 0.5 * np.exp(-0.5 * ((x - lam) / 0.01) ** 2)
    spec = _spectrum(resol=45000.0, y=y)
    fitter = VoigtFitterV2(spec, _systs(None, 30000.0))
    systs, *_ = fitter.fit(max_nfev=20)
    assert [c.resol for c in systs.components] == [None, 30000.0]


# --- Operations ---------------------------------------------------------------

def test_rebin_returns_resolving_power():
    spec = _spectrum(resol_col=np.full(N, 6.0), resol_unit=au.km / au.s)
    reb = spec.rebin(None, None, 5.0 * au.km / au.s, None, np.nan)
    assert reb._data.aux_cols["resol"].unit == au.dimensionless_unscaled
    r = reb.resolution()
    assert np.nanmedian(r) == pytest.approx(C_KMS / 6.0, rel=1e-3)


@pytest.mark.parametrize("op", [
    lambda s: s.apply_expression("y", "y * 2"),
    lambda s: s.smooth_column("y", 10.0),
    lambda s: s.x_convert(4.1, au.km / au.s),
])
def test_operations_keep_the_resolution(op):
    # F22: these methods dropped the attribute; the column survives
    variable = np.linspace(40000.0, 50000.0, N)
    out = op(_spectrum(resol_col=variable))
    assert out.resolution() == pytest.approx(variable)


def test_dataclass_replace_keeps_the_column():
    spec = _spectrum(resol=45000.0)
    data = dataclasses.replace(spec._data, resol=0.0)
    assert SpectrumV2(data).resolution() == pytest.approx(np.full(N, 45000.0))


# --- GUI on the same resolver (PR 2) ------------------------------------------------

def test_resolution_summary():
    from astrocook.core.resolution import resolution_summary
    assert resolution_summary(_spectrum()) == (0.0, False)
    assert resolution_summary(_spectrum(resol=45000.0)) == (45000.0, False)
    median, variable = resolution_summary(_spectrum(resol_col=np.linspace(40000.0, 50000.0, N)))
    assert median == pytest.approx(45000.0, rel=1e-3) and variable


@pytest.mark.parametrize("col, unit, expected", [
    (np.full(N, 45000.0), au.dimensionless_unscaled, 45000.0),
    (np.full(N, 6.0), au.km / au.s, C_KMS / 6.0),
])
def test_system_inspector_reads_the_column_by_its_unit(qtbot, col, unit, expected):
    # Before #526 the inspector always read the column as km/s
    from types import SimpleNamespace
    from astrocook.core.session import SessionV2
    from astrocook.gui.system_inspector import VelocityPlotWidget
    holder = SimpleNamespace(_current_session=SessionV2(
        "s", None, spec=_spectrum(resol_col=col, resol_unit=unit)))
    r, unit_out = VelocityPlotWidget._get_resolution_at(holder, 616.0)
    assert unit_out == 'R' and r == pytest.approx(expected)
