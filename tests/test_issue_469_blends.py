"""Issue #469: lines of different systems that are blended must be fitted together.

Scenario from the report: the CIV 1550 line of one system is blended with
AlII 1670 of another system at a different redshift. Selecting the AlII
component must bring the CIV component into the fitted group.
"""
import numpy as np
import astropy.units as au

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import (
    ComponentDataV2, DataColumnV2, SpectrumDataV2, SystemListDataV2,
)
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2

C_KMS = 299792.458
Z_ALII = 2.0
# CIV 1550 placed 5 km/s away from AlII 1670 at Z_ALII
_W_BLEND = ATOM_DATA["AlII_1670"]["wave"] * (1 + Z_ALII)
Z_CIV = _W_BLEND / ATOM_DATA["CIV_1550"]["wave"] * (1 + 5.0 / C_KMS) - 1


def _components(civ_logN=13.5):
    return [
        ComponentDataV2(id=1, z=Z_ALII, dz=None, logN=12.8, dlogN=None, b=8.0,
                        db=None, series="AlII_1670"),
        ComponentDataV2(id=2, z=Z_CIV, dz=None, logN=civ_logN, dlogN=None,
                        b=10.0, db=None, series="CIV"),
        # Unrelated line far away: must stay out of the group
        ComponentDataV2(id=3, z=3.5, dz=None, logN=13.0, dlogN=None, b=20.0,
                        db=None, series="Ly_a"),
    ]


def _systs(civ_logN=13.5):
    return SystemListV2(SystemListDataV2(components=_components(civ_logN)))


def test_blended_components_are_grouped():
    systs = _systs()
    alii, civ, lya = systs.components
    assert systs.get_connected_group([alii.uuid]) == {alii.uuid, civ.uuid}
    assert systs.get_connected_group([civ.uuid]) == {alii.uuid, civ.uuid}


def test_fitting_context_activates_blended_component():
    systs = _systs()
    alii, civ, lya = systs.components
    with systs.fitting_context([alii.uuid]):
        active = set(systs.constraint_model._active_uuids)
    assert active == {alii.uuid, civ.uuid}


def test_fitting_alii_also_fits_blended_civ():
    # Synthetic data from the true model; the CIV starting value is wrong.
    # If CIV is fitted together with AlII, its logN moves to the true value.
    x = np.arange(4990.0, 5025.0, 0.04) / 10.0  # nm, covers both CIV lines
    one = np.ones_like(x)

    def spectrum(y):
        return SpectrumV2(SpectrumDataV2(
            x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - 0.002, au.nm),
            xmax=DataColumnV2(x + 0.002, au.nm),
            y=DataColumnV2(y, au.dimensionless_unscaled),
            dy=DataColumnV2(0.01 * one, au.dimensionless_unscaled),
            aux_cols={"cont": DataColumnV2(one, au.dimensionless_unscaled)},
            z_em=4.0, resol=45000.0))

    _, true_model = VoigtFitterV2(spectrum(one), _systs()).compute_model_flux()
    data = spectrum(true_model + np.random.default_rng(1).normal(0, 0.01, len(x)))

    systs = _systs(civ_logN=13.0)  # true value 13.5
    alii = systs.components[0]
    with systs.fitting_context([alii.uuid]):
        fitted, _, result = VoigtFitterV2(data, systs).fit(max_nfev=2000)

    assert result.success
    civ_fit = fitted.components[1]
    assert abs(civ_fit.logN - 13.5) < 0.05


def test_inspector_highlights_blended_component(qapp):
    # The reporter's symptom: selecting AlII must highlight CIV (the group
    # members that are not selected get a yellow-orange background)
    from types import SimpleNamespace
    from PySide6.QtCore import Qt
    from astrocook.gui.system_inspector import SystemInspector

    systs = _systs()
    alii, civ, lya = systs.components
    inspector = SystemInspector(None)
    inspector.current_session = SimpleNamespace(systs=systs)
    inspector.table_model.update_data(systs.components, {})
    inspector._update_group_definition([alii])

    def background(row):
        model = inspector.table_model
        return model.data(model.index(row, 0), Qt.BackgroundRole)

    assert background(1) is not None   # CIV highlighted
    assert background(2) is None       # unrelated Ly-alpha not highlighted


def _pair_at(dv_kms):
    """AlII 1670 and CIV 1550 (b = 8 and 10 km/s) separated by dv_kms."""
    z_civ = _W_BLEND / ATOM_DATA["CIV_1550"]["wave"] * (1 + dv_kms / C_KMS) - 1
    comps = [ComponentDataV2(id=1, z=Z_ALII, dz=None, logN=13.5, dlogN=None,
                             b=8.0, db=None, series="AlII_1670"),
             ComponentDataV2(id=2, z=z_civ, dz=None, logN=14.0, dlogN=None,
                             b=10.0, db=None, series="CIV")]
    return SystemListV2(SystemListDataV2(components=comps))


def _flat_spectrum(resol, unit=au.dimensionless_unscaled):
    """Flat spectrum over both CIV lines, with a flat resol column."""
    x = np.arange(4990.0, 5025.0, 0.04) / 10.0  # nm
    one = np.ones_like(x)
    return SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - 0.002, au.nm),
        xmax=DataColumnV2(x + 0.002, au.nm),
        y=DataColumnV2(one, au.dimensionless_unscaled),
        dy=DataColumnV2(0.01 * one, au.dimensionless_unscaled),
        aux_cols={"resol": DataColumnV2(np.full_like(x, resol), unit)}))


def test_blend_threshold_without_resolution():
    # No spectrum (or no resolution): the threshold is max(4 b, 30 km/s)
    # (here 4 b = 40 km/s),
    # as before F14
    for dv, grouped in [(39.0, True), (41.0, False), (60.0, False)]:
        systs = _pair_at(dv)
        alii, civ = systs.components
        assert (civ.uuid in systs.get_connected_group([alii.uuid])) is grouped
        no_resol = _flat_spectrum(np.nan)
        assert (civ.uuid in systs.get_connected_group(
            [alii.uuid], spec=no_resol)) is grouped


def test_blend_threshold_follows_resolution():
    # F14: lines closer than the instrumental FWHM are blended. At R = 5000
    # (FWHM 60 km/s) the pair is grouped at 50 km/s, not at 70 km/s
    spec = _flat_spectrum(5000.0)
    for dv, grouped in [(50.0, True), (70.0, False)]:
        systs = _pair_at(dv)
        alii, civ = systs.components
        assert (civ.uuid in systs.get_connected_group(
            [alii.uuid], spec=spec)) is grouped
        with systs.fitting_context([alii.uuid], spec=spec):
            active = set(systs.constraint_model._active_uuids)
        assert (civ.uuid in active) is grouped


def test_blend_threshold_fwhm_column_and_floor():
    # A resol column in km/s is a FWHM; a FWHM below 30 km/s (high
    # resolution) leaves the floor max(4 b, 30 km/s) = 40 km/s in place
    systs = _pair_at(50.0)
    alii, civ = systs.components
    assert civ.uuid in systs.get_connected_group(
        [alii.uuid], spec=_flat_spectrum(60.0, au.km / au.s))
    for dv, grouped in [(39.0, True), (41.0, False)]:
        systs = _pair_at(dv)
        alii, civ = systs.components
        assert (civ.uuid in systs.get_connected_group(
            [alii.uuid], spec=_flat_spectrum(45000.0))) is grouped


def test_blend_threshold_component_resolution_override():
    # A component resolution (an override, R) wins over the column, as in
    # the fitter; the larger FWHM of the two lines is used
    spec = _flat_spectrum(45000.0)
    systs = _pair_at(50.0)
    alii, civ = systs.components
    assert civ.uuid not in systs.get_connected_group([alii.uuid], spec=spec)
    systs = systs.update_component(civ.uuid, resol=5000.0)
    assert civ.uuid in systs.get_connected_group([alii.uuid], spec=spec)


def test_blend_threshold_same_series_follows_resolution():
    # Kinematic proximity of the same species uses the same threshold
    comps = [ComponentDataV2(id=i, z=2.0 + dz, dz=None, logN=13.0, dlogN=None,
                             b=8.0, db=None, series="Ly_a")
             for i, dz in [(1, 0.0), (2, 3.0 * 50.0 / C_KMS)]]
    systs = SystemListV2(SystemListDataV2(components=comps))
    a, b = systs.components
    x = np.linspace(3600.0, 3660.0, 2000) / 10.0
    one = np.ones_like(x)
    spec = SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x, au.nm),
        xmax=DataColumnV2(x, au.nm), y=DataColumnV2(one, au.dimensionless_unscaled),
        dy=DataColumnV2(one, au.dimensionless_unscaled),
        aux_cols={"resol": DataColumnV2(np.full_like(x, 5000.0),
                                        au.dimensionless_unscaled)}))
    assert b.uuid not in systs.get_connected_group([a.uuid])
    assert b.uuid in systs.get_connected_group([a.uuid], spec=spec)


def test_fit_includes_component_blended_at_low_resolution():
    # F14 at fit level (the mock benchmark cannot see it: it fits only the
    # sparse metal doublets): at R = 5000 (FWHM 60 km/s) CIV 1550 lies
    # 50 km/s from AlII 1670, beyond max(4 b, 30 km/s) = 40 km/s but within
    # the FWHM. Fitting AlII with the spectrum passed to the grouping also
    # fits CIV (its wrong starting logN moves to the true value)
    x = np.arange(4990.0, 5025.0, 0.04) / 10.0  # nm
    one = np.ones_like(x)

    def spectrum(y):
        return SpectrumV2(SpectrumDataV2(
            x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - 0.002, au.nm),
            xmax=DataColumnV2(x + 0.002, au.nm),
            y=DataColumnV2(y, au.dimensionless_unscaled),
            dy=DataColumnV2(0.01 * one, au.dimensionless_unscaled),
            aux_cols={"cont": DataColumnV2(one, au.dimensionless_unscaled),
                      "resol": DataColumnV2(np.full_like(x, 5000.0),
                                            au.dimensionless_unscaled)},
            z_em=4.0))

    _, true_model = VoigtFitterV2(spectrum(one), _pair_at(50.0)).compute_model_flux()
    data = spectrum(true_model + np.random.default_rng(1).normal(0, 0.01, len(x)))

    for spec, civ_fitted in [(None, False), (data, True)]:
        systs = _pair_at(50.0)
        alii, civ = systs.components
        systs = systs.update_component(civ.uuid, logN=13.5)  # true value 14.0
        with systs.fitting_context([alii.uuid], spec=spec):
            fitted, _, result = VoigtFitterV2(data, systs).fit(max_nfev=2000)
        assert result.success
        assert bool(abs(fitted.components[1].logN - 13.5) > 0.1) is civ_fitted
