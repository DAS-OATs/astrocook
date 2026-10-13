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
