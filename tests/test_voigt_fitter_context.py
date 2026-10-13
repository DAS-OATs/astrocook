"""Tests for the preparation of the Voigt fit context."""
import numpy as np
import astropy.units as au

from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import (
    ComponentDataV2, DataColumnV2, SpectrumDataV2, SystemListDataV2
)
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2


def test_prepare_fit_context_with_no_valid_pixels():
    # With no valid pixel the fit mask stays empty; the debug branch that
    # lists the components used an unassigned variable and raised
    # UnboundLocalError.
    x = np.linspace(4630.0, 4660.0, 300)
    one = np.ones_like(x)
    spec = SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.Angstrom),
        xmin=DataColumnV2(x - 0.05, au.Angstrom),
        xmax=DataColumnV2(x + 0.05, au.Angstrom),
        y=DataColumnV2(np.full_like(x, np.nan), au.dimensionless_unscaled),
        dy=DataColumnV2(0.1 * one, au.dimensionless_unscaled),
        aux_cols={'cont': DataColumnV2(one, au.dimensionless_unscaled)},
    ))
    comp = ComponentDataV2(id=1, z=2.0, dz=None, logN=13.5, dlogN=None,
                           b=20.0, db=None, series='CIV')
    fitter = VoigtFitterV2(spec, SystemListV2(SystemListDataV2(components=[comp])))

    lower, upper = fitter.prepare_fit_context()

    assert len(lower) == len(upper) == 3  # z, logN, b (btur fixed)
    assert np.all(lower < upper)
