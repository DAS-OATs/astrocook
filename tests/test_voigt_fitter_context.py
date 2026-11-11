"""Tests for the preparation of the Voigt fit context."""
import numpy as np
import pytest
import astropy.units as au

from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import (
    ComponentData, DataColumn, SpectrumData, SystemListData
)
from astrocook.core.system_list import SystemList
from astrocook.fitting.voigt_fitter import NoPixelsToFit, VoigtFitter


def test_prepare_fit_context_with_no_valid_pixels():
    # With no valid pixel the fit mask stays empty; the debug branch that
    # lists the components used an unassigned variable and raised
    # UnboundLocalError. Then the whole spectrum was fitted; since F111
    # (phase 6-ter, PR 7-quater) the fit is refused with a clear message.
    x = np.linspace(4630.0, 4660.0, 300)
    one = np.ones_like(x)
    spec = Spectrum(SpectrumData(
        x=DataColumn(x, au.Angstrom),
        xmin=DataColumn(x - 0.05, au.Angstrom),
        xmax=DataColumn(x + 0.05, au.Angstrom),
        y=DataColumn(np.full_like(x, np.nan), au.dimensionless_unscaled),
        dy=DataColumn(0.1 * one, au.dimensionless_unscaled),
        aux_cols={'cont': DataColumn(one, au.dimensionless_unscaled)},
    ))
    comp = ComponentData(id=1, z=2.0, dz=None, logN=13.5, dlogN=None,
                           b=20.0, db=None, series='CIV')
    fitter = VoigtFitter(spec, SystemList(SystemListData(components=[comp])))

    with pytest.raises(NoPixelsToFit, match=r"\(CIV at z = 2\.0000\): .* or on pixels "
                                             "that cannot be used"):
        fitter.prepare_fit_context()
