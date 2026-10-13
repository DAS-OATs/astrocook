"""Issue #455: fit only a subset of the transitions of a large multiplet.

The inspector offers "Add <lines> (Subset)" for multiplets with more than
two transitions; it calls the ``add_component`` recipe with a
comma-separated series of the selected transitions. The component must
model those lines only, and be fitted.
"""
import numpy as np
import astropy.units as au
import pytest

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS
from astrocook.core.session import SessionV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import (
    ComponentDataV2, DataColumnV2, SpectrumDataV2, SystemListDataV2,
)
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.recipes.absorbers import RecipeAbsorbersV2

Z = 1.0
SUBSET = "FeII_2344,FeII_2382"


def _spectrum(x_nm, y):
    one = np.ones_like(x_nm)
    return SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x_nm, au.nm), xmin=DataColumnV2(x_nm - 0.0015, au.nm),
        xmax=DataColumnV2(x_nm + 0.0015, au.nm),
        y=DataColumnV2(y, au.dimensionless_unscaled),
        dy=DataColumnV2(0.01 * one, au.dimensionless_unscaled),
        aux_cols={"cont": DataColumnV2(one, au.dimensionless_unscaled)},
        z_em=2.0, resol=45000.0))


def _line_depth(x_nm, model, trans):
    w = ATOM_DATA[trans]["wave"] * (1 + Z) / 10.0
    near = np.abs(x_nm - w) < 0.05
    return 1.0 - model[near].min()


def test_niII_offers_subset():
    # Reported in the issue: NiII did not get the subset option
    assert len(STANDARD_MULTIPLETS["NiII"]) > 2


@pytest.fixture
def restore_multiplets():
    # add_component registers custom series in the global STANDARD_MULTIPLETS
    saved = dict(STANDARD_MULTIPLETS)
    yield
    STANDARD_MULTIPLETS.clear()
    STANDARD_MULTIPLETS.update(saved)


def test_add_subset_component_models_and_fits_only_selected_lines(restore_multiplets):
    x = np.arange(466.0, 522.0, 0.003)  # FeII 2344-2600 at z = 1
    true = ComponentDataV2(id=1, z=Z, dz=None, logN=13.8, dlogN=None, b=8.0,
                           db=None, series=SUBSET)
    _, model = VoigtFitterV2(
        _spectrum(x, np.ones_like(x)),
        SystemListV2(SystemListDataV2(components=[true]))).compute_model_flux()
    data = model + np.random.default_rng(3).normal(0, 0.01, len(x))
    session = SessionV2(name="feii", gui=None, spec=_spectrum(x, data))

    # Start from a wrong logN, as a user clicking on the plot would
    result = RecipeAbsorbersV2(session).add_component(
        series=SUBSET, z=f"{Z}", logN="13.3", b="10.0")
    assert result != 0

    comps = result.systs.components
    assert len(comps) == 1 and comps[0].series == SUBSET
    assert comps[0].logN == pytest.approx(13.8, abs=0.05)

    fitted = result.spec.get_column("model").value
    assert _line_depth(x, fitted, "FeII_2344") > 0.1
    assert _line_depth(x, fitted, "FeII_2382") > 0.1
    # Transitions of the multiplet that were not selected are not modelled
    for trans in ("FeII_2374", "FeII_2586", "FeII_2600"):
        assert _line_depth(x, fitted, trans) < 1e-3
