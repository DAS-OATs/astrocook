"""Errors of the linked parameters (issue #630).

A parameter linked to another (e.g. the z of SiIV to that of CIV, or its b
to the b of CIV times the ratio of the thermal widths) is not free in the
fit, and was given error 0, which the ASCII list of the systems exported.
Its error is now propagated from the free parameters through the
expression of the link. Data from the mock generator.
"""
import dataclasses

import numpy as np
import pytest

from astrocook.core.structures import ComponentDataV2, SystemListDataV2
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.mocks import MockComponent, MockTruth, generate

Z = 2.0
RATIO = 0.654  # b of SiIV over b of CIV (thermal)


@pytest.fixture(scope="module")
def spec():
    truth = MockTruth(
        components=[MockComponent("CIV", z=Z, logN=13.5, b=15.0),
                      MockComponent("SiIV", z=Z, logN=13.0, b=15.0 * RATIO)],
        lam_min_nm=415.0, lam_max_nm=470.0, dv_kms=2.5, resol=45000.0,
        snr=50.0, z_em=2.2, seed=3)
    return generate(truth)


def _linked_systs(b_expr="{}"):
    comps = [ComponentDataV2(id=1, z=Z + 2e-5, dz=None, logN=13.3, dlogN=None,
                             b=17.0, db=None, series="CIV"),
             ComponentDataV2(id=2, z=Z - 2e-5, dz=None, logN=12.8, dlogN=None,
                             b=8.0, db=None, series="SiIV")]
    systs = SystemListV2(SystemListDataV2(components=comps))
    civ, siiv = (c.uuid for c in systs.components)
    target = f"p['{civ}']"
    return (systs
            .update_constraint(siiv, "z", is_free=False, expression=f"{target}.z",
                               target_uuid=civ)
            .update_constraint(siiv, "b", is_free=False,
                               expression=b_expr.format(f"{target}.b"),
                               target_uuid=civ))


def test_linked_parameters_carry_the_error_of_their_target(spec):
    systs = _linked_systs(b_expr="{}*" + str(RATIO))
    fitted, _, _ = VoigtFitterV2(spec, systs).fit(max_nfev=2000, z_window_kms=20.0)
    civ, siiv = fitted.components
    # The links hold...
    assert siiv.z == civ.z
    assert siiv.b == pytest.approx(civ.b * RATIO, rel=1e-10)
    # ...and the errors follow them (before: 0)
    assert civ.dz > 0 and civ.db > 0
    assert siiv.dz == pytest.approx(civ.dz, rel=1e-6)
    assert siiv.db == pytest.approx(civ.db * RATIO, rel=1e-6)
    # The free parameter of the linked component keeps its own error
    assert siiv.dlogN > 0
    # The frozen btur still has error 0, and nothing is flagged
    assert siiv.dbtur == 0.0 and civ.dbtur == 0.0
    assert "nan_errors" not in siiv.flags
    # The table (used by the ASCII export) shows them
    t = fitted.t
    assert t["dz"][1] == siiv.dz and t["db"][1] == siiv.db


def test_nonlinear_link(spec):
    # A link that is not linear: the gradient of the expression, sqrt here
    systs = _linked_systs(b_expr="np.sqrt({})*2.5")
    fitted, _, _ = VoigtFitterV2(spec, systs).fit(max_nfev=2000, z_window_kms=20.0)
    civ, siiv = fitted.components
    assert siiv.b == pytest.approx(2.5 * np.sqrt(civ.b), rel=1e-10)
    assert siiv.db == pytest.approx(civ.db * 1.25 / np.sqrt(civ.b), rel=1e-5)


def test_a_parameter_linked_to_a_frozen_one_has_no_error(spec):
    systs = _linked_systs()
    civ_uuid = systs.components[0].uuid
    systs = systs.update_constraint(civ_uuid, "b", is_free=False)
    fitted, _, _ = VoigtFitterV2(spec, systs).fit(max_nfev=2000, z_window_kms=20.0)
    civ, siiv = fitted.components
    assert civ.db == 0.0 and siiv.db == 0.0
    assert siiv.dz == pytest.approx(civ.dz, rel=1e-6)


def test_evaluate_uncertainties_also_gives_the_linked_errors(spec):
    # The errors recomputed at the end of the automatic pipelines (#534)
    systs = _linked_systs(b_expr="{}*" + str(RATIO))
    fitted, _, _ = VoigtFitterV2(spec, systs).fit(max_nfev=2000, z_window_kms=20.0)
    zeroed = SystemListV2(dataclasses.replace(
        fitted._data, components=[dataclasses.replace(c, dz=0.0, db=0.0)
                                  for c in fitted.components]))
    civ, siiv = VoigtFitterV2(spec, zeroed).evaluate_uncertainties().components
    assert siiv.dz == pytest.approx(civ.dz, rel=1e-6) and siiv.dz > 0
    assert siiv.db == pytest.approx(civ.db * RATIO, rel=1e-6)
