"""Sessions whose pixel edges have no unit (#623): since #526 the fit
converts the edges to Angstrom, and a refit failed with "'' (dimensionless)
and 'Angstrom' (length) are not convertible"."""
import dataclasses

import astropy.units as au
import numpy as np

from astrocook.core.session import SessionV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import (
    ComponentDataV2, DataColumnV2, SystemListDataV2,
)
from astrocook.core.system_list import SystemListV2
from astrocook.mocks import MockComponent, MockTruth, generate


def _session(scale=1.0):
    """A mock with a CIV doublet, the edges stored without a unit (and
    multiplied by ``scale``), as in some sessions of earlier versions."""
    spec = generate(MockTruth(
        components=[MockComponent("CIV", 2.97, 13.5, 10.0)], lam_min_nm=612.0,
        lam_max_nm=618.0, dv_kms=2.5, resol=45000.0, snr=30.0, z_em=3.2, seed=1))
    d = spec._data
    data = dataclasses.replace(
        d, xmin=DataColumnV2(d.xmin.values * scale, au.dimensionless_unscaled),
        xmax=DataColumnV2(d.xmax.values * scale, au.dimensionless_unscaled))
    systs = SystemListV2(SystemListDataV2(components=[ComponentDataV2(
        id=1, z=2.97, dz=None, logN=13.3, dlogN=None, b=12.0, db=None,
        series="CIV")]))
    return SessionV2("old", None, spec=SpectrumV2(data), systs=systs), spec


def test_edges_without_a_unit_take_the_unit_of_x():
    session, original = _session()
    assert session.spec.xmin.unit == au.nm and session.spec.xmax.unit == au.nm
    np.testing.assert_array_equal(session.spec.xmin.value, original.xmin.value)


def test_refit_of_a_session_with_edges_without_a_unit():
    session, _ = _session()
    uuid = session.systs.components[0].uuid
    out = session.absorbers.refit_all(selected_uuids=[uuid])
    assert out != 0
    assert abs(out.systs.components[0].logN - 13.5) < 0.1


def test_edges_that_are_not_around_x_are_recomputed():
    # e.g. in Angstrom, while x is in nm
    session, original = _session(scale=10.0)
    assert session.spec.xmin.unit == au.nm
    np.testing.assert_allclose(session.spec.xmin.value, original.xmin.value, rtol=1e-9)
    np.testing.assert_allclose(session.spec.xmax.value, original.xmax.value, rtol=1e-9)


def test_acs2_round_trip(tmp_path):
    session, _ = _session()
    path = tmp_path / "old.acs2"
    session.save(str(path))
    reloaded = SessionV2.open_new(str(path), "old", None, "auto")
    assert reloaded.spec.xmin.unit == au.nm
    np.testing.assert_array_equal(reloaded.spec.xmin.value, session.spec.xmin.value)
