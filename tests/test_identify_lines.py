"""``SpectrumV2.identify_lines`` on mocks with a known answer."""
import json

import numpy as np
import pytest

from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.mocks.physics import C_KMS

Z_EM = 3.0

#: Isolated metal components on the red side (Ly-a at 486 nm), with a
#: two-component CIV system (40 km/s apart) to check that real close
#: components are kept apart
COMPONENTS = [
    MockComponent("CIV", 2.50, 14.0, 10.0),
    MockComponent("CIV", 2.70, 13.8, 12.0),
    MockComponent("CIV", 2.70 + 40.0 / C_KMS * 3.70, 13.6, 10.0),
    MockComponent("SiIV", 2.65, 13.5, 10.0),
    MockComponent("MgII", 1.20, 13.5, 8.0),
]


@pytest.fixture(scope="module")
def identified():
    truth = MockTruth(components=COMPONENTS, lam_min_nm=490.0, lam_max_nm=700.0,
                      dv_kms=2.5, resol=45000.0, snr=100.0, z_em=Z_EM, seed=1)
    return generate(truth).identify_lines(["CIV", "SiIV", "MgII"])


def _systems(spec):
    series = json.loads(spec.meta["series_map_json"])
    z = json.loads(spec.meta["z_map_json"])
    return [(series[k], z[k]) for k in series]


@pytest.mark.parametrize("series", ["CIV", "SiIV", "MgII"])
def test_one_candidate_per_component(identified, series):
    # F49: a doublet gave one candidate per line, ~1 km/s apart
    found = sorted(z for s, z in _systems(identified) if s == series)
    true = sorted(c.z for c in COMPONENTS if c.series == series)
    assert len(found) == len(true)
    for z_found, z_true in zip(found, true):
        assert abs(z_found - z_true) / (1 + z_true) * C_KMS < 5.0


def test_both_lines_of_a_doublet_are_labelled(identified):
    # Both regions of a doublet keep their label (the GUI shows them); only
    # the region of the stronger line gives the candidate systems
    lines = {line for entries in
             json.loads(identified.meta["region_identifications"]).values()
             for line, _ in entries}
    assert {"CIV_1548", "CIV_1550", "SiIV_1393", "SiIV_1402",
            "MgII_2796", "MgII_2803"} <= lines
    assert np.count_nonzero(identified.get_column("abs_ids").value) > 0
