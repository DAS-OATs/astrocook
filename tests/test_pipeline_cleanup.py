"""Output of the automatic pipelines (#517, F21): no Ly-a narrower than
min_b after the last refit of lya_auto, no components at the lower bound
of logN (they absorb nothing measurable)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import full_chain as fc  # noqa: E402


@pytest.fixture(scope="module")
def forest():
    fc.seed_uuids(0)
    session = fc.mock_session().continuum.estimate_auto()
    return session.absorbers.lya_auto()


def test_no_lya_narrower_than_min_b(forest):
    # Before #517 the filter on min_b ran before the last refit, which
    # could narrow other lines below it (39 on J0103-1305)
    narrow = [c for c in forest.systs.components if c.series == "Ly_a" and c.b < 10.0]
    assert not narrow


def test_no_component_at_the_lower_bound_of_logN(forest):
    assert not [c for c in forest.systs.components if "logN_at_lower_bound" in c.flags]
    assert all(c.logN > 10.01 for c in forest.systs.components)


def test_forest_auto_applies_min_b():
    # F65: forest_auto took min_b and did not use it. A short stretch of
    # forest (the whole mock of the full chain takes minutes with the
    # Lyman series)
    from astrocook.core.session import SessionV2
    from astrocook.mocks import MockComponent, MockTruth, generate
    z_em = 3.0
    comps = [MockComponent("Ly_ab", z, logN, b) for z, logN, b in (
        (2.80, 14.0, 25.0), (2.8015, 13.5, 6.0), (2.82, 13.8, 30.0), (2.85, 13.2, 4.0))]
    spec = generate(MockTruth(components=comps, lam_min_nm=385.0, lam_max_nm=470.0,
                              dv_kms=5.0, resol=45000.0, snr=30.0, z_em=z_em, seed=5))
    session = SessionV2("forest", None, spec=spec).continuum.estimate_auto()
    out = session.absorbers.forest_auto(min_b="10.0")
    assert out != 0 and out.systs.components
    assert not [c for c in out.systs.components if c.series == "Lyman" and c.b < 10.0]
