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
