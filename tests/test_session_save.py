"""Saving a session: only the V2 ``.acs2`` format (phase 3).

Saving in the V1 ``.acs`` format was removed (agreed with Guido): V1
archives can still be opened, and Astrocook 1.9.0 can still write them.
"""
import os

import numpy as np
import astropy.units as au
import pytest

from astrocook.core.session import SessionV2, load_session_from_file
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, HistoryLogV2, SpectrumDataV2


def _make_session():
    x = np.linspace(400.0, 410.0, 50)
    one = np.ones_like(x)
    data = SpectrumDataV2(
        x=DataColumnV2(x, au.nm),
        xmin=DataColumnV2(x - 0.1, au.nm),
        xmax=DataColumnV2(x + 0.1, au.nm),
        y=DataColumnV2(one, au.dimensionless_unscaled),
        dy=DataColumnV2(0.1 * one, au.dimensionless_unscaled),
    )
    session = SessionV2(name="save_test", gui=None, spec=SpectrumV2(data))
    log = HistoryLogV2()
    log.add_entry("smooth", {"std": 10.0})
    session.log_manager = log
    return session


@pytest.mark.parametrize("name", ["session.acs", "SESSION.ACS"])
def test_saving_to_acs_is_refused(tmp_path, name):
    path = tmp_path / name
    with pytest.raises(ValueError, match=r"\.acs2.*1\.9\.0"):
        _make_session().save(str(path))
    assert os.listdir(tmp_path) == []


def test_save_adds_the_acs2_extension(tmp_path):
    session = _make_session()
    session.save(str(tmp_path / "session"), initial_session=session)
    assert os.listdir(tmp_path) == ["session.acs2"]
    loaded = load_session_from_file(str(tmp_path / "session.acs2"), "s", None, None)
    assert [e.recipe_name for e in loaded.log_manager.entries] == ["smooth"]
