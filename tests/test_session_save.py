"""Saving a session: only the V2 ``.acs2`` format (phase 3).

Saving in the V1 ``.acs`` format was removed (agreed with Guido): V1
archives can still be opened, and Astrocook 1.9.0 can still write them.
"""
import os

import numpy as np
import astropy.units as au
import pytest

from astrocook.session.session import Session, load_session_from_file
from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, HistoryLog, SpectrumData


def _make_session():
    x = np.linspace(400.0, 410.0, 50)
    one = np.ones_like(x)
    data = SpectrumData(
        x=DataColumn(x, au.nm),
        xmin=DataColumn(x - 0.1, au.nm),
        xmax=DataColumn(x + 0.1, au.nm),
        y=DataColumn(one, au.dimensionless_unscaled),
        dy=DataColumn(0.1 * one, au.dimensionless_unscaled),
    )
    session = Session(name="save_test", gui=None, spec=Spectrum(data))
    log = HistoryLog()
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
