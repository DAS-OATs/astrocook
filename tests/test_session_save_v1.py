"""Tests for saving a session to the V1 .acs format."""
import json
import logging
import tarfile

import numpy as np
import astropy.units as au

from astrocook.core.session import SessionV2
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
    session = SessionV2(name="v1_save_test", gui=None, spec=SpectrumV2(data))
    log = HistoryLogV2()
    log.add_entry("smooth", {"std": 10.0})
    session.log_manager = log
    return session


def _read_log_from_archive(path):
    with tarfile.open(path) as tar:
        names = [m.name for m in tar.getmembers() if m.name.endswith('.json')]
        assert names, "no JSON log in the archive"
        return json.loads(tar.extractfile(names[0]).read())


def test_v2_log_is_not_stored_in_v1_archive(tmp_path, caplog):
    session = _make_session()
    path = tmp_path / "session.acs"

    assert session.log_lost_on_save(str(path)) is True
    with caplog.at_level(logging.WARNING):
        session.save(str(path), initial_session=session)

    assert _read_log_from_archive(path) == {"set_menu": []}
    assert any(".acs2" in r.getMessage() for r in caplog.records
               if r.levelno == logging.WARNING)


def test_v2_log_is_kept_in_acs2():
    session = _make_session()
    assert session.log_lost_on_save("session.acs2") is False
    assert session.log_lost_on_save("SESSION.ACS2") is False
