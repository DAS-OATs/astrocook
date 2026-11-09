"""A save that fails leaves the previous archive in place (F94, phase 6).

The ``.acs2`` archive used to be written directly over the target file:
``tarfile.open(path, 'w:gz')`` truncates it before its members are added,
so a failure while writing left a partial archive (or an empty file) in
place of the previous one. A failure in the original data
(``initial_session``) was logged and the save went on, failing later with
``FileNotFoundError`` after the target had been truncated. Now the archive
is written to a temporary file in the same folder and replaces the target
only when complete; any failure leaves the target as it was and no
temporary file behind.
"""
import os
import stat
import tarfile

import pytest

import astrocook.io.adapter as adapter
from astrocook.core.session import SessionV2

from test_recipes_absorbers import _systs
from test_recipes_edit_flux import _spec


@pytest.fixture
def session():
    return SessionV2("s", None, spec=_spec(),
                     systs=_systs(("CIV", 2.97, 13.6, 12.0)))


def _names(path):
    with tarfile.open(path) as tar:
        return sorted(tar.getnames())


def _fail_on_call(monkeypatch, obj, name, n, exc=OSError("simulated failure")):
    """Make the ``n``-th call of ``obj.name`` raise ``exc``."""
    original = getattr(obj, name)
    calls = {"n": 0}

    def wrapper(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == n:
            raise exc
        return original(*args, **kwargs)

    monkeypatch.setattr(obj, name, wrapper)


def test_save_writes_the_archive(session, tmp_path):
    path = tmp_path / "s.acs2"
    session.save(str(path), initial_session=session)
    # The metadata of the original data since F109
    assert _names(path) == ["s_meta.json", "s_meta_orig.json", "s_spec.fits",
                            "s_spec_orig.fits", "s_systs.fits", "s_systs_orig.fits"]
    assert os.listdir(tmp_path) == ["s.acs2"]


def test_failure_while_writing_keeps_the_previous_archive(session, tmp_path,
                                                          monkeypatch):
    path = tmp_path / "s.acs2"
    session.save(str(path))
    before = path.read_bytes()
    # The second member fails while the archive is being written
    _fail_on_call(monkeypatch, tarfile.TarFile, "add", 2)
    with pytest.raises(OSError, match="simulated failure"):
        session.save(str(path), initial_session=session)
    assert path.read_bytes() == before
    assert os.listdir(tmp_path) == ["s.acs2"]


def test_failure_of_the_original_data_keeps_the_previous_archive(
        session, tmp_path, monkeypatch):
    path = tmp_path / "s.acs2"
    session.save(str(path))
    before = path.read_bytes()
    # The second conversion of a spectrum is that of the original data
    _fail_on_call(monkeypatch, adapter, "_convert_spec_data_to_table", 2,
                  RuntimeError("original data"))
    with pytest.raises(RuntimeError, match="original data"):
        session.save(str(path), initial_session=session)
    assert path.read_bytes() == before
    assert os.listdir(tmp_path) == ["s.acs2"]


def test_failed_save_of_a_new_file_leaves_nothing(session, tmp_path,
                                                  monkeypatch):
    _fail_on_call(monkeypatch, tarfile.TarFile, "add", 1)
    with pytest.raises(OSError):
        session.save(str(tmp_path / "s.acs2"))
    assert os.listdir(tmp_path) == []


def test_new_archive_has_the_default_permissions(session, tmp_path):
    # As a file created with open(): 0o666 less the umask (the temporary
    # file is created 0o600)
    umask = os.umask(0)
    os.umask(umask)
    path = tmp_path / "s.acs2"
    session.save(str(path))
    assert stat.S_IMODE(path.stat().st_mode) == 0o666 & ~umask


def test_overwritten_archive_keeps_its_permissions(session, tmp_path):
    path = tmp_path / "s.acs2"
    session.save(str(path))
    path.chmod(0o640)
    session.save(str(path))
    assert stat.S_IMODE(path.stat().st_mode) == 0o640


@pytest.mark.skipif(not hasattr(os, "symlink"), reason="no symbolic links")
def test_save_through_a_symbolic_link(session, tmp_path):
    # The archive is written where the link points, the link is kept
    target = tmp_path / "real.acs2"
    session.save(str(target))
    link = tmp_path / "link.acs2"
    link.symlink_to(target)
    session.save(str(link), initial_session=session)
    assert link.is_symlink()
    assert "link_spec_orig.fits" in _names(target)


def test_name_without_extension(session, tmp_path):
    session.save(str(tmp_path / "s"))
    assert os.listdir(tmp_path) == ["s.acs2"]
