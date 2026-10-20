"""Archives are unpacked without writing outside the temporary folder (F35).

A `.acs2` or `.acs` archive is a tar.gz; a member named `../file` (or with
an absolute path, or a link pointing outside) would be written outside the
temporary folder by a plain `extractall`. The loaders use the "data"
extraction filter, which refuses such members.
"""
import io
import tarfile

from astrocook.core.session import load_session_from_file
from astrocook.io.v1_archive import V1ArchiveManager


def _malicious_archive(path):
    data = b"written outside the temporary folder"
    with tarfile.open(path, "w:gz") as tar:
        info = tarfile.TarInfo("../f35_escape.txt")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))


def _unpack_dir(tmp_path, monkeypatch):
    # The archive is unpacked in tmp_path/unpack, so an escaping member
    # would land in tmp_path
    unpack = tmp_path / "unpack"
    unpack.mkdir()
    monkeypatch.setattr("tempfile.mkdtemp", lambda *a, **k: str(unpack))


def test_acs2_member_outside_the_folder_is_refused(tmp_path, monkeypatch):
    archive = tmp_path / "evil.acs2"
    _malicious_archive(archive)
    _unpack_dir(tmp_path, monkeypatch)

    result = load_session_from_file(str(archive), "evil", None, "auto")

    assert result == 0
    assert not (tmp_path / "f35_escape.txt").exists()


def test_acs_member_outside_the_folder_is_refused(tmp_path, monkeypatch):
    archive = tmp_path / "evil.acs"
    _malicious_archive(archive)
    _unpack_dir(tmp_path, monkeypatch)

    assert V1ArchiveManager(str(archive)).unpack() is None
    assert not (tmp_path / "f35_escape.txt").exists()


def test_regular_acs_archive_is_unpacked(tmp_path, monkeypatch):
    # No V1 sample archive yet: regular files in a subfolder, as V1 writes
    archive = tmp_path / "good.acs"
    data = b"x"
    with tarfile.open(archive, "w:gz") as tar:
        for name in ("good_spec.fits", "sub/good_systs.fits"):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    _unpack_dir(tmp_path, monkeypatch)

    root = V1ArchiveManager(str(archive)).unpack()

    assert (tmp_path / "unpack" / "sub" / "good_systs.fits").exists()
    assert root == str(tmp_path / "unpack")
