"""Software Update (phase 7-quater, PR U1): version comparison, release parsing,
recognition of the installation, checksum, and the dialog. No network: the
requests go through a fake ``urlopen``."""
import hashlib
import io
import json
from pathlib import Path

import pytest

from astrocook.gui import software_update as su
from astrocook.gui.software_update import (
    ReleaseInfo, SoftwareUpdateDialog, is_newer, parse_release, pip_command)


@pytest.mark.parametrize("tag,current,expected", [
    ("2.8.0", "2.7.0", True),
    ("2.10.0", "2.8.0", True),          # not a string comparison
    ("2.8.0", "2.10.0", False),
    ("v2.8.0", "2.7.1", True),          # a tag with v
    ("2.8.0", "2.8.0", False),
    ("2.8.0", "2.8.0rc1", True),        # a pre-release is older than its release
    ("2.8.0rc1", "2.8.0", False),
    ("2.8.0rc2", "2.8.0rc1", True),
    ("nightly", "2.8.0", False),        # not a version
    ("2.8.0", "dev", False),
])
def test_is_newer(tag, current, expected):
    assert is_newer(tag, current) is expected


RELEASE_JSON = {
    "tag_name": "v2.9.0", "body": "Notes **here**", "html_url": "https://example/r",
    "assets": [{"name": "other.zip", "browser_download_url": "https://example/o"},
               {"name": "Astrocook.dmg", "browser_download_url": "https://example/a.dmg",
                "digest": "sha256:abc", "size": 5}]}


def test_parse_release():
    r = parse_release(RELEASE_JSON)
    assert (r.tag, r.notes, r.dmg_url, r.dmg_digest, r.dmg_size) == (
        "v2.9.0", "Notes **here**", "https://example/a.dmg", "sha256:abc", 5)
    bare = parse_release({"tag_name": "2.9.0", "body": None})
    assert bare.dmg_url is None and bare.dmg_digest is None and bare.notes == ""


class FakeResponse(io.BytesIO):
    def __init__(self, data: bytes):
        super().__init__(data)
        self.headers = {"Content-Length": str(len(data))}

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()


def fake_urlopen(payloads):
    calls = []

    def urlopen(req, timeout=None):
        calls.append(req)
        data = payloads[req.full_url]
        if isinstance(data, Exception):
            raise data
        return FakeResponse(data)
    urlopen.calls = calls
    return urlopen


def test_fetch_sends_only_a_generic_user_agent():
    uo = fake_urlopen({su.RELEASES_API: json.dumps(RELEASE_JSON).encode()})
    assert su.fetch_latest_release(uo).tag == "v2.9.0"
    headers = {k.lower(): v for k, v in uo.calls[0].header_items()}
    assert headers["user-agent"] == su.USER_AGENT and "authorization" not in headers


def test_install_kind(tmp_path):
    assert su.install_kind(frozen=True, platform="darwin") == "frozen_macos"
    pkg = tmp_path / "astrocook"
    pkg.mkdir()
    assert su.install_kind(frozen=False, platform="darwin", package_dir=pkg) == "pip"
    (tmp_path / ".git").mkdir()
    assert su.install_kind(frozen=False, platform="linux", package_dir=pkg) == "editable"
    # frozen on another platform is not the macOS app
    assert su.install_kind(frozen=True, platform="win32", package_dir=pkg) == "editable"


def test_pip_command():
    assert pip_command("v2.9.0") == (
        "pip install -U git+https://github.com/DAS-OATs/astrocook@v2.9.0")


def test_verify_digest(tmp_path):
    f = tmp_path / "x"
    f.write_bytes(b"data")
    good = "sha256:" + hashlib.sha256(b"data").hexdigest().upper()
    assert su.verify_digest(str(f), good) is True
    assert su.verify_digest(str(f), "sha256:" + "0" * 64) is False
    assert su.verify_digest(str(f), None) is None
    assert su.verify_digest(str(f), "md5:abc") is None


def _release(digest=None):
    return ReleaseInfo("v2.9.0", "Notes", "https://example/r", "https://example/a.dmg", digest)


def test_dialog_up_to_date(qtbot):
    dlg = SoftwareUpdateDialog(current="2.9.0", kind="pip", release=_release())
    qtbot.addWidget(dlg)
    assert "up to date" in dlg.status.text() and not dlg.action_button.isVisibleTo(dlg)


def test_dialog_pip_shows_command_and_copies(qtbot):
    from PySide6.QtGui import QGuiApplication
    dlg = SoftwareUpdateDialog(current="2.8.0", kind="pip", release=_release())
    qtbot.addWidget(dlg)
    assert pip_command("v2.9.0") in dlg.command.text()
    assert dlg.notes.toPlainText() == "Notes"
    dlg.action_button.click()
    assert QGuiApplication.clipboard().text() == pip_command("v2.9.0")


def test_dialog_editable_only_warns(qtbot):
    dlg = SoftwareUpdateDialog(current="2.8.0", kind="editable", release=_release())
    qtbot.addWidget(dlg)
    assert "git pull" in dlg.status.text() and not dlg.action_button.isVisibleTo(dlg)


def test_dialog_check_failure_is_a_message(qtbot):
    uo = fake_urlopen({su.RELEASES_API: OSError("offline")})
    dlg = SoftwareUpdateDialog(current="2.8.0", kind="pip", urlopen=uo)
    qtbot.addWidget(dlg)
    qtbot.waitUntil(lambda: "Could not check" in dlg.status.text(), timeout=5000)


def test_dialog_check_through_thread(qtbot):
    uo = fake_urlopen({su.RELEASES_API: json.dumps(RELEASE_JSON).encode()})
    dlg = SoftwareUpdateDialog(current="2.8.0", kind="pip", urlopen=uo)
    qtbot.addWidget(dlg)
    qtbot.waitUntil(lambda: dlg.release is not None, timeout=5000)
    assert "v2.9.0" in dlg.status.text()


@pytest.mark.parametrize("digest_ok,expect", [(True, "Checksum verified"),
                                              (None, "could not be verified")])
def test_dialog_download_opens_dmg(qtbot, tmp_path, digest_ok, expect):
    data = b"dmg-bytes" * 1000
    digest = "sha256:" + hashlib.sha256(data).hexdigest() if digest_ok else None
    opened = []
    uo = fake_urlopen({"https://example/a.dmg": data})
    dlg = SoftwareUpdateDialog(current="2.8.0", kind="frozen_macos", urlopen=uo,
                               opener=opened.append, download_dir=str(tmp_path),
                               release=_release(digest))
    qtbot.addWidget(dlg)
    dlg.action_button.click()
    qtbot.waitUntil(lambda: bool(opened), timeout=5000)
    assert Path(opened[0]).read_bytes() == data
    assert expect in dlg.status.text()


def test_dialog_download_bad_checksum_deletes_file(qtbot, tmp_path):
    uo = fake_urlopen({"https://example/a.dmg": b"abc"})
    opened = []
    dlg = SoftwareUpdateDialog(current="2.8.0", kind="frozen_macos", urlopen=uo,
                               opener=opened.append, download_dir=str(tmp_path),
                               release=_release("sha256:" + "0" * 64))
    qtbot.addWidget(dlg)
    dlg.action_button.click()
    qtbot.waitUntil(lambda: "does not match" in dlg.status.text(), timeout=5000)
    assert not opened and not (tmp_path / "Astrocook.dmg").exists()


def test_dialog_development_install_compares_nothing(qtbot, monkeypatch):
    """A checkout has no meaningful version: no network, no comparison, a
    neutral message with the commit when it can be read."""
    monkeypatch.setattr(su, "git_commit", lambda *a, **k: "abc1234")
    uo = fake_urlopen({})  # any call would raise KeyError
    dlg = SoftwareUpdateDialog(current="2.7.0", kind="editable", urlopen=uo)
    qtbot.addWidget(dlg)
    qtbot.wait(100)
    assert uo.calls == [] and dlg._worker is None and dlg.release is None
    text = dlg.status.text()
    assert "development installation" in text and "version 2.7.0, commit abc1234" in text
    assert "new version" not in text and not dlg.action_button.isVisibleTo(dlg)


def test_dialog_development_install_without_a_readable_commit(qtbot, monkeypatch):
    monkeypatch.setattr(su, "git_commit", lambda *a, **k: None)
    dlg = SoftwareUpdateDialog(current="2.7.0", kind="editable", urlopen=fake_urlopen({}))
    qtbot.addWidget(dlg)
    assert "version 2.7.0)" in dlg.status.text() and "commit" not in dlg.status.text()


def test_git_commit_is_none_outside_a_repository(tmp_path):
    assert su.git_commit(tmp_path / "pkg") is None
