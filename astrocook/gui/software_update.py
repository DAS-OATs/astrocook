"""Software Update: check for a newer release on GitHub and offer the update.

Phase 7-quater, PR U1 (level 1 only: no automatic replacement of the running
application). The pieces:

* pure functions (version comparison, parsing of the release, recognition of
  the installation, SHA-256 check): no Qt, tested without network;
* two small workers (``QThread``) that bring their result back through Qt
  signals: the check of the latest release and the download of the ``.dmg``;
* ``SoftwareUpdateDialog``, whose offer depends on the installation: the
  frozen macOS app downloads the ``.dmg``, an editable/git install is told to
  ``git pull``, any other (pip) install is shown the command to copy.

The request carries a generic User-Agent and nothing of the user.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import urllib.request
from http.client import HTTPException
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Literal, Optional

from packaging.version import InvalidVersion, Version
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QLabel, QPlainTextEdit, QProgressBar,
    QPushButton, QVBoxLayout)

RELEASES_API = "https://api.github.com/repos/DAS-OATs/astrocook/releases/latest"
REPO_URL = "https://github.com/DAS-OATs/astrocook"
#: Fixed name of the asset of the macOS app (``docs/installation.md``)
DMG_ASSET = "Astrocook.dmg"
USER_AGENT = "Astrocook-update-check"
TIMEOUT = 10  # seconds

#: What a request or the reading of its answer can raise (JSON errors are ValueErrors)
_NETWORK_ERRORS = (OSError, ValueError, HTTPException, AttributeError)

InstallKind = Literal["frozen_macos", "editable", "pip"]


# --- Pure functions ---------------------------------------------------------

def parse_version(tag: str) -> Optional[Version]:
    """The version of a tag (a leading ``v`` is accepted), None if it is not one."""
    try:
        return Version(tag.strip().lstrip("vV"))
    except InvalidVersion:
        return None


def is_newer(tag: str, current: str) -> bool:
    """True if the release ``tag`` is newer than the running ``current`` version.

    Compared with ``packaging.version``, never as strings (``2.10.0`` > ``2.8.0``;
    a pre-release is older than its release). An unparsable tag is not newer.
    """
    latest, running = parse_version(tag), parse_version(current)
    return latest is not None and running is not None and latest > running


@dataclass(frozen=True)
class ReleaseInfo:
    tag: str
    notes: str
    html_url: str
    dmg_url: Optional[str] = None
    dmg_digest: Optional[str] = None  # "sha256:<hex>" as the API gives it
    dmg_size: Optional[int] = None


def parse_release(data: dict) -> ReleaseInfo:
    """Build a ``ReleaseInfo`` from the JSON of ``releases/latest``."""
    tag = str(data.get("tag_name", ""))
    dmg = next((a for a in data.get("assets") or [] if a.get("name") == DMG_ASSET), None)
    return ReleaseInfo(
        tag=tag,
        notes=str(data.get("body") or ""),
        html_url=str(data.get("html_url") or f"{REPO_URL}/releases/tag/{tag}"),
        dmg_url=dmg.get("browser_download_url") if dmg else None,
        dmg_digest=dmg.get("digest") if dmg else None,
        dmg_size=dmg.get("size") if dmg else None)


def fetch_latest_release(urlopen: Callable = urllib.request.urlopen) -> ReleaseInfo:
    """Ask GitHub for the latest release (blocking: call it from a thread)."""
    req = urllib.request.Request(RELEASES_API, headers={
        "User-Agent": USER_AGENT, "Accept": "application/vnd.github+json"})
    with urlopen(req, timeout=TIMEOUT) as resp:
        return parse_release(json.loads(resp.read().decode("utf-8")))


def install_kind(frozen: Optional[bool] = None, platform: Optional[str] = None,
                 package_dir: Optional[Path] = None) -> InstallKind:
    """How Astrocook is installed: the frozen macOS app, an editable install or
    git checkout, or any other (pip) install. The arguments are for the tests."""
    frozen = getattr(sys, "frozen", False) if frozen is None else frozen
    platform = sys.platform if platform is None else platform
    if frozen and platform == "darwin":
        return "frozen_macos"
    if package_dir is None:
        package_dir = Path(__file__).resolve().parent.parent
    # A git checkout (editable or not) has the .git folder above the package;
    # an editable install from a checkout is the same thing.
    if (package_dir.parent / ".git").exists():
        return "editable"
    return "pip"


def pip_command(tag: str) -> str:
    """The command to update by pip (no PyPI yet: from the tag on GitHub)."""
    return f"pip install -U git+{REPO_URL}@{tag}"


def sha256_of(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def verify_digest(path: str, digest: Optional[str]) -> Optional[bool]:
    """True/False if the file matches the API ``digest`` (``sha256:<hex>``);
    None if there is no (usable) digest: the file is then not verified."""
    if not digest or not digest.lower().startswith("sha256:"):
        return None
    return sha256_of(path) == digest.split(":", 1)[1].strip().lower()


def download_file(url: str, dest: str, progress: Callable[[int, int], None] = lambda a, b: None,
                  cancelled: Callable[[], bool] = lambda: False,
                  urlopen: Callable = urllib.request.urlopen, chunk: int = 1 << 16) -> bool:
    """Download ``url`` to ``dest`` (blocking). False if cancelled (the partial
    file is removed)."""
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    done = 0
    with urlopen(req, timeout=TIMEOUT) as resp, open(dest, "wb") as out:
        total = int(resp.headers.get("Content-Length") or 0)
        while block := resp.read(chunk):
            if cancelled():
                out.close()
                os.remove(dest)
                return False
            out.write(block)
            done += len(block)
            progress(done, total)
    return True


# --- Workers ----------------------------------------------------------------

class _CheckWorker(QThread):
    found = Signal(object)   # ReleaseInfo
    failed = Signal(str)

    def __init__(self, parent=None, urlopen: Callable = urllib.request.urlopen):
        super().__init__(parent)
        self._urlopen = urlopen

    def run(self):
        try:
            self.found.emit(fetch_latest_release(self._urlopen))
        except _NETWORK_ERRORS as exc:  # network, JSON, HTTP: all end as a message
            self.failed.emit(str(exc))


class _DownloadWorker(QThread):
    progress = Signal(int, int)
    finished_ok = Signal(str, object)  # path, verified (True/False/None)
    failed = Signal(str)

    def __init__(self, url: str, dest: str, digest: Optional[str], parent=None,
                 urlopen: Callable = urllib.request.urlopen):
        super().__init__(parent)
        self._url, self._dest, self._digest, self._urlopen = url, dest, digest, urlopen
        self._cancel = False

    def cancel(self):
        self._cancel = True

    def run(self):
        try:
            ok = download_file(self._url, self._dest, self.progress.emit,
                               lambda: self._cancel, self._urlopen)
            if ok:
                self.finished_ok.emit(self._dest, verify_digest(self._dest, self._digest))
        except _NETWORK_ERRORS as exc:
            self.failed.emit(str(exc))


# --- Dialog -----------------------------------------------------------------

class SoftwareUpdateDialog(QDialog):
    """Checks for a newer release when shown and offers what the installation
    allows. ``urlopen``, ``kind`` and ``opener`` are injectable for the tests
    and the screenshot (no network)."""

    def __init__(self, parent=None, current: Optional[str] = None,
                 kind: Optional[InstallKind] = None,
                 urlopen: Callable = urllib.request.urlopen,
                 opener: Optional[Callable[[str], None]] = None,
                 download_dir: Optional[str] = None,
                 release: Optional[ReleaseInfo] = None):
        super().__init__(parent)
        from astrocook import __version__
        self.current = current or __version__
        self.kind: InstallKind = kind or install_kind()
        self._urlopen = urlopen
        self._opener = opener or self._open_with_system
        self._download_dir = download_dir
        self.release: Optional[ReleaseInfo] = None
        self._worker: Optional[QThread] = None
        self.setWindowTitle("Software Update")
        self.setMinimumWidth(480)

        layout = QVBoxLayout(self)
        self.status = QLabel("Checking for updates...")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.notes = QPlainTextEdit()
        self.notes.setReadOnly(True)
        self.notes.setVisible(False)
        layout.addWidget(self.notes)
        self.command = QLabel()
        self.command.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.command.setVisible(False)
        layout.addWidget(self.command)
        self.progress = QProgressBar()
        self.progress.setVisible(False)
        layout.addWidget(self.progress)
        self.action_button = QPushButton()
        self.action_button.setVisible(False)
        layout.addWidget(self.action_button)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        if release is not None:  # already known (screenshot)
            self._on_found(release)
        else:
            worker = _CheckWorker(self, self._urlopen)
            worker.found.connect(self._on_found)
            worker.failed.connect(self._on_failed)
            self._worker = worker
            worker.start()

    # -- check
    def _on_failed(self, message: str):
        self.status.setText("Could not check for updates (no network connection, or "
                            f"GitHub is not reachable).\n\n{message}")

    def _on_found(self, release: ReleaseInfo):
        self.release = release
        if not is_newer(release.tag, self.current):
            self.status.setText(f"Astrocook {self.current} is up to date "
                                f"(latest release: {release.tag}).")
            return
        self.status.setText(f"A new version is available: <b>{release.tag}</b> "
                            f"(you have {self.current}).")
        self.notes.setPlainText(release.notes)
        self.notes.setVisible(True)
        if self.kind == "frozen_macos":
            if release.dmg_url:
                self.action_button.setText("Download")
                self.action_button.clicked.connect(self._start_download)
                self.action_button.setVisible(True)
            else:
                self.status.setText(self.status.text() + f"<br>The release has no {DMG_ASSET}; "
                                    f"see {release.html_url}")
        elif self.kind == "editable":
            self.status.setText(self.status.text() + "<br>This is a git checkout or an editable "
                                "install: update it with <code>git pull</code>. Nothing is changed "
                                "automatically.")
        else:
            cmd = pip_command(release.tag)
            self.status.setText(self.status.text() + "<br>To update, run in the same environment:")
            self.command.setText(f"<code>{cmd}</code>")
            self.command.setVisible(True)
            self.action_button.setText("Copy")
            self.action_button.clicked.connect(lambda: QGuiApplication.clipboard().setText(cmd))
            self.action_button.setVisible(True)

    # -- download (frozen macOS app)
    def _start_download(self):
        rel = self.release
        assert rel is not None and rel.dmg_url
        folder = self._download_dir or str(Path.home() / "Downloads")
        os.makedirs(folder, exist_ok=True)
        dest = os.path.join(folder, DMG_ASSET)
        self.action_button.setEnabled(False)
        self.progress.setVisible(True)
        worker = _DownloadWorker(rel.dmg_url, dest, rel.dmg_digest, self, self._urlopen)
        worker.progress.connect(self._on_progress)
        worker.finished_ok.connect(self._on_downloaded)
        worker.failed.connect(self._on_download_failed)
        self._worker = worker
        worker.start()

    def _on_progress(self, done: int, total: int):
        self.progress.setRange(0, total or 0)
        self.progress.setValue(done)

    def _on_download_failed(self, message: str):
        self.progress.setVisible(False)
        self.action_button.setEnabled(True)
        self.status.setText(f"Download failed: {message}")

    def _on_downloaded(self, path: str, verified):
        self.progress.setVisible(False)
        if verified is False:
            os.remove(path)
            self.action_button.setEnabled(True)
            self.status.setText("The downloaded file does not match its SHA-256 checksum: it was "
                                "deleted. Try again.")
            return
        note = ("Checksum verified." if verified else
                "The checksum could not be verified (the release gives none).")
        self.status.setText(f"Downloaded {path}. {note}<br>Close Astrocook, then drag the new "
                            "app from the disk image into Applications.")
        self.action_button.setVisible(False)
        self._opener(path)

    @staticmethod
    def _open_with_system(path: str):
        subprocess.Popen(["open", path])

    def done(self, result):
        w = self._worker
        if w is not None and w.isRunning():
            if isinstance(w, _DownloadWorker):
                w.cancel()
            w.wait(3000)
        super().done(result)
