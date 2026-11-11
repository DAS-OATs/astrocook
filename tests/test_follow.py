"""The shared part of a session followed by the GUI (``astrocook/follow.py``,
phase 6-ter PR 5-bis): the marker, the archive, whether the GUI can run and
the command that launches it."""
import json
import os
import subprocess
import sys

import pytest

from astrocook import __version__, follow


def test_the_marker_is_next_to_the_archive(tmp_path):
    archive = tmp_path / "s.acs2"
    assert follow.marker_path(archive) == tmp_path / "s.acs2.follow.json"
    assert follow.read_marker(archive) is None
    path = follow.write_marker(archive, "s-1")
    assert json.loads(path.read_text()) == {"session_id": "s-1", "pid": os.getpid(),
                                            "version": __version__}
    marker = follow.read_marker(archive)
    assert marker == follow.FollowMarker("s-1", os.getpid(), __version__)
    assert marker.writer_alive()
    # Written atomically: no temporary file is left
    assert sorted(p.name for p in tmp_path.iterdir()) == ["s.acs2.follow.json"]


@pytest.mark.parametrize("text", ["not json", "{}", '{"session_id": "s", "pid": "x"}'])
def test_a_marker_that_cannot_be_read_is_none(tmp_path, text):
    archive = tmp_path / "s.acs2"
    follow.marker_path(archive).write_text(text)
    assert follow.read_marker(archive) is None


@pytest.mark.skipif(os.name != "posix", reason="checked on POSIX only")
def test_a_dead_writer():
    process = subprocess.Popen([sys.executable, "-c", "pass"])
    process.wait()
    assert not follow.FollowMarker("s", process.pid, "").writer_alive()


def test_remove(tmp_path):
    archive = tmp_path / "s.acs2"
    archive.write_bytes(b"x")
    follow.write_marker(archive, "s")
    follow.remove(archive)
    assert list(tmp_path.iterdir()) == []
    follow.remove(archive)  # nothing left: no error


def test_the_archive_is_the_session_as_saved(tmp_path):
    from astrocook.session.session import load_session_from_file
    from astrocook.session.history import SessionHistory
    from astrocook.core.structures import HistoryLogV2
    from test_gui_smoke import SESSION_FILE
    session = load_session_from_file(str(SESSION_FILE), "s", None, "auto")
    history = SessionHistory(session, HistoryLogV2())
    history.log_manager.add_note("A note.")
    follow.write_archive(history, tmp_path / "s.acs2")
    saved = load_session_from_file(str(tmp_path / "s.acs2"), "s", None, "auto")
    assert len(saved.systs.components) == len(session.systs.components)
    assert [n.text for n in saved.log_manager.notes] == ["A note."]


# --- Whether the GUI can run here, and how it is launched ------------------------------------

def test_gui_unavailable_without_pyside(monkeypatch):
    real = follow.importlib.util.find_spec
    monkeypatch.setattr(follow.importlib.util, "find_spec",
                        lambda name: None if name == "PySide6" else real(name))
    assert "PySide6 is not installed" in follow.gui_unavailable_reason()


@pytest.mark.parametrize("platform, env, available", [
    ("linux", {}, False),
    ("linux", {"DISPLAY": ":0"}, True),
    ("linux", {"WAYLAND_DISPLAY": "wayland-0"}, True),
    ("darwin", {}, True),
    ("win32", {}, True),
])
def test_gui_unavailable_without_a_display(monkeypatch, platform, env, available):
    monkeypatch.setattr(follow.sys, "platform", platform)
    for name in ("DISPLAY", "WAYLAND_DISPLAY"):
        monkeypatch.delenv(name, raising=False)
    for name, value in env.items():
        monkeypatch.setenv(name, value)
    reason = follow.gui_unavailable_reason()
    assert (reason is None) == available
    if not available:
        assert "no display" in reason


def test_the_command_on_macos(monkeypatch, tmp_path):
    monkeypatch.setattr(follow.sys, "platform", "darwin")
    app = tmp_path / "Astrocook.app"
    app.mkdir()
    monkeypatch.setattr(follow, "MACOS_APPS", (tmp_path / "missing.app", app))
    assert follow.gui_command("/x/s.acs2") == ["open", "-a", "Astrocook", "/x/s.acs2"]
    # Without the application: the command, as elsewhere
    monkeypatch.setattr(follow, "MACOS_APPS", (tmp_path / "missing.app",))
    monkeypatch.setattr(follow.shutil, "which", lambda name: "/bin/astrocook")
    assert follow.gui_command("/x/s.acs2") == ["/bin/astrocook", "/x/s.acs2"]


def test_the_command_elsewhere(monkeypatch):
    monkeypatch.setattr(follow.sys, "platform", "linux")
    monkeypatch.setattr(follow.shutil, "which", lambda name: "/usr/bin/astrocook")
    assert follow.gui_command("/x/s.acs2") == ["/usr/bin/astrocook", "/x/s.acs2"]
    monkeypatch.setattr(follow.shutil, "which", lambda name: None)
    assert follow.gui_command("/x/s.acs2") == [sys.executable, "-m", "astrocook", "/x/s.acs2"]


def test_launch_is_detached_and_silent(monkeypatch):
    launched = []
    monkeypatch.setattr(follow.subprocess, "Popen",
                        lambda command, **kwargs: launched.append((command, kwargs)))
    monkeypatch.setattr(follow, "gui_command", lambda path: ["astrocook", str(path)])
    assert follow.launch_gui("/x/s.acs2") == ["astrocook", "/x/s.acs2"]
    ((command, kwargs),) = launched
    # Nothing reaches the standard output of the server (its channel to the
    # client), and the GUI outlives it
    for stream in ("stdin", "stdout", "stderr"):
        assert kwargs[stream] is subprocess.DEVNULL
    if os.name == "posix":
        assert kwargs["start_new_session"]
