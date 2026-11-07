"""A session followed by the GUI (phase 6-ter, PR 5-bis; decision 7 of the
phase, option b).

A program that changes a session (the MCP server, a script) saves it to an
archive after every change; the GUI, launched on that archive, reloads it
each time it is replaced, and shows the session read only. This module is
the part shared by the two sides, with no import of Qt or of ``mcp``:
writing the archive, the *marker* that says it is followed, whether the GUI
can run here, and the command that launches it.

The marker is a small JSON file next to the archive,
``<archive>.follow.json``, with the identifier of the session, the process
that writes the archive and the version of Astrocook. The GUI looks for it
whenever it opens an archive, whichever way the file arrives (the command
line, the macOS ``open`` event, drag and drop): a flag of the command line
would not reach an application already running on macOS. When the writer
stops following the session it removes the marker, and the GUI says that the
session is no longer followed.

The archive is never read half written: :func:`write_archive` saves it as
``session_save`` does, to a temporary file in the same folder that replaces
the archive only when complete (``os.replace``, atomic on the same file
system; ``_write_replacing`` of :mod:`astrocook.io.adapter`, F94). The
marker is written in the same way.
"""
import importlib.util
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union

from astrocook import __version__

logger = logging.getLogger(__name__)

#: Appended to the path of a followed archive to make the path of its marker
MARKER_SUFFIX = ".follow.json"

#: Where the macOS application is looked for, in this order
MACOS_APPS = (Path("/Applications/Astrocook.app"),
              Path.home() / "Applications" / "Astrocook.app")

PathLike = Union[str, os.PathLike]


def marker_path(archive: PathLike) -> Path:
    """The path of the marker of ``archive``."""
    return Path(os.fspath(archive) + MARKER_SUFFIX)


@dataclass(frozen=True)
class FollowMarker:
    """What the marker of a followed archive says."""
    #: The identifier of the session in the program that writes the archive
    session_id: str
    #: The process that writes the archive
    pid: int
    #: The version of Astrocook of that program
    version: str

    def writer_alive(self) -> bool:
        """Whether the process that writes the archive is still running.

        Checked on macOS and Linux only (a signal 0 to the process); on
        Windows it is taken as running. A marker left behind by a writer that
        did not stop cleanly (killed) belongs to a dead process.
        """
        if os.name != "posix":
            return True
        try:
            os.kill(self.pid, 0)
        except ProcessLookupError:
            return False
        except OSError:
            return True  # e.g. running under another user (PermissionError)
        return True


def write_marker(archive: PathLike, session_id: str) -> Path:
    """Write the marker of ``archive`` (atomically) and return its path."""
    path = marker_path(archive)
    text = json.dumps({"session_id": session_id, "pid": os.getpid(),
                       "version": __version__})
    fd, temp = tempfile.mkstemp(dir=path.parent, prefix=".", suffix=".follow.tmp")
    try:
        with os.fdopen(fd, "w") as f:
            f.write(text)
        os.replace(temp, path)
    finally:
        if os.path.exists(temp):
            os.remove(temp)
    return path


def read_marker(archive: PathLike) -> Optional[FollowMarker]:
    """The marker of ``archive``, or None if there is none (or it cannot be
    read: the archive then opens as any other)."""
    path = marker_path(archive)
    try:
        data = json.loads(path.read_text())
        return FollowMarker(session_id=str(data["session_id"]), pid=int(data["pid"]),
                            version=str(data.get("version", "")))
    except FileNotFoundError:
        return None
    except (OSError, ValueError, KeyError, TypeError) as e:
        logger.debug(f"The marker {path} could not be read ({e}): the archive "
                     "is not followed.")
        return None


def write_archive(history, archive: PathLike) -> None:
    """Save the current state of a session history to ``archive``, with its
    first state as the original data and its log, as ``session_save`` does
    (atomically)."""
    history.current_state.save(os.fspath(archive), initial_session=history.states[0])


def remove(archive: PathLike) -> None:
    """Stop following ``archive``: remove its marker, then the archive."""
    for path in (marker_path(archive), Path(archive)):
        try:
            path.unlink()
        except FileNotFoundError:
            pass


def gui_unavailable_reason() -> Optional[str]:
    """Why the GUI cannot be shown here, or None if it can.

    It cannot without PySide6 (the Docker image of the server installs the
    package without it), or on Linux without a display (``DISPLAY`` or
    ``WAYLAND_DISPLAY``).
    """
    if importlib.util.find_spec("PySide6") is None:
        return ("PySide6 is not installed here (e.g. the Docker image of the "
                "server), so the GUI cannot run.")
    if sys.platform.startswith("linux") and not (
            os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")):
        return "There is no display here (DISPLAY and WAYLAND_DISPLAY are not set)."
    return None


def gui_command(archive: PathLike) -> List[str]:
    """The command that opens ``archive`` in the GUI.

    On macOS the application (``open -a Astrocook``), if it is installed in
    ``/Applications`` or ``~/Applications``; else, and on the other systems,
    the ``astrocook`` command, or this Python with ``-m astrocook`` if the
    command is not on the path.
    """
    archive = os.fspath(archive)
    if sys.platform == "darwin" and any(app.is_dir() for app in MACOS_APPS):
        return ["open", "-a", "Astrocook", archive]
    exe = shutil.which("astrocook")
    if exe:
        return [exe, archive]
    return [sys.executable, "-m", "astrocook", archive]


def launch_gui(archive: PathLike) -> List[str]:
    """Launch the GUI on ``archive``, detached from this process, and return
    the command used.

    The GUI reads nothing from this process and writes nothing to it: its
    standard streams go to the null device (the standard output of the MCP
    server is its channel to the client), and it keeps running when this
    process ends.

    Raises
    ------
    OSError
        If the command cannot be started.
    """
    command = gui_command(archive)
    kwargs = dict(stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                  stderr=subprocess.DEVNULL, close_fds=True)
    if os.name == "posix":
        kwargs["start_new_session"] = True
    else:
        kwargs["creationflags"] = (subprocess.DETACHED_PROCESS
                                   | subprocess.CREATE_NEW_PROCESS_GROUP)
    subprocess.Popen(command, **kwargs)
    logger.info(f"GUI launched on {archive}: {' '.join(command)}")
    return command
