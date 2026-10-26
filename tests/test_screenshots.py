"""Smoke test of the script that rebuilds the screenshots of the tutorials.

The scenes drive the GUI through its internals (methods of the main window,
menus, dialogs): a change of the GUI can break them. Two quick scenes run on
a small test archive in place of the spectra of the tutorials, so that such
a change shows up in the pull request that makes it. The images themselves
are not checked.
"""
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "docs" / "screenshots" / "make_screenshots.py"
ARCHIVE = Path(__file__).parent / "data" / "J0103-1305_metallines_lya_new_3.acs2"
SCENES = ["getting_started_main_view", "getting_started_session_list"]


def test_scenes_run(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    for name in ("My_QSO", "My_other_QSO"):
        shutil.copyfile(ARCHIVE, data / f"{name}.acs2")
    out = tmp_path / "out"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--data", str(data), "--out", str(out),
         "--only", *SCENES],
        capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, result.stderr[-3000:]
    for scene in SCENES:
        with Image.open(out / f"{scene}.png") as image:
            # Twice the logical size of the window, or a crop of it
            assert image.width > 1000 and image.height > 700
