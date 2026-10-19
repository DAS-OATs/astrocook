"""The GUI launcher is reachable from an installed package (phase 2).

`python -m astrocook` and the `astrocook` console script both run
`astrocook._run_astrocook.main`; `run_astrocook.py` in the repository root
is a thin wrapper around it (also used by PyInstaller).
"""
import os
import subprocess
import sys
from importlib.metadata import entry_points

import astrocook._run_astrocook as launcher


def test_console_script_points_to_the_launcher():
    scripts = entry_points(group="console_scripts")
    (ep,) = [ep for ep in scripts if ep.name == "astrocook"]
    assert ep.value == "astrocook._run_astrocook:main"
    assert ep.load() is launcher.main


def test_python_m_astrocook_parses_arguments():
    # --help exits in argparse, before any window is created
    result = subprocess.run(
        [sys.executable, "-m", "astrocook", "--help"],
        capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stderr
    assert "session_files" in result.stdout


def test_assets_resolved_from_the_repository_root():
    path = launcher.resource_path("assets/logo_3d_LR.png")
    assert path.endswith("assets/logo_3d_LR.png")
    assert os.path.exists(path)
