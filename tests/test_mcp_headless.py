"""The MCP server without Qt (phase 6-ter, PR 6): the headless Docker image
installs the package without PySide6, so nothing the server imports may need
it. ``tests/test_fit_groups_plotting.py`` does the same for
``astrocook.plotting``; here it is the whole server, in a process of its own
(the modules of the suite would hide an import)."""
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("mcp")

CHECK = Path(__file__).parent / "mcp_stdio_check.py"
#: The server in a process where PySide6 cannot be imported (``None`` in
#: ``sys.modules`` makes ``import PySide6`` fail and ``find_spec`` answer
#: None), as in the image
WITHOUT_PYSIDE6 = ("import sys; sys.modules['PySide6'] = None; "
                   "from astrocook.mcp_server.server import main; main()")


def test_building_the_server_imports_no_qt():
    code = (
        "import sys\n"
        "from astrocook.mcp_server.server import build_server\n"
        "server = build_server()\n"
        "qt = sorted(m for m in sys.modules if m.split('.')[0] in\n"
        "            ('PySide6', 'shiboken6', 'qtawesome')\n"
        "            or m == 'astrocook.gui' or m.startswith('astrocook.gui.'))\n"
        "assert not qt, qt\n"
        "assert 'matplotlib.pyplot' not in sys.modules, 'pyplot'\n")
    subprocess.run([sys.executable, "-I", "-c", code], check=True)


def test_the_server_runs_over_stdio_without_pyside6():
    # What the CI job does in the image, here with this interpreter: open a
    # mock, a recipe, save, plot, and gui_open refused because PySide6 is
    # not there (not for the lack of a display)
    done = subprocess.run(
        [sys.executable, str(CHECK), "--", sys.executable, "-c", WITHOUT_PYSIDE6],
        capture_output=True, text=True, timeout=300)
    assert done.returncode == 0, done.stdout + done.stderr
    assert "all steps passed" in done.stdout
    assert "PySide6 is not installed here" in done.stdout
