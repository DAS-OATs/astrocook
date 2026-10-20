"""Logging is configured by the application, not by the package (phase 3).

Like other libraries, Astrocook leaves the logging configuration to the
script that uses it (by default Python shows warnings and errors only);
the GUI launcher configures it. Before phase 3 the import of
``astrocook.legacy`` called ``logging.basicConfig`` as a side effect.
"""
import subprocess
import sys

import astrocook._run_astrocook as launcher


def _run(code):
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def test_import_does_not_configure_logging():
    code = ("import logging, astrocook, astrocook.core.session, "
            "astrocook.recipes.absorbers, astrocook.io.loaders\n"
            "print(len(logging.getLogger().handlers))")
    assert _run(code) == "0"


def test_script_can_choose_the_level():
    code = ("import logging\n"
            "logging.basicConfig(level=logging.INFO, format='%(levelname)s')\n"
            "import astrocook\n"
            "print(logging.getLevelName(logging.getLogger().level))")
    assert _run(code) == "INFO"


def test_gui_launcher_configures_logging(monkeypatch):
    calls = []
    monkeypatch.setattr(launcher.logging, "basicConfig",
                        lambda **kwargs: calls.append(kwargs))
    monkeypatch.setattr(sys, "argv", ["astrocook", "--help"])
    try:
        launcher.main()
    except SystemExit:
        pass  # --help exits before any window is created
    assert calls and calls[0]["level"] == launcher.logging.INFO
