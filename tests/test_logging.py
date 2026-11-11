"""Logging is configured by the application, not by the package (phase 3).

Like other libraries, Astrocook leaves the logging configuration to the
script that uses it (by default Python shows warnings and errors only);
the GUI launcher configures it. Before phase 3 the import of
``astrocook.legacy`` called ``logging.basicConfig`` as a side effect.
"""
import logging
import subprocess
import sys
from pathlib import Path

import pytest

import astrocook._run_astrocook as launcher


def _run(code):
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def test_import_does_not_configure_logging():
    code = ("import logging, astrocook, astrocook.session.session, "
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


# --- Where the messages of the package go (phase 6, PR 3) --------------------------
# Characterization tests written before the move to loggers per module
# (#566): the GUI and the MCP server attach their handlers to the root
# logger, the launcher configures the root logger with the name of the
# module in the format, a script that configures nothing still sees the
# warnings. A recipe that logs a warning: rebin with norm (F37).

TESTS_DIR = str(Path(__file__).parent)
NORM_WARNING = "'norm' is ignored"

_REBIN = ("from test_recipes_edit_flux import _spec\n"
          "from astrocook.session.session import SessionV2\n"
          "SessionV2('s', None, spec=_spec()).flux.rebin(norm='True')\n")


def _run_stderr(code):
    code = f"import sys; sys.path.insert(0, {TESTS_DIR!r})\n" + code
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    return result.stderr


class _Records(logging.Handler):
    def __init__(self):
        super().__init__(level=logging.DEBUG)
        self.records = []

    def emit(self, record):
        self.records.append(record)


@pytest.fixture
def root_records():
    """Records reaching a handler on the root logger, as in the GUI."""
    root = logging.getLogger()
    handler, level = _Records(), root.level
    root.addHandler(handler)
    root.setLevel(logging.DEBUG)
    try:
        yield handler.records
    finally:
        root.removeHandler(handler)
        root.setLevel(level)


def test_recipe_messages_reach_a_root_handler(root_records):
    from test_recipes_edit_flux import _spec
    from astrocook.session.session import SessionV2
    SessionV2("s", None, spec=_spec()).flux.rebin(norm="True")
    warnings = [r for r in root_records if r.levelno == logging.WARNING]
    assert [r.module for r in warnings] == ["flux"]
    assert NORM_WARNING in warnings[0].getMessage()
    # The debug lines of the core reach it too
    assert {"session", "spectrum_combine"} <= {r.module for r in root_records
                                               if r.levelno == logging.DEBUG}


def test_gui_warning_capture_sees_recipe_warnings():
    from astrocook.gui.recipe_progress import WarningCaptureHandler
    from test_recipes_edit_flux import _spec
    from astrocook.session.session import SessionV2
    handler = WarningCaptureHandler()
    logging.getLogger().addHandler(handler)
    try:
        SessionV2("s", None, spec=_spec()).flux.rebin(norm="True")
    finally:
        logging.getLogger().removeHandler(handler)
    assert len(handler.warnings) == 1 and NORM_WARNING in handler.warnings[0]


def test_launcher_format_names_the_module():
    stderr = _run_stderr(
        "import logging\n"
        "logging.basicConfig(level=logging.INFO, "
        "format='[%(levelname)s] %(module)s: %(message)s')\n" + _REBIN)
    assert f"[WARNING] flux: The parameter {NORM_WARNING}" in stderr


def test_script_without_logging_configuration_sees_warnings():
    stderr = _run_stderr(_REBIN)
    assert NORM_WARNING in stderr
    assert "Rebin: Interpolating" not in stderr  # DEBUG: not shown


def test_script_can_configure_logging_after_a_recipe():
    # F96: the module-level functions (logging.warning(...)) call
    # logging.basicConfig() when the root logger has no handler, so a
    # recipe run before the script configured logging made the script's
    # own basicConfig a no-op: its level and format were ignored. The
    # loggers per module do not configure the root logger.
    stderr = _run_stderr(
        _REBIN +
        "import logging\n"
        "logging.basicConfig(level=logging.INFO, "
        "format='[%(levelname)s] %(module)s: %(message)s')\n"
        + _REBIN)
    assert f"[WARNING] flux: The parameter {NORM_WARNING}" in stderr


def test_records_name_their_module(root_records):
    from test_recipes_edit_flux import _spec
    from astrocook.session.session import SessionV2
    SessionV2("s", None, spec=_spec()).flux.rebin(norm="True")
    names = {r.name for r in root_records if r.name.startswith("astrocook")}
    assert "astrocook.recipes.flux" in names
    assert "astrocook.core.spectrum_combine" in names
    assert "root" not in {r.name for r in root_records}


def _direct_root_calls():
    """Calls of the module-level logging functions in the package."""
    import ast
    import astrocook
    funcs = {"debug", "info", "warning", "error", "exception", "critical",
             "log", "warn", "fatal"}
    package = Path(astrocook.__file__).parent
    found = []
    for path in sorted(package.rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                    and node.func.attr in funcs
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id == "logging"):
                found.append(f"{path.relative_to(package)}:{node.lineno}")
    return found


def test_package_logs_through_module_loggers():
    # Each module logs with logger = logging.getLogger(__name__) (#566)
    assert _direct_root_calls() == []
