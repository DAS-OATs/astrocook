"""The Dockerfile of the headless image of the MCP server (phase 6-ter, PR 6)
stays aligned with ``pyproject.toml``.

The image installs the package with ``pip install --no-deps`` and lists the
other dependencies by hand, to leave out PySide6 without changing the
dependencies of the package. Nothing can build the image in the suite (no
Docker), so these tests read the file: the list must be that of
``pyproject.toml`` (the dependencies and the ``[mcp]`` extra) but for what
only the GUI needs, and what only the GUI needs must really be used only by
it. The CI job ``docker`` builds the image and runs the server in it.
"""
import ast
import fnmatch
import re
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DOCKERFILE = ROOT / "Dockerfile"

#: Required dependencies that the image leaves out, each with its reason
GUI_ONLY = {
    "pyside6": "the GUI (astrocook/gui); useless without a display, ~0.7 GB with Qt",
    "qtawesome": "the icons of the GUI (astrocook/gui)",
}


def _name(spec):
    return re.split(r"[<>=!~ \[;]", spec.strip(), maxsplit=1)[0].lower().replace("_", "-")


def _normalized(spec):
    return spec.replace(" ", "").lower()


def pyproject_dependencies():
    data = tomllib.loads((ROOT / "pyproject.toml").read_text())
    return (data["project"]["dependencies"]
            + data["project"]["optional-dependencies"]["mcp"])


def dockerfile_dependencies():
    """The specifiers between the markers ``deps:begin`` and ``deps:end``."""
    text = DOCKERFILE.read_text()
    block = re.search(r"^# deps:begin\n(.*?)^# deps:end", text, re.S | re.M)
    assert block, "the Dockerfile has no deps:begin / deps:end markers"
    return re.findall(r'"([^"]+)"', block.group(1))


def test_the_dockerfile_lists_the_dependencies_of_the_package():
    docker = {_name(s): _normalized(s) for s in dockerfile_dependencies()}
    expected = {_name(s): _normalized(s) for s in pyproject_dependencies()
                if _name(s) not in GUI_ONLY}
    assert docker == expected, (
        "the dependencies of the Dockerfile differ from pyproject.toml "
        f"(without {sorted(GUI_ONLY)}): only in the Dockerfile "
        f"{sorted(set(docker.items()) - set(expected.items()))}, only in "
        f"pyproject.toml {sorted(set(expected.items()) - set(docker.items()))}")


def test_what_the_image_leaves_out_is_still_a_dependency():
    # A stale entry of GUI_ONLY would hide a dependency that was removed
    # (or renamed) in pyproject.toml
    names = {_name(s) for s in pyproject_dependencies()}
    assert set(GUI_ONLY) <= names


def test_only_the_gui_imports_what_the_image_leaves_out():
    # Without PySide6 (and qtawesome) everything but astrocook/gui must work
    modules = {"pyside6": "PySide6", "qtawesome": "qtawesome",
               "shiboken6": "shiboken6"}
    #: The GUI: its package and the launcher of the ``astrocook`` command
    gui = ("gui", "_run_astrocook.py")
    offenders = []
    for path in (ROOT / "astrocook").rglob("*.py"):
        if path.relative_to(ROOT / "astrocook").parts[0] in gui:
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            names = ([a.name for a in node.names] if isinstance(node, ast.Import)
                     else [node.module or ""] if isinstance(node, ast.ImportFrom)
                     and node.level == 0 else [])
            offenders += [f"{path.relative_to(ROOT)}: {n}" for n in names
                          if n.split(".")[0] in modules.values()]
    assert not offenders, offenders


def test_nothing_imports_scienceplots():
    # F110: it was a dependency that no code imported, and was removed from
    # pyproject.toml (and from the macOS build): an import of it would fail
    # in a new installation
    for path in (ROOT / "astrocook").rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            names = ([a.name for a in node.names] if isinstance(node, ast.Import)
                     else [node.module or ""] if isinstance(node, ast.ImportFrom)
                     else [])
            assert not [n for n in names if n.split(".")[0] == "scienceplots"], path
    assert "scienceplots" not in (ROOT / "Astrocook.spec").read_text()


def test_the_dockerfile_installs_the_package_without_the_gui():
    text = DOCKERFILE.read_text()
    assert re.search(r"^RUN pip install --no-deps \.", text, re.M)
    # PySide6 and qtawesome only appear in comments
    code = "\n".join(line for line in text.splitlines()
                     if not line.lstrip().startswith("#"))
    assert "pyside6" not in code.lower() and "qtawesome" not in code.lower()
    # stdio, the command a client starts, and not as root
    assert re.search(r'^ENTRYPOINT \["astrocook-mcp"\]', text, re.M)
    assert re.search(r"^USER (?!root)\w+", text, re.M)


def test_the_python_of_the_image_is_supported_by_the_package():
    data = tomllib.loads((ROOT / "pyproject.toml").read_text())
    minimum = tuple(int(n) for n in
                    re.search(r">=\s*([\d.]+)", data["project"]["requires-python"])
                    .group(1).split("."))
    tag = re.search(r"^FROM python:(\d+)\.(\d+)", DOCKERFILE.read_text(), re.M)
    assert tag and (int(tag.group(1)), int(tag.group(2))) >= minimum


def test_the_build_context_has_what_the_dockerfile_copies():
    copied = [name for line in DOCKERFILE.read_text().splitlines()
              if line.startswith("COPY ")
              for name in line.split()[1:-1]]
    assert copied
    for name in copied:
        assert (ROOT / name).exists(), name
    # ... and the .dockerignore does not leave them out
    ignored = [line.strip() for line in (ROOT / ".dockerignore").read_text().splitlines()
               if line.strip() and not line.startswith("#")]
    for name in copied:
        assert not any(fnmatch.fnmatch(name, pat) for pat in ignored), name


@pytest.mark.parametrize("pattern", ["tests", "docs", ".venv", ".git"])
def test_dockerignore_leaves_out_the_heavy_folders(pattern):
    lines = [line.strip() for line in (ROOT / ".dockerignore").read_text().splitlines()]
    assert pattern in lines
