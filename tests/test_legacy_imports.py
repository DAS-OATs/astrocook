"""No new dependencies of the V2 code on ``astrocook.legacy`` (phase 3).

Phase 3 removes the imports from ``astrocook.legacy`` one at a time
(``ROADMAP.md``, section 3). This test lists the ones still allowed: an
import that is not in the list fails the test, and so does an entry of the
list that is no longer needed, so that the list only shrinks. When a PR
removes a dependency, it removes its entries here.

Imports inside functions and ``try`` blocks count as well.
"""
import ast
from pathlib import Path

import astrocook

PACKAGE_DIR = Path(astrocook.__file__).parent

# (module, imported legacy module) pairs still allowed: none since the
# native reader of the V1 archives (io/v1_archive.py)
ALLOWED = set()


def _imported_modules(path):
    """Absolute names of the modules imported by a file."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    package = ".".join(path.relative_to(PACKAGE_DIR.parent).parent.parts)
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                # Relative import: resolve it against the file's package
                base = package.split(".")
                base = base[:len(base) - node.level + 1]
                name = ".".join(base + ([node.module] if node.module else []))
            else:
                name = node.module
            yield name
            # `from astrocook import legacy` and `from astrocook.legacy import x`
            for alias in node.names:
                yield f"{name}.{alias.name}"


def _legacy_imports():
    found = set()
    for path in sorted(PACKAGE_DIR.rglob("*.py")):
        rel = path.relative_to(PACKAGE_DIR.parent)
        if rel.parts[1] == "legacy":
            continue
        for name in _imported_modules(path):
            parts = name.split(".")
            if parts[:2] == ["astrocook", "legacy"] and len(parts) >= 3:
                found.add((rel.as_posix(), ".".join(parts[:3])))
            elif parts == ["astrocook", "legacy"]:
                found.add((rel.as_posix(), name))
    return found


def test_no_new_legacy_imports():
    new = _legacy_imports() - ALLOWED
    assert not new, (
        "New imports from astrocook.legacy (phase 3 removes them; use or "
        f"write a V2 equivalent): {sorted(new)}")


def test_tests_use_the_v1_reference():
    # The parity tests compare V2 with the V1 copies in tests/reference/v1,
    # which stay when astrocook/legacy goes
    tests_dir = Path(__file__).parent
    found = []
    for path in sorted(tests_dir.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            names = ([a.name for a in node.names] if isinstance(node, ast.Import)
                     else [node.module or ""] if isinstance(node, ast.ImportFrom)
                     else [])
            found += [(path.name, n) for n in names
                      if n == "astrocook.legacy" or n.startswith("astrocook.legacy.")]
    assert not found, f"Use reference.v1 instead of astrocook.legacy: {found}"


def test_allowed_list_is_up_to_date():
    stale = ALLOWED - _legacy_imports()
    assert not stale, (
        "These imports from astrocook.legacy are gone: remove them from "
        f"ALLOWED in {Path(__file__).name}: {sorted(stale)}")
