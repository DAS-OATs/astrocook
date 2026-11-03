"""The broad ``except`` clauses of the package can only decrease (phase 6).

Phase 6 reviews the clauses that catch every exception (``except:``,
``except Exception``, ``except BaseException``) one by one: each is
narrowed, removed, kept at a boundary with an explicit log, or discussed
(``ROADMAP.md``, phase 6, "Inventory of the broad except clauses"). This
test holds the number of such clauses per file at the count of the
inventory, so that no new one is added in the meantime and the inventory
stays up to date: a PR that removes or narrows a clause lowers the count
here and marks the row of the inventory as done.

Ruff's ``BLE001`` is not used for this: it does not report a clause that
logs with ``exc_info=True`` or re-raises, so it would miss about a third
of the inventory (the recipe entry points among them).
"""
import ast
from pathlib import Path

import astrocook

PACKAGE = Path(astrocook.__file__).parent

#: Broad except clauses per file, as in the inventory of phase 6.
BROAD_EXCEPTS = {
    "astrocook/_run_astrocook.py": 4,
    "astrocook/core/session.py": 2,
    "astrocook/gui/about_dialog.py": 1,
    "astrocook/gui/continuum_controls.py": 2,
    "astrocook/gui/continuum_editor.py": 3,
    "astrocook/gui/data_inspector.py": 1,
    "astrocook/gui/identification_viewer_dialog.py": 3,
    "astrocook/gui/inspector_table.py": 1,
    "astrocook/gui/log_scripter_dialog.py": 4,
    "astrocook/gui/main_panels.py": 3,
    "astrocook/gui/main_recipes.py": 9,
    "astrocook/gui/main_scripts.py": 3,
    "astrocook/gui/main_sessions.py": 10,
    "astrocook/gui/main_window.py": 7,
    "astrocook/gui/plot_canvas.py": 7,
    "astrocook/gui/plot_mouse.py": 6,
    "astrocook/gui/plot_style.py": 5,
    "astrocook/gui/pyside_plot.py": 15,
    "astrocook/gui/qt_workers.py": 3,
    "astrocook/gui/recipe_dialog.py": 5,
    "astrocook/gui/script_dialog.py": 1,
    "astrocook/mcp_server/errors.py": 1,
    "astrocook/recipes/absorbers.py": 15,
    "astrocook/recipes/continuum.py": 4,
    "astrocook/recipes/edit.py": 6,
    "astrocook/recipes/file.py": 2,
    "astrocook/recipes/flux.py": 5,
}


def _is_broad(handler: ast.ExceptHandler) -> bool:
    if handler.type is None:
        return True
    types = (handler.type.elts if isinstance(handler.type, ast.Tuple)
             else [handler.type])
    return any(isinstance(t, ast.Name) and t.id in ("Exception", "BaseException")
               for t in types)


def _count_broad_excepts() -> dict:
    counts = {}
    for path in sorted(PACKAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        n = sum(1 for node in ast.walk(tree)
                if isinstance(node, ast.ExceptHandler) and _is_broad(node))
        if n:
            counts[f"astrocook/{path.relative_to(PACKAGE).as_posix()}"] = n
    return counts


def test_is_broad():
    tree = ast.parse(
        "try: pass\nexcept: pass\n"
        "try: pass\nexcept Exception: pass\n"
        "try: pass\nexcept (ValueError, BaseException): pass\n"
        "try: pass\nexcept ValueError: pass\n"
        "try: pass\nexcept (KeyError, TypeError): pass\n")
    handlers = [n for n in ast.walk(tree) if isinstance(n, ast.ExceptHandler)]
    assert [_is_broad(h) for h in handlers] == [True, True, True, False, False]


def test_broad_excepts_match_the_inventory():
    counts = _count_broad_excepts()
    added = {f: (BROAD_EXCEPTS.get(f, 0), n) for f, n in counts.items()
             if n > BROAD_EXCEPTS.get(f, 0)}
    assert not added, (
        "New broad except clauses (file: (inventory, now)); catch the "
        f"specific exceptions instead: {added}")
    removed = {f: (n, counts.get(f, 0)) for f, n in BROAD_EXCEPTS.items()
               if counts.get(f, 0) < n}
    assert not removed, (
        "Broad except clauses removed (file: (inventory, now)): lower the "
        "count in BROAD_EXCEPTS and mark the rows of the inventory in "
        f"ROADMAP.md as done: {removed}")

