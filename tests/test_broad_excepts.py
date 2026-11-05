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
    "astrocook/_run_astrocook.py": 2,
    "astrocook/core/session.py": 2,
    "astrocook/gui/main_recipes.py": 1,
    "astrocook/gui/main_sessions.py": 2,
    "astrocook/gui/main_window.py": 1,
    "astrocook/gui/plot_canvas.py": 1,
    "astrocook/gui/pyside_plot.py": 1,
    "astrocook/gui/qt_workers.py": 3,
    "astrocook/gui/script_dialog.py": 1,
    "astrocook/mcp_server/errors.py": 1,
    # The boundaries of the cleaning of the narrow lines of lya_auto and
    # forest_auto (phase 6, Guido's choice: the result is kept, with a
    # WARNING); the 13 others went in phase 6-bis (PR 3c)
    "astrocook/recipes/absorbers.py": 2,
    # The boundary of the recipes (typed_recipe): it names the recipe in
    # a note and re-raises (phase 6-bis, PR 3, which removed the 11 of
    # continuum.py, file.py and flux.py in 3a, the 6 of edit.py in 3b, 13 of
    # absorbers.py in 3c)
    "astrocook/recipes/params.py": 1,
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

