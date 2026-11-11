"""``SessionV2`` and ``SessionHistory`` live in ``astrocook.session`` (phase
6-quater, PR 5); the old paths in ``core/`` are aliases that give the same
objects, with no warning (the scripts and the logs use them)."""
import importlib


def test_the_old_and_the_new_path_give_the_same_objects():
    from astrocook.core import SessionV2 as from_core_package
    from astrocook.core.session import SessionV2 as old_session
    from astrocook.core.session import load_session_from_file as old_load
    from astrocook.core.session_manager import SessionHistory as old_history
    from astrocook.session import SessionHistory, SessionV2
    from astrocook.session.history import SessionHistory as history_module
    from astrocook.session.session import SessionV2 as session_module
    from astrocook.session.session import load_session_from_file

    assert old_session is SessionV2 is session_module is from_core_package
    assert old_history is SessionHistory is history_module
    assert old_load is load_session_from_file
    assert SessionV2.__module__ == "astrocook.session.session"
    assert SessionHistory.__module__ == "astrocook.session.history"
    assert importlib.import_module("astrocook").SessionV2 is SessionV2


def test_an_unknown_name_of_the_old_modules_is_an_attribute_error():
    import astrocook.core
    import astrocook.core.session
    import astrocook.core.session_manager
    for module in (astrocook.core, astrocook.core.session,
                   astrocook.core.session_manager):
        assert not hasattr(module, "NoSuchName")


def test_core_imports_the_session_package_only_lazily():
    """No module-level import of ``astrocook.session`` in ``core/``: the
    aliases import it inside a function, so there is no cycle."""
    import ast
    from pathlib import Path

    import astrocook.core
    for path in Path(astrocook.core.__file__).parent.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in tree.body:
            names = []
            if isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            elif isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            assert not any(n.startswith("astrocook.session") for n in names), path
