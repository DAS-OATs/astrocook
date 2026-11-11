"""Old home of ``SessionHistory``, now ``astrocook.session.history`` (phase
6-quater, PR 5); see ``astrocook/core/session.py``."""
from importlib import import_module
from typing import Any


def __getattr__(name: str) -> Any:
    if name.startswith("__"):
        raise AttributeError(name)
    return getattr(import_module("astrocook.session.history"), name)
