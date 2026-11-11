"""Old home of ``SessionV2``, now ``astrocook.session.session`` (phase
6-quater, PR 5). Every name of the new module is served from here, lazily
(PEP 562): ``core/`` imports nothing of ``astrocook.session`` until a name
is asked, so there is no import cycle. The objects are the same, and there
is no warning: the scripts use this path."""
from importlib import import_module
from typing import Any


def __getattr__(name: str) -> Any:
    if name.startswith("__"):
        raise AttributeError(name)
    return getattr(import_module("astrocook.session.session"), name)
