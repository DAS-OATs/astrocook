"""Sessions: ``Session`` (the data of one analysis) and ``SessionHistory``
(its undo/redo states).

They lived in ``astrocook/core/`` until phase 6-quater (PR 5). They moved
here because a session runs the recipes, which sit above the core: the core
is the bottom layer and does not depend on this package. The old import
paths (``astrocook.session.session``, ``astrocook.session.history``) still
work, by lazy aliases in those modules, with no warning: the scripts and the
logs of the sessions use them.
"""
from astrocook.session.history import SessionHistory
from astrocook._renamed import renamed
from astrocook.session.session import Session, load_session_from_file

__all__ = ["Session", "SessionHistory", "load_session_from_file"]

# Old name (phase 6-quater, PR 6): see astrocook/_renamed.py
__getattr__ = renamed(__name__, {"SessionV2": "Session"})
