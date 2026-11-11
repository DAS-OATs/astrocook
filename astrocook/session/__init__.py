"""Sessions: ``SessionV2`` (the data of one analysis) and ``SessionHistory``
(its undo/redo states).

They lived in ``astrocook/core/`` until phase 6-quater (PR 5). They moved
here because a session runs the recipes, which sit above the core: the core
is the bottom layer and does not depend on this package. The old import
paths (``astrocook.session.session``, ``astrocook.session.history``) still
work, by lazy aliases in those modules, with no warning: the scripts and the
logs of the sessions use them.
"""
from astrocook.session.history import SessionHistory
from astrocook.session.session import SessionV2, load_session_from_file

__all__ = ["SessionV2", "SessionHistory", "load_session_from_file"]
