"""The other sessions that a recipe reads (phase 6-quater, PR 9).

One rule for the parameters of kind ``session`` (one open session) and
``sessions`` (several, comma-separated) of the registry, agreed with Guido:

- a session is named by its name among the open sessions: in the GUI the
  session list, in the MCP server its store (both are
  ``session._gui.session_histories``, by ``display_name``); the name of
  the session the recipe runs on is that session, also without a GUI;
- in a script, a parameter takes the sessions themselves too (``Session``
  objects), and the log records their names (:func:`session_names`);
- a name not found makes the recipe fail, saying which (before, each recipe
  did its own: skipped without a word, skipped with a WARNING, or an
  ``AttributeError`` without a GUI, F135, F136).
"""
from __future__ import annotations

import ast
from typing import TYPE_CHECKING, Any, List

from astrocook.recipes.errors import InvalidParameter, RecipeFailed
from astrocook.recipes.messages import msg_session_not_found

if TYPE_CHECKING:
    from astrocook.session.session import Session

#: Said when a name is not found and there is no GUI (a script)
SCRIPT_HINT = " In a script, pass the session itself."


def is_session(value: Any) -> bool:
    """Whether ``value`` is a session (a ``Session`` object)."""
    from astrocook.session.session import Session
    return isinstance(value, Session)


def split_names(value: Any) -> List[Any]:
    """The items of a parameter of kind ``sessions``: the names of a
    comma-separated string (or of a string that writes a list, as the logs
    of ``stitch`` before PR 9 held it), or the items of a list."""
    if isinstance(value, str):
        text = value.strip()
        if text.startswith(('[', '(')):
            try:
                parsed = ast.literal_eval(text)
            except (ValueError, SyntaxError):
                parsed = None
            if isinstance(parsed, (list, tuple)):
                return [str(v).strip() for v in parsed if str(v).strip()]
        return [n.strip() for n in text.split(',') if n.strip()]
    if isinstance(value, (list, tuple)):
        return [v.strip() if isinstance(v, str) else v for v in value
                if not (isinstance(v, str) and not v.strip())]
    return [value]


def session_names(value: Any) -> Any:
    """``value`` with the sessions replaced by their names, for the log: a
    session for a parameter of kind ``session``, a comma-separated string
    for one of kind ``sessions`` given as a list."""
    if is_session(value):
        return value.name
    if isinstance(value, (list, tuple)):
        return ", ".join(v.name if is_session(v) else str(v) for v in value)
    return value


def open_session(session: 'Session', item: Any, param: str) -> 'Session':
    """The session named ``item`` (or ``item`` itself, a session), for the
    recipe that runs on ``session``; RecipeFailed if no open session has
    that name."""
    if is_session(item):
        return item
    if not isinstance(item, str):
        raise InvalidParameter(f"{param}: {item!r} is not the name of a session")
    name = item.strip()
    gui = getattr(session, '_gui', None)
    for hist in getattr(gui, 'session_histories', None) or []:
        if hist.display_name == name:
            return hist.current_state
    if name == session.name:
        return session
    hint = "" if gui is not None else SCRIPT_HINT
    raise RecipeFailed(msg_session_not_found.format(name=name) + hint)


def open_sessions(session: 'Session', value: Any, param: str) -> List['Session']:
    """The sessions named by ``value`` (see :func:`split_names`), in their
    order; InvalidParameter if it names none."""
    items = split_names(value)
    if not items:
        raise InvalidParameter(f"{param}: name at least one session")
    return [open_session(session, item, param) for item in items]
