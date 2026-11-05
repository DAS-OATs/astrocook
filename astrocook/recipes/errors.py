"""The exceptions of the recipes (phase 6-bis, PR 3).

A recipe that cannot do what was asked raises a :class:`RecipeError`,
instead of logging the reason and returning 0. It does not log it: the
error is logged once, where it is caught (the GUI shows it in a dialog and
the status bar, the MCP server passes it to the agent, the worker of the
scripts logs it with the line of the script; a Python script gets the
exception). Guido's decisions at the start of PR 3.

Every error says in which recipe it happened and why. The decorator of
the recipes (:func:`astrocook.recipes.params.typed_recipe`) records each
recipe the error goes through, so that an error of a step of a pipeline
names the pipeline and the step::

    Auto-fit metal doublets → Refit all systems: the region 'nope' is unknown.

The same decorator turns the errors of the core that mean "this data or
this value cannot be used" (``ValueError``, ``KeyError``, ``OSError``) into
:class:`RecipeFailed`, with the original error as its cause; the others,
usually bugs, go on as they are, with a note naming the recipe.
"""
from __future__ import annotations

from typing import List, Optional, Tuple


class RecipeError(Exception):
    """A recipe could not do what was asked.

    ``reason`` is why; ``chain`` the recipes the error went through, as
    ``(category, name, title)`` from the outermost to the one where it
    happened (filled by the decorator of the recipes). ``str()`` gives the
    titles and the reason, the message shown to the user.
    """

    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason
        self.chain: List[Tuple[str, str, str]] = []

    @property
    def recipe(self) -> Optional[str]:
        """The recipe where the error happened, as ``category.name``."""
        if not self.chain:
            return None
        category, name, _ = self.chain[-1]
        return f"{category}.{name}"

    def enter(self, category: str, name: str, title: str) -> None:
        """Record a recipe the error goes through (called from the inside
        out, so each new one is the outermost)."""
        self.chain.insert(0, (category, name, title))

    def __str__(self) -> str:
        reason = self.reason if self.reason.endswith(('.', '!', '?')) else self.reason + '.'
        if not self.chain:
            return reason
        return " → ".join(title for _, _, title in self.chain) + ": " + reason


class InvalidParameter(RecipeError, ValueError):
    """A parameter cannot be converted to its type, or is out of its
    bounds; the reason names the parameter and the value."""


class RecipeFailed(RecipeError):
    """The recipe ran and could not give a result: the data or the session
    do not allow it (a session that is not open, a file that cannot be
    read, an expression that cannot be evaluated, ...)."""


class MissingPrerequisite(RecipeError):
    """The session lacks what the recipe needs (z_em, the resolution, the
    continuum, ...); phase 6-bis, PR 5, fills these in ("Prerequisites of
    the recipes")."""


def reason_of(error: BaseException) -> str:
    """The text of an error of the core, as a reason (a ``KeyError`` gives
    its key quoted, not the repr of its arguments)."""
    if isinstance(error, KeyError) and error.args:
        return f"{error.args[0]!r} not found"
    return str(error) or type(error).__name__
