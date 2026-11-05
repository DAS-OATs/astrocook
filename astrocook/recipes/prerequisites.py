"""What a recipe needs from the session, and who provides it (phase 6-bis, PR 5).

The registry declares what each recipe *requires* (``z_em``, the
``resolution``, the ``continuum``) and what it *provides*; this module
says what is missing for a recipe (:func:`missing`, for the GUI and the
MCP server, which ask before the recipe runs) and fills the gaps at the
boundary of the recipe (:func:`fill`, called by ``typed_recipe``),
following a policy chosen by the script:

- ``"provide"`` (the default): what a recipe can provide with its
  defaults is provided first, with a WARNING (the continuum, by
  *Auto-estimate Continuum*); what no recipe can provide by itself (the
  emission redshift, the resolution) is left to the recipe, which goes on
  as before with its own WARNING (Guido: the numbers do not change);
- ``"raise"``: a :class:`~astrocook.recipes.errors.MissingPrerequisite`
  that says what is missing and which recipe provides it.

The policy is set for a block of a script with :func:`policy`, or for one
call with ``run(..., on_missing="raise")``.
"""
from __future__ import annotations

import contextvars
import logging
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, Callable, Iterator, List, Mapping, Optional, Tuple

from astrocook.recipes.errors import MissingPrerequisite

logger = logging.getLogger(__name__)

#: The policies of the scripts
POLICIES = ("provide", "raise")

_POLICY = contextvars.ContextVar('astrocook_prerequisites_policy', default="provide")
_PROVIDED = contextvars.ContextVar('astrocook_prerequisites_provided', default=None)


def _has_z_em(spec) -> bool:
    return spec._data.z_em != 0.0


def _has_resolution(spec) -> bool:
    from astrocook.core.resolution import resolution_summary
    return resolution_summary(spec)[0] > 0


def _has_continuum(spec) -> bool:
    # A spectrum normalized already (flux ~1) has a continuum of 1
    return spec.norm is not None


@dataclass(frozen=True)
class Prerequisite:
    """Something a recipe can require from the session."""
    name: str
    #: Said when it is missing, e.g. "the spectrum has no continuum"
    missing: str
    #: Whether the session has it (from its spectrum)
    check: Callable[[Any], bool]
    #: The recipe that provides it with its defaults, as ``category.name``;
    #: None when a value must be given (e.g. z_em)
    automatic: Optional[str] = None


PREREQUISITES: Mapping[str, Prerequisite] = {
    "z_em": Prerequisite("z_em", "the emission redshift (z_em) is not set", _has_z_em),
    "resolution": Prerequisite("resolution", "the spectrum has no resolution",
                               _has_resolution),
    "continuum": Prerequisite("continuum", "the spectrum has no continuum", _has_continuum,
                              automatic="continuum.estimate_auto"),
}


@dataclass(frozen=True)
class Missing:
    """A prerequisite that a recipe needs and the session does not have."""
    #: ``z_em``, ``resolution`` or ``continuum``
    name: str
    #: e.g. "the spectrum has no continuum"
    reason: str
    #: The recipes that provide it, as ``category.name``
    providers: Tuple[str, ...]
    #: The one that provides it with its defaults, or None
    automatic: Optional[str]

    @property
    def message(self) -> str:
        """The reason and how to fill the gap, e.g. "The spectrum has no
        continuum: Auto-estimate continuum (single-click)
        (continuum.estimate_auto) provides it with its defaults"; all the
        recipes that provide it are in :attr:`providers`."""
        from astrocook.recipes.registry import find_recipe

        def named(p):
            return f"{find_recipe(p.split('.')[1]).brief.rstrip('.')} ({p})"
        if self.automatic is not None:
            how = f"{named(self.automatic)} provides it with its defaults"
        else:
            how = "set it with " + " or ".join(named(p) for p in self.providers)
        return f"{self.reason[0].upper()}{self.reason[1:]}: {how}"


@dataclass(frozen=True)
class ProvidedStep:
    """A recipe run first to provide a prerequisite (policy "provide")."""
    #: The recipe, as ``category.name``
    recipe: str
    #: The values of its parameters (its defaults), for the session log
    params: Mapping[str, Any]
    #: Why it was run, e.g. "the spectrum has no continuum"
    reason: str
    #: The session it gave, on which the recipe then ran
    session: Any


def providers_of(name: str) -> Tuple[str, ...]:
    """The recipes that declare that they provide ``name``, as
    ``category.name``, the automatic one first."""
    from astrocook.recipes.registry import all_recipes
    found = [f"{s.category}.{s.name}" for s in all_recipes() if name in s.provides]
    automatic = PREREQUISITES[name].automatic
    return tuple(sorted(found, key=lambda p: p != automatic))


def missing(session, recipe) -> Tuple[Missing, ...]:
    """What ``recipe`` (a :class:`~astrocook.recipes.registry.RecipeSpec`)
    needs and ``session`` does not have; pure, nothing is run."""
    spec = getattr(session, 'spec', None)
    if spec is None:
        return ()
    out = []
    for name in sorted(recipe.requires, key=list(PREREQUISITES).index):
        prerequisite = PREREQUISITES[name]
        if not prerequisite.check(spec):
            out.append(Missing(name, prerequisite.missing, providers_of(name),
                               prerequisite.automatic))
    return tuple(out)


@contextmanager
def policy(name: str) -> Iterator[None]:
    """The policy for the recipes run in the block: ``"provide"`` or
    ``"raise"``.

    >>> with policy("raise"):
    ...     session = session.absorbers.refit_all()   # doctest: +SKIP
    """
    if name not in POLICIES:
        raise ValueError(f"Unknown policy {name!r}: use one of {', '.join(POLICIES)}.")
    token = _POLICY.set(name)
    try:
        yield
    finally:
        _POLICY.reset(token)


def current_policy() -> str:
    """The policy in force."""
    return _POLICY.get()


@contextmanager
def collect_provided() -> Iterator[List[ProvidedStep]]:
    """The steps run to provide the prerequisites of the outermost recipe
    called in the block (used by ``run`` and by the callers that write the
    session log)."""
    provided: List[ProvidedStep] = []
    token = _PROVIDED.set(provided)
    try:
        yield provided
    finally:
        _PROVIDED.reset(token)


def fill(session, recipe, outermost: bool = True):
    """The session with the prerequisites of ``recipe`` filled, as the
    policy says (see the module); called by ``typed_recipe`` before the
    recipe runs. Raises :class:`MissingPrerequisite` with the policy
    ``"raise"``."""
    gaps = missing(session, recipe)
    if not gaps:
        return session
    if _POLICY.get() == "raise":
        raise MissingPrerequisite("; ".join(g.message for g in gaps) + ".")
    from astrocook.recipes.params import resolve_params
    from astrocook.recipes.registry import find_recipe
    for gap in gaps:
        if gap.automatic is None:
            continue  # the recipe goes on, with its own WARNING (as before)
        category, name = gap.automatic.split('.')
        title = find_recipe(name).brief.rstrip('.')
        logger.warning(f"{gap.reason[0].upper()}{gap.reason[1:]}: it is estimated first "
                       f"with {title}, with its defaults.")
        session = getattr(getattr(session, category), name)()
        provided = _PROVIDED.get()
        if outermost and provided is not None:
            provided.append(ProvidedStep(recipe=gap.automatic,
                                         params=resolve_params(name, {}),
                                         reason=gap.reason, session=session))
    return session

