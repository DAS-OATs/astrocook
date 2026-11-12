"""The checks of the parameters that name components, series and links
(phase 6-quater, PR 10).

The manual fits of Absorbers take the UUIDs of components, the names of
series and the expressions of links. Before PR 10 a UUID that names no
component was taken silently (``update_components``, ``delete_component``,
``fit_component``, ``optimize_system``: a state that changes nothing;
``update_component``: the same, with a WARNING; ``update_constraint``: a
constraint stored for no component; ``fit_bayesian``: "cond is not defined
on empty arrays"), and an unknown series was stored (a component without
lines, which no fit can use): F142, F143. Here they are refused, with
:class:`~astrocook.recipes.errors.InvalidParameter`, before the recipe runs.
"""
from __future__ import annotations

import re
from typing import Iterable, List, Optional

from astrocook.core.atomic_data import ATOM_DATA, STANDARD_MULTIPLETS, series_lines
from astrocook.recipes.errors import InvalidParameter

#: The components named by an expression of a link: ``p['<uuid>']``
_NAMED = re.compile(r"""p\[\s*(['"])(.*?)\1\s*\]""")


def check_uuids(session, uuids: Iterable[str], param: str) -> List[str]:
    """``uuids`` as a list; :class:`InvalidParameter` naming ``param`` if
    one of them names no component of the session."""
    uuids = list(uuids)
    known = {c.uuid for c in session.systs.components} if session.systs else set()
    unknown = [u for u in uuids if u not in known]
    if unknown:
        listed = ", ".join(repr(u) for u in unknown[:5])
        if len(unknown) > 5:
            listed += f" and {len(unknown) - 5} more"
        raise InvalidParameter(f"{param}: there is no component {listed} in the session")
    return uuids


def check_series(series: str, param: str, several: str = "") -> str:
    """``series`` if every transition of it is known (a transition, a
    standard multiplet, or a custom one written as its transitions
    separated by commas); :class:`InvalidParameter` naming ``param``
    otherwise. ``several``: the recipe to name when the series is a list
    of series (e.g. ``'CIV,SiIV'``)."""
    if not series.strip():
        raise InvalidParameter(f"{param}: give a transition or a multiplet")
    unknown = [t for t in series_lines(series) if t not in ATOM_DATA]
    if not unknown:
        return series
    if several and ',' in series and all(t in STANDARD_MULTIPLETS for t in unknown):
        raise InvalidParameter(f"{param}: {series!r} is a list of series, not a multiplet; "
                               f"to add several series at the same redshift use {several}")
    raise InvalidParameter(f"{param}: {', '.join(repr(t) for t in unknown)} is not a "
                           f"known transition or multiplet")


def link_target(uuid: str, expression: str, target_uuid: Optional[str]) -> Optional[str]:
    """The component a link of ``uuid`` points to (phase 6-quater, PR 10,
    F146): ``target_uuid`` if it is named by ``expression``, the first
    component the expression names if ``target_uuid`` is not given (the
    fit puts the target in the group of the component, so a link made
    without it, by a script or an agent, fitted its component without its
    target), None if it names none (a constant). :class:`InvalidParameter`
    if the expression names the component itself, or if ``target_uuid`` is
    not named by it."""
    named = [m.group(2) for m in _NAMED.finditer(expression)]
    if uuid in named:
        raise InvalidParameter(f"expression: {expression!r} links the component to itself")
    if target_uuid is None:
        return named[0] if named else None
    if target_uuid not in named:
        raise InvalidParameter(f"target_uuid: {target_uuid!r} is not the component of the "
                               f"expression {expression!r}")
    return target_uuid
