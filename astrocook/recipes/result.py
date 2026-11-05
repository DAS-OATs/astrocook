"""The structured result of a recipe (phase 6-bis, PR 4).

:func:`run` runs a recipe of the registry and returns a :class:`RecipeResult`:
the new session (or sessions), what changed in it, listed and not only
counted (the components added, removed and changed, the columns, the
properties), the warnings the recipe logged, the recipes it ran as steps,
and the values of the parameters it used. The GUI, the MCP server and the
scripts that want the whole result call it; the methods of the recipes
keep returning the session (Guido, decision 1 of the phase).

What changed is found by comparing the session before and after, in one
place for every recipe (Guido's choice at the start of PR 4: a generic
comparison, not a description written by each recipe). What a recipe
knows and the comparison cannot see (the score of an identification, the
groups a fit could not fit) goes in ``extra``, PR 4b; the fields of what
the algorithms do not provide yet (the uncertainty of the continuum, the
credibility of an identification and its alternatives) are there, empty,
so that phase 7 fills them without changing the interface.
"""
from __future__ import annotations

import contextvars
import logging
import threading
from collections import Counter
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass, field
from typing import Any, Dict, Iterator, List, Mapping, Optional, Tuple

import numpy as np

#: The parameters of a component compared by :func:`compare_sessions`
COMPONENT_FIELDS = ('series', 'z', 'logN', 'b', 'btur', 'resol', 'flags')


@dataclass(frozen=True)
class ComponentRef:
    """A component, as a change names it."""
    uuid: str
    series: str
    z: float
    logN: float
    b: float


@dataclass(frozen=True)
class ComponentChange:
    """A component whose parameters changed: ``changed`` maps a parameter
    to its values before and after."""
    component: ComponentRef
    changed: Mapping[str, Tuple[Any, Any]]


@dataclass(frozen=True)
class Changes:
    """What a recipe changed in a session (:func:`compare_sessions`)."""
    components_added: Tuple[ComponentRef, ...] = ()
    components_removed: Tuple[ComponentRef, ...] = ()
    components_changed: Tuple[ComponentChange, ...] = ()
    columns_added: Tuple[str, ...] = ()
    columns_removed: Tuple[str, ...] = ()
    columns_changed: Tuple[str, ...] = ()
    #: name, object, pixels, z_em, resolution (median R), units: before and after
    properties_changed: Mapping[str, Tuple[Any, Any]] = field(default_factory=dict)

    @property
    def empty(self) -> bool:
        """True if nothing changed."""
        return not (self.components_added or self.components_removed
                    or self.components_changed or self.columns_added
                    or self.columns_removed or self.columns_changed
                    or self.properties_changed)

    def describe(self) -> List[str]:
        """The changes in a few words each, e.g. ``['3 components added',
        "column 'cont' changed"]``."""
        out = []
        for what, items in (("added", self.components_added),
                            ("removed", self.components_removed),
                            ("changed", self.components_changed)):
            if items:
                out.append(f"{len(items)} component{'s' if len(items) > 1 else ''} {what}")
        for what, names in (("added", self.columns_added), ("removed", self.columns_removed),
                            ("changed", self.columns_changed)):
            if names:
                cols = ", ".join(f"'{n}'" for n in names)
                out.append(f"column{'s' if len(names) > 1 else ''} {cols} {what}")
        for name, (old, new) in self.properties_changed.items():
            out.append(f"{name} {old!r} → {new!r}")
        return out


@dataclass(frozen=True)
class RecipeResult:
    """The result of a recipe run by :func:`run`."""
    #: The recipe, as ``category.name``
    recipe: str
    #: Its title (the brief of the registry)
    title: str
    #: The new session; a tuple for a recipe that opens several (e.g.
    #: ``extract_preset`` with two regions); None for a recipe whose result
    #: is not a session (``value``)
    session: Any
    #: The result opens a new session, instead of being the next state of
    #: the same one
    branching: bool
    #: What changed, compared with the session the recipe ran on; one per
    #: session of ``session`` when it is a tuple
    changes: Any
    #: The values of the parameters the recipe used, defaults included
    params: Mapping[str, Any]
    #: The warnings the recipe logged
    warnings: Tuple[str, ...] = ()
    #: The recipes it ran as steps, as ``category.name``, in their order
    steps: Tuple[str, ...] = ()
    #: The result when it is not a session (the redshift of
    #: ``detect_anchor``, whether ``export_ascii`` wrote a file)
    value: Any = None
    #: What the recipe knows and the comparison cannot see (PR 4b):
    #: ``removed`` (uuid -> why a component of the session was removed),
    #: ``discarded`` (uuid -> why a component the recipe made was dropped,
    #: e.g. a narrow line of a pipeline), ``not_fitted``
    #: (the components whose group could not be fitted), ``candidates`` and
    #: ``regions`` (the identification of ``identify_lines``)
    extra: Mapping[str, Any] = field(default_factory=dict)
    #: The recipes run first to provide what the session lacked (PR 5),
    #: as :class:`~astrocook.recipes.prerequisites.ProvidedStep`: each with
    #: its session, for the callers that keep the history and the log
    provided: Tuple[Any, ...] = ()
    #: Reserved for phase 7 (the algorithms do not provide them yet): the
    #: uncertainty of the continuum, the credibility of an identification
    #: and its alternatives
    continuum_uncertainty: Any = None
    credibility: Any = None
    alternatives: Tuple[Any, ...] = ()

    @property
    def sessions(self) -> Tuple[Any, ...]:
        """The new sessions, always as a tuple."""
        if self.session is None:
            return ()
        return self.session if isinstance(self.session, tuple) else (self.session,)

    @property
    def notes(self) -> Tuple[str, ...]:
        """What the user should know of the result besides the changes
        (phase 6-bis, PR 4c, Guido): the components removed and why, those
        a pipeline made and dropped, those whose fit failed, the candidate
        systems found; one line each, e.g. ``'2 components removed: logN <
        min_logN = 14.0 (too weak)'``."""
        from astrocook.recipes.registry import find_recipe
        out = [f"{find_recipe(p.recipe.split('.')[1]).brief.rstrip('.')} run first: "
               f"{p.reason}" for p in self.provided]
        for key, what in (("removed", "removed"),
                          ("discarded", "made and discarded by the recipe")):
            reasons = Counter(self.extra.get(key, {}).values())
            for reason, n in reasons.items():
                out.append(f"{_count(n, 'component')} {what}: {reason}")
        not_fitted = self.extra.get("not_fitted", ())
        if not_fitted:
            out.append(f"{_count(len(not_fitted), 'component')} not fitted: the fit of "
                       f"{'its' if len(not_fitted) == 1 else 'their'} group failed, "
                       f"the values are those from before")
        candidates = self.extra.get("candidates")
        if candidates is not None:
            out.append(f"{_count(len(candidates), 'candidate system')} found")
        return tuple(out)

    @property
    def counts(self) -> Tuple[str, ...]:
        """The notes in a word each, for a status bar, e.g. ``('2
        removed', '1 not fitted')``."""
        out = []
        for key, what in (("removed", "removed"), ("discarded", "discarded"),
                          ("not_fitted", "not fitted"), ("candidates", "candidates")):
            items = self.extra.get(key)
            if items is not None and (len(items) or key == "candidates"):
                out.append(f"{len(items)} {what}")
        return tuple(out)

    @property
    def summary(self) -> str:
        """One line: the recipe, what changed, the warnings."""
        all_changes = self.changes if isinstance(self.changes, tuple) else (self.changes,)
        parts = [p for c in all_changes if c is not None for p in c.describe()]
        if isinstance(self.session, tuple):
            parts.insert(0, f"{len(self.session)} sessions")
        text = f"{self.title}: " + ("; ".join(parts) if parts else "nothing changed")
        if self.warnings:
            n = len(self.warnings)
            text += f", with {n} warning{'s' if n > 1 else ''}"
        return text + "."


def _count(n: int, noun: str) -> str:
    return f"{n} {noun}{'s' if n != 1 else ''}"


# --- The comparison --------------------------------------------------------

def _plain(value):
    """A NumPy scalar as a Python one (the result is also read as JSON)."""
    return value.item() if isinstance(value, np.generic) else value


def _ref(c) -> ComponentRef:
    return ComponentRef(uuid=c.uuid, series=c.series, z=_plain(c.z),
                        logN=_plain(c.logN), b=_plain(c.b))


def _components(session) -> Dict[str, Any]:
    systs = getattr(session, 'systs', None)
    return {c.uuid: c for c in (systs.components if systs else [])}


def _columns(spec) -> Dict[str, Any]:
    if spec is None:
        return {}
    columns = {'x': spec.x, 'y': spec.y, 'dy': spec.dy}
    columns.update(spec._data.aux_cols)
    return columns


def _values(column) -> np.ndarray:
    values = getattr(column, 'values', None)
    if values is None:
        values = getattr(column, 'value', column)
    return np.asarray(values)


def _same_column(a, b) -> bool:
    if str(getattr(a, 'unit', '')) != str(getattr(b, 'unit', '')):
        return False
    va, vb = _values(a), _values(b)
    if va.shape != vb.shape:
        return False
    if va.dtype.kind in 'fc' or vb.dtype.kind in 'fc':
        return bool(np.array_equal(va, vb, equal_nan=True))
    return bool(np.array_equal(va, vb))


def _properties(session) -> Dict[str, Any]:
    from astrocook.core.resolution import resolution_summary
    spec = getattr(session, 'spec', None)
    if spec is None:
        return {'name': getattr(session, 'name', None)}
    return {'name': session.name, 'object': spec.meta.get('OBJECT'),
            'pixels': len(spec.x), 'z_em': float(spec._data.z_em),
            'resolution': float(resolution_summary(spec)[0]),
            'x unit': str(spec.x.unit), 'y unit': str(spec.y.unit)}


def compare_sessions(before, after) -> Changes:
    """What changed from the session ``before`` to ``after``."""
    old, new = _components(before), _components(after)
    added = tuple(_ref(new[u]) for u in new if u not in old)
    removed = tuple(_ref(old[u]) for u in old if u not in new)
    changed = []
    for u in new:
        if u in old:
            diff = {f: (_plain(getattr(old[u], f, None)), _plain(getattr(new[u], f, None)))
                    for f in COMPONENT_FIELDS
                    if getattr(old[u], f, None) != getattr(new[u], f, None)}
            if diff:
                changed.append(ComponentChange(_ref(new[u]), diff))
    cols_old = _columns(getattr(before, 'spec', None))
    cols_new = _columns(getattr(after, 'spec', None))
    p_old, p_new = _properties(before), _properties(after)
    # On another grid (a split, a rebin) every column differs: the change is
    # in the pixels, not in the columns that both have
    same_grid = p_old.get('pixels') == p_new.get('pixels')
    return Changes(
        components_added=added, components_removed=removed,
        components_changed=tuple(changed),
        columns_added=tuple(n for n in cols_new if n not in cols_old),
        columns_removed=tuple(n for n in cols_old if n not in cols_new),
        columns_changed=tuple(n for n in cols_new if same_grid and n in cols_old
                              and not _same_column(cols_old[n], cols_new[n])),
        properties_changed={k: (p_old.get(k), v) for k, v in p_new.items()
                            if p_old.get(k) != v})


# --- What a recipe knows ------------------------------------------------------

_EXTRA: contextvars.ContextVar = contextvars.ContextVar('astrocook_recipe_extra',
                                                        default=None)


def report(key: str, value: Any) -> None:
    """Put ``value`` in the ``extra`` of the result, under ``key``: what a
    recipe knows and the comparison of the sessions cannot see (PR 4b). A
    list value is extended by the steps that report it again. Nothing
    happens outside :func:`run`."""
    extra = _EXTRA.get()
    if extra is None:
        return
    if isinstance(value, list) and isinstance(extra.get(key), list):
        extra[key].extend(value)
    else:
        extra[key] = value


def report_removed(uuids, reason: str) -> None:
    """Say why the components ``uuids`` were removed: ``extra['removed']``
    maps each to its reason (the comparison lists them, not why)."""
    extra = _EXTRA.get()
    if extra is None or not uuids:
        return
    removed = extra.setdefault('removed', {})
    for uuid in uuids:
        removed.setdefault(uuid, reason)


@contextmanager
def _collect_extra() -> Iterator[Dict[str, Any]]:
    extra: Dict[str, Any] = {}
    token = _EXTRA.set(extra)
    try:
        yield extra
    finally:
        _EXTRA.reset(token)


# --- Running a recipe --------------------------------------------------------

class _Warnings(logging.Handler):
    """The warnings logged by Astrocook in one thread."""

    def __init__(self, thread_id: int):
        super().__init__(level=logging.WARNING)
        self.thread_id = thread_id
        self.messages: List[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        if record.thread == self.thread_id and record.levelno == logging.WARNING:
            message = record.getMessage()
            if message not in self.messages:
                self.messages.append(message)


@contextmanager
def _collect_warnings() -> Iterator[List[str]]:
    handler = _Warnings(threading.get_ident())
    logger = logging.getLogger('astrocook')
    logger.addHandler(handler)
    try:
        yield handler.messages
    finally:
        logger.removeHandler(handler)


def run(session, recipe: str, *, on_missing: Optional[str] = None,
        **params) -> RecipeResult:
    """Run the recipe ``recipe`` (``'category.name'``, or the name alone)
    on ``session`` and return its :class:`RecipeResult`.

    The parameters are given as to the method of the recipe (typed, or as
    strings); a recipe that cannot do what was asked raises its
    :class:`~astrocook.recipes.errors.RecipeError`, as the method does.
    ``on_missing`` is the policy for a prerequisite the session lacks,
    ``"provide"`` or ``"raise"`` (:mod:`astrocook.recipes.prerequisites`);
    by default, the one in force.
    """
    from astrocook.recipes.params import collect_steps, resolve_params
    from astrocook.recipes.prerequisites import collect_provided, policy
    from astrocook.recipes.registry import find_recipe, get_recipe
    if '.' in recipe:
        spec = get_recipe(*recipe.split('.', 1))
    else:
        spec = find_recipe(recipe)
        if spec is None:
            raise KeyError(f"There is no recipe {recipe}")
    method = getattr(getattr(session, spec.category), spec.name)
    chosen = policy(on_missing) if on_missing is not None else nullcontext()
    with chosen, _collect_warnings() as warnings, collect_steps() as steps, \
            _collect_extra() as extra, collect_provided() as provided:
        out = method(**params)
    used = resolve_params(spec.name, params)
    if 'removed' in extra:
        # The components of the session the recipe ran on, and those it made
        # and dropped itself (e.g. the narrow lines of a pipeline)
        before = _components(session)
        reasons = extra.pop('removed')
        removed = {u: r for u, r in reasons.items() if u in before}
        discarded = {u: r for u, r in reasons.items() if u not in before}
        if removed:
            extra['removed'] = removed
        if discarded:
            extra['discarded'] = discarded
    if isinstance(out, list):
        out = tuple(out)
    if isinstance(out, tuple):
        new, changes, value = out, tuple(compare_sessions(session, s) for s in out), None
    elif hasattr(out, 'spec'):
        new, changes, value = out, compare_sessions(session, out), None
    else:
        new, changes, value = None, None, out
    return RecipeResult(recipe=f"{spec.category}.{spec.name}", title=spec.brief.rstrip('.'),
                        session=new, branching=spec.branching, changes=changes,
                        params=used, warnings=tuple(warnings), steps=tuple(steps),
                        value=value, extra=dict(extra), provided=tuple(provided))
