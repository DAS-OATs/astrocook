"""Typed parameters of the recipes: one conversion, at the boundary.

The GUI dialogs, the session log, the log window and the scripts pass the
parameters of a recipe as strings (``'5000.0'``, ``'True'``, ``'None'``);
the MCP server and the scripts may pass them typed. The decorator
:func:`typed_recipe` converts them once, to the types declared in the
recipe registry (:mod:`astrocook.recipes.registry`), and checks them
against its bounds before the recipe runs, so that the body of a recipe
receives typed values (phase 6-bis, PR 2).

:func:`resolve_params` gives the values a recipe actually used (the
parameters passed, converted, and the defaults of the others), which the
session log records (#573, item 3).
"""
from __future__ import annotations

import ast
import contextlib
import contextvars
import functools
import inspect
import logging
import math
from typing import Any, Callable, Dict, Mapping

# InvalidParameter is defined with the other exceptions of the recipes
# (phase 6-bis, PR 3) and raised here
from astrocook.recipes.errors import InvalidParameter, RecipeError, RecipeFailed, reason_of

logger = logging.getLogger(__name__)

#: The string that means "no value" for an optional parameter: what the
#: dialogs and the log write for None. An empty field of a dialog (``''``)
#: means None only for a parameter whose default is None (e.g. the start of
#: *Rebin*: the range of the data); elsewhere it is an error, so that a field
#: emptied by mistake does not switch an option off (Guido, phase 6-bis, PR
#: 2a: an empty ``kappa`` of *Rebin* is not "no clipping")
NONE_STRING = 'None'



def _number(param, value, kind):
    if isinstance(value, bool):
        # bool is an int: True for a number is a mistake, not 1
        raise InvalidParameter(f"{param.name}: {value!r} is not a number")
    if isinstance(value, str):
        text = value.strip()
        try:
            return kind(text)
        except ValueError:
            raise InvalidParameter(f"{param.name}: {value!r} is not "
                                   f"{'an integer' if kind is int else 'a number'}") from None
    if kind is int:
        if isinstance(value, float) and not value.is_integer():
            raise InvalidParameter(f"{param.name}: {value!r} is not an integer")
        try:
            return int(value)
        except (TypeError, ValueError):
            raise InvalidParameter(f"{param.name}: {value!r} is not an integer") from None
    try:
        return float(value)
    except (TypeError, ValueError):
        raise InvalidParameter(f"{param.name}: {value!r} is not a number") from None


def _boolean(param, value):
    if isinstance(value, bool):
        return value
    # 1 and 0 too, which scripts often write (Guido, phase 6-bis, PR 2a)
    if isinstance(value, int) and value in (0, 1):
        return bool(value)
    if isinstance(value, str) and value.strip().lower() in ('true', 'false', '1', '0'):
        # The dialogs and the log write 'True' and 'False'; some recipes
        # also read 'true' and 'false'
        return value.strip().lower() in ('true', '1')
    raise InvalidParameter(f"{param.name}: {value!r} is not True or False")


def _list(param, value):
    if isinstance(value, (list, tuple)):
        return list(value)
    if isinstance(value, str):
        try:
            parsed = ast.literal_eval(value)
        except (ValueError, SyntaxError):
            parsed = None
        if isinstance(parsed, (list, tuple)):
            return list(parsed)
    raise InvalidParameter(f"{param.name}: {value!r} is not a list")


def _choice(param, text: str) -> str:
    """One of the choices of ``param``, in its spelling (the case is not
    checked: ``'NESTED'`` is ``'nested'``, as the recipes read it before);
    :class:`InvalidParameter` otherwise (phase 6-bis, PR 6c: an unknown
    value was taken silently, e.g. ``engine='foo'`` ran MCMC)."""
    for choice in param.choices:
        if text.strip().lower() == choice.lower():
            return choice
    allowed = ", ".join(repr(c) for c in param.choices)
    raise InvalidParameter(f"{param.name}: {text!r} is not one of {allowed}")


def convert(param, value: Any) -> Any:
    """``value`` converted to the type of ``param`` (a
    :class:`~astrocook.recipes.registry.ParamSpec`) and checked against its
    bounds; :class:`InvalidParameter` if it cannot be."""
    if isinstance(value, str) and value in param.keep:
        return value  # a special value, such as '_current_'
    if param.optional and param.type is not str and (
            value is None or (isinstance(value, str) and (
                value.strip() == NONE_STRING
                or (value.strip() == '' and param.default is None)))):
        return None
    if value is None:
        if param.optional:
            return None  # e.g. a string that may be absent (phase 6-bis, PR 2c)
        raise InvalidParameter(f"{param.name}: a value is needed")
    if param.type is float:
        result = _number(param, value, float)
    elif param.type is int:
        result = _number(param, value, int)
    elif param.type is bool:
        return _boolean(param, value)
    elif param.type is list:
        return _list(param, value)
    elif param.type is str:
        if param.optional and isinstance(value, str) and value == 'None':
            return None
        text = value if isinstance(value, str) else str(value)
        return _choice(param, text) if param.choices else text
    else:
        return value
    if math.isnan(result):
        if param.allow_nan:
            return result
        raise InvalidParameter(f"{param.name}: {value!r} is not a number")
    if math.isinf(result):
        # Not a value of any parameter so far (Guido, phase 6-bis, PR 2b)
        raise InvalidParameter(f"{param.name}: {value!r} is not a finite number")
    if param.gt is not None and not result > param.gt:
        raise InvalidParameter(f"{param.name}: {result} is not above {param.gt}")
    if param.min is not None and result < param.min:
        raise InvalidParameter(f"{param.name}: {result} is below {param.min}")
    if param.max is not None and result > param.max:
        raise InvalidParameter(f"{param.name}: {result} is above {param.max}")
    return result


def convert_params(spec, params: Mapping[str, Any]) -> Dict[str, Any]:
    """The parameters of ``spec`` in ``params`` converted to their types;
    the others (e.g. ``alias_map``, deprecated parameters) unchanged."""
    declared = {p.name: p for p in spec.params}
    return {name: convert(declared[name], value) if name in declared else value
            for name, value in params.items()}


def resolve_params(name: str, params: Mapping[str, Any]) -> Dict[str, Any]:
    """The values that the recipe ``name`` uses with ``params``: those
    passed, converted, and the defaults of the others, in the order of the
    registry; for the session log (#573, item 3). The parameters passed that
    the registry does not declare (e.g. ``alias_map``), and those of a name
    that is not a recipe, are returned as passed. A parameter passed by its
    old name is returned by the new one (phase 6-quater, PR 2)."""
    from astrocook.recipes.registry import find_recipe
    spec = find_recipe(name)
    if spec is None:
        return dict(params)
    params = spec.canonical(params)
    resolved = {}
    for p in spec.params:
        resolved[p.name] = convert(p, params[p.name]) if p.name in params else p.default
    resolved.update({k: v for k, v in params.items() if k not in resolved})
    return resolved


#: The errors of the core that mean "this data or this value cannot be
#: used": the decorator of the recipes turns them into RecipeFailed (Guido,
#: phase 6-bis, PR 3). The others, usually bugs, go on as they are
DATA_ERRORS = (ValueError, KeyError, OSError)


class _Steps:
    """The recipes run inside the one that :func:`collect_steps` watches."""

    def __init__(self):
        self.depth = 0
        self.names: list = []


_STEPS: contextvars.ContextVar = contextvars.ContextVar('astrocook_recipe_steps',
                                                        default=None)


@contextlib.contextmanager
def collect_steps():
    """Collect, as ``category.name``, the recipes that the next recipe runs
    as its steps (phase 6-bis, PR 4: the ``steps`` of a RecipeResult)."""
    steps = _Steps()
    token = _STEPS.set(steps)
    try:
        yield steps.names
    finally:
        _STEPS.reset(token)


#: How many recipes are running, one inside the other: the prerequisites
#: filled by the outermost are the steps of the session log (PR 5)
_DEPTH: contextvars.ContextVar = contextvars.ContextVar('astrocook_recipe_depth', default=0)


def typed_recipe(method: Callable) -> Callable:
    """Decorator of a recipe method, its boundary.

    The parameters passed are converted to the types of the registry
    before the method runs (:class:`InvalidParameter` if one cannot be, or
    is out of its bounds). Every :class:`RecipeError` that leaves the
    method records the recipe (so that the message names the recipe, and
    the steps of a pipeline); the errors of the core in
    :data:`DATA_ERRORS` become :class:`RecipeFailed`, with the original
    error as their cause; any other exception (a bug) goes on unchanged,
    with a note naming the recipe. Nothing is logged here: the error is
    logged where it is caught (phase 6-bis, PR 3). Before the method runs,
    the prerequisites that the registry declares for it are filled, or
    their absence raised, as the policy of
    :mod:`astrocook.recipes.prerequisites` says (PR 5). A parameter
    passed by its old name (an alias in the registry) is passed to the
    method by its new name (phase 6-quater, PR 2).
    """
    signature = inspect.signature(method)

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        from astrocook.recipes.registry import find_recipe
        spec = find_recipe(method.__name__)
        where = (spec.category, spec.name, spec.brief.rstrip('.'))
        steps = _STEPS.get()
        if steps is not None:
            if steps.depth:
                steps.names.append(f"{spec.category}.{spec.name}")
            steps.depth += 1
        depth = _DEPTH.set(_DEPTH.get() + 1)
        try:
            bound = signature.bind(self, *args, **spec.canonical(kwargs))
            passed = {k: v for k, v in bound.arguments.items() if k != 'self'}
            converted = convert_params(spec, passed)
            from astrocook.recipes.prerequisites import fill
            session = fill(self._session, spec, outermost=_DEPTH.get() == 1)
            if session is not self._session:
                self = type(self)(session)
            return method(self, **converted)
        except RecipeError as e:
            e.enter(*where)
            raise
        except DATA_ERRORS as e:
            error = RecipeFailed(reason_of(e))
            error.enter(*where)
            raise error from e
        except Exception as e:
            e.add_note(f"In the recipe {spec.category}.{spec.name} ({where[2]}).")
            raise
        finally:
            _DEPTH.reset(depth)
            if steps is not None:
                steps.depth -= 1

    wrapper.typed_recipe = True  # checked by tests/test_recipe_params.py
    return wrapper
