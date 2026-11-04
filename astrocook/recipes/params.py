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
import functools
import inspect
import logging
import math
from typing import Any, Callable, Dict, Mapping

from astrocook.recipes.messages import msg_param_fail

logger = logging.getLogger(__name__)

#: The strings that mean "no value" for an optional parameter (``'None'``
#: is what the dialogs and the log write for None; an empty field of a
#: dialog is ``''``)
NONE_STRINGS = ('None', '')

#: The names of the recipes whose parameters are typed (filled by
#: :func:`typed_recipe`): the others still receive and log what is passed,
#: until phase 6-bis types them all (PR 2b, 2c)
TYPED_RECIPES: set = set()


class InvalidParameter(ValueError):
    """A parameter of a recipe cannot be converted to its type, or is out
    of its bounds; the message names the parameter and the value."""


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
    if isinstance(value, str) and value.strip().lower() in ('true', 'false'):
        # The dialogs and the log write 'True' and 'False'; some recipes
        # also read 'true' and 'false'
        return value.strip().lower() == 'true'
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


def convert(param, value: Any) -> Any:
    """``value`` converted to the type of ``param`` (a
    :class:`~astrocook.recipes.registry.ParamSpec`) and checked against its
    bounds; :class:`InvalidParameter` if it cannot be."""
    if isinstance(value, str) and value in param.keep:
        return value  # a special value, such as '_current_'
    if value is None or (isinstance(value, str) and value.strip() in NONE_STRINGS
                         and param.type is not str):
        if param.optional:
            return None
        if value is None:
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
        return value if isinstance(value, str) else str(value)
    else:
        return value
    if math.isnan(result) and not param.allow_nan:
        raise InvalidParameter(f"{param.name}: {value!r} is not a number")
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
    registry; for the session log (#573, item 3). The parameters of a recipe
    that is not typed yet, and the ones passed that the registry does not
    declare (e.g. ``alias_map``), are returned as passed."""
    from astrocook.recipes.registry import find_recipe
    spec = find_recipe(name)
    if spec is None or name not in TYPED_RECIPES:
        return dict(params)
    resolved = {}
    for p in spec.params:
        resolved[p.name] = convert(p, params[p.name]) if p.name in params else p.default
    resolved.update({k: v for k, v in params.items() if k not in resolved})
    return resolved


def typed_recipe(method: Callable) -> Callable:
    """Decorator of a recipe method: the parameters passed are converted to
    the types of the registry before the method runs. A value that cannot
    be converted, or is out of bounds, stops the recipe as a failure (the
    reason logged at ERROR, 0 returned; an exception in PR 3 of the phase).
    """
    TYPED_RECIPES.add(method.__name__)
    signature = inspect.signature(method)

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        from astrocook.recipes.registry import find_recipe
        spec = find_recipe(method.__name__)
        bound = signature.bind(self, *args, **kwargs)
        passed = {k: v for k, v in bound.arguments.items() if k != 'self'}
        try:
            converted = convert_params(spec, passed)
        except InvalidParameter as e:
            logger.error(f"{msg_param_fail} {e}.")
            return 0
        return method(self, **converted)

    return wrapper
