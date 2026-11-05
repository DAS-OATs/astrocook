"""The JSON schema of the parameters of a recipe, from the registry.

The tools of the full MCP server are to be generated from the recipe
registry (``ROADMAP.md``, section 4): this is the part that does not depend
on the server, written and tested on every recipe in phase 6-bis (PR 8) to
show that the registry is enough. The pilot keeps its five tools written by
hand; :func:`parameters_schema` checks that they offer what the registry
says.

What the schema says, for each parameter offered to agents (``agents``,
phase 6-bis PR 6a): its type, its description, its default (a value
that cannot be written in JSON, such as NaN, is left out), its bounds,
the values allowed (``enum``, from ``choices``) and, as
``x-astrocook-kind``, what it names (an open session, a column of the
spectrum). A parameter without a default is required.
"""
import math
from typing import Any, Dict

#: JSON types of the Python types of the registry
JSON_TYPES = {float: "number", int: "integer", bool: "boolean", str: "string",
              list: "array"}


def _json_default(value: Any) -> Any:
    """The default as JSON, or ``...`` if it cannot be written."""
    if isinstance(value, float) and not math.isfinite(value):
        return ...
    return value


def parameter_schema(param) -> Dict[str, Any]:
    """The schema of one parameter (a
    :class:`~astrocook.recipes.registry.ParamSpec`)."""
    schema: Dict[str, Any] = {"type": JSON_TYPES[param.type], "description": param.doc}
    if param.optional:
        schema["type"] = [schema["type"], "null"]
    if param.keep:
        # Special strings (e.g. '_current_', the value of the session) of a
        # parameter that is not a string
        types = schema["type"] if isinstance(schema["type"], list) else [schema["type"]]
        if "string" not in types:
            schema["type"] = types + ["string"]
    if param.choices:
        schema["enum"] = list(param.choices)
    if param.min is not None:
        schema["minimum"] = param.min
    if param.max is not None:
        schema["maximum"] = param.max
    if param.gt is not None:
        schema["exclusiveMinimum"] = param.gt
    if param.kind:
        schema["x-astrocook-kind"] = param.kind
    default = _json_default(param.default)
    if default is not ... and not _required(param):
        schema["default"] = default
    return schema


def _required(param) -> bool:
    """No default: the dialogs fill it, a script or an agent must give it."""
    return param.default in ("", None) and not param.optional


def parameters_schema(spec) -> Dict[str, Any]:
    """The JSON schema (an object) of the parameters of a recipe (a
    :class:`~astrocook.recipes.registry.RecipeSpec`) offered to agents."""
    offered = [p for p in spec.params if p.agents]
    return {
        "type": "object",
        "description": f"{spec.brief} {spec.details}".strip(),
        "properties": {p.name: parameter_schema(p) for p in offered},
        "required": [p.name for p in offered if _required(p)],
        "additionalProperties": False,
    }
