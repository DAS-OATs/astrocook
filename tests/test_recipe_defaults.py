"""The defaults of the recipes: one value for the GUI and the API (#573, F40).

Each recipe defines the default of a parameter twice: in its schema
(the registry, ``astrocook.recipes.registry``, which fills the GUI
dialogs and the tools of the MCP server) and in the signature of its method (scripts, and the GUI
actions that pass only some parameters, such as "Refit" and "Fit
component" in the System Inspector). The schema is the source of truth
(agreed with Guido): every optional parameter of a schema must have the
same default in the signature. The differences found in phase 3 are
strict xfails until Guido decides each value; a decided value turns its
xfail into a pass (and strict makes it fail until it is removed from
``KNOWN_DIFFERENCES``).

Parameters that the API requires and the dialog pre-fills (e.g. the
expression of ``apply_expression``) are not compared.
"""
import inspect
import math

import pytest

from astrocook.recipes.absorbers import RecipeAbsorbersV2
from astrocook.recipes.continuum import RecipeContinuumV2
from astrocook.recipes.edit import RecipeEditV2
from astrocook.recipes.file import RecipeFileV2
from astrocook.recipes.flux import RecipeFluxV2
from astrocook.recipes.registry import RECIPES

CLASSES = {"absorbers": RecipeAbsorbersV2, "continuum": RecipeContinuumV2,
           "edit": RecipeEditV2, "file": RecipeFileV2, "flux": RecipeFluxV2}

#: (menu, recipe, parameter): schema default, signature default, as found
#: in phase 3 (issue #573). To be decided by Guido, one by one.
KNOWN_DIFFERENCES = {}
# The 15 differences found in phase 3 were decided by Guido and aligned in
# phase 4 (#573): see ROADMAP_HISTORY.md, F40, and CHANGELOG.md


def _schema_parameters():
    for menu, specs in RECIPES.items():
        for recipe, spec in specs.items():
            for param in spec.params:
                yield menu, recipe, param, CLASSES[menu]


def _cases():
    for menu, recipe, param, cls in _schema_parameters():
        key = (menu, recipe, param.name)
        marks = ([pytest.mark.xfail(strict=True, reason="#573: to be decided by Guido")]
                 if key in KNOWN_DIFFERENCES else [])
        yield pytest.param(menu, recipe, param, cls, marks=marks,
                           id=".".join(key))


@pytest.mark.parametrize("menu, recipe, param, cls", list(_cases()))
def test_schema_default_is_the_signature_default(menu, recipe, param, cls):
    method = getattr(cls, recipe)
    signature = inspect.signature(method).parameters
    assert param.name in signature, "not in the signature"
    default = signature[param.name].default
    if default is inspect.Parameter.empty:
        return  # required by the API, pre-filled by the dialog
    # Typed signatures (phase 6-bis, PR 2): the same value, not only the
    # same string
    both_nan = (isinstance(default, float) and isinstance(param.default, float)
                and math.isnan(default) and math.isnan(param.default))
    assert both_nan or (default == param.default
                        and type(default) is type(param.default)), (
        f"GUI {param.default!r}, API {default!r}")


def test_every_schema_recipe_has_a_method():
    for menu, specs in RECIPES.items():
        for recipe in specs:
            assert callable(getattr(CLASSES[menu], recipe, None)), f"{menu}.{recipe}"


def test_known_differences_are_schema_parameters():
    # A row left here after its recipe or parameter is renamed would hide
    # nothing: catch it
    present = {(m, r, p.name) for m, r, p, _ in _schema_parameters()}
    assert set(KNOWN_DIFFERENCES) <= present
