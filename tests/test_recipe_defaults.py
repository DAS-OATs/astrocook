"""The defaults of the recipes: one value for the GUI and the API (#573, F40).

Each recipe defines the default of a parameter twice: in its schema
(``*_RECIPES_SCHEMAS``, which fills the GUI dialogs and the tools of the
MCP server) and in the signature of its method (scripts, and the GUI
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

import pytest

from astrocook.recipes.absorbers import ABSORBERS_RECIPES_SCHEMAS, RecipeAbsorbersV2
from astrocook.recipes.continuum import CONTINUUM_RECIPES_SCHEMAS, RecipeContinuumV2
from astrocook.recipes.edit import EDIT_RECIPES_SCHEMAS, RecipeEditV2
from astrocook.recipes.file import FILE_RECIPES_SCHEMAS, RecipeFileV2
from astrocook.recipes.flux import FLUX_RECIPES_SCHEMAS, RecipeFluxV2

MENUS = {
    "absorbers": (ABSORBERS_RECIPES_SCHEMAS, RecipeAbsorbersV2),
    "continuum": (CONTINUUM_RECIPES_SCHEMAS, RecipeContinuumV2),
    "edit": (EDIT_RECIPES_SCHEMAS, RecipeEditV2),
    "file": (FILE_RECIPES_SCHEMAS, RecipeFileV2),
    "flux": (FLUX_RECIPES_SCHEMAS, RecipeFluxV2),
}

#: (menu, recipe, parameter): schema default, signature default, as found
#: in phase 3 (issue #573). To be decided by Guido, one by one.
KNOWN_DIFFERENCES = {
    ("absorbers", "refit_all", "group_depth"): (2, "0"),
    ("absorbers", "refit_all", "max_nfev"): (200, "100"),
    ("absorbers", "fit_component", "max_nfev"): (200, "2000"),
    ("absorbers", "optimize_system", "max_components"): (3, "5"),
    ("absorbers", "optimize_system", "threshold_sigma"): (3.0, "2.5"),
    ("absorbers", "optimize_system", "aic_penalty"): (5.0, "0.0"),
    ("absorbers", "populate_from_identification", "series_map_json"): ("", None),
    ("absorbers", "populate_from_identification", "z_map_json"): ("", None),
    ("absorbers", "fit_bayesian", "uuid"): ("", "not in the signature"),
    ("continuum", "estimate_auto", "fudge"): ("auto", "1.0"),
    ("edit", "x_convert", "z_ref"): ("_current_", "0.0"),
    ("edit", "coadd", "dx"): (10.0, "0.01"),
    ("edit", "coadd", "xunit"): ("km/s", "nm"),
    ("file", "import_ascii", "path"): ("", "not in the signature"),
    ("flux", "rebin", "kappa"): ("None", "5.0"),
}


def _schema_parameters():
    for menu, (schemas, cls) in MENUS.items():
        for recipe, schema in schemas.items():
            for param in schema["params"]:
                yield menu, recipe, param, cls


def _same(schema_default, signature_default, kind):
    """The two defaults as the recipe would see them: the API takes
    strings ('2.0', 'True'), the schema typed values."""
    if kind in (int, float):
        try:
            return float(schema_default) == float(signature_default)
        except (TypeError, ValueError):
            pass
    return str(schema_default) == str(signature_default)


def _cases():
    for menu, recipe, param, cls in _schema_parameters():
        key = (menu, recipe, param["name"])
        marks = ([pytest.mark.xfail(strict=True, reason="#573: to be decided by Guido")]
                 if key in KNOWN_DIFFERENCES else [])
        yield pytest.param(menu, recipe, param, cls, marks=marks,
                           id=".".join(key))


@pytest.mark.parametrize("menu, recipe, param, cls", list(_cases()))
def test_schema_default_is_the_signature_default(menu, recipe, param, cls):
    method = getattr(cls, recipe)
    signature = inspect.signature(method).parameters
    assert param["name"] in signature, "not in the signature"
    default = signature[param["name"]].default
    if default is inspect.Parameter.empty:
        return  # required by the API, pre-filled by the dialog
    assert _same(param["default"], default, param["type"]), (
        f"GUI {param['default']!r}, API {default!r}")


def test_every_schema_recipe_has_a_method():
    for menu, (schemas, cls) in MENUS.items():
        for recipe in schemas:
            assert callable(getattr(cls, recipe, None)), f"{menu}.{recipe}"


def test_known_differences_are_schema_parameters():
    # A row left here after its recipe or parameter is renamed would hide
    # nothing: catch it
    present = {(m, r, p["name"]) for m, r, p, _ in _schema_parameters()}
    assert set(KNOWN_DIFFERENCES) <= present
