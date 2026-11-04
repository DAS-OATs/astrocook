"""The recipe registry (phase 6-bis, PR 1): one description of each recipe.

The registry (``astrocook.recipes.registry``) is built from the
``*_RECIPES_SCHEMAS`` dictionaries of the recipe modules; these tests
check that it agrees with the recipe classes, that it is read only, and
that the flags moved into it (branching, requirements) are the ones the
code had before (``BRANCHING_RECIPES`` in ``core/utils.py`` and the two
lists of ``gui/main_recipes.py``, as of 2.4.0).
"""
import inspect

import pytest

from astrocook.recipes import registry
from astrocook.recipes.absorbers import RecipeAbsorbersV2
from astrocook.recipes.continuum import RecipeContinuumV2
from astrocook.recipes.edit import RecipeEditV2
from astrocook.recipes.file import RecipeFileV2
from astrocook.recipes.flux import RecipeFluxV2
from astrocook.recipes.registry import RECIPES, build_spec, find_recipe, get_recipe

CLASSES = {"absorbers": RecipeAbsorbersV2, "continuum": RecipeContinuumV2,
           "edit": RecipeEditV2, "file": RecipeFileV2, "flux": RecipeFluxV2}

#: Public methods of the recipe classes that are not recipes
NOT_RECIPES = {("file", "export_ascii_outputs")}  # helper of the GUI export (F82)

#: As in ``BRANCHING_RECIPES`` (core/utils.py) before the registry
BRANCHING = {"split", "stitch", "equalize_and_stitch", "extract_preset", "coadd",
             "resample", "identify_lines"}
#: As in ``MainWindowV2._RECIPES_REQUIRING_Z_EM`` before the registry
REQUIRING_Z_EM = {"fit_powerlaw", "estimate_auto", "find_absorbed", "extract_preset",
                  "doublets_auto", "lya_auto", "forest_auto", "metals_auto",
                  "identify_lines"}
#: As in ``MainWindowV2._RECIPES_REQUIRING_RESOL`` before the registry
REQUIRING_RESOLUTION = {"fit_component", "add_component", "doublets_auto", "lya_auto",
                        "forest_auto", "metals_auto", "identify_lines"}


def test_categories_in_the_order_of_the_menus():
    assert tuple(RECIPES) == registry.CATEGORIES
    assert len(list(registry.all_recipes())) == 48


@pytest.mark.parametrize("category", registry.CATEGORIES)
def test_every_method_has_its_spec_and_every_spec_its_method(category):
    cls = CLASSES[category]
    methods = {name for name, _ in inspect.getmembers(cls, inspect.isfunction)
               if not name.startswith("_") and (category, name) not in NOT_RECIPES}
    assert methods == set(RECIPES[category])


@pytest.mark.parametrize("spec", list(registry.all_recipes()),
                         ids=lambda s: f"{s.category}.{s.name}")
def test_every_parameter_is_in_the_signature(spec):
    signature = inspect.signature(getattr(CLASSES[spec.category], spec.name)).parameters
    for p in spec.params:
        assert p.name in signature, p.name
    assert spec.brief


def test_recipe_names_are_unique_and_found_by_name():
    names = [s.name for s in registry.all_recipes()]
    assert len(names) == len(set(names))
    for spec in registry.all_recipes():
        assert find_recipe(spec.name) is spec
        assert get_recipe(spec.category, spec.name) is spec
    assert find_recipe("no_such_recipe") is None
    with pytest.raises(KeyError, match="no recipe edit.no_such_recipe"):
        get_recipe("edit", "no_such_recipe")


def test_flags_as_before_the_registry():
    specs = list(registry.all_recipes())
    assert {s.name for s in specs if s.branching} == BRANCHING
    assert {s.name for s in specs if "z_em" in s.requires} == REQUIRING_Z_EM
    assert {s.name for s in specs if "resolution" in s.requires} == REQUIRING_RESOLUTION


def test_registry_is_read_only():
    spec = get_recipe("edit", "equalize_and_stitch")
    with pytest.raises(TypeError):
        RECIPES["edit"]["new"] = spec
    with pytest.raises(AttributeError):
        spec.brief = "changed"
    with pytest.raises(AttributeError):
        spec.params[0].doc = "changed"
    # A list default is given as a copy
    knots = get_recipe("continuum", "update_from_knots").param("knots_x")
    knots.default.append(1.0)
    assert knots.default == []


def test_param_lookup_and_defaults():
    spec = get_recipe("continuum", "estimate_auto")
    assert spec.param("kappa").default == 2.0
    assert spec.param("kappa").type is float
    assert spec.defaults()["fudge"] == "auto"
    with pytest.raises(KeyError, match="no parameter 'nope'"):
        spec.param("nope")
    # A parameter without a description is described by its name
    assert get_recipe("absorbers", "delete_component").param("uuid").doc == "uuid"


def test_build_spec_rejects_unknown_keys_and_requirements():
    entry = {"brief": "A recipe.", "params": []}
    assert build_spec("edit", "r", entry).requires == frozenset()
    with pytest.raises(ValueError, match="unknown keys"):
        build_spec("edit", "r", {**entry, "branchng": True})
    with pytest.raises(ValueError, match="unknown keys"):
        build_spec("edit", "r", {**entry, "params": [
            {"name": "x", "type": int, "default": 0, "dok": "x"}]})
    with pytest.raises(ValueError, match="unknown requirements"):
        build_spec("edit", "r", {**entry, "requires": ("continuum",)})


def test_schema_dictionaries_read_only_by_the_registry():
    """Only the registry reads the ``*_RECIPES_SCHEMAS`` dictionaries: a
    reader that changed them would change the description of a recipe for
    the whole process (F12, F79)."""
    from pathlib import Path
    import astrocook
    root = Path(astrocook.__file__).parent
    readers = sorted(str(path.relative_to(root)) for path in root.rglob("*.py")
                     if "_RECIPES_SCHEMAS" in path.read_text(encoding="utf-8")
                     and path.parent.name != "recipes")
    assert readers == []
