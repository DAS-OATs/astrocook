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


#: Parameters of the methods of the recipes that are not in the registry,
#: each with its reason. F106 (phase 6-ter): ``prefix`` of ``export_ascii``
#: was missing, and nothing checked this direction, only that the registry
#: names parameters of the method
NOT_IN_REGISTRY = {
    ("edit", "apply_expression", "alias_map"): "filled by the GUI from the names of the sessions",
    ("edit", "mask_expression", "alias_map"): "filled by the GUI from the names of the sessions",
    ("edit", "split", "alias_map"): "filled by the GUI from the names of the sessions",
    ("flux", "rebin", "norm"): "no longer offered (#560, F37)",
    ("continuum", "estimate_auto", "smooth_std"): "deprecated, kept for old scripts",
    ("continuum", "fit_continuum", "smooth_std"): "deprecated, kept for old scripts",
    ("continuum", "fit_continuum", "template"): "not offered since 6-bis; to review in 6-quater",
    ("absorbers", "update_components", "series"): "kept a string, 'None' for unchanged",
}


@pytest.mark.parametrize("spec", list(registry.all_recipes()),
                         ids=lambda s: f"{s.category}.{s.name}")
def test_every_parameter_of_the_method_is_in_the_registry(spec):
    signature = inspect.signature(getattr(CLASSES[spec.category], spec.name)).parameters
    declared = {p.name for p in spec.params}
    for name, p in signature.items():
        if name == "self" or p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD):
            continue
        key = (spec.category, spec.name, name)
        assert name in declared or key in NOT_IN_REGISTRY, (
            f"{spec.category}.{spec.name}: the parameter {name} is not in the registry")
        # The exceptions have a default: the registry alone can run the recipe
        if key in NOT_IN_REGISTRY:
            assert p.default is not inspect.Parameter.empty, key


def test_the_exceptions_are_still_parameters():
    for (category, name, param), reason in NOT_IN_REGISTRY.items():
        assert param in inspect.signature(getattr(CLASSES[category], name)).parameters
        assert reason


def test_recipe_names_are_unique_and_found_by_name():
    names = [s.name for s in registry.all_recipes()]
    assert len(names) == len(set(names))
    for spec in registry.all_recipes():
        assert find_recipe(spec.name) is spec
        assert get_recipe(spec.category, spec.name) is spec
    assert find_recipe("no_such_recipe") is None
    with pytest.raises(KeyError, match="no recipe edit.no_such_recipe"):
        get_recipe("edit", "no_such_recipe")


#: Added in phase 6-bis (PR 5): every recipe that fits needs the resolution
#: (F99), components_auto z_em as the pipelines that wrap it
ADDED_Z_EM = {"components_auto"}
ADDED_RESOLUTION = {"components_auto", "refit_all", "optimize_system", "fit_bayesian",
                    "add_linked_system"}


def test_flags_as_before_the_registry():
    specs = list(registry.all_recipes())
    assert {s.name for s in specs if s.branching} == BRANCHING
    assert {s.name for s in specs if "z_em" in s.requires} == REQUIRING_Z_EM | ADDED_Z_EM
    assert ({s.name for s in specs if "resolution" in s.requires}
            == REQUIRING_RESOLUTION | ADDED_RESOLUTION)


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
    # A parameter without a description is described by its name (none of
    # the registry is, since phase 6-bis PR 6: see the test below)
    bare = build_spec("edit", "r", {"brief": "B.", "params": [
        {"name": "uuid", "type": str, "default": ""}]})
    assert bare.param("uuid").doc == "uuid"


def test_every_recipe_and_parameter_is_described():
    """The texts are what the dialogs show and what a model reads (phase
    6-bis, PR 6, F48): a recipe has its details, a parameter a sentence
    that is not its name; the parameters not offered to agents are those
    agreed with Guido (debug, questions of the GUI)."""
    for spec in registry.all_recipes():
        assert spec.brief.endswith(".") and spec.details, spec.name
        for p in spec.params:
            assert p.doc != p.name and p.doc.endswith((".", ")")), (spec.name, p.name)
    hidden = {(s.name, p.name) for s in registry.all_recipes() for p in s.params
              if not p.agents}
    assert hidden == {("identify_lines", "debug_rating"), ("identify_lines", "bypass_scoring"),
                      ("estimate_auto", "renorm_model"), ("fit_continuum", "renorm_model"),
                      ("update_from_knots", "renorm_model"), ("smooth_column", "renorm_model")}
    # F48: the mask is True for the absorbed pixels, as find_absorbed makes it
    for name in ("identify_lines", "fit_continuum"):
        assert "True for the absorbed pixels" in get_recipe(
            "absorbers" if name == "identify_lines" else "continuum", name
        ).param("mask_col").doc


def test_build_spec_rejects_unknown_keys_and_requirements():
    entry = {"brief": "A recipe.", "params": []}
    assert build_spec("edit", "r", entry).requires == frozenset()
    with pytest.raises(ValueError, match=r"^edit\.r: unknown keys \['branchng'\]"):
        build_spec("edit", "r", {**entry, "branchng": True})
    with pytest.raises(ValueError, match=r"^edit\.r: unknown keys \['dok'\]"):
        build_spec("edit", "r", {**entry, "params": [
            {"name": "x", "type": int, "default": 0, "dok": "x"}]})
    with pytest.raises(ValueError, match="unknown requirements"):
        build_spec("edit", "r", {**entry, "requires": ("mask",)})
    with pytest.raises(ValueError, match="unknown kind 'file'"):
        build_spec("edit", "r", {**entry, "params": [
            {"name": "x", "type": str, "default": "", "kind": "file"}]})
    with pytest.raises(ValueError, match="the default 'c' is not a choice"):
        build_spec("edit", "r", {**entry, "params": [
            {"name": "x", "type": str, "default": "c", "choices": ("a", "b")}]})
    with pytest.raises(ValueError, match="unknown provisions"):
        build_spec("edit", "r", {**entry, "provides": ("mask",)})


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


def test_the_registry_built_again_is_the_same():
    # The registry is built at import; built again here, every key of every
    # entry must give the same spec (mutation testing, phase 6-bis PR 4a:
    # the code that runs only at import was not checked by any test)
    rebuilt = registry._build()
    assert list(rebuilt) == list(RECIPES)
    for category in RECIPES:
        assert list(rebuilt[category]) == list(RECIPES[category])
        for name, spec in RECIPES[category].items():
            assert repr(rebuilt[category][name]) == repr(spec), f"{category}.{name}"
            assert ([repr(p.default) for p in rebuilt[category][name].params]
                    == [repr(p.default) for p in spec.params])


def test_build_spec_reads_every_key():
    entry = {"brief": "A recipe.", "details": "More.", "url": "page.html",
             "gui_hidden": True, "branching": True, "requires": ("z_em", "resolution"),
             "provides": ("continuum",),
             "params": [{"name": "x", "type": float, "default": 1.0, "doc": "The x.",
                         "gui_hidden": True, "optional": True, "keep": ("_current_",),
                         "min": 0.0, "max": 2.0, "gt": -1.0, "allow_nan": True,
                         "agents": False},
                        {"name": "s", "type": str, "default": "a", "choices": ("a", "b"),
                         "kind": "session"},
                        {"name": "y", "type": int, "default": 3}]}
    spec = build_spec("edit", "r", entry)
    assert (spec.category, spec.name, spec.brief, spec.details, spec.url) == \
        ("edit", "r", "A recipe.", "More.", "page.html")
    assert spec.gui_hidden and spec.branching
    assert spec.requires == frozenset({"z_em", "resolution"})
    assert spec.provides == frozenset({"continuum"})
    x, z, y = spec.params
    assert (z.choices, z.kind) == (("a", "b"), "session")
    assert (y.choices, y.kind) == ((), "")
    assert (x.name, x.type, x.default, x.doc, x.gui_hidden, x.optional, x.keep,
            x.min, x.max, x.gt, x.allow_nan, x.agents) == \
        ("x", float, 1.0, "The x.", True, True, ("_current_",), 0.0, 2.0, -1.0, True, False)
    # The defaults of the optional keys
    assert (y.doc, y.gui_hidden, y.optional, y.keep, y.min, y.max, y.gt, y.allow_nan,
            y.agents) == ("y", False, False, (), None, None, None, False, True)
    bare = build_spec("flux", "s", {"brief": "B."})
    assert (bare.details, bare.url, bare.params, bare.gui_hidden, bare.branching,
            bare.requires, bare.provides) == ("", None, (), False, False, frozenset(),
                                              frozenset())
    # The text of the methods (phase 6-ter, PR 5c): none, a generic sentence
    assert (bare.methods, bare.references) == (None, ())
    told = build_spec("flux", "s", {"brief": "B.", "methods": "Done{menu}",
                                    "references": ["C26"]})
    assert (told.methods, told.references) == ("Done{menu}", ("C26",))
