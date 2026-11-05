"""The JSON schemas of the parameters of every recipe, from the registry
(phase 6-bis, PR 8): the test that the registry is enough to generate the
tools of the full MCP server (``ROADMAP.md``, section 4)."""
import math

import pytest

jsonschema = pytest.importorskip("jsonschema")

from astrocook.mcp_server.schemas import parameters_schema  # noqa: E402
from astrocook.recipes.params import convert  # noqa: E402
from astrocook.recipes.registry import all_recipes, get_recipe  # noqa: E402

SPECS = list(all_recipes())


@pytest.mark.parametrize("spec", SPECS, ids=lambda s: f"{s.category}.{s.name}")
def test_every_recipe_has_a_valid_schema(spec):
    schema = parameters_schema(spec)
    jsonschema.Draft202012Validator.check_schema(schema)
    offered = {p.name for p in spec.params if p.agents}
    assert set(schema["properties"]) == offered
    # The defaults are valid values of the schema, and of the recipe
    defaults = {name: prop["default"] for name, prop in schema["properties"].items()
                if "default" in prop}
    jsonschema.validate(defaults, {**schema, "required": []})
    for name, value in defaults.items():
        converted = convert(spec.param(name), value)
        assert converted == spec.param(name).default or (
            isinstance(converted, float) and math.isnan(converted))
    # Every value allowed is a value of the recipe
    for name, prop in schema["properties"].items():
        for value in prop.get("enum", []):
            assert convert(spec.param(name), value) == value


def test_what_a_schema_says():
    schema = parameters_schema(get_recipe("absorbers", "refit_all"))
    region = schema["properties"]["region_limit"]
    assert region["enum"] == ["None", "lya_forest", "all_ly_forest", "red_side"]
    assert region["default"] == "None"
    max_b = schema["properties"]["max_b"]
    assert (max_b["type"], max_b["exclusiveMinimum"], max_b["maximum"]) == ("number", 1, 300)
    # Not offered to agents: the options of debug and the questions of the GUI
    identify = parameters_schema(get_recipe("absorbers", "identify_lines"))["properties"]
    assert "debug_rating" not in identify and "bypass_scoring" not in identify
    assert identify["mask_col"]["x-astrocook-kind"] == "column"
    assert "renorm_model" not in parameters_schema(
        get_recipe("continuum", "estimate_auto"))["properties"]
    # Required: no default; NaN is not JSON (rebin.filling) and is left out
    resample = parameters_schema(get_recipe("flux", "resample"))
    assert resample["properties"]["target_session"]["x-astrocook-kind"] == "session"
    assert "default" not in parameters_schema(
        get_recipe("flux", "rebin"))["properties"]["filling"]
    assert parameters_schema(get_recipe("file", "export_ascii"))["required"] == [
        "path"]
    # The current value of the session: a number or '_current_'
    z_em = parameters_schema(get_recipe("edit", "set_properties"))["properties"]["z_em"]
    assert z_em["type"] == ["number", "string"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate({"region_limit": "nowhere"}, schema)
