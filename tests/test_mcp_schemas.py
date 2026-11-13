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
    # F106: the prefix of the file names is required, as the folder
    # ('folder' was 'path' until phase 6-quater, PR 7)
    assert parameters_schema(get_recipe("file", "export_ascii"))["required"] == [
        "folder", "prefix"]
    # The current value of the session: a number or '_current_'
    z_em = parameters_schema(get_recipe("edit", "set_properties"))["properties"]["z_em"]
    assert z_em["type"] == ["number", "string"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate({"region_limit": "nowhere"}, schema)


# --- The descriptions paid in every tool (phase 6-quater, PR 12) -------------------------------

#: The longest description of a shared parameter, and of any parameter
SHARED_MAX = 160
ANY_MAX = 300
#: The total of the shared parameters over the tools (21 360 characters in the 48 recipe
#: tools when they were shortened, from 32 352)
SHARED_TOTAL_MAX = 24_000


def _tools():
    import asyncio
    pytest.importorskip("mcp")
    from astrocook.mcp_server.server import build_server
    return asyncio.run(build_server().list_tools())


def test_the_shared_parameters_have_one_text_and_are_short():
    from astrocook.mcp_server.schemas import SHARED_PARAMS
    tools = _tools()
    total = 0
    for tool in tools:
        props = tool.input_schema["properties"]
        # One text in every tool that has the parameter, the one of the
        # registry of shared parameters (the reason and the maximum listed of
        # session_undo, session_history and the like are their own)
        if "session_id" in props:
            assert props["session_id"]["description"] == \
                SHARED_PARAMS["session_id"]["description"], tool.name
        if "provide_missing" in props:  # a recipe tool
            for name, shared in SHARED_PARAMS.items():
                assert props[name] == {**shared}, (tool.name, name)
                total += len(props[name]["description"])
        for name, prop in props.items():
            assert len(prop.get("description", "")) <= ANY_MAX, (tool.name, name)
    for name, shared in SHARED_PARAMS.items():
        assert len(shared["description"]) <= SHARED_MAX, name
    assert 0 < total <= SHARED_TOTAL_MAX, total
