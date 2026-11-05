"""The recipes report the same way the parameters they cannot read.

Characterization test for phase 3: the message came from
``astrocook.legacy.message`` in V1 and was kept when it moved into the V2
code; phase 6 rewrote it by the rule of the log messages (``ROADMAP.md``);
phase 6-bis (PR 3) turned it into an InvalidParameter that names the recipe
and the parameter. One recipe per module, each with a non-numeric parameter.
"""
import logging

import pytest

from astrocook.recipes.errors import InvalidParameter
from test_recipes_edit_flux import _session

@pytest.mark.parametrize("recipe, kwargs, title", [
    ("flux.calculate_running_std", dict(window_pix="abc"),
     "Calculate running standard deviation"),
    ("edit.set_properties", dict(z_em="abc"), "Set session properties"),
    ("continuum.find_absorbed", dict(smooth_len_lya="abc"), "Find absorbed regions (V1 logic)"),
    ("absorbers.identify_lines", dict(min_pix_region="abc"),
     "Identify likely absorption systems"),
])
def test_bad_parameter_raises_and_names_it(recipe, kwargs, title, caplog):
    # Since phase 6-bis (PR 3) an exception, logged where it is caught: the
    # message names the recipe, the parameter and its value
    group, name = recipe.split(".")
    method = getattr(getattr(_session(), group), name)
    (param, value), = kwargs.items()
    with caplog.at_level(logging.ERROR), pytest.raises(InvalidParameter) as info:
        method(**kwargs)
    assert str(info.value).startswith(f"{title}: {param}: {value!r} is not")
    assert info.value.recipe == recipe
    assert caplog.text == ""
