"""The recipes log the same message on parameters they cannot parse.

Characterization test for phase 3: the message comes from
``astrocook.legacy.message`` in V1 and must not change when it moves into
the V2 code. One recipe per module, each with a non-numeric parameter.
"""
import logging

import pytest

from test_recipes_edit_flux import _session

MSG = "I can't understand the values of the parameters. Please try again."


@pytest.mark.parametrize("recipe, kwargs", [
    ("flux.calculate_running_std", dict(window_pix="abc")),
    ("edit.set_properties", dict(z_em="abc")),
    ("continuum.find_absorbed", dict(smooth_len_lya="abc")),
    ("absorbers.identify_lines", dict(min_pix_region="abc")),
])
def test_bad_parameter_logs_message_and_returns_0(recipe, kwargs, caplog):
    group, name = recipe.split(".")
    method = getattr(getattr(_session(), group), name)
    with caplog.at_level(logging.ERROR):
        assert method(**kwargs) == 0
    assert MSG in [r.getMessage() for r in caplog.records
                   if r.levelno == logging.ERROR]
