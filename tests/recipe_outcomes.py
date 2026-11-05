"""The outcomes of a recipe that does not give a result, as the batteries of
``test_recipe_inputs*.py`` record them (phase 6-bis, PR 3).

- REFUSED: a parameter cannot be used; the recipe raises InvalidParameter,
  whose message names the recipe and the parameter.
- FAILS: the value is converted, then the recipe cannot give a result; it
  raises another RecipeError, whose message names the recipe.
- RAISES: another exception leaves the recipe (a bug, or an error of a
  library that the boundary does not convert); it says which recipe in a
  note.
"""
import pytest

from astrocook.recipes.errors import InvalidParameter, RecipeError
from astrocook.recipes.registry import find_recipe

REFUSED = "refused"
FAILS = "fails"
RAISES = "raises"
FAILURES = (REFUSED, FAILS, RAISES)


def check_failure(run, name, param, expected):
    """Run ``run()`` (the recipe ``name`` with a value of ``param``) and
    check its failure ``expected``, one of :data:`FAILURES`."""
    title = find_recipe(name).brief.rstrip('.')
    if expected == REFUSED:
        with pytest.raises(InvalidParameter) as info:
            run()
        assert f"{param}: " in str(info.value), str(info.value)
    elif expected == FAILS:
        with pytest.raises(RecipeError) as info:
            run()
        assert not isinstance(info.value, InvalidParameter), str(info.value)
    else:
        with pytest.raises(Exception) as info:
            run()
        assert not isinstance(info.value, RecipeError), str(info.value)
        assert any(name in note for note in getattr(info.value, "__notes__", [])), \
            "the recipe is not named"
        return
    # The message names the recipe, and the recipe is recorded
    assert title in str(info.value), str(info.value)
    assert info.value.chain and info.value.chain[0][1] == name
