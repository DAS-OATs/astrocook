"""The old names of renamed parameters (phase 6-quater, PR 2).

A parameter renamed in the registry keeps its old names as ``aliases``: the
scripts, the logs saved in the archives (and their replay) and the MCP
tools that use an old name still run, with the same result; the log, the
methods text, the Log Scripter and the dialogs show the new name only. An
alias is silent and permanent; an old name given with the new one is an
error (decisions agreed with Guido).

No recipe of the package is renamed here: the fixture ``renamed_parameter``
(``conftest.py``) declares ``sigma_kms`` of ``flux.smooth`` as renamed from
``sigma_old``. The tests of Run All are in ``test_gui_smoke.py``, those of
the MCP tools in ``test_mcp_tools.py``.
"""
import inspect
import warnings

import astropy.units as au
import numpy as np
import pytest

from astrocook.session.session import Session, load_session_from_file
from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, HistoryLog, SpectrumData
from astrocook.recipes.errors import InvalidParameter
from astrocook.recipes.params import resolve_params
from astrocook.recipes.registry import all_recipes, build_spec, canonical_params, get_recipe
from astrocook.recipes.result import run

FLUX = au.erg / au.cm**2 / au.s / au.nm


def _session(name="main"):
    """A noisy power law, z_em = 2.5."""
    rng = np.random.default_rng(0)
    x = np.arange(400.0, 440.0, 0.01)
    cont = (x / 420.0) ** -1.5
    y = cont + 0.02 * cont * rng.standard_normal(len(x))
    spec = Spectrum(SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x - 0.005, au.nm),
        xmax=DataColumn(x + 0.005, au.nm), y=DataColumn(y, FLUX),
        dy=DataColumn(0.02 * cont, FLUX), aux_cols={}, z_em=2.5))
    return Session(name=name, gui=None, spec=spec)


def _same(a, b):
    return np.array_equal(a.spec.y.value, b.spec.y.value, equal_nan=True)


# --- The registry ------------------------------------------------------------

def _entry(*params):
    return {"brief": "R.", "params": list(params)}


def test_an_alias_is_read_with_its_parameter():
    spec = build_spec("flux", "r", _entry(
        {"name": "width", "type": float, "default": 1.0, "aliases": ["w", "wid"]},
        {"name": "n", "type": int, "default": 3}))
    width, n = spec.params
    assert width.aliases == ("w", "wid") and n.aliases == ()
    assert spec.renamed() == {"w": "width", "wid": "width"}
    assert spec.param("width") is width
    with pytest.raises(KeyError):
        spec.param("w")  # the old name is not a parameter


@pytest.mark.parametrize("params", [
    # The name of another parameter
    [{"name": "a", "type": float, "default": 1.0, "aliases": ["b"]},
     {"name": "b", "type": float, "default": 1.0}],
    # Its own name
    [{"name": "a", "type": float, "default": 1.0, "aliases": ["a"]}],
    # The old name of two parameters
    [{"name": "a", "type": float, "default": 1.0, "aliases": ["c"]},
     {"name": "b", "type": float, "default": 1.0, "aliases": ["c"]}],
    # A parameter declared twice
    [{"name": "a", "type": float, "default": 1.0},
     {"name": "a", "type": float, "default": 2.0}],
    # A string, not a list of names; not a name
    [{"name": "a", "type": float, "default": 1.0, "aliases": "old"}],
    [{"name": "a", "type": float, "default": 1.0, "aliases": ["not a name"]}],
])
def test_an_alias_that_is_not_one_name_of_one_parameter_is_refused(params):
    with pytest.raises(ValueError):
        build_spec("flux", "r", _entry(*params))


def test_canonical_renames_the_old_names_in_their_place():
    spec = build_spec("flux", "r", _entry(
        {"name": "a", "type": float, "default": 1.0},
        {"name": "width", "type": float, "default": 1.0, "aliases": ["w"]},
        {"name": "n", "type": int, "default": 3}))
    assert list(spec.canonical({"n": 1, "w": 2.0, "a": 3.0}).items()) == [
        ("n", 1), ("width", 2.0), ("a", 3.0)]
    # Nothing to rename, or not a parameter (e.g. alias_map): as given
    assert spec.canonical({"width": 2.0, "other": 1}) == {"width": 2.0, "other": 1}
    given = {"width": 2.0}
    assert spec.canonical(given) == given and spec.canonical(given) is not given


@pytest.mark.parametrize("params", [{"w": 2.0, "width": 2.0}, {"width": 2.0, "w": 3.0}])
def test_the_old_and_the_new_name_together_are_refused(params):
    spec = build_spec("flux", "r", _entry(
        {"name": "width", "type": float, "default": 1.0, "aliases": ["w"]}))
    with pytest.raises(InvalidParameter, match="w is the old name of width: give only one"):
        spec.canonical(params)


#: The parameters of the package renamed, as in the table of the renamed
#: parameters of ROADMAP.md: (recipe, new name) -> old names
RENAMED = {("export_ascii", "folder"): ("path",)}  # phase 6-quater, PR 7


def test_the_renamed_parameters_of_the_package():
    # A renaming adds its row here and in the table of ROADMAP.md
    assert {(s.name, p.name): p.aliases for s in all_recipes() for p in s.params
            if p.aliases} == RENAMED


def test_the_fixture_renames_one_parameter(renamed_parameter):
    assert get_recipe("flux", "smooth") is renamed_parameter
    assert renamed_parameter.renamed() == {"sigma_old": "sigma_kms"}


# --- The scripts ---------------------------------------------------------------

def test_a_script_with_the_old_name_gives_the_same_result(renamed_parameter):
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # silent: no warning
        old = _session().flux.smooth(sigma_old=50.0)
    assert _same(old, _session().flux.smooth(sigma_kms=50.0))
    assert not _same(old, _session().flux.smooth())
    # Strings, as from the dialogs and the logs
    assert _same(_session().flux.smooth(sigma_old="50.0"), old)


def test_the_old_name_is_still_checked(renamed_parameter):
    with pytest.raises(InvalidParameter, match="sigma_kms: -1.0 is not above 0"):
        _session().flux.smooth(sigma_old=-1.0)


def test_a_script_with_both_names_is_refused(renamed_parameter):
    with pytest.raises(InvalidParameter) as error:
        _session().flux.smooth(sigma_old=50.0, sigma_kms=50.0)
    assert error.value.recipe == "flux.smooth"
    assert "sigma_old is the old name of sigma_kms" in str(error.value)


def test_without_the_alias_the_old_name_is_refused():
    with pytest.raises(TypeError):
        _session().flux.smooth(sigma_old=50.0)


def test_the_result_and_the_log_have_the_new_name(renamed_parameter):
    result = run(_session(), "flux.smooth", sigma_old=50.0)
    assert result.params == {"sigma_kms": 50.0}
    assert resolve_params("smooth", {"sigma_old": "50.0"}) == {"sigma_kms": 50.0}
    assert _same(result.session, _session().flux.smooth(sigma_kms=50.0))


def test_the_signature_and_the_help_have_the_new_name(renamed_parameter):
    from astrocook.recipes.flux import RecipeFlux
    names = inspect.signature(RecipeFlux.smooth).parameters
    assert "sigma_kms" in names and "sigma_old" not in names


# --- The logs of the archives ---------------------------------------------------

def _archive_with_old_names(tmp_path):
    """An archive saved with the old name in its log (as by the version
    before the renaming)."""
    session = _session().flux.smooth(sigma_kms=50.0)
    log = HistoryLog()
    log.add_entry("smooth", {"sigma_old": 50.0}, reason="Too noisy.")
    log.add_entry("not_a_recipe", {"sigma_old": 1.0})
    session.log_manager = log
    path = tmp_path / "old.acs2"
    session.save(str(path))
    return path


def test_the_log_of_an_archive_is_read_with_the_new_names(renamed_parameter, tmp_path):
    path = _archive_with_old_names(tmp_path)
    loaded = load_session_from_file(str(path), "old", None, "auto")
    smooth, other = loaded.log_manager.entries
    assert smooth.params == {"sigma_kms": 50.0} and smooth.reason == "Too noisy."
    assert other.params == {"sigma_old": 1.0}  # not a recipe: as read
    # Its replay gives the state saved
    replayed = _session().flux.smooth(**smooth.params)
    assert _same(replayed, loaded)


def test_the_archive_itself_is_not_changed(renamed_parameter, tmp_path):
    import json
    import tarfile
    path = _archive_with_old_names(tmp_path)
    load_session_from_file(str(path), "old", None, "auto")
    with tarfile.open(path) as tar:
        member = next(m for m in tar.getmembers() if m.name.endswith(".json"))
        saved = json.load(tar.extractfile(member))["log_history_json"]
    assert saved["entries"][0]["params"] == {"sigma_old": 50.0}


def test_a_log_entry_with_both_names_is_kept_as_read(renamed_parameter):
    # Only by hand: the replay says what is wrong
    both = {"sigma_old": 1.0, "sigma_kms": 2.0}
    assert canonical_params("smooth", both) == both
    assert canonical_params("smooth", {"sigma_old": 1.0}) == {"sigma_kms": 1.0}


def test_the_methods_text_and_the_log_scripter_show_the_new_name(renamed_parameter,
                                                                  tmp_path):
    from astrocook.gui.log_script import log_to_script
    from astrocook.methods_text import methods_text
    loaded = load_session_from_file(str(_archive_with_old_names(tmp_path)), "old", None,
                                    "auto")
    (smooth,) = [p.text for p in methods_text(loaded).paragraphs if p.recipe == "flux.smooth"]
    assert "standard deviation 50 km/s" in smooth and "sigma" not in smooth
    lines = log_to_script(loaded.log_manager)
    assert "smooth(sigma_kms=50.0)" in lines
    assert not any("sigma_old" in line for line in lines if line.startswith("smooth"))


def test_the_dialog_shows_the_new_name(renamed_parameter, qtbot):
    from astrocook.gui.recipe_dialog import RecipeDialog
    dialog = RecipeDialog("flux", "smooth", _session())
    qtbot.addWidget(dialog)
    assert list(dialog.input_widgets) == ["sigma_kms"]
