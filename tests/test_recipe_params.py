"""Typed parameters of the recipes (phase 6-bis, PR 2).

The parameters are converted once, by ``typed_recipe``, to the types of the
registry: a recipe gives the same result with the strings of the dialogs
and of the session log as with typed values, and the log records the values
actually used (#573, item 3). The conversion itself; what each recipe does
with each kind of input is in ``test_recipe_inputs*.py``.
"""
import logging
import math

import numpy as np
import pytest

import astropy.units as au

from astrocook.core.session import SessionV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, HistoryLogV2, SpectrumDataV2
from astrocook.recipes.params import (InvalidParameter, convert,
                                      resolve_params)
from astrocook.recipes.registry import ParamSpec, get_recipe

FLUX = au.erg / au.cm**2 / au.s / au.nm


def _session():
    """A noisy power law with a Gaussian absorption line, z_em = 2.5."""
    rng = np.random.default_rng(0)
    x = np.arange(380.0, 460.0, 0.01)
    cont = (x / 420.0) ** -1.5
    y = cont * (1 - 0.5 * np.exp(-0.5 * ((x - 420.0) / 0.05) ** 2))
    y = y + 0.02 * cont * rng.standard_normal(len(x))
    spec = SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - 0.005, au.nm),
        xmax=DataColumnV2(x + 0.005, au.nm), y=DataColumnV2(y, FLUX),
        dy=DataColumnV2(0.02 * cont, FLUX), aux_cols={}, z_em=2.5))
    return SessionV2(name="main", gui=None, spec=spec)


def _p(type_, **kw):
    return ParamSpec(name="x", type=type_, doc="x", **kw)


# --- The conversion -----------------------------------------------------

@pytest.mark.parametrize("type_, value, expected", [
    (float, "2.5", 2.5), (float, 2, 2.0), (float, " 3 ", 3.0), (float, 1e-3, 1e-3),
    (int, "21", 21), (int, 21, 21), (int, 3.0, 3),
    (bool, "True", True), (bool, "False", False), (bool, "true", True), (bool, False, False),
    (str, "auto", "auto"), (str, 1.0, "1.0"),
    (list, [1, 2], [1, 2]), (list, (1, 2), [1, 2]), (list, "[1.0, 2.0]", [1.0, 2.0]),
])
def test_values_are_converted_to_their_type(type_, value, expected):
    result = convert(_p(type_), value)
    assert result == expected and type(result) is type(expected)


@pytest.mark.parametrize("type_, value, message", [
    (float, "abc", "'abc' is not a number"),
    (float, True, "True is not a number"),
    (float, "", "'' is not a number"),
    (float, "nan", "'nan' is not a number"),
    (float, None, "a value is needed"),
    (int, "3.5", "'3.5' is not an integer"),
    (int, 3.5, "3.5 is not an integer"),
    (bool, "yes", "'yes' is not True or False"),
    (bool, 2, "2 is not True or False"),
    (list, "abc", "'abc' is not a list"),
])
def test_values_that_cannot_be_converted_are_refused(type_, value, message):
    with pytest.raises(InvalidParameter, match=f"x: {message}"):
        convert(_p(type_), value)


def test_optional_none_nan_kept_and_bounds():
    optional = _p(float, optional=True)
    assert convert(optional, None) is None
    assert convert(optional, "None") is None
    assert convert(optional, "") is None  # an empty field, default None
    with pytest.raises(InvalidParameter, match="x: '' is not a number"):
        # an empty field where the default is not None: an error (Guido)
        convert(ParamSpec(name="x", type=float, doc="x", optional=True, _default=5.0), "")
    assert convert(_p(bool), 1) is True and convert(_p(bool), "0") is False
    assert convert(optional, "4") == 4.0
    assert math.isnan(convert(_p(float, allow_nan=True), "nan"))
    assert convert(_p(float, keep=("_current_",)), "_current_") == "_current_"
    assert convert(_p(str, optional=True), "None") is None
    assert convert(_p(str, optional=True), None) is None  # PR 2c
    assert convert(_p(str), "None") == "None"  # a name, e.g. of a session
    bounded = _p(float, min=0.0, max=1.0)
    assert convert(bounded, "0.5") == 0.5
    with pytest.raises(InvalidParameter, match="x: -1.0 is below 0.0"):
        convert(bounded, -1)
    with pytest.raises(InvalidParameter, match="x: 2.0 is above 1.0"):
        convert(bounded, 2)


def test_registry_defaults_convert_to_themselves():
    """Every default of a typed recipe is a valid value of its parameter."""
    from astrocook.recipes.registry import all_recipes
    for spec in all_recipes():
        for p in spec.params:
            if p.default in ("", None) and not p.optional:
                continue  # required by the API, pre-filled by the dialog
            value = convert(p, p.default)
            assert value == p.default or (isinstance(value, float) and math.isnan(value)), \
                (spec.name, p.name)


# --- The recipes ------------------------------------------------------------

def test_every_recipe_is_typed():
    """Every recipe of the registry converts its parameters (phase 6-bis,
    PR 2a-2c): a new recipe must be decorated with ``typed_recipe``."""
    from astrocook.recipes.registry import all_recipes
    from astrocook.core.session import SessionV2
    for spec in all_recipes():
        method = getattr(type(getattr(SessionV2("s", None), spec.category)), spec.name)
        assert getattr(method, "typed_recipe", False), f"{spec.category}.{spec.name}"


def test_strings_and_typed_values_give_the_same_result():
    session = _session()
    strings = session.flux.rebin(xstart="None", xend="None", dx="20.0", xunit="km/s",
                                 kappa="None", filling="nan")
    typed = session.flux.rebin(xstart=None, xend=None, dx=20.0, kappa=None,
                               filling=np.nan)
    assert np.array_equal(strings.spec.y.value, typed.spec.y.value, equal_nan=True)
    smooth_s = session.flux.smooth(sigma_kms="50.0")
    smooth_t = session.flux.smooth(sigma_kms=50)
    assert np.array_equal(smooth_s.spec.y.value, smooth_t.spec.y.value)
    std_s = session.flux.calculate_running_std(window_pix="20")
    std_t = session.flux.calculate_running_std(window_pix=20)
    assert np.array_equal(std_s.spec.t["running_std"], std_t.spec.t["running_std"])


def test_a_value_that_cannot_be_read_stops_the_recipe(caplog):
    session = _session()
    with caplog.at_level(logging.ERROR):
        assert session.flux.smooth(sigma_kms="fifty") == 0
    assert "The values of the parameters cannot be read" in caplog.text
    assert "sigma_kms: 'fifty' is not a number" in caplog.text


def test_an_unknown_unit_stops_rebin(caplog):
    with caplog.at_level(logging.ERROR):
        assert _session().flux.rebin(xunit="furlongs") == 0
    assert "xunit: 'furlongs' is not a unit" in caplog.text


def test_parameters_outside_the_registry_pass_unchanged(caplog):
    # smooth_std: deprecated, not in the registry (F16); a string as before
    session = _session().continuum.estimate_auto()
    assert session != 0


# --- The values in the log -------------------------------------------------

def test_resolve_params_gives_the_values_used():
    resolved = resolve_params("rebin", {"dx": "20", "kappa": "None"})
    assert list(resolved) == [p.name for p in get_recipe("flux", "rebin").params]
    assert resolved["dx"] == 20.0 and resolved["kappa"] is None
    assert resolved["xstart"] is None and math.isnan(resolved["filling"])
    assert resolved["xunit"] == "km/s"
    # Passed but not declared: kept as passed
    assert resolve_params("rebin", {"norm": "False"})["norm"] == "False"
    # A name that is not a recipe: as passed
    assert resolve_params("no_such_recipe", {"max_b": "50"}) == {"max_b": "50"}
    assert resolve_params("doublets_auto", {"max_b": "50"})["max_b"] == 50.0


def test_the_gui_logs_the_values_used():
    from astrocook.gui.main_recipes import _add_log_entry

    class History:
        log_manager = HistoryLogV2()

    history = History()
    _add_log_entry(history, "smooth", {"sigma_kms": "50.0"})
    (entry,) = history.log_manager.entries
    assert entry.params == {"sigma_kms": 50.0}
    _add_log_entry(history, "estimate_auto", {"kappa": "3"})
    entry = history.log_manager.entries[-1]
    assert entry.params["kappa"] == 3.0 and entry.params["template"] is False
    assert entry.params["smooth_len_lya"] == 5000.0  # a default, written out


def test_a_log_of_strings_still_replays(tmp_path):
    """The logs written before the typed parameters (strings, as by the
    dialogs of 2.4.0) replay with the same results: saved in an archive,
    read back and run again."""
    from astrocook.core.session import load_session_from_file
    session = _session()
    log = HistoryLogV2()
    log.add_entry("smooth", {"sigma_kms": "50.0"})
    log.add_entry("rebin", {"xstart": "None", "xend": "None", "dx": "20.0",
                            "xunit": "km/s", "kappa": "5.0", "filling": "nan"})
    session.log_manager = log
    path = tmp_path / "old_log.acs2"
    session.save(str(path))
    loaded = load_session_from_file(str(path), "old", None, "auto")
    entries = loaded.log_manager.to_dict()["entries"]
    assert [e["params"]["dx"] for e in entries if e["recipe_name"] == "rebin"] == ["20.0"]
    replayed = _session()
    for e in entries:
        replayed = getattr(replayed.flux, e["recipe_name"])(**e["params"])
    expected = _session().flux.smooth(sigma_kms=50.0).flux.rebin(dx=20.0)
    assert np.array_equal(replayed.spec.y.value, expected.spec.y.value, equal_nan=True)
