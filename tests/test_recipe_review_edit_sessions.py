"""Review of the recipes of edit on several spectra (phase 6-quater, PR 9).

``delete``, ``split``, ``extract_preset``, ``trim_common``, ``stitch``,
``equalize_and_stitch``, ``coadd``, ``import_systems``,
``import_systems_ascii``: the behaviour kept and the behaviour changed,
each test saying what it was. The fingerprints of the valid results were
taken on master before the change (the numbers of a valid call do not move,
except where a finding says so: F132, F133, F134).
"""
import hashlib
import warnings
from types import SimpleNamespace

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, SpectrumData
from astrocook.core.units import ARBITRARY_FLUX
from astrocook.recipes.edit import RecipeEdit
from astrocook.recipes.errors import InvalidParameter, RecipeFailed
from astrocook.session.session import Session

PER_A = au.erg / (au.cm**2 * au.s * au.AA)
PER_NM = au.erg / (au.cm**2 * au.s * au.nm)


def _spec(x0, x1, step, yunit=PER_A, scale=1.0, xunit=au.nm, seed=0, z_em=2.0):
    """A smooth power law with noise and one absorption line, in ``xunit``
    (``x0``, ``x1``, ``step`` in nm), its flux times ``scale`` in ``yunit``."""
    rng = np.random.default_rng(seed)
    x = np.arange(x0, x1, step)
    cont = (x / 450.0) ** -1.5
    y = cont * (1 - 0.5 * np.exp(-0.5 * ((x - 0.5 * (x0 + x1)) / 0.1) ** 2))
    y = scale * (y + 0.01 * cont * rng.standard_normal(len(x)))
    dy = scale * 0.01 * cont
    f = (1 * au.nm).to_value(xunit, equivalencies=au.spectral())
    return Spectrum(SpectrumData(
        x=DataColumn(x * f, xunit), xmin=DataColumn((x - step / 2) * f, xunit),
        xmax=DataColumn((x + step / 2) * f, xunit), y=DataColumn(y, yunit),
        dy=DataColumn(dy, yunit), aux_cols={"cont": DataColumn(scale * cont, yunit)},
        z_em=z_em))


def _session(main, others, name="A", gui=True):
    """``main`` as the session ``name``, open in a stand-in of the GUI with
    ``others`` ({name: spectrum}), as the GUI and the MCP server hold them."""
    s = Session(name=name, gui=None, spec=main)
    if gui:
        histories = [SimpleNamespace(display_name=name, current_state=s)]
        histories += [SimpleNamespace(display_name=k,
                                      current_state=Session(name=k, gui=None, spec=v))
                      for k, v in others.items()]
        s._gui = SimpleNamespace(session_histories=histories)
    return s


def _fingerprint(session):
    spec = session.spec
    digest = hashlib.md5()
    for name in spec.t.colnames:
        values = np.asarray(spec.t[name], dtype=float)
        digest.update(name.encode())
        digest.update(np.round(np.nan_to_num(values, nan=-9e9), 9).tobytes())
    digest.update(f"{spec.x.unit}|{spec.y.unit}|{session.name}".encode())
    return digest.hexdigest()


@pytest.fixture
def arms():
    """Three arms in nm, the same unit, the redder ones miscalibrated."""
    blue = _spec(380.0, 460.0, 0.02, seed=1)
    mid = _spec(450.0, 560.0, 0.03, scale=0.8, seed=2)
    red = _spec(550.0, 700.0, 0.05, scale=1.3, seed=3)
    return _session(blue, {"mid": mid, "red": red})


def _quiet(run):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run()


# --- Kept: the numbers of the valid calls (fingerprints taken on master) ---

KEPT = {
    "stitch": lambda s: RecipeEdit(s).stitch(["mid", "red"]),
    "stitch_unsorted": lambda s: RecipeEdit(s).stitch(["red"], sort=False),
    "eqst_auto": lambda s: RecipeEdit(s).equalize_and_stitch("mid, red", "455, 555"),
    "eqst_ranges": lambda s: RecipeEdit(s).equalize_and_stitch(
        "mid, red", "455, 555", equalize_ranges="451-459, 551-559"),
    "eqst_narrow_range": lambda s: RecipeEdit(s).equalize_and_stitch(
        "mid, red", "455, 555", equalize_ranges="451-451.05, 551-559"),
    "eqst_factors": lambda s: RecipeEdit(s).equalize_and_stitch(
        "red, mid", "455, 555", manual_factors="1.25, 0.77"),
    "coadd_eq0": lambda s: RecipeEdit(s).coadd("A, mid"),
    "coadd_eq1": lambda s: RecipeEdit(s).coadd("A, mid", equalize_order=1),
    "coadd_none": lambda s: RecipeEdit(s).coadd("A, mid", equalize_order=-1),
    "coadd_nm": lambda s: RecipeEdit(s).coadd("A, mid, red", xstart=400, xend=600,
                                              dx=0.05, xunit="nm", kappa=None),
    # F138, added to PR 9: flux.equalize and the core of coadd
    "equalize_0": lambda s: s.flux.equalize("mid", order=0),
    "equalize_1": lambda s: s.flux.equalize("mid", order=1),
    "equalize_-1": lambda s: s.flux.equalize("mid", order=-1),
    "core_coadd": lambda s: s.with_new_spectrum(s.spec.coadd(
        [s._gui.session_histories[1].current_state.spec], None, None,
        au.Quantity(10.0, "km/s"))),
}


#: Taken on master before the change (3198345)
FINGERPRINTS = {
    "coadd_eq0": "f5660f5d7c8a4f2f5be591f325135feb",
    "coadd_eq1": "66ab85a4fd02a3f38c3a426320e328e2",
    "coadd_nm": "3b56e6f6b6f1460d93db93268541753e",
    "coadd_none": "79db46278518322b449d1f6ee5f95cfc",
    "eqst_auto": "070c9e43c5b5d05239aeab95877aa6ea",
    "eqst_factors": "f4856a51342ef4a2591aa9fb4c10f931",
    "eqst_narrow_range": "e76215238dd0f53677d8667c7ba57dae",
    "eqst_ranges": "f8f1e0b4e9b95e55bec773bac4da5fa8",
    "stitch": "c99259d1886fc72b5bba9a1c18167dbf",
    "stitch_unsorted": "03cec5a59dae2dfae99f306fdef6f688",
    "equalize_0": "e0ae22a5cc94dae585baba4afc612664",
    "equalize_1": "c47b53aa0af3e0a64dcbf724c1ba0a79",
    "equalize_-1": "2a287b990ba5d8da07b43608fb0e7520",
    "core_coadd": "f5660f5d7c8a4f2f5be591f325135feb",  # the same as coadd_eq0
}


@pytest.mark.parametrize("case", sorted(KEPT))
def test_valid_results_unchanged(arms, case):
    """The arms in nm and in one unit of the flux: the refactoring of
    equalize_and_stitch, the conversions of F132-F134 and the new checks
    leave these results as they were, to the bit."""
    result = _quiet(lambda: KEPT[case](arms))
    assert _fingerprint(result) == FINGERPRINTS[case]


# --- Changed: the units of the flux and of the axis (F132, F133, F134) ---

def _level(session, lo=None, hi=None):
    """The median flux over the continuum of the pixels between lo and hi
    (nm): 1 where the flux has the level of the reference arm."""
    spec = session.spec
    x = spec.x.to_value(au.nm)
    keep = np.ones(len(x), dtype=bool)
    if lo is not None:
        keep &= x > lo
    if hi is not None:
        keep &= x < hi
    cont = (x[keep] / 450.0) ** -1.5
    return float(np.nanmedian(spec.y.value[keep] / cont))


def test_equalize_and_stitch_converts_the_flux_before_scaling():
    # F132: an arm per nm (the same physical flux, 10 times the values) was
    # equalized on its values, then converted by the stitch: 10 times too
    # faint (0.1)
    s = _session(_spec(380.0, 460.0, 0.02), {"red": _spec(450.0, 560.0, 0.03, PER_NM, 10.0)})
    out = _quiet(lambda: RecipeEdit(s).equalize_and_stitch("red", "455"))
    assert out.spec.y.unit == PER_A
    assert _level(out, 470.0) == pytest.approx(1.0, abs=0.01)
    # A factor multiplies the flux in its unit, as before: 1 keeps it
    out = _quiet(lambda: RecipeEdit(s).equalize_and_stitch("red", "455", manual_factors="1"))
    assert _level(out, 470.0) == pytest.approx(1.0, abs=0.01)


def test_coadd_converts_the_flux():
    # F133: the values were co-added as they were: without scaling, an arm
    # per nm 10 times too bright (with the scaling, right by chance)
    s = _session(_spec(380.0, 460.0, 0.02), {"red": _spec(450.0, 560.0, 0.03, PER_NM, 10.0)})
    for order in (None, 0):
        out = _quiet(lambda o=order: RecipeEdit(s).coadd("A, red", equalize_order=o))
        assert out.spec.y.unit == PER_A
        assert _level(out, 470.0) == pytest.approx(1.0, abs=0.01), order


@pytest.mark.parametrize("recipe", ["coadd", "equalize_and_stitch"])
def test_an_arm_in_angstrom_is_equalized(recipe):
    # F134: the overlap was looked for in the values of the axis, nm against
    # Å: none, so no scaling (the arm stayed twice as bright); the cuts of
    # equalize_and_stitch were applied in Å (the whole arm kept, pixels twice)
    s = _session(_spec(380.0, 460.0, 0.02),
                 {"red": _spec(450.0, 560.0, 0.03, scale=2.0, xunit=au.AA)})
    run = (lambda: RecipeEdit(s).coadd("A, red")) if recipe == "coadd" else \
        (lambda: RecipeEdit(s).equalize_and_stitch("red", "455"))
    out = _quiet(run)
    assert out.spec.x.unit == au.nm
    assert _level(out, 470.0) == pytest.approx(1.0, abs=0.01)
    assert np.all(np.diff(out.spec.x.value) > 0)


def test_equalize_and_stitch_orders_the_arms_in_one_unit():
    # F134: the arms were ordered by the values of their axis: an arm in Å
    # (4500) came after one in nm (450) even when bluer
    blue = _spec(380.0, 460.0, 0.02, xunit=au.AA)
    s = _session(blue, {"red": _spec(450.0, 560.0, 0.03, scale=2.0)})
    out = _quiet(lambda: RecipeEdit(s).equalize_and_stitch("red", "455"))
    assert out.spec.x.unit == au.AA
    x = out.spec.x.to_value(au.nm)
    assert x.min() == pytest.approx(380.0) and x.max() == pytest.approx(560.0, abs=0.05)
    assert np.all(np.diff(x) > 0)
    assert _level(out, 470.0) == pytest.approx(1.0, abs=0.01)
    assert out.name == "Stitched_A_red"


def test_an_arbitrary_flux_takes_the_unit_of_the_reference_when_scaled():
    # Until PR 9: refused by equalize_and_stitch (in its last stitch), and
    # co-added by coadd with the values as they were in either case
    s = _session(_spec(380.0, 460.0, 0.02),
                 {"red": _spec(450.0, 560.0, 0.03, ARBITRARY_FLUX, 5.0)})
    out = _quiet(lambda: RecipeEdit(s).equalize_and_stitch("red", "455"))
    assert out.spec.y.unit == PER_A and _level(out, 470.0) == pytest.approx(1.0, abs=0.01)
    out = _quiet(lambda: RecipeEdit(s).coadd("A, red"))
    assert out.spec.y.unit == PER_A and _level(out, 470.0) == pytest.approx(1.0, abs=0.01)
    # Not scaled: refused
    with pytest.raises(RecipeFailed, match="arbitrary units.*cannot be combined"):
        RecipeEdit(s).coadd("A, red", equalize_order=None)
    with pytest.raises(RecipeFailed, match="arbitrary units.*cannot be combined"):
        RecipeEdit(s).stitch("red")


@pytest.mark.parametrize("recipe", ["stitch", "coadd", "equalize_and_stitch"])
def test_a_normalized_flux_is_not_combined_with_a_physical_one(recipe):
    # coadd co-added them with the values as they were (the others refused)
    s = _session(_spec(380.0, 460.0, 0.02),
                 {"red": _spec(450.0, 560.0, 0.03, au.dimensionless_unscaled)})
    runs = {"stitch": lambda: RecipeEdit(s).stitch("red"),
            "coadd": lambda: RecipeEdit(s).coadd("A, red"),
            "equalize_and_stitch": lambda: RecipeEdit(s).equalize_and_stitch("red", "455")}
    with pytest.raises(RecipeFailed, match="normalized"):
        _quiet(runs[recipe])


def test_an_axis_in_velocity_is_refused():
    red = _spec(450.0, 560.0, 0.03).x_convert(2.0, au.km / au.s)
    s = _session(_spec(380.0, 460.0, 0.02), {"red": red})
    for run in (lambda: RecipeEdit(s).stitch("red"),
                lambda: RecipeEdit(s).coadd("A, red"),
                lambda: RecipeEdit(s).equalize_and_stitch("red", "455")):
        with pytest.raises(RecipeFailed, match="x_convert"):
            _quiet(run)


# --- Changed: one rule for the sessions (F135, F136) ---

def test_without_a_gui_a_name_fails_with_a_hint():
    # F135: an AttributeError (no session_histories) for coadd,
    # equalize_and_stitch, import_systems, trim_common
    s = _session(_spec(380.0, 460.0, 0.02), {}, gui=False)
    runs = [lambda: RecipeEdit(s).coadd("A, other"),
            lambda: RecipeEdit(s).equalize_and_stitch("other", "455"),
            lambda: RecipeEdit(s).import_systems("other"),
            lambda: RecipeEdit(s).trim_common("other", 2.0, "Ly_a", "Ly_a"),
            lambda: RecipeEdit(s).stitch("other"),
            lambda: s.flux.equalize("other")]
    for run in runs:
        with pytest.raises(RecipeFailed, match="In a script, pass the session itself"):
            run()


def test_in_a_script_the_sessions_themselves(arms):
    # stitch took them already; now every recipe that reads other sessions
    blue = arms.spec
    mid = Session(name="mid", gui=None, spec=_spec(450.0, 560.0, 0.03, scale=0.8, seed=2))
    red = Session(name="red", gui=None, spec=_spec(550.0, 700.0, 0.05, scale=1.3, seed=3))
    s = Session(name="A", gui=None, spec=blue)
    assert _fingerprint(_quiet(lambda: RecipeEdit(s).stitch([mid, red]))) \
        == FINGERPRINTS["stitch"]
    assert _fingerprint(_quiet(lambda: RecipeEdit(s).equalize_and_stitch(
        [mid, red], "455, 555"))) == FINGERPRINTS["eqst_auto"]
    assert _fingerprint(_quiet(lambda: RecipeEdit(s).coadd([s, mid]))) \
        == FINGERPRINTS["coadd_eq0"]
    # This session by its name, also without a GUI
    assert _fingerprint(_quiet(lambda: RecipeEdit(s).coadd(["A", mid]))) \
        == FINGERPRINTS["coadd_eq0"]


def test_the_log_holds_names():
    from astrocook.recipes.params import resolve_params
    mid = Session(name="mid", gui=None, spec=_spec(450.0, 560.0, 0.03))
    assert resolve_params("stitch", {"other_sessions": [mid, "red"]})["other_sessions"] \
        == "mid, red"
    assert resolve_params("stitch", {"other_sessions": ["mid", "red"]})["other_sessions"] \
        == "mid, red"
    assert resolve_params("import_systems", {"source_session": mid})["source_session"] == "mid"


def test_the_old_logs_of_stitch_replay(arms):
    # The GUI passed a list, the logs held it (other_sessions was a list)
    for value in (["mid", "red"], "['mid', 'red']", "mid, red"):
        assert _fingerprint(_quiet(lambda v=value: RecipeEdit(arms).stitch(v))) \
            == FINGERPRINTS["stitch"]


@pytest.mark.parametrize("run", [
    lambda s: RecipeEdit(s).coadd("A, nope"),  # skipped without a word
    lambda s: RecipeEdit(s).stitch(["nope", "mid"]),  # skipped with a WARNING
    lambda s: RecipeEdit(s).equalize_and_stitch("nope, red", "455, 555"),  # failed so
    lambda s: RecipeEdit(s).trim_common("nope", 2.0, "Ly_a", "Ly_a"),  # "no range in common"
])
def test_a_name_not_found_fails_naming_it(arms, run):
    with pytest.raises(RecipeFailed, match="no open session named 'nope'"):
        _quiet(lambda: run(arms))


def test_trim_common_takes_its_old_name(arms):
    a = _quiet(lambda: arms.edit.trim_common(others_names="mid", z_target=2.7,
                                             trans_self="Ly_a", trans_others="Ly_a"))
    b = _quiet(lambda: arms.edit.trim_common(other_sessions="mid", z_target=2.7,
                                             trans_self="Ly_a", trans_others="Ly_a"))
    assert _fingerprint(a) == _fingerprint(b)


# --- Changed: refused before the recipe runs ---

@pytest.mark.parametrize("kwargs, param", [
    (dict(session_names="mid, red"), "session_names"),  # without this session
    (dict(session_names="A, mid", xstart=500, xend=450), "xstart"),  # 0 pixels (F124)
    (dict(session_names="A, mid", xstart=450, xend=450), "xstart"),
    (dict(session_names="A, mid", xunit="s"), "xunit"),  # a failed conversion
    (dict(session_names="A, mid", equalize_order=-2), "equalize_order"),  # no scaling
    (dict(session_names="A, mid", equalize_order=6), "equalize_order"),  # accepted
])
def test_coadd_refused(arms, kwargs, param):
    with pytest.raises(InvalidParameter, match=f"{param}: "):
        arms.edit.coadd(**kwargs)


def test_coadd_minus_one_is_none(arms):
    # F123: -1 (no scaling, as the old logs hold it) is None
    none = _quiet(lambda: arms.edit.coadd("A, mid", equalize_order=None))
    assert _fingerprint(none) == FINGERPRINTS["coadd_none"]


@pytest.mark.parametrize("kwargs, param", [
    (dict(stitch_wavelengths="555, 455"), "stitch_wavelengths"),  # an empty split
    (dict(stitch_wavelengths="455, 455"), "stitch_wavelengths"),
    (dict(stitch_wavelengths="455, 750"), "stitch_wavelengths"),  # outside red
    (dict(stitch_wavelengths="455, 555", manual_factors="0, 1"), "manual_factors"),  # taken
    (dict(stitch_wavelengths="455, 555", manual_factors="-1, 1"), "manual_factors"),
    (dict(stitch_wavelengths="455, 555", manual_factors="nan, 1"), "manual_factors"),
    (dict(stitch_wavelengths="455, 555", equalize_ranges="459-451, 551-559"),
     "equalize_ranges"),  # an empty split
])
def test_equalize_and_stitch_refused(arms, kwargs, param):
    with pytest.raises(InvalidParameter, match=f"{param}: "):
        _quiet(lambda: arms.edit.equalize_and_stitch("mid, red", **kwargs))


@pytest.mark.parametrize("run, param", [
    (lambda s: s.edit.delete(""), "targets"),  # the session as it was
    (lambda s: s.edit.delete("nope"), "targets"),  # a WARNING, nothing done
    (lambda s: s.edit.delete("dF"), "targets"),  # a WARNING, nothing done
    (lambda s: s.edit.split(""), "expression"),  # RecipeFailed
    (lambda s: s.edit.extract_preset("lya_forest, nope"), "region"),  # nope skipped
    (lambda s: s.edit.trim_common("mid", 2.7, "Nope", "Ly_a"), "trans_self"),  # no range
    (lambda s: s.edit.trim_common("mid", 2.7, "Ly_a", "Nope"), "trans_others"),  # ignored
    (lambda s: s.edit.trim_common("mid, red", 2.7, "Ly_a", "Ly_a"), "trans_others"),
])
def test_refused(arms, run, param):
    with pytest.raises(InvalidParameter, match=f"{param}"):
        run(arms)


def test_delete_still_deletes(arms):
    out = arms.edit.delete("cont")
    assert not out.spec.has_aux_column("cont")


def test_required_parameters():
    from astrocook.mcp_server.schemas import parameters_schema
    from astrocook.recipes.registry import get_recipe
    split = parameters_schema(get_recipe("edit", "split"))
    assert "expression" in split["required"]
    imp = parameters_schema(get_recipe("edit", "import_systems"))
    assert "source_session" in imp["required"]  # it was 'None'
    stitch = parameters_schema(get_recipe("edit", "stitch"))
    assert stitch["properties"]["other_sessions"]["type"] == "string"  # an array
    assert stitch["properties"]["other_sessions"]["x-astrocook-kind"] == "sessions"
    coadd = parameters_schema(get_recipe("edit", "coadd"))["properties"]["equalize_order"]
    assert coadd["maximum"] == 5


# --- Changed: F138, flux.equalize and the core (added to PR 9) ---

def _to_equalize(main, ref):
    """``main`` as the session A, to be equalized to the session ref."""
    return _session(main, {"ref": ref})


def test_equalize_compares_the_flux_in_one_unit():
    # F138: a session per nm (10 times the values of the same flux)
    # equalized to one per Å kept its unit with the values of the other:
    # 10 times too faint
    s = _to_equalize(_spec(450.0, 560.0, 0.03, PER_NM, 10.0), _spec(380.0, 460.0, 0.02))
    out = _quiet(lambda: s.flux.equalize("ref"))
    assert out.spec.y.unit == PER_NM
    level = _level(out) * PER_NM
    assert level.to_value(PER_A) == pytest.approx(1.0, abs=0.01)


def test_equalize_compares_the_axis_in_one_unit():
    # F138: an axis in Å against one in nm: no overlap of the values, so
    # no scaling (the flux stayed twice as bright)
    s = _to_equalize(_spec(450.0, 560.0, 0.03, scale=2.0, xunit=au.AA),
                     _spec(380.0, 460.0, 0.02))
    out = _quiet(lambda: s.flux.equalize("ref"))
    assert out.spec.x.unit == au.AA
    assert _level(out) == pytest.approx(1.0, abs=0.01)


def test_equalize_an_arbitrary_flux_takes_the_unit_of_the_reference():
    # The values were scaled and the unit kept arbitrary
    s = _to_equalize(_spec(450.0, 560.0, 0.03, ARBITRARY_FLUX, 5.0), _spec(380.0, 460.0, 0.02))
    out = _quiet(lambda: s.flux.equalize("ref"))
    assert out.spec.y.unit == PER_A and _level(out) == pytest.approx(1.0, abs=0.01)


@pytest.mark.parametrize("main, match", [
    (_spec(450.0, 560.0, 0.03, au.dimensionless_unscaled), "normalized"),  # scaled
    (_spec(450.0, 560.0, 0.03).x_convert(2.0, au.km / au.s), "x_convert"),  # not scaled
])
def test_equalize_refused(main, match):
    s = _to_equalize(main, _spec(380.0, 460.0, 0.02))
    with pytest.raises(RecipeFailed, match=match):
        _quiet(lambda: s.flux.equalize("ref"))


def test_the_core_in_a_script():
    # Spectrum.coadd and Spectrum.equalize_to_reference, called by a
    # script, compare in the units of the reference too (F133, F134, F138)
    blue = _spec(380.0, 460.0, 0.02)
    red = _spec(450.0, 560.0, 0.03, PER_NM, 10.0)
    for order in (-1, 0):
        out = _quiet(lambda o=order: blue.coadd([red], None, None,
                                                au.Quantity(10.0, "km/s"), equalize_order=o))
        assert _level(Session(name="c", gui=None, spec=out), 470.0) \
            == pytest.approx(1.0, abs=0.01), order
    red_aa = _spec(450.0, 560.0, 0.03, scale=2.0, xunit=au.AA)
    out = _quiet(lambda: red_aa.equalize_to_reference(blue))
    assert _level(Session(name="e", gui=None, spec=out)) == pytest.approx(1.0, abs=0.01)
    with pytest.raises(ValueError, match="x_convert"):
        _spec(450.0, 560.0, 0.03).x_convert(2.0, au.km / au.s).equalize_to_reference(blue)
