"""Characterization tests of the recipe layer of ``recipes/edit.py`` and
``recipes/flux.py``.

The SpectrumV2 methods behind several of these recipes are covered by the
snapshots; here the target is the recipe itself: parameter parsing (string
inputs as from the GUI dialogs), the lookup of other sessions by display
name, the chaining with the core methods, and what is returned on bad
input (a session, the same session, or 0). Numerical outputs not covered
elsewhere are stored as snapshots.

These tests record the current behaviour, including two questionable
points listed in the findings register (``ROADMAP.md``, section 3.2): F22
(some operations drop the global resolution attribute) and F23
(``flux.rebin`` raises on a non-numeric parameter instead of returning 0).
"""
from types import SimpleNamespace

import numpy as np
import astropy.units as au
import pytest

from astrocook.core.session import SessionV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import (
    ComponentDataV2, DataColumnV2, SpectrumDataV2, SystemListDataV2,
)
from astrocook.core.system_list import SystemListV2

FLUX = au.erg / (au.cm**2 * au.s * au.nm)
Z_EM = 4.1
C_KMS = 299792.458
N = len(np.arange(612.0, 620.0, 0.005))  # pixels of the default _spec()


def _spec(x0=612.0, x1=620.0, step=0.005, scale=1.0, seed=0, resol=0.0,
          z_em=Z_EM, extra=None):
    """Noisy power-law spectrum with ``cont`` and ``model`` columns."""
    rng = np.random.default_rng(seed)
    x = np.arange(x0, x1, step)
    n = len(x)
    cont = scale * (x / 616.0) ** -1.5
    model = cont * (1.0 - 0.5 * np.exp(-0.5 * ((x - 616.0) / 0.02) ** 2))
    y = model + 0.02 * cont * rng.standard_normal(n)
    aux = {"cont": DataColumnV2(cont, FLUX), "model": DataColumnV2(model, FLUX)}
    aux.update(extra or {})
    return SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - step / 2, au.nm),
        xmax=DataColumnV2(x + step / 2, au.nm), y=DataColumnV2(y, FLUX),
        dy=DataColumnV2(0.02 * cont, FLUX), aux_cols=aux, z_em=z_em,
        resol=resol))


def _systs(*zs, series="CIV"):
    comps = [ComponentDataV2(id=i + 1, z=z, dz=None, logN=13.5, dlogN=None,
                             b=10.0, db=None, series=series)
             for i, z in enumerate(zs)]
    return SystemListV2(SystemListDataV2(components=comps))


def _gui(**sessions):
    """Minimal stand-in for the GUI: the recipes look up other sessions by
    display name among the open session histories."""
    return SimpleNamespace(session_histories=[
        SimpleNamespace(display_name=name, current_state=state)
        for name, state in sessions.items()])


def _session(spec=None, gui=None, systs=None, name="main"):
    return SessionV2(name=name, gui=gui, spec=spec if spec is not None else _spec(),
                     systs=systs)


# --- Edit > Set properties ------------------------------------------------

def test_set_properties_numeric():
    session = _session(_spec(extra={"resol": DataColumnV2(
        np.full(N, 30000.0), au.dimensionless_unscaled)}))
    out = session.edit.set_properties(name="new", object="QSO", z_em="3.5",
                                      resol="45000")
    spec = out.spec
    assert out.name == "new"
    assert spec._data.z_em == 3.5
    assert spec._data.resol == 45000.0
    assert spec.meta["OBJECT"] == "QSO"
    assert spec.meta["resol"] == 45000.0
    assert spec.meta["resol_fwhm"] == pytest.approx(C_KMS / 45000.0)
    # The resolution column is overwritten with the new value
    np.testing.assert_array_equal(spec.get_column("resol").value, 45000.0)
    # The input session is unchanged
    assert session.spec._data.z_em == Z_EM


def test_set_properties_current_keeps_values():
    session = _session(_spec(resol=20000.0))
    out = session.edit.set_properties()
    assert out.name == "main"
    assert out.spec._data.z_em == Z_EM
    assert out.spec._data.resol == 20000.0


def test_set_properties_pixel_fwhm():
    # "Npx": R = median(lambda) / (N * median(d lambda))
    session = _session()
    out = session.edit.set_properties(resol="3px")
    x = session.spec.x.value
    expected = np.median(x) / (3 * np.median(np.diff(x)))
    assert out.spec._data.resol == pytest.approx(expected)


def test_set_properties_keeps_column_if_asked():
    column = np.linspace(40000.0, 50000.0, N)
    session = _session(_spec(extra={"resol": DataColumnV2(
        column, au.dimensionless_unscaled)}))
    out = session.edit.set_properties(resol="45000", overwrite_resol_col=False)
    np.testing.assert_array_equal(out.spec.get_column("resol").value, column)
    assert out.spec._data.resol == 45000.0


@pytest.mark.parametrize("kwargs", [dict(z_em="abc"), dict(resol="abc"),
                                    dict(resol="xpx")])
def test_set_properties_bad_input_returns_0(kwargs):
    assert _session().edit.set_properties(**kwargs) == 0


# --- Edit > Apply / Mask expression ---------------------------------------

def test_apply_expression_ui_names():
    session = _session()
    out = session.edit.apply_expression("F", "F * 2")
    np.testing.assert_allclose(out.spec.y.value, 2 * session.spec.y.value)
    assert out.spec.y.unit == FLUX  # the target keeps its unit
    snr = session.edit.apply_expression("snr", "F / dF").spec.get_column("snr")
    np.testing.assert_allclose(snr.value, session.spec.y.value / session.spec.dy.value)
    assert snr.unit == au.dimensionless_unscaled  # a new column has no unit


def test_apply_expression_other_session():
    # "s1.F": the other session is resampled on this grid
    other = _spec(scale=2.0, seed=1)
    session = _session(gui=_gui(other=SimpleNamespace(spec=other)))
    out = session.edit.apply_expression("diff", "F - s1.F",
                                        alias_map={"s1": "other"})
    np.testing.assert_allclose(out.spec.get_column("diff").value,
                               session.spec.y.value - other.y.value)


@pytest.mark.parametrize("expression", ["", "F +", "unknown * 2"])
def test_apply_expression_bad_input_returns_0(expression):
    assert _session().edit.apply_expression("F", expression) == 0


def test_mask_expression():
    session = _session()
    out = session.edit.mask_expression("F", "λ < 614")
    masked = session.spec.x.value < 614
    assert np.all(np.isnan(out.spec.y.value[masked]))
    np.testing.assert_array_equal(out.spec.y.value[~masked],
                                  session.spec.y.value[~masked])


@pytest.mark.parametrize("target, expression", [
    ("missing", "λ < 614"), ("F", "F * 2"), ("F", "")])
def test_mask_expression_bad_input_returns_0(target, expression):
    assert _session().edit.mask_expression(target, expression) == 0


@pytest.mark.parametrize("operation", [
    lambda s: s.edit.apply_expression("F", "F * 1"),
    lambda s: s.edit.mask_expression("F", "λ < 0"),
    lambda s: s.edit.smooth_column("dF", "10.0"),
    lambda s: s.edit.split("λ > 0"),
])
def test_some_operations_drop_the_resolution_attribute(operation):
    # Characterization of F22: these operations rebuild SpectrumDataV2
    # without the resol attribute (reset to 0). The Voigt model changes
    # only when the resolution is in the attribute alone: set_properties
    # and the loaders also store it in the metadata, which is kept.
    out = operation(_session(_spec(resol=45000.0)))
    assert out.spec._data.resol == 0.0


# --- Edit > Smooth column, Delete -----------------------------------------

def test_smooth_column(snapshot):
    session = _session()
    out = session.edit.smooth_column("F", "10.0")
    np.testing.assert_array_equal(out.spec.dy.value, session.spec.dy.value)
    snapshot.check("recipe_smooth_column_F_10kms", {"y": out.spec.y})


def test_smooth_column_cont_renormalises_model():
    session = _session()
    out = session.edit.smooth_column("cont", "2000.0", renorm_model="True")
    old = session.spec
    new = out.spec
    assert not np.allclose(new.get_column("cont").value, old.get_column("cont").value)
    # The normalised model profile is kept
    np.testing.assert_allclose(
        new.get_column("model").value / new.get_column("cont").value,
        old.get_column("model").value / old.get_column("cont").value)
    # Without the flag the model is left as it is
    out = session.edit.smooth_column("cont", "2000.0")
    np.testing.assert_array_equal(out.spec.get_column("model").value,
                                  old.get_column("model").value)


@pytest.mark.parametrize("target, sigma", [("missing", "10.0"), ("F", "abc")])
def test_smooth_column_bad_input_returns_0(target, sigma):
    assert _session().edit.smooth_column(target, sigma) == 0


def test_delete_columns_and_systems():
    session = _session(systs=_systs(2.97))
    out = session.edit.delete("cont, systems")
    assert not out.spec.has_aux_column("cont")
    assert out.spec.has_aux_column("model")
    assert out.systs.components == []
    assert session.edit.delete("") is session


# --- Edit > Split, Extract preset, Trim common ----------------------------

def test_split_keeps_rows_and_systems_in_range():
    # CIV 1548 at z = 2.97 falls at 614.7 nm, at z = 2.99 at 617.8 nm
    session = _session(systs=_systs(2.97, 2.99))
    out = session.edit.split("(λ > 612.5) & (λ < 616)")
    x = out.spec.x.value
    assert x.min() > 612.5 and x.max() < 616
    assert [c.z for c in out.systs.components] == [2.97]


def test_split_empty_expression_returns_0():
    assert _session().edit.split("") == 0


def test_extract_preset():
    # z_em = 4.1: Lyman-beta at 523.1 nm, Lyman-alpha at 620.0 nm
    session = _session(_spec(500.0, 700.0, 0.05))
    forest = session.edit.extract_preset("lya_forest")
    x = forest.spec.x.value
    assert x.min() == pytest.approx(102.5722 * (1 + Z_EM), abs=0.05)
    assert x.max() == pytest.approx(121.567 * (1 + Z_EM), abs=0.05)
    both = session.edit.extract_preset("lya_forest, red_side")
    assert isinstance(both, list) and len(both) == 2
    assert both[1].spec.x.value.min() == pytest.approx(x.max(), abs=0.1)


@pytest.mark.parametrize("z_em, region", [(Z_EM, "unknown"), (0.0, "lya_forest")])
def test_extract_preset_bad_input_returns_0(z_em, region):
    assert _session(_spec(z_em=z_em)).edit.extract_preset(region) == 0


def test_trim_common():
    # CIV 1548 at z = 2.97 (614.7 nm) in this session, Ly-alpha at the same z
    # (482.6 nm) in the other: the common coverage is +-200 km/s
    ly = _spec(482.0, 483.0, 0.005)
    session = _session(gui=_gui(ly=SimpleNamespace(spec=ly)))
    out = session.edit.trim_common("ly", "2.97", "CIV_1548", "Ly_a",
                                   window_kms="200.0")
    lam = 154.8204 * 3.97
    x = out.spec.x.value
    assert x.min() == pytest.approx(lam * (1 - 200 / C_KMS), abs=0.006)
    assert x.max() == pytest.approx(lam * (1 + 200 / C_KMS), abs=0.006)


def test_trim_common_mismatch_returns_0():
    # Two transitions for one session found
    session = _session(gui=_gui(ly=SimpleNamespace(spec=_spec())))
    assert session.edit.trim_common("ly", "2.97", "CIV_1548", "Ly_a, Ly_b") == 0


# --- Edit > Import systems ------------------------------------------------

@pytest.mark.parametrize("append, expected", [(True, [2.95, 2.97]),
                                              (False, [2.97])])
def test_import_systems(append, expected):
    # Of the source systems, only z = 2.97 falls in this spectrum (612-620 nm)
    source = SimpleNamespace(systs=_systs(2.97, 3.5))
    session = _session(gui=_gui(source=source), systs=_systs(2.95))
    out = session.edit.import_systems("source", append=append)
    assert sorted(c.z for c in out.systs.components) == expected


def test_import_systems_nothing_to_import():
    session = _session(gui=_gui(empty=SimpleNamespace(systs=_systs()),
                                far=SimpleNamespace(systs=_systs(3.5))))
    for name in ("missing", "empty", "far"):
        assert session.edit.import_systems(name) is session


# --- Edit > Import systems from ASCII (#480) ------------------------------

def test_import_systems_ascii_round_trip(tmp_path):
    # The _systems.csv written by Export ASCII can be read back
    source = _session(systs=_systs(2.95, 2.97))
    assert source.file.export_ascii("systems", str(tmp_path), "src")
    session = _session(systs=_systs(2.96))
    out = session.edit.import_systems_ascii(str(tmp_path / "src_systems.csv"))
    assert sorted(c.z for c in out.systs.components) == [2.95, 2.96, 2.97]
    replaced = session.edit.import_systems_ascii(str(tmp_path / "src_systems.csv"),
                                                 append="False")
    assert sorted(c.z for c in replaced.systs.components) == [2.95, 2.97]


def test_import_systems_ascii_minimal_columns(tmp_path):
    # Whitespace-separated, only series and z; out-of-range rows skipped
    path = tmp_path / "list.txt"
    path.write_text("series z b\nCIV 2.97 15.0\nCIV 3.5 10.0\nSiIV 4.0 8\n")
    out = _session().edit.import_systems_ascii(str(path))
    comps = out.systs.components
    assert [(c.series, c.z, c.b, c.logN) for c in comps] == [("CIV", 2.97, 15.0, 13.0)]


def test_import_systems_ascii_bad_input(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("name,redshift\nCIV,2.97\n")
    session = _session()
    assert session.edit.import_systems_ascii(str(path)) == 0
    assert session.edit.import_systems_ascii(str(tmp_path / "missing.csv")) == 0


# --- Edit > Stitch, Co-add ------------------------------------------------

def test_stitch():
    blue, red = _spec(612.0, 616.0), _spec(616.0, 620.0, seed=1)
    session = _session(blue, gui=_gui(red=SimpleNamespace(spec=red)))
    for others in (["red"], [_session(red)]):
        out = session.edit.stitch(others)
        np.testing.assert_allclose(out.spec.x.value,
                                   np.concatenate([blue.x.value, red.x.value]))
    # Unknown names are skipped; nothing to stitch returns the same session
    assert session.edit.stitch(["missing"]) is session
    assert session.edit.stitch([]) is session


def test_stitch_overlap_asks_for_scatter_plot():
    # Issue #442: overlapping regions give interleaved pixels, drawn as a
    # zigzag by the default step plot; stitch asks the GUI for a scatter plot
    blue, red = _spec(612.0, 616.5), _spec(616.0, 620.0, seed=1)
    session = _session(blue, gui=_gui(red=SimpleNamespace(spec=red)))
    assert session.edit.stitch(["red"]).spec.meta.get("plot_style") == "scatter"
    # Adjacent regions: no overlap, default style
    adjacent = _session(_spec(612.0, 616.0), gui=_gui(red=SimpleNamespace(
        spec=_spec(616.0, 620.0, seed=1))))
    assert "plot_style" not in adjacent.edit.stitch(["red"]).spec.meta
    # Rebinning the stitched spectrum goes back to the default style
    rebinned = session.edit.stitch(["red"]).flux.rebin(dx="5.0")
    assert "plot_style" not in rebinned.spec.meta


def test_coadd(snapshot):
    a, b = _spec(seed=1), _spec(seed=2, scale=1.1)
    gui = _gui(a=SimpleNamespace(spec=a), b=SimpleNamespace(spec=b))
    out = _session(a, gui=gui).edit.coadd("a, b", dx="0.01", xunit="nm")
    x = out.spec.x.value
    np.testing.assert_allclose(np.diff(x), 0.01, rtol=1e-6)
    snapshot.check("recipe_coadd_nm", {"x": out.spec.x, "y": out.spec.y,
                                       "dy": out.spec.dy}, rtol=1e-8)


@pytest.mark.parametrize("names, kwargs", [("", {}), ("missing", {}),
                                           ("a", dict(dx="abc"))])
def test_coadd_bad_input_returns_0(names, kwargs):
    session = _session(gui=_gui(a=SimpleNamespace(spec=_spec())))
    assert session.edit.coadd(names, **kwargs) == 0


# --- Flux recipes ---------------------------------------------------------

def test_rebin(snapshot):
    out = _session().flux.rebin(dx="5.0", xunit="km/s")
    dv = np.diff(np.log(out.spec.x.value)) * C_KMS
    np.testing.assert_allclose(dv, 5.0, rtol=1e-3)
    snapshot.check("recipe_rebin_5kms", {"x": out.spec.x, "y": out.spec.y,
                                         "dy": out.spec.dy}, rtol=1e-8)


@pytest.mark.parametrize("params", [{"dx": "abc"}, {"xunit": "not_a_unit"},
                                    {"kappa": "x"}, {"xstart": "?"}])
def test_rebin_bad_input_returns_0(params):
    # Regression test for F23: the parsing error was logged but the recipe
    # went on and raised UnboundLocalError instead of returning 0
    assert _session().flux.rebin(**params) == 0


def test_rebin_does_not_offer_norm():
    # F37 (#560): 'norm' (V1: return the normalized spectrum) was offered by
    # the dialog but never applied; it is no longer offered
    from astrocook.recipes.flux import FLUX_RECIPES_SCHEMAS
    params = [p["name"] for p in FLUX_RECIPES_SCHEMAS["rebin"]["params"]]
    assert "norm" not in params


def test_rebin_norm_is_ignored_with_a_warning(caplog):
    session = _session()
    plain = session.flux.rebin(dx="5.0")
    with caplog.at_level("WARNING"):
        normed = session.flux.rebin(dx="5.0", norm="True")
    assert "'norm' is ignored" in caplog.text
    np.testing.assert_array_equal(normed.spec.y.value, plain.spec.y.value)
    np.testing.assert_array_equal(normed.spec.cont.value, plain.spec.cont.value)
    # The alternative named in the warning
    out = plain.edit.apply_expression("y_norm", "y / cont")
    np.testing.assert_allclose(out.spec.get_column("y_norm").value,
                               plain.spec.y.value / plain.spec.cont.value)


def test_smooth_bad_input_returns_0():
    assert _session().flux.smooth("abc") == 0


@pytest.mark.parametrize("window, used", [("20", 21), ("1", 3)])
def test_running_std_window_is_made_odd(window, used):
    session = _session()
    out = session.flux.calculate_running_std(window_pix=window)
    ref = session.spec.calculate_running_std("y", "running_std", used)
    np.testing.assert_array_equal(out.spec.get_column("running_std").value,
                                  ref.get_column("running_std").value)
    assert session.flux.calculate_running_std(window_pix="abc") == 0


def test_equalize():
    ref, faint = _spec(seed=1), _spec(seed=2, scale=0.5)
    out = _session(faint, gui=_gui(ref=SimpleNamespace(spec=ref))).flux.equalize("ref")
    ratio = np.median(out.spec.y.value / ref.y.value)
    assert ratio == pytest.approx(1.0, abs=0.01)


@pytest.mark.parametrize("name, order", [("None", "0"), ("missing", "0"),
                                         ("ref", "abc")])
def test_equalize_bad_input_returns_0(name, order):
    session = _session(gui=_gui(ref=SimpleNamespace(spec=_spec())))
    assert session.flux.equalize(name, order) == 0


def test_resample():
    target = _spec(613.0, 619.0, 0.007, seed=3)
    out = _session(gui=_gui(target=SimpleNamespace(spec=target))).flux.resample("target")
    np.testing.assert_allclose(out.spec.x.value, target.x.value)


@pytest.mark.parametrize("name", ["None", "", "missing"])
def test_resample_bad_input_returns_0(name):
    assert _session(gui=_gui()).flux.resample(name) == 0
