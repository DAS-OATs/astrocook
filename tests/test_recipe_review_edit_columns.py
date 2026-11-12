"""The review of the parameters of the recipes of edit on columns and
expressions (phase 6-quater, PR 8, agreed with Guido, 12 November 2026),
with F120 (one set of names of the columns).

Characterization first: the behaviour that the review keeps, and (before
the change) the behaviour that it changes, each test saying which. The
numbers of a valid call do not change.
"""
from types import SimpleNamespace

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, SpectrumData
from astrocook.mcp_server.schemas import parameters_schema
from astrocook.recipes.errors import InvalidParameter, RecipeFailed
from astrocook.recipes.registry import get_recipe
from astrocook.session.session import Session

FLUX = au.erg / au.cm**2 / au.s / au.nm


def _spec(scale=1.0, x0=380.0, x1=460.0):
    rng = np.random.default_rng(0)
    x = np.arange(x0, x1, 0.02)
    cont = scale * (x / 420.0) ** -1.5
    y = cont * (1 + 0.02 * rng.standard_normal(len(x)))
    dy = 0.02 * cont
    dy[10] = 0.0
    return Spectrum(SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x - 0.01, au.nm),
        xmax=DataColumn(x + 0.01, au.nm), y=DataColumn(y, FLUX),
        dy=DataColumn(dy, FLUX),
        aux_cols={"cont": DataColumn(cont, FLUX)}, z_em=2.5))


@pytest.fixture
def session():
    other = Session(name="other", gui=None, spec=_spec(2.0, 390.0, 450.0))
    main = Session(name="main", gui=None, spec=_spec())
    main._gui = SimpleNamespace(session_histories=[
        SimpleNamespace(display_name="other", current_state=other)])
    return main


# --- Kept: the numbers of a valid call ----------------------------------------

@pytest.mark.parametrize("target,expression", [
    ("snr", "F / dF"), ("snr", "y / dy"), ("snr", "F / dy")])
def test_apply_expression_both_sets_of_names(session, target, expression):
    out = session.edit.apply_expression(target_col=target, expression=expression)
    s = session.spec
    with np.errstate(divide="ignore"):
        assert np.array_equal(out.spec.t["snr"].value, s.y.value / s.dy.value)


@pytest.mark.parametrize("target", ["F", "y"])
def test_apply_expression_writes_the_flux(session, target):
    out = session.edit.apply_expression(target_col=target, expression="F * 2")
    assert np.array_equal(out.spec.y.value, 2 * session.spec.y.value)
    assert out.spec.y.unit == FLUX


def test_apply_expression_with_another_session(session):
    out = session.edit.apply_expression(target_col="d", expression="F - s1.F",
                                        alias_map={"s1": "other"})
    assert np.isfinite(out.spec.t["d"].value).any()


@pytest.mark.parametrize("target,expression", [
    ("F", "dF <= 0"), ("y", "dy <= 0"), ("F", "(λ < 300) | (dF <= 0)")])
def test_mask_expression_both_sets_of_names(session, target, expression):
    out = session.edit.mask_expression(target_col=target, expression=expression)
    y = out.spec.y.value
    assert np.isnan(y[10]) and np.isfinite(np.delete(y, 10)).all()


@pytest.mark.parametrize("target", ["dF", "dy"])
def test_smooth_column_both_sets_of_names(session, target):
    out = session.edit.smooth_column(target_col=target, sigma_kms=50.0)
    ref = session.spec.smooth_column(target_col="dy", sigma_kms=50.0)
    assert np.array_equal(out.spec.dy.value, ref.dy.value)


def test_x_convert_round_trip(session):
    kms = session.edit.x_convert()
    assert kms.spec.x.unit == au.km / au.s
    back = kms.edit.x_convert(xunit="nm")
    assert np.allclose(back.spec.x.value, session.spec.x.value)
    ang = session.edit.x_convert(xunit="Angstrom")
    assert np.allclose(ang.spec.x.value, 10 * session.spec.x.value)


def test_y_convert_per_angstrom(session):
    out = session.edit.y_convert(yunit="erg/(Angstrom s cm^2)")
    assert np.allclose(out.spec.y.value, session.spec.y.value / 10)
    assert out.spec.t["cont"].unit == au.erg / au.cm**2 / au.s / au.Angstrom


def test_set_properties(session):
    out = session.edit.set_properties(name="new", z_em=3.0, resol="40000")
    assert out.name == "new" and out.spec._data.z_em == 3.0
    assert np.all(out.spec.t["resol"].value == 40000)


def test_export_ascii_both_sets_of_names(session, tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    session.file.export_ascii(targets="F, dF", folder=str(tmp_path / "a"), prefix="p")
    session.file.export_ascii(targets="y, dy", folder=str(tmp_path / "b"), prefix="p")
    a = (tmp_path / "a" / "p_spec.csv").read_text()
    assert a == (tmp_path / "b" / "p_spec.csv").read_text()
    # The file keeps the names of the data
    assert a.splitlines()[0] == "x,y,dy"


def test_delete_and_split_with_the_names_of_the_gui(session):
    assert "cont" not in session.edit.delete(targets="cont").spec.t.colnames
    part = session.edit.split(expression="(λ > 400) & (λ < 410)")
    assert part.spec.x.value.min() > 400 and part.spec.x.value.max() < 410


# --- Changed (each test says what it was) -------------------------------------

def test_apply_expression_target_and_expression_are_required():
    # Was: defaults 'F' and 'F / dF', which replaced the flux with its S/N
    schema = parameters_schema(get_recipe("edit", "apply_expression"))
    assert {"target_col", "expression"} <= set(schema["required"])
    assert "default" not in schema["properties"]["target_col"]
    assert "default" not in schema["properties"]["expression"]


@pytest.mark.parametrize("name", ["", " ", "a b", "1x", "s1.F"])
def test_apply_expression_refuses_a_name_no_expression_can_use(session, name):
    # Was: columns '' and 'a b' made, that no expression could name
    with pytest.raises(InvalidParameter, match="not a name of a column"):
        session.edit.apply_expression(target_col=name, expression="F * 2")


@pytest.mark.parametrize("target", ["λ", "λmin", "λmax", "x", "lambda"])
@pytest.mark.parametrize("recipe,kw", [
    ("apply_expression", {"expression": "λ * 0 + 1"}),
    ("mask_expression", {"expression": "λ > 400"}),
    ("smooth_column", {})])
def test_the_axis_is_not_written(session, recipe, kw, target):
    # Was: the axis replaced, masked (NaN) or smoothed (F127)
    with pytest.raises(InvalidParameter, match="is the axis"):
        getattr(session.edit, recipe)(target_col=target, **kw)


@pytest.mark.parametrize("recipe,kw", [
    ("mask_expression", {"expression": "F > 0"}), ("smooth_column", {})])
def test_a_missing_column_is_named_with_the_columns(session, recipe, kw):
    # Was: a RecipeFailed that blamed the expression (mask) or the core
    with pytest.raises(InvalidParameter, match="no column zz. The columns: λ, λmin, "
                       "λmax, F, dF, cont"):
        getattr(session.edit, recipe)(target_col="zz", **kw)


@pytest.mark.parametrize("name", ["", "  "])
def test_set_properties_refuses_an_empty_name(session, name):
    # Was: the name kept without a word
    with pytest.raises(InvalidParameter, match="cannot be empty"):
        session.edit.set_properties(name=name)


@pytest.mark.parametrize("unit", ["Hz", "Jy", "s"])
def test_x_convert_refuses_a_unit_that_is_not_an_axis(session, unit):
    # Was: a RecipeFailed of the conversion
    with pytest.raises(InvalidParameter, match="neither a wavelength"):
        session.edit.x_convert(xunit=unit)


@pytest.mark.parametrize("unit", ["Jy", "erg/(Hz s cm^2)", "km"])
def test_y_convert_refuses_a_unit_the_flux_cannot_take(session, unit):
    # Was: a RecipeFailed of the conversion
    with pytest.raises(InvalidParameter, match="cannot be converted"):
        session.edit.y_convert(yunit=unit)


def test_running_std_with_the_names_of_the_gui(session):
    # Was: 'F' (the button of the dialog) not found; 'y' the default
    by_gui = session.flux.calculate_running_std(input_col="F")
    by_data = session.flux.calculate_running_std(input_col="y")
    assert np.array_equal(by_gui.spec.t["running_std"].value,
                          by_data.spec.t["running_std"].value, equal_nan=True)
    assert get_recipe("flux", "calculate_running_std").param("input_col").default == "F"
    default = session.flux.calculate_running_std()
    assert np.array_equal(default.spec.t["running_std"].value,
                          by_data.spec.t["running_std"].value, equal_nan=True)


@pytest.mark.parametrize("name", ["F", "dF", "λ", "lambda", "y"])
def test_running_std_output_refuses_the_names_of_the_gui(session, name):
    # Was: 'F' made a column that the expressions read as y
    with pytest.raises(InvalidParameter, match="not a name for a new column"):
        session.flux.calculate_running_std(output_col=name)


def test_function_that_numexpr_does_not_have(session):
    # Was: "'VariableNode' object is not callable" (F128)
    with pytest.raises(RecipeFailed, match="a function that the expressions do not "
                       "have"):
        session.edit.apply_expression(target_col="m", expression="median(F)")


def test_columns_in_range_gives_the_axis_as_asked(session):
    # Was: always under 'x'
    assert set(session.spec.columns_in_range(["F"])["columns"]) == {"λ", "F"}
    data = session.spec.columns_in_range(["x", "y"])
    assert set(data["columns"]) == {"x", "y"} and data["axis"] == "x"
    data = session.spec.columns_in_range(["lambda", "dF"])
    assert data["axis"] == "lambda"
    assert np.array_equal(data["columns"]["lambda"], session.spec.x.value)


# --- One set of names, the old one accepted (F120) -----------------------------

def test_lambda_is_the_ascii_spelling_of_the_axis(session):
    a = session.edit.mask_expression(target_col="F", expression="lambda > 400")
    b = session.edit.mask_expression(target_col="F", expression="λ > 400")
    assert np.array_equal(a.spec.y.value, b.spec.y.value, equal_nan=True)
    part = session.edit.split(expression="(lambdamin > 400) & (lambdamax < 410)")
    assert part.spec.x.value.min() > 400


def test_names_of_the_column_module():
    from astrocook.core.column_names import (display_column, internal_column,
                                             internal_expression)
    assert internal_expression("(λ > 1) & (λmin < F) | dF_old + s1.F + lambda") == \
        "(x > 1) & (xmin < y) | dF_old + s1.y + x"
    assert [internal_column(n) for n in ("F", "dF", "y", "cont", "λmax")] == \
        ["y", "dy", "y", "cont", "xmax"]
    assert [display_column(n) for n in ("x", "dy", "cont")] == ["λ", "dF", "cont"]


def test_registry_documents_the_names_of_the_gui():
    # No default nor description of a column in the names of the data
    for name, param in [("edit.apply_expression", "target_col"),
                        ("edit.mask_expression", "target_col"),
                        ("edit.smooth_column", "target_col"),
                        ("flux.calculate_running_std", "input_col")]:
        spec = get_recipe(*name.split("."))
        p = spec.param(param)
        assert p.default not in ("x", "y", "dy")
        assert "'y'" not in p.doc and "'dy'" not in p.doc
