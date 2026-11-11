"""The review of the parameters of the recipes of File and Flux (phase
6-quater, PR 7, agreed with Guido, 11 November 2026).

Characterization first: the behaviour that the review keeps, and (before
the change) the behaviour that it changes, each test saying which. The
numbers of a valid call do not change: the cases changed are inputs that
gave an empty or hidden result without a word.
"""
import logging
from types import SimpleNamespace

import astropy.units as au
import numpy as np
import pytest

from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, SpectrumData
from astrocook.mcp_server.schemas import parameters_schema
from astrocook.methods_text import _recipe_sentence
from astrocook.recipes.errors import InvalidParameter, RecipeFailed
from astrocook.recipes.registry import get_recipe
from astrocook.session.session import Session

FLUX = au.erg / au.cm**2 / au.s / au.nm


def _spec(scale=1.0, x0=380.0, x1=460.0):
    rng = np.random.default_rng(0)
    x = np.arange(x0, x1, 0.02)
    cont = scale * (x / 420.0) ** -1.5
    y = cont * (1 + 0.02 * rng.standard_normal(len(x)))
    return Spectrum(SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x - 0.01, au.nm),
        xmax=DataColumn(x + 0.01, au.nm), y=DataColumn(y, FLUX),
        dy=DataColumn(0.02 * cont, FLUX),
        aux_cols={"cont": DataColumn(cont, FLUX)}, z_em=2.5))


@pytest.fixture
def session():
    other = Session(name="other", gui=None, spec=_spec(2.0, 390.0, 450.0))
    main = Session(name="main", gui=None, spec=_spec())
    main._gui = SimpleNamespace(session_histories=[
        SimpleNamespace(display_name="other", current_state=other)])
    return main


# --- Kept ------------------------------------------------------------------

@pytest.mark.parametrize("value", ["None", ""])
def test_equalize_without_a_session_fails_as_before(session, value):
    with pytest.raises(RecipeFailed, match="No reference session is selected"):
        session.flux.equalize(reference_session=value)


@pytest.mark.parametrize("value", ["None", ""])
def test_resample_without_a_session_fails_as_before(session, value):
    with pytest.raises(RecipeFailed, match="No target session is selected"):
        session.flux.resample(target_session=value)


def test_equalize_and_resample_without_arguments_fail_as_before(session):
    with pytest.raises(RecipeFailed, match="No reference session"):
        session.flux.equalize()
    with pytest.raises(RecipeFailed, match="No target session"):
        session.flux.resample()


def test_equalize_order_minus_one_is_a_spline(session):
    spline = session.flux.equalize(reference_session="other", order=-1)
    const = session.flux.equalize(reference_session="other", order=0)
    assert not np.allclose(spline.spec.y.value, const.spec.y.value)


@pytest.mark.parametrize("value", ["False", False])
def test_rebin_norm_false_is_silent(session, caplog, value):
    with caplog.at_level(logging.WARNING):
        out = session.flux.rebin(dx=20.0, norm=value)
    assert "norm" not in caplog.text
    ref = session.flux.rebin(dx=20.0)
    assert np.array_equal(out.spec.y.value, ref.spec.y.value, equal_nan=True)


def test_rebin_norm_true_warns_and_changes_nothing(session, caplog):
    with caplog.at_level(logging.WARNING):
        out = session.flux.rebin(dx=20.0, norm="True")
    assert "The parameter 'norm' is ignored" in caplog.text
    ref = session.flux.rebin(dx=20.0)
    assert np.array_equal(out.spec.y.value, ref.spec.y.value, equal_nan=True)


def test_rebin_with_a_valid_range(session):
    out = session.flux.rebin(xstart=400.0, xend=420.0, dx=0.1, xunit="nm")
    assert 190 <= len(out.spec.x) <= 210


@pytest.mark.parametrize("magnitudes", ["SDSS_r=17.0", " SDSS_g = 17.5 , SDSS_r=17.0 "])
def test_calibrate_with_valid_magnitudes(session, magnitudes):
    out = session.flux.calibrate_from_magnitudes(magnitudes=magnitudes)
    assert not np.allclose(out.spec.y.value, session.spec.y.value)


def test_calibrate_outside_the_spectrum_fails(session):
    with pytest.raises(RecipeFailed, match="does not overlap"):
        session.flux.calibrate_from_magnitudes(magnitudes="SDSS_z=17.0")


def test_running_std_window_rule(session):
    # An even window is increased by 1, one below 3 is 3
    assert np.array_equal(
        np.asarray(session.flux.calculate_running_std(window_pix=20).spec.t["running_std"]),
        np.asarray(session.flux.calculate_running_std(window_pix=21).spec.t["running_std"]))
    assert np.array_equal(
        np.asarray(session.flux.calculate_running_std(window_pix=1).spec.t["running_std"]),
        np.asarray(session.flux.calculate_running_std(window_pix=3).spec.t["running_std"]))


def test_export_ascii_by_position_and_keywords(session, tmp_path):
    assert session.file.export_ascii("F", str(tmp_path), "a")
    assert (tmp_path / "a_spec.csv").exists()


# --- Changed by the review -----------------------------------------------------
# Each test says what the recipe did before (characterized before the change)

def test_equalize_order_below_minus_one_is_refused(session):
    # Before: a spline, as -1
    with pytest.raises(InvalidParameter, match="order: -2 is below -1"):
        session.flux.equalize(reference_session="other", order=-2)


def test_rebin_with_start_not_below_end_is_refused(session):
    # Before: a spectrum without pixels
    for xstart, xend in [(420.0, 400.0), (410.0, 410.0)]:
        with pytest.raises(InvalidParameter, match="xstart: .* is not below xend"):
            session.flux.rebin(xstart=xstart, xend=xend, dx=20.0)


@pytest.mark.parametrize("name", ["y", "x", "xmin", "xmax", "dy", "", "  "])
def test_running_std_into_a_core_or_empty_name_is_refused(session, name):
    # Before: the result went to a hidden column, the core one unchanged
    with pytest.raises(InvalidParameter, match="output_col"):
        session.flux.calculate_running_std(output_col=name)


def test_running_std_can_replace_an_auxiliary_column(session):
    out = session.flux.calculate_running_std(output_col="cont")
    assert not np.allclose(np.asarray(out.spec.t["cont"]),
                           np.asarray(session.spec.t["cont"]))


@pytest.mark.parametrize("magnitudes, match", [
    ("Foo=17", "the filter 'Foo' is not known; the filters are Johnson_U"),
    ("SDSS_r=17,", "is not a list of 'filter=magnitude' pairs"),
    ("", "is not a list"), ("SDSS_r", "is not a list"), ("SDSS_r=abc", "is not a list"),
    ("SDSS_r=17=18", "is not a list")])
def test_calibrate_with_wrong_magnitudes_is_refused(session, magnitudes, match):
    # Before: RecipeFailed, from the core
    with pytest.raises(InvalidParameter, match=match):
        session.flux.calibrate_from_magnitudes(magnitudes=magnitudes)


@pytest.mark.parametrize("value", ["true", " TRUE ", "1", True, 1])
def test_rebin_norm_true_in_any_spelling_warns(session, caplog, value):
    # Before: only 'True' (and True) warned
    with caplog.at_level(logging.WARNING):
        session.flux.rebin(dx=20.0, norm=value)
    assert "The parameter 'norm' is ignored" in caplog.text


def test_export_ascii_items_not_found_are_named(session, tmp_path, caplog):
    # Before: skipped without a word
    with caplog.at_level(logging.WARNING):
        assert session.file.export_ascii("F, model, sytems", str(tmp_path), "a")
    assert "Not exported, not found in the session: model, sytems." in caplog.text
    with caplog.at_level(logging.WARNING):
        caplog.clear()
        assert session.file.export_ascii("spec, λ, F, dF", str(tmp_path), "b")
    assert "Not exported" not in caplog.text


def test_export_ascii_folder_and_its_old_name(session, tmp_path):
    # 'path' until this PR: an alias of the registry
    assert session.file.export_ascii(targets="F", folder=str(tmp_path / ""), prefix="new")
    assert session.file.export_ascii(targets="F", path=str(tmp_path), prefix="old")
    assert ((tmp_path / "new_spec.csv").read_text()
            == (tmp_path / "old_spec.csv").read_text())
    with pytest.raises(InvalidParameter, match="path is the old name of folder"):
        session.file.export_ascii(targets="F", folder=str(tmp_path), path=str(tmp_path),
                                  prefix="both")


def test_session_parameters_are_required_in_the_schema():
    # Before: 'None' as a default, not required
    for category, name, param in [("flux", "equalize", "reference_session"),
                                  ("flux", "resample", "target_session")]:
        schema = parameters_schema(get_recipe(category, name))
        assert "default" not in schema["properties"][param]
        assert param in schema["required"]
        assert schema["properties"][param]["x-astrocook-kind"] == "session"


@pytest.mark.parametrize("passed, used", [(20, 21), (21, 21), (1, 3)])
def test_running_std_methods_text_names_the_window_used(passed, used):
    # Before: the window passed (20 for a window of 21)
    text, _ = _recipe_sentence(get_recipe("flux", "calculate_running_std"),
                               {"input_col": "y", "output_col": "s", "window_pix": passed},
                               menus=False)
    assert f"{used} pixels" in text
    assert "Other parameters" not in text


def test_running_std_methods_text_of_a_log_without_the_window():
    text, _ = _recipe_sentence(get_recipe("flux", "calculate_running_std"),
                               {"input_col": "y", "output_col": "s"}, menus=False)
    assert "running window of" in text
