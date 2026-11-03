"""Errors of core, fitting and io reported once, with their context (phase 6).

Phase 6, PR 4a: the clauses that logged an error and re-raised it (the
loaders, the expressions, the power-law and the calibration, the save of
an archive) no longer log: the error is logged once, by the caller that
catches it (a recipe, ``load_session_from_file``), with what the user
needs to know (the file and its format, the expression). The clauses that
caught every exception catch only the expected ones.
"""
import logging

import pytest

from astrocook.core.session import SessionV2, load_session_from_file
from astrocook.io.loaders import detect_file_format

from test_recipes_edit_flux import _spec


def _errors(caplog):
    return [r for r in caplog.records if r.levelno >= logging.ERROR]


def test_unreadable_file_is_reported_once_with_its_format(tmp_path, caplog):
    bad = tmp_path / "bad.fits"
    bad.write_bytes(b"not a FITS file")
    with caplog.at_level(logging.ERROR):
        assert load_session_from_file(str(bad), "bad", None, "astrocook_fits") == 0
    errors = _errors(caplog)
    assert len(errors) == 1
    assert errors[0].getMessage().startswith("Could not read bad.fits as astrocook_fits:")


def test_open_new_returns_0_for_an_unreadable_file(tmp_path, caplog):
    bad = tmp_path / "bad.fits"
    bad.write_bytes(b"not a FITS file")
    with caplog.at_level(logging.ERROR):
        assert SessionV2.open_new(str(bad), "bad", None, "astrocook_fits") == 0
    assert len(_errors(caplog)) == 1


def test_corrupted_v1_archive(tmp_path, caplog):
    bad = tmp_path / "bad.acs"
    bad.write_bytes(b"not a tar.gz archive")
    with caplog.at_level(logging.ERROR):
        assert load_session_from_file(str(bad), "bad", None, "auto") == 0
    assert "Could not unpack" in caplog.text


def test_short_csv_is_read_as_a_generic_spectrum(tmp_path):
    # F97 (characterization): the detection reads the first 10 lines, so a
    # CSV file with fewer lines is not recognized as CSV
    short = tmp_path / "short.csv"
    short.write_text("500.0,1.0,0.1\n501.0,1.0,0.1\n")
    assert detect_file_format(str(short)) == "generic_spectrum"
    long = tmp_path / "long.csv"
    long.write_text("".join(f"{500 + i}.0,1.0,0.1\n" for i in range(12)))
    assert detect_file_format(str(long)) == "simple_csv"


def test_failed_expression_is_logged_once_with_the_expression(caplog):
    session = SessionV2("s", None, spec=_spec())
    with caplog.at_level(logging.ERROR):
        assert session.edit.apply_expression(target_col="y", expression="nope * 2") == 0
    errors = _errors(caplog)
    assert len(errors) == 1
    assert "Could not apply the expression 'nope * 2'" in errors[0].getMessage()


def test_split_keeps_the_type_of_the_error():
    # numexpr raises KeyError for an unknown name: the type reaches the caller
    with pytest.raises(KeyError):
        _spec().split("nope > 1")


# --- PR 4b: the rows D of core, fitting and io ------------------------------------
# Each case reproduces the failure: the fallbacks went on with a result
# that looked valid (decided with Guido; ROADMAP.md, phase 6, "Rows D of
# core, fitting and io").

import dataclasses  # noqa: E402
import io  # noqa: E402
import tarfile  # noqa: E402

import astropy.units as au  # noqa: E402
import numpy as np  # noqa: E402

import astrocook.core.spectrum_operations as ops  # noqa: E402
import astrocook.core.system_list as system_list_module  # noqa: E402
from astrocook.core.structures import ParameterConstraintV2  # noqa: E402
from astrocook.core.system_list import SystemListV2  # noqa: E402

from conftest import DATA_DIR  # noqa: E402
from test_recipes_absorbers import _systs  # noqa: E402

SECONDS = np.arange(20.0) * au.s  # an x axis that cannot become km/s


def test_smooth_spectrum_does_not_return_the_unsmoothed_array():
    with pytest.raises(au.UnitConversionError):
        ops.smooth_spectrum(SECONDS, np.ones(20), 10.0)


def test_merge_regions_does_not_return_them_unmerged():
    regions = np.zeros(20, dtype=int)
    regions[3:6], regions[8:10] = 1, 2
    with pytest.raises(au.UnitConversionError):
        ops.merge_regions_by_velocity(regions, SECONDS, 0.0, 50.0)


def test_kinematic_check_does_not_return_no_candidates():
    with pytest.raises(au.UnitConversionError):
        ops._find_kinematic_doublet_candidates(np.zeros(20, dtype=int), SECONDS,
                                               0.0, "CIV", 3.0)


def _rate(spec):
    n = spec.x.size
    red, blue = np.zeros(n, dtype=bool), np.zeros(n, dtype=bool)
    red[60:70], blue[40:50] = True, True
    return ops.rate_doublet_candidate(spec, red, blue, "CIV", 2.25)


def test_a_candidate_that_cannot_be_rated_scores_0(monkeypatch, caplog):
    def bad_interp(*args, **kwargs):
        raise ValueError("too few points")
    monkeypatch.setattr(ops, "interp1d", bad_interp)
    with caplog.at_level(logging.DEBUG, logger="astrocook.core.spectrum_operations"):
        assert _rate(_spec()) == 0.0
    assert [r.levelno for r in caplog.records
            if "cannot be rated" in r.getMessage()] == [logging.DEBUG]


def test_an_unexpected_error_in_the_rating_reaches_the_caller(monkeypatch):
    def broken(*args, **kwargs):
        raise KeyError("bug")
    monkeypatch.setattr(ops, "interp1d", broken)
    with pytest.raises(KeyError):
        _rate(_spec())


def test_renormalize_model_does_not_hide_a_bug():
    spec = _spec()
    model = spec._data.aux_cols["model"]
    with pytest.raises(ValueError):
        spec._renormalize_model(np.ones(3), np.ones(4), model)


def test_no_fallback_constraint_model(monkeypatch):
    # With an empty fallback model a save lost the links and the frozen
    # parameters; with None it crashed
    def broken(systs):
        raise RuntimeError("constraint model")
    monkeypatch.setattr(system_list_module, "VoigtModelConstraintV2", broken)
    with pytest.raises(RuntimeError, match="constraint model"):
        _systs(("CIV", 2.0, 13.0, 10.0))


def test_unreadable_v1_system_list_fails_the_opening(tmp_path, caplog):
    # It opened the session without its systems, and a save in .acs2 lost
    # them; now the opening fails, as for a V2 archive (Guido's choice)
    broken = tmp_path / "v1_systems.acs"
    with tarfile.open(DATA_DIR / "v1_systems.acs", "r:gz") as src, \
            tarfile.open(broken, "w:gz") as dst:
        for member in src.getmembers():
            data = src.extractfile(member).read()
            if member.name.endswith("_systs.fits"):
                data = b"not a FITS table"
            member.size = len(data)
            dst.addfile(member, io.BytesIO(data))
    with caplog.at_level(logging.ERROR):
        assert load_session_from_file(str(broken), "v1", None, "auto") == 0
    assert "Could not open v1_systems.acs" in caplog.text
    # The intact archive still opens with its systems
    session = load_session_from_file(str(DATA_DIR / "v1_systems.acs"), "v1", None, "auto")
    assert session != 0 and len(session.systs.components) > 0


# --- links ---

def _two():
    systs = _systs(("CIV", 2.0, 13.0, 10.0), ("SiIV", 2.0001, 12.5, 8.0))
    return systs, [c.uuid for c in systs.components]


@pytest.mark.parametrize("expression, error", [
    ("p['{a}'].z +", "SyntaxError"),
    ("p['missing'].z", "KeyError"),
    ("p['{a}'].nope", "AttributeError"),
    ("foo(p['{a}'].z)", "NameError"),
])
def test_update_constraint_refuses_a_wrong_link(expression, error):
    systs, (a, b) = _two()
    with pytest.raises(ValueError, match=f"cannot be evaluated \\({error}"):
        systs.update_constraint(b, "z", is_free=False,
                                expression=expression.format(a=a), target_uuid=a)


def test_update_constraint_accepts_a_good_link():
    systs, (a, b) = _two()
    new = systs.update_constraint(b, "z", is_free=False,
                                  expression=f"p['{a}'].z + 1e-4", target_uuid=a)
    assert new.constraint_model._compiled_links


def _with_link(expression):
    """A system list with a link written directly, as read from an archive
    (``{a}`` in the expression: the UUID of the first component)."""
    systs, (a, b) = _two()
    link = ParameterConstraintV2(is_free=False, expression=expression.format(a=a),
                                 target_uuid=a)
    data = dataclasses.replace(systs._data, v2_constraints_map={b: {"z": link}})
    return SystemListV2(data), a, b


def test_a_wrong_link_from_an_archive_is_kept_frozen(caplog):
    with caplog.at_level(logging.ERROR):
        systs, a, b = _with_link("p['missing'].z")
    model = systs.constraint_model
    assert model._compiled_links == []
    assert "cannot be evaluated (KeyError" in caplog.text
    assert "kept frozen at its value" in caplog.text
    # Frozen: the value of z of b does not change in the full vector
    p_full = model.map_p_free_to_full(model._p_free_vector)
    assert p_full[model._uuid_to_idx[b]] == pytest.approx(2.0001)


def test_a_good_link_from_an_archive_is_applied():
    systs, a, b = _with_link("p['{a}'].z + 1e-3")
    model = systs.constraint_model
    p_full = model.map_p_free_to_full(model._p_free_vector)
    assert p_full[model._uuid_to_idx[b]] == pytest.approx(2.001)


def test_a_link_that_fails_during_the_fit_is_not_silent():
    systs, a, b = _with_link("p['{a}'].z")
    model = systs.constraint_model
    target_idx, _ = model._compiled_links[0]
    model._compiled_links.append((target_idx, compile("p['gone'].z", "<link>", "eval")))
    with pytest.raises(KeyError):
        model.map_p_free_to_full(model._p_free_vector)
