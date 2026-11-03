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
