"""The warnings of the libraries, silenced where they are raised (phase 6, PR 7).

The rule of the log messages: a cosmetic warning of a library is silenced
where it is raised (``warnings.catch_warnings``), never with a global filter;
a warning that tells of a loss becomes a WARNING of Astrocook that says what
is lost. Each test records every warning (``simplefilter("always")``).
"""
import logging
import warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table

from astrocook.session.session import Session, load_session_from_file
from astrocook.core.spectrum_operations import compute_auto_fudge
from astrocook.io.adapter import _fits_meta, _write_fits_table

from conftest import DATA_DIR
from test_recipes_absorbers import _systs
from test_recipes_edit_flux import _spec


def _recorded(fn):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fn()
    return result, [f"{w.category.__name__}: {w.message}" for w in caught]


def test_save_of_a_v1_session_has_no_warnings(tmp_path):
    # The metadata of a V1 archive keep the first lines of its header, with
    # the keywords of the FITS structure: astropy ignored them with one
    # warning each (TTYPEn, TFORMn, NAXIS...), at every save
    session = load_session_from_file(str(DATA_DIR / "v1_systems.acs"), "v1", None, "auto")
    assert any(k.startswith("TTYPE") for k in session.spec.meta)
    _, caught = _recorded(lambda: session.save(str(tmp_path / "v1.acs2")))
    assert caught == []


def test_the_header_written_is_the_same(tmp_path):
    # Leaving out the structural keywords does not change the file: astropy
    # ignored them anyway
    t = Table({"x": [1.0, 2.0], "y": [3.0, 4.0]})
    meta = {"NAXIS": 2, "TTYPE1": "wave", "TFORM1": "D", "OBJECT": "QSO", "Z_EM": 3.0}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t.meta.update(meta)
        t.write(tmp_path / "raw.fits")
    t.meta.clear()
    t.meta.update(_fits_meta(meta))
    _write_fits_table(t, str(tmp_path / "filtered.fits"))
    raw = fits.getheader(tmp_path / "raw.fits", 1)
    filtered = fits.getheader(tmp_path / "filtered.fits", 1)
    assert list(raw.items()) == list(filtered.items())


def test_long_keywords_are_written_without_warnings(tmp_path):
    # series_map_json and the like become HIERARCH cards, the standard way
    # to hold them: astropy warned at every save of an identification
    session = Session("s", None, spec=_spec(), systs=_systs(("CIV", 2.97, 13.6, 12.0)))
    session.systs._data.meta["series_map_json"] = '{"1": "CIV"}'
    _, caught = _recorded(lambda: session.save(str(tmp_path / "s.acs2")))
    assert caught == []
    reloaded = load_session_from_file(str(tmp_path / "s.acs2"), "s", None, "auto")
    assert reloaded.systs._data.meta["series_map_json"] == '{"1": "CIV"}'


def test_bayesian_samples_not_saved_is_said_by_astrocook(tmp_path, caplog):
    # F28: astropy skipped them with a warning of its own; now a WARNING of
    # Astrocook says what is lost (they are still not saved: the strict xfail
    # test_bayesian_results_are_saved)
    session = Session("s", None, spec=_spec(), systs=_systs(("CIV", 2.97, 13.6, 12.0)))
    session.systs._data.meta["bayesian_results"] = {"samples": [[1.0, 2.0]]}
    with caplog.at_level(logging.WARNING):
        _, caught = _recorded(lambda: session.save(str(tmp_path / "s.acs2")))
    assert caught == []
    assert "The samples of the Bayesian fit are not saved in the archive" in caplog.text


def test_auto_fudge_of_a_nan_model_has_no_warning():
    result, caught = _recorded(lambda: compute_auto_fudge(
        np.full(20, 1.1), np.zeros(20, dtype=bool), np.full(20, np.nan)))
    assert result == 1.0 and caught == []

