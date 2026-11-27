"""The statistics of the group fits (phase 7, PR 7d.9): the record of each
fit (``FitRecord``), held by the components it fitted, made by the fitter,
named after the recipe, kept by the archive, cleared by an edit.

On the synthetic spectrum of ``test_recipes_absorbers.py`` (a blended CIV
doublet and a Ly-alpha line, R = 45000).
"""
import dataclasses
import json
import tarfile

import numpy as np
import pytest

from astrocook.core.structures import FitRecord
from astrocook.fitting.fit_record import (
    current_recipe, information_criteria, records_of, running_recipe)
from astrocook.fitting.voigt_fitter import VoigtFitter
from astrocook.session.session import Session

from test_recipes_absorbers import CIV_Z, _session, _systs, data  # noqa: F401 (fixture)


@pytest.fixture(scope="module")
def pair(data):  # noqa: F811
    """Two blended CIV components, fitted together by the recipe."""
    systs = _systs(("CIV", CIV_Z, 13.5, 10.0), ("CIV", CIV_Z + 1e-4, 13.0, 8.0))
    session = _session(data, systs)
    uuid = session.systs.components[0].uuid
    return session.absorbers.fit_component(uuid)


def test_a_fit_leaves_one_record_shared_by_its_group(pair):
    a, b = pair.systs.components
    assert a.fit is not None and a.fit is b.fit
    rec = a.fit
    assert set(rec.members) == {a.uuid, b.uuid}
    assert rec.n_free == 6  # z, logN, b of each; btur frozen
    assert rec.n_pix > rec.n_free
    assert rec.n_constraining <= rec.n_pix
    # The reduced chi2 is the one of the components (the one of the errors)
    assert rec.chi2_reduced == a.chi2
    aic, _, bic = information_criteria(rec.chi2, rec.n_pix, rec.n_free)
    assert rec.aic == pytest.approx(aic) and rec.bic == pytest.approx(bic)
    assert rec.aic == pytest.approx(rec.chi2 + 12)
    assert rec.converged is True and 0 < rec.nfev <= rec.max_nfev == 2000
    assert rec.recipe == "absorbers.fit_component"
    assert rec.time.endswith("+00:00") and not rec.final
    # The bounds of the free parameters: z in the window of the fit (20 km/s)
    lo, hi = rec.bounds_of(a.uuid, "z")
    assert lo < a.z < hi and (hi - lo) / (1 + a.z) * 299792.458 == pytest.approx(40, rel=0.01)
    assert rec.bounds_of(a.uuid, "logN") is not None
    assert rec.bounds_of(a.uuid, "btur") is None  # frozen
    assert {u for u, *_ in rec.bounds} == {a.uuid, b.uuid}


def test_the_chi2_is_that_of_the_pixels_of_the_fit(data):  # noqa: F811
    systs = _systs(("CIV", CIV_Z, 13.5, 10.0))
    fitter = VoigtFitter(data, systs)
    out, model, res = fitter.fit()
    rec = out.components[0].fit
    assert rec.chi2 == pytest.approx(float(np.sum(res.fun ** 2)))
    assert rec.n_pix == len(res.fun)
    chi2, n = fitter.chi2_on_fit_pixels(model)
    assert (rec.chi2, rec.n_pix) == (pytest.approx(chi2, rel=1e-6), n)
    # Outside a recipe there is no recipe to name
    assert rec.recipe is None


def test_a_fit_that_stops_at_max_nfev_did_not_converge(data):  # noqa: F811
    systs = _systs(("CIV", CIV_Z, 13.5, 10.0))
    out, _, res = VoigtFitter(data, systs).fit(max_nfev=2, model=False)
    rec = out.components[0].fit
    assert not res.success
    assert rec.converged is False and rec.nfev == res.nfev and rec.max_nfev == 2


def test_a_spectator_keeps_its_record(data):  # noqa: F811
    session = _session(data, _systs(("CIV", CIV_Z, 13.5, 10.0), ("Ly_a", 4.059, 13.0, 20.0)))
    civ, lya = session.systs.components
    session = session.absorbers.fit_component(civ.uuid, group_depth="0")
    first = {c.uuid: c.fit for c in session.systs.components}
    assert first[lya.uuid] is None and first[civ.uuid] is not None
    session = session.absorbers.fit_component(lya.uuid, group_depth="0")
    after = {c.uuid: c.fit for c in session.systs.components}
    assert after[civ.uuid] is first[civ.uuid]
    assert after[lya.uuid].members == (lya.uuid,)


def test_the_uncertainties_at_the_final_values_keep_the_convergence(pair, data):  # noqa: F811
    systs = pair.systs
    rec = systs.components[0].fit
    with systs.fitting_context([c.uuid for c in systs.components], group_depth=0, spec=data):
        out = VoigtFitter(data, systs).evaluate_uncertainties()
    new = out.components[0].fit
    assert new is out.components[1].fit and new.id != rec.id
    assert new.final and (new.converged, new.nfev, new.max_nfev) == (
        rec.converged, rec.nfev, rec.max_nfev)
    assert new.chi2_reduced == out.components[0].chi2
    # On the pixels of the window at the final values (the fit took them at
    # the initial values)
    assert new.n_pix > new.n_free and np.isfinite(new.chi2)


def test_refit_all_leaves_the_records_of_its_final_values(data):  # noqa: F811
    session = _session(data, _systs(("CIV", CIV_Z, 13.5, 10.0), ("Ly_a", 4.059, 13.0, 20.0)))
    out = session.absorbers.refit_all()
    recs = [c.fit for c in out.systs.components]
    assert all(r is not None and r.final and r.recipe == "absorbers.refit_all"
               and r.converged is not None for r in recs)
    assert recs[0] is not recs[1]  # two groups


def test_an_edit_clears_the_record(pair):
    uuid = pair.systs.components[0].uuid
    edited = pair.absorbers.update_component(uuid, logN="13.4")
    comp = next(c for c in edited.systs.components if c.uuid == uuid)
    assert comp.fit is None and comp.flags == ()
    other = next(c for c in edited.systs.components if c.uuid != uuid)
    assert other.fit is not None


def test_running_recipe_names_the_outermost():
    assert current_recipe() is None
    with running_recipe("absorbers.doublets_auto"):
        with running_recipe("absorbers.refit_all"):
            assert current_recipe() == "absorbers.doublets_auto"
    assert current_recipe() is None


def test_the_record_is_not_compared_with_the_component(pair):
    a = pair.systs.components[0]
    assert dataclasses.replace(a, fit=None) == a


def test_record_to_and_from_dict():
    rec = FitRecord(chi2=10.0, n_pix=20, n_free=3, aic=16.0, bic=None, chi2_reduced=0.6,
                    n_constraining=18, members=("a", "b"), converged=True, nfev=7,
                    max_nfev=100, message="ok", recipe="absorbers.fit_component",
                    time="2026-11-27T10:00:00+00:00", bounds=(("a", "z", 2.9, 3.1),))
    d = json.loads(json.dumps(rec.to_dict()))
    assert FitRecord.from_dict(d) == rec
    # A key of a later version is ignored
    assert FitRecord.from_dict({**d, "later": 1}) == rec


def test_records_of(pair):
    assert records_of([]) == {}
    rec = pair.systs.components[0].fit
    assert records_of(pair.systs.components) == {rec.id: rec}


# --- The archive -----------------------------------------------------------

def _round_trip(session, tmp_path):
    path = tmp_path / "fits.acs2"
    session.save(str(path))
    return Session.open_new(str(path), "fits", None, "auto"), path


def test_the_records_are_saved_once_and_shared_again(pair, tmp_path):
    reloaded, path = _round_trip(pair, tmp_path)
    a, b = reloaded.systs.components
    assert a.fit is not None and a.fit is b.fit
    assert a.fit == pair.systs.components[0].fit
    with tarfile.open(path) as tar:
        names = tar.getnames()
        assert "fits_fit_records.json" in names
        member = tar.extractfile("fits_fit_records.json")
        stored = json.load(member)
    assert stored["schema"] == "astrocook.fit_records"
    assert len(stored["records"]) == 1
    assert set(stored["components"]) == {a.uuid, b.uuid}


def test_an_archive_without_records_opens_without(data, tmp_path):  # noqa: F811
    session = _session(data, _systs(("CIV", CIV_Z, 13.5, 10.0)))
    reloaded, path = _round_trip(session, tmp_path)
    assert reloaded.systs.components[0].fit is None
    with tarfile.open(path) as tar:
        assert not any(n.endswith("_fit_records.json") for n in tar.getnames())


def test_a_later_schema_is_not_read(pair, tmp_path):
    from astrocook.io import archive_data
    path = tmp_path / "x_fit_records.json"
    value = archive_data._get_fit_records(pair.systs)
    archive_data._write_fit_records(value, str(path))
    d = json.loads(path.read_text())
    d["schema_version"] = 99
    path.write_text(json.dumps(d))
    with pytest.raises(ValueError, match="version 99"):
        archive_data._read_fit_records(str(path))
