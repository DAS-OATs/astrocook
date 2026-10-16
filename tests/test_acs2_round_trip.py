"""Round trips of the V2 session archive (``.acs2``): save, reload, compare.

A session with everything the format is meant to keep (units, auxiliary
columns of different types, z_em, resolution, object name, components
with uncertainties, turbulent b, chi2 and resolution, a custom multiplet,
frozen and linked parameters, identification metadata), the real
J0103-1305 archives, and a second round trip (save, load, save, load) to
check that the format is stable.

What is lost or changed is recorded (findings register, ``ROADMAP.md``,
section 3.2): ``region_identifications`` reloaded as a dict instead of a JSON string
(F27), the Bayesian samples (F28, strict xfail), and the metadata key
``resol`` reloaded as ``RESOL`` (FITS header, see F8). The V2 session log
(F26) is now kept, without the entries undone at save time.
"""
import dataclasses
import json
import tarfile

import numpy as np
import astropy.units as au
import pytest

from astrocook.core.session import SessionV2
from astrocook.core.structures import DataColumnV2, HistoryLogV2

from conftest import DATA_DIR, J0103_SESSION
from test_recipes_absorbers import _systs
from test_recipes_edit_flux import N, _spec

SMALL_ARCHIVE = DATA_DIR / "J0103-1305_metallines_lya_new_3.acs2"


def _round_trip(session, tmp_path, name="session", initial=None):
    path = tmp_path / f"{name}.acs2"
    session.save(str(path), initial_session=initial)
    reloaded = SessionV2.open_new(str(path), name, None, "auto")
    assert reloaded != 0
    return reloaded, path


def _constraints(systs):
    return {uuid: {p: (c.is_free, c.expression, c.target_uuid)
                   for p, c in params.items()}
            for uuid, params in systs.constraint_model.v2_constraints_by_uuid.items()}


def _assert_same(a, b, meta=True):
    """Spectrum, components and constraints of ``b`` equal to those of ``a``."""
    da, db = a.spec._data, b.spec._data
    assert (db.z_em, db.resol) == (da.z_em, da.resol)
    for key in ("x", "xmin", "xmax", "y", "dy"):
        ca, cb = getattr(da, key), getattr(db, key)
        assert cb.unit == ca.unit, key
        np.testing.assert_array_equal(cb.values, ca.values, err_msg=key)
    assert sorted(db.aux_cols) == sorted(da.aux_cols)
    for key, ca in da.aux_cols.items():
        cb = db.aux_cols[key]
        assert cb.unit == ca.unit, key
        np.testing.assert_array_equal(cb.values, ca.values, err_msg=key)
    if meta:
        assert db.meta == da.meta
    assert ([dataclasses.astuple(c) for c in b.systs.components]
            == [dataclasses.astuple(c) for c in a.systs.components])
    assert [c.uuid for c in b.systs.components] == [c.uuid for c in a.systs.components]
    assert _constraints(b.systs) == _constraints(a.systs)


@pytest.fixture
def rich_session():
    spec = _spec(extra={
        "abs_mask": DataColumnV2(np.arange(N) % 3 == 0, au.dimensionless_unscaled),
        "abs_ids": DataColumnV2((np.arange(N) // 100).astype(float),
                                au.dimensionless_unscaled),
        "resol": DataColumnV2(np.linspace(40000.0, 50000.0, N),
                              au.dimensionless_unscaled),
    })
    systs = _systs(("CIV", 2.97, 13.6, 12.0), ("SiIV", 2.98, 13.0, 8.0),
                   ("CIV_1548,SiIV_1393", 2.99, 13.0, 9.0))
    session = SessionV2("rich", None, spec=spec, systs=systs)
    session = session.edit.set_properties(object="QSO J0000-0000", z_em="4.1",
                                          resol="45000", overwrite_resol_col=False)
    civ, siiv, _ = session.systs.components
    session = session.absorbers.update_component(civ.uuid, btur="3.0")
    session = session.absorbers.update_constraint(civ.uuid, "b", is_free=False)
    session = session.absorbers.update_constraint(
        siiv.uuid, "z", is_free=False, expression=f"p['{civ.uuid}'].z",
        target_uuid=civ.uuid)
    comps = [dataclasses.replace(c, dz=1e-5, dlogN=0.05, db=1.0, dbtur=0.5,
                                 chi2=1.1, resol=40000.0)
             for c in session.systs.components]
    session = session.with_new_system_list(type(session.systs)(
        dataclasses.replace(session.systs._data, components=comps)))
    meta = session.spec._data.meta
    meta["series_map_json"] = json.dumps({str(i): "CIV" for i in range(300)})
    meta["region_identifications"] = json.dumps({"1": [["CIV", 0.9]]})
    return session


def test_rich_session_round_trip(rich_session, tmp_path):
    reloaded, _ = _round_trip(rich_session, tmp_path)
    _assert_same(rich_session, reloaded, meta=False)
    before, after = rich_session.spec.meta, reloaded.spec.meta
    assert after["OBJECT"] == "QSO J0000-0000"
    assert after["series_map_json"] == before["series_map_json"]
    assert after["resol_fwhm"] == before["resol_fwhm"]
    # The resolution comes back from resol_fwhm; the key 'resol' comes back
    # upper case (short FITS keyword), see F8
    assert "resol" not in after and after["RESOL"] == before["resol"]
    # F27: the identifications come back as a dict with integer keys
    assert after["region_identifications"] == {1: [["CIV", 0.9]]}


def test_integer_aux_column_round_trip(rich_session, tmp_path):
    # Integer region ids, as written by identify_lines (this case was in
    # tests/test_v2_io_integrity.py, removed: findings register F33)
    ids = np.zeros(N, dtype=int)
    ids[N // 3: N // 2] = 1
    spec = rich_session.spec.update_column("abs_ids", ids)
    reloaded, _ = _round_trip(rich_session.with_new_spectrum(spec), tmp_path)
    np.testing.assert_array_equal(reloaded.spec.get_column("abs_ids").value, ids)


def test_second_round_trip_is_identical(rich_session, tmp_path):
    first, _ = _round_trip(rich_session, tmp_path, "first")
    second, _ = _round_trip(first, tmp_path, "second")
    _assert_same(first, second)


def test_archive_contents(rich_session, tmp_path):
    # The state before the changes is stored too (for reproducibility), but
    # it is not loaded back
    initial = SessionV2("rich", None, spec=_spec())
    _, path = _round_trip(rich_session, tmp_path, initial=initial)
    with tarfile.open(path) as tar:
        names = sorted(tar.getnames())
    assert names == ["session_meta.json", "session_spec.fits",
                     "session_spec_orig.fits", "session_systs.fits"]


def test_small_archive_round_trip(tmp_path):
    original = SessionV2.open_new(str(SMALL_ARCHIVE), "small", None, "auto")
    reloaded, _ = _round_trip(original, tmp_path)
    _assert_same(original, reloaded)
    assert len(reloaded.systs.components) == 764


@pytest.mark.slow
def test_large_archive_round_trip(tmp_path):
    # 15 MB archive, ~20 s
    original = SessionV2.open_new(str(J0103_SESSION), "large", None, "auto")
    reloaded, _ = _round_trip(original, tmp_path)
    _assert_same(original, reloaded)


def test_session_log_is_saved(rich_session, tmp_path):
    # Regression test for F26: the archive stored dataclasses.asdict of the
    # log, which was {}, and the loader did not read the log back.
    log = HistoryLogV2()
    log.add_entry("set_properties", {"z_em": "4.1"})
    log.add_entry("smooth", {"sigma_kms": "20"})
    rich_session.log_manager = log
    reloaded, path = _round_trip(rich_session, tmp_path)
    with tarfile.open(path) as tar:
        meta = json.load(tar.extractfile("session_meta.json"))
    assert meta["log_type"] == "v2"
    assert "smooth" in json.dumps(meta["log_history_json"])
    restored = reloaded.log_manager
    assert isinstance(restored, HistoryLogV2)
    assert [(e.recipe_name, e.params) for e in restored.entries] == [
        ("set_properties", {"z_em": "4.1"}), ("smooth", {"sigma_kms": "20"})]
    assert restored.current_index == 1


def test_undone_log_entries_are_not_saved(rich_session, tmp_path):
    log = HistoryLogV2()
    log.add_entry("set_properties", {"z_em": "4.1"})
    log.add_entry("smooth", {"sigma_kms": "20"})
    log.undo()
    rich_session.log_manager = log
    reloaded, _ = _round_trip(rich_session, tmp_path)
    assert [e.recipe_name for e in reloaded.log_manager.entries] == ["set_properties"]
    assert reloaded.log_manager.current_index == 0


def test_session_without_log_reloads_without_log(rich_session, tmp_path):
    reloaded, _ = _round_trip(rich_session, tmp_path)
    assert reloaded.log_manager is None


@pytest.mark.xfail(strict=True, reason=(
    "The Bayesian samples (systs meta 'bayesian_results', used by the corner "
    "plot) are not written to the archive: the FITS header of the system "
    "table cannot hold them and they are dropped without a warning. Findings "
    "register F28 (ROADMAP.md, section 3.2)."))
def test_bayesian_results_are_saved(rich_session, tmp_path):
    uuid = rich_session.systs.components[0].uuid
    samples = np.zeros((10, 3))
    meta = {"bayesian_results": {uuid: {"samples": samples, "labels": ["z", "logN", "b"],
                                        "engine": "mcmc"}}}
    systs = type(rich_session.systs)(
        dataclasses.replace(rich_session.systs._data, meta=meta))
    reloaded, _ = _round_trip(rich_session.with_new_system_list(systs), tmp_path)
    assert uuid in reloaded.systs._data.meta.get("bayesian_results", {})


@pytest.mark.parametrize("round_trip", [False, True])
def test_identification_viewer_after_round_trip(qtbot, rich_session, tmp_path,
                                                round_trip):
    # Regression test for F27: after an .acs2 round trip
    # region_identifications is a dict (the loader parses the JSON and makes
    # the keys integers), and the viewer, which called json.loads on it,
    # showed "Error parsing metadata"
    from astrocook.gui.identification_viewer_dialog import IdentificationViewerDialog
    session = _round_trip(rich_session, tmp_path)[0] if round_trip else rich_session
    dialog = IdentificationViewerDialog(session.spec, "rich")
    qtbot.addWidget(dialog)
    assert "Error" not in dialog.summary_label.text()
    assert "Found <b>" in dialog.summary_label.text()
