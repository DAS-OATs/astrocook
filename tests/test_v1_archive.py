"""Opening V1 ``.acs`` archives (``astrocook/io/v1_archive.py``, phase 3).

Two small archives cut from real V1 sessions provided by Guido:
``v1_systems.acs`` (a spectrum with all its systems fitted) and
``v1_metals_constraints.acs`` (metal systems with V1 constraints). The
spectrum is cut to 395-412 nm, the system tables are complete, and the
headers are anonymised (object, coordinates, dates, file names); the V1 log
is replaced with an empty one.

The snapshots were recorded with the former reader, which went through
the V1 classes of ``astrocook/legacy`` (before phase 3, 8b), and are
compared exactly: the native reader gives the same session.
"""
import io
import json
import tarfile
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from astrocook.core.session import load_session_from_file

DATA = Path(__file__).parent / "data"
ARCHIVES = ["v1_systems", "v1_metals_constraints"]
PARAMS = ["z", "dz", "logN", "dlogN", "b", "db", "btur", "resol"]


def _load(name, fmt="auto"):
    session = load_session_from_file(str(DATA / f"{name}.acs"), name, None, fmt)
    assert session, f"{name}.acs did not load"
    return session


@pytest.fixture(scope="module", params=ARCHIVES)
def archive(request):
    return request.param, _load(request.param)


def test_spectrum(snapshot, archive):
    name, session = archive
    data = session.spec._data
    values = {}
    for key in ("x", "xmin", "xmax", "y", "dy"):
        column = getattr(data, key)
        values[key] = column.values
        values[f"{key}__unit_str"] = np.array(str(column.unit))
    for key, column in data.aux_cols.items():
        values[f"aux_{key}"] = column.values
        values[f"aux_{key}__unit_str"] = np.array(str(column.unit))
    values["aux_order"] = np.array(list(data.aux_cols))
    values["meta_keys"] = np.array(list(data.meta))
    values["meta_values"] = np.array([str(v) for v in data.meta.values()])
    values["z_em"] = np.float64(data.z_em)
    snapshot.check(f"v1_archive_{name}_spectrum", values, rtol=0.0)


def test_systems(snapshot, archive):
    name, session = archive
    comps = session.systs.components
    values = {p: np.array([np.nan if getattr(c, p) is None else getattr(c, p)
                           for c in comps], dtype=float) for p in PARAMS}
    values["id"] = np.array([c.id for c in comps])
    values["series"] = np.array([c.series for c in comps])
    values["func"] = np.array([c.func for c in comps])
    values["dbtur_is_none"] = np.array([c.dbtur is None for c in comps])
    values["chi2_is_none"] = np.array([c.chi2 is None for c in comps])
    parsed = session.systs._data.parsed_constraints
    keys = sorted(parsed, key=str)
    values["constraint_keys"] = np.array([str(k) for k in keys])
    values["constraint_values"] = np.array([str(parsed[k]) for k in keys])
    snapshot.check(f"v1_archive_{name}_systems", values, rtol=0.0)


def test_v1_archive_contents():
    # What the tests cover, independently of the snapshots
    systems = _load("v1_systems")
    assert len(systems.spec.x) == 6316 and len(systems.systs.components) == 425
    ids = [c.id for c in systems.systs.components]
    assert ids == sorted(ids)  # V1 order of the components
    metals = _load("v1_metals_constraints")
    assert len(metals.systs.components) == 185
    parsed = metals.systs._data.parsed_constraints
    assert len(parsed) == 226
    assert parsed[(147, "z")].expression == "lines_voigt_14_z"


def test_format_name_does_not_matter():
    # The GUI passes 'auto', an .acs is always read as a V1 archive
    a = _load("v1_systems", "auto")
    b = _load("v1_systems", "archive")
    np.testing.assert_array_equal(a.spec.y.value, b.spec.y.value)


def _repack(name, tmp_path, edit):
    """Copy of an archive with one member changed by ``edit(name, bytes)``."""
    out = tmp_path / f"{name}.acs"
    with tarfile.open(DATA / f"{name}.acs") as src, tarfile.open(out, "w:gz") as dst:
        for member in src.getmembers():
            data = edit(member.name, src.extractfile(member).read())
            if data is None:
                continue
            member.size = len(data)
            dst.addfile(member, io.BytesIO(data))
    return out


def test_broken_spectrum_is_an_error(tmp_path):
    # Before phase 3 a spectrum the V1 reader could not read was replaced
    # with mock data (flux 1 between 100 and 500 nm), with only a warning
    def drop_xmin(name, data):
        if not name.endswith("_spec.fits"):
            return data
        with fits.open(io.BytesIO(data)) as hdul:
            table = hdul[1].data
            cols = fits.ColDefs([c for c in table.columns if c.name != "xmin"])
            hdu = fits.BinTableHDU.from_columns(cols, header=hdul[1].header)
            buf = io.BytesIO()
            fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(buf)
            return buf.getvalue()

    path = _repack("v1_systems", tmp_path, drop_xmin)
    assert load_session_from_file(str(path), "broken", None, "auto") == 0


def test_acs2_without_metadata_reads_the_data(tmp_path):
    # Before phase 3 an .acs2 archive without _meta.json went to the V1
    # reader and opened with mock data and no systems
    name = "J0103-1305_metallines_lya_new_3"
    full = load_session_from_file(str(DATA / f"{name}.acs2"), "full", None, "auto")
    out = tmp_path / "nometa.acs2"
    with tarfile.open(DATA / f"{name}.acs2") as src, tarfile.open(out, "w:gz") as dst:
        for member in src.getmembers():
            if not member.name.endswith("_meta.json"):
                dst.addfile(member, src.extractfile(member))
    session = load_session_from_file(str(out), "nometa", None, "auto")
    np.testing.assert_array_equal(session.spec.y.value, full.spec.y.value)
    assert len(session.systs.components) == len(full.systs.components)
    assert session.log_manager is None


def test_unknown_raw_format_is_an_error(tmp_path):
    # Before phase 3 a format without a V2 loader went to the V1 loaders
    path = tmp_path / "spectrum.fits"
    with tarfile.open(DATA / "v1_systems.acs") as tar:
        path.write_bytes(tar.extractfile("v1_systems_spec.fits").read())
    assert load_session_from_file(str(path), "raw", None, "eso_adp") == 0


def test_v1_log_is_not_read():
    # The V1 log (<name>.json) stays in the archive; the session has none
    with tarfile.open(DATA / "v1_systems.acs") as tar:
        assert json.loads(tar.extractfile("v1_systems.json").read()) == {"set_menu": []}
    assert _load("v1_systems").log_manager is None
