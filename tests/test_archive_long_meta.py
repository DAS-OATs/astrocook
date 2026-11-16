"""A long string in the metadata (e.g. ``region_map_json`` after the automatic
identification of doublets, issue #781) made the save of an archive fail: a
HIERARCH card cannot continue on other cards. It is now split in chunks and
joined back at load."""
import json

import pytest
from astropy.table import Table

from astrocook.io.adapter import (_fits_meta, _join_meta_chunks,
                                  _write_fits_table)
from astrocook.session.session import load_session_from_file

from test_gui_smoke import SESSION_FILE


def _long_json(n=400):
    return json.dumps({str(i): i % 7 + 0.5 for i in range(n)})


@pytest.mark.parametrize("value", [
    _long_json(), _long_json(3), "it's " * 100 + "end", json.dumps({"a'b": "c'd" * 40}),
])
def test_long_meta_string_round_trip(tmp_path, value):
    t = Table({"x": [1.0, 2.0]})
    t.meta.update(_fits_meta({"region_map_json": value, "short_json": "{}"}))
    path = str(tmp_path / "t.fits")
    _write_fits_table(t, path)
    back = Table.read(path)
    _join_meta_chunks(back.meta)
    assert back.meta["region_map_json"] == value
    assert back.meta["short_json"] == "{}"


def test_session_with_long_identification_meta_saves(tmp_path):
    session = load_session_from_file(str(SESSION_FILE), "s", None, "auto")
    value = _long_json()
    session.spec._data.meta["region_map_json"] = value
    path = tmp_path / "s.acs2"
    session.save(str(path), initial_session=session)
    loaded = load_session_from_file(str(path), "s2", None, "auto")
    assert loaded.spec.meta["region_map_json"] == value
