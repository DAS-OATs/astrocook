"""Numerical snapshots of the line identification on a real spectrum.

Input: J0103-1305 window with the continuum and absorbed-region mask saved
in the session (frozen inputs, see ``conftest.py``). Parameters are the
defaults of the ``identify_lines`` recipe.
"""
import json

import numpy as np
import pytest

# Defaults of recipes/absorbers.py (identify_lines)
MULTIPLETS = ["CIV", "SiIV", "MgII", "Ly_ab"]


@pytest.fixture(scope="module")
def identified(j0103_window_with_continuum):
    return j0103_window_with_continuum.identify_lines(MULTIPLETS)


def _region_maps(meta):
    """Flatten the per-region JSON maps into arrays sorted by region id."""
    series = json.loads(meta["series_map_json"])
    ids = sorted(series, key=int)
    out = {"region": np.array([int(i) for i in ids]),
           "series": np.array([series[i] for i in ids])}
    for key in ("z", "logN", "b"):
        values = json.loads(meta[f"{key}_map_json"])
        assert sorted(values, key=int) == ids
        out[key] = np.array([values[i] for i in ids], dtype=float)
    return out


def _identifications(meta):
    """Flatten region_identifications into (region, line, score) arrays."""
    rows = sorted(
        (int(rid), line, score)
        for rid, entries in json.loads(meta["region_identifications"]).items()
        for line, score in entries)
    return {"region": np.array([r[0] for r in rows]),
            "line": np.array([r[1] for r in rows]),
            "score": np.array([r[2] for r in rows], dtype=float)}


def test_identify_lines_regions(snapshot, identified):
    abs_ids = identified.get_column("abs_ids").value
    snapshot.check("identification_j0103_abs_ids",
                   {"abs_ids": abs_ids.astype(np.int64)})


def test_identify_lines_identifications(snapshot, identified):
    snapshot.check("identification_j0103_identifications",
                   _identifications(identified.meta), rtol=1e-8)


def test_identify_lines_region_maps(snapshot, identified):
    snapshot.check("identification_j0103_region_maps",
                   _region_maps(identified.meta), rtol=1e-8)
