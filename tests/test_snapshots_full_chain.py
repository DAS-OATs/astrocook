"""Numerical snapshots of the full analysis chain.

Continuum, doublets on the red side and Lyman-alpha forest, run with the
recipe defaults as from the GUI menus (see ``full_chain.py``), on two
inputs:

- the raw J0103-1305 spectrum, 500-700 nm (Lyman-beta emission, Lyman-alpha
  forest, Lyman-alpha emission and red side), z_em = 4.072, R = 45000. It
  takes about 15 minutes, so it is marked ``slow`` and runs only with
  ``pytest --run-slow``;
- a mock spectrum (``astrocook.mocks``), 355-460 nm at z_em = 2.5, with a
  Lyman-alpha + Lyman-beta forest and CIV and SiIV systems: about 40 s.

Unlike the stage snapshots, each stage here takes the output of the
previous one, so a change anywhere in the chain shows up. The chain runs in
a child process with a fixed hash seed and seeded UUIDs: without them the
result changes from run to run (findings register F19, ``ROADMAP.md``).

These are characterization tests: they record what the chain does now,
including the questionable results listed in the register (F19-F21), not
what is correct.
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

RUNNER = Path(__file__).parent / "full_chain.py"
SNAPSHOT_DIR = Path(__file__).parent / "snapshots"


def _run_chain(case, tmp_path):
    out = tmp_path / f"{case}.npz"
    env = dict(os.environ, PYTHONHASHSEED="0")
    subprocess.run([sys.executable, str(RUNNER), case, str(out)], env=env,
                   check=True, cwd=RUNNER.parent.parent)
    with np.load(out, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def _groups(case, arrays):
    """Snapshot name -> arrays, one snapshot per stage."""
    groups = {f"full_chain_{case}_continuum":
              {k: arrays[k] for k in ("cont", "abs_mask")}}
    for stage in ("doublets", "forest"):
        groups[f"full_chain_{case}_{stage}"] = {
            k: v for k, v in arrays.items() if k.startswith(stage)}
    groups[f"full_chain_{case}_model"] = {"model_norm": arrays["model_norm"]}
    return groups


def _report_differences(case, arrays):
    """Print every difference from the stored snapshots (shown by pytest
    on failure), so that one run tells how far the whole chain moved."""
    for name, values in _groups(case, arrays).items():
        path = SNAPSHOT_DIR / f"{name}.npz"
        if not path.exists():
            continue
        with np.load(path, allow_pickle=False) as ref:
            for key in ref.files:
                expected, actual = ref[key], np.asarray(values.get(key))
                if expected.shape != actual.shape:
                    print(f"{name}[{key}]: shape {actual.shape}, expected {expected.shape}")
                elif np.issubdtype(expected.dtype, np.inexact):
                    diff = np.abs(actual.astype(float) - expected)
                    rel = diff / np.maximum(np.abs(expected), 1e-300)
                    print(f"{name}[{key}]: max abs {np.nanmax(diff, initial=0):.3g}, "
                          f"max rel {np.nanmax(rel, initial=0):.3g}, "
                          f"> 1e-6 rel: {np.sum(rel > 1e-6)}/{rel.size}")
                else:
                    print(f"{name}[{key}]: {np.sum(actual != expected)}/{expected.size} differ")


def _check(snapshot, case, arrays):
    _report_differences(case, arrays)
    groups = _groups(case, arrays)
    snapshot.check(f"full_chain_{case}_continuum", groups[f"full_chain_{case}_continuum"])
    for stage in ("doublets", "forest"):
        snapshot.check(f"full_chain_{case}_{stage}",
                       groups[f"full_chain_{case}_{stage}"], rtol=1e-6)
    # float32 keeps the file small; the tolerance is set accordingly
    snapshot.check(f"full_chain_{case}_model", groups[f"full_chain_{case}_model"],
                   rtol=1e-6, atol=1e-6)


@pytest.mark.slow
def test_full_chain_j0103(snapshot, tmp_path):
    _check(snapshot, "j0103", _run_chain("j0103", tmp_path))


def test_full_chain_mock(snapshot, tmp_path):
    _check(snapshot, "mock", _run_chain("mock", tmp_path))
