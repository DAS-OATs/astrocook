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

Even so, the automatic pipelines amplify floating-point differences: on
another CPU, the J0103-1305 chain gives one Lyman-alpha component less and
a model that differs by up to 0.8 on single pixels, although the continuum
agrees to 1e-15. An exact comparison would also fail on any harmless
refactoring that changes the order of operations. So only the continuum is
compared exactly; the fitted stages are compared statistically
(:func:`_compare_components`, :func:`_compare_model`), with thresholds a
few times wider than the differences seen between runs with different
groupings (F19) and between machines, and much narrower than a real
regression (e.g. a wrong resolution or continuum).

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
C_KMS = 299792.458

# Thresholds of the statistical comparison. Observed so far (J0103-1305,
# runs with different groupings and on another machine): counts within 1
# component, >= 94.7% of the components matched, 90th percentile of
# |dlogN| <= 0.005 and of |db| <= 0.11 km/s, < 2% of the pixels of the
# normalised model differing by more than 0.01.
MATCH_DV_KMS = 5.0       # a component is matched within this velocity
COUNT_RTOL = 0.02        # number of components per series
MIN_MATCHED = 0.90       # fraction of reference components matched
LOGN_P50, LOGN_P90 = 0.005, 0.05   # |dlogN| of matched components (dex)
B_P50, B_P90 = 0.05, 0.5           # |db| of matched components (km/s)
ERR_RTOL_P50 = 0.05      # median relative difference of the uncertainties
MODEL_DIFF = 0.01        # a model pixel differs if |d model_norm| > this
MODEL_MAX_FRAC = 0.05    # ... on at most this fraction of the pixels
MODEL_P99 = 0.1          # 99th percentile of |d model_norm|


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


def _load(name):
    path = SNAPSHOT_DIR / f"{name}.npz"
    if not path.exists():
        pytest.fail(f"Missing snapshot {path.name}: run pytest --update-snapshots")
    with np.load(path, allow_pickle=False) as ref:
        return {key: ref[key] for key in ref.files}


def _compare_components(ref, act, stage):
    """Statistics of the components of ``act`` against ``ref``; returns the
    list of failed checks (empty if all pass)."""
    def col(d, key):
        return d[f"{stage}_{key}"]

    failures = []
    for series in sorted(set(col(ref, "series")) | set(col(act, "series"))):
        n_ref = int(np.sum(col(ref, "series") == series))
        n_act = int(np.sum(col(act, "series") == series))
        if abs(n_act - n_ref) > max(1, COUNT_RTOL * n_ref):
            failures.append(f"{series}: {n_act} components, expected {n_ref}")

    # Match each reference component with the nearest one of the same series
    pairs = []
    for i, (z, series) in enumerate(zip(col(ref, "z"), col(ref, "series"))):
        same = np.flatnonzero(col(act, "series") == series)
        if len(same):
            j = same[np.argmin(np.abs(col(act, "z")[same] - z))]
            if abs(col(act, "z")[j] - z) / (1 + z) * C_KMS < MATCH_DV_KMS:
                pairs.append((i, j))
    n_ref = len(col(ref, "z"))
    matched = len(pairs) / n_ref if n_ref else 1.0
    print(f"{stage}: {len(col(act, 'z'))} components (expected {n_ref}), "
          f"{matched:.1%} matched")
    if matched < MIN_MATCHED:
        failures.append(f"only {matched:.1%} of the components matched")
    if not pairs:
        return failures

    i, j = np.array(pairs).T
    for key, p50, p90 in (("logN", LOGN_P50, LOGN_P90), ("b", B_P50, B_P90)):
        diff = np.abs(col(act, key)[j] - col(ref, key)[i])
        q50, q90 = np.percentile(diff, [50, 90])
        print(f"  |d{key}|: median {q50:.3g}, 90th percentile {q90:.3g}")
        if q50 > p50 or q90 > p90:
            failures.append(f"|d{key}| median {q50:.3g} (max {p50}), "
                            f"90th percentile {q90:.3g} (max {p90})")
    for key in ("dz", "dlogN", "db"):
        expected, actual = col(ref, key)[i], col(act, key)[j]
        ok = np.isfinite(expected) & np.isfinite(actual) & (expected != 0)
        rel = np.median(np.abs(actual[ok] / expected[ok] - 1)) if ok.any() else 0.0
        print(f"  {key}: median relative difference {rel:.3g}")
        if rel > ERR_RTOL_P50:
            failures.append(f"{key}: median relative difference {rel:.3g} "
                            f"(max {ERR_RTOL_P50})")
    return failures


def _compare_model(ref, act):
    diff = np.abs(act["model_norm"].astype(float) - ref["model_norm"])
    frac = np.mean(diff > MODEL_DIFF)
    p99 = np.percentile(diff, 99)
    print(f"model: {frac:.2%} of the pixels differ by > {MODEL_DIFF}, "
          f"99th percentile {p99:.3g}")
    failures = []
    if frac > MODEL_MAX_FRAC or p99 > MODEL_P99:
        failures.append(f"model: {frac:.2%} of the pixels differ by > {MODEL_DIFF} "
                        f"(max {MODEL_MAX_FRAC:.0%}), 99th percentile {p99:.3g} "
                        f"(max {MODEL_P99})")
    return failures


def _check(snapshot, case, arrays):
    groups = _groups(case, arrays)
    if snapshot.update:
        for name, values in groups.items():
            snapshot.check(name, values)
        return

    # The continuum does not go through any fit: compared exactly
    snapshot.check(f"full_chain_{case}_continuum",
                   groups[f"full_chain_{case}_continuum"], rtol=1e-8)
    failures = []
    for stage in ("doublets", "forest"):
        failures += _compare_components(_load(f"full_chain_{case}_{stage}"),
                                        groups[f"full_chain_{case}_{stage}"],
                                        stage)
    failures += _compare_model(_load(f"full_chain_{case}_model"),
                               groups[f"full_chain_{case}_model"])
    assert not failures, "\n".join(failures)


@pytest.mark.slow
def test_full_chain_j0103(snapshot, tmp_path):
    _check(snapshot, "j0103", _run_chain("j0103", tmp_path))


def test_full_chain_mock(snapshot, tmp_path):
    _check(snapshot, "mock", _run_chain("mock", tmp_path))
