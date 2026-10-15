"""Timing benchmark of the Voigt fitting, to compare before and after a change.

Not a test (not collected by pytest, not run in CI): run it by hand on the
same machine before and after a change that touches the model, the
convolution or the fitter, and compare the numbers. Three cases, on mock
spectra (``astrocook.mocks``) plus, optionally, the J0103-1305 data:

- ``single``: fit of one Ly-alpha component (R = 45000), repeated;
- ``group``: fit of a blended group of forest components;
- ``chain_mock``: the full chain of ``tests/full_chain.py`` on the mock
  (continuum, doublets, Ly-alpha forest), about 40 s;
- ``chain_j0103`` (``--j0103``): the same chain on J0103-1305, 500-700 nm,
  about 10 minutes.

Usage::

    python benchmarks/bench_fitting.py [--j0103] [--repeat N] [--json out.json]

The times depend on the machine: only compare runs made on the same one.
"""
import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))


def _systs(components):
    from astrocook.core.structures import ComponentDataV2, SystemListDataV2
    from astrocook.core.system_list import SystemListV2
    comps = [ComponentDataV2(id=i + 1, z=c[1], dz=None, logN=c[2], dlogN=None,
                             b=c[3], db=None, series=c[0])
             for i, c in enumerate(components)]
    return SystemListV2(SystemListDataV2(components=comps))


def bench_single(repeat):
    from astrocook.fitting.voigt_fitter import VoigtFitterV2
    from astrocook.mocks import MockComponent, MockTruth, generate
    spec = generate(MockTruth(
        components=[MockComponent("Ly_a", 4.0, 13.8, 25.0)],
        lam_min_nm=605.0, lam_max_nm=611.0, dv_kms=2.5, resol=45000.0,
        snr=50.0, z_em=4.1, seed=1))
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        VoigtFitterV2(spec, _systs([("Ly_a", 4.0, 13.5, 20.0)])).fit(
            max_nfev=2000, z_window_kms=20.0)
        times.append(time.perf_counter() - t0)
    return float(np.median(times))


def bench_group(repeat):
    from astrocook.fitting.voigt_fitter import VoigtFitterV2
    from astrocook.mocks import MockComponent, MockTruth, generate
    z0 = 3.0
    dz = 40.0 / 299792.458 * (1 + z0)
    truth = [MockComponent("Ly_a", z0 + k * dz, 13.0 + 0.1 * (k % 3), 20.0 + 3 * k)
             for k in range(8)]
    spec = generate(MockTruth(components=truth, lam_min_nm=484.0,
                              lam_max_nm=490.0, dv_kms=2.5, resol=45000.0,
                              snr=30.0, z_em=3.2, seed=2))
    guess = [("Ly_a", c.z + 2e-5, c.logN - 0.2, c.b * 0.8) for c in truth]
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        VoigtFitterV2(spec, _systs(guess)).fit(max_nfev=2000, z_window_kms=20.0)
        times.append(time.perf_counter() - t0)
    return float(np.median(times))


def bench_chain(case):
    import full_chain as fc
    fc.seed_uuids(0)
    session = {"mock": fc.mock_session, "j0103": fc.j0103_session}[case]()
    t0 = time.perf_counter()
    fc.run_chain(session)
    return time.perf_counter() - t0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--j0103", action="store_true",
                        help="also run the J0103-1305 chain (~10 min)")
    parser.add_argument("--repeat", type=int, default=10,
                        help="repetitions of the single and group fits")
    parser.add_argument("--json", help="write the results to this file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.WARNING)

    results = {
        "single_s": bench_single(args.repeat),
        "group_s": bench_group(args.repeat),
        "chain_mock_s": bench_chain("mock"),
    }
    if args.j0103:
        results["chain_j0103_s"] = bench_chain("j0103")
    for key, value in results.items():
        print(f"{key:15s} {value:10.3f}")
    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
