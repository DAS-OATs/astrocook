"""Build the telluric template of the mocks (phase 7, PR 7a.1).

Reads ``sky_telluric.fits``, the ESO SkyCalc sky model of the pre-release
``data-telluric-2026`` of the repository (also in the git history, commit
9740b65c, ``astrocook/data/sky_telluric.fits``), and writes
``astrocook/mocks/data/telluric_paranal.npz``: the optical depth
tau = -ln(trans_ma) of the molecular absorption, between 550 and 1000 nm
on the native grid of the model (0.005 nm, vacuum wavelengths), in
float32, with tau below ``TAU_FLOOR`` set to zero (a transmission changed
by less than 1e-5, for the compression). Usage::

    python benchmarks/make_telluric_template.py path/to/sky_telluric.fits

The parameters of the SkyCalc run are in the release notes of
``data-telluric-2026`` and in the header of the file: SkyCalc 2.0.9,
Paranal, airmass 1.1, PWV 2.5 mm, no line spread function.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "astrocook" / "mocks" / "data" / "telluric_paranal.npz"
LAM_RANGE_NM = (550.0, 1000.0)
#: Optical depths below this are set to zero
TAU_FLOOR = 1e-5
SOURCE_SHA256 = "aa700893aac1f6c44541a91a24928b09292581b200fc13a234ac07ee9d4d0a6a"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source", help="sky_telluric.fits")
    args = parser.parse_args()
    sha = hashlib.sha256(Path(args.source).read_bytes()).hexdigest()
    if sha != SOURCE_SHA256:
        raise SystemExit(f"Unexpected SHA-256 of {args.source}: {sha}")
    with fits.open(args.source) as hdul:
        data = hdul[1].data
        lam = np.asarray(data["lam"], dtype=float)
        trans = np.asarray(data["trans_ma"], dtype=float)
    sel = (lam >= LAM_RANGE_NM[0] - 1e-9) & (lam <= LAM_RANGE_NM[1] + 1e-9)
    lam, trans = lam[sel], trans[sel]
    step = float(np.round(np.median(np.diff(lam)), 6))
    if not np.allclose(np.diff(lam), step, rtol=0, atol=1e-6):
        raise SystemExit("The wavelength grid of the source is not uniform")
    tau = -np.log(trans)
    tau[tau < TAU_FLOOR] = 0.0
    meta = dict(
        source="ESO SkyCalc sky model, column trans_ma (molecular absorption)",
        release="data-telluric-2026 (DAS-OATs/astrocook), sky_telluric.fits",
        source_sha256=SOURCE_SHA256, skycalc_version="2.0.9",
        references="Noll et al. 2012, A&A 543, A92; Jones et al. 2013, A&A 560, A91",
        observatory="Paranal", airmass=1.1, pwv_mm=2.5, lsf="none",
        wavelengths="vacuum, nm", tau_floor=TAU_FLOOR)
    OUT.parent.mkdir(exist_ok=True)
    np.savez_compressed(OUT, lam_start_nm=np.float64(lam[0]),
                        lam_step_nm=np.float64(step),
                        tau=tau.astype(np.float32), meta=json.dumps(meta))
    print(f"{OUT}: {len(tau)} values, {OUT.stat().st_size / 1024:.0f} kB")


if __name__ == "__main__":
    main()
