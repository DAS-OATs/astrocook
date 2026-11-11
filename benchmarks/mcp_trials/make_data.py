"""The input spectra of the trials of the MCP server with fresh agents.

Phase 6-ter, PR 7 (``ROADMAP.md``, section 4.3). Each trial gets one file,
written to the folder given (the agent of the trial must not read this
script: it holds the true parameters of the mocks)::

    python benchmarks/mcp_trials/make_data.py <trial> <folder>

with ``<trial>`` one of:

- ``mock``: the mock of the full-chain snapshots (``tests/full_chain.py``),
  without the true continuum and the noiseless flux (in the MCP pilot the
  agent fitted on the true continuum); z_em and R are set, as the generator
  sets them;
- ``j0103``: the J0103-1305 spectrum of the tests, raw: wavelength, flux,
  error and resolution only (no continuum, model, telluric model or
  components) and z_em not set, as an observed spectrum before the analysis;
- ``edit``: a mock for the editing trial, with a CIV system of two blended
  components (the true parameters are in :func:`edit_truth`).
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests"))

TRIALS = ("mock", "j0103", "edit")


def edit_truth():
    """CIV at z = 2.2 in two components 9.4 km/s apart, a single CIV at
    z = 2.25 for comparison, SiIV of the first system; z_em = 2.6, red side
    of the Lyman-alpha emission only."""
    from astrocook.mocks import MockComponent, MockTruth
    from astrocook.mocks.truth import QSO_EMISSION_LINES
    comps = [MockComponent("CIV", 2.20000, 13.6, 9.0),
             MockComponent("CIV", 2.20010, 13.3, 7.0),
             MockComponent("SiIV", 2.20000, 12.9, 9.0),
             MockComponent("CIV", 2.25000, 13.5, 11.0)]
    return MockTruth(components=comps, lam_min_nm=440.0, lam_max_nm=520.0,
                     dv_kms=2.5, resol=45000.0, snr=50.0, z_em=2.6, seed=11,
                     emission_lines=QSO_EMISSION_LINES, cont_slope=-1.0)


def make(trial: str, folder: Path) -> Path:
    from astrocook.core.session import Session
    from astrocook.core.spectrum import Spectrum
    from astrocook.core.structures import SpectrumData
    from astrocook.mocks import generate
    folder.mkdir(parents=True, exist_ok=True)
    if trial == "mock":
        import full_chain
        spec = generate(full_chain.mock_truth(), truth_columns=False)
        path = folder / "spectrum.acs2"
        Session("spectrum", None, spec=spec).save(str(path))
    elif trial == "edit":
        spec = generate(edit_truth(), truth_columns=False)
        path = folder / "spectrum.acs2"
        Session("spectrum", None, spec=spec).save(str(path))
    elif trial == "j0103":
        from conftest import J0103_SESSION
        full = Session.open_new(str(J0103_SESSION), "J0103-1305", None, "auto")
        d = full.spec._data
        data = SpectrumData(x=d.x, xmin=d.xmin, xmax=d.xmax, y=d.y, dy=d.dy,
                              aux_cols={"resol": d.aux_cols["resol"]})
        path = folder / "J0103-1305.acs2"
        Session("J0103-1305", None, spec=Spectrum(data)).save(str(path))
    else:
        raise SystemExit(f"unknown trial {trial!r}: one of {TRIALS}")
    return path


if __name__ == "__main__":
    print(make(sys.argv[1], Path(sys.argv[2])))
