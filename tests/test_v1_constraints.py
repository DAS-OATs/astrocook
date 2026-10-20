"""V1 constraints of ``.acs`` archives in the V2 model (F43, phase 3).

V1 stores a link as an lmfit expression on the names of the V1 model,
``lines_voigt_<id>_<par>`` (e.g. ``lines_voigt_14_b*0.654``): it is
translated into the V2 syntax, ``p['<uuid>'].<par>``, with the target
component. Before, the linked parameters were frozen at their V1 values
with no target (the System Inspector could not show the links) and the
expression failed silently at each evaluation; constraints on ``logN``
were lost (the parameter name was lower-cased).
"""
import io
import tarfile
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from astrocook.core.session import load_session_from_file
from astrocook.io.v1_archive import parse_v1_fits_constraints
from astrocook.mocks import MockComponent, MockTruth, generate

DATA = Path(__file__).parent / "data"
RATIO = "0.65395596200518"  # SiIV b / CIV b in the V1 archive (thermal)


@pytest.fixture(scope="module")
def metals():
    return load_session_from_file(str(DATA / "v1_metals_constraints.acs"),
                                  "m", None, "auto")


def _by_id(session):
    return {int(c.id): c for c in session.systs.components}


def test_links_are_translated(metals):
    comps = _by_id(metals)
    constraints = metals.systs.constraint_model.v2_constraints_by_uuid
    target = comps[14].uuid
    z = constraints[comps[147].uuid]["z"]
    b = constraints[comps[147].uuid]["b"]
    assert (z.target_uuid, z.expression) == (target, f"p['{target}'].z")
    assert (b.target_uuid, b.expression) == (target, f"p['{target}'].b*{RATIO}")
    assert not z.is_free and not b.is_free
    # A frozen V1 parameter stays frozen, with no link
    frozen = constraints[comps[14].uuid]["b"]
    assert not frozen.is_free and frozen.expression is None


def test_links_follow_their_target(metals):
    # Move the free z of a CIV component: the linked SiIV z follows it, and
    # the SiIV b is the CIV b times the V1 ratio
    model = metals.systs.constraint_model
    comps = _by_id(metals)
    idx = model._uuid_to_idx
    civ, siiv = idx[comps[14].uuid], idx[comps[147].uuid]
    free = model._param_map["is_free"]
    assert free[civ] and not free[siiv] and not free[siiv + 2]
    p_free = model.p_free_vector.copy()
    p_free[np.flatnonzero(free).tolist().index(civ)] += 1e-4
    p = model.map_p_free_to_full(p_free)
    assert p[civ] == pytest.approx(comps[14].z + 1e-4, rel=1e-12)
    assert p[siiv] == p[civ]
    assert p[siiv + 2] == pytest.approx(p[civ + 2] * float(RATIO), rel=1e-12)


def test_constraints_survive_an_acs2_round_trip(metals, tmp_path):
    # Translated, they are V2 constraints: saved and restored by .acs2
    path = tmp_path / "metals.acs2"
    metals.save(str(path), initial_session=metals)
    reloaded = load_session_from_file(str(path), "r", None, "auto")
    comps = _by_id(reloaded)
    target = comps[14].uuid
    b = reloaded.systs.constraint_model.v2_constraints_by_uuid[comps[147].uuid]["b"]
    assert (b.target_uuid, b.expression) == (target, f"p['{target}'].b*{RATIO}")


def test_logN_constraints_are_kept():
    hdr = {"AC CONSTR ID 0": 2, "AC CONSTR PAR 0": "logN",
           "AC CONSTR VAL 0": "lines_voigt_1_logN-0.3"}
    assert list(parse_v1_fits_constraints(hdr)) == [(2, "logN")]


# --- Fit of a V1 session with linked components ---------------------------

Z, B_CIV, LOGN = 2.0, 12.0, (13.5, 13.0)
C_KMS = 299792.458


def _write_v1_archive(path, spectrum, rows, constraints):
    """A V1 .acs archive (spectrum and system tables) for ``spectrum``."""
    data = spectrum._data
    spec = fits.BinTableHDU.from_columns([
        fits.Column(name=n, format="D", array=a) for n, a in [
            ("x", data.x.values), ("xmin", data.xmin.values),
            ("xmax", data.xmax.values), ("y", data.y.values),
            ("dy", data.dy.values), ("cont", data.aux_cols["cont"].values)]])
    cols = {k: [r[k] for r in rows] for k in rows[0]}
    systs = fits.BinTableHDU.from_columns(
        [fits.Column(name=k, format="5A" if k == "func" else
                     "100A" if k == "series" else "K" if k == "id" else "D",
                     array=np.array(v)) for k, v in cols.items()])
    for n, (cid, par, val) in enumerate(constraints):
        systs.header[f"HIERARCH AC CONSTR ID {n}"] = cid
        systs.header[f"HIERARCH AC CONSTR PAR {n}"] = par
        systs.header[f"HIERARCH AC CONSTR VAL {n}"] = val
    with tarfile.open(path, "w:gz") as tar:
        for name, hdu in [("mock_spec.fits", spec), ("mock_systs.fits", systs)]:
            buf = io.BytesIO()
            fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(buf)
            info = tarfile.TarInfo(name)
            info.size = buf.tell()
            buf.seek(0)
            tar.addfile(info, buf)


def test_refit_follows_the_v1_links(tmp_path):
    ratio = float(RATIO)
    truth = MockTruth(
        components=[MockComponent("CIV", z=Z, logN=LOGN[0], b=B_CIV),
                    MockComponent("SiIV", z=Z, logN=LOGN[1], b=B_CIV * ratio)],
        lam_min_nm=416.0, lam_max_nm=466.5, dv_kms=2.5,
        resol=45000.0, snr=100.0, z_em=2.2, seed=3)
    spectrum = generate(truth)

    def row(cid, series, z, logN, b):
        return dict(func="voigt", series=series, z0=z, z=z, dz=np.nan,
                    logN=logN, dlogN=np.nan, b=b, db=np.nan, btur=0.0,
                    dbtur=np.nan, resol=45000.0, chi2r=np.nan, snr=np.nan,
                    id=cid)

    # Initial guesses off the truth; the SiIV z and b are linked to CIV
    dz = 8.0 / C_KMS * (1 + Z)
    rows = [row(1, "CIV", Z + dz, 13.3, 16.0),
            row(2, "SiIV", Z - dz, 12.8, 5.0)]
    links = [(2, "z", "lines_voigt_1_z"), (2, "b", f"lines_voigt_1_b*{RATIO}")]
    path = tmp_path / "mock.acs"
    _write_v1_archive(path, spectrum, rows, links)

    session = load_session_from_file(str(path), "mock", None, "auto")
    session = session.edit.set_properties(z_em="2.2", resol="45000")
    fitted = _by_id(session.absorbers.refit_all())
    civ, siiv = fitted[1], fitted[2]
    # The links hold after the fit...
    assert siiv.z == civ.z
    assert siiv.b == pytest.approx(civ.b * ratio, rel=1e-10)
    # ...and the fit recovers the truth (SiIV followed CIV from its guess)
    assert abs(civ.z - Z) / (1 + Z) * C_KMS < 1.0
    assert civ.b == pytest.approx(B_CIV, rel=0.1)
    assert siiv.logN == pytest.approx(LOGN[1], abs=0.1)


def test_system_inspector_shows_the_links(qtbot, metals):
    # Before F43 the tooltip read "Target: Unknown" with the V1 expression
    from PySide6.QtCore import Qt
    from astrocook.gui.system_inspector import COLUMNS, SystemTableModel
    comps = list(metals.systs.components)
    model = SystemTableModel()
    model.update_data(comps, metals.systs.constraint_model.v2_constraints_by_uuid)
    row = [int(c.id) for c in comps].index(147)
    tip = model.data(model.index(row, COLUMNS.index("b")), Qt.ToolTipRole)
    assert "<b>Target:</b> CIV (z=1.94004)" in tip
    assert f"CIV.b*{RATIO}" in tip
