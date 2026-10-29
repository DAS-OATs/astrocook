"""Tests for the file recipes (ASCII export)."""
import numpy as np
import astropy.units as au
from astropy.io import ascii

from astrocook.core.session import SessionV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, SpectrumDataV2


def _make_session(n=20):
    x = np.linspace(400.0, 410.0, n)
    y = np.linspace(1.0, 2.0, n)
    data = SpectrumDataV2(
        x=DataColumnV2(x, au.nm),
        xmin=DataColumnV2(x - 0.25, au.nm),
        xmax=DataColumnV2(x + 0.25, au.nm),
        y=DataColumnV2(y, au.dimensionless_unscaled),
        dy=DataColumnV2(0.1 * y, au.dimensionless_unscaled),
    )
    return SessionV2(name="export_test", gui=None, spec=SpectrumV2(data))


def test_export_ascii_writes_spectrum_columns(tmp_path):
    # Regression test for #476 and #481 ("name 're' is not defined").
    session = _make_session()

    saved = session.file.export_ascii(targets="F, dF", path=str(tmp_path),
                                      prefix="qso")

    assert saved is True
    table = ascii.read(tmp_path / "qso_spec.csv", format="csv")
    assert table.colnames == ["x", "y", "dy"]
    np.testing.assert_allclose(table["x"], session.spec.x.value)
    np.testing.assert_allclose(table["y"], session.spec.y.value)
    np.testing.assert_allclose(table["dy"], session.spec.dy.value)


def test_export_ascii_maps_ui_names_on_whole_words_only(tmp_path):
    # 'λmin' must map to 'xmin', not to 'x' + 'min'.
    session = _make_session()

    saved = session.file.export_ascii(targets="λmin, λmax", path=str(tmp_path),
                                      prefix="qso")

    assert saved is True
    table = ascii.read(tmp_path / "qso_spec.csv", format="csv")
    assert table.colnames == ["x", "xmin", "xmax"]


def test_export_ascii_outputs(tmp_path):
    # The files written, as the GUI asks before exporting (F82)
    session = _make_session()
    for targets, outputs in (("F, dF", ["spec"]), ("λmin", ["spec"]), ("λ", []),
                             ("spec", ["spec"]), ("cont", []), ("systems", ["systems"]),
                             ("F, systems", ["spec", "systems"]), ("", [])):
        assert session.file.export_ascii_outputs(targets) == outputs, targets
        for f in tmp_path.iterdir():
            f.unlink()
        session.file.export_ascii(targets=targets, path=str(tmp_path), prefix="qso")
        written = sorted(f.name for f in tmp_path.iterdir())
        assert written == [f"qso_{o}.csv" for o in outputs], targets


def test_export_ascii_whole_spectrum(tmp_path):
    # 'spec': all the columns of the spectrum (as documented; it wrote
    # nothing before F82 was fixed)
    session = _make_session()
    assert session.file.export_ascii(targets="spec", path=str(tmp_path), prefix="qso")
    table = ascii.read(tmp_path / "qso_spec.csv", format="csv")
    assert table.colnames == ["x", "xmin", "xmax", "y", "dy"]
