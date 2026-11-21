"""External models in a spectrum (phase 7-ter, PR T2).

``import_model`` reads a model from a file or another session, brings it
onto the grid of the spectrum (a copy on the same grid, otherwise an
average over the pixels) and writes it in the column of its role
(``telluric``: 1 beyond the model). The reader of the tables by name
(``io/tables.py``) is shared with the loader of the ASCII spectra with
named columns (FT2), which refuses a telluric model without a flux (FT1).
``import_ascii`` refuses a table whose axis is not that of the spectrum,
and says which columns it replaced (FT3, Guido, 20 November 2026). The
role ``cont`` (T3): a continuum in the unit of the flux, NaN beyond the
model, a normalized one refused, the ``model`` following it on request.
"""
import logging

import astropy.units as au
import numpy as np
import pytest
from astropy.io import fits

from astrocook.core import model_grid
from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, SpectrumData
from astrocook.io import loaders, tables
from astrocook.recipes.errors import InvalidParameter, RecipeFailed
from astrocook.session.session import Session, load_session_from_file

C_KMS = 299792.458
N = 600


def _grid(n=N, start=500.0, dv=1.0):
    """A log-linear grid of ``dv`` km/s: centres and edges (nm)."""
    i = np.arange(n)
    x = start * np.exp(i * dv / C_KMS)
    return x, x * np.exp(-0.5 * dv / C_KMS), x * np.exp(0.5 * dv / C_KMS)


def _session(n=N, dv=1.0, **aux):
    x, lo, hi = _grid(n, dv=dv)
    cols = {"cont": DataColumn(np.full(n, 2.0), au.dimensionless_unscaled)}
    cols.update(aux)
    data = SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(lo, au.nm), xmax=DataColumn(hi, au.nm),
        y=DataColumn(np.full(n, 1.5), au.dimensionless_unscaled),
        dy=DataColumn(np.full(n, 0.05), au.dimensionless_unscaled), aux_cols=cols)
    return Session("qso", None, spec=Spectrum(data))


def _line(n):
    """A telluric line: a transmission with a Gaussian dip."""
    i = np.arange(n)
    return 1.0 - 0.6 * np.exp(-0.5 * ((i - n / 2) / 4.0) ** 2)


def _molecfit_file(path, x, lo, hi, trans, decimals=6):
    """The format of the request (V. D'Odorico): spaces, a commented
    header, the dy of the science spectrum."""
    with open(path, "w") as f:
        f.write("# x xmin xmax trans_molecfit dy\n")
        for a, b, c, t in zip(x, lo, hi, trans, strict=True):
            f.write(f"{a:.{decimals}f} {b:.{decimals}f} {c:.{decimals}f} {t:.6f} 12.3\n")
    return str(path)


# --- The grid ----------------------------------------------------------------------

def test_same_grid_is_copied():
    x, lo, hi = _grid(50)
    values = np.linspace(0.5, 1.0, 50)
    out = model_grid.model_on_grid(x + 0.05 * (hi - lo), lo, hi, values, x, lo, hi, 1.0)
    assert out.method == "copy" and out.covered.all()
    np.testing.assert_array_equal(out.values, values)


def test_average_conserves_the_integral():
    # A model on a grid three times finer: the equivalent width is the same
    x, lo, hi = _grid(200, dv=3.0)
    mx, mlo, mhi = _grid(600, dv=1.0)
    trans = _line(600)
    out = model_grid.model_on_grid(mx, mlo, mhi, trans, x, lo, hi, 1.0)
    assert out.method == "average"
    ew_model = np.sum((1 - trans) * (mhi - mlo))
    ew_spec = np.sum((1 - out.values) * (hi - lo))
    assert ew_spec == pytest.approx(ew_model, rel=1e-9)


def test_a_model_on_a_part_of_the_grid_or_a_wider_one_is_copied():
    x, lo, hi = _grid(100)
    values = np.linspace(0.5, 1.0, 100)
    # A part (rounded as in a text file), and a grid wider than the spectrum
    out = model_grid.model_on_grid(np.round(x[20:70], 6), lo[20:70], hi[20:70],
                                   values[20:70], x, lo, hi, 1.0)
    assert out.method == "copy" and out.covered.sum() == 50
    np.testing.assert_array_equal(out.values[20:70], values[20:70])
    out = model_grid.model_on_grid(x, lo, hi, values, x[10:40], lo[10:40], hi[10:40], 1.0)
    assert out.method == "copy" and out.covered.all()
    np.testing.assert_array_equal(out.values, values[10:40])


def test_pixels_beyond_the_model_get_the_fill():
    x, lo, hi = _grid(100)
    sub = slice(30, 60)
    trans = np.full(30, 0.5)
    out = model_grid.model_on_grid(x[sub], lo[sub], hi[sub], trans, x, lo, hi, 1.0)
    assert out.covered.sum() == 30
    np.testing.assert_array_equal(out.values[sub], 0.5)
    np.testing.assert_allclose(out.values[:30], 1.0, rtol=1e-9)
    np.testing.assert_allclose(out.values[60:], 1.0, rtol=1e-9)
    # NaN beyond (a continuum, T3)
    out = model_grid.model_on_grid(x[sub], lo[sub], hi[sub], trans, x, lo, hi, np.nan)
    assert np.isnan(out.values[:30]).all() and np.isnan(out.values[60:]).all()


def test_a_pixel_covered_in_part_mixes_the_model_and_the_fill():
    x, lo, hi = _grid(10)
    # The model covers half of pixel 5 and all the pixels after it
    m_lo = np.array([0.5 * (lo[5] + hi[5])])
    out = model_grid.average_over_pixels(m_lo, np.array([hi[-1]]), np.array([0.0]),
                                         lo, hi, 1.0)
    assert out.values[5] == pytest.approx(0.5, rel=1e-6)
    assert not out.covered[5] and out.covered[6:].all()
    np.testing.assert_array_equal(out.values[:5], 1.0)


def test_non_finite_values_of_the_model_are_not_used():
    x, lo, hi = _grid(20)
    values = np.full(20, 0.8)
    values[7] = np.nan
    out = model_grid.model_on_grid(x, lo, hi, values, x, lo, hi, 1.0)
    assert out.values[7] == 1.0 and not out.covered[7]


# --- The reader of the tables -------------------------------------------------------

@pytest.mark.parametrize("text", [
    "# x xmin xmax trans_molecfit dy\n500.0 499.9 500.1 0.9 3\n500.2 500.1 500.3 0.8 3\n",
    "x\txmin\txmax\ttrans_molecfit\tdy\n500.0\t499.9\t500.1\t0.9\t3\n500.2\t500.1\t500.3\t0.8\t3\n",
    "#x,xmin,xmax,trans_molecfit,dy\n500.0,499.9,500.1,0.9,3\n500.2,500.1,500.3,0.8,3\n",
    "# a comment\nx xmin xmax trans_molecfit dy\n500.0 499.9 500.1 0.9 3\n500.2 500.1 500.3 0.8 3\n",
])
def test_text_tables_by_name(tmp_path, text):
    path = tmp_path / "m.dat"
    path.write_text(text)
    table = tables.read_table(str(path))
    assert table.colnames == ["x", "xmin", "xmax", "trans_molecfit", "dy"]
    x, xmin, xmax = tables.axis_of(table)
    np.testing.assert_array_equal(x.to_value(au.nm), [500.0, 500.2])
    np.testing.assert_array_equal(xmin.to_value(au.nm), [499.9, 500.1])
    np.testing.assert_array_equal(table["trans_molecfit"], [0.9, 0.8])


def test_text_table_without_a_header(tmp_path):
    path = tmp_path / "m.txt"
    path.write_text("500.0 0.9\n500.2 0.8\n")
    table = tables.read_table(str(path))
    assert table.colnames == ["col1", "col2"]
    x, xmin, xmax = tables.axis_of(table)
    np.testing.assert_allclose(x.to_value(au.nm), [500.0, 500.2])
    np.testing.assert_allclose(xmin.to_value(au.nm), [499.9, 500.1])


def test_find_column_prefers_the_exact_name():
    from astropy.table import Table
    table = Table({"Wave": [1.0], "x": [2.0], "MTRANS": [0.5]})
    assert tables.find_column(table, tables.AXIS_NAMES) == "x"
    assert tables.find_column(table, ("mtrans",)) == "MTRANS"
    assert tables.find_column(table, ("flux",)) is None


def test_fits_of_one_row_with_the_axis_in_another_extension(tmp_path):
    # As an X-Shooter product: SPECTRUM (WAVE in nm) and TELLURIC_DATA
    # (mtrans), both tables of one row of arrays
    wave = np.linspace(600.0, 601.0, 30)
    mtrans = np.linspace(0.7, 1.0, 30)
    spec = fits.BinTableHDU.from_columns(
        [fits.Column(name="WAVE", format="30D", unit="nm", array=wave[None]),
         fits.Column(name="FLUX", format="30D", array=np.ones((1, 30)))], name="SPECTRUM")
    tell = fits.BinTableHDU.from_columns(
        [fits.Column(name="mtrans", format="30D", array=mtrans[None])], name="TELLURIC_DATA")
    path = tmp_path / "xsh.fits"
    fits.HDUList([fits.PrimaryHDU(), spec, tell]).writeto(path)
    table = tables.read_table(str(path), wanted=("mtrans",))
    np.testing.assert_array_equal(table["mtrans"], mtrans)
    x, _, _ = tables.axis_of(table)
    np.testing.assert_array_equal(x.to_value(au.nm), wave)


def test_fits_axis_in_micron(tmp_path):
    # The Molecfit outputs give lambda in micron
    lam = np.linspace(0.6, 0.601, 20)
    hdu = fits.BinTableHDU.from_columns(
        [fits.Column(name="lambda", format="D", unit="um", array=lam),
         fits.Column(name="mtrans", format="D", array=np.full(20, 0.9))])
    path = tmp_path / "tac.fits"
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)
    x, _, _ = tables.axis_of(tables.read_table(str(path)))
    np.testing.assert_allclose(x.to_value(au.nm), lam * 1000)


# --- The recipe ---------------------------------------------------------------------

def test_molecfit_file_of_the_request(tmp_path, caplog):
    # A sub-range of the grid of the spectrum, written with six decimals:
    # the model goes to its pixels, 1 beyond
    session = _session()
    x, lo, hi = _grid()
    sub = slice(100, 400)
    trans = _line(300)
    path = _molecfit_file(tmp_path / "m.dat", x[sub], lo[sub], hi[sub], trans)
    with caplog.at_level(logging.INFO):
        out = session.file.import_model(file_path=path)
    tel = out.spec.get_column("telluric")
    assert tel.unit == au.dimensionless_unscaled
    np.testing.assert_allclose(tel.value[sub], trans, atol=1e-5)
    np.testing.assert_allclose(tel.value[:100], 1.0, rtol=1e-9)
    np.testing.assert_allclose(tel.value[400:], 1.0, rtol=1e-9)
    # The normalization is cont × telluric; the error of the spectrum is
    # not the dy of the file
    np.testing.assert_allclose(out.spec.norm.value, 2.0 * tel.value)
    np.testing.assert_array_equal(out.spec.dy.value, session.spec.dy.value)
    assert "copied" in caplog.text and "covers 300 of 600 pixels; beyond, 1" in caplog.text
    # The input session is unchanged
    assert not session.spec.has_aux_column("telluric")


def test_the_same_grid_is_copied(tmp_path, caplog):
    session = _session()
    x, lo, hi = _grid()
    trans = _line(N)
    path = tmp_path / "m.dat"
    path.write_text("# x trans_molecfit\n"
                    + "".join(f"{a:.4f} {t:.6f}\n" for a, t in zip(x, trans, strict=True)))
    with caplog.at_level(logging.INFO):
        out = session.file.import_model(file_path=str(path))
    np.testing.assert_allclose(out.spec.get_column("telluric").value, trans, atol=1e-6)
    assert "copied (the model is on the grid of the spectrum)" in caplog.text


def test_a_finer_model_is_averaged(tmp_path):
    session = _session(n=200, dv=3.0)
    mx, mlo, mhi = _grid(600, dv=1.0)
    trans = _line(600)
    path = _molecfit_file(tmp_path / "m.dat", mx, mlo, mhi, trans, decimals=8)
    out = session.file.import_model(file_path=path)
    _, lo, hi = _grid(200, dv=3.0)
    tel = out.spec.get_column("telluric").value
    assert np.sum((1 - tel) * (hi - lo)) == pytest.approx(
        np.sum((1 - trans) * (mhi - mlo)), rel=1e-5)


def test_a_model_that_does_not_overlap_is_refused(tmp_path):
    session = _session()
    x, lo, hi = _grid(50, start=700.0)
    path = _molecfit_file(tmp_path / "m.dat", x, lo, hi, np.full(50, 0.9))
    with pytest.raises(RecipeFailed, match="does not overlap the spectrum"):
        session.file.import_model(file_path=path)


@pytest.mark.parametrize("params", [{}, {"file_path": "m.dat", "source_session": "other"}])
def test_one_source_and_only_one(params):
    with pytest.raises(InvalidParameter, match="file_path or source_session"):
        _session().file.import_model(**params)


def test_a_missing_file(tmp_path):
    with pytest.raises(RecipeFailed, match="does not exist"):
        _session().file.import_model(file_path=str(tmp_path / "none.dat"))


def test_a_file_without_a_column_of_the_model(tmp_path):
    path = tmp_path / "m.dat"
    path.write_text("# x xmin xmax a b\n500 499.9 500.1 1 2\n")
    with pytest.raises(RecipeFailed, match="no column of a telluric model.*a, b"):
        _session().file.import_model(file_path=str(path))
    # Named: found
    x, lo, hi = _grid()
    path.write_text("# x a\n" + "".join(f"{v:.7f} 0.5\n" for v in x))
    out = _session().file.import_model(file_path=str(path), column="a")
    np.testing.assert_allclose(out.spec.get_column("telluric").value, 0.5)
    with pytest.raises(RecipeFailed, match="no column c"):
        _session().file.import_model(file_path=str(path), column="c")


def test_a_role_that_is_not_known():
    with pytest.raises(InvalidParameter, match="role"):
        _session().file.import_model(file_path="m.dat", role="model")


def test_a_transmission_with_a_unit_is_refused(tmp_path):
    from astropy.table import Table
    x, lo, hi = _grid()
    table = Table({"lam": x * au.nm, "telluric": np.full(N, 0.5) * au.erg})
    path = tmp_path / "m.ecsv"
    table.write(path, format="ascii.ecsv")
    with pytest.raises(RecipeFailed, match="dimensionless"):
        _session().file.import_model(file_path=str(path))


def test_values_outside_0_1_are_a_warning(tmp_path, caplog):
    x, lo, hi = _grid()
    path = _molecfit_file(tmp_path / "m.dat", x, lo, hi, np.full(N, 3.0))
    with caplog.at_level(logging.WARNING):
        _session().file.import_model(file_path=path)
    assert "a transmission is between 0 and 1" in caplog.text


def test_a_flux_flagged_as_corrected_is_a_warning(tmp_path, caplog):
    x, lo, hi = _grid()
    path = _molecfit_file(tmp_path / "m.dat", x, lo, hi, np.full(N, 0.5))
    session = _session().edit.set_properties(telluric_corrected="yes")
    with caplog.at_level(logging.WARNING):
        out = session.file.import_model(file_path=path)
    assert "does not enter the normalization" in caplog.text
    np.testing.assert_array_equal(out.spec.norm.value, 2.0)


def test_from_another_session():
    trans = _line(N)
    model = _session(telluric=DataColumn(trans, au.dimensionless_unscaled))
    out = _session().file.import_model(source_session=model)
    np.testing.assert_array_equal(out.spec.get_column("telluric").value, trans)
    # A column given by name, also the flux (a model opened as a spectrum)
    out = _session().file.import_model(source_session=model, column="F")
    np.testing.assert_array_equal(out.spec.get_column("telluric").value, 1.5)
    with pytest.raises(RecipeFailed, match="no column of a telluric model"):
        _session().file.import_model(source_session=_session())


def test_the_session_by_its_name_and_the_log(tmp_path):
    # As the GUI and the MCP server: a name among the open sessions, and
    # the log records the parameters by name
    from types import SimpleNamespace
    from astrocook.recipes.params import resolve_params
    model = _session(telluric=DataColumn(_line(N), au.dimensionless_unscaled))
    model.name = "model"
    session = _session()
    session._gui = SimpleNamespace(session_histories=[
        SimpleNamespace(display_name="model", current_state=model)])
    out = session.file.import_model(source_session="model")
    np.testing.assert_array_equal(out.spec.get_column("telluric").value, _line(N))
    assert resolve_params("import_model", {"source_session": model}) == {
        "file_path": "", "source_session": "model", "role": "telluric", "column": "",
        "renorm_model": True}


def test_an_axis_in_velocity_is_refused(tmp_path):
    x, lo, hi = _grid()
    path = _molecfit_file(tmp_path / "m.dat", x, lo, hi, np.full(N, 0.5))
    vel = _session().edit.x_convert(z_ref=0.0, xunit="km/s")
    with pytest.raises(RecipeFailed, match="cannot be compared"):
        vel.file.import_model(file_path=path)


def test_the_archive_keeps_the_column(tmp_path):
    x, lo, hi = _grid()
    path = _molecfit_file(tmp_path / "m.dat", x, lo, hi, _line(N))
    out = _session().file.import_model(file_path=path)
    out.save(str(tmp_path / "s.acs2"))
    back = load_session_from_file(str(tmp_path / "s.acs2"), "s", None, "auto")
    np.testing.assert_allclose(back.spec.get_column("telluric").value, _line(N), atol=1e-6)


# --- The loaders (FT1, FT2) ---------------------------------------------------------

def test_a_telluric_model_opened_as_a_spectrum_is_refused(tmp_path, caplog):
    x, lo, hi = _grid(50)
    path = _molecfit_file(tmp_path / "m.dat", x, lo, hi, np.full(50, 0.9))
    assert loaders.detect_file_format(path) == "ascii_table"
    with pytest.raises(ValueError, match="File > Import Model"):
        loaders.get_loader("ascii_table")(path)
    # The opening fails with that message in the log
    with caplog.at_level(logging.ERROR):
        assert load_session_from_file(path, "m", None, "auto") == 0
    assert "Import Model" in caplog.text


def test_the_old_telluric_format_is_a_spectrum_with_its_transmission(tmp_path):
    path = tmp_path / "s.dat"
    path.write_text("#wave\tflux\terror\ttrans_molecfit\n"
                    + "".join(f"{w}\t1.0\t0.1\t0.9\n" for w in np.linspace(500, 501, 20)))
    assert loaders.detect_file_format(str(path)) == "ascii_table"
    data = loaders.get_loader("ascii_table")(str(path))
    np.testing.assert_array_equal(data.aux_cols["telluric"].values, 0.9)
    # FT1: the transmission is not the absorption model any more
    assert "model" not in data.aux_cols
    # The old format name still works in a script
    assert loaders.get_loader("ascii_telluric_model") is loaders.get_loader("ascii_table")


def test_an_ascii_spectrum_with_named_columns_opens(tmp_path):
    # FT2: it went to the generic loader ("No SIMPLE card found")
    path = tmp_path / "s.dat"
    path.write_text("# x xmin xmax y dy cont\n"
                    "500.0 499.9 500.1 1.0 0.1 2.0\n500.2 500.1 500.3 1.5 0.2 2.0\n")
    assert loaders.detect_file_format(str(path)) == "ascii_table"
    session = load_session_from_file(str(path), "s", None, "auto")
    np.testing.assert_array_equal(session.spec.x.value, [500.0, 500.2])
    np.testing.assert_array_equal(session.spec.xmin.value, [499.9, 500.1])
    np.testing.assert_array_equal(session.spec.y.value, [1.0, 1.5])
    np.testing.assert_array_equal(session.spec.dy.value, [0.1, 0.2])
    np.testing.assert_array_equal(session.spec.get_column("cont").value, 2.0)
    assert str(session.spec.x.unit) == "nm"


def test_an_ascii_spectrum_without_a_header_opens_by_position(tmp_path):
    path = tmp_path / "s.txt"
    path.write_text("500.0 1.0 0.1\n500.2 1.5 0.2\n500.4 2.0 0.3\n")
    assert loaders.detect_file_format(str(path)) == "ascii_table"
    data = loaders.get_loader("ascii_table")(str(path))
    np.testing.assert_array_equal(data.y.values, [1.0, 1.5, 2.0])


def test_a_csv_spectrum_is_still_read_by_position(tmp_path):
    # Unchanged (F100): a CSV file whose names are not known by name
    path = tmp_path / "s.csv"
    path.write_text("wave,flux,err\n500.0,1.0,0.1\n501.0,1.5,0.2\n")
    assert loaders.detect_file_format(str(path)) == "simple_csv"


# --- import_ascii (FT3) -------------------------------------------------------------

def test_import_ascii_refuses_another_axis(tmp_path):
    session = _session(n=20)
    x, _, _ = _grid(20)
    path = tmp_path / "c.csv"
    path.write_text("x,cont\n" + "".join(f"{v + 0.5},3.0\n" for v in x))
    with pytest.raises(RecipeFailed, match="Import Model"):
        session.file.import_ascii(file_path=str(path))


def test_import_ascii_same_axis_and_the_columns_replaced(tmp_path, caplog):
    session = _session(n=20)
    x, _, _ = _grid(20)
    path = tmp_path / "c.csv"
    path.write_text("x,cont,snr\n" + "".join(f"{v:.6f},3.0,9.0\n" for v in x))
    with caplog.at_level(logging.WARNING):
        out = session.file.import_ascii(file_path=str(path))
    np.testing.assert_array_equal(out.spec.get_column("cont").value, 3.0)
    assert "replaced by those of the file: cont." in caplog.text


def test_import_ascii_without_an_axis_is_unchanged(tmp_path, caplog):
    session = _session(n=20)
    path = tmp_path / "c.csv"
    path.write_text("snr\n" + "9.0\n" * 20)
    with caplog.at_level(logging.WARNING):
        out = session.file.import_ascii(file_path=str(path))
    np.testing.assert_array_equal(out.spec.get_column("snr").value, 9.0)
    assert "replaced" not in caplog.text


# --- The dialog ---------------------------------------------------------------------

def test_the_dialog_has_a_browse_button(qapp, monkeypatch, tmp_path):
    from astrocook.gui import recipe_dialog
    from astrocook.gui.recipe_dialog import RecipeDialog
    dialog = RecipeDialog("file", "import_model", _session())
    assert set(dialog.input_widgets) == {"file_path", "source_session", "role", "column"}
    chosen = str(tmp_path / "m.dat")
    monkeypatch.setattr(recipe_dialog.QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (chosen, "")))
    dialog.browse_button.click()
    sent = []
    dialog.recipe_requested.connect(lambda cat, rec, params, aliases: sent.append(params))
    dialog.accept()
    assert sent == [{"file_path": chosen, "source_session": "", "role": "telluric",
                     "column": ""}]


def test_a_transmission_a_little_above_1_is_not_a_warning(tmp_path, caplog):
    # The Molecfit model of the request reaches 1 plus a little (Guido, T2)
    x, lo, hi = _grid()
    trans = np.full(N, 0.5)
    trans[10] = 1.0000004
    path = _molecfit_file(tmp_path / "m.dat", x, lo, hi, trans)
    with caplog.at_level(logging.WARNING):
        _session().file.import_model(file_path=path)
    assert "between 0 and 1" not in caplog.text


def test_the_dialog_comes_back_in_front_after_browse(qapp, monkeypatch, tmp_path):
    from astrocook.gui import recipe_dialog
    from astrocook.gui.recipe_dialog import RecipeDialog
    dialog = RecipeDialog("file", "import_model", _session())
    raised = []
    monkeypatch.setattr(dialog, "raise_", lambda: raised.append("raise"))
    monkeypatch.setattr(dialog, "activateWindow", lambda: raised.append("activate"))
    monkeypatch.setattr(recipe_dialog.QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (str(tmp_path / "m.dat"), "")))
    dialog.browse_button.click()
    assert raised == ["raise", "activate"]


# --- The role cont (T3) -------------------------------------------------------------

def _physical_session(n=N, **aux):
    """A spectrum with a physical flux (~1e-16 erg/(s cm2 Angstrom)),
    without a continuum."""
    x, lo, hi = _grid(n)
    flam = au.erg / (au.s * au.cm ** 2 * au.Angstrom)
    data = SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(lo, au.nm), xmax=DataColumn(hi, au.nm),
        y=DataColumn(np.full(n, 1.5e-16), flam), dy=DataColumn(np.full(n, 5e-18), flam),
        aux_cols=dict(aux))
    return Session("qso", None, spec=Spectrum(data))


def _cont_file(path, x, cont, name="cont"):
    with open(path, "w") as f:
        f.write(f"# x {name}\n")
        for a, c in zip(x, cont, strict=True):
            f.write(f"{a:.7f} {c:.6e}\n")
    return str(path)


def test_a_continuum_beyond_the_model_is_nan(tmp_path, caplog):
    session = _physical_session()
    x, _, _ = _grid()
    sub = slice(100, 400)
    cont = np.linspace(2.0e-16, 3.0e-16, 300)
    path = _cont_file(tmp_path / "c.dat", x[sub], cont)
    with caplog.at_level(logging.INFO):
        out = session.file.import_model(file_path=path, role="cont")
    got = out.spec.get_column("cont")
    assert got.unit == session.spec.y.unit
    np.testing.assert_allclose(got.value[sub], cont, rtol=1e-6)
    assert np.all(np.isnan(got.value[:100])) and np.all(np.isnan(got.value[400:]))
    assert "The continuum of c.dat (column cont) was copied" in caplog.text
    assert "covers 300 of 600 pixels; beyond, NaN" in caplog.text
    assert "300 of 600 pixels have no continuum" in caplog.text
    # The normalization is the new continuum (NaN beyond it)
    np.testing.assert_allclose(out.spec.norm.value[sub], cont, rtol=1e-6)


def test_a_continuum_averaged_keeps_only_the_pixels_covered_in_full(tmp_path):
    session = _physical_session(n=200)
    mx, mlo, mhi = _grid(600, dv=1.0 / 3.0)
    path = tmp_path / "c.dat"
    path.write_text("# x xmin xmax continuum\n" + "".join(
        f"{a:.9f} {b:.9f} {c:.9f} 2.0e-16\n" for a, b, c in zip(mx, mlo, mhi, strict=True)))
    out = session.file.import_model(file_path=str(path), role="cont")
    got = out.spec.get_column("cont").value
    covered = np.isfinite(got)
    assert 0 < np.count_nonzero(covered) < 200
    np.testing.assert_allclose(got[covered], 2.0e-16, rtol=1e-6)


def test_a_continuum_with_a_unit_is_converted(tmp_path, caplog):
    from astropy.table import Table
    session = _physical_session()
    x, _, _ = _grid()
    flam = au.erg / (au.s * au.cm ** 2 * au.Angstrom)
    cont = np.full(N, 2.0e-16) * flam
    table = Table({"lam": x * au.nm,
                   "cont": cont.to(au.W / (au.m ** 2 * au.nm))})
    path = tmp_path / "c.ecsv"
    table.write(path, format="ascii.ecsv")
    with caplog.at_level(logging.INFO):
        out = session.file.import_model(file_path=str(path), role="cont")
    np.testing.assert_allclose(out.spec.get_column("cont").value, 2.0e-16, rtol=1e-9)
    assert "converted from" in caplog.text
    # Per unit frequency: converted at the wavelength of each pixel
    table["cont"] = cont.to(au.erg / (au.s * au.cm ** 2 * au.Hz),
                            equivalencies=au.spectral_density(x * au.nm))
    table.write(path, format="ascii.ecsv", overwrite=True)
    out = session.file.import_model(file_path=str(path), role="cont")
    np.testing.assert_allclose(out.spec.get_column("cont").value, 2.0e-16, rtol=1e-9)


def test_a_continuum_in_a_unit_that_is_not_of_the_flux_is_refused(tmp_path):
    from astropy.table import Table
    x, _, _ = _grid()
    path = tmp_path / "c.ecsv"
    Table({"lam": x * au.nm, "cont": np.full(N, 2.0) * au.km}).write(
        path, format="ascii.ecsv")
    with pytest.raises(RecipeFailed, match="a continuum is in the unit of the flux"):
        _physical_session().file.import_model(file_path=str(path), role="cont")


def test_a_normalized_continuum_is_refused(tmp_path):
    x, _, _ = _grid()
    path = _cont_file(tmp_path / "c.dat", x, np.full(N, 1.0))
    with pytest.raises(RecipeFailed, match="a normalized continuum is not one"):
        _physical_session().file.import_model(file_path=path, role="cont")
    # On the scale of the flux (a factor 10 at most): taken
    path = _cont_file(tmp_path / "c.dat", x, np.full(N, 1.2e-15))
    out = _physical_session().file.import_model(file_path=path, role="cont")
    np.testing.assert_allclose(out.spec.get_column("cont").value, 1.2e-15, rtol=1e-6)


def test_the_model_follows_the_new_continuum(tmp_path, caplog):
    # The model of a line on the old continuum (2): rescaled to the new
    # one (4) where it is finite, kept beyond it
    x, _, _ = _grid()
    profile = _line(N)
    model = DataColumn(2.0 * profile, au.dimensionless_unscaled, "the fit")
    sub = slice(0, 300)
    path = _cont_file(tmp_path / "c.dat", x[sub], np.full(300, 4.0))
    with caplog.at_level(logging.INFO):
        out = _session(model=model).file.import_model(file_path=path, role="cont")
    got = out.spec.get_column("model").value
    np.testing.assert_allclose(got[sub], 4.0 * profile[sub], rtol=1e-6)
    np.testing.assert_array_equal(got[300:], 2.0 * profile[300:])
    assert "kept as it was on 300 pixels without a continuum" in caplog.text
    # Not asked: the model as it was
    out = _session(model=model).file.import_model(file_path=path, role="cont",
                                                  renorm_model=False)
    np.testing.assert_array_equal(out.spec.get_column("model").value, 2.0 * profile)
    # A telluric model does not touch the model
    path = _molecfit_file(tmp_path / "m.dat", *_grid(), np.full(N, 0.5))
    out = _session(model=model).file.import_model(file_path=path)
    np.testing.assert_array_equal(out.spec.get_column("model").value, 2.0 * profile)


def test_a_continuum_from_another_session():
    x, _, _ = _grid()
    other = _session()            # cont = 2, in the unit of its flux
    out = _session().file.import_model(source_session=other, role="cont")
    np.testing.assert_array_equal(out.spec.get_column("cont").value, 2.0)
    # A physical continuum into an arbitrary flux: refused
    phys = _physical_session(cont=DataColumn(np.full(N, 2e-16), _physical_session().spec.y.unit))
    arb = _session()
    with pytest.raises(RecipeFailed, match="a continuum is in the unit of the flux"):
        arb.edit.set_properties(yunit="arb_flux").file.import_model(
            source_session=phys, role="cont")


@pytest.mark.parametrize("role, asked", [("cont", True), ("telluric", False)])
def test_the_gui_asks_to_renormalize_for_a_continuum(role, asked):
    # As for fit_continuum: the window asks, the answer goes in renorm_model
    from types import SimpleNamespace
    from PySide6.QtWidgets import QMessageBox
    from astrocook.gui.main_recipes import MainRecipesMixin
    session = _session(model=DataColumn(np.ones(N), au.dimensionless_unscaled))
    calls = []
    window = SimpleNamespace(
        active_history=SimpleNamespace(current_state=session),
        _show_custom_message=lambda **k: calls.append(k) or QMessageBox.StandardButton.No)
    params = {"file_path": "c.dat", "role": role}
    MainRecipesMixin._ask_to_renormalize_model(window, "import_model", params)
    assert bool(calls) == asked
    assert params.get("renorm_model") == ("False" if asked else None)
