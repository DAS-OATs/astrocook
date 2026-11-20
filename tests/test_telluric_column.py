"""The telluric transmission of a spectrum (phase 7-ter, PR T1).

The column of the transmission is ``telluric`` (dimensionless, 1 where
there is no absorption); ``telluric_model``, the name read by the core
before 2.8.0, is renamed when a spectrum is built (from an archive, a file
or a script). The flag ``TELLURIC_CORRECTED`` of the metadata says that the
flux is already corrected: then the normalization is the continuum alone.

The characterization tests (written on master before the change) fix what
must not change: the normalization of a spectrum with ``telluric_model``
(cont × transmission), of one without a telluric column (cont), and of a
V1 archive whose ``telluric`` column is the boolean mask of V1's
``telluric_mask`` (cont: the mask is not a transmission).
"""
import astropy.units as au
import numpy as np
import pytest
from astropy.io import fits

from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, SpectrumData
from astrocook.io import loaders
from astrocook.io.v1_archive import read_v1_spectrum
from astrocook.session.session import Session

FLUX = au.erg / au.cm**2 / au.s / au.nm
N = 200


def _data(**aux):
    x = np.linspace(600.0, 610.0, N)
    cont = np.full(N, 2.0)
    cols = {"cont": DataColumn(cont, FLUX)}
    cols.update(aux)
    return SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x - 0.025, au.nm),
        xmax=DataColumn(x + 0.025, au.nm), y=DataColumn(1.5 * cont, FLUX),
        dy=DataColumn(np.full(N, 0.05), FLUX), aux_cols=cols, z_em=3.0)


TRANS = np.linspace(0.5, 1.0, N)


def _trans():
    return DataColumn(TRANS.copy(), au.dimensionless_unscaled)


# --- Characterization: what must not change -----------------------------------

def test_norm_without_a_telluric_column_is_the_continuum():
    spec = Spectrum(_data())
    np.testing.assert_array_equal(spec.norm.value, np.full(N, 2.0))
    assert spec.norm.unit == FLUX


def test_norm_with_telluric_model_is_continuum_times_transmission():
    spec = Spectrum(_data(telluric_model=_trans()))
    np.testing.assert_array_equal(spec.norm.value, 2.0 * TRANS)
    assert spec.norm.unit == FLUX


def test_norm_of_an_archive_with_telluric_model_survives_the_round_trip(tmp_path):
    session = Session("tell", None, spec=Spectrum(_data(telluric_model=_trans())))
    path = tmp_path / "tell.acs2"
    session.save(str(path))
    reloaded = Session.open_new(str(path), "tell", None, "auto")
    np.testing.assert_array_equal(reloaded.spec.norm.value, 2.0 * TRANS)


def _v1_spec_fits(path, **extra):
    x = np.linspace(600.0, 610.0, N)
    cols = [fits.Column(name=n, format='D', array=a) for n, a in (
        ("x", x), ("xmin", x - 0.025), ("xmax", x + 0.025),
        ("y", np.full(N, 3.0)), ("dy", np.full(N, 0.05)), ("cont", np.full(N, 2.0)))]
    for name, (fmt, arr) in extra.items():
        cols.append(fits.Column(name=name, format=fmt, array=arr))
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols)]).writeto(path)


V1_MASK = (np.arange(N) % 7) == 0  # V1 telluric_mask: True in the bands


def test_v1_boolean_telluric_mask_does_not_enter_the_normalization(tmp_path):
    path = tmp_path / "v1_spec.fits"
    _v1_spec_fits(path, telluric=("L", V1_MASK))
    spec = Spectrum(read_v1_spectrum(str(path)))
    np.testing.assert_array_equal(spec.norm.value, np.full(N, 2.0))


def test_smooth_smooths_the_transmission_of_telluric_model():
    spec = Spectrum(_data(telluric_model=_trans()))
    new = spec.smooth(20.0)
    smoothed = new.get_column("telluric") if new.has_aux_column("telluric") \
        else new.get_column("telluric_model")
    assert not np.array_equal(smoothed.value, TRANS)
    assert smoothed.unit == au.dimensionless_unscaled


# --- The new column and the flag ------------------------------------------------

def test_telluric_model_is_renamed_telluric_at_construction():
    spec = Spectrum(_data(telluric_model=_trans()))
    assert spec.has_aux_column("telluric") and not spec.has_aux_column("telluric_model")
    np.testing.assert_array_equal(spec.get_column("telluric").value, TRANS)
    # In place: the order of the columns is kept
    assert list(spec._data.aux_cols) == ["cont", "telluric"]


def test_norm_reads_the_telluric_column():
    spec = Spectrum(_data(telluric=_trans()))
    np.testing.assert_array_equal(spec.norm.value, 2.0 * TRANS)


def test_the_flag_leaves_the_continuum_alone():
    data = _data(telluric=_trans())
    data.meta["TELLURIC_CORRECTED"] = True
    spec = Spectrum(data)
    np.testing.assert_array_equal(spec.norm.value, np.full(N, 2.0))
    # The column stays, for the identification (design/identification.md, 7)
    assert spec.has_aux_column("telluric")


@pytest.mark.parametrize("value, corrected", [
    (True, True), (False, False), ("T", True), ("F", False), (1, True), (0, False),
    ("yes", True), ("no", False)])
def test_values_of_the_flag(value, corrected):
    from astrocook.core.telluric import is_telluric_corrected
    assert is_telluric_corrected({"TELLURIC_CORRECTED": value}) is corrected
    assert is_telluric_corrected({}) is False


def test_v1_boolean_telluric_mask_is_renamed_telluric_mask(tmp_path):
    path = tmp_path / "v1_spec.fits"
    _v1_spec_fits(path, telluric=("L", V1_MASK))
    spec = Spectrum(read_v1_spectrum(str(path)))
    assert not spec.has_aux_column("telluric")
    np.testing.assert_array_equal(spec.get_column("telluric_mask").value, V1_MASK)


def test_both_names_keep_the_one_the_core_read():
    """``telluric_model`` (read by the core before) becomes ``telluric``; a
    ``telluric`` already there (written by the ASCII telluric loader, read by
    nothing) is kept as ``telluric_unused``, not lost."""
    other = DataColumn(np.full(N, 0.7), au.dimensionless_unscaled)
    spec = Spectrum(_data(telluric=other, telluric_model=_trans()))
    np.testing.assert_array_equal(spec.get_column("telluric").value, TRANS)
    np.testing.assert_array_equal(spec.get_column("telluric_unused").value, 0.7)
    np.testing.assert_array_equal(spec.norm.value, 2.0 * TRANS)


def test_the_flag_survives_the_archive(tmp_path):
    data = _data(telluric=_trans())
    data.meta["TELLURIC_CORRECTED"] = True
    session = Session("tell", None, spec=Spectrum(data))
    path = tmp_path / "tell.acs2"
    session.save(str(path))
    reloaded = Session.open_new(str(path), "tell", None, "auto")
    np.testing.assert_array_equal(reloaded.spec.norm.value, np.full(N, 2.0))
    from astrocook.core.telluric import is_telluric_corrected
    assert is_telluric_corrected(reloaded.spec.meta)


def test_fit_continuum_renormalizes_the_model_without_the_tellurics_when_corrected():
    x = np.linspace(500.0, 510.0, 201)
    mask = (np.arange(x.size) >= 90) & (np.arange(x.size) < 110)
    trans = np.linspace(0.9, 1.0, x.size)
    data = SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x - 0.025, au.nm),
        xmax=DataColumn(x + 0.025, au.nm), y=DataColumn(np.full(x.size, 2.0), FLUX),
        dy=DataColumn(np.full(x.size, 0.05), FLUX), z_em=3.0,
        aux_cols={"abs_mask": DataColumn(mask, au.dimensionless_unscaled),
                  "cont": DataColumn(np.full(x.size, 2.0), FLUX),
                  "model": DataColumn(np.where(mask, 1.0, 2.0), FLUX),
                  "telluric": DataColumn(trans, au.dimensionless_unscaled)},
        meta={"TELLURIC_CORRECTED": True})
    spec = Spectrum(data)
    new = spec.fit_continuum(1.0, 0.0, renorm_model=True)
    cont = new.get_column("cont").value
    expected = cont * spec.get_column("model").value / 2.0
    np.testing.assert_allclose(new.get_column("model").value, expected, rtol=1e-14)


# --- The loaders (FT4) ----------------------------------------------------------

WAVE_NM = np.linspace(600.0, 610.0, 50)


def _xshooter(path, mtrans=None, flux_reduced=False):
    cols = [fits.Column(name='WAVE', format='50D', array=[WAVE_NM], unit='nm'),
            fits.Column(name='FLUX', format='50D', array=[np.ones(50)],
                        unit='erg cm**(-2) s**(-1) angstrom**(-1)'),
            fits.Column(name='ERR', format='50D', array=[np.full(50, 0.1)],
                        unit='erg cm**(-2) s**(-1) angstrom**(-1)')]
    if flux_reduced:
        cols.append(fits.Column(name='FLUX_REDUCED', format='50D', array=[np.ones(50)]))
    prim = fits.PrimaryHDU()
    prim.header['INSTRUME'] = 'XSHOOTER'
    hdus = [prim, fits.BinTableHDU.from_columns(cols, name='SPECTRUM')]
    if mtrans is not None:
        hdus.append(fits.BinTableHDU.from_columns(
            [fits.Column(name='mtrans', format='50D', array=[mtrans])],
            name='TELLURIC_DATA'))
    fits.HDUList(hdus).writeto(path)


def _load(path, fmt):
    assert loaders.detect_file_format(str(path)) == fmt
    return Spectrum(loaders.get_loader(fmt)(str(path)))


def test_xshooter_transmission_goes_to_telluric_with_the_flag(tmp_path):
    mtrans = np.linspace(0.6, 1.0, 50)
    path = tmp_path / "xsh.fits"
    _xshooter(path, mtrans=mtrans, flux_reduced=True)
    spec = _load(path, "xshooter_fits")
    np.testing.assert_array_equal(spec.get_column("telluric").value, mtrans)
    assert not spec.has_aux_column("telluric_transmission")
    from astrocook.core.telluric import is_telluric_corrected
    assert is_telluric_corrected(spec.meta)
    # The normalization is unchanged: the continuum alone
    spec = spec.update_column("cont", np.full(50, 2.0), unit=spec.y.unit)
    np.testing.assert_array_equal(spec.norm.value, np.full(50, 2.0))


def test_xshooter_without_a_telluric_model_has_no_flag(tmp_path):
    path = tmp_path / "xsh.fits"
    _xshooter(path)
    spec = _load(path, "xshooter_fits")
    assert not spec.has_aux_column("telluric")
    assert "TELLURIC_CORRECTED" not in spec.meta


def _pypeit_tellmodel(path, telluric):
    cols = [fits.Column(name='WAVE', format='50D', array=[WAVE_NM * 10]),
            fits.Column(name='TELLURIC', format='50D', array=[telluric]),
            fits.Column(name='OBJ_MODEL', format='50D', array=[np.full(50, 3.0)])]
    fits.HDUList([fits.PrimaryHDU(),
                  fits.BinTableHDU.from_columns(cols, name='MODEL')]).writeto(path)


def test_pypeit_telluric_goes_to_telluric(tmp_path):
    telluric = np.linspace(0.6, 1.0, 50)
    path = tmp_path / "tell.fits"
    _pypeit_tellmodel(path, telluric)
    spec = _load(path, "pypeit_tellmodel")
    np.testing.assert_array_equal(spec.get_column("telluric").value, telluric)
    # The flux is still the transmission (a model, not a spectrum: T2)
    np.testing.assert_array_equal(spec.y.value, telluric)
    assert "TELLURIC_CORRECTED" not in spec.meta


def test_ascii_telluric_model_column_now_enters_the_normalization(tmp_path):
    path = tmp_path / "s.dat"
    path.write_text("#wave\tflux\terror\ttrans_molecfit\n"
                    + "".join(f"{w}\t1.0\t0.1\t0.9\n" for w in WAVE_NM))
    # Read by name since T2 (FT1), the transmission still in telluric
    spec = _load(path, "ascii_table")
    spec = spec.update_column("cont", np.full(50, 2.0), unit=spec.y.unit)
    np.testing.assert_allclose(spec.norm.value, np.full(50, 1.8))


# --- set_properties --------------------------------------------------------------

def test_set_properties_sets_and_clears_the_flag():
    from astrocook.core.telluric import is_telluric_corrected
    from astrocook.recipes.errors import InvalidParameter
    session = Session("tell", None, spec=Spectrum(_data(telluric=_trans())))
    corrected = session.edit.set_properties(telluric_corrected="yes")
    assert is_telluric_corrected(corrected.spec.meta)
    np.testing.assert_array_equal(corrected.spec.norm.value, np.full(N, 2.0))
    # '_current_' keeps it
    kept = corrected.edit.set_properties(z_em="3.1")
    assert is_telluric_corrected(kept.spec.meta)
    cleared = kept.edit.set_properties(telluric_corrected="no")
    assert not is_telluric_corrected(cleared.spec.meta)
    np.testing.assert_array_equal(cleared.spec.norm.value, 2.0 * TRANS)
    # The input session is unchanged
    assert "TELLURIC_CORRECTED" not in session.spec.meta
    with pytest.raises(InvalidParameter):
        session.edit.set_properties(telluric_corrected="maybe")
