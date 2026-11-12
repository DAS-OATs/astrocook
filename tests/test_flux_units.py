"""The units of the flux and the CSV files (phase 6-quater, PR 8-bis, with
F129 and F130; agreed with Guido, 12 November 2026).

Characterization first (the previous commit): the unit that each loader
assigned to the flux, the conversion of the flux, the CSV files of
``export_ascii`` and ``import_ascii``; each test says what it was. No value
changes, only units, except the continuum converted with the flux (F130)
and the columns of a CSV file read by name (F129).
"""
import astropy.units as au
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, SpectrumData
from astrocook.io import loaders
from astrocook.core.units import ARBITRARY_FLUX, parse_flux_unit
from astrocook.recipes.errors import InvalidParameter, RecipeFailed
from astrocook.session.session import Session, load_session_from_file

from conftest import DATA_DIR, J0103_SESSION

FLAM = au.erg / (au.cm**2 * au.s * au.Angstrom)
FLAM_NM = au.erg / (au.cm**2 * au.s * au.nm)
WAVE_AA = np.linspace(4000.0, 4010.0, 50)


# --- Synthetic files, one per loader -----------------------------------------

def _astrocook_fits(path, y_unit=FLAM_NM):
    t = Table()
    t['x'] = WAVE_AA / 10 * au.nm
    t['y'] = np.ones(50) * y_unit
    t['dy'] = np.full(50, 0.1) * y_unit
    t['cont'] = np.ones(50) * y_unit
    t.write(path, format='fits')
    with fits.open(path, mode='update') as hdul:
        hdul[0].header['ORIGIN'] = 'Astrocook'


def _pypeit_coadd(path):
    t = Table({'wave': WAVE_AA, 'flux': np.ones(50), 'ivar': np.ones(50),
               'mask': np.ones(50, dtype=int)})
    fits.HDUList([fits.PrimaryHDU(), fits.table_to_hdu(t)]).writeto(path)
    with fits.open(path, mode='update') as hdul:
        hdul[1].name = 'SPECTRUM'


def _pypeit_spec1d(path):
    t = Table({'OPT_WAVE': WAVE_AA, 'OPT_FLAM': np.ones(50),
               'OPT_FLAM_SIG': np.full(50, 0.1)})
    hdu = fits.table_to_hdu(t)
    hdu.name = 'SPAT0100-SLIT0001-DET01-ORDER0030'
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)


def _pypeit_stack(path):
    t = Table({'wave_stack': [WAVE_AA], 'flux_stack': [np.ones(50)],
               'ivar_stack': [np.ones(50)]})
    hdu = fits.table_to_hdu(t)
    hdu.name = 'SPECTRUM'
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)


def _pypeit_tellmodel(path):
    t = Table({'WAVE': [WAVE_AA], 'TELLURIC': [np.full(50, 0.9)],
               'OBJ_MODEL': [np.ones(50)]})
    hdu = fits.table_to_hdu(t)
    hdu.name = 'MODEL'
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)


def _espresso(path, bunit=None):
    wave = np.vstack([WAVE_AA, WAVE_AA + 20])
    prim = fits.PrimaryHDU()
    sci = fits.ImageHDU(np.ones((2, 50)), name='SCIDATA')
    if bunit:
        sci.header['BUNIT'] = bunit
    fits.HDUList([prim, sci, fits.ImageHDU(np.full((2, 50), 0.1), name='ERRDATA'),
                  fits.ImageHDU(wave, name='WAVEDATA_VAC_BARY')]).writeto(path)


def _xshooter(path, tunit='erg cm**(-2) s**(-1) angstrom**(-1)', fluxcal=None):
    cols = [fits.Column(name='WAVE', format='50D', array=[WAVE_AA / 10], unit='nm'),
            fits.Column(name='FLUX', format='50D', array=[np.ones(50)], unit=tunit),
            fits.Column(name='ERR', format='50D', array=[np.full(50, 0.1)], unit=tunit)]
    prim = fits.PrimaryHDU()
    prim.header['INSTRUME'] = 'XSHOOTER'
    if fluxcal:
        prim.header['FLUXCAL'] = fluxcal
    fits.HDUList([prim, fits.BinTableHDU.from_columns(cols)]).writeto(path)


def _generic(path, tunit=None):
    cols = [fits.Column(name='LAMBDA', format='D', array=WAVE_AA),
            fits.Column(name='FLUX', format='D', array=np.ones(50), unit=tunit),
            fits.Column(name='ERR', format='D', array=np.full(50, 0.1), unit=tunit)]
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols)]).writeto(path)


def _simple_csv(path):
    path.write_text("".join(f"{w / 10},1.0,0.1\n" for w in WAVE_AA))


def _telluric_ascii(path):
    path.write_text("#wave\tflux\terror\ttrans_molecfit\n"
                    + "".join(f"{w / 10}\t1.0\t0.1\t0.9\n" for w in WAVE_AA))


def _resvel(path):
    path.write_text("RESVEL = 6.6\n" + "".join(f"{w} 1.0 0.1\n" for w in WAVE_AA))


FORMATS = {
    'astrocook_fits': ('s.fits', _astrocook_fits),
    'pypeit_coadd': ('s.fits', _pypeit_coadd),
    'pypeit_spec1d_multi': ('s.fits', _pypeit_spec1d),
    'pypeit_stack': ('s.fits', _pypeit_stack),
    'pypeit_tellmodel': ('s.fits', _pypeit_tellmodel),
    'espresso_s2d_fits': ('s.fits', _espresso),
    'xshooter_fits': ('s.fits', _xshooter),
    'generic_spectrum': ('s.fits', _generic),
    'simple_csv': ('s.csv', _simple_csv),
    'ascii_telluric_model': ('s.dat', _telluric_ascii),
    'ascii_resvel_header': ('s.dat', _resvel),
}


def _load(tmp_path, fmt, *args, **kwargs):
    name, make = FORMATS[fmt]
    path = tmp_path / name
    make(path, *args, **kwargs)
    assert loaders.detect_file_format(str(path)) == fmt
    return loaders.get_loader(fmt)(str(path))


ARB = ARBITRARY_FLUX

# The unit of the flux of each loader, and before PR 8-bis
UNITS = {
    'astrocook_fits': (FLAM_NM, FLAM),          # TUNIT read; it was not
    'pypeit_coadd': (ARB, au.dimensionless_unscaled),
    'pypeit_spec1d_multi': (ARB, au.dimensionless_unscaled),
    'pypeit_stack': (ARB, au.dimensionless_unscaled),
    # The transmission stays dimensionless (its object model is arb_flux)
    'pypeit_tellmodel': (au.dimensionless_unscaled, au.dimensionless_unscaled),
    'espresso_s2d_fits': (ARB, au.dimensionless_unscaled),
    'xshooter_fits': (FLAM, FLAM),              # TUNIT read; it was assumed
    'generic_spectrum': (ARB, au.dimensionless_unscaled),
    'simple_csv': (ARB, FLAM),                  # FLAM was assumed
    'ascii_telluric_model': (ARB, FLAM),        # FLAM was assumed
    'ascii_resvel_header': (ARB, au.dimensionless_unscaled),
}


@pytest.mark.parametrize("fmt", sorted(FORMATS))
def test_unit_of_the_flux_of_each_loader(tmp_path, fmt):
    data = _load(tmp_path, fmt)
    assert data.y.unit == UNITS[fmt][0] and data.dy.unit == UNITS[fmt][0]


def test_object_model_of_the_telluric_file(tmp_path):
    assert _load(tmp_path, 'pypeit_tellmodel').aux_cols['obj_model'].unit == ARB


def test_astrocook_fits_without_unit_and_its_continuum(tmp_path):
    data = _load(tmp_path, 'astrocook_fits', y_unit=au.dimensionless_unscaled)
    assert data.y.unit == ARB and data.aux_cols['cont'].unit == ARB


@pytest.mark.parametrize("bunit,unit", [
    ("adu", au.adu / au.Angstrom), (None, ARB), ("not a unit", ARB)])
def test_espresso_bunit(tmp_path, bunit, unit):
    # Per pixel, divided by the width of the pixel in Angstrom
    assert _load(tmp_path, 'espresso_s2d_fits', bunit=bunit).y.unit == unit


@pytest.mark.parametrize("tunit,fluxcal,unit", [
    ("erg cm**(-2) s**(-1) angstrom**(-1)", "ABSOLUTE", FLAM),
    ("erg cm**(-2) s**(-1) angstrom**(-1)", "UNCALIBRATED", ARB),
    ("adu", None, au.adu), ("", None, ARB)])
def test_xshooter_tunit_and_fluxcal(tmp_path, tunit, fluxcal, unit):
    assert _load(tmp_path, 'xshooter_fits', tunit=tunit, fluxcal=fluxcal).y.unit == unit


@pytest.mark.parametrize("tunit,unit", [
    ("1e-17 erg/(s cm2 Angstrom)", 1e-17 * FLAM), ("Jy", au.Jy), (None, ARB)])
def test_generic_tunit(tmp_path, tunit, unit):
    assert _load(tmp_path, 'generic_spectrum', tunit=tunit).y.unit == unit


def test_v1_archive_is_in_arbitrary_units():
    # Before: erg/(Angstrom s cm2), which V1 assumed; the V1 archives write no TUNIT
    session = Session.open_new(str(DATA_DIR / "v1_systems.acs"), "v1", None, "auto")
    assert session.spec.y.unit == ARB and session.spec.get_column('cont').unit == ARB


# --- The unit arb_flux --------------------------------------------------------

def test_arb_flux_is_a_unit_of_its_own():
    assert au.Unit('arb_flux') is ARB
    assert ARB != au.dimensionless_unscaled and not ARB.is_equivalent(FLAM)


@pytest.mark.parametrize("text,unit", [
    ("arb_flux", ARB), ("erg/(s cm2 Angstrom)", FLAM), ("10**(-17) erg/s/cm**2/Angstrom", 1e-17 * FLAM),
    ("", None), (None, None), ("not a unit", None)])
def test_parse_flux_unit(text, unit):
    assert parse_flux_unit(text) == unit


def _arb_session():
    x = WAVE_AA / 10
    spec = Spectrum(SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x - 0.01, au.nm),
        xmax=DataColumn(x + 0.01, au.nm), y=DataColumn(np.ones(50), ARB),
        dy=DataColumn(np.full(50, 0.1), ARB),
        aux_cols={'cont': DataColumn(np.full(50, 2.0), ARB),
                  'quality': DataColumn(np.zeros(50), au.dimensionless_unscaled)},
        z_em=2.5))
    return Session(name="s", gui=None, spec=spec)


def test_archive_round_trip_of_arb_flux(tmp_path):
    # No warning (they are errors in the tests): the unit is written in
    # TUNITn by hand and read back by its name
    path = tmp_path / "a.acs2"
    _arb_session().save(str(path))
    back = Session.open_new(str(path), "a", None, "auto")
    assert back.spec.y.unit == ARB and back.spec.dy.unit == ARB
    assert back.spec.get_column('cont').unit == ARB
    assert back.spec.get_column('quality').unit == au.dimensionless_unscaled


def test_archive_header_holds_arb_flux(tmp_path):
    import tarfile
    path = tmp_path / "a.acs2"
    _arb_session().save(str(path))
    with tarfile.open(path) as tar:
        tar.extractall(tmp_path / "x", filter="data")
    spec_file = next((tmp_path / "x").glob("*_spec.fits"))
    header = fits.getheader(spec_file, 1)
    units = {header[f"TTYPE{i}"]: header.get(f"TUNIT{i}") for i in range(1, header["TFIELDS"] + 1)}
    assert units['y'] == units['dy'] == units['cont'] == 'arb_flux'
    assert units['x'] == 'nm'


def test_old_archive_with_a_flux_without_unit(tmp_path):
    # An archive of a loader that gave a dimensionless flux: read in arb_flux
    # (the archive is not changed)
    session = _arb_session()
    session = session.with_new_spectrum(session.spec.with_flux_unit(au.dimensionless_unscaled))
    path = tmp_path / "a.acs2"
    session.save(str(path))
    back = Session.open_new(str(path), "a", None, "auto")
    assert back.spec.y.unit == ARB and back.spec.get_column('cont').unit == ARB
    assert back.spec.get_column('quality').unit == au.dimensionless_unscaled


# --- Conversion of the flux (F130) --------------------------------------------

def test_archive_continuum_without_unit_takes_the_unit_of_the_flux():
    # F130: cont and model without a unit, in an archive whose flux has one
    session = Session.open_new(str(J0103_SESSION), "J0103-1305", None, "auto")
    assert session.spec.get_column('cont').unit == FLAM
    assert session.spec.get_column('model').unit == FLAM


def test_y_convert_converts_the_continuum():
    # F130, before: cont without a unit was not converted with the flux, and
    # the normalized flux changed by a factor 10
    session = Session.open_new(str(J0103_SESSION), "J0103-1305", None, "auto")
    spec = session.spec
    out = session.edit.y_convert(yunit='erg/(nm s cm^2)').spec
    assert np.allclose(out.y.value, spec.y.value * 10)
    assert np.allclose(out.get_column('cont').value, spec.get_column('cont').value * 10)
    with np.errstate(divide='ignore', invalid='ignore'):
        assert np.allclose(out.y.value / out.get_column('cont').value,
                           spec.y.value / spec.get_column('cont').value, equal_nan=True)


def test_y_convert_of_a_continuum_without_unit_in_memory():
    spec = _session().spec
    data = spec._data
    aux = dict(data.aux_cols, cont=DataColumn(data.aux_cols['cont'].values,
                                              au.dimensionless_unscaled))
    import dataclasses
    spec = Spectrum(dataclasses.replace(data, aux_cols=aux))
    out = spec.y_convert(FLAM)
    assert out.get_column('cont').unit == FLAM
    assert np.allclose(out.get_column('cont').value, 0.25)


def test_y_convert_refuses_an_arbitrary_flux():
    with pytest.raises(InvalidParameter, match="arbitrary units.*set_properties"):
        _arb_session().edit.y_convert(yunit='erg/(nm s cm^2)')


# --- Declaring the unit (set_properties) --------------------------------------

def test_set_properties_declares_the_unit():
    session = _arb_session()
    out = session.edit.set_properties(yunit='1e-17 erg/(Angstrom s cm^2)')
    unit = 1e-17 * FLAM
    assert out.spec.y.unit == unit and out.spec.dy.unit == unit
    assert out.spec.get_column('cont').unit == unit
    assert out.spec.get_column('quality').unit == au.dimensionless_unscaled
    for name in ('y', 'dy'):
        assert np.array_equal(getattr(out.spec, name).value, getattr(session.spec, name).value)
    assert np.array_equal(out.spec.get_column('cont').value, session.spec.get_column('cont').value)
    # And back to arbitrary units
    back = out.edit.set_properties(yunit='arb_flux')
    assert back.spec.y.unit == ARB and back.spec.get_column('cont').unit == ARB


def test_set_properties_declares_columns_in_the_unit_of_the_flux():
    session = _session()
    out = session.edit.set_properties(yunit='erg/(Angstrom s cm^2)')
    assert out.spec.get_column('cont').unit == FLAM
    assert np.array_equal(out.spec.y.value, session.spec.y.value)


def test_set_properties_keeps_the_unit_by_default():
    session = _arb_session()
    assert session.edit.set_properties(z_em=3.0).spec.y.unit == ARB


@pytest.mark.parametrize("yunit,match", [
    ("bogus", "not a unit"), ("km/s", "not a flux density"), ("Jy", "not a flux density")])
def test_set_properties_refuses_other_units(yunit, match):
    with pytest.raises(InvalidParameter, match=match):
        _arb_session().edit.set_properties(yunit=yunit)


def test_set_properties_methods_text():
    from astrocook.core.structures import HistoryLog
    from astrocook.methods_text import methods_text
    from astrocook.session.history import SessionHistory
    from test_methods_text import step
    history = SessionHistory(_arb_session(), HistoryLog())
    step(history, "edit.set_properties", yunit='erg/(Angstrom s cm^2)')
    step(history, "edit.set_properties", yunit='arb_flux')
    text = methods_text(history).text
    assert "unit of the flux was declared as erg/(Angstrom s cm^2), the values unchanged" in text
    assert "The flux was declared in arbitrary units" in text


# --- CSV files (F129) ----------------------------------------------------------

def _session():
    x = WAVE_AA / 10
    y = np.linspace(1.0, 2.0, 50)
    spec = Spectrum(SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x - 0.01, au.nm),
        xmax=DataColumn(x + 0.01, au.nm), y=DataColumn(y, FLAM_NM),
        dy=DataColumn(y / 10, FLAM_NM),
        aux_cols={'cont': DataColumn(np.full(50, 2.5), FLAM_NM)}, z_em=2.5))
    return Session(name="s", gui=None, spec=spec)


def test_export_ascii_writes_the_names_of_the_gui_and_the_units(tmp_path):
    # Before: x,y,dy,cont and no units
    _session().file.export_ascii(targets='F, dF, cont', folder=str(tmp_path), prefix='p')
    text = (tmp_path / 'p_spec.csv').read_text()
    lines = text.splitlines()
    assert lines[0] == '# %ECSV 1.0'
    header = [line for line in lines if not line.startswith('#')][0]
    assert header == 'lam,F,dF,cont'
    assert '{name: F, unit: erg / (nm s cm2), datatype: float64}' in text
    assert 'z_em: 2.5' in text


@pytest.mark.parametrize("targets", ['F, dF, cont', 'spec'])
@pytest.mark.parametrize("yunit", ['erg/(nm s cm^2)', 'erg/(Angstrom s cm^2)', 'arb_flux'])
def test_csv_round_trip(tmp_path, yunit, targets):
    # F129, before: the flux came back in erg/(Angstrom s cm2), cont was
    # lost, and with 'spec' xmin came back as the flux
    session = _session().edit.set_properties(yunit=yunit)
    session.file.export_ascii(targets=targets, folder=str(tmp_path), prefix='p')
    back = load_session_from_file(str(tmp_path / 'p_spec.csv'), 'b', None, 'auto')
    assert back.spec.y.unit == session.spec.y.unit
    for name in ('x', 'y', 'dy'):
        assert np.array_equal(getattr(back.spec, name).value, getattr(session.spec, name).value)
    assert back.spec.get_column('cont').unit == session.spec.y.unit
    assert np.array_equal(back.spec.get_column('cont').value, session.spec.get_column('cont').value)
    assert back.spec._data.z_em == 2.5
    if targets == 'spec':
        assert np.array_equal(back.spec._data.xmin.values, session.spec._data.xmin.values)


def test_csv_round_trip_with_the_axis_in_velocity(tmp_path):
    session = _session().edit.x_convert(xunit='km/s')
    session.file.export_ascii(targets='spec', folder=str(tmp_path), prefix='p')
    back = load_session_from_file(str(tmp_path / 'p_spec.csv'), 'b', None, 'auto')
    assert back.spec.x.unit == au.km / au.s
    assert np.array_equal(back.spec.x.value, session.spec.x.value)


def test_old_csv_of_the_whole_spectrum(tmp_path):
    # Written by export_ascii before PR 8-bis: by name (before, by position,
    # so xmin came back as the flux); the flux in arbitrary units
    path = tmp_path / 'old.csv'
    path.write_text("x,xmin,xmax,y,dy,cont\n"
                    + "".join(f"{w / 10},{w / 10 - 0.01},{w / 10 + 0.01},1.5,0.1,2.0\n"
                              for w in WAVE_AA))
    back = load_session_from_file(str(path), 'b', None, 'auto')
    assert np.all(back.spec.y.value == 1.5) and back.spec.y.unit == ARB
    assert np.all(back.spec.get_column('cont').value == 2.0)
    assert back.spec.get_column('cont').unit == ARB
    assert back.spec.x.unit == au.nm


def test_import_ascii_reads_the_names_of_the_gui(tmp_path):
    # Before: a column named F was imported as an auxiliary column F (and lam
    # as a column lam)
    session = _session()
    (tmp_path / 'i.csv').write_text("lam,F\n" + "".join(f"{w / 10},3.0\n" for w in WAVE_AA))
    out = session.file.import_ascii(file_path=str(tmp_path / 'i.csv'))
    assert 'F' not in out.spec._data.aux_cols and 'lam' not in out.spec._data.aux_cols
    assert np.all(out.spec.y.value == 3.0) and out.spec.y.unit == FLAM_NM


def test_import_ascii_of_an_exported_file(tmp_path):
    session = _session()
    session.file.export_ascii(targets='spec', folder=str(tmp_path), prefix='p')
    out = session.file.import_ascii(file_path=str(tmp_path / 'p_spec.csv'))
    assert np.array_equal(out.spec.y.value, session.spec.y.value)
    assert out.spec.get_column('cont').unit == FLAM_NM


def test_import_ascii_refuses_a_flux_in_another_unit(tmp_path):
    session = _session()
    other = session.edit.y_convert(yunit='erg/(Angstrom s cm^2)')
    other.file.export_ascii(targets='F', folder=str(tmp_path), prefix='p')
    with pytest.raises(RecipeFailed, match="must be the same"):
        session.file.import_ascii(file_path=str(tmp_path / 'p_spec.csv'))


# --- Calibration --------------------------------------------------------------

def _wide_session(unit):
    x = np.linspace(400.0, 800.0, 2000)
    spec = Spectrum(SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x - 0.1, au.nm),
        xmax=DataColumn(x + 0.1, au.nm), y=DataColumn(np.ones(2000), unit),
        dy=DataColumn(np.full(2000, 0.1), unit)))
    return Session(name="s", gui=None, spec=spec)


@pytest.mark.parametrize("unit", [au.dimensionless_unscaled, ARB])
def test_calibration_reads_a_flux_without_unit_as_flam(caplog, unit):
    # An arbitrary flux as a dimensionless one (PR 8-bis: it failed in the
    # conversion to erg/(Angstrom s cm2))
    before = _wide_session(au.dimensionless_unscaled).flux.calibrate_from_magnitudes(
        magnitudes='SDSS_r=17.0')
    out = _wide_session(unit).flux.calibrate_from_magnitudes(magnitudes='SDSS_r=17.0')
    assert out.spec.y.unit == FLAM
    assert np.array_equal(out.spec.y.value, before.spec.y.value)
    assert "taken in erg/cm2/s" in caplog.text


# --- GUI and MCP --------------------------------------------------------------

def test_label_of_the_unit():
    from astrocook.gui.plot_style import flux_unit_label
    assert flux_unit_label(ARB) == "arbitrary units"
    assert flux_unit_label(FLAM) == "erg s\u207b\u00b9 cm\u207b\u00b2 \u00c5\u207b\u00b9"
