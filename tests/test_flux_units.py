"""The units of the flux and the CSV files (phase 6-quater, PR 8-bis, with
F129 and F130; agreed with Guido, 12 November 2026).

Characterization first: the unit that each loader assigns to the flux, the
conversion of the flux, the CSV files of ``export_ascii`` and
``import_ascii``, before the change; each test says what it was. No value
changes: only units.
"""
import astropy.units as au
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, SpectrumData
from astrocook.io import loaders
from astrocook.session.session import Session, load_session_from_file

from conftest import J0103_SESSION

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


# Before PR 8-bis: the unit of the flux of each loader, whatever the file said
BEFORE = {
    'astrocook_fits': FLAM,             # TUNIT not read
    'pypeit_coadd': au.dimensionless_unscaled,
    'pypeit_spec1d_multi': au.dimensionless_unscaled,
    'pypeit_stack': au.dimensionless_unscaled,
    'pypeit_tellmodel': au.dimensionless_unscaled,
    'espresso_s2d_fits': au.dimensionless_unscaled,
    'xshooter_fits': FLAM,              # TUNIT not read
    'generic_spectrum': au.dimensionless_unscaled,
    'simple_csv': FLAM,                 # assumed
    'ascii_telluric_model': FLAM,       # assumed
    'ascii_resvel_header': au.dimensionless_unscaled,
}


@pytest.mark.parametrize("fmt", sorted(FORMATS))
def test_unit_of_the_flux_of_each_loader(tmp_path, fmt):
    data = _load(tmp_path, fmt)
    assert data.y.unit == BEFORE[fmt] and data.dy.unit == BEFORE[fmt]


# --- Conversion of the flux (F130) --------------------------------------------

def test_y_convert_leaves_a_continuum_without_unit():
    # F130, before: cont without a unit in an archive whose flux has one is
    # not converted with it
    session = Session.open_new(str(J0103_SESSION), "J0103-1305", None, "auto")
    spec = session.spec
    assert spec.get_column('cont').unit == au.dimensionless_unscaled
    out = session.edit.y_convert(yunit='erg/(nm s cm^2)').spec
    assert np.allclose(out.y.value, spec.y.value * 10)
    assert np.array_equal(out.get_column('cont').value, spec.get_column('cont').value)


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


def test_export_ascii_writes_the_names_of_the_data(tmp_path):
    # Before: x, y, dy and no units
    _session().file.export_ascii(targets='F, dF, cont', folder=str(tmp_path), prefix='p')
    lines = (tmp_path / 'p_spec.csv').read_text().splitlines()
    assert lines[0] == 'x,y,dy,cont'
    assert float(lines[1].split(',')[0]) == pytest.approx(400.0)


def test_csv_round_trip_loses_units_and_cont(tmp_path):
    # F129, before: the flux comes back in erg/(Angstrom s cm2), cont lost
    session = _session().edit.y_convert(yunit='erg/(nm s cm^2)')
    session.file.export_ascii(targets='F, dF, cont', folder=str(tmp_path), prefix='p')
    back = load_session_from_file(str(tmp_path / 'p_spec.csv'), 'b', None, 'auto')
    assert back.spec.y.unit == FLAM
    assert np.allclose(back.spec.y.value, session.spec.y.value)
    assert 'cont' not in back.spec._data.aux_cols


def test_csv_of_the_whole_spectrum_comes_back_wrong(tmp_path):
    # F129, before: 'spec' writes x, xmin, xmax, y, dy, ...; the loader takes
    # the first three columns, so xmin comes back as the flux
    session = _session()
    session.file.export_ascii(targets='spec', folder=str(tmp_path), prefix='p')
    header = (tmp_path / 'p_spec.csv').read_text().splitlines()[0]
    assert header == 'x,xmin,xmax,y,dy,cont'
    back = load_session_from_file(str(tmp_path / 'p_spec.csv'), 'b', None, 'auto')
    assert np.allclose(back.spec.y.value, session.spec._data.xmin.values)


def test_import_ascii_makes_a_column_F(tmp_path):
    # Before: a column named F is imported as an auxiliary column F
    session = _session()
    (tmp_path / 'i.csv').write_text("lam,F\n" + "".join(f"{w / 10},3.0\n" for w in WAVE_AA))
    out = session.file.import_ascii(file_path=str(tmp_path / 'i.csv'))
    assert 'F' in out.spec._data.aux_cols and 'lam' in out.spec._data.aux_cols
    assert np.array_equal(out.spec.y.value, session.spec.y.value)


# --- Calibration --------------------------------------------------------------

def _wide_session(unit):
    x = np.linspace(400.0, 800.0, 2000)
    spec = Spectrum(SpectrumData(
        x=DataColumn(x, au.nm), xmin=DataColumn(x - 0.1, au.nm),
        xmax=DataColumn(x + 0.1, au.nm), y=DataColumn(np.ones(2000), unit),
        dy=DataColumn(np.full(2000, 0.1), unit)))
    return Session(name="s", gui=None, spec=spec)


def test_calibration_reads_a_flux_without_unit_as_flam(caplog):
    out = _wide_session(au.dimensionless_unscaled).flux.calibrate_from_magnitudes(
        magnitudes='SDSS_r=17.0')
    assert out.spec.y.unit == FLAM
    assert "taken in erg/cm2/s" in caplog.text
