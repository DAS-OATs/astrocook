from .spectrum import Spectrum
from .line_list import LineList
from .syst_list import SystList
from astropy import units as au
import numpy as np

prefix = "Format:"

class Format(object):
    """ Class for file formats. """

    def __init__(self):
        pass

    def _create_xmin_xmax(self, x):
        mean = 0.5*(x[1:]+x[:-1])
        xmin = np.append(x[0], mean)
        xmax = np.append(mean, x[-1])
        return xmin, xmax

    def astrocook(self, hdul, struct):
        hdr = hdul[1].header
        data = hdul[1].data

        if struct in ['spec', 'nodes', 'lines']:

            x = data['x']
            xmin = data['xmin']
            xmax = data['xmax']
            y = data['y']
            dy = data['dy']
            xunit = au.nm
            yunit = au.erg/au.cm**2/au.s/au.Angstrom
            meta = hdr
            try:
                meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
            except:
                meta['object'] = ''
                print(prefix, "HIERARCH ESO OBS TARG NAME not defined.")
            if struct in ['spec', 'nodes']:
                return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)
            if struct in ['lines']:
                return LineList(x, xmin, xmax, y, dy, xunit, yunit, meta)

        if struct in ['systs']:
            series = data['series']
            func = data['func']
            z = data['z']
            logN = data['logN']
            b = data['b']
            chi2r = data['chi2r']
            id = data['id']
            return SystList(series, func, z, logN, b, chi2r, id)


    def eso_midas(self, hdul):
        """ ESO-MIDAS Table """

        hdr = hdul[1].header
        data = hdul[1].data
        x = data[hdr['TTYPE1']]
        xmin = data[hdr['TTYPE1']]
        xmax = data[hdr['TTYPE1']]
        y = data[hdr['TTYPE2']]
        dy = data[hdr['TTYPE3']]
        xunit = au.Angstrom
        yunit = au.erg/au.cm**2/au.s/au.Angstrom
        meta = {'instr': ''}
        try:
            meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
        except:
            meta['object'] = ''
            print(prefix, "HIERARCH ESO OBS TARG NAME not defined.")
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

    def espresso_das_spectrum(self, hdul):
        """ ESPRESSO DAS FSPEC/RSPEC format """

        hdr = hdul[0].header
        data = hdul[1].data
        x = data['WAVEL']
        xmin = x-data['PIXSIZE']*0.5
        xmax = x+data['PIXSIZE']*0.5
        y = data['FLUX']
        dy = data['FLUXERR']
        xunit = au.nm
        yunit = au.electron/au.nm #erg/au.cm**2/au.s/au.nm
        meta = {'instr': 'ESPRESSO'}
        try:
            meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
        except:
            meta['object'] = ''
            print(prefix, "HIERARCH ESO OBS TARG NAME not defined.")
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

    def espresso_drs_spectrum(self, hdul):
        """ ESPRESSO DRS S1D format """

        hdr = hdul[0].header
        data = hdul[1].data
        x = data['wavelength']
        xmin, xmax = self._create_xmin_xmax(x)
        y = data['flux']/(xmax-xmin)#*10#au.nm/au.Angstrom
        dy = data['error']/(xmax-xmin)#*10#au.nm/au.Angstrom
        xunit = au.Angstrom
        yunit = au.electron/au.Angstrom #erg/au.cm**2/au.s/au.nm
        meta = {'instr': 'ESPRESSO'}
        try:
            meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
        except:
            meta['object'] = ''
            print(prefix, "HIERARCH ESO OBS TARG NAME not defined.")
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

    def espresso_spectrum_format(self, data):
        """ ESPRESSO spectrum format """

        x = data['col2']
        xmin = data['col8']
        xmax = data['col9']
        y = np.zeros(len(x))
        dy = np.zeros(len(x))
        xunit = au.nm
        yunit = None
        meta = {}
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

    def uves_popler_spectrum(self, hdul):
        """ UVES POPLER format """

        hdr = hdul[0].header
        crval1 = hdr['CRVAL1']
        cdelt1 = hdr['CDELT1']
        naxis1 = hdr['NAXIS1']
        data = hdul[0].data
        y = data[:][0]*data[:][3]
        dy = data[:][1]*data[:][3]
        x = 10**np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)]
        xmin, xmax = self._create_xmin_xmax(x)
        xunit = au.Angstrom
        yunit = au.electron/au.Angstrom
        meta = {'instr': 'UVES'}
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
            print(prefix, "OBJECT not defined.")
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

    def xshooter_reduce_spectrum(self, hdul, hdul_e):
        hdr = hdul[0].header
        crval1 = hdr['CRVAL1']
        cdelt1 = hdr['CDELT1']
        naxis1 = hdr['NAXIS1']
        data = hdul[0].data
        data_e = hdul_e[0].data
        y = data
        dy = data_e
        x = 10**np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)]
        xmin, xmax = self._create_xmin_xmax(x)
        xunit = au.Angstrom
        yunit = au.electron/au.Angstrom
        meta = {'instr': 'X-shooter'}
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
            print(prefix, "OBJECT not defined.")
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)
