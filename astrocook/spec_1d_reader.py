from . import Spec1D, Spec1DCont
from astropy import units as u
from astropy.io import fits
import numpy as np


class Spec1DReader:
    def __init__(self, verbose=0):
        ''' Constructor for the Spec1DReader class. '''
        self._method = None
        self._source = None
        self._verbose = int(verbose)

    @property
    def source(self):
        """Source providing the spectrum."""
        return self._source
        
    def cont(self, filename):
        """Read a spectrum with continuum from an astrocook FITS file"""

        hdulist = fits.open(filename)
        data = hdulist[1].data
        meta = {}
        for (key, val) in hdulist[0].header.items():
            meta[key] = val
            
            
        names = np.array(hdulist[1].columns.names)
        units = np.array(hdulist[1].columns.units)
        x_unit = units[np.where(names == 'X')][0]
        y_unit = units[np.where(names == 'Y')][0]
        #print(np.where(names == 'X'), np.where(names == 'X'))
        #print(x_unit, y_unit)
        xmin = data['XMIN']
        xmax = data['XMAX']
        x = data['X']
        y = data['Y']            
        dy = data['DY']
        group = data['GROUP']
        resol = data['RESOL']
        abs_fit = data['ABS_FIT']
        em_fit = data['EM_FIT']
        abs_rem = data['ABS_REM']
        em_rem = data['EM_REM']
        cont = data['CONT']

        dx = 0.5 * (xmax - xmin)
        
        c1 = np.argwhere(y > 0)
        c2 = np.argwhere(dy > 0)
        igood = np.intersect1d(c1, c2)

        good = np.repeat(-1, len(x))
        good[igood] = 1

        gen = Spec1D(x, y, dy=dy, xmin=xmin, xmax=xmax, xUnit=x_unit, 
                     yUnit=y_unit, group=good, resol=resol, meta=meta)
        cont = Spec1DCont(gen, abs_fit=abs_fit, em_fit=em_fit, abs_rem=abs_rem, 
                          em_rem=em_rem, cont=cont)
                          
        return cont
        

    def sdss_dr10(self, filename):
        """Read a spectrum from a SDSS-DR10 FITS file."""
        hdulist = fits.open(filename)

        if (self._verbose > 0):
            print("spec1reader.sdss_dr10: reading from " + filename)
            hdulist.info()

        #Convert header from first extension into a dict
        meta = {}
        for (key, val) in hdulist[0].header.items():
            meta[key] = val

        x  = 10.**hdulist[1].data['loglam']
        y  = hdulist[1].data['flux']
        dy = np.repeat(float('nan'), len(x))

        c1 = np.argwhere(hdulist[1].data['and_mask'] == 0)
        c2 = np.argwhere(hdulist[1].data['ivar'] > 0)
        c3 = np.argwhere(hdulist[1].data['flux'] > 0)
        igood = np.intersect1d(c1, c2)
        igood = np.intersect1d(igood, c3)

        tmp = hdulist[1].data['ivar']
        dy[igood] = 1. / np.sqrt(tmp[igood])

        good = np.repeat(-1, len(x))
        good[igood] = 1

        s = Spec1D(x, y, dy=dy, 
                   xUnit=u.Angstrom, 
                   yUnit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom,
                   group=good, 
                   meta=meta)
        hdulist.close()

        #self._method = sdss_dr10 #TODO: check this
        self._source = filename
        return s


    def uves(self, filename):
        """Read a spectrum from a UVES file."""
        hdulist = fits.open(filename)

        if (self._verbose > 0):
            print("spec1reader.uves: reading from " + filename)
            hdulist.info()

        #Convert header from first extension into a dict
        meta = {}
        for (key, val) in hdulist[0].header.items():
            meta[key] = val

        data = hdulist[1].data
        x = data.field('WAVEL')
        dx = data.field('PIXSIZE')
        y = data.field('FLUX')            
        dy = data.field('FLUXERR')
        resol = [60000.] * len(data)
        resol_e = [1000.] * len(data)

        c1 = np.argwhere(hdulist[1].data['FLUX'] > 0)
        c2 = np.argwhere(hdulist[1].data['FLUXERR'] > 0)
        igood = np.intersect1d(c1, c2)

        good = np.repeat(-1, len(x))
        good[igood] = 1

        s = Spec1D(x, y, dy=dy, 
                   xUnit=u.nm, 
                   yUnit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom,
                   group=good, 
                   resol=resol,
                   meta=meta)
        hdulist.close()

        #self._method = sdss_dr10 #TODO: check this
        self._source = filename
        return s

    def simul(self, filename):
        """Read a simulated spectrum."""

        hdulist = fits.open(filename)

        if (self._verbose > 0):
            print("spec1reader.simul: reading from " + filename)
            hdulist.info()

        #Convert header from first extension into a dict
        meta = {}
        for (key, val) in hdulist[0].header.items():
            meta[key] = val

        data = hdulist[1].data
        x = data.field('WAVE') * 0.1
        xmin = np.append(
                   data.field('WAVE')[0],
                   0.5 * (data.field('WAVE')[:-1] + data.field('WAVE')[:-1]))
        xmax = np.append(
                   0.5 * (data.field('WAVE')[:-1] + data.field('WAVE')[:-1]),
                   data.field('WAVE')[len(data) - 1])
        xmin = xmin * 0.1
        xmax = xmax * 0.1
        y = data.field('NORMFLUX')            
        dy = data.field('STDEV')
        resol = data.field('WAVE') / data.field('FWHM')

        c1 = np.argwhere(hdulist[1].data['NORMFLUX'] > 0)
        c2 = np.argwhere(hdulist[1].data['STDEV'] > 0)
        igood = np.intersect1d(c1, c2)

        good = np.repeat(-1, len(x))
        good[igood] = 1

        s = Spec1D(x, y, dy=dy, 
                   xmin=xmin,
                   xmax=xmax,
                   xUnit=u.nm, 
                   yUnit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom,
                   group=good, 
                   resol=resol,
                   meta=meta)
        hdulist.close()

        #self._method = sdss_dr10 #TODO: check this
        self._source = filename
        return s
        
