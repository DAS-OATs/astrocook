import numpy as np
from astropy.io import fits
from astropy import units as u
from . import spec1d


class spec1dreader:
    def __init__(self, verbose=0):
        ''' Constructor for the spec1dreader class. '''
        self._method = None
        self._source = None
        self._verbose = int(verbose)

    @property
    def source(self):
        """Source providing the spectrum."""
        return self._source

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

        s = spec1d(x, y, dy=dy, 
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


        data = fits.open(file)[1].data
        x = data.field('WAVEL')
        dx = data.field('PIXSIZE')
        y = data.field('FLUX')            
        dy = data.field('FLUXERR')
        resol = [60000.] * len(data)
        resol_e = [1000.] * len(data)



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

        s = spec1d(x, y, dy=dy, 
                   xUnit=u.Angstrom, 
                   yUnit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom,
                   group=good, 
                   meta=meta)
        hdulist.close()

        #self._method = sdss_dr10 #TODO: check this
        self._source = filename
        return s






