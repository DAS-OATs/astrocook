import numpy as np
from astropy.io import fits
from astropy import units as u
from . import List


class ListReader:
    def __init__(self, verbose=0):
        ''' Constructor for the Spec1DReader class. '''
        self._method = None
        self._source = None
        self._verbose = int(verbose)

    @property
    def source(self):
        """Source providing the spectrum."""
        return self._source

    def gen(self, filename):
        """Read a simulated spectrum."""

        hdulist = fits.open(filename)

        if (self._verbose > 0):
            print("spec1reader.ac: reading from " + filename)
            hdulist.info()

        #Convert header from first extension into a dict
        meta = {}
        for (key, val) in hdulist[0].header.items():
            meta[key] = val

        data = hdulist[1].data
        x = data.field('X')
        xmin = data.field('XMIN')
        xmax = data.field('XMAX')
        y = data.field('Y')            
        dy = data.field('DY')

        """
        c1 = np.argwhere(hdulist[1].data['Y'] > 0)
        c2 = np.argwhere(hdulist[1].data['DY'] > 0)
        igood = np.intersect1d(c1, c2)
        """
        
        good = np.repeat(1, len(x))
        #good[igood] = 1

        gen = List(x, y, dy=dy, 
                 xmin=xmin,
                 xmax=xmax,
                 xUnit=u.nm, 
                 yUnit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom,
                 group=good)
        hdulist.close()

        #self._method = sdss_dr10 #TODO: check this
        self._source = filename
        return gen