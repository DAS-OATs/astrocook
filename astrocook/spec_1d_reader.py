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
    
    def ac(self, filename):
        """Read a spectrum from an astrocook FITS file"""

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

        dx = 0.5 * (xmax - xmin)
        
        c1 = np.argwhere(y > 0)
        c2 = np.argwhere(dy > 0)
        igood = np.intersect1d(c1, c2)

        good = np.repeat(-1, len(x))
        good[igood] = 1

        gen = Spec1D(x, y, dy=dy, xmin=xmin, xmax=xmax, xunit=x_unit, 
                     yunit=y_unit, group=good, resol=resol, meta=meta)
        return gen
    
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

        gen = Spec1D(x, y, dy=dy, xmin=xmin, xmax=xmax, xunit=x_unit, 
                     yunit=y_unit, group=good, resol=resol, meta=meta)
        cont = Spec1DCont(gen, abs_fit=abs_fit, em_fit=em_fit, abs_rem=abs_rem, 
                          em_rem=em_rem, cont=cont)
                          
        return cont
    
    def syst(self, name):

        hdulist = fits.open(name + '_syst_spec.fits')
        data = hdulist[1].data
        meta = {}
        for (key, val) in hdulist[0].header.items():
            meta[key] = val
            
        names = np.array(hdulist[1].columns.names)
        units = np.array(hdulist[1].columns.units)
        x_unit = units[np.where(names == 'X')][0]
        y_unit = units[np.where(names == 'Y')][0]

        # Currently units are hardcoded
        x_unit = u.nm
        y_unit = 1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom
        #print(np.where(names == 'X'), np.where(names == 'X'))
        #print(x_unit, y_unit)
        xmin = data['XMIN']
        xmax = data['XMAX']
        x = data['X']
        y = data['Y']            
        y_fit = data['Y_FIT']
        y_rem = data['Y_REM']
        dy = data['DY']
        group = data['GROUP']
        resol = data['RESOL']
        dx = 0.5 * (xmax - xmin)
        
        c1 = np.argwhere(y > 0)
        c2 = np.argwhere(dy > 0)
        igood = np.intersect1d(c1, c2)

        good = np.repeat(-1, len(x))
        good[igood] = 1

        spec = Spec1D(x, y, dy=dy, xmin=xmin, xmax=xmax, xunit=x_unit, 
                     yunit=y_unit, group=good, resol=resol, meta=meta)
        fit = Spec1D(x, y=y_fit, dy=dy, xmin=xmin, xmax=xmax, xunit=x_unit, 
                     yunit=y_unit, group=good, resol=resol, meta=meta)
        rem = Spec1D(x, y=y_rem, dy=dy, xmin=xmin, xmax=xmax, xunit=x_unit, 
                     yunit=y_unit, group=good, resol=resol, meta=meta)
        
        return (spec, fit, rem)
        

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
                   xunit=u.Angstrom, 
                   yunit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom,
                   group=good, 
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
                   xunit=u.nm, 
                   yunit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom,
                   group=good, 
                   resol=resol,
                   meta=meta)
        hdulist.close()

        #self._method = sdss_dr10 #TODO: check this
        self._source = filename
        return s

    def uves(self, filename, resol=None):
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

        
        
        try:
            cond = ~np.isnan(data.field('WAVEL'))
            x = data.field('WAVEL')[cond]
            x_name = 'WAVEL'
        except:
            cond = ~np.isnan(data.field('wave'))            
            x = data.field('wave')[cond] * 0.1
            x_name = 'wave'
        try:
            dx = data.field('PIXSIZE')[cond]
        except:
            dx = data.field('wpix')[cond] * 0.1
        try:
            y = data.field('FLUX')[cond]
            y_name = 'FLUX'
        except:
            y = data.field('flux')[cond]
            y_name = 'flux'
        try:
            dy = data.field('FLUXERR')[cond]
            dy_name = 'FLUXERR'
        except:
            dy = data.field('sigma')[cond]
            dy_name = 'sigma'
        #resol = [60000.] * len(data)
        #resol_e = [1000.] * len(data)
        resol_arr = [resol] * len(x)
        
        c1 = np.argwhere(hdulist[1].data[y_name] > 0)
        c2 = np.argwhere(hdulist[1].data[dy_name] > 0)
        igood = np.intersect1d(c1, c2)
        
        good = np.repeat(-1, len(x))
        good[igood] = 1

        s = Spec1D(x, y, dy=dy, 
                   xunit=u.nm, 
                   yunit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom,
                   group=good, 
                   resol=resol_arr,
                   meta=meta)
        hdulist.close()

        #self._method = sdss_dr10 #TODO: check this
        self._source = filename
        return s

    def xsh(self, filename, resol=None):
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
            
        try:
            x = data.field('WAVEL')
        except:
            x = data.field('wave') * 0.1
        try:
            dx = data.field('PIXSIZE')
        except:
            dx = data.field('pixsize') * 0.1
        try:
            y = data.field('FLUX')
            y_name = 'FLUX'
        except:
            y = data.field('flux')
            y_name = 'flux'
        try:
            dy = data.field('FLUXERR')
            dy_name = 'FLUXERR'
        except:
            dy = data.field('err')
            dy_name = 'err'
        try:
            resol_arr = data.field('resol')
            resol_name = 'resol'
        except:
            resol_arr = [resol] * len(data)

        c1 = np.argwhere(hdulist[1].data[y_name] > 0)
        c2 = np.argwhere(hdulist[1].data[dy_name] > 0)
        igood = np.intersect1d(c1, c2)

        good = np.repeat(-1, len(x))
        good[igood] = 1

        s = Spec1D(x, y, dy=dy, 
                   xunit=u.nm, 
                   yunit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom,
                   group=good, 
                   resol=resol_arr,
                   meta=meta)
        hdulist.close()

        self._source = filename
        return s
    
    def xq(self, filename, resol=None):
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
        
        x = data.field('wave')*0.1
        dx = data.field('wpix')*0.1
       
        y = data.field('flux')
        
        y_name = 'flux'
        dy = data.field('err')
        dy_name = 'err'
        try:
            resol_arr = data.field('resolution')
        except:
            resol_arr = data.field('resol')            
        resol_name = 'resol'


        resol_arr[resol_arr==5100] = 5500
        resol_arr[resol_arr==8800] = 8900

        c1 = np.argwhere(hdulist[1].data[y_name] > 0)
        c2 = np.argwhere(hdulist[1].data[dy_name] > 0)
        igood = np.intersect1d(c1, c2)

        good = np.repeat(-1, len(x))
        good[igood] = 1

        s = Spec1D(x, y, dy=dy, 
                   xunit=u.nm, 
                   yunit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom,
                   group=good, 
                   resol=resol_arr,
                   meta=meta)
        hdulist.close()

        self._source = filename
        return s
    
