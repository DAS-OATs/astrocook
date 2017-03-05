import numpy as np
from astropy.io import fits as fits
from astropy import units as u
from astropy.table import Column, QTable, Table
from specutils import extinction
import copy


class spec1d():
    """Class for generic spectra
    
    A generic spectrum is a QTable with the following columns: 
        -# @exp: exposure;
        -# @order: spectral order;
        -# @x: channels;
        -# @xmin: lower limit for each channel; 
        -# @xmax: upper limit for each channel; 
        -# @y: flux density in the channel;
        -# @dy: error on @flux.
        -# @group: quality/grouping of the channel;
        -# @resol: spectral resolution in the channel; 
        
    The following metadata are also expected:
        -# @order: spectrum order;
        -# @exptime: integration time of all exposures;
        -# @meta: miscellaneous information (TBD).        
    """
    
    def __init__(self, x, y, 
                 xmin=None,
                 xmax=None,
                 dx=None,
                 dy=None,
                 group=None,
                 resol=None,
                 xUnit=u.dimensionless_unscaled, 
                 yUnit=u.dimensionless_unscaled, 
                 exptime=float('nan'), 
                 order=-1,
                 meta=None,
                 dtype=float):
        ''' Constructor for the spec1d class. '''

        col_x  = Column(np.asarray(copy.deepcopy(x) , dtype=float), name='X')
        col_y  = Column(np.asarray(copy.deepcopy(y) , dtype=float), name='Y')

        if (xmin is None   and   xmax is not None):
            raise Exception('When XMAX is used also XMIN must be given')
        if (xmin is not None   and   xmax is None):
            raise Exception('When XMIN is used also XMAX must be given')
        if (xmin is None   and   xmax is None):
            if (dx is None): 
                dx = np.roll((x - np.roll(x, 1))/2., -1)
                dx[-1] = dx[-2]
            xmin = x - dx
            xmax = x + dx
        col_xmin = Column(np.asarray(copy.deepcopy(xmin), dtype=float), name='XMIN')
        col_xmax = Column(np.asarray(copy.deepcopy(xmax), dtype=float), name='XMAX')

        if (dy is None):
            dy = np.repeat(float('nan'), len(col_x))
        col_dy = Column(np.asarray(copy.deepcopy(dy), dtype=float), name='DY')
        
        if (group is None):
            group = np.ones(len(col_x))
        col_g = Column(np.asarray(copy.deepcopy(group), dtype=int), name='GROUP')
        
        if (resol is None):
            resol = np.repeat(float('nan'), len(col_x))
        col_r = Column(np.asarray(copy.deepcopy(resol), dtype=float), name='RESOL')
        
        # Auxiliary data and meta
        self._exptime = float(exptime)
        self._order   = int(order)
        if (meta is None):
            meta = {}

        # Table creation
        self._t = Table(data=(col_xmin, col_xmax, col_x, col_y, col_dy, col_g, col_r), masked=True, meta=meta)
        self._t['XMIN'].unit = xUnit
        self._t['XMAX'].unit = xUnit
        self._t['X'].unit    = xUnit
        self._t['Y'].unit    = yUnit
        self._t['DY'].unit   = yUnit

        self._t['XMIN'].mask  = np.isnan(self._t['XMIN'].quantity.value)
        self._t['XMAX'].mask  = np.isnan(self._t['XMAX'].quantity.value)
        self._t['X'].mask     = np.isnan(self._t['X'].quantity.value)
        self._t['Y'].mask     = np.isnan(self._t['Y'].quantity.value)
        self._t['DY'].mask    = np.isnan(self._t['DY'].quantity.value)
        self._t['RESOL'].mask = np.isnan(self._t['RESOL'].quantity.value)

        self._useGood = False


    def _getWithMask(self,colName):
        if self._useGood:
            ret = self._t[colName].quantity[self._igood]
        else:
            ret = self._t[colName].quantity
        null = np.argwhere(self._t[colName].mask)
        if null.size > 0:
            ret[null] = float('nan')
        return ret

    @property
    def t(self):
        if self._useGood:
            return self._t[self._igood]
        else:
            return self._t

    @property
    def useGood(self):
        """Tells whether x, y, etc. getters return only data from channels flagged as good."""
        return self._useGood

    @useGood.setter
    def useGood(self, value):
        self._useGood = value
        if self._useGood:
            self._igood = np.argwhere(self._t['GROUP'] >= 0)

    @property
    def exptime(self):
        """Spectrum exposure time (in seconds)."""
        return self._exptime

    @property
    def order(self):
        """Spectrum order (-1 if not specified)."""
        return self._order

    @property
    def meta(self):
        """Meta information for the spectrum."""
        return self._t.meta

    def nchan(self):
        """Number of channels in the spectrum."""
        if self._useGood:
            return len(self._igood)
        else:
            return len(self._t)

    @property
    def x(self):
        """Quantities associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self._getWithMask('X')

    @x.setter
    def x(self, value):
        if self._useGood:
            self._t['X'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['X'] = np.asarray(value, dtype='float')
        self._t['X'].unit = self._t['XMIN'].unit

    @property
    def xmin(self):
        """Lower limit of quantity associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self._getWithMask('XMIN')

    @xmin.setter
    def xmin(self, value):
        if self._useGood:
            self._t['XMIN'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['XMIN'] = np.asarray(value, dtype='float')
        self._t['XMIN'].unit = self._t['X'].unit

    @property
    def xmax(self):
        """Upper limit of quantity associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self._getWithMask('XMAX')

    @xmax.setter
    def xmax(self, value):
        if self._useGood:
            self._t['XMAX'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['XMAX'] = np.asarray(value, dtype='float')
        self._t['XMAX'].unit = self._t['X'].unit

    @property
    def dx(self):
        """Widths of spectrum channels, in the same units as the spectrum channels."""
        xmin = self._getWithMask('XMIN')
        xmax = self._getWithMask('XMAX')
        return xmax - xmin

    @property
    def y(self):
        """Quantities associated to spectrum intensities (e.g. flux density, luminosity density, nuF_nu, lambdaF_lambda, etc.)."""
        return self._getWithMask('Y')

    @y.setter
    def y(self, value):
        if self._useGood:
            self._t['Y'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['Y'] = np.asarray(value, dtype='float')
        self._t['Y'].unit = self._t['DY'].unit

    @property
    def dy(self):
        """Uncertainties associated to y values."""
        return self._getWithMask('DY')

    @dy.setter
    def dy(self, value):
        if self._useGood:
            self._t['DY'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['DY'] = np.asarray(value, dtype='float')
        self._t['DY'].unit = self._t['Y'].unit

    @property
    def group(self):
        """Return group flag for each spectrum channel."""
        if self._useGood:
            return self._t['GROUP'].quantity.value[self._igood]
        else:
            return self._t['GROUP'].quantity.value

    @group.setter
    def group(self, value):
        if self._useGood:
            self._t['GROUP'][self._iGood] = np.asarray(value, dtype='int')
        else:
            self._t['GROUP'] = np.asarray(value, dtype='int')

    @property
    def resol(self):
        """Return resolution for each spectrum channel."""
        return self._getWithMask('RESOL')

    @property
    def xUnit(self):
        """Physical unit for the x property, to be expressed as an astropy unit."""
        return self._t['X'].unit

    @property
    def yUnit(self):
        """Physical unit for the y property, to be expressed as an astropy unit."""
        return self._t['Y'].unit

    def convert(self, xUnit=None, yUnit=None):
        """Convert x and/or y values into equivalent quantities."""
        if not (xUnit is None):
            mask = self._t['X'].mask
            q = self._t['X']
            p = q.to(xUnit, equivalencies=u.spectral())
            self._t['X'] = p
            self._t['X'].mask = mask

            mask = self._t['XMIN'].mask
            q = self._t['XMIN']
            p = q.to(xUnit, equivalencies=u.spectral())
            self._t['XMIN'] = p
            self._t['XMIN'].mask = mask

            mask = self._t['XMAX'].mask
            q = self._t['XMAX']
            p = q.to(xUnit, equivalencies=u.spectral())
            self._t['XMAX'] = p
            self._t['XMAX'].mask = mask

        if not (yUnit is None):
            mask = self._t['Y'].mask
            q = self._t['Y']
            p = q.to(yUnit, equivalencies=u.spectral_density(self._t['X']))
            self._t['Y'] = p
            self._t['Y'].mask = mask

            mask = self._t['DY']
            q = self._t['DY']
            p = q.to(yUnit, equivalencies=u.spectral_density(self._t['X']))
            self._t['DY'] = p
            self._t['DY'].mask = mask


    def deredden(self, A_v, model='od94'):
        extFactor = extinction.reddening(self._t['X'], A_v, model=model)
        self._t['Y']  *= extFactor
        self._t['DY'] *= extFactor
