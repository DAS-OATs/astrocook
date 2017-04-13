from astropy import units as u
from astropy.io import fits
from astropy.table import Column, Table
import copy
import numpy as np
import sys

class List():
    """ Class for generic line lists
    
    A generic line list is a QTable with the following columns: 
        -# @x: channel;
        -# @y: flux density in the channel;
        -# @dy: error on @y.
        -# @group: quality/grouping of the channel;

    The following metadata are also expected:
        -# @meta: miscellaneous information (TBD).        
    """

    def __init__(self, x, y, 
                 xmin=None,
                 xmax=None,
                 dx=None,
                 dy=None,
                 group=None,
                 xUnit=u.dimensionless_unscaled, 
                 yUnit=u.dimensionless_unscaled, 
                 meta=None,
                 dtype=float):
        """ Constructor for the LineList class.
        """

        col_x = Column(np.asarray(copy.deepcopy(x) , dtype=dtype), name='X')
        col_y = Column(np.asarray(copy.deepcopy(y) , dtype=dtype), name='Y')

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
        col_xmin = Column(np.asarray(copy.deepcopy(xmin), dtype=dtype), 
            name='XMIN')
        col_xmax = Column(np.asarray(copy.deepcopy(xmax), dtype=dtype), 
            name='XMAX')
        
        if (dy is None):
            dy = np.repeat(float('nan'), len(col_x))
        col_dy = Column(np.asarray(copy.deepcopy(dy), dtype=dtype), 
            name='DY')
        
        if (group is None):
            group = np.ones(len(col_x))
        col_group = Column(np.asarray(copy.deepcopy(group), dtype=int), 
            name='GROUP')

        # Auxiliary data and meta
        if (meta is None):
            meta = {}

        # Table creation
        self._t = Table(
            data=(col_xmin, col_xmax, col_x, col_y, col_dy, col_group), 
            masked=True, meta=meta)

        self._t['X'].unit = xUnit
        self._t['XMIN'].unit = xUnit
        self._t['XMAX'].unit = xUnit
        self._t['Y'].unit = yUnit
        self._t['DY'].unit = yUnit

        self._t['XMIN'].mask = np.isnan(self._t['XMIN'].quantity.value)
        self._t['XMAX'].mask = np.isnan(self._t['XMAX'].quantity.value)
        self._t['X'].mask = np.isnan(self._t['X'].quantity.value)
        self._t['Y'].mask = np.isnan(self._t['Y'].quantity.value)
        self._t['DY'].mask = np.isnan(self._t['DY'].quantity.value)

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
    def useGood(self):
        """Tells whether x, y, etc. getters return only data from channels flagged as good."""
        return self._useGood

    @useGood.setter
    def useGood(self, value):
        self._useGood = value
        if self._useGood:
            self._igood = np.argwhere(self._t['GROUP'] >= 0)

    @property
    def t(self):
        if self._useGood:
            return self._t[self._igood]
        else:
            return self._t

    @property
    def meta(self):
        """Meta information for the spectrum."""
        return self._t.meta
        
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
            
    def from_table(self, table, meta = {}):
        """Read a line list from a (list-like) table"""
        
        xmin = table['XMIN']
        xmax = table['XMAX']
        x = table['X']
        y = table['Y']            
        group = table['GROUP']
        
        igood = np.argwhere(y > 0)

        good = np.repeat(-1, len(x))
        good[igood] = 1

        list = List(x, y, xmin=xmin, xmax=xmax, xUnit=x.unit, yUnit=y.unit, 
                    group=good, meta=meta)
        return list

    def group(self):
        """Group the lines by fitting range
        
        Lines with overlapping fitting ranges [XMIN-XMAX] must be fitted 
        together. This function scans a list of lines and updates it so that 
        lines to be fitted together have the same GROUP (from 1 onward).
        """

        # Sort the line list by X, just in case
        self.t.sort('X')

        # Check whether the intervals [XMIN-XMAX] in two consecutive rows are 
        # disjoint
        disjoint = self.xmax[:-1] < self.xmin[1:] 

        # GROUP is updated by cumulatively summing its elements each time the
        # interval [XMIN-XMAX] in a row is disjoint from the previous one
        self.t['GROUP'][1:] = self.t['GROUP'][0] \
            + np.cumsum(self.t['GROUP'][1:] * disjoint)
        
    def guess_logN(self, cont):
        """Guess the column densities of absorption lines
        
        The column densities are guessed from line peaks, assuming a HI curve of
        growth with thermal broadening of 20 km/s. The curve of growth is 
        parametrized by a polynomial. 
        """
        
        a_0 =  14.096;
        a_1 = - 4.6251;
        a_2 =  18.657;
        a_3 = -46.299;
        a_4 =  53.301;
        a_5 = -23.442;
        return a_0 + a_1 *    (self.y.value / cont)    \
                   + a_2 * pow(self.y.value / cont, 2) \
                   + a_3 * pow(self.y.value / cont, 3) \
                   + a_4 * pow(self.y.value / cont, 4) \
                   + a_5 * pow(self.y.value / cont, 5)
        
   # return colden;

        

    def save(self, filename):
        hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='XMIN', format='E', array=self.xmin),
             fits.Column(name='XMAX', format='E', array=self.xmax),
             fits.Column(name='X', format='E', array=self.x),
             fits.Column(name='Y', format='E', array=self.y),
             fits.Column(name='DY', format='E', array=self.dy)])
        hdu.writeto(filename)