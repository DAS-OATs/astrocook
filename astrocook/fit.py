from astropy.table import Column, Table
import copy
import numpy as np

lya_x = 121.567

class Fit():

    def __init__(self, spec, lines,
                 id=None,
                 z=None,
                 logN=None,
                 b=None,
                 btur=None):

        self._xmin = copy.deepcopy(lines.xmin)
        self._xmax = copy.deepcopy(lines.xmax)
        self._x = copy.deepcopy(lines.x)
        self._y = copy.deepcopy(lines.y)
        self._t = copy.deepcopy(lines.t)
        self._spec = copy.deepcopy(spec.t)
        
        # Array with line IDs
        if (id == None):
            self._id = np.full(len(self._x), 'Ly_a', dtype=object)
        elif (len(id) == len(self._x)):
            self._id = np.array(id, dtype=object)
        else:
            print("IDs not recognized!")
        col_id = Column(self._id, name='ID')

        # Array with redshift
        if (z == None):
            lya_z = self._x / lya_x - 1.
            self._z = np.full(len(self._x), lya_z, dtype=float)
        elif (len(z) == len(self._x)):
            self._id = np.array(z, dtype=object)
        else:
            print("Redshifts not recognized!")
        col_z = Column(self._z, name='Z')
        
        # Array with column densities
        if (logN == None):
            self._logN = np.full(len(self._x), 14.0, dtype=float)
        elif (len(logN) == len(self._x)):
            self._logN = np.array(logN, dtype=object)
        else:
            print("Column densities not recognized!")
        col_logN = Column(self._logN, name='LOGN')

        # Array with Doppler broadenings
        if (b == None):
            self._b = np.full(len(self._x), 20.0, dtype=float)
        elif (len(b) == len(self._x)):
            self._b = np.array(b, dtype=object)
        else:
            print("Column densities not recognized!")
        col_b = Column(self._b, name='B')

        # Array with turbulence broadenings
        if (btur == None):
            self._btur = np.full(len(self._x), 0.0, dtype=float)
        elif (len(btur) == len(self._x)):
            self._btur = np.array(btur, dtype=object)
        else:
            print("Column densities not recognized!")
        col_btur = Column(self._btur, name='BTUR')
        
        col_add = Table(data=(col_id, col_z, col_logN, col_b, col_btur),
                        masked=True)
        self._t.add_columns(col_add.columns.values())
        
        self._use_good = lines.use_good

    def _mask(self, prop):
        _prop = getattr(self, prop)
        if self._use_good:
            ret = _prop[self._igood]
        else:
            ret = _prop
        return ret

    def _mask_col(self, col):
        if self._use_good:
            ret = self._t[col].quantity[self._igood]
        else:
            ret = self._t[col].quantity
        null = np.argwhere(self._t[col].mask)
        if null.size > 0:
            ret[null] = float('nan')
        return ret

    def comp(self, line=-1):
        """ Create an array of companion lines for each line in the list """
        iter = range(len(self._t))
        if line == -1:
            ret = np.array([self._t[self.comp_sel(l)] for l in iter])
        else:
            ret = self._t[self.comp_sel(line)]
        return ret
            
    def comp_sel(self, line):
        """ Define the selection indexes for companion lines """
        iter = range(len(self._t))
        self._t.sort('X')
        groups = np.append([0], np.cumsum(self._t['XMAX'][:-1] <
                                          self._t['XMIN'][1:]))
        return np.array([np.logical_and(groups[l] == groups[line],
                                        l != line) for l in iter])

    def group(self, line=-1):
        """ Create a group of line for each line in the list """
        iter = range(len(self._t))
        if line == -1:
            ret = np.array([self._t[self.group_sel(l)] for l in iter])
        else:
            ret = self._t[self.group_sel(line)]
        return ret

    def group_sel(self, line):
        """ Define the selection indexes for group lines """
        iter = range(len(self._t))
        self._t.sort('X')
        groups = np.append([0], np.cumsum(self._t['XMAX'][:-1] <
                                          self._t['XMIN'][1:]))
        return np.array([groups[l] == groups[line] for l in iter])

    def range(self, line=-1):
        """ Extract the spectral range for each line in the list, taking into
        account its companion lines """
        iter = range(len(self._t))
        if line == -1:
            ret = np.array([self._spec[self.range_sel(l)] for l in iter])
        else:
            ret = self._spec[self.range_sel(line)]
        return ret
            
    def range_sel(self, line):
        """ Define the selection indexes for spectral ranges """
        return np.logical_and(
            self._spec['X'] > min(self._t[self.group_sel(line)]['XMIN']),
            self._spec['X'] < max(self._t[self.group_sel(line)]['XMAX']))
    
    @property
    def lines(self):
        if self._use_good:
            return self._t[self._igood]
        else:
            return self._t
        
    @property
    def spec(self):
        if self._use_good:
            return self._spec[self._igood]
        else:
            return self._spec

    @property
    def t(self):
        if self._use_good:
            return self._t[self._igood]
        else:
            return self._t

    @property
    def x(self):
        """Quantities associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self._mask_col('X')

    @x.setter
    def x(self, value):
        if self._use_good:
            self._t['X'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['X'] = np.asarray(value, dtype='float')
        self._t['X'].unit = self._t['XMIN'].unit
    
    @property
    def y(self):
        """Quantities associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self._mask_col('Y')

    @y.setter
    def y(self, value):
        if self._use_good:
            self._t['Y'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['Y'] = np.asarray(value, dtype='float')
        self._t['Y'].unit = self._t['YMIN'].unit
