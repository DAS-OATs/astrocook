from . import List
from astropy import units as u
from astropy.io import fits
from astropy.table import Column, Table
import copy
import numpy as np
import sys

class ListSyst(List):
    """ Class for line lists of absorption systems
    
    A line list of absorption systems is a generic line list with the following 
    additional columns: 
        -# @id: identification of the system;
        -# @name: name of the transition;
        -# @z: redshift
    """

    def __init__(self, list=None, 
                 syst=None,
                 z=None,
                 name=None, 
                 meta=None,
                 dtype=float):

        # Column definition
        if (syst is None):
            syst = np.repeat(0, len(list.x))
        col_syst = Column(np.asarray(copy.deepcopy(syst), dtype=int), 
                          name='SYST')

        if (name is None):
            name = np.repeat(None, len(list.x))
        col_name = Column(np.asarray(copy.deepcopy(name), dtype=str), 
                          name='NAME')

        if (z is None):
            z = np.repeat(float('nan'), len(list.x))
        col_z = Column(np.asarray(copy.deepcopy(z), dtype=dtype), name='Z')

        # Auxiliary data and meta
        if (meta is None):
            meta = {}

        # Table creation
        self._t = copy.deepcopy(list.t)
        self._t.add_columns([col_syst, col_name, col_z])

        self._useGood = False
        
    @property
    def syst(self):
        """Quantities associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self._getWithMask('SYST')

    @syst.setter
    def syst(self, value):
        if self._useGood:
            self._t['SYST'][self._iGood] = np.asarray(value, dtype='int')
        else:
            self._t['SYST'] = np.asarray(value, dtype='int')
        self._t['SYST'].unit = self._t['SYST'].unit

    @property
    def name(self):
        """Quantities associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self._getWithMask('NAME')

    @name.setter
    def name(self, value):
        if self._useGood:
            self._t['NAME'][self._iGood] = np.asarray(value, dtype='string')
        else:
            self._t['NAME'] = np.asarray(value, dtype='str')
        self._t['NAME'].unit = self._t['NAME'].unit
        
    @property
    def z(self):
        """Quantities associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self._getWithMask('SYST')

    @z.setter
    def z(self, value):
        if self._useGood:
            self._t['Z'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['Z'] = np.asarray(value, dtype='float')
        self._t['Z'].unit = self._t['Z'].unit
        
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
            
        for i in range(0, len(self.t)):
            syst_coin = self.t['SYST'] == self.t['SYST'][i] 
            self.t['GROUP'][syst_coin] = self.t['GROUP'][i]

        # Sort the line list by X, just in case
        self.t.sort('GROUP')

        #print self.t
