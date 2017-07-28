from . import Line
from astropy.table import Column, Table
from copy import deepcopy as dc
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
import numpy as np
import warnings
import sys

class Syst(Line):

    def __init__(self,
                 line=None,
                 spec=None,
                 x=[],
                 y=[],
                 xmin=[],
                 xmax=[],
                 dy=[],
                 ion=[],
                 yunit=None,
                 meta=None,
                 dtype=float):
        """ Constructor for the Syst class """ 
        
        # Exceptions
        #if ((ion is not []) and (line is not None) \
        #    and (len(ion) != len(line.t))):
        #    raise Exception("Ion list and line list must have the same length.")
        if ((x is []) != (y is [])):
            raise Exception("X and Y must be provided together.")
        if ((xmin is []) != (xmax is [])):
            raise Exception("XMIN and XMAX must be provided together.")
        sumlen = len(x) + len(y) + len(xmin) + len(xmax) + len(dy) + len(ion)
        if ((x != []) and (sumlen % len(x) != 0)):
            raise Exception("Data arrays must have the same length.")
        if ((xmin != []) and (sumlen % len(xmin) != 0)):
            raise Exception("Data arrays must have the same length.")

        # Warnings
        if ((spec is None) and (line is None) and (x is [])):
            warnings.warn("No spectrum, line list, or data provided.")

        if ((line is not None) and (ion is [])):
            ion = np.full(len(line._t), 'Ly_a')           

        if ((x != []) and (ion is [])):
            ion = np.full(len(x), 'Ly_a')           

        # Line list
        self._line = None
        if (line is not None):
            self._line = dc(line)

        # Spectrum
        self._spec = None
        if (spec is not None):
            self._spec = dc(spec)
            
        # Redshift list
        self._redsh = None
        if ((ion is not []) and (line is not None) \
            and (len(ion) == len(line.t))):
            self._redsh = dc(line)
            self._redsh.to_z(ion)
            self._redsh.t.add_column(Column(np.asarray(dc(ion)), name='ION'))

        # System list
        
        data = ()
        if (x != []):
            col_x = Column(np.asarray(dc(x), dtype=dtype), name='X')
            col_y = Column(np.asarray(dc(y)), name='Y')
            data = (col_x, col_y)
            if yunit is None:
                try:
                    yunit = y.unit
                except:
                    raise Exception("Y unit not provided.")
        if (xmin != []):
            col_xmin = Column(np.asarray(dc(xmin), dtype=dtype), name='XMIN')
            col_xmax = Column(np.asarray(dc(xmax), dtype=dtype), name='XMAX')
            data += (col_xmin, col_xmax)
        if (dy != []):
            col_dy = Column(np.asarray(dc(dy), dtype=dtype), name='DY')
            data += (col_dy,)
        if ((x != []) and (ion != [])):
            col_ion = Column(np.asarray(dc(ion)), name='ION')
            data += (col_ion,)
        if data is ():
            data = None
            
        if (meta is None):
            meta = {}

        self._t = Table(data=data, masked=True, meta=meta)
        if (y != []):
            self._t['Y'].unit = yunit
        if (dy != []):
            self._t['DY'].unit = yunit

        self._use_good = False

        if (self._redsh is not None):
            if self._use_good:
                self._redsh.ion = np.asarray(self._redsh.t['ION'][self._igood])
            else:
                self._redsh.ion = np.asarray(self._redsh.t['ION'])

                
# Properties

    @property
    def ion(self):
        if self._use_good:
            ret = np.asarray(self._t['ION'][self._igood])
        else:
            ret = np.asarray(self._t['ION'])
        return ret 

    @ion.setter
    def ion(self, value):
        if self._use_good:
            self._t['ION'][self._igood] = np.asarray(value)
        else:
            self._t['ION'] = np.asarray(value)
        
    @property
    def line(self):
        return self._line

    @line.setter
    def line(self, value):
        if isinstance(value, Line):
            self._line = value
        else:
            raise Exception("Line list has a wrong format.")
    
    @property
    def redsh(self):
        return self._redsh

    @redsh.setter
    def redsh(self, value):
        if isinstance(value, Syst):
            self._redsh = value
        else:
            raise Exception("Redshift list has a wrong format.")


# Methods

    def plot(self, figsize=(6,6), split=False, block=True, **kwargs):
        ion = np.unique(self.ion)
        n = len(ion)
        z = self.x
        zmin = np.min(self.xmin)
        zmax = np.max(self.xmax)
        fig = plt.figure(figsize=figsize)
        if split == True:
            r = min(n,4)
            c = int(np.ceil(n/4))
            figsize = (n*2, c*4)
            grid = gs(r,c)            
            for p in range(n):
                ax = fig.add_subplot(grid[p%4, int(np.floor(p/4))])
                ax.set_ylabel("Flux [" + str(self._spec.y.unit) + "]")
                spec = dc(self._spec)
                line = dc(self._line)
                spec.to_z([ion[p]])
                line.to_z([ion[p]])
                chunk = spec.t[np.logical_and(spec.x > zmin, spec.x < zmax)]
                memb = line.t[np.logical_and(line.x > zmin, line.x < zmax)]
                ax.set_xlim(zmin, zmax)
                ax.plot(chunk['X'], chunk['Y'])
                ax.scatter(memb['X'], memb['Y'])
                for comp in z:
                    ax.axvline(x=comp, ymin=0.65, ymax=0.85, color='lightgray')
                if (p+1) % r != 0:
                    ax.set_xticks([], [])
                else:
                    ax.set_xlabel("Redshift")
                ax.text(0.5, 0.92, ion[p], horizontalalignment="center",
                        verticalalignment="center", transform=ax.transAxes,
                        fontsize=12)
        else:
            grid = gs(1,1)
            ax = fig.add_subplot(grid[:,:])
            ax.set_xlabel("Redshift")
            ax.set_ylabel("Flux [" + str(self._spec.y.unit) + "]")
            for p in range(n):
                spec = dc(self._spec)
                line = dc(self._line)
                spec.to_z([ion[p]])
                line.to_z([ion[p]])
                chunk = spec.t[np.logical_and(spec.x > zmin, spec.x < zmax)]
                memb = line.t[np.logical_and(line.x > zmin, line.x < zmax)]
                ax.set_xlim(zmin, zmax)
                ax.plot(chunk['X'], chunk['Y'])
                ax.scatter(memb['X'], memb['Y'])
            text = ', '.join(str(p) for p in ion)
            for comp in z:
                ax.axvline(x=comp, ymin=0.65, ymax=0.85, color='lightgray')

            ax.text(0.5, 0.92, text, horizontalalignment="center",
                    verticalalignment="center", transform=ax.transAxes,
                    fontsize=12)
           
                    
        grid.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])
        grid.update(hspace=0.0)
        if block is False:
            plt.ion()
            plt.draw()
        else:
            plt.show()

    def zmatch(self, ztol=1e-5):
        """ Match redshifts in a list, to define systems """

        # Flatten arrays
        # N.B. Y and DY don't need flattening, as they have only one entry
        # per row. They will be repeated to have the shape of the other arrays.
        z = np.ravel(self._redsh.x)
        y = np.repeat(self._redsh.y, self._redsh.x.shape[1])        
        zmin = np.ravel(self._redsh.xmin)
        zmax = np.ravel(self._redsh.xmax)
        dy = np.repeat(self._redsh.dy, self._redsh.x.shape[1])
        ion = np.ravel(self._redsh.ion)

        # Sort arrays
        argsort = np.argsort(z, kind='mergesort')
        z_sort = np.asarray(z[argsort])
        y_sort = np.asarray(y[argsort])      
        zmin_sort = np.asarray(zmin[argsort])
        zmax_sort = np.asarray(zmax[argsort])
        dy_sort = np.asarray(dy[argsort])
        ion_sort = np.asarray(ion[argsort])
        
        # Find coincidences
        coinc = np.isclose(z_sort[1:], z_sort[:-1], atol=ztol) 
        coupl = np.core.defchararray.not_equal(ion_sort[1:], ion_sort[:-1])
        coinc = np.logical_and(coinc, coupl)

        new = True
        z_coinc = []
        y_coinc = []
        zmin_coinc = []        
        zmax_coinc = []
        ion_coinc = []
        dy_coinc = []
        for i in range(len(z_sort) - 2):
            if ((new == True) and (coinc[i] == True)):
                new = False
                z_coinc_row = z_sort[i]
                y_coinc_row = (y_sort[i],)
                zmin_coinc_row = zmin_sort[i]
                zmax_coinc_row = zmax_sort[i]                
                dy_coinc_row = (dy_sort[i],)
                ion_coinc_row = (ion_sort[i],)
            if (coinc[i] == True):
                z_coinc_row = np.append(z_coinc_row, z_sort[i+1])
                y_coinc_row = y_coinc_row + (y_sort[i+1],)
                zmin_coinc_row = np.append(zmin_coinc_row, zmin_sort[i+1])
                zmax_coinc_row = np.append(zmax_coinc_row, zmax_sort[i+1])
                dy_coinc_row = dy_coinc_row + (dy_sort[i+1],)
                ion_coinc_row = ion_coinc_row + (ion_sort[i+1],)
            if (coinc[i+1] == False):
                if (new == False):
                    z_coinc.append(np.mean(z_coinc_row))
                    y_coinc.append(y_coinc_row)
                    zmin_coinc.append(np.mean(zmin_coinc_row))
                    zmax_coinc.append(np.mean(zmax_coinc_row))
                    dy_coinc.append(dy_coinc_row)
                    ion_coinc.append(ion_coinc_row)
                new = True
        syst = Syst(self.line, self.spec, x=z_coinc, y=y_coinc,
                    xmin=zmin_coinc, xmax=zmax_coinc, dy=dy_coinc,
                    ion=ion_coinc, yunit=y.unit)
        self.__dict__.update(syst.__dict__)

