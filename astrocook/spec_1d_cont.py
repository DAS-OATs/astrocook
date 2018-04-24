from . import Spec1D
from .utils import many_gauss, many_voigt, savitzky_golay
from astropy.constants import c
from astropy.io import fits as fits
from astropy.table import Column, Table
import copy
import numpy as np
import signal
from scipy.optimize import curve_fit
#TODO: uncomment this import statsmodels.api as sm
#import matplotlib.pyplot as plt

# From http://stackoverflow.com/questions/25027122/break-the-function-after-certain-time
class TimeoutException(Exception):   # Custom exception class
    pass

def handler(signum, frame):   # Custom signal handler
    raise Exception("stop")
    

class Spec1DCont(Spec1D):
    """Class for continuum fitting in spectra 
    
    A spectrum with fitted continuum is a generic spectrum with the following 
    additional columns: 
        -# @abs_fit: fitted profile of absorption lines;
        -# @em_fit: fitted profile of emission lines;      
        -# @abs_rem: flux density after removal of the absorption lines;
        -# @em_rem: flux density after removal of the emission lines;
        -# @cont: fiducial fitting of continuum flux density.
    """

    def __init__(self, spec=None,
                 abs_fit=None,
                 em_fit=None,
                 abs_rem=None,
                 em_rem=None,
                 cont=None,
                 meta=None,
                 dtype=float):

        # Column definition
        if (abs_fit is None):
            abs_fit = np.repeat(float('nan'), len(spec.x))
        col_abs_fit = Column(np.asarray(copy.deepcopy(abs_fit), dtype=dtype), 
                                        name='ABS_FIT')

        if (em_fit is None):
            em_fit = np.repeat(float('nan'), len(spec.x))
        col_em_fit = Column(np.asarray(copy.deepcopy(em_fit), dtype=dtype), 
                                        name='EM_FIT')

        if (abs_rem is None):
            abs_rem = np.repeat(float('nan'), len(spec.x))
        col_abs_rem = Column(np.asarray(copy.deepcopy(abs_rem), dtype=dtype), 
                                        name='ABS_REM')

        if (em_rem is None):
            em_rem = np.repeat(float('nan'), len(spec.x))
        col_em_rem = Column(np.asarray(copy.deepcopy(em_rem), dtype=dtype), 
                                        name='EM_REM')

        if (cont is None):
            cont = np.repeat(float('nan'), len(spec.x))
        col_cont = Column(np.asarray(copy.deepcopy(cont), dtype=dtype), 
                                        name='CONT')

        # Auxiliary data and meta
        if (meta is None):
            meta = {}

        # Table expansion
        self._t = copy.deepcopy(spec.t)
        self._t.add_columns([col_abs_fit, col_em_fit, col_abs_rem, col_em_rem, 
                             col_cont])
        self._t['ABS_FIT'].unit = self._t['Y'].unit
        self._t['EM_FIT'].unit = self._t['Y'].unit
        self._t['ABS_REM'].unit = self._t['Y'].unit
        self._t['EM_REM'].unit = self._t['Y'].unit
        self._t['CONT'].unit = self._t['Y'].unit
        self._use_good = False
        
    @property
    def abs_fit(self):
        return self._getWithMask('ABS_FIT')

    @abs_fit.setter
    def abs_fit(self, value):
        if self._use_good:
            self._t['ABS_FIT'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['ABS_FIT'] = np.asarray(value, dtype='float')
        self._t['ABS_FIT'].unit = self._t['ABS_FIT'].unit

    @property
    def em_fit(self):
        return self._getWithMask('EM_FIT')

    @em_fit.setter
    def em_fit(self, value):
        if self._use_good:
            self._t['EM_FIT'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['EM_FIT'] = np.asarray(value, dtype='float')
        self._t['EM_FIT'].unit = self._t['EM_FIT'].unit

    @property
    def abs_rem(self):
        return self._getWithMask('ABS_REM')

    @abs_rem.setter
    def abs_rem(self, value):
        if self._use_good:
            self._t['ABS_REM'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['ABS_REM'] = np.asarray(value, dtype='float')
        self._t['ABS_REM'].unit = self._t['ABS_REM'].unit
        
    @property
    def em_rem(self):
        return self._getWithMask('EM_REM')

    @em_rem.setter
    def em_rem(self, value):
        if self._use_good:
            self._t['EM_REM'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['EM_REM'] = np.asarray(value, dtype='float')
        self._t['EM_REM'].unit = self._t['EM_REM'].unit

    @property
    def cont(self):
        return self._getWithMask('CONT')

    @cont.setter
    def cont(self, value):
        if self._use_good:
            self._t['CONT'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['CONT'] = np.asarray(value, dtype='float')
        self._t['CONT'].unit = self._t['CONT'].unit

    def fit_lines_2(self, list, zem='0.0', col='y', timer=5):
        """Fit lines in a spectrum with Voigt profiles"""
    
        self_col = getattr(self, col)        
        y_fit = np.zeros(len(self_col)) * self_col.unit
        y_fit.fill(np.nan)
        y_rem = copy.deepcopy(self_col)
        
        rchisq = 10
        rchisq_thres = 5
        i = 0
        i_thres = 30
#        while rchisq > rchisq_thres:
        lyaf_where = self.t['X'] < (1 + zem) * 121.567 * self.x.unit 
        lyaf_spec = self.from_table(self.t[lyaf_where])
        fail = 0
        while i < i_thres:
        
            # Find the position of the maximum residual
            #res_idx = np.argmax(np.abs(1 - lyaf_spec.y / lyaf_spec.cont))
            res_idx = np.argmax(1 - lyaf_spec.y / lyaf_spec.cont)
            res_x = lyaf_spec.x[res_idx]
            print(res_x)
            # Find the line group where the maximum residual is located
            list_idx = (np.abs(list.x - res_x)).argmin()
            #list_where = list.t['GROUP'] == list.t['GROUP'][list_idx]
            list_where = list.t['X'] == list.t['X'][list_idx]
            group_list = list.from_table(list.t[list_where])
            group_xmin = np.amin(group_list.xmin)
            group_xmax = np.amax(group_list.xmax)
            spec_where = np.logical_and(self.x >= group_xmin, 
                                         self.x <= group_xmax)
            spec_sel = lyaf_spec.from_table(lyaf_spec.t[spec_where])
            spec_sel_col = getattr(spec_sel, col)
            #print(res_x, group_xmin, group_xmax)
            
            # Define the guess Voigt parameters for the residual
            mean = res_x.value
            log_N = 14.0
            b = 20.0
            p0_voigt = np.stack([mean, log_N, b])
            p0_voigt = np.ndarray.flatten(p0_voigt, order='F')
            #p0_voigt = np.tan(p0_voigt / 3600 - 0.0)
            #print(p0_voigt)
            
            try:

            
                # Fit a Voigt profile to the residual
                coeff, var_matrix = curve_fit(many_voigt, spec_sel.x.value, 
                                              spec_sel.y/spec_sel.cont, 
                                              p0=p0_voigt)
                #coeff = p0_voigt
                #print(coeff)
                #print(*coeff)
                prof = many_voigt(spec_sel.x.value, *coeff)
                #print(spec_sel.x.value)
                #print(prof)
                y_fit[spec_where] = spec_sel.cont.value * prof * self_col.unit                              
                y_rem[spec_where] = spec_sel_col + spec_sel.cont.value \
                                        * (1 - prof) * self_col.unit
                lyaf_spec.t['Y'][spec_where] = y_rem[spec_where]
                rchisq = np.sum(
                    (spec_sel.y - y_fit[spec_where]) ** 2 \
                    / spec_sel.dy ** 2) \
                    / (len(spec_sel.x) - len(group_list.x) * 3)

            except (RuntimeError, Exception):
                    
                # Use the guess local continuum
                print("!")
                y_fit[spec_where] = spec_sel_col                           
                y_rem[spec_where] = spec_sel.cont.value * self_col.unit                           
                lyaf_spec.t['Y'][spec_where] = y_rem[spec_where]
                fail = fail + 1
            finally:
                signal.alarm(0)


            i += 1

        self._t['ABS_FIT'] = y_fit
        self._t['ABS_REM'] = y_rem
    

    def fit_lines(self, list, col='y', timer=5):
        """Fit lines in a spectrum"""
        
        self_col = getattr(self, col)
        
        y_fit = np.zeros(len(self_col)) * self_col.unit
        y_fit.fill(np.nan)
        y_rem = copy.deepcopy(self_col)

        if list is not None:

            groups = np.unique(list.t['GROUP'])
            #groups = np.array([2])

            # Not working with Python 2.7
            #end=""
            #end="\n   "
            #print("  %i line groups (v = Voigt fit, X = failure): " \
            #      % (len(np.unique(list.t['GROUP']))), end=end, flush=True)
            fail = 0
            for group in groups:

                # Extract the lines of the group
                where_list = list.t['GROUP'] == group
                group_list = list.from_table(list.t[where_list])

                # Extract the spectral region of the group
                xmin = np.amin(group_list.xmin)
                xmax = np.amax(group_list.xmax)
                where_spec = np.logical_and(
                    self.x.value >= xmin.value, 
                    self.x.value <= xmax.value)
                group_spec = self.from_table(self.t[where_spec])
                group_spec_col = getattr(group_spec, col)

                fact = 3600

                # Compute the guess means of the profiles
                mean = group_list.x.value
                mean_min = group_list.x.value * (1 - 20000 / c.value)
                mean_max = group_list.x.value * (1 + 20000 / c.value)
                mean_norm = [500 * fact] * len(group_list.y)

                # Compute the guess local continuum by interpolating across the 
                # contiguous maxima
                x_i = group_spec.x[0].value
                x_f = group_spec.x[len(group_spec.x) - 1].value
                y_i = group_spec_col[0].value
                y_f = group_spec_col[len(group_spec_col) - 1].value
                
                cont = group_spec.cont.value
                if np.any(np.isnan(cont)):
                    cont = y_i + (group_spec.x.value - x_i) * (y_f - y_i) / (x_f - x_i)
            
                group_list_cont = np.interp(group_list.x.value,
                                            group_spec.x.value, cont)
                
                # Compute the guess column densities 
                #log_N = [13.0] * len(group_list.y)
                #print(group_list.guess_logN(group_list_cont))
                log_N = group_list.guess_logN(group_list_cont)
                log_N_min = [10.0] * len(group_list.y)
                log_N_max = [20.0] * len(group_list.y)
                log_N_norm = [14.0 * fact] * len(group_list.y)

                # Compute the guess broadening parameters
                b = [20.0] * len(group_list.y)
                b_min = [2.0] * len(group_list.y)
                b_max = [100.0] * len(group_list.y)
                b_norm = [20.0 * fact] * len(group_list.y)

                # Sort the parameters by descending column density
                """
                argsort = np.fliplr([np.argsort(log_N)])[0]
                mean = mean[argsort]
                mean_min = mean_min[argsort]
                mean_max = mean_max[argsort]
                log_N = log_N[argsort]
                """
                                
                p0_voigt = np.stack([mean, log_N, b])
                p0_voigt = np.ndarray.flatten(p0_voigt, order='F')
                min_voigt = np.stack([mean_min, log_N_min, b_min])
                min_voigt = np.ndarray.flatten(min_voigt, order='F')
                max_voigt = np.stack([mean_max, log_N_max, b_max])
                max_voigt = np.ndarray.flatten(max_voigt, order='F')
                norm_voigt = np.stack([mean_norm, log_N_norm, b_norm])
                norm_voigt = np.ndarray.flatten(norm_voigt, order='F')

                # Transform Voigt parameters
                p0_voigt = np.tan(p0_voigt / norm_voigt - 0.0)
                min_voigt = np.tan(min_voigt / norm_voigt - 0.0) 
                max_voigt = np.tan(max_voigt / norm_voigt - 0.0) 

                bounds_voigt = (min_voigt, max_voigt)

                norm = group_spec_col.value / cont
                dnorm = group_spec.dy.value / cont
                
                signal.signal(signal.SIGALRM, handler)
                signal.setitimer(signal.ITIMER_REAL, timer)


                norm_add = np.array([500 * fact, 14 * fact, 20.0 * fact])


                try:

                    coeff = []
                    
                    rchisq = 10
                    i = 0
                    
                    while rchisq > 5:
                    #for i in range(3, len(p0_voigt) + 1, 3):
                    
                        i = i + 1

                        #p0_ran = range(i - 3, i)
                        #p0_sel = np.append(coeff, p0_voigt[p0_ran])
                        #print(i, p0_voigt, p0_sel)
                        #bounds_ran = range(0, i)
                        #bounds_voigt = (min_voigt, max_voigt)
                        #bounds_sel = (bounds_voigt[0][bounds_ran], 
                        #              bounds_voigt[1][bounds_ran])

                        # Fit a sum of Voigt profiles
                        #print(i, p0_voigt)
                        coeff, var_matrix = curve_fit(
                                                many_voigt, group_spec.x.value, 
                                                norm, 
                                                p0=p0_voigt, sigma=dnorm, 
                                                #bounds=bounds_sel,
                                                bounds=bounds_voigt,
                                                method='trf')
                        #print(i, *coeff)
                        prof = many_voigt(group_spec.x.value, *coeff)
                        y_fit[where_spec] = cont * prof * self_col.unit                              
                        y_rem[where_spec] = group_spec_col + cont \
                                                * (1 - prof) * self_col.unit
                        rchisq = np.sum(
                            (group_spec.y - y_fit[where_spec]) ** 2 \
                            / group_spec.dy ** 2) \
                            / (len(group_spec.x) - len(group_list.x) * 3)
                        #print(i, rchisq)
                            
                        """
                        max_res = np.argmax(group_spec_col + cont \
                                            * (1 - prof) * self_col.unit)
                        p0_add = np.array([group_spec.x[max_res].value, 12, 20.0])
                        #min_add = np.stack([group_spec.x[max_res].value * (1 - 20000 / c.value), 10.0, 20.0])
                        #max_add = np.stack([group_spec.x[max_res].value * (1 + 20000 / c.value), 10.0, 20.0])

                        p0_add = np.tan(p0_add / norm_add)
                        #min_add = np.tan(min_add / norm_add - 0.0) 
                        #max_add = np.tan(max_add / norm_add - 0.0) 
                        p0_voigt = np.append(p0_voigt, p0_add)
                        #min_voigt = np.append(min_voigt, min_add)
                        #max_voigt = np.append(max_voigt, max_add)
                        """
                    # Not working with Python 2.7
                    #if end == "":
                    #    print("v ", end=end, flush=True)
                    #else:
                    #    print("[%f-%f]: v (%f) " 
                    #          % (xmin.value, xmax.value, rchisq), 
                    #          end=end, flush=True)

                #"""
                except (RuntimeError, Exception):
                    
                    # Use the guess local continuum
                    y_fit[where_spec] = group_spec_col                           
                    y_rem[where_spec] = cont * self_col.unit                           
                    fail = fail + 1
                    # Not working with Python 2.7
                    #if end == "":
                    #    print("X ", end=end, flush=True)
                    #else:
                    #    print("[%f-%f]:  X" % (xmin.value, xmax.value), 
                    #          end=end, flush=True)
                finally:
                    signal.alarm(0)
                #"""

            if end == "":
                print(" ")
            print("  Failure rate: %.2f%%" % (fail / len(groups) * 100))
            
        
        # Updating the properties em_fit, em_rem, etc. instead of the columns
        # has the effect of erasing the unit from the columns
        self._t['ABS_FIT'] = y_fit
        self._t['ABS_REM'] = y_rem
        
    def from_table(self, table, meta = {}):
        """Read a spectrum from a (spectrum-like) table"""
        
        xmin = table['XMIN']
        xmax = table['XMAX']
        x = table['X']
        y = table['Y']            
        dy = table['DY']
        group = table['GROUP']
        resol = table['RESOL']
        abs_fit = table['ABS_FIT']
        em_fit = table['EM_FIT']
        abs_rem = table['ABS_REM']
        em_rem = table['EM_REM']
        cont = table['CONT']
        dx = 0.5 * (xmax - xmin)
        
        c1 = np.argwhere(y > 0)
        c2 = np.argwhere(dy > 0)
        igood = np.intersect1d(c1, c2)

        good = np.repeat(-1, len(x))
        good[igood] = 1

        gen = Spec1D(x, y, dy=dy, xmin=xmin, xmax=xmax, xunit=x.unit, 
                     yunit=y.unit, group=good, resol=resol, meta=meta)
        cont = Spec1DCont(gen, abs_fit=abs_fit, em_fit=em_fit, abs_rem=abs_rem, 
                          em_rem=em_rem, cont=cont)
        return cont

    def lowess(self, x, y, frac):
        """Fit the continuum with a LOWESS algorithm"""

        self.cont = sm.nonparametric.lowess(y, x, frac=frac)[:, 1] * y.unit        

        if mode is 'em':
            self.em_rem = y_rem
        else:
            self.abs_rem = y_rem
        self.y_rem = y_rem        

    def rem_lines(self, list, mode='abs', col='y', size=101, order=3, timer=5):
        """Remove lines from a spectrum with Savitzky-Golay profiles"""
        
        self_col = getattr(self, col)
        y_fit = np.zeros(len(self_col)) * self_col.unit
        y_fit.fill(np.nan)
        y_rem = copy.deepcopy(self_col)

        if list is not None:

            groups = np.unique(list.t['GROUP'])
            #groups = [0]

            # Not working with Python 2.7
            #end=""
            #end="\n   "
            #print("  %i line groups (sg = Savitzky-Golay fit, " \
            #      "X = failure): " \
            #      % (len(np.unique(list.t['GROUP']))), end=end, flush=True)

            # Loop over groups
            fail = 0
            for group in groups:

                # Extract the lines of the group
                where_list = list.t['GROUP'] == group
                group_list = list.from_table(list.t[where_list])

                # Extract the spectral region of the group
                xmin = np.amin(group_list.xmin)
                xmax = np.amax(group_list.xmax)
                where_spec = np.logical_and(
                    self.x.value >= xmin.value, 
                    self.x.value <= xmax.value)
                group_spec = self.from_table(self.t[where_spec])
                group_spec_col = getattr(group_spec, col)

                # Compute the guess local continuum by interpolating across the 
                # contiguous maxima
                x_i = group_spec.x[0].value
                x_f = group_spec.x[len(group_spec.x) - 1].value
                y_i = group_spec_col[0].value
                y_f = group_spec_col[len(group_spec_col) - 1].value
                cont = group_spec.cont.value
                if np.any(np.isnan(cont)):
                    cont = y_i + (group_spec.x.value - x_i) \
                           * (y_f - y_i) / (x_f - x_i)
                
                group_list_cont = np.interp(group_list.x.value,
                                            group_spec.x.value, cont)
                # Set timer
                signal.signal(signal.SIGALRM, handler)
                signal.setitimer(signal.ITIMER_REAL, timer)

                try:

                    # Fit a Savitzky-Golay profile
                    group_spec.sg(cont * self_col.unit - group_spec_col, 
                                  size, order)
                    prof = group_spec.cont.value
                    y_fit[where_spec] = cont * self_col.unit \
                                        - prof * self_col.unit                              
                    y_rem[where_spec] = group_spec_col \
                                        + prof * self_col.unit  
                    # Not working with Python 2.7
                    #if end == "":
                    #    print("sg ", end=end, flush=True)
                    #else:
                    #    print("[%f-%f]: sg " % (xmin.value, xmax.value), 
                    #          end=end, flush=True)

                except:
                    
                    # Use the guess local continuum
                    y_fit[where_spec] = group_spec_col                           
                    y_rem[where_spec] = cont * self_col.unit                           

                    fail = fail + 1
                    # Not working with Python 2.7
                    #if end == "":
                    #    print("X ", end=end, flush=True)
                    #else:
                    #    print("[%f-%f]:  X" % (xmin.value, xmax.value), 
                    #          end=end, flush=True)
                    
                finally:
                    signal.alarm(0)
            if end == "":
                print(" ")

            if end == "":
                print(" ")
            print("  Failure rate: %.2f%%" % (fail / len(groups) * 100))

        # Updating the properties em_fit, em_rem, etc. instead of the columns
        # has the effect of erasing the unit from the columns
        if mode is 'em':
            self._t['EM_FIT'] = y_fit
            self._t['EM_REM'] = y_rem
        else:
            self._t['ABS_FIT'] = y_fit
            self._t['ABS_REM'] = y_rem


    def save(self, filename):
        x_unit = "{0}".format(self.x.unit)
        y_unit = "{0}".format(self.y.unit)
        hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='XMIN', format='E', unit=x_unit, array=self.xmin),
             fits.Column(name='XMAX', format='E', unit=x_unit, array=self.xmax),
             fits.Column(name='X', format='E', unit=x_unit, array=self.x),
             fits.Column(name='Y', format='E', unit=y_unit, array=self.y),
             fits.Column(name='DY', format='E', unit=y_unit, array=self.dy),
             fits.Column(name='GROUP', format='I', array=self.group),
             fits.Column(name='RESOL', format='E', array=self.resol),
             fits.Column(name='ABS_FIT', format='E', unit=y_unit, array=self.abs_fit),
             fits.Column(name='EM_FIT', format='E', unit=y_unit, array=self.em_fit),
             fits.Column(name='ABS_REM', format='E', unit=y_unit, array=self.abs_rem),
             fits.Column(name='EM_REM', format='E', unit=y_unit, array=self.em_rem),
             fits.Column(name='CONT', format='E', unit=y_unit, array=self.cont)])
        hdu.writeto(filename, overwrite=True)

    """
    def save(self, filename):
#        print(self._t['XMIN'].quantity.value)
        hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='XMIN', format='E', unit=self.xmin.unit, array=self.xmin),
             fits.Column(name='XMAX', format='E', unit=self.xmax.unit, array=self.xmax),
             fits.Column(name='X', format='E', unit=self.x.unit, array=self.x),
             fits.Column(name='Y', format='E', unit=self.y.unit, array=self.y),
             fits.Column(name='DY', format='E', unit=self.dy.unit, array=self.dy),
             fits.Column(name='GROUP', format='I', array=self.group),
             fits.Column(name='RESOL', format='E', array=self.resol),
             fits.Column(name='ABS_FIT', format='E', unit=self.abs_fit.unit, array=self.abs_fit),
             fits.Column(name='EM_FIT', format='E', unit=self.em_fit.unit, array=self.em_fit),
             fits.Column(name='ABS_REM', format='E', unit=self.abs_rem.unit, array=self.abs_rem),
             fits.Column(name='EM_REM', format='E', unit=self.em_rem.unit, array=self.em_rem),
             fits.Column(name='CONT', format='E', unit=self.cont.unit, array=self.cont)])
            [fits.Column(name='XMIN', format='E', array=self._t['XMIN']),
             fits.Column(name='XMAX', format='E', array=self._t['XMAX']),
             fits.Column(name='X', format='E', array=self._t['X']),
             fits.Column(name='Y', format='E', array=self._t['Y']),
             fits.Column(name='DY', format='E', array=self._t['DY']),
             fits.Column(name='GROUP', format='I', array=self._t['GROUP']),
             fits.Column(name='RESOL', format='E', array=self._t['RESOL']),
             fits.Column(name='ABS_FIT', format='E', array=self._t['ABS_FIT']),
             fits.Column(name='EM_FIT', format='E', array=self._t['EM_FIT']),
             fits.Column(name='ABS_REM', format='E', array=self._t['ABS_REM']),
             fits.Column(name='EM_REM', format='E', array=self._t['EM_REM']),
             fits.Column(name='CONT', format='E', array=self._t['CONT'])])
        hdu.writeto(filename, overwrite=True)
    """

    def sg(self, y, width, deg):
        """Fit the continuum with a Savitzky-Golay filter"""
    
        self._t['CONT'] = savitzky_golay(y, width, deg) * y.unit


    def findStretchable(self,
                        spec,
                        stiff=None,
                        pnratio=None,
                        minStep=None,
                        tol=None):
        """Find absorbed continuum.

        """

        if (stiff is None):
            stiff = 0.9

        if (pnratio is None):
            pnratio = 20.

        if (minStep is None):
            minStep = 5e-4

        if (tol is None):
            tol = 0.01

        if (abs(stiff) > 1):
            raise ValueError('Stiff must be a number between 0 and 1')

        if (abs(stiff) <= 0):
            raise ValueError('pnRatio must be a number > 0')

        cStr = stiff
        cUp  = 1.-stiff
        cDn  = cUp / float(pnratio)
        print("Settings: ", cStr, cUp, cDn)

        offsetX = spec.x.min().value
        dataX   = spec.x.value - offsetX
        scaleX  = dataX.max()
        dataX  /= scaleX

        offsetY = spec.y.min().value
        dataY   = spec.y.value - offsetY
        scaleY  = dataY.max()
        dataY  /= scaleY

        deltaXsq = ( dataX[1:] - dataX[0:-1] )**2.
        deltaXsq = np.insert(deltaXsq, 0, deltaXsq[0])

        around = 1
        npar = 2
        modelY = None
        lastTotScore = None
        while True:
            npar += (npar-1)

            if (npar > len(dataX)):
                print('Algorithm did not converged')
                break

            index = np.int_(np.round(np.linspace(0, len(dataX)-1, npar)))
            if (modelY is not None):
                newModelX = dataX[index]
                modelY = np.interp(newModelX, modelX, modelY)
                modelX = copy.deepcopy(newModelX)
            else:
                modelX = dataX[index]
                modelY = np.interp(modelX, dataX, dataY)

            go = np.ones(npar)
            cur = -1
            while go.sum() > 0:
                cur = (cur+1) % npar
                if (go[cur] == 0):
                    cur += 1
                    continue

                i0 = (cur-around)
                if (i0 < 0):
                    i0 = 0

                i1 = (cur+around)
                if (i1>(npar-1)):
                    i1 = npar-1

                jj = (dataX >= modelX[i0])  &  (dataX <= modelX[i1])
                valAtStart = modelY[cur]
                step = 0.1
                while step > minStep:
                    for sign in [1, -1]:

                        lastScore = -1
                        while True:
                            if (lastScore >= 0):
                                modelY[cur] += sign * step

                            modelCmp = np.interp(dataX[jj], modelX[i0:i1+1], modelY[i0:i1+1])
                            deltaModelSq = (modelCmp[1:] - modelCmp[0:-1])**2.
                            deltaModelSq = np.insert(deltaModelSq, 0, 0)
                            aa = cStr * np.ones(len(deltaXsq[jj]))
                            aa[(dataX[jj] < 0.3)] *= 10

                            tmp = dataY[jj]-modelCmp
                            iUp = (tmp > 0)
                            iDn = (tmp < 0)
                            score = \
                                    (aa * np.sqrt(deltaModelSq + deltaXsq[jj])).sum() + \
                                    cUp  * tmp[iUp].sum()                             - \
                                    cDn  * tmp[iDn].sum()

                            if (lastScore != -1):
                                if (score >= lastScore):
                                    modelY[cur] -= sign * step #revert
                                    break

                            lastScore = score
                    step /= 20.

                diff = abs(modelY[cur] - valAtStart)
                if (diff < minStep):
                    go[cur] = 0
                else:
                    go[i0:i1+1] = 1

            modelCmp = np.interp(dataX, modelX, modelY)
            deltaModelSq = (modelCmp[1:] - modelCmp[0:-1])**2.
            deltaModelSq = np.insert(deltaModelSq, 0, 0)
            tmp = dataY-modelCmp
            iUp = (tmp > 0)
            iDn = (tmp < 0)
            score1 = np.sqrt(deltaModelSq + deltaXsq).sum()
            score2 = tmp[iUp].sum()
            score3 = tmp[iDn].sum()
            totScore = cStr * score1 + \
                       cUp  * score2 - \
                       cDn  * score3
            print('Scores: ', score1, score2, score3, totScore, len(modelX))
            if (lastTotScore is not None):
                if (abs(lastTotScore - totScore)/totScore < tol):
                    break

            lastTotScore = totScore

            #plt.plot(dataX, dataY, marker='.', linestyle='None', markersize=1)
            #plt.plot(dataX, modelCmp)
            #plt.show()

        i = np.argwhere(modelCmp > dataY)
        print('Length:   ', score1)
        print('PN ratio: ', len(i)/float(len(dataX) - len(i)), ' ~ ', pnratio)

        modelCmp = np.interp(dataX, modelX, modelY)
        modelCmp *= scaleY
        modelCmp += offsetY

        return modelCmp
