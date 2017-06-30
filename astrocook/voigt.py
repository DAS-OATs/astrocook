from astropy.constants import c, e, m_e
from astropy.table import Column, Table
import copy
from lmfit import Model, Parameters
from lmfit.models import ConstantModel, LinearModel, VoigtModel, ExpressionModel
import numpy as np
import random
from scipy.signal import fftconvolve
from scipy.special import wofz

lya_x = 121.567

    
def func_norm(x, y_norm):
    return y_norm

def func_voigt(x, z, logN, b, btur):
    wave = 121.567
    f = 0.416
    Gamma = 6.265e8
    #Why sqrt(2)? Apparently it works
    amplitude = np.sqrt(2) * np.pi * f * e.esu.value ** 2 \
                / (m_e.value * c.value) * 1e-18 * np.power(10, logN)
    center = wave * (1 + z)
    sigma = np.sqrt(b ** 2 + btur ** 2) * 1e-1 / (np.sqrt(2) * wave)
    gamma = Gamma * 1e-10 / (4 * np.pi)
    return np.exp(-amplitude * np.real(wofz(
        (x - center + 1j * gamma) /\
        (sigma * np.sqrt(2)))) /\
        (sigma * np.sqrt(2 * np.pi)))

def func_voigt2(x, z, N, b, btur, data=None, eps=None):
    wave = 121.567
    f = 0.416
    Gamma = 6.265e8
    #Why sqrt(2)? Apparently it works
    amplitude = np.sqrt(2) * np.pi * f * e.esu.value ** 2 \
                / (m_e.value * c.value) * 1e-18 * N
    center = wave * (1 + z)
    sigma = np.sqrt(b ** 2 + btur ** 2) * 1e-1 / (np.sqrt(2) * wave)
    gamma = Gamma * 1e-10 / (4 * np.pi)
    model = np.exp(-amplitude * np.real(wofz(
        (x - center + 1j * gamma) /\
        (sigma * np.sqrt(2)))) /\
        (sigma * np.sqrt(2 * np.pi)))
    #print(np.array(model))
    #print(np.array(data))
    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model - data)/eps

class Voigt():

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

    @t.setter
    def t(self, value):
        self._t = value

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
        """Quantities associated to spectrum channels"""
        return self._mask_col('Y')

    @y.setter
    def y(self, value):
        if self._use_good:
            self._t['Y'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['Y'] = np.asarray(value, dtype='float')
        self._t['Y'].unit = self._t['YMIN'].unit

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


    def line_add(self, xmin, xmax, x, y):
        z = x / 121.567 - 1
        if np.isclose(self._t['Z'], z, atol=1e-7).any() == False or 1 == 1:
            self.t.add_row([xmin, xmax , x, y, None, 1, 'Ly_a', z, 14.0, 20.0,
                            0.0])
            ret = True
        else:
            ret = False
        return ret
        
    def fit(self, line, iter=0, prev=None):
        self.prep(line)

        #cont_model, cont_param =
        self.model_cont(line, prev)        
        #trasm_model, trasm_param = self.model_trasm(line)
        trasm_model = self.model_trasm(line)

        ftol = 1e-7 # max(pow(10, -iter), 1e-7)
        out = trasm_model.fit(self._y_rect, self._param, x=self._x_ran,
                              fit_kws={
                                  #'ftol': ftol,
                                  'maxfev': 1000
                              },
                              weights=1/self._dy_rect)

        self._y_fit = out.best_fit
        self._y_resid = self._y_rect - self._y_fit
        self._y_norm = np.full(len(self._x_ran), out.params['y_norm'].value)
        self._y_cont = self._y_slope0 * out.params['y_norm'].value #self._y_norm
        self._param = out.params
        self._out_redchi = np.sum(
                ((self._y_rect - self._y_fit) / self._dy_rect) ** 2) / out.nfree

        self.fit_save(line)
        return out        

    def fit_auto(self, line):
        add = True


        out = self.fit(line)

        stop = self._out_redchi >= self._redchi
        i = 0
        while stop == False: 
            i += 1
            print("[Iteration %2i] %4i function calls, " \
                  "reduced chi-squared: %f, AIC: %f" \
                  % (i, out.nfev, self._out_redchi, out.aic))

            xmin = np.min(self.group(line)['XMIN'])
            xmax = np.max(self.group(line)['XMAX'])        

            par = np.stack([[np.mean(self._x_ran)], [0.2], [1]])

            prof = np.exp(-(self._x_ran-np.mean(self._x_ran))**2/(2.*0.01**2))
            prof = prof / np.sum(prof)
            self._y_conv = fftconvolve(self._y_resid, prof, mode='same')
            
            x = self._x_ran[np.argmin(self._y_conv)]
            y = np.interp(x, self._x_ran, self._y_ran)

            add = self.line_add(xmin, xmax, x, y)
            ret = out
            out = self.fit(line, iter=i, prev=out)

            stop = (add == False or self._out_redchi >= self._redchi)
            if stop == True:
                self._y_fit = ret.best_fit
                self._y_resid = self._y_rect - self._y_fit
                self._y_norm = np.full(len(self._x_ran),
                                       ret.params['y_norm'].value)
                self._y_cont = self._y_slope0 * ret.params['y_norm'].value 
                self._param = ret.params
                self._out_redchi = self._redchi
            self._redchi = self._out_redchi
            """
            print(out.redchi)
            print(np.sum(out.residual ** 2) / out.nfree)
            print(np.sum(((self._y_rect - self._y_fit) / self._dy_rect) ** 2) / out.nfree)
            print(np.array(out.residual[0:6]))
            print(np.array((self._y_rect - self._y_fit) / self._dy_rect)[0:6])
            print(np.array(self._dy_rect[0:6]))
            """
            #out = self.fit(line, iter=i, prev=out)
        return ret
        

    def fit_save(self, line):
        param_copy = [value for (key,value) in sorted(self._param.items())]
        g = 0
        for l in range(len(self.t)):
            if self.group_sel(line)[l] == True:
                self.t['Z'][l] = param_copy[g + 4]
                self.t['LOGN'][l] = np.log10(param_copy[g + 1]) 
                self.t['B'][l] = param_copy[g + 2] 
                self.t['BTUR'][l] = param_copy[g + 3]         
                g += 4
                
    def guess(self, y_norm):
        a_0 =  14.096;
        a_1 = - 4.6251;
        a_2 =  18.657;
        a_3 = -46.299;
        a_4 =  53.301;
        a_5 = -23.442;
        return a_0 + a_1 * (y_norm) + a_2 * pow(y_norm, 2) \
               + a_3 * pow(y_norm, 3) + a_4 * pow(y_norm, 4) \
               + a_5 * pow(y_norm, 5);
    
    def model_cont(self, line, prev=None):
        y_thres = 1e-3
        if prev == None:
            cont = self._y_extr[0] + (self._x_ran - self._x_extr[0]) \
                   * (self._y_extr[1] - self._y_extr[0]) \
                   / (self._x_extr[1] - self._x_extr[0])
            self._y_cont = cont
            self._y_norm = np.full(len(self._x_ran), np.mean(self._y_cont))

            self._y_cont0 = self._y_cont
            self._y_norm0 = self._y_norm
            self._y_slope0 = self._y_cont0 / self._y_norm0

        self._y_slope = self._y_cont / self._y_norm
        self._y_rect = self._y_ran / self._y_slope
        self._dy_rect = self._dy_ran / self._y_slope        

        #return model, param

    def model_trasm(self, line):
        wave = 121.567
        f = 0.416
        Gamma = 6.265e8

        group = self.group(line)
        norm_model = Model(func_norm)
        param = norm_model.make_params()
        y_norm = self._y_norm[0]
        y_norm0 = self._y_norm0[0]        
        y_span = 1 + np.max(self._dy_rect) / np.max(self._y_rect)
        y_min = np.max(self._y_rect) / y_span
        y_max = np.max(self._y_rect) * y_span        
        param['y_norm'].set(y_norm, #vary=False,
                            min=y_min, max=y_max)

        
        z_diff = 1e-4
        logN_min = 10
        logN_max = 22
        N_min = 1e10
        N_max = 1e22
        b_min = 1
        b_max = 100
        btur_min = 0
        btur_max = 100
        model = norm_model
        pref_list = []
        for g in range(len(group)):
            pref = 'z' + str(group['Z'][g]).replace('.', '') + '_'
            if pref in pref_list:
                pref += str(g) + '_'
            pref_list.append(pref)
            voigt_model = Model(func_voigt2, prefix=pref)#, data=self._y_rect)#, eps=self._dy_rect)
            param.update(voigt_model.make_params())            

            z = group['Z'][g]
            #if self._improve == False or group['LOGN'][g] == 14.0:
            if group['LOGN'][g] == 14.0:
                try:
                    """
                    if self._improve == True:
                        y_interp = np.interp(group['X'][g], self._x_ran,
                                             1 - self._y_resid / self._y_cont)
                    else:
                        y_interp = np.interp(group['X'][g], self._x_ran,
                                             self._y_ran / self._y_cont)
                    """
                    y_interp = np.interp(group['X'][g], self._x_ran,
                                         self._y_ran / self._y_cont)
                except:
                    y_interp = np.interp(group['X'][g], self._x_ran,
                                         self._y_ran / self._y_cont)
                logN = self.guess(y_interp)
            else:
                logN = group['LOGN'][g]
            N = np.power(10, logN)
            b = group['B'][g]
            btur = group['BTUR'][g]
            param[pref + 'z'].set(z, min=z-z_diff, max=z+z_diff)
            param[pref + 'N'].set(N, min=N_min, max=N_max)
            param[pref + 'b'].set(b, min=b_min, max=b_max)
            param[pref + 'btur'].set(btur, vary=False, min=btur_min,
                                     max=btur_max)
            model *= voigt_model

        self._y_trasm = model.eval(param, x=self._x_ran)
        #self._y_trasm = (self._y_rect + model.eval(param, x=self._x_ran)) \
        #                * self._dy_rect
        #self._y_norm = np.full(len(self._x_ran), param['y_norm'].value)
        self._param = param
        return model #, param
    
    """
    def param(self, line):
        # Fit parameters
        param = Parameters()
        if line == -1:
            list = range(len(self._x))
        else:
            list = range(line, line + 1)
        for l in list:
            pref = 'z' + str(self._z[l]).replace('.', '') + '_'
            param.add(pref + 'z', value=self._z[l])
            param.add(pref + 'logN', value=self._logN[l])
            param.add(pref + 'b', value=self._b[l])
            param.add(pref + 'btur', value=self._btur[l])
        return param
    """
    def prep(self, line=-1):
        ran = self.range(line)
        self._x_ran = ran['X']
        self._y_ran = ran['Y']
        self._dy_ran = ran['DY']
        self._x_extr = np.array([ran['X'][0], ran['X'][len(ran) - 1]]) 
        self._y_extr = np.array([ran['Y'][0], ran['Y'][len(ran) - 1]]) 
        
        
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

