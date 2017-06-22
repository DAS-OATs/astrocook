from astropy.constants import c, e, m_e
from astropy.table import Column, Table
import copy
from lmfit import Parameters
from lmfit.models import LinearModel, VoigtModel
import numpy as np
import random
from scipy.signal import fftconvolve

lya_x = 121.567

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
        """Quantities associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
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
        """
        if np.isclose(self._t['Z'], z, atol=1e-7).any() == True:
            z += 1e-7 * (random.randint(0, 1) * 2 - 1)
        """
        if np.isclose(self._t['Z'], z, atol=1e-7).any() == False:
            self.t.add_row([xmin, xmax , x, y, None, 1, 'Ly_a', z, 14.0, 20.0,
                            0.0])
            ret = True
        else:
            ret = False

        return ret
        
    def fit(self, line, prev=None):
        self.prep(line)
        
        cont_model, cont_param = self.model_cont(line)        
        line_model, line_param = self.model_line(line)
        model = line_model
        std_param = line_param
        if prev != None:
            for param in prev.params:
                std_param[param] = prev.params[param]
        
        out = model.fit(self._tau_ran, std_param, x=self._x_ran,
                        weights=1/self._dtau_ran)

        self._tau_fit = out.best_fit
        self._tau_resid = out.residual * self._dtau_ran
        self._y_fit = np.exp(-self._tau_fit) * self._y_cont
        self._y_resid = self._y_ran - self._y_fit

        """
        cont_out = cont_model.fit(np.exp(-self._tau_resid) * self._y_cont,
                                  cont_param, x=self._x_ran,
                                  weights=1/self._dy_ran) 
        
        print(out.fit_report())
        print(cont_out.fit_report())
        """
        return out

    def fit_auto(self, line):
        out = self.fit(line)
        add = True
        i = 0
        while out.redchi > 1 and add == True:
            i += 1
            print("Iteration %i, reduced chi-squared: %f" % (i, out.redchi))
            xmin = np.min(self.group(line)['XMIN'])
            xmax = np.max(self.group(line)['XMAX'])        

            par = np.stack([[np.mean(self._x_ran)], [0.2], [1]])

            prof = np.exp(-(self._x_ran-np.mean(self._x_ran))**2/(2.*0.01**2))
            prof = prof / np.sum(prof)
            self._y_conv = fftconvolve(self._y_resid, prof, mode='same')
            
            x = self._x_ran[np.argmin(self._y_conv)]
            y = np.interp(x, self._x_ran, self._y_ran)
            add = self.line_add(xmin, xmax, x, y)
            #print(self.group(line))
            out = self.fit(line, prev=out)

    
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

    def guess_line(self, z, y_norm):
        wave = 121.567
        f = 0.416
        Gamma = 6.265e8

        a_0 =  14.096;
        a_1 = - 4.6251;
        a_2 =  18.657;
        a_3 = -46.299;
        a_4 =  53.301;
        a_5 = -23.442;

        logN = a_0 + a_1 * (y_norm) + a_2 * pow(y_norm, 2) \
               + a_3 * pow(y_norm, 3) + a_4 * pow(y_norm, 4) \
               + a_5 * pow(y_norm, 5)
        b = 20.0
        btur = 0.0
        
        line_model = VoigtModel()
        std_param = line_model.make_params()
        """
        z = group['Z'][g]
        logN = group['logN'][g]
        b = group['B'][g]
        btur = group['BTUR'][g]
        """
        ampl = np.sqrt(2) * np.pi * f * e.esu.value ** 2 \
               / (m_e.value * c.value) * 1e-18 * np.power(10, logN) #1e-18
        center = wave * (1 + z)
        sigma = np.sqrt(b ** 2 + btur ** 2) * 1e-1 / (np.sqrt(2) * wave) #1e-1
        gamma = Gamma * 1e-10 / (4 * np.pi) #1e-10
        std_param['amplitude'].set(ampl, min=0)
        std_param['center'].set(center)
        std_param['sigma'].set(sigma)
        std_param['gamma'].set(gamma)

        self._tau_pre += line_model.eval(std_param, x=self._x_ran)
        self._y_pre = np.exp(-self._tau_pre) * self._y_cont

        out = line_model.fit(self._tau_ran, std_param, x=self._x_ran,
                             weights=self._dtau_ran)
        return out
        
    def model_cont(self, line):
        model = LinearModel(prefix='cont_')
        std_param = model.guess(self._y_extr, x=self._x_extr)
        intercept = std_param['cont_intercept'].value
        slope = std_param['cont_slope'].value
        std_param['cont_intercept'].set(intercept)
        std_param['cont_slope'].set(slope)#, vary=False)            

        self._y_cont = model.eval(std_param, x=self._x_ran)
        self._y_norm = self._y_ran / self._y_cont
        self._dy_norm = self._dy_ran / self._y_cont
        self._tau_ran = -np.log(np.absolute(self._y_norm))
        y_thres = 1e-3
        self._tau_ran[self._y_norm < y_thres] = \
            -self._y_norm[self._y_norm < y_thres] / y_thres + 1 -np.log(y_thres)
        self._dtau_ran = np.absolute(self._dy_norm / self._y_norm)
        return model, std_param

    
    def model_line(self, line):
        wave = 121.567
        f = 0.416
        Gamma = 6.265e8

        group = self.group(line)
        self._tau_pre = np.zeros(len(self._x_ran))
        for g in range(len(group)):

            pref = 'z' + str(group['Z'][g]).replace('.', '') + '_'
            line_model = VoigtModel(prefix=pref)
            if g == 0:
                std_param = line_model.make_params()
                #std_param = line_model.guess(self._y_ran, x=self._x_ran)
            else:
                std_param.update(line_model.make_params())
                #std_param.update(line_model.guess(self._y_ran, x=self._x_ran))
            y_norm = np.interp(group['X'][g], self._x_ran, self._y_norm)
            #"""
            z = group['Z'][g]
            logN = self.guess(y_norm)
            #print(logN)
            b = group['B'][g]
            btur = group['BTUR'][g]
            #Why sqrt(2)? Apparently it works
            ampl = np.sqrt(2) * np.pi * f * e.esu.value ** 2 \
                   / (m_e.value * c.value) * 1e-18 * np.power(10, logN)
            center = wave * (1 + z)
            sigma = np.sqrt(b ** 2 + btur ** 2) * 1e-1 / (np.sqrt(2) * wave)
            gamma = Gamma * 1e-10 / (4 * np.pi)

            """
            out = self.guess_line(group['Z'][g], y_norm)
            out.params.pretty_print()
            ampl = out.params['amplitude'].value
            center = out.params['center'].value
            sigma = out.params['sigma'].value
            gamma = out.params['gamma'].value
            """
            std_param[pref + 'amplitude'].set(ampl, min=0)
            std_param[pref + 'center'].set(center)#, min=center-0.01, max=center+0.01)
            std_param[pref + 'sigma'].set(sigma)
            std_param[pref + 'gamma'].set(gamma)

            if g == 0:
                model = line_model
            else:
                model += line_model

        #model2 = np.exp(-model) * self._y_cont
                    
        self._tau_line = model.eval(std_param, x=self._x_ran)
        self._y_line = np.exp(-self._tau_line) * self._y_cont

        return model, std_param

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

