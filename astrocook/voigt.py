from astropy.constants import c, e, m_e
from astropy.table import Column, Table
from collections import OrderedDict
import copy
from lmfit import CompositeModel, Model, Parameters
from lmfit.models import ConstantModel, LinearModel, VoigtModel, ExpressionModel
from lmfit.lineshapes import gaussian
import matplotlib.pyplot as plt
import os
import random
from scipy.interpolate import RectBivariateSpline, bisplrep, bisplev, interp2d
from scipy.signal import fftconvolve
from scipy.special import wofz
from scipy.stats import linregress
import sys

#wave = 121.567
#f = 0.416
#Gamma = 6.265e8
redchi_thr = 0.9

ion_dict = {'Ly_a': [121.567, 0.416, 6.265e8],
            'CIV_1548': [154.8204, 0.1899, 2.643e8],
            'CIV_1550': [155.0781, 0.09475, 2.628e8]} 

resol = 3.0e4
def pause():
    programPause = input("Press the <ENTER> key to continue...")
    
def convolve(arr, kernel):
    """ Convolve an array with a kernel """
    npts = min(len(arr), len(kernel))
    pad  = np.ones(npts)
    tmp  = np.concatenate((pad*arr[0], arr, pad*arr[-1]))
    out  = np.convolve(tmp, kernel, mode='valid')
    noff = int((len(out) - npts)/2)
    return out[noff:noff+npts] / np.sum(kernel)

def func_fadd(a, u):
    """ Compute the real part of the Faddeeva function """
    return np.real(wofz(u + 1j * a))

def func_lin(x, y_norm, y_slope):
    #print(np.array(y_norm + (x - np.mean(x)) * y_slope)[0])
    return y_norm + (x - np.mean(x)) * y_slope

def func_norm(x, y_norm):
    return y_norm

def func_voigt(x, z, N, b, btur, ion='Ly_a', tab=None):
    """ Compute the Voigt function """
    wave = ion_dict[ion][0] * 1e-9
    f = ion_dict[ion][1]
    Gamma = ion_dict[ion][2]
    x_si = x * 1e-9
    N_si = N * 1e4
    b_si = np.sqrt(b**2 + btur**2) * 1e3
    tau0 = N_si * np.sqrt(np.pi) * f * e.esu.value**2 / (m_e.value * c.value) \
           * 1e-9 \
           * wave / b_si
    a = 0.25 * Gamma * wave / (np.pi * b_si)
    u = c.value / b_si * (x_si / (wave * (1 + z)) - 1)
    if tab == None:
        model = np.exp(-tau0 * func_fadd(a, u))
    else:
        model = np.exp(-tau0 * tab(a, u).flatten())
    return model
        
class Voigt():

    def __init__(self, spec, lines,
                 chosen=-1,
                 syst=None,
                 logN=None,
                 b=None,
                 btur=None):

        self._xmin = copy.deepcopy(lines.xmin)
        self._xmax = copy.deepcopy(lines.xmax)
        self._x = copy.deepcopy(lines.x)
        self._y = copy.deepcopy(lines.y)
        self._t = copy.deepcopy(lines.t)
        self._spec = copy.deepcopy(spec.t)
        self._lines = copy.deepcopy(lines)
        self._chosen = chosen

        # Arrays with system IDs and redshifts
        if (syst == None):
            lya_z = self._x / ion_dict['Ly_a'][0] - 1.
            self._syst = np.full((len(self._x), 1), 'Ly_a', dtype=object)
            self._z = np.full(len(self._x), lya_z, dtype=float)
        elif (len(syst) == len(self._x)):
            self._syst = np.array(syst, dtype=object)
            wave = []
            for syst in self._syst:
                wave.append(ion_dict[syst[0]][0])
            self._z = np.array(self._x / wave - 1., dtype=float)
        else:
            print("System IDs not recognized!")
        col_syst = Column(self._syst, name='SYST')
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
        
        col_add = Table(data=(col_syst, col_z, col_logN, col_b, col_btur),
                        masked=True)
        
        # To handle table with preexisting Voigt columns
        try:
            self._t.add_columns(col_add.columns.values())
        except:
            pass
            
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
        if self._chosen == -1:
            if self._use_good:
                ret = self._t[self._igood]
            else:
                ret = self._t
        else:
            line = self._lines.from_table(self._t[self._chosen: self._chosen+1])
            if self._use_good:
                ret = line._t[self._igood]
            else:
                ret = line._t
        return ret
                
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
        """Quantities associated to spectrum channels"""
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

    def fit(self, line, iter=0, maxfev=500, ax=None):
        """ Fit a composite continuum + Voigt model to a system """

        # Fit the model
        out = self._trasm_model.fit(self._y_ran, self._param, x=self._x_ran,
                                    fit_kws={'maxfev': maxfev},
                                    weights=1/self._dy_ran)
        #self._param.pretty_print()
        self._fit = out
        #out.params.pretty_print()
        
        # Save the results

        # Plot the results
            
        #return out        

    def fit_auto(self, line, cont=None, redchi=float('inf'), aic=float('inf'),
                 ax=None):
        """ Incrementally fit a system by adding components for residuals """

        i = 0
        self._out_redchi = redchi
        self._out_aic = aic
        stop = False #self._out_redchi < 1
        prev = None
        add = True
        old_t = copy.deepcopy(self._t)
        while stop == False:
            i += 1
            old_redchi = self._out_redchi
            old_aic = self._out_aic
            self.prep(line, prev)
            self.model_cont(line, prev, cont)
            self._trasm_model = self.model_trasm(line)
            # Plot fit
            if ax != None:
                comp = self.group(line)
                for x in comp['X']:
                    ax.axvline(x=x, ymin=0.75, ymax=0.95, color='lightgray')
                ax.plot(self._x_ran, self._y_ran / self._y_cont, c='b')
                ax.plot(self._x_ran, self._dy_ran / self._y_cont, c='b',
                        linestyle=':')
                ax.plot(self._x_ran, -self._dy_ran / self._y_cont, c='b',
                        linestyle=':')
                ax.plot(self._x_ran, self._y_trasm / self._y_cont, c='r',
                        linestyle=':')
                ax.plot(self._x_ran, self._y_cont0 / self._y_cont, c='y',
                        linestyle=':')
                try:
                    ax.plot(self._x_ran, self._y_fit / self._y_cont, c='g')
                    ax.plot(self._x_ran, self._y_cont / self._y_cont, c='y')
                    ax.plot(self._x_ran, self._y_resid / self._y_cont, c='g',
                            linestyle='--')
                except:
                    pass
                plt.draw()
                plt.pause(0.1)

            self.fit(line, ax=ax)
            

            stop = (add == False or (self._fit.redchi < redchi_thr and i > 1) \
                    or (self._fit.aic > old_aic))
            
            if stop == False:

                prev = self._fit
                norm = self._y_cont

                print("[Iteration %2i] %3i components, %4i function calls, " \
                  "reduced chi-squared: %f, AIC: %f" \
                      % (i, len(self.group(line)), self._fit.nfev,
                         self._fit.redchi, self._fit.aic))

                # Save fit
                if self._fit.redchi < old_redchi: #or 1 == 1:
                    self.fit_save(line)
                    old_redchi = self._out_redchi
                    old_aic = self._out_aic

                # Add line where the minimum residual is located
                x_add = self._x_ran[np.argmin(self._y_resid / self._dy_ran)]
                y_add = np.interp(x_add, self._x_ran, self._y_ran)
                y_resid_add = np.interp(x_add, self._x_ran,
                                        self._y_resid / self._y_cont)
                old_t = copy.deepcopy(self._t)
                syst_add = self._t['SYST'][line]  # Check this, because it
                                                  # happens that the largest
                                                  # residual is the other member
                """
                xmin = np.min(self.group(line)['XMIN'])
                xmax = np.max(self.group(line)['XMAX'])        
                """
                x_where= np.where(self.group(line)['SYST'] == syst_add)[0]
                xmin = np.min(self.group(line)['XMIN'][x_where])
                xmax = np.max(self.group(line)['XMAX'][x_where])        
                add = self.line_add(xmin, xmax, x_add, y_add, syst=syst_add)

                # Merge close lines (not so stable)
                tol = 1 / 3 * np.mean(self._x_ran) / resol
                #self.line_merge(tol=tol)

                if ax != None:# and add == True:
                    ax.cla()
                    ax.scatter(x_add, y_resid_add)

                # Delete weak lines
                delete = self.line_delete(self._t['LOGN'] < 10)


            else:
                self._fit = prev
                self._t = copy.deepcopy(old_t)
                
        if stop == True and i == 0:
            self.prep(line, prev)
            self.model_cont(line, prev, cont)        
            trasm_model = self.model_trasm(line)

            
    def fit_save(self, line):

        self._y_fit = self._fit.best_fit
        self._y_resid = self._y_ran - self._y_fit
        self._y_cont = np.full(len(self._x_ran),
                               self._fit.params['y_norm'].value) \
                               + (self._x_ran - np.mean(self._x_ran)) \
                               * self._fit.params['y_slope'].value
        self._param = self._fit.params
        self._out_redchi = np.sum(
                ((self._y_ran - self._y_fit) / self._dy_ran) ** 2) \
                / self._fit.nfree
        self._out_aic = self._fit.aic

        #print(self._param.items())
        #print([value for value in self._param.valuesdict().values()])
        sorted_keys = sorted(self._param.valuesdict().keys())
        param_sort = OrderedDict((key, self._param.valuesdict()[key]) \
                                 for key in sorted_keys)
        param_copy = [value for value in param_sort.values()]
        g = 0
        #s = 5  # To skip over continuum and PSF parameters
        s = 2  # To skip over continuum parameters
        n = len(self.group(line)['SYST'][line])  # Number of ions in the system
        for l in np.where(self.group_sel(line) == True)[0]: 
            self._t['X'][l] = ion_dict[self._t['SYST'][l][0]][0] \
                              * (param_copy[s + g + 3] + 1)
            self._t['Z'][l] = param_copy[s + g + 3]
            self._t['LOGN'][l] = np.log10(param_copy[s + g]) 
            self._t['B'][l] = param_copy[s + g + 1] 
            self._t['BTUR'][l] = param_copy[s + g + 2]
            g += 4 * n
        
    def group(self, line=-1):
        """ Create a group of line for each line in the list """
        iter = range(len(self._t))
        #self._t.sort('X')
        if line == -1:
            sel = np.array([self._t[self.group_sel(l)] for l in iter])
        else:
            sel = self._t[self.group_sel(line)]
        ret = copy.deepcopy(sel)
        for row in sel:
            if (len(row['SYST']) > 1):
                ref = row['SYST'][0]
                for i in range(1, len(row['SYST'])):
                    ion = row['SYST'][i]
                    xfact = ion_dict[ion][0] / ion_dict[ref][0]
                    yfact = ion_dict[ion][1] / ion_dict[ref][1]
                    xmin = row['XMIN'] * xfact 
                    xmax = row['XMAX'] * xfact
                    x = row['X'] * xfact
                    try:
                        y = np.interp(x, self._x_ran, self._y_ran)
                    except:
                        try:
                            y = self._cont - (self._cont - row['Y']) * yfact
                        except:
                            y = row['Y']
                    syst = np.append(ion, np.delete(row['SYST'][:].data, i))
                    z = row['Z']
                    logN = row['LOGN']
                    b = row['B']
                    btur = row['BTUR']
                    ret.add_row([xmin, xmax, x, y, None, 1, syst, z, logN, b,
                                 btur])
        return ret
    
    def group_sel(self, line):
        """ Define the selection indexes for group lines """
        iter = range(len(self._t))
        self._t.sort('X')
        groups = np.append([0], np.cumsum(self._t['XMAX'][:-1] <
                                          self._t['XMIN'][1:]))
        return np.array([groups[l] == groups[line] for l in iter])


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
    
    def line_add(self, xmin, xmax, x, y, syst=['Ly_a']):
        """ Add a line to a line list """
        z = x / ion_dict[syst[0]][0] - 1
        #print(self._t['Z'], z)
        #if np.isclose(self._t['Z'], z, atol=1e-7).any() == False or 1 == 1:
        self.t.add_row([xmin, xmax , x, y, None, 1, syst, z, 14.0, 20.0,
                        0.0])
        ret = True
        #else:
        #    ret = False
        return ret

    def line_delete(self, cond):
        """ Delete a line from a line list """
        where = np.where(cond)[0]
        for w in where:
            self.t.remove_row(w)
        return len(where)

    def line_merge(self, tol):
        """ Merge close lines """
        x = self.t['X'][:-1]
        x_shift = self.t['X'][1:]
        where = np.isclose(x, x_shift, atol=tol)
        l = 0
        for w in where:
            if w == True:
                #print(l, w)
                self.t['LOGN'][l] = np.log(np.power(10, self.t['LOGN'][l]) \
                                    + np.power(10, self.t['LOGN'][l + 1]))
                self.t['B'][l] = max(self.t['B'][l], self.t['B'][l + 1])
                self.t.remove_row(l + 1)    
            l += 1
            
    def model_cont(self, line, prev=None, cont=None):
        if prev == None:
            if cont == None:
                cont = self._y_extr[0] + (self._x_ran - self._x_extr[0]) \
                       * (self._y_extr[1] - self._y_extr[0]) \
                       / (self._x_extr[1] - self._x_extr[0])
            
            self._y_cont = cont
            self._y_cont0 = self._y_cont

            
    def model_trasm(self, line):
        """ Create the composite continuum + Voigt model for a system """
        
        group = self.group(line)

        # Continuum model
        norm_model = Model(func_lin)
        param = norm_model.make_params()
        y_norm = np.mean(self._y_cont)
        y_norm_fact = 1.02
        y_norm_max = np.mean(self._y_cont0) * y_norm_fact
        y_norm_min = np.mean(self._y_cont0) / y_norm_fact        
        y_slope = (self._y_cont[len(self._y_cont) - 1] - self._y_cont[0]) \
                  / (self._x_ran[len(self._x_ran) - 1] - self._x_ran[0])
        y_slope_add = 0.5
        y_slope_max = (self._y_cont0[len(self._y_cont0) - 1]
                       - self._y_cont0[0]) \
                       / (self._x_ran[len(self._x_ran) - 1] - self._x_ran[0]) \
                       * (1 + y_slope_add)
        y_slope_min = (self._y_cont0[len(self._y_cont0) - 1]
                       - self._y_cont0[0]) \
                       / (self._x_ran[len(self._x_ran) - 1] - self._x_ran[0]) \
                       * (1 - y_slope_add)
        param['y_norm'].set(y_norm, max=y_norm_max, min=y_norm_min)
        param['y_slope'].set(y_slope, max=y_slope_max, min=y_slope_min)
        model = norm_model

        z_diff = 5e-4
        logN_min = 0
        logN_max = 20
        N_min = 1
        N_max = 1e20
        b_min = 1
        b_max = 100
        btur_min = 0
        btur_max = 100
        ztag_list = []
        pref_list = []
        y_interp = np.ones(len(group))

        # First loop to identify newly added lines
        g_interp = len(group)
        for g in range(len(group)):
            if group['LOGN'][g] == 14.0:
                g_interp = g
                try:
                    y_interp[g] = np.interp(group['X'][g], self._x_ran,
                                            1 + self._y_resid / self._y_cont)
                except:
                    y_interp[g] = np.interp(group['X'][g], self._x_ran,
                                            self._y_ran / self._y_cont)
                    
        # Second loop to initialize parameters
        for g in range(len(group)):
            ion = group['SYST'][g][0]
            ref = group['SYST'][g][1]
            ztag = 'z' + str(group['Z'][g]).replace('.', '')
            pref = ztag + '_' + ion + '_'
            if pref in pref_list:
                pref += str(g) + '_'
            pref_list.append(pref)
            voigt_model = Model(func_voigt, prefix=pref, ion=ion,
                                tab=self.splrep)
            param.update(voigt_model.make_params())            

            z = group['Z'][g]
            #print(y_interp, g)
            #if len(y_interp) > 1:
            if g_interp == len(group) - 1:
                N = np.power(10, self.guess(y_interp[g]))
            else:

                #if g_interp == len(group):
                #    N = np.power(10, group['LOGN'][g])
                                        
                # Newly added lines: guess values of N
                if (y_interp[g] != 1):
                #elif g == g_interp:
                    N = np.power(10, self.guess(y_interp[g]))

                # Lines adjacent to newly added lines: N is decreased
                #elif g == g_interp - 1 or g == g_interp + 1:
                    #N = max(np.power(10, group['LOGN'][g]) \
                     #       - np.power(10, self.guess(y_interp[0])) * 0.5, 1e10)

                # Other lines: previous value
                else:
                    N = np.power(10, group['LOGN'][g])
                #N = 1e10

            b = group['B'][g] #* 0.8
            btur = group['BTUR'][g]

            # Set constraints among parameter
            z_expr = ''
            N_expr = ''
            b_expr = ''
            btur_expr = ''
            if ztag in ztag_list:
                z_expr = ztag + '_' + ref + '_z'
                N_expr = ztag + '_' + ref + '_N'
                b_expr = ztag + '_' + ref + '_b'
                btur_expr = ztag + '_' + ref + '_btur'                
            
            ztag_list.append(ztag)

            param[pref + 'z'].set(z, min=z-z_diff, max=z+z_diff, expr=z_expr)
            param[pref + 'N'].set(N, min=N_min, max=N_max, expr=N_expr)
            param[pref + 'b'].set(b, min=b_min, max=b_max, expr=b_expr)
            param[pref + 'btur'].set(btur, vary=False, expr=btur_expr)
            model *= voigt_model

        #param.pretty_print()
            
        # Convolve PSF
        """
        psf = Model(gaussian)
        conv_model = CompositeModel(model, psf, convolve)
        center = np.mean(self._x_ran)
        sigma = 1/3 * center / resol
        param.update(psf.make_params())
        param['amplitude'].set(1, vary=False)
        param['center'].set(center, vary=False)
        param['sigma'].set(sigma, vary=False)        
        """
        conv_model = model
        
        # Save results
        self._y_trasm = conv_model.eval(param, x=self._x_ran)
        self._param = param
        return conv_model

    def prep(self, line, prev=None):
        """ Prepare data structure for modeling and fitting """
        if prev == None:
            ran = self.range(line)
            self._x_ran = ran['X']
            self._y_ran = ran['Y']
            self._dy_ran = ran['DY'] #* 1.5  # Errors may be underestimated
            self._x_extr = np.array([ran['X'][0], ran['X'][len(ran) - 1]]) 
            x_max = np.array(
                np.append(self.group(line)['XMIN'], self.group(line)['XMAX']))
            y_max = np.array(np.interp(x_max, self._x_ran, self._y_ran))
            if (len(x_max) > 3): 
                x_top = np.extract(y_max > np.mean(y_max), x_max)
                y_top = np.array(np.interp(x_top, self._x_ran, self._y_ran))
            else:
                x_top = x_max
                y_top = y_max
            m, q, r, p, e = linregress(x_top, y_top)
            self._y_extr = m * self._x_extr + q
        
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
        """
        return np.logical_and(
                self._spec['X'] > min(self._t[self.group_sel(line)]['XMIN']),
                self._spec['X'] < max(self._t[self.group_sel(line)]['XMAX']))
        """
        ret = self._spec['X'] < 0.0
        for row in self.group(line):
            ret = np.logical_or(ret,
                                np.logical_and(
                                    self._spec['X'] > row['XMIN'],
                                    self._spec['X'] < row['XMAX']))
        return ret

    def tabulate(self, a_min=1e-6, a_max=1e-2, da=1e-3,
                 u_min=-1e3, u_max=1e3, du=0.01, pref=None):
        """ Tabulate the Voigt function for faster fitting """
        a_range = np.arange(a_min, a_max, da)
        u_range = np.arange(u_min, u_max, du)        
        ret = np.empty([len(u_range), len(a_range)])
        #gamma = Gamma * 1e-10 / (4 * np.pi)
        a_i = -1
        for a in a_range:
            a_i += 1
            u_i = -1
            for u in u_range:
                u_i += 1
                ret[u_i, a_i] = func_fadd(a, u)                
        self.tab = ret
        self.splrep = interp2d(a_range, u_range, ret)
        if pref != None:
            np.savetxt(pref + 'tab.dat', self.tab)
            np.savetxt(pref + 'a_range.dat', a_range)            
            np.savetxt(pref + 'u_range.dat', u_range)
