from .utils import *
from astropy import units as u
from astropy.constants import c, e, m_e
from copy import deepcopy as dc
from lmfit import Model as lmm
from lmfit import Parameters as lmp
from lmfit.lineshapes import gaussian
import numpy as np
from scipy.special import wofz

yunit = u.erg / (u.Angstrom * u.cm**2 * u.s)

def fadd_func(a, u):
    """ Compute the real part of the Faddeeva function """
    return np.real(wofz(u + 1j * a))

def linear_func(x, slope, norm):
    return norm + (x - np.mean(x)) * slope

def linear_step_func(x, xmin, xmax, slope, norm):
    where = np.where(np.logical_and(x>=xmin, x<=xmax))[0]
    ret = x * 0.0
    ret[where] = norm + (x[where] - np.mean(x[where])) * slope
    return ret

def norm_func(x, norm):
    ret = norm
    return ret

def norm_step_func(x, xmin, xmax, norm):
    where = np.where(np.logical_and(x>=xmin, x<=xmax))[0]
    ret = x * 0.0
    ret[where] = norm
    return ret

def voigt_func(x, z, N, b, btur, ion='Ly_a', tab=None):
    """ Compute the Voigt function """
    wave = dict_wave[ion].value * 1e-9
    f = dict_f[ion]
    gamma = dict_gamma[ion]
    x_si = x * 1e-9
    N_si = N * 1e4
    b_si = np.sqrt(b**2 + btur**2) * 1e3
    tau0 = N_si * np.sqrt(np.pi) * f * e.esu**2 / (m_e * c) \
           * 1e-9 * wave / b_si
    a = 0.25 * gamma * wave / (np.pi * b_si)
    u = c.value / b_si * (x_si / (wave * (1 + z)) - 1)
    if tab == None:
        model = np.exp(-tau0.value * fadd_func(a, u)) #* yunit
    else:
        model = np.exp(-tau0.value * tab(a, u).flatten()) #* yunit
    return model


class Model():

    def __init__(self,
                 spec=None,
                 line=None,
                 syst=None,
                 group=None,
                 chunk=None):
        """ Constructor for the Fit class """

        if ((syst is None) and (line is None)):
            raise Exception("syst or line must be provided.")
        
        self._spec = spec
        if (syst is None):
            self._syst = line
        else:
            self._syst = syst
        if (hasattr(self._syst, '_cont')):
            self._precont = self._syst.cont
            self._cont = self._syst._cont
        self._group = group
        self._chunk = chunk
        
# Methods

    def N_guess(self, unabs, ion='Ly_a'):
        """ Guess the column density, given the line peak """
        
        logN_arr = np.arange(20.0, 10.0, -0.1)

        logN = []
        syst = self._syst.t
        if (self._group is not None):
            syst = self._syst.t[self._group[1]]

        for r in syst: #self._syst.t[self._group[1]]:
            if (hasattr(self._syst, '_ion')):
                try:
                    ion = r['ION'][0]
                except:
                    ion = r['ION'] 
                    y = r['Y']                   
                x = r['X']
            else:
                ion = 'Ly_a'
            try:
                y = r['Y'][0]
            except:
                y = r['Y']                   
            x = r['X']

            cont = dc(unabs)
            cont.to_z([ion])
            norm = np.mean(y / np.interp(x, cont.x, cont.y))
            voigt_arr = voigt_func(dict_wave[ion].value, 0,
                                   np.power(10, logN_arr),
                                   voigt_def['b'],
                                   voigt_def['btur'], ion=ion)
            logN_interp = np.interp(norm, voigt_arr, logN_arr)

            if (logN_interp > 14.5):
                logN_interp = 14.5
            if logN == []:
                logN = np.array([logN_interp])
            else:
                logN = np.append(logN, logN_interp)

        ret = np.power(10, logN) / u.cm**2
                
        return ret

    def norm(self, value=1.0, vary=False):

        #if (self._chunk is None):
        #    raise Exception("Chunk must be provided.")

        if (self._cont is None):
            raise Exception("Continuum be provided.")


    
        """
        for c in range(1, len(self._chunk)):
            pref = 'cont' + str(c) + '_'
            model = lmm(norm_step_func, prefix=pref)
            param = model.make_params()
            xmin = np.min(self._spec.x[self._chunk[c]].value)
            xmax = np.max(self._spec.x[self._chunk[c]].value)
            cont = self._cont.y.value
            norm = value #np.mean(cont)
            
            param[pref+'xmin'].set(xmin, vary=False)
            param[pref+'xmax'].set(xmax, vary=False)
            if (vary == False):
                param[pref+'norm'].set(norm, vary=False)
            else:
                param[pref+'norm'].set(norm, min=0.95)#, max=1.05, min=0.95)

            if (c == 1):
                ret = (model, param)
            else:
                ret += (model, param)
        """

        pref = 'cont_'
        model = lmm(norm_func, prefix=pref)
        param = model.make_params()
        cont = self._cont.y.value
        norm = value #np.mean(cont)
            
        #param[pref+'chunk_sum'].set(chunk_sum, vary=False)
        if (vary == False):
            param[pref+'norm'].set(norm, vary=False)
        else:
            param[pref+'norm'].set(norm, min=0.95)#, max=1.05, min=0.95)

        ret = (model, param)

        return ret
            
    
    def psf(self, resol): #, center, sigma):

        if (self._chunk is None):
            raise Exception("Chunk must be provided.")

        """
        for c in range(1, len(self._chunk)):
            pref = 'psf' + str(c) + '_'
            model = lmm(gaussian, prefix=pref)
            param = model.make_params()
            #print(np.sum(self._chunk[c]))
            mean = np.mean(self._spec.x[self._chunk[c]].value)
            center = mean
            sigma = mean / resol * 4.246609001e-1  # Factor to convert FWHM into
                                                   # standard deviation
            param[pref+'amplitude'].set(1.0, vary=False)
            #print("center", center)
            param[pref+'center'].set(center, vary=False)
            param[pref+'sigma'].set(sigma, vary=False)        
        
            if (c == 1):
                ret = (model, param)
            else:
                ret += (model, param)
        """
        for c in range(1, len(self._chunk)):
            if (c == 1):
                chunk_sum = dc(self._chunk[c])
            else: 
                chunk_sum += self._chunk[c]    
                   
        pref = 'psf_'
        model = lmm(gaussian, prefix=pref)
        param = model.make_params()

        # Check this: it only works if the center of the gaussian is inside
        # the chunk of interest
        center = np.median(self._spec.x[chunk_sum].value)
        sigma = center / resol * 4.246609001e-1  # Factor to convert FWHM into
                                                   # standard deviation
        param[pref+'amplitude'].set(1.0, vary=False)
        param[pref+'center'].set(center, vary=False)
        param[pref+'sigma'].set(sigma, vary=False)        
        
        ret = (model, param)    
        #"""
        return ret

    def unabs(self):
        """ Create a model of the unabsorbed continuum level for a line """
        
        if (self._chunk is None):
            raise Exception("Chunk must be provided.")

        for c in range(1, len(self._chunk)):
            pref = 'cont' + str(c) + '_'
            model = lmm(linear_step_func, prefix=pref)
            param = model.make_params()
            xmin = np.min(self._spec.x[self._chunk[c]].value)
            xmax = np.max(self._spec.x[self._chunk[c]].value)

            where = np.where(np.logical_and(
                self._syst._maxima.x.value >= xmin - 1e-7,
                self._syst._maxima.x.value <= xmax + 1e-7))
            maxima = self._syst._maxima._t[where]
            argsort = np.argsort(maxima['Y'])
            if (len(maxima) > 1):
                x_2best = maxima['X'][argsort][-2:]
                y_2best = maxima['Y'][argsort][-2:]
            elif (len(maxima) == 1):
                x_2best = np.full(2, maxima['X'][argsort][-1])
                y_2best = np.full(2, maxima['Y'][argsort][-1])
            else:
                x_2best = np.array([self._syst.xmin[self._group[c]][0].value,
                                    self._syst.xmax[self._group[c]][-1].value])
                y_2best = np.array([np.interp(x_2best[0], self._spec.x.value,
                                              self._spec.y.value),
                                    np.interp(x_2best[1], self._spec.x.value,
                                              self._spec.y.value)])
            bestsort = np.argsort(x_2best)
            x_2sort = x_2best[bestsort]
            y_2sort = y_2best[bestsort]
            
            xi = self._spec.x[self._chunk[c]][0].value
            xf = self._spec.x[self._chunk[c]][-1].value
            #xi = self._syst.xmin[self._group[c]][0].value
            #xf = self._syst.xmax[self._group[c]][-1].value            
            xm = np.mean(self._spec.x[self._chunk[c]].value)
            yi = self._spec.y[self._chunk[c]][0].value
            yf = self._spec.y[self._chunk[c]][-1].value
            #xi = self._syst.xmin[self._group[c]][0].value
            #xf = self._syst.xmax[self._group[c]][-1].value            
            yi = np.interp(xi, self._spec.x.value,
                           self._spec.y.value)
            yf = np.interp(xf, self._spec.x.value,
                           self._spec.y.value)

            slope = (yf - yi) / (xf - xi)
            norm = yi + (xm - xi) * slope

            #slope = (y_2sort[1] - y_2sort[0]) / (x_2sort[1] - x_2sort[0])
            #norm = y_2sort[0] + (xm - x_2sort[0]) * slope
            param[pref+'xmin'].set(xmin, vary=False)
            param[pref+'xmax'].set(xmax, vary=False)
            param[pref+'slope'].set(slope, vary=False)
                                    # min=slope/unabs_fact['slope'],
                                    #max=slope*unabs_fact['slope'])
            param[pref+'norm'].set(norm)#, vary=False)
                                   #min=norm/unabs_fact['norm'],
                                   #max=norm*unabs_fact['norm'])
            if (c == 1):
                ret = (model, param)
            else:
                ret += (model, param)

        return ret

    def voigt(self, z, N, b, btur, ion):
        """ Create a Voigt model for a line """

        self._z = dc(z)
        self._N = dc(N)
        self._b = dc(b)
        self._btur = dc(btur)        
        
        z_list = []
        pref_list = []
        expr_dict = {}
        i = 0
        
        """
        for c in range(1, len(self._chunk)):
            for l in range(len(self._syst.t[self._group[1]])):
                pref = 'voigt' + str(c) + '_z' + str(z[l]).replace('.', '') \
                       + '_'
                #print(c, l, pref, self._syst.t[self._group[1]][l])
                if (z[l] in z_list):
                    expr = np.asarray(pref_list)[
                        np.where(np.asarray(z_list)==z[l])][0]
                    i += 1
                    pref = pref + str(i) + '_'
                    z_expr = expr + 'z'
                    N_expr = expr + 'N'
                    b_expr = expr + 'b'
                    btur_expr = expr + 'btur'
                else:
                    z_list.append(z[l])
                    pref_list.append(pref)                    
                    z_expr = ''
                    N_expr = ''
                    b_expr = ''
                    btur_expr = ''
                expr_dict[pref+'z'] = z_expr 
                expr_dict[pref+'N'] = N_expr 
                expr_dict[pref+'b'] = b_expr 
                expr_dict[pref+'btur'] = btur_expr 

                if (type(ion) is str):
                    model_l = lmm(voigt_func, prefix=pref, ion=ion)
                else:
                    model_l = lmm(voigt_func, prefix=pref, ion=ion[c-1])
                if (l == 0):
                    model = model_l
                    param = model_l.make_params()
                else:
                    model *= model_l
                    param.update(model_l.make_params())
                param[pref+'z'].set(z[l].value, min=z[l].value/z_fact,
                                    max=z[l].value*z_fact)#, expr=str(z_expr))
                param[pref+'N'].set(N[l].value, min=voigt_min['N'],
                                    max=voigt_max['N'])#, expr=N_expr)
                param[pref+'b'].set(b[l].value, min=voigt_min['b'],
                                    max=voigt_max['b'])#, expr=b_expr)
                param[pref+'btur'].set(btur[l].value, min=voigt_min['btur'],
                                       max=voigt_max['btur'])#, expr=btur_expr)
            if (c == 1):
                ret = (model, param)
            else:
                ret += (model, param)

        #ret = (model, param)
        """
        z_list = []
        pref_list = []
        expr_dict = {}
        mult_old = ''
        i = 0
        for l in range(len(self._syst.t[self._group[1]])):
            ion = np.sort(self._syst.ion[self._group[1]][l])
            for c in range(len(ion)):
                mult = ion[c].split('_')[0]
                pref = 'voigt' + str(c) + '_z' + str(z[l]).replace('.', '') \
                       + '_'
                if (z[l] in z_list):
                    expr = np.asarray(pref_list)[
                        np.where(np.asarray(z_list)==z[l])][0]
                    i += 1
                    pref = pref + str(i) + '_'
                    z_expr = expr + 'z'
                    N_expr = expr + 'N'
                    b_expr = expr + 'b'
                    btur_expr = expr + 'btur'
                else:
                    z_list.append(z[l])
                    pref_list.append(pref)                    
                    z_expr = ''
                    N_expr = ''
                    b_expr = ''
                    btur_expr = ''
                expr_dict[pref+'z'] = z_expr 
                expr_dict[pref+'N'] = N_expr 
                expr_dict[pref+'b'] = b_expr 
                expr_dict[pref+'btur'] = btur_expr 
                if (type(ion) is str):
                    model_l = lmm(voigt_func, prefix=pref, ion=ion)
                else:
                    model_l = lmm(voigt_func, prefix=pref, ion=ion[c-1])
                if ((l == 0) and (c == 0)):
                    model = model_l
                    param = model_l.make_params()
                else:
                    model *= model_l
                    param.update(model_l.make_params())
                param[pref+'z'].set(z[l].value, min=z[l].value/z_fact,
                                    max=z[l].value*z_fact)#, expr=str(z_expr))
                param[pref+'N'].set(N[l].value, min=voigt_min['N'],
                                    max=voigt_max['N'])#, expr=N_expr)
                param[pref+'b'].set(b[l].value, min=voigt_min['b'],
                                    max=voigt_max['b'])#, expr=b_expr)
                param[pref+'btur'].set(btur[l].value, min=voigt_min['btur'],
                                       max=voigt_max['btur'])#, expr=btur_expr)
                
                if (mult == mult_old):
                    for k in expr_dict:
                        param[k].set(expr=expr_dict[k])    
                #param.pretty_print()
                mult_old = mult

        ret2 = (model, param)


        #if (hasattr(self._syst, 'ion')):
        #    ret += (ret, expr_dict)


        return ret2
