from .utils import *
from astropy import units as u
from astropy.constants import c, e, m_e
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
        self._group = group
        self._chunk = chunk
        
# Methods

    def N_guess(self, unabs, ion='Ly_a'):
        """ Guess the column density, given the line peak """
        
        logN_arr = np.arange(20.0, 10.0, -0.1)

        logN = []
        for r in self._syst.t[self._group[1]]:
            if (hasattr(self._syst, '_ion')):
                try:
                    ion = r['ION'][0]
                except:
                    ion = r['ION']                    
            else:
                ion = 'Ly_a'
            norm = np.mean(np.asarray(r['Y']) \
                           / np.interp(r['X'], self._spec.x, unabs.y))
            #norm = norm * 
            voigt_arr = voigt_func(dict_wave[ion].value, 0,
                                   np.power(10, logN_arr),
                                   voigt_def['b'],
                                   voigt_def['btur'], ion=ion)
            if logN == []:
                logN = np.array([np.interp(norm, voigt_arr, logN_arr)])
            else:
                logN = np.append(logN, np.interp(norm, voigt_arr, logN_arr))

        return np.power(10, logN) / u.cm**2

    def psf(self, resol): #, center, sigma):

        if (self._chunk is None):
            raise Exception("Chunk must be provided.")

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
            xi = self._spec.x[self._chunk[c]][0].value
            xf = self._spec.x[self._chunk[c]][-1].value
            xm = np.mean(self._spec.x[self._chunk[c]].value)
            yi = self._spec.y[self._chunk[c]][0].value
            yf = self._spec.y[self._chunk[c]][-1].value
            slope = (yf - yi) / (xf - xi)
            norm = yi + (xm - xi) * slope
            param[pref+'xmin'].set(xmin, vary=False)
            param[pref+'xmax'].set(xmax, vary=False)
            param[pref+'slope'].set(slope) #, vary=False)
                                    # min=slope/unabs_fact['slope'],
                                    #max=slope*unabs_fact['slope'])
            param[pref+'norm'].set(norm) #, vary=False)
                                   #min=norm/unabs_fact['norm'],
                                   #max=norm*unabs_fact['norm'])
            if (c == 1):
                ret = (model, param)
            else:
                ret += (model, param)

        return ret

    def voigt(self, z, N, b, btur, ion):
        """ Create a Voigt model for a line """

        z_list = []
        pref_list = []
        expr_dict = {}
        i = 0
        for c in range(1, len(self._chunk)):
            for l in range(len(self._syst.t[self._group[1]])):
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

        if (hasattr(self._syst, 'ion')):
            ret += (ret, expr_dict)

        return ret
