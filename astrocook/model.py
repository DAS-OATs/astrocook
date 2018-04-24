from . import Spec1D
from .utils import *
from astropy import units as u
from astropy.constants import c, e, m_e
from copy import deepcopy as dc
from lmfit import Model as lmm
from lmfit import Parameters as lmp
from lmfit.lineshapes import gaussian
import numpy as np
from scipy.special import wofz
import sys

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
    ret = np.full(len(x), norm)
    return ret

def norm_step_func(x, xmin, xmax, norm):
    where = np.where(np.logical_and(x>=xmin, x<=xmax))[0]
    ret = x * 0.0
    ret[where] = norm
    return ret

def psf_func(x, c_min, c_max, center, resol):
    sigma = center / resol * 4.246609001e-1
    psf = np.exp(-(0.5 * (x-center) / sigma)**2)
    psf[np.where(psf < 1e-4)] = 0.0
    ret = [np.array(psf)[c_min:c_max]]
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
    ret = np.array(model)
    return ret

def voigt_params(syst, **kwargs):
    """ Read voigt parameters, provided as kwargs, and format them as arrays
    of the same length of the system to be fitted """

    try:
        z = kwargs['z']
        try:
            z_arr = [z.value] * len(syst.t) * z.unit
        except:
            z_arr = [z] * len(syst.t) * u.nm/u.nm            
    except:
        z_arr = []
    try:
        N = kwargs['N']
        try:
            N_arr = [N.value] * len(syst.t) * N.unit
        except:
            N_arr = [N] * len(syst.t) / u.cm**2            
    except:
        N_arr = []
    try:
        b = kwargs['b']
        try:
            b_arr = [b.value] * len(syst.t) * b.unit
        except:
            b_arr = [b] * len(syst.t) * u.km/u.s            
    except:
        b_arr = []
    try:
        btur = kwargs['btur']
        try:
            btur_arr = [btur.value] * len(syst.t) * btur.unit
        except:
            btur_arr = [btur] * len(syst.t) * u.km/u.s            
    except:
        btur_arr = []

    ret = {'z': z_arr, 'N': N_arr, 'b': b_arr, 'btur': btur_arr}
    return ret

class Model():

    def __init__(self,
                 spec=None,
                 line=None,
                 syst=None,
                 group=None,
                 chunk=None):
        """ Constructor for the Fit class """

        #if ((syst is None) and (line is None)):
        #    raise Exception("syst or line must be provided.")
        
        self._spec = spec
        if ((syst is None) and (line is not None)):
            self._syst = line
        elif (syst is not None):
            self._syst = syst
        #print(self._syst.t)    
        if (hasattr(self, '_syst')):
            
            if (hasattr(self._syst, '_cont')):
                self._precont = self._syst.cont
                self._cont = self._syst._cont
            if (hasattr(self._syst, '_group')):            
                self._group = group
            if (hasattr(self._syst, '_chunk')):            
                self._chunk = chunk
        else:
            if (hasattr(self._spec, '_cont')):
                self._precont = self._spec.cont
                self._cont = self._spec._cont
            if (hasattr(self._spec, '_group')):            
                self._group = group
            if (hasattr(self._spec, '_chunk')):            
                self._chunk = chunk
        
# Methods

    def N_guess(self, unabs, ion='Ly_a'):
        """ Guess the column density, given the line peak """
        
        logN_arr = np.arange(20.0, 10.0, -0.1)

        logN = []
        
        syst = self._syst.t
        if (self._group is not None):
            syst = self._syst.t[self._group[1]]

        r_i = 0
        for r in syst: #self._syst.t[self._group[1]]:
            if (hasattr(self._syst, '_ion')):
                try:
                    for i in range(len(r['Y']) - len(r['Y']) + 1):
                        ion = r['ION'][i]
                except:
                    try:
                        ion = r['ION']
                    except:
                        ion = self._syst._ion[r_i]
                    #y = r['Y']                   
                #x = r['X']
            else:
                ion = 'Ly_a'
            try:
                y = r['Y'][0]
            except:
                y = r['Y']                   
            x = r['X']

            r_i += 1
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
            raise Exception("Continuum must be provided.")
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

    def norm_new(self, level=1.0, vary=[True], expr=[None], pref='norm'):
        """ Normalization of the continuum """
        
        self._norm_fun = lmm(norm_func, prefix=pref+'_')
        self._norm_par = self._norm_fun.make_params()
        self._norm_par[pref+'_norm'].set(level, vary=vary[0], expr=expr[0])
        
            
    def prof(self, spec, ion, wave_step=1e-3*u.nm, width=0.2 * u.nm, #0.03*u.nm,
             prof='voigt', ion_mask=None, **kwargs):
        """Create a window with a model profile of a line or multiplet"""

        temp = dc(spec)
        temp._norm_guess = self.norm()
        if (prof == 'voigt'):
            try:
                z = kwargs['z']
            except:
                z = u.Quantity(0.0)
            try:
                N = kwargs['N']
            except:
                N = voigt_def['N'] / u.cm**2
            try:
                b = kwargs['b']
            except:
                b = voigt_def['b'] * u.km/u.s 
            try:
                btur = kwargs['btur'] 
            except:
                btur = voigt_def['btur'] * u.km/u.s
                
            max_size = 0
            max_size = np.size(z) if np.size(z) > max_size else max_size
            max_size = np.size(N) if np.size(N) > max_size else max_size
            max_size = np.size(b) if np.size(b) > max_size else max_size
            max_size = np.size(btur) if np.size(btur) > max_size else max_size
            
            z = np.full(max_size, z) * z.unit if np.size(z)==1 else z
            N = np.full(max_size, N) * N.unit if np.size(N)==1 else N
            b = np.full(max_size, b) * b.unit if np.size(b)==1 else b   
            btur = np.full(max_size, btur) * btur.unit if np.size(btur)==1 \
                   else btur

            # Redshifts are translated to be around zero
            z = z - np.mean(z)

            ion_prof = ion
            if (ion_mask is not None):
                if (np.sum(ion_mask) == 0):
                    N = N * 0.0
                else:
                    ion_prof = np.array(ion)[ion_mask]

            temp._prof_guess = self.voigt(z, N, b, btur, ion=ion_prof)    
            voigt = temp.model(psf=False)
            
            x = np.array([])
            for p in range(len(ion)):
                wave_min = dict_wave[ion[p]] * (1 + z[0]) - width
                wave_max = dict_wave[ion[p]] * (1 + z[-1]) + width
            
                x = np.append(x, np.arange(wave_min.value, wave_max.value,
                                           wave_step.value))
            y = voigt[0].eval(voigt[1], x=x)
            ret = Spec1D(x, y, xunit=temp.x.unit, yunit=temp.y.unit)
        else:
            raise Exception("Only Voigt profile is supported.")
        
        return ret
    
    def prof_mult(self, spec, ion, plot=False, **kwargs):
        ion_mask = np.full(len(ion), False)
        prof = []
        prof.append(self.prof(spec, ion, **kwargs))       
        f_list = np.array([dict_f[ion[i]] for i in range(len(ion))]) 
        f_sort = np.argsort(f_list)[::-1]
        for i in range(len(ion)):
            ion_mask_temp = dc(ion_mask)
            ion_mask_temp[f_sort[i]] = True
            #print(i, ion, ion_mask_temp)
            #ion_mask_temp[i] = True
            prof.append(self.prof(spec, ion, ion_mask=ion_mask_temp, **kwargs))
        prof.append(self.prof(spec, ion, ion_mask=ion_mask, **kwargs))        

        if (plot == True):
            for p in range(len(prof)):
                prof[p].plot()        

        return prof
    
    def psf(self): #, center, sigma):

        """
        for c in range(1, len(self._chunk)):
            pref = 'psf' + str(c) + '_'
            model = lmm(gaussian, prefix=pref)
            param = model.make_params()
            #print(np.sum(self._chunk[c]))
            mean = np.mean(self._spec.x[self._chunk[c]].value)
            center = mean
            print(center)
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
        x = self._spec.x[chunk_sum]
        if (len(x) % 2 == 0):
            x = x[:-1]
        resol = self._spec.t['RESOL'][int((len(x)-1)*0.5)]
        #resol = 20000000
        center = np.median(x.value)
        sigma = center / resol * 4.246609001e-1  # Factor to convert FWHM into
                                                   # standard deviation
        param[pref+'amplitude'].set(1.0, vary=False)
        param[pref+'center'].set(center, vary=False)
        param[pref+'sigma'].set(sigma, vary=False)        
        
        ret = (model, param)    
        #"""
        return ret

    def psf_new(self, center=None, resol=None, vary=[False, False],
                expr=[None, None], pref='psf'):
        """ Instrumental PSF """

        sigma = center / resol * 4.246609001e-1  # Convert FWHM into st. dev.
        self._psf_fun = lmm(gaussian, prefix=pref+'_')
        self._psf_par = self._psf_fun.make_params()
        self._psf_par[pref+'_amplitude'].set(1.0, vary=False)
        self._psf_par[pref+'_center'].set(center, vary=vary[0], expr=expr[0])
        self._psf_par[pref+'_sigma'].set(sigma, vary=vary[1], expr=expr[1])
        
    def psf_new2(self, c_min, c_max, center, resol, vary=False, expr=None,
                 pref='psf'):
        self._psf_fun = lmm(psf_func, prefix=pref+'_')
        self._psf_par = self._psf_fun.make_params()
        self._psf_par[pref+'_c_min'].set(c_min, vary=False)
        self._psf_par[pref+'_c_max'].set(c_max, vary=False)
        self._psf_par[pref+'_center'].set(center, vary=False)
        self._psf_par[pref+'_resol'].set(resol, vary=vary,
                                         min=resol/1.1, max=resol*1.1,
                                         expr=expr)
   
    def psf2(self):
        pref = 'psf_'
        model = lmm(psf_func, prefix=pref)
        param = model.make_params() 
        #param[pref+'resol'].set(resol, vary=False)        
        ret = (model, param)    

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

    def voigt_new(self, ion, z, N, b, btur,
                  vary=[True, True, True, False],
                  expr=[None, None, None, None], pref='voigt'):
        """ Voigt profile """

        self._prof_fun = lmm(voigt_func, prefix=pref+'_', ion=ion)
        self._prof_par = self._prof_fun.make_params()
        #self._prof_par[pref+'_chunk_lim'].set(chunk_lim, vary=False)
        self._prof_par[pref+'_z'].set(z, vary=vary[0], min=0, expr=expr[0])
        self._prof_par[pref+'_N'].set(N, vary=vary[1], min=0, expr=expr[1])
        self._prof_par[pref+'_b'].set(b, vary=vary[2], min=0, expr=expr[2])
        self._prof_par[pref+'_btur'].set(btur, vary=vary[3], min=0,
                                         expr=expr[3])

    def voigt(self, z, N, b, btur, z_vary=True, N_vary=True, b_vary=True,
              btur_vary=False, ion=['Ly_a']):
        """ Create a Voigt model for a line """

        self._z = dc(z)
        self._N = dc(N)
        self._b = dc(b)
        self._btur = dc(btur)        
        
        z_list = []
        pref_list = []
        expr_dict = {}
        mult_old = ''
        i = 0
        ran = 1
        #if (hasattr(self, '_syst')):
        if (len(ion) > 1 or 1==1):
            #ran = len(self._syst.t[self._group[1]])
            ran = len(z)
        for l in range(ran):
            try:
                ion = np.sort(self._syst.ion[self._group[1]][l])
            except:
                pass
            for c in range(len(ion)):
                mult = ion[c].split('_')[0]
                pref = 'voigt' + str(c) + '_z' \
                       + str(z[l]).replace('.', '').replace('-', 'm') + '_'
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
                N_min = 0 if N[l].value == 0 else voigt_min['N']
                param[pref+'z'].set(z[l].value, min=z[l].value/z_fact,
                                    max=z[l].value*z_fact)#, vary=z_vary)#, expr=str(z_expr))
                param[pref+'N'].set(N[l].value, min=N_min, #min=voigt_min['N'],
                                    max=voigt_max['N'])#, vary=N_vary)#, expr=N_expr)
                param[pref+'b'].set(b[l].value, min=voigt_min['b'],
                                    max=voigt_max['b'])#, vary=b_vary)#, expr=b_expr)
                param[pref+'btur'].set(btur[l].value, min=voigt_min['btur'],
                                       max=voigt_max['btur'], vary=False)
                #, expr=btur_expr)
                
                if (mult == mult_old):
                    for k in expr_dict:
                        param[k].set(expr=expr_dict[k])    
                mult_old = mult

        ret = (model, param)


        return ret

