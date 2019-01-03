from . import Spec1D
from .utils import *
from astropy import units as u
from astropy.constants import c, e, m_e
from astropy.table import Table, Column
from copy import deepcopy as dc
from lmfit import Model as lmm
from lmfit import Parameters as lmp
from lmfit.lineshapes import gaussian
import numpy as np
from scipy.special import wofz
import sys

yunit = u.erg / (u.Angstrom * u.cm**2 * u.s)

def fadd_f(a, u):
    """ @brief Real part of the Faddeeva function Re(F)
    @param a First abstract variable
    @param u Second abstrac variable
    @return Re(F(a, u))
    """
    
    return np.real(wofz(u + 1j * a))

def gauss_psf(x, c_min, c_max, center, resol):
    """ @brief Gaussian PSF

    The function returns a gaussian array for each element of a selected region
    in the wavelength domain

    @param x Wavelength domain (in nm)
    @param c_min Starting pixel of the region
    @param c_max Ending pixel of the region
    @param center Center wavelength of the region
    @param resol Resolution
    @return Gaussian PSF over x
    """

    sigma = center / resol * 4.246609001e-1
    psf = np.exp(-(0.5 * (x-center) / sigma)**2)
    psf[np.where(psf < 1e-4)] = 0.0
    ret = [np.array(psf)[c_min:c_max]]
    return ret
    

def linear_f(x, norm, slope):
    """ @brief Linear function
    
    @param x Wavelength domain 
    @param norm Normalization
    @param slope Slope
    @return Linear function over x, with value norm at the center
    """
    
    return np.array(norm + (x - np.mean(x)) * slope)

def voigt_f(x, z, N, b, btur, ion='Ly_a', wave=0.0, tab=None):
    """ @brief Voigt function (real part of the Faddeeva function, after a 
    change of variables)

    @param x Wavelength domain (in nm)
    @param z Redshift
    @param N Column density (in cm^-2)
    @param b Doppler broadening (in km s^-1)
    @param btur Turbulent broadening (in km s^-1)
    @param ion Ionic transition
    @param wave Wavelength of the line (in nm)
    @param tab Table with the Faddeeva function
    @return Voigt function over x
    """
    if wave == 0.0:
        wave = dict_wave[ion].value*1e-9
        wave_z = wave*(1+z)
    else:  # For unknown species, Ly-alpha rest-frame wavelength is assumed
        wave_z = wave*1e-9
        wave = dict_wave['Ly_a'].value*1e-9
    f = dict_f[ion]
    gamma = dict_gamma[ion]
    x_si = x * 1e-9
    N_si = N * 1e4
    b_si = np.sqrt(b**2 + btur**2) * 1e3
    tau0 = N_si * np.sqrt(np.pi) * f * e.esu**2 / (m_e * c) \
           * 1e-9 * wave / b_si
    a = 0.25 * gamma * wave / (np.pi * b_si)
    u = c.value / b_si * (x_si / wave_z - 1)
    if tab == None:
        model = np.exp(-tau0.value * fadd_f(a, u)) #* yunit
    else:
        model = np.exp(-tau0.value * tab(a, u).flatten()) #* yunit
    ret = np.array(model)
    return ret


"""
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
"""
def voigt_params(syst, **kwargs):

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


    
class Model(Spec1D):

    def __init__(self, x=None, xmin=None, xmax=None, y=None, dy=None,
                 yresid=None, yadj=None, xunit=xunit_def, yunit=yunit_def):
                 #spec=None,
                 #line=None,
                 #syst=None,
                 #group=None,
                 #chunk=None):
        """ @brief Constructor for the Model class 
        """

        self._t = self.create_t(x, xmin, xmax, y, dy, yresid, yadj)
        """
        #if ((syst is None) and (line is None)):
        #    raise Exception("syst or line must be provided.")
        
        self._spec = spec
        if ((syst is None) and (line is not None)):
            self._syst = line
        elif (syst is not None):
            self._syst = syst
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
        """
        #if chunk is not None:
        #    self._model = dc(chunk)
        self._use_good = False
        
# Methods

    def create_t(self, x=None, xmin=None, xmax=None, y=None, dy=None,
                 yresid=None, yadj=None, xunit=xunit_def, yunit=yunit_def,
                 mask=None, dtype=float):
        """ @brief Create a model table

        @param x Domain
        @param y Model
        @param norm Normalization component of the model
        @return Table
        """

        t = Table(masked=True)
        if x is not None:
            t['X'] = Column(np.array(x, ndmin=1), dtype=dtype, unit=xunit)
        if xmin is not None:
            t['XMIN'] = Column(np.array(xmin, ndmin=1), dtype=dtype, unit=xunit)
        if xmax is not None:
            t['XMAX'] = Column(np.array(xmax, ndmin=1), dtype=dtype, unit=xunit)
        if y is not None:
            t['Y'] = Column(np.array(y, ndmin=1), dtype=dtype, unit=yunit)
        if dy is not None:
            t['DY'] = Column(np.array(dy, ndmin=1), dtype=dtype, unit=yunit)
        if yresid is not None:
            t['YRESID'] = Column(np.array(yresid, ndmin=1), dtype=dtype,
                                 unit=yunit)
        if yadj is not None:
            t['YADJ'] = Column(np.array(yadj, ndmin=1), dtype=dtype, unit=yunit)
        if mask is not None:
            t['X'].mask = mask

        return t
        
    def linear(self, value=[1.0, 0.0], vary=[True, False], min=[0.9, 0.0],
               max=[1.1, 0.0], expr=[None, None], pref='adj'):
        """ @brief Linear adjustment to continuum 

        The model has two parameters, norm and slope, and is multiplied to the
        continuum to correct it.

        @param value Array of guess values for the parameters
        @param vary Array of constraint on variability for the parameters
        @param min Array of minimum values for the parameters
        @param max Array of maximum values for the parameters
        @param expr Array of constraining expressions for the parameters
        @param pref Prefix to be used in composite models
        """
        
        self._linear_fun = lmm(linear_f, prefix=pref+'_')
        self._linear_par = self._linear_fun.make_params()
        for i, k in enumerate(self._linear_par.keys()):
            self._linear_par[k].set(value=value[i], vary=vary[i], min=min[i],
                                    max=max[i], expr=expr[i])
        

    def psf_gauss(self, value, vary=[False, False, False, False],
                  min=[None, None, None, None], max=[None, None, None, None],
                  expr=[None, None, None, None], pref='psf_gauss'):
        """ @brief Instrument PSF with gaussian profile

        The instrument PSF is defined over a spectral region and has four 
        parameters: c_min (starting pixel of the region), c_max (ending pixel of
        the region), center (center wavelength of the region), and resol

        @param value Array of guess values for the parameters
        @param vary Array of constraint on variability for the parameters
        @param min Array of minimum values for the parameters
        @param max Array of maximum values for the parameters
        @param expr Array of constraining expressions for the parameters
        @param pref Prefix to be used in composite models        
        """
        self._psf_gauss_fun = lmm(gauss_psf, prefix=pref+'_')
        self._psf_gauss_par = self._psf_gauss_fun.make_params()
        for i, k in enumerate(self._psf_gauss_par.keys()):
            self._psf_gauss_par[k].set(value=value[i], vary=vary[i], min=min[i],
                                       max=max[i], expr=expr[i])
                

    def voigt(self, ion, wave, value=[z_def, N_def, b_def, btur_def],
              vary=[True, True, True, False], min=[None, 1e10, 0.0, None],
              max=[None, 1e23, 200.0, None], expr=[None, None, None, None],
              pref='voigt'):
        """ @brief Voigt profile for a single line 

        The Voigt profile has four parameters: z, N, b, btur

        @param ion Name of the ionic transition
        @param wave Wavelength of the ionic transition
        @param value Array of guess values for the parameters
        @param vary Array of constraint on variability for the parameters
        @param min Array of minimum values for the parameters
        @param max Array of maximum values for the parameters
        @param expr Array of constraining expressions for the parameters
        @param pref Prefix to be used in composite models        
        """

        self._voigt_fun = lmm(voigt_f, prefix=pref+'_', ion=ion)
        self._voigt_par = self._voigt_fun.make_params()
        for i, k in enumerate(self._voigt_par.keys()):
            if i < 4:
                self._voigt_par[k].set(value=value[i], vary=vary[i], min=min[i],
                                       max=max[i], expr=expr[i])
            else:
                self._voigt_par[k].set(value=wave, vary=False)


# Deprecated


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
            mean = np.mean(self._spec.x[self._chunk[c]].value)
            center = mean
            sigma = mean / resol * 4.246609001e-1  # Factor to convert FWHM into
                                                   # standard deviation
            param[pref+'amplitude'].set(1.0, vary=False)
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



