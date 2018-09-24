from astropy import units as u
from collections import OrderedDict as od


# Description
rec_descr = {'cont_line_rem': "Find Continuum by Removing Lines",
             'cont_max_smooth': "Find Continuum by Smoothing the Flux Maxima",
             'line_find': "Find Lines",
             'spec_extract': "Extract Spectral Region",
             'syst_def': "Define System",
             'syst_find': "Find Systems",
             'syst_fit': "Fit Systems"}

class Recipe():
    def __init__(self, obj, name):

        self.obj = obj
        self.name = name
        self.descr = rec_descr[name]        

        if name == 'cont_line_rem':
            self.params = {'frac': 0.03}
            self.dialog = od([('Fraction:', 'frac')])
            self.procs = ['line_rem']
                             
        if name == 'cont_max_smooth':
            self.params = {'smooth': 4.0,
                           'flux_corr': 1.0,
                           'kappa_low': 3.0,
                           'kappa_high': 3.0}
            self.dialog = od([('Smoothing:', 'smooth'),
                              ('Flux correction:', 'flux_corr'),
                              ('Lower thresh. (sigma):', 'kappa_low'),
                              ('Upper thresh. (sigma):', 'kappa_high')])
            self.procs = ['max_smooth']

        if name == 'line_find':
            self.params = {'mode': 'abs',
                           'diff': 'min',
                           'kappa': 5.0,
                           'sigma_min': 5.0,
                           'sigma_max': 100.0}
            self.dialog = od([('Mode:', 'mode'),
                              ('Difference:', 'diff'),
                              ('Threshold (sigma):', 'kappa'),
                              ('Min. sigma (sigma):', 'sigma_min'),
                              ('Max. sigma (sigma):', 'sigma_max')])
            self.procs = ['find_special']


        if name == 'spec_extract':
            self.params = {'forest': True,
                           'ion': 'Ly_a',
                           'zem': 0.0,
                           'prox_vel': 0.0,
                           'reg': False,
                           'xmin': 0.0,
                           'xmax': 0.0}
            self.dialog = od([('Use forest:', 'forest'),
                              ('Ion:', 'ion'),
                              ('Emission redshift:', 'zem'),
                              ('Prox. velocity:', 'prox_vel'),
                              ('Use region:', 'reg'), 
                              ('Min. wavelength:', 'xmin'),
                              ('Max. wavelength:', 'xmax')])
            self.procs = ['extract']

        if name == 'syst_def':
            self.params = {'forest': 'Ly',
                           'zem': 0.0,
                           'prox_vel': 0.0,
                           'xmin': 0.0,
                           'xmax': 0.0}
            self.dialog = od([('Ion:', 'forest'),
                              ('Emission redshift:', 'zem'),
                              ('Prox. velocity:', 'prox_vel'),
                              ('Min. wavelength:', 'xmin'),
                              ('Max. wavelength:', 'xmax')])
            self.procs = ['syst_def']
            
        if name == 'syst_find':
            self.params = {'series': 'CIV',
                           'ztol': 3e-4}
            self.dialog = od([('Series:', 'series'),
                              ('Redshift tolerance:', 'ztol')])
            self.procs = ['find']

    def execute(self):
        for p in zip(self.procs):
            ret = getattr(self.obj, p)(**self.params)

        return ret
