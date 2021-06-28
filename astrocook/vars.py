import ast
from astropy import units as au
from astropy import constants as aconst
from astropy.io import ascii, fits
import numpy as np
import os
import operator as op
import pathlib
 #c, e, m_e

"""
ac_ops = {'>': operator.gt,
          '<': operator.lt,
          '==': operator.eq,
          '+': operator.add,
          '-': operator.sub,
          '*': operator.mul,
          '/': operator.truediv,  # use operator.div for Python 2
          '%': operator.mod,
          '^': operator.xor}
"""
py_ops = {ast.Add: op.add,
          ast.Sub: op.sub,
          ast.Mult: op.mul,
          ast.Div: op.truediv,
          ast.Pow: op.pow,
          ast.Lt: op.lt,
          ast.Gt: op.gt,
          ast.BitXor: op.xor,
          ast.USub: op.neg}

np_ops = {ast.Add: op.add,
          ast.Sub: op.sub,
          ast.Mult: op.mul,
          ast.Div: op.truediv,
          ast.Pow: op.pow,
          ast.Lt: np.less,
          ast.Gt: op.gt,
          ast.BitXor: op.xor,
          ast.USub: op.neg}

xunit_def = au.nm
yunit_def = au.erg / (au.Angstrom * au.cm**2 * au.s)
zunit_def = au.nm / au.nm
Nunit_def = 1 / au.cm**2
bunit_def = au.km / au.s

equiv_w_v = [(au.nm, au.km/au.s,
              lambda x: np.log(x/121.567)*aconst.c.to(au.km/au.s),
              lambda x: np.exp(x/aconst.c.to(au.km/au.s).value)*121.567)]


logN_def = 14
b_def = 10

resol_def = None
max_nfev_def = 100

hwin_def = 250.0

seq = ['spec', 'nodes', 'lines', 'systs', 'mods']
seq_menu = seq + ['y_conv', 'cont', 'z0']
graph_sel = [#'spec_x_y',
             #'spec_x_y_det',
             #'lines_x_y', 'spec_x_cont', 'spec_x_model', 'spec_x_yfitmask',
             #'systs_z_series',
             'spec_h2o_reg'
             ]
graph_cols_sel = ''

graph_elem="spec,x,y,None,step,-,1,C0,1\n"\
           "spec,x,dy,None,step,-,1,C1,0.5\n"\
           "lines,x,y,None,scatter,+,1.5,C2,1\n"\
           "nodes,x,y,None,scatter,o,1,C3,1\n"\
           "spec,x,cont,None,plot,-,1,C8,1\n"\
           "spec,x,model,None,plot,-,1,C9,1\n"\
           "spec,x,model,fit_mask,plot,-,3,C9,0.5\n"\
           "systs,z,None,None,axvline,--,0.8,C2,1.0"

pars_std_d =  {
    'z': 0.0, 'logN': 13, 'b': 10.0, 'btur': 0.0, 'resol': 35000,
    'z_vary': True, 'logN_vary': True, 'b_vary': True, 'btur_vary': False, 'resol_vary': False,
    'z_min': 1e-3, 'logN_min': 10, 'b_min': 1.0, 'btur_min': 0.0, 'resol_min': 0,
    'z_max': 1e-3, 'logN_max': 18, 'b_max': 100.0, 'btur_max': 100.0, 'resol_max': 1e6,
#    'z_max': 1e-3, 'logN_max': 18, 'b_max': 200.0, 'btur_max': 200.0, 'resol_max': 1e6,
#    'z_max': 1e-3, 'logN_max': 22, 'b_max': 1000.0, 'btur_max': 200.0, 'resol_max': 1e6,
    'z_expr': None, 'logN_expr': None, 'b_expr': None, 'btur_expr': None, 'resol_expr': None}



# Default values for continuum adjustment parameters
adj_gauss_d = {
    'z': 0.0, 'ampl': 0.0, 'sigma': 0.01,
    'z_vary': False, 'ampl_vary': False, 'sigma_vary': False,
    'z_min': 0.0, 'ampl_min': -0.05, 'sigma_min': 0.0,
    'z_max': 10.0, 'ampl_max': 0.05, 'sigma_max': 1.0,
    'z_expr': None, 'ampl_expr': None, 'sigma_expr': None}


# Default values for line Voigt parameters
lines_voigt_d = {
    'z': 0.0, 'N': 1.e13, 'b': 5.0, 'btur': 0.0,
    'z_vary': True, 'N_vary': True, 'b_vary': True, 'btur_vary': False,
    'z_min': 0.0, 'N_min': 1.e11, 'b_min': 1.0, 'btur_min': 0.0,
    'z_max': 10.0, 'N_max': 1.e22, 'b_max': 100.0, 'btur_max': 100.0,
    'z_expr': None, 'N_expr': None, 'b_expr': None, 'btur_expr': None}

# Default values for PSF gaussian Parameters
psf_gauss_d = {
    'z': 0.0, 'resol': 35000,
    'z_vary': False, 'resol_vary': False,
    'z_min': 0.0, 'resol_min': 0,
    'z_max': 10.0, 'resol_max': 1e6,
    'z_expr': None, 'resol_expr': None}

forbidden_keywords = ['XTENSION', 'BITPIX', 'PCOUNT', 'GCOUNT', 'TFIELDS',
                      'NAXIS', 'TTYPE', 'TFORM', 'TUNIT', 'TDISP']

x_col_names = np.array(['x', 'wave', 'WAVE', 'col1', 'lambda'])
y_col_names = np.array(['y', 'flux', 'FLUX', 'col2'])
dy_col_names = np.array(['dy', 'err', 'ERR', 'fluxerr', 'FLUXERR', 'error', 'ERROR', 'col3'])

h2o_reg = np.array([[1350, 1450], [1800, 1950], [2500, 3400]])

filt_x_skymap = {'u': 348.0, 'v': 382.5, 'g': 493.0, 'r': 629.4, 'i': 770.2,
                 'z': 923.6}
zero_point_skymap = {'u': 29.005687, 'v': 28.481306, 'g': 29.55393973,
                     'r': 29.0143769, 'i': 28.4342476, 'z': 27.9522006}

p = '/'.join(pathlib.PurePath(os.path.realpath(__file__)).parts[0:-1]) + '/../'
atom_par = ascii.read(pathlib.Path(p+'/atom_par.dat'))

xem_d = {k: v*au.nm for (k, v) in atom_par['col1', 'col2']}
fosc_d = {k: v for (k, v) in atom_par['col1', 'col3']}
gamma_d = {k: v for (k, v) in atom_par['col1', 'col4']}

telluric = fits.open(pathlib.Path(p+'/telluric.fits'))[1].data

pars_d = {'lines_voigt_d': lines_voigt_d,
          'psf_gauss_d': psf_gauss_d}

trans_d = np.array([t for t in atom_par['col1'] if '-' not in t])
series_d = {k: None for k in np.unique([a.split('_')[0] for a in atom_par['col1']])}
for s in series_d:
    series_d[s] = [a for a in atom_par['col1'] if a.split('_')[0]==s]

#trans_d_short = ['SiIV_1393', 'SiIV_1402', 'SiII_1526', 'CIV_1548', 'CIV_1550', 'AlII_1670', 'NiII_1741', 'NiII_1751', 'AlIII_1854', 'AlIII_1862', 'FeII_2344', 'FeII_2374', 'FeII_2382', 'MnII_2576', 'FeII_2586', 'MnII_2594', 'FeII_2600', 'MnII_2606', 'MgII_2796', 'MgII_2803']
trans_d_short = ['CIV_1548', 'CIV_1550', 'MgII_2796', 'MgII_2803', 'SiIV_1393', 'SiIV_1402', 'AlIII_1854', 'AlIII_1862', 'FeII_2586', 'FeII_2600']
trans_d_short = ['CIV_1548', 'CIV_1550']
series_d_short = {}
for (k,v) in series_d.items():
    for vi in v:
        if vi in trans_d_short:
            if k not in series_d_short:
                series_d_short[k] = []
            series_d_short[k].append(vi)

series_d['Ly-a'] = ['Ly_a']
series_d['Ly-ab'] = ['Ly_b', 'Ly_a']
series_d['Ly-abg'] = ['Ly_g', 'Ly_b', 'Ly_a']

series_d_old = {'Ly': ['Ly_15', 'Ly_14', 'Ly_13', 'Ly_12', 'Ly_11', 'Ly_10',
                       'Ly_9', 'Ly_8', 'Ly_7', 'Ly_6', 'Ly_e', 'Ly_d', 'Ly_g',
                       'Ly_b', 'Ly_a'],
               'Ly-abg': ['Ly_g', 'Ly_b', 'Ly_a'],
               'Ly-ab': ['Ly_b', 'Ly_a'],
               'Ly-a': ['Ly_a'],
               'NV': ['NV_1238', 'NV_1242'],
               'SiII': ['SiII_1260', 'SiII_1304', 'SiII_1526'],
               'OI': ['OI_1302'],
               'SiIV': ['SiIV_1393', 'SiIV_1402'],
               'CIV': ['CIV_1548', 'CIV_1550'],
               'CIV-1548': ['CIV_1548'],
               'CIV-1550': ['CIV_1550'],
               'FeII': ['FeII_2344', 'FeII_2374', 'FeII_2382', 'FeII_2586', 'FeII_2600'],
               'FeII-2382': ['FeII_2344', 'FeII_2374', 'FeII_2382'],
               'FeII-2586': ['FeII_2586', 'FeII_2600'],
               'MgII': ['MgII_2796', 'MgII_2803'],
               'CaII': ['CaII_3934', 'CaII_3969'],
               'unknown': ['unknown']}

xem_d_old = {'Ly_a': 121.567 * au.nm,
          'Ly_b': 102.5722200 * au.nm,
          'Ly_g': 97.2536700 * au.nm,
             'Ly_d': 94.9743000 * au.nm,
             'Ly_e': 93.7803400 * au.nm,
             'Ly_6': 93.0748200 * au.nm,
             'Ly_7': 92.6225600 * au.nm,
             'Ly_8': 92.3150300 * au.nm,
             'Ly_9': 92.0963000 * au.nm,
             'Ly_10': 91.9351300 * au.nm,
             'Ly_11': 91.8129300 * au.nm,
             'Ly_12': 91.7180500 * au.nm,
             'Ly_13': 91.6429100 * au.nm,
             'Ly_14': 91.5823800 * au.nm,
             'Ly_15': 91.5328900 * au.nm,
             'Ly_lim': 91.18 * au.nm,
             'NV_1238': 123.8821 * au.nm,
             'NV_1242': 124.2804 * au.nm,
             'SiII_1260': 126.04221 * au.nm,
             'SiII_1304': 130.43702 * au.nm,
             'SiII_1526': 152.6707 * au.nm,
             'OI_1302': 130.21685 * au.nm,
             'CII_1334': 133.45323 * au.nm,
             'SiIV_1393': 139.37602 * au.nm,
             'SiIV_1402': 140.27729 * au.nm,
             'CIV': 154.94925 * au.nm,
             'CIV_1548': 154.8204 * au.nm,
             'CIV_1550': 155.0781 * au.nm,
             'AlII_1670': 167.078861 * au.nm,
             'FeII_2344': 234.42139 * au.nm,
             'FeII_2374': 237.44612 * au.nm,
             'FeII_2382': 238.27652 * au.nm,
             'FeII_2586': 258.665 * au.nm,
             'FeII_2600': 260.01729 * au.nm,
             'MgII_2796': 279.63543 * au.nm,
             'MgII_2803': 280.35315 * au.nm,
             'CaII_3934': 393.4775 * au.nm,
             'CaII_3969': 396.95901 * au.nm}

# Ionic oscillator strengths
fosc_d_old = {'Ly_a': 0.416,
          'Ly_b': 0.0791000,
          'Ly_g': 0.0290100,
          'Ly_d': 0.0139000,
          'Ly_e': 0.0078000,
          'Ly_6': 0.0048100,
          'Ly_7': 0.0031850,
          'Ly_8': 0.0022170,
          'Ly_9': 0.0016060,
          'Ly_10': 0.0012010,
          'Ly_11': 0.0009219,
          'Ly_12': 0.0007231,
          'Ly_13': 0.0005777,
          'Ly_14': 0.0004689,
          'Ly_15': 0.0003858,
          'NV_1238': 0.156,
          'NV_1242': 0.777,
          'SiII_1260': 1.1799999,
          'SiII_1304': 0.0863,
          'CII_1334': 0.128,
          'SiII_1526': 0.133,
          'OI_1302': 0.048,
          'SiIV_1393': 0.513,
          'SiIV_1402': 0.254,
          'CIV_1548': 0.1899,
          'CIV_1550': 0.09475,
          'AlII_1670': 1.74,
          'FeII_2344': 0.114,
          'FeII_2374': 0.0313,
          'FeII_2382': 0.32,
          'FeII_2586': 0.0691,
          'FeII_2600': 0.239,
          'MgII_2796': 0.6155,
          'MgII_2803': 0.3058,
          'CaII_3934': 0.6267,
          'CaII_3969': 0.3116,
          'neb': 0.1,
          'unknown': 0.416}

# Ionic damping lengths
gamma_d_old = {'Ly_a': 6.265e+08,
              'Ly_b': 1.8970e+08,
              'Ly_g': 8.1260e+07,
              'Ly_d': 4.2040e+07,
              'Ly_e': 2.4500e+07,
              'Ly_6': 1.2360e+07,
              'Ly_7': 8.2550e+06,
              'Ly_8': 5.7850e+06,
              'Ly_9': 4.2100e+06,
              'Ly_10': 3.1600e+06,
              'Ly_11': 2.4320e+06,
              'Ly_12': 1.9110e+06,
              'Ly_13': 1.5290e+06,
              'Ly_14': 1.2430e+06,
              'Ly_15': 1.0240e+06,
              'NV_1238': 3.391e+08,
              'NV_1242': 3.356e+08,
              'SiII_1260': 2.95e+09,
              'SiII_1304': 1.01e+09,
              'SiII_1526': 1.13e+09,
              'OI_1302': 5.65e+08,
              'CII_1334': 2.88e+08,
              'SiIV_1393': 8.8e+08,
              'SiIV_1402': 8.62e+08,
              'CIV_1548': 2.643e+08,
              'CIV_1550': 2.628e+08,
              'AlII_1670': 1.39e+09,
              'FeII_2344': 2.68e+08,
              'FeII_2374': 3.090e+08,
              'FeII_2382': 3.13e+08,
              'FeII_2586': 2.72e+08,
              'FeII_2600': 2.70e+08,
              'MgII_2796': 2.625e+08,
              'MgII_2803': 2.595e+08,
              'CaII_3934': 1.444e+08,
              'CaII_3969': 1.409e+08,
              'neb': 5e8,
              'unknown': 6.265e+08}


log_seed = {'set_menu': []}
