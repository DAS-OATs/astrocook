from astropy import units as au
 #c, e, m_e

xunit_def = au.nm
yunit_def = au.erg / (au.Angstrom * au.cm**2 * au.s)
zunit_def = au.nm / au.nm
Nunit_def = 1 / au.cm**2
bunit_def = au.km / au.s

pars_std_d =  {
    'z': 0.0, 'logN': 13, 'b': 10.0, 'btur': 0.0, 'resol': 35000,
    'z_vary': True, 'logN_vary': True, 'b_vary': True, 'btur_vary': False, 'resol_vary': False,
    'z_min': None, 'logN_min': 10, 'b_min': 1.0, 'btur_min': None, 'resol_min': None,
    'z_max': None, 'logN_max': 17, 'b_max': 100.0, 'btur_max': None, 'resol_max': None,
    'z_expr': None, 'logN_expr': None, 'b_expr': None, 'btur_expr': None, 'resol_expr': None}



# Default values for continuum adjustment parameters
adj_gauss_d = {
    'z': 0.0, 'ampl': 0.0, 'sigma': 0.01,
    'z_vary': False, 'ampl_vary': False, 'sigma_vary': False,
    'z_min': None, 'ampl_min': -0.05, 'sigma_min': 0.0,
    'z_max': None, 'ampl_max': 0.05, 'sigma_max': 1.0,
    'z_expr': None, 'ampl_expr': None, 'sigma_expr': None}


# Default values for line Voigt parameters
lines_voigt_d = {
    'z': 0.0, 'N': 1.e13, 'b': 10.0, 'btur': 0.0,
    'z_vary': True, 'N_vary': True, 'b_vary': True, 'btur_vary': False,
    'z_min': None, 'N_min': 1.e11, 'b_min': 5.0, 'btur_min': None,
    'z_max': None, 'N_max': 1.e17, 'b_max': 100.0, 'btur_max': None,
    'z_expr': None, 'N_expr': None, 'b_expr': None, 'btur_expr': None}

# Default values for PSF gaussian Parameters
psf_gauss_d = {
    'z': 0.0, 'resol': 35000,
    'z_vary': False, 'resol_vary': False,
    'z_min': None, 'resol_min': None,
    'z_max': None, 'resol_max': None,
    'z_expr': None, 'resol_expr': None}


pars_d = {'lines_voigt_d': lines_voigt_d,
          'psf_gauss_d': psf_gauss_d}

series_d = {'Ly': ['Ly_15', 'Ly_14', 'Ly_13', 'Ly_12', 'Ly_11', 'Ly_10',
                       'Ly_9', 'Ly_8', 'Ly_7', 'Ly_6', 'Ly_e', 'Ly_d', 'Ly_g',
                       'Ly_b', 'Ly_a'],
               'Ly_abg': ['Ly_g', 'Ly_b', 'Ly_a'],
               'Ly_ab': ['Ly_b', 'Ly_a'],
               'Ly_a': ['Ly_a'],
               'NV': ['NV_1238', 'NV_1242'],
               'SiII': ['SiII_1260', 'SiII_1304', 'SiII_1526'],
               'OI': ['OI_1302'],
               'SiIV': ['SiIV_1393', 'SiIV_1402'],
               'CIV': ['CIV_1548', 'CIV_1550'],
               'CIV_1548': ['CIV_1548'],
               'CIV_1550': ['CIV_1550'],
               'FeII': ['FeII_2344', 'FeII_2382', 'FeII_2586', 'FeII_2600'],
               'MgII': ['MgII_2796', 'MgII_2803'],
               'CaII': ['CaII_3934', 'CaII_3969'],
               'unknown': ['unknown']}

xem_d = {'Ly_a': 121.567 * au.nm,
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
             'SiIV_1393': 139.37602 * au.nm,
             'SiIV_1402': 140.27729 * au.nm,
             'CIV_1548': 154.8204 * au.nm,
             'CIV_1550': 155.0781 * au.nm,
             'FeII_2344': 234.42139 * au.nm,
             'FeII_2382': 238.27652 * au.nm,
             'FeII_2586': 258.665 * au.nm,
             'FeII_2600': 260.01729 * au.nm,
             'MgII_2796': 279.63543 * au.nm,
             'MgII_2803': 280.35315 * au.nm,
             'CaII_3934': 393.4775 * au.nm,
             'CaII_3969': 396.95901 * au.nm}

# Ionic oscillator strengths
fosc_d = {'Ly_a': 0.416,
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
          'SiII_1526': 0.133,
          'OI_1302': 0.048,
          'SiIV_1393': 0.513,
          'SiIV_1402': 0.254,
          'CIV_1548': 0.1899,
          'CIV_1550': 0.09475,
          'FeII_2344': 0.114,
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
gamma_d = {'Ly_a': 6.265e+08,
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
              'SiIV_1393': 8.8e+08,
              'SiIV_1402': 8.62e+08,
              'CIV_1548': 2.643e+08,
              'CIV_1550': 2.628e+08,
              'FeII_2344': 2.68e+08,
              'FeII_2382': 3.13e+08,
              'FeII_2586': 2.72e+08,
              'FeII_2600': 2.70e+08,
              'MgII_2796': 2.625e+08,
              'MgII_2803': 2.595e+08,
              'CaII_3934': 1.444e+08,
              'CaII_3969': 1.409e+08,
              'neb': 5e8,
              'unknown': 6.265e+08}
