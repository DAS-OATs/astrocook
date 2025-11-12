from astropy import constants as const
from astropy import units as au
import numpy as np
from typing import Dict


# --- 1. Master Atomic Data Registry ---
# Key: Unique ID from atom_par.dat
# Data: {wave (Angstrom), f (oscillator strength), gamma (s^-1), mass (amu)}
ATOM_DATA = {
    # --- Hydrogen (HI Lyman Series) ---
    'Ly_a': {'wave': 1215.6700, 'f': 0.4164, 'gamma': 6.265e8, 'mass': 1.00794},
    'Ly_b': {'wave': 1025.7222, 'f': 0.07912, 'gamma': 1.897e8, 'mass': 1.00794},
    'Ly_g': {'wave': 972.5367, 'f': 0.02900, 'gamma': 8.127e7, 'mass': 1.00794},
    'Ly_d': {'wave': 949.7430, 'f': 0.01394, 'gamma': 4.204e7, 'mass': 1.00794},
    'Ly_e': {'wave': 937.8034, 'f': 0.007799, 'gamma': 2.450e7, 'mass': 1.00794},
    'Ly_6': {'wave': 930.7482, 'f': 0.00481, 'gamma': 1.236e7, 'mass': 1.00794},
    'Ly_7': {'wave': 926.2256, 'f': 0.003185, 'gamma': 8.255e6, 'mass': 1.00794},
    'Ly_8': {'wave': 923.1503, 'f': 0.002216, 'gamma': 5.785e6, 'mass': 1.00794},
    'Ly_9': {'wave': 920.9630, 'f': 0.001605, 'gamma': 4.210e6, 'mass': 1.00794},
    'Ly_10': {'wave': 919.3513, 'f': 0.001201, 'gamma': 3.160e6, 'mass': 1.00794},
    'Ly_lim': {'wave': 911.75, 'f': 0.0, 'gamma': 0.0, 'mass': 1.00794},

    # --- Carbon ---
    'CII_1036': {'wave': 1036.3367, 'f': 0.118, 'gamma': 2.200e9, 'mass': 12.011},
    'CII_1334': {'wave': 1334.5323, 'f': 0.128, 'gamma': 2.880e8, 'mass': 12.011},
    'CII*_1335': {'wave': 1335.7077, 'f': 0.1149, 'gamma': 2.880e8, 'mass': 12.011}, 
    'CIII_977': {'wave': 977.0201, 'f': 0.757, 'gamma': 1.760e9, 'mass': 12.011}, 
    'CIV_1548': {'wave': 1548.2040, 'f': 0.1899, 'gamma': 2.643e8, 'mass': 12.011},
    'CIV_1550': {'wave': 1550.7810, 'f': 0.09475, 'gamma': 2.628e8, 'mass': 12.011},

    # --- Nitrogen ---
    'NI_1134a': {'wave': 1134.1653, 'f': 0.0146, 'gamma': 1.510e8, 'mass': 14.0067},
    'NI_1134b': {'wave': 1134.4149, 'f': 0.0287, 'gamma': 1.490e8, 'mass': 14.0067},
    'NI_1134c': {'wave': 1134.9803, 'f': 0.0416, 'gamma': 1.440e8, 'mass': 14.0067},
    'NI_1199': {'wave': 1199.5496, 'f': 0.132, 'gamma': 4.070e8, 'mass': 14.0067},
    'NI_1200': {'wave': 1200.2233, 'f': 0.0869, 'gamma': 4.020e8, 'mass': 14.0067},
    'NI_1201': {'wave': 1200.7098, 'f': 0.0432, 'gamma': 4.000e8, 'mass': 14.0067}, 
    'NII_1083': {'wave': 1083.9937, 'f': 0.111, 'gamma': 3.740e8, 'mass': 14.0067}, # <<< *** ADDED ***
    'NIII_989': {'wave': 989.799, 'f': 0.123, 'gamma': 5.000e8, 'mass': 14.0067}, # <<< *** ADDED ***
    'NV_1238': {'wave': 1238.821, 'f': 0.156, 'gamma': 3.391e8, 'mass': 14.0067},
    'NV_1242': {'wave': 1242.804, 'f': 0.0777, 'gamma': 3.356e8, 'mass': 14.0067},

    # --- Oxygen ---
    'OI_988a': {'wave': 988.7734, 'f': 0.0465, 'gamma': 2.260e8, 'mass': 15.9994},
    'OI_1039': {'wave': 1039.2304, 'f': 0.00907, 'gamma': 1.870e8, 'mass': 15.9994},
    'OI_1302': {'wave': 1302.1685, 'f': 0.0480, 'gamma': 5.650e8, 'mass': 15.9994},
    'OVI_1031': {'wave': 1031.9261, 'f': 0.1325, 'gamma': 4.149e8, 'mass': 15.9994},
    'OVI_1037': {'wave': 1037.6167, 'f': 0.0658, 'gamma': 4.076e8, 'mass': 15.9994},

    # --- Magnesium ---
    'MgI_1827': {'wave': 1827.9351, 'f': 0.0242, 'gamma': 1.610e7, 'mass': 24.305}, # <<< *** ADDED ***
    'MgI_2026': {'wave': 2026.4768, 'f': 0.113, 'gamma': 7.250e7, 'mass': 24.305}, # <<< *** ADDED ***
    'MgI_2852': {'wave': 2852.9631, 'f': 1.83, 'gamma': 5.000e8, 'mass': 24.305}, # <<< *** ADDED ***
    'MgII_1239': {'wave': 1239.9253, 'f': 0.000632, 'gamma': 1.370e6, 'mass': 24.305}, # <<< *** ADDED ***
    'MgII_1240': {'wave': 1240.3947, 'f': 0.000356, 'gamma': 1.540e6, 'mass': 24.305}, # <<< *** ADDED ***
    'MgII_2796': {'wave': 2796.3543, 'f': 0.6155, 'gamma': 2.625e8, 'mass': 24.305},
    'MgII_2803': {'wave': 2803.5315, 'f': 0.3058, 'gamma': 2.595e8, 'mass': 24.305},

    # --- Aluminum ---
    'AlII_1670': {'wave': 1670.7874, 'f': 1.74, 'gamma': 1.390e9, 'mass': 26.9815},
    'AlIII_1854': {'wave': 1854.7184, 'f': 0.559, 'gamma': 5.420e8, 'mass': 26.9815},
    'AlIII_1862': {'wave': 1862.7910, 'f': 0.278, 'gamma': 5.340e8, 'mass': 26.9815},

    # --- Silicon ---
    'SiII_1190': {'wave': 1190.4158, 'f': 0.292, 'gamma': 4.080e9, 'mass': 28.0855},
    'SiII_1193': {'wave': 1193.2897, 'f': 0.582, 'gamma': 4.070e9, 'mass': 28.0855},
    'SiII_1260': {'wave': 1260.4221, 'f': 1.18, 'gamma': 2.950e9, 'mass': 28.0855},
    'SiII*_1264': {'wave': 1264.7377, 'f': 1.05, 'gamma': 2.930e9, 'mass': 28.0855}, # <<< *** ADDED ***
    'SiII_1304': {'wave': 1304.3702, 'f': 0.0863, 'gamma': 1.010e9, 'mass': 28.0855},
    'SiII_1526': {'wave': 1526.7066, 'f': 0.133, 'gamma': 1.130e9, 'mass': 28.0855},
    'SiII_1808': {'wave': 1808.0129, 'f': 0.00208, 'gamma': 2.380e6, 'mass': 28.0855},
    'SiIII_1206': {'wave': 1206.500, 'f': 1.63, 'gamma': 2.480e9, 'mass': 28.0855},
    'SiIV_1393': {'wave': 1393.7602, 'f': 0.513, 'gamma': 8.800e8, 'mass': 28.0855},
    'SiIV_1402': {'wave': 1402.7729, 'f': 0.254, 'gamma': 8.620e8, 'mass': 28.0855},

    # --- Phosphorous ---
    'PII_961': {'wave': 961.0412, 'f': 0.349, 'gamma': 5.090e9, 'mass': 30.9736}, # <<< *** ADDED ***
    'PII_963': {'wave': 963.8005, 'f': 1.46, 'gamma': 4.140e9, 'mass': 30.9736}, # <<< *** ADDED ***
    'PII_1152': {'wave': 1152.818, 'f': 0.245, 'gamma': 1.230e9, 'mass': 30.9736}, # <<< *** ADDED ***

    # --- Sulfur ---
    'SII_1250': {'wave': 1250.578, 'f': 0.00543, 'gamma': 4.630e7, 'mass': 32.066}, # <<< *** ADDED ***
    'SII_1253': {'wave': 1253.805, 'f': 0.0109, 'gamma': 4.620e7, 'mass': 32.066}, # <<< *** ADDED ***
    'SII_1259': {'wave': 1259.518, 'f': 0.0166, 'gamma': 4.650e7, 'mass': 32.066}, # <<< *** ADDED ***
    'SIII_1012': {'wave': 1012.495, 'f': 0.0438, 'gamma': 2.810e8, 'mass': 32.066}, # <<< *** ADDED ***
    'SIII_1190': {'wave': 1190.203, 'f': 0.0237, 'gamma': 6.650e7, 'mass': 32.066}, # <<< *** ADDED ***
    'SIV_1062': {'wave': 1062.664, 'f': 0.0494, 'gamma': 1.690e8, 'mass': 32.066}, # <<< *** ADDED ***
    'SVI_933': {'wave': 933.378, 'f': 0.437, 'gamma': 1.670e9, 'mass': 32.066}, # <<< *** ADDED ***
    'SVI_944': {'wave': 944.523, 'f': 0.215, 'gamma': 1.610e9, 'mass': 32.066}, # <<< *** ADDED ***

    # --- Argon ---
    'ArI_1048': {'wave': 1048.2199, 'f': 0.263, 'gamma': 3.180e7, 'mass': 39.948}, # <<< *** ADDED ***
    'ArI_1066': {'wave': 1066.6598, 'f': 0.0675, 'gamma': 1.320e7, 'mass': 39.948}, # <<< *** ADDED ***

    # --- Chromium ---
    'CrII_2056': {'wave': 205.62569, 'f': 0.103, 'gamma': 4.070e8, 'mass': 51.9961}, # <<< *** ADDED ***
    'CrII_2062': {'wave': 206.22361, 'f': 0.0759, 'gamma': 4.060e8, 'mass': 51.9961}, # <<< *** ADDED ***
    'CrII_2066': {'wave': 206.6164, 'f': 0.0512, 'gamma': 4.170e8, 'mass': 51.9961}, # <<< *** ADDED ***

    # --- Manganese ---
    'MnII_1197': {'wave': 1197.184, 'f': 0.217, 'gamma': 7.850e8, 'mass': 54.9309}, # <<< *** ADDED ***
    'MnII_1199': {'wave': 1199.391, 'f': 0.169, 'gamma': 7.850e8, 'mass': 54.9309}, # <<< *** ADDED ***
    'MnII_2576': {'wave': 257.6877, 'f': 0.361, 'gamma': 2.820e8, 'mass': 54.9308}, # <<< *** ADDED ***
    'MnII_2594': {'wave': 259.4499, 'f': 0.28, 'gamma': 2.780e8, 'mass': 54.9308}, # <<< *** ADDED ***
    'MnII_2606': {'wave': 260.6462, 'f': 0.198, 'gamma': 2.720e8, 'mass': 54.9308}, # <<< *** ADDED ***

    # --- Iron ---
    'FeII_1081': {'wave': 1081.8748, 'f': 0.0126, 'gamma': 5.980e7, 'mass': 55.847}, # <<< *** ADDED ***
    'FeII_1096': {'wave': 1096.8769, 'f': 0.0327, 'gamma': 2.260e8, 'mass': 55.847}, # <<< *** ADDED ***
    'FeII_1122': {'wave': 1121.9748, 'f': 0.029, 'gamma': 1.920e8, 'mass': 55.847}, # <<< *** ADDED ***
    'FeII_1133': {'wave': 1133.6654, 'f': 0.00472, 'gamma': 3.060e7, 'mass': 55.847}, # <<< *** ADDED ***
    'FeII_1142': {'wave': 1142.3656, 'f': 0.00401, 'gamma': 2.560e7, 'mass': 55.847}, # <<< *** ADDED ***
    'FeII_1143': {'wave': 1143.226, 'f': 0.0192, 'gamma': 9.810e7, 'mass': 55.847}, # <<< *** ADDED ***
    'FeII_1145': {'wave': 1144.9379, 'f': 0.083, 'gamma': 3.520e8, 'mass': 55.847}, # <<< *** ADDED ***
    'FeII_1260': {'wave': 1260.533, 'f': 0.024, 'gamma': 1.260e8, 'mass': 55.847}, # <<< *** ADDED ***
    'FeII_1608': {'wave': 1608.4511, 'f': 0.0577, 'gamma': 2.740e8, 'mass': 55.845},
    'FeII_1611': {'wave': 1611.2005, 'f': 0.00138, 'gamma': 2.860e8, 'mass': 55.845}, # <<< *** ADDED ***
    'FeII_2249': {'wave': 2249.8768, 'f': 0.00182, 'gamma': 3.310e8, 'mass': 55.847}, 
    'FeII_2260': {'wave': 2260.7805, 'f': 0.00244, 'gamma': 2.580e8, 'mass': 55.847}, 
    'FeII_2344': {'wave': 2344.214, 'f': 0.114, 'gamma': 2.68e8, 'mass': 55.845},
    'FeII_2367': {'wave': 2367.5905, 'f': 2.16e-05, 'gamma': 3.070e8, 'mass': 55.847}, # <<< *** ADDED ***
    'FeII_2374': {'wave': 2374.461, 'f': 0.0313, 'gamma': 3.09e8, 'mass': 55.845},
    'FeII_2382': {'wave': 2382.765, 'f': 0.320, 'gamma': 3.13e8, 'mass': 55.845},
    'FeII_2586': {'wave': 2586.650, 'f': 0.0691, 'gamma': 2.72e8, 'mass': 55.845},
    'FeII_2600': {'wave': 2600.173, 'f': 0.239, 'gamma': 2.70e8, 'mass': 55.845},
    'FeIII_1122': {'wave': 1122.524, 'f': 0.0544, 'gamma': 3.700e8, 'mass': 55.847}, # <<< *** ADDED ***

    # --- Nickel ---
    'NiII_1317': {'wave': 1317.217, 'f': 0.07786, 'gamma': 4.2038e8, 'mass': 58.6934}, # <<< *** ADDED ***
    'NiII_1370': {'wave': 1370.132, 'f': 0.0769, 'gamma': 4.100e8, 'mass': 58.6934}, # <<< *** ADDED ***
    'NiII_1454': {'wave': 1454.842, 'f': 0.0323, 'gamma': 1.0179e8, 'mass': 58.6934}, # <<< *** ADDED ***
    'NiII_1709': {'wave': 1709.6042, 'f': 0.0324, 'gamma': 4.350e8, 'mass': 58.6934}, # <<< *** ADDED ***
    'NiII_1741': {'wave': 1741.5531, 'f': 0.0427, 'gamma': 5.000e8, 'mass': 58.6934}, # <<< *** ADDED ***
    'NiII_1751': {'wave': 1751.9157, 'f': 0.0277, 'gamma': 3.700e8, 'mass': 58.6934}, # <<< *** ADDED ***

    # --- Zinc ---
    'ZnII_2026': {'wave': 202.6137, 'f': 0.501, 'gamma': 4.070e8, 'mass': 65.39}, 
    'ZnII_2062': {'wave': 206.26604, 'f': 0.246, 'gamma': 3.860e8, 'mass': 65.39}, 
}

# --- 2. Standard Multiplet Definitions ---
# Gives a clean way to refer to groups without ambiguity.
STANDARD_MULTIPLETS = {
    # Doublets (strongest first for primary ID)
    'CIV': ['CIV_1548', 'CIV_1550'],
    'SiIV': ['SiIV_1393', 'SiIV_1402'],
    'NV': ['NV_1238', 'NV_1242'],
    'OVI': ['OVI_1031', 'OVI_1037'],
    'MgII': ['MgII_2796', 'MgII_2803'],
    'AlIII': ['AlIII_1854', 'AlIII_1862'],
    'CII': ['CII_1334', 'CII*_1335'],
    'ZnII': ['ZnII_2026', 'ZnII_2062'], # <<< *** ADDED ***
    
    # Complex multiplets
    'SiII': ['SiII_1260', 'SiII_1526', 'SiII_1193', 'SiII_1190', 'SiII_1304', 'SiII_1808'],
    'FeII': ['FeII_2382', 'FeII_2600', 'FeII_2344', 'FeII_2586', 'FeII_1608', 'FeII_2374'],
    'NI': ['NI_1200', 'NI_1199', 'NI_1201'],
    'CrII': ['CrII_2056', 'CrII_2062', 'CrII_2066'], # <<< *** ADDED ***

    # Series
    'Lyman': ['Ly_a', 'Ly_b', 'Ly_g', 'Ly_d', 'Ly_e', 'Ly_6', 'Ly_7', 'Ly_8'],
    'Ly_a': ['Ly_a'], 
    'Ly_ab': ['Ly_a', 'Ly_b'], 
}

# Built from ATOM_DATA to be used by V2 functions.
# This replaces the dependency on v1.vars.xem_d
xem_d = {
    key: data['wave'] * au.Angstrom
    for key, data in ATOM_DATA.items()
}

def is_hydrogen_line(transition_name: str) -> bool:
    """
    Checks if a transition is a Hydrogen line based on atomic mass.
    """
    if transition_name not in ATOM_DATA:
        return False
    # Check if mass is < 2.0 (i.e., H or D)
    return ATOM_DATA[transition_name]['mass'] < 2.0

def get_multiplet_velocity_lags(multiplet_name: str) -> Dict[str, float]:
    """
    Calculates the velocity lags (in km/s) of lines in a multiplet 
    relative to the primary (strongest) line.
    """
    if multiplet_name not in STANDARD_MULTIPLETS:
        raise ValueError(f"Unknown multiplet '{multiplet_name}'")
        
    lines = STANDARD_MULTIPLETS[multiplet_name]
    primary_wave = ATOM_DATA[lines[0]]['wave']
    c_kms = const.c.to(au.km/au.s).value
    
    lags = {}
    for line_id in lines:
        wave = ATOM_DATA[line_id]['wave']
        # dv = c * ln(lambda / lambda_primary)
        dv = c_kms * np.log(wave / primary_wave)
        lags[line_id] = dv
        
    return lags