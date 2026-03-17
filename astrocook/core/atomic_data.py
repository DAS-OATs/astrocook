from astropy import constants as const
from astropy import units as au
import numpy as np
from typing import Dict


# --- 1. Master Atomic Data Registry ---
# Key: Unique ID from atom_par.dat
# Data: {wave (Angstrom), f (oscillator strength), gamma (s^-1), mass (amu)}
ATOM_DATA = {
    # --- Hydrogen (Lyman Series) ---
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
    'Ly_11': {'wave': 918.1293, 'f': 0.0009219, 'gamma': 2.432e6, 'mass': 1.00794},
    'Ly_12': {'wave': 917.1805, 'f': 0.0007231, 'gamma': 1.911e6, 'mass': 1.00794},
    'Ly_13': {'wave': 916.4291, 'f': 0.0005777, 'gamma': 1.529e6, 'mass': 1.00794},
    'Ly_14': {'wave': 915.8238, 'f': 0.0004689, 'gamma': 1.243e6, 'mass': 1.00794},
    'Ly_15': {'wave': 915.3289, 'f': 0.0003858, 'gamma': 1.024e6, 'mass': 1.00794},
    'Ly_lim': {'wave': 911.75, 'f': 0.0, 'gamma': 0.0, 'mass': 1.00794},

    # --- Hydrogen (Balmer Series) ---
    'H_a': {'wave': 6562.83, 'f': 0.64108, 'gamma': 6.265e8, 'mass': 1.00794},
    'H_b': {'wave': 4861.34, 'f': 0.11938, 'gamma': 1.897e8, 'mass': 1.00794},
    'H_g': {'wave': 4340.47, 'f': 0.044694, 'gamma': 8.126e7, 'mass': 1.00794},
    'H_d': {'wave': 4101.74, 'f': 0.022105, 'gamma': 4.204e7, 'mass': 1.00794},
    'H_e': {'wave': 3970.08, 'f': 0.012711, 'gamma': 2.450e7, 'mass': 1.00794},

    # --- Deuterium (Lyman Series) ---
    'D_a': {'wave': 1215.3394, 'f': 0.416, 'gamma': 6.270e8, 'mass': 2.0079},
    'D_b': {'wave': 1025.4432, 'f': 0.0791, 'gamma': 1.897e8, 'mass': 2.0079},
    'D_g': {'wave': 972.2722, 'f': 0.02899, 'gamma': 8.127e7, 'mass': 2.0079},
    'D_d': {'wave': 949.4846, 'f': 0.01394, 'gamma': 4.203e7, 'mass': 2.0079},
    'D_e': {'wave': 937.5483, 'f': 0.0078, 'gamma': 2.450e7, 'mass': 2.0079},
    'D_6': {'wave': 930.4950, 'f': 0.00482, 'gamma': 1.237e7, 'mass': 2.0079},
    'D_7': {'wave': 925.9737, 'f': 0.00318, 'gamma': 8.261e6, 'mass': 2.0079},
    'D_8': {'wave': 922.8992, 'f': 0.00222, 'gamma': 5.789e6, 'mass': 2.0079},
    'D_9': {'wave': 920.7125, 'f': 0.0016, 'gamma': 4.214e6, 'mass': 2.0079},
    'D_10': {'wave': 919.1013, 'f': 0.0012, 'gamma': 3.162e6, 'mass': 2.0079},
    'D_11': {'wave': 917.8796, 'f': 0.00092, 'gamma': 2.434e6, 'mass': 2.0079},
    'D_12': {'wave': 916.9310, 'f': 0.00072, 'gamma': 1.913e6, 'mass': 2.0079},

    # --- Helium ---
    'HeI_3889': {'wave': 3889.0, 'f': 0.0, 'gamma': 0.0, 'mass': 4.0026},
    'HeI_5876': {'wave': 5876.0, 'f': 0.0, 'gamma': 0.0, 'mass': 4.0026},
    'HeII_1640': {'wave': 1640.4, 'f': 0.0, 'gamma': 0.0, 'mass': 4.0026},
    'HeII_4686': {'wave': 4686.0, 'f': 0.0, 'gamma': 0.0, 'mass': 4.0026},

    # --- Lithium ---
    'LiI_6709': {'wave': 6709.6276, 'f': 0.498, 'gamma': 3.689e7, 'mass': 6.941},

    # --- Carbon ---
    'CI_945': {'wave': 945.191, 'f': 0.152, 'gamma': 3.790e8, 'mass': 12.011},
    'CI_1277': {'wave': 1277.2452, 'f': 0.0853, 'gamma': 2.320e8, 'mass': 12.011},
    'CI_1560': {'wave': 1560.3092, 'f': 0.0774, 'gamma': 1.170e8, 'mass': 12.011},
    'CI_1656': {'wave': 1656.9284, 'f': 0.149, 'gamma': 3.600e8, 'mass': 12.011},
    'CII_903': {'wave': 903.6234, 'f': 0.168, 'gamma': 6.862e8, 'mass': 12.011},
    'CII_904': {'wave': 903.9616, 'f': 0.336, 'gamma': 2.743e9, 'mass': 12.011},
    'CII_1036': {'wave': 1036.3367, 'f': 0.118, 'gamma': 2.200e9, 'mass': 12.011},
    'CII_1334': {'wave': 1334.5323, 'f': 0.128, 'gamma': 2.880e8, 'mass': 12.011},
    'CII*_1335': {'wave': 1335.7077, 'f': 0.1149, 'gamma': 2.880e8, 'mass': 12.011}, 
    'CII_2326': {'wave': 2326.0, 'f': 0.0, 'gamma': 0.0, 'mass': 12.011},
    'CIII_977': {'wave': 977.0201, 'f': 0.757, 'gamma': 1.760e9, 'mass': 12.011}, 
    'CIII_1908': {'wave': 1908.734, 'f': 0.0, 'gamma': 0.0, 'mass': 12.011},
    'CIV_1548': {'wave': 1548.2040, 'f': 0.1899, 'gamma': 2.643e8, 'mass': 12.011},
    'CIV_1550': {'wave': 1550.7810, 'f': 0.09475, 'gamma': 2.628e8, 'mass': 12.011},

    # --- Nitrogen ---
    'NI_963': {'wave': 963.9903, 'f': 0.0124, 'gamma': 8.550e7, 'mass': 14.0067},
    'NI_964': {'wave': 964.6256, 'f': 0.0079, 'gamma': 8.310e7, 'mass': 14.0067},
    'NI_965': {'wave': 965.0413, 'f': 0.00386, 'gamma': 8.160e7, 'mass': 14.0067},
    'NI_1134a': {'wave': 1134.1653, 'f': 0.0146, 'gamma': 1.510e8, 'mass': 14.0067},
    'NI_1134b': {'wave': 1134.4149, 'f': 0.0287, 'gamma': 1.490e8, 'mass': 14.0067},
    'NI_1134c': {'wave': 1134.9803, 'f': 0.0416, 'gamma': 1.440e8, 'mass': 14.0067},
    'NI_1199': {'wave': 1199.5496, 'f': 0.132, 'gamma': 4.070e8, 'mass': 14.0067},
    'NI_1200': {'wave': 1200.2233, 'f': 0.0869, 'gamma': 4.020e8, 'mass': 14.0067},
    'NI_1201': {'wave': 1200.7098, 'f': 0.0432, 'gamma': 4.000e8, 'mass': 14.0067}, 
    'NI_6529': {'wave': 6529.03, 'f': 0.0, 'gamma': 0.0, 'mass': 14.0067},
    'NII_915': {'wave': 915.6131, 'f': 0.159, 'gamma': 1.270e9, 'mass': 14.0067},
    'NII_1083': {'wave': 1083.9937, 'f': 0.111, 'gamma': 3.740e8, 'mass': 14.0067},
    'NII_6549': {'wave': 6549.86, 'f': 0.0, 'gamma': 0.0, 'mass': 14.0067},
    'NII_6585': {'wave': 6585.27, 'f': 0.0, 'gamma': 0.0, 'mass': 14.0067},
    'NIII_989': {'wave': 989.799, 'f': 0.123, 'gamma': 5.000e8, 'mass': 14.0067},
    'NV_1238': {'wave': 1238.821, 'f': 0.156, 'gamma': 3.391e8, 'mass': 14.0067},
    'NV_1242': {'wave': 1242.804, 'f': 0.0777, 'gamma': 3.356e8, 'mass': 14.0067},

    # --- Oxygen ---
    'OI_936': {'wave': 936.6295, 'f': 0.00365, 'gamma': 0.0, 'mass': 15.9994},
    'OI_948': {'wave': 948.6855, 'f': 0.00631, 'gamma': 0.0, 'mass': 15.9994},
    'OI_971': {'wave': 971.7382, 'f': 0.0116, 'gamma': 5.850e7, 'mass': 15.9994},
    'OI_976': {'wave': 976.4481, 'f': 0.0033, 'gamma': 3.860e7, 'mass': 15.9994},
    'OI_988a': {'wave': 988.7734, 'f': 0.0465, 'gamma': 2.260e8, 'mass': 15.9994},
    'OI_988b': {'wave': 988.6549, 'f': 0.0083, 'gamma': 5.660e7, 'mass': 15.9994},
    'OI_1025': {'wave': 1025.7616, 'f': 0.0163, 'gamma': 1.020e8, 'mass': 15.9994},
    'OI_1039': {'wave': 1039.2304, 'f': 0.00907, 'gamma': 1.870e8, 'mass': 15.9994},
    'OI_1302': {'wave': 1302.1685, 'f': 0.0480, 'gamma': 5.650e8, 'mass': 15.9994},
    'OI_6302': {'wave': 6302.046, 'f': 0.0, 'gamma': 0.0, 'mass': 15.9994},
    'OI_6365': {'wave': 6365.536, 'f': 0.0, 'gamma': 0.0, 'mass': 15.9994},
    'OII_832': {'wave': 832.7572, 'f': 0.0444, 'gamma': 2.135e8, 'mass': 15.9994},
    'OII_833': {'wave': 833.3294, 'f': 0.0886, 'gamma': 8.510e8, 'mass': 15.9994},
    'OII_834': {'wave': 834.4655, 'f': 0.132, 'gamma': 8.430e8, 'mass': 15.9994},
    'OII_3727': {'wave': 3727.092, 'f': 0.0, 'gamma': 0.0, 'mass': 15.9994},
    'OII_3729': {'wave': 3729.875, 'f': 0.0, 'gamma': 0.0, 'mass': 15.9994},
    'OIII_832': {'wave': 832.2927, 'f': 0.107, 'gamma': 3.429e8, 'mass': 15.9994},
    'OIII_1665': {'wave': 1665.85, 'f': 0.0, 'gamma': 0.0, 'mass': 15.9994},
    'OIII_4364': {'wave': 4364.436, 'f': 0.0, 'gamma': 0.0, 'mass': 15.9994},
    'OIII_4932': {'wave': 4932.603, 'f': 0.0, 'gamma': 0.0, 'mass': 15.9994},
    'OIII_4960': {'wave': 4960.295, 'f': 0.0, 'gamma': 0.0, 'mass': 15.9994},
    'OIII_5008': {'wave': 5008.24, 'f': 0.0, 'gamma': 0.0, 'mass': 15.9994},
    'OVI_1031': {'wave': 1031.9261, 'f': 0.1325, 'gamma': 4.149e8, 'mass': 15.9994},
    'OVI_1037': {'wave': 1037.6167, 'f': 0.0658, 'gamma': 4.076e8, 'mass': 15.9994},

    # --- Neon ---
    'NeIII_3869': {'wave': 3868.760, 'f': 3.9e-10, 'gamma': 1.74e-1, 'mass': 19.992},
    'NeV_3346': {'wave': 3345.821, 'f': 2.1e-10, 'gamma': 7.60e-2, 'mass': 19.992},

    # --- Sodium ---
    'NaI_5891': {'wave': 5891.5833, 'f': 0.6408, 'gamma': 61570000.0, 'mass': 22.9898},
    'NaI_5897': {'wave': 5897.5581, 'f': 0.3201, 'gamma': 61390000.0, 'mass': 22.9898},

    # --- Magnesium ---
    'MgI_1827': {'wave': 1827.9351, 'f': 0.0242, 'gamma': 1.610e7, 'mass': 24.305},
    'MgI_2026': {'wave': 2026.4768, 'f': 0.113, 'gamma': 7.250e7, 'mass': 24.305},
    'MgI_2852': {'wave': 2852.9631, 'f': 1.83, 'gamma': 5.000e8, 'mass': 24.305},
    'MgII_1025': {'wave': 1025.9681, 'f': 0.000743, 'gamma': 2.350e6, 'mass': 24.305},
    'MgII_1026': {'wave': 1026.1134, 'f': 0.000392, 'gamma': 2.480e6, 'mass': 24.305},
    'MgII_1239': {'wave': 1239.9253, 'f': 0.000632, 'gamma': 1.370e6, 'mass': 24.305},
    'MgII_1240': {'wave': 1240.3947, 'f': 0.000356, 'gamma': 1.540e6, 'mass': 24.305},
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
    'SiII*_1264': {'wave': 1264.7377, 'f': 1.05, 'gamma': 2.930e9, 'mass': 28.0855},
    'SiII_1304': {'wave': 1304.3702, 'f': 0.0863, 'gamma': 1.010e9, 'mass': 28.0855},
    'SiII_1526': {'wave': 1526.7066, 'f': 0.133, 'gamma': 1.130e9, 'mass': 28.0855},
    'SiII_1808': {'wave': 1808.0129, 'f': 0.00208, 'gamma': 2.380e6, 'mass': 28.0855},
    'SiIII_1206': {'wave': 1206.500, 'f': 1.63, 'gamma': 2.480e9, 'mass': 28.0855},
    'SiIV_1393': {'wave': 1393.7602, 'f': 0.513, 'gamma': 8.800e8, 'mass': 28.0855},
    'SiIV_1402': {'wave': 1402.7729, 'f': 0.254, 'gamma': 8.620e8, 'mass': 28.0855},

    # --- Phosphorous ---
    'PII_961': {'wave': 961.0412, 'f': 0.349, 'gamma': 5.090e9, 'mass': 30.9736},
    'PII_963': {'wave': 963.8005, 'f': 1.46, 'gamma': 4.140e9, 'mass': 30.9736},
    'PII_1152': {'wave': 1152.818, 'f': 0.245, 'gamma': 1.230e9, 'mass': 30.9736},
    'PIII_913': {'wave': 913.39683, 'f': 0.203, 'gamma': 4.790e9, 'mass': 30.9736},
    'PIII_917': {'wave': 917.1178, 'f': 0.404, 'gamma': 4.780e9, 'mass': 30.9736},
    'PV_1117': {'wave': 1117.79774, 'f': 0.472, 'gamma': 1.260e9, 'mass': 30.9736},
    'PV_1128': {'wave': 1128.0078, 'f': 0.233, 'gamma': 1.220e9, 'mass': 30.9736},

    # --- Sulfur ---
    'SII_1250': {'wave': 1250.578, 'f': 0.00543, 'gamma': 4.630e7, 'mass': 32.066},
    'SII_1253': {'wave': 1253.805, 'f': 0.0109, 'gamma': 4.620e7, 'mass': 32.066},
    'SII_1259': {'wave': 1259.518, 'f': 0.0166, 'gamma': 4.650e7, 'mass': 32.066},
    'SII_4072': {'wave': 4072.3, 'f': 0.0, 'gamma': 0.0, 'mass': 32.066},
    'SII_6718': {'wave': 6718.29, 'f': 0.0, 'gamma': 0.0, 'mass': 32.066},
    'SII_6732': {'wave': 6732.67, 'f': 0.0, 'gamma': 0.0, 'mass': 32.066},
    'SIII_1012': {'wave': 1012.495, 'f': 0.0438, 'gamma': 2.810e8, 'mass': 32.066},
    'SIII_1190': {'wave': 1190.203, 'f': 0.0237, 'gamma': 6.650e7, 'mass': 32.066},
    'SIII_9069': {'wave': 9069.0, 'f': 0.0, 'gamma': 0.0, 'mass': 32.066},
    'SIII_9532': {'wave': 9532.0, 'f': 0.0, 'gamma': 0.0, 'mass': 32.066},
    'SIV_1062': {'wave': 1062.664, 'f': 0.0494, 'gamma': 1.690e8, 'mass': 32.066},
    'SVI_933': {'wave': 933.378, 'f': 0.437, 'gamma': 1.670e9, 'mass': 32.066},
    'SVI_944': {'wave': 944.523, 'f': 0.215, 'gamma': 1.610e9, 'mass': 32.066},

    # --- Argon ---
    'ArI_1048': {'wave': 1048.2199, 'f': 0.263, 'gamma': 3.180e7, 'mass': 39.948},
    'ArI_1066': {'wave': 1066.6598, 'f': 0.0675, 'gamma': 1.320e7, 'mass': 39.948},

    # --- Potassium ---
    'KI_7701': {'wave': 7701.0835, 'f': 0.3327, 'gamma': 3.742e7, 'mass': 39.0983},
    
    # --- Calcium ---
    'CaI_4227': {'wave': 4227.918, 'f': 1.77, 'gamma': 2.200e8, 'mass': 40.078},
    'CaII_3934': {'wave': 3934.775, 'f': 0.6267, 'gamma': 144400000.0, 'mass': 40.078},
    'CaII_3969': {'wave': 3969.5901, 'f': 0.3116, 'gamma': 140900000.0, 'mass': 40.078},
    'CaII_8500': {'wave': 8500.36, 'f': 0.0, 'gamma': 0.0, 'mass': 40.078},
    'CaII_8544': {'wave': 8544.44, 'f': 0.0, 'gamma': 0.0, 'mass': 40.078},
    'CaII_8664': {'wave': 8664.52, 'f': 0.0, 'gamma': 0.0, 'mass': 40.078},

    # --- Titanium ---
    'TiII_3066':  {'wave': 3066.362, 'f': 0.107, 'gamma': 2.10e8,  'mass': 47.867},
    'TiII_3384': {'wave': 3384.7304, 'f': 0.358, 'gamma': 1.750e8, 'mass': 47.867},
    'TiIII_1298': {'wave': 1298.697, 'f': 0.0964, 'gamma': 6.350e8, 'mass': 47.867},

    # --- Chromium ---
    'CrII_2056': {'wave': 2056.2569, 'f': 0.103, 'gamma': 4.070e8, 'mass': 51.9961},
    'CrII_2062': {'wave': 2062.2361, 'f': 0.0759, 'gamma': 4.060e8, 'mass': 51.9961},
    'CrII_2066': {'wave': 2066.164, 'f': 0.0512, 'gamma': 4.170e8, 'mass': 51.9961},

    # --- Manganese ---
    'MnII_1197': {'wave': 1197.184, 'f': 0.217, 'gamma': 7.850e8, 'mass': 54.9309},
    'MnII_1199': {'wave': 1199.391, 'f': 0.169, 'gamma': 7.850e8, 'mass': 54.9309},
    'MnII_2576': {'wave': 2576.877, 'f': 0.361, 'gamma': 2.820e8, 'mass': 54.9308},
    'MnII_2594': {'wave': 2594.499, 'f': 0.28, 'gamma': 2.780e8, 'mass': 54.9308},
    'MnII_2606': {'wave': 2606.462, 'f': 0.198, 'gamma': 2.720e8, 'mass': 54.9308},

    # --- Iron ---
    'FeI_1851': {'wave': 1851.6902, 'f': 0.02222, 'gamma': 5.558e7, 'mass': 55.847},
    'FeI_2167': {'wave': 2167.4534, 'f': 0.15, 'gamma': 2.736e8, 'mass': 55.847},
    'FeI_2298': {'wave': 2298.8769, 'f': 0.0245, 'gamma': 1.350e8, 'mass': 55.847},
    'FeI_2463': {'wave': 2463.3922, 'f': 0.0532, 'gamma': 5.000e8, 'mass': 55.847},
    'FeI_2484': {'wave': 2484.0209, 'f': 0.544, 'gamma': 5.000e8, 'mass': 55.847},
    'FeI_2501': {'wave': 2501.18858, 'f': 0.0493, 'gamma': 3.850e8, 'mass': 55.847},
    'FeI_2523': {'wave': 2523.6083, 'f': 0.203, 'gamma': 3.850e8, 'mass': 55.847},
    'FeI_2719': {'wave': 2719.8329, 'f': 0.122, 'gamma': 1.750e8, 'mass': 55.847},
    'FeI_2967': {'wave': 2967.7646, 'f': 0.0438, 'gamma': 1.270e8, 'mass': 55.847},
    'FeI_2984': {'wave': 2984.4402, 'f': 0.029, 'gamma': 1.750e8, 'mass': 55.847},
    'FeI_3021': {'wave': 3021.5187, 'f': 0.104, 'gamma': 1.690e8, 'mass': 55.847},
    'FeI_3441': {'wave': 3441.5918, 'f': 0.0236, 'gamma': 1.710e7, 'mass': 55.847},
    'FeI_3720': {'wave': 3720.9928, 'f': 0.0411, 'gamma': 1.630e7, 'mass': 55.847},
    'FeI_3861': {'wave': 3861.0058, 'f': 0.0217, 'gamma': 1.280e7, 'mass': 55.847},
    'FeII_1081': {'wave': 1081.8748, 'f': 0.0126, 'gamma': 5.980e7, 'mass': 55.847},
    'FeII_1096': {'wave': 1096.8769, 'f': 0.0327, 'gamma': 2.260e8, 'mass': 55.847},
    'FeII_1122': {'wave': 1121.9748, 'f': 0.029, 'gamma': 1.920e8, 'mass': 55.847},
    'FeII_1133': {'wave': 1133.6654, 'f': 0.00472, 'gamma': 3.060e7, 'mass': 55.847},
    'FeII_1142': {'wave': 1142.3656, 'f': 0.00401, 'gamma': 2.560e7, 'mass': 55.847},
    'FeII_1143': {'wave': 1143.226, 'f': 0.0192, 'gamma': 9.810e7, 'mass': 55.847},
    'FeII_1145': {'wave': 1144.9379, 'f': 0.083, 'gamma': 3.520e8, 'mass': 55.847},
    'FeII_1260': {'wave': 1260.533, 'f': 0.024, 'gamma': 1.260e8, 'mass': 55.847},
    'FeII_1608': {'wave': 1608.4511, 'f': 0.0577, 'gamma': 2.740e8, 'mass': 55.845},
    'FeII_1611': {'wave': 1611.2005, 'f': 0.00138, 'gamma': 2.860e8, 'mass': 55.845},
    'FeII_2249': {'wave': 2249.8768, 'f': 0.00182, 'gamma': 3.310e8, 'mass': 55.847}, 
    'FeII_2260': {'wave': 2260.7805, 'f': 0.00244, 'gamma': 2.580e8, 'mass': 55.847}, 
    'FeII_2344': {'wave': 2344.214, 'f': 0.114, 'gamma': 2.68e8, 'mass': 55.845},
    'FeII_2367': {'wave': 2367.5905, 'f': 2.16e-05, 'gamma': 3.070e8, 'mass': 55.847},
    'FeII_2374': {'wave': 2374.461, 'f': 0.0313, 'gamma': 3.09e8, 'mass': 55.845},
    'FeII_2382': {'wave': 2382.765, 'f': 0.320, 'gamma': 3.13e8, 'mass': 55.845},
    'FeII_2586': {'wave': 2586.650, 'f': 0.0691, 'gamma': 2.72e8, 'mass': 55.845},
    'FeII_2600': {'wave': 2600.173, 'f': 0.239, 'gamma': 2.70e8, 'mass': 55.845},
    'FeIII_859': {'wave': 859.723, 'f': 0.115, 'gamma': 8.491e8, 'mass': 55.847},
    'FeIII_1122': {'wave': 1122.524, 'f': 0.0544, 'gamma': 3.700e8, 'mass': 55.847},

    # --- Cobalt ---
    'CoII_1466': {'wave': 1466.211, 'f': 0.031, 'gamma': 7.870e7, 'mass': 58.9332},
    'CoII_1574': {'wave': 1574.5508, 'f': 0.025, 'gamma': 6.730e7, 'mass': 58.9332},
    'CoII_1941': {'wave': 1941.2852, 'f': 0.034, 'gamma': 4.350e8, 'mass': 58.9332},
    'CoII_2012': {'wave': 2012.1664, 'f': 0.0368, 'gamma': 3.450e8, 'mass': 58.9332},

    # --- Nickel ---
    'NiII_1317': {'wave': 1317.217, 'f': 0.07786, 'gamma': 4.2038e8, 'mass': 58.6934},
    'NiII_1370': {'wave': 1370.132, 'f': 0.0769, 'gamma': 4.100e8, 'mass': 58.6934},
    'NiII_1454': {'wave': 1454.842, 'f': 0.0323, 'gamma': 1.0179e8, 'mass': 58.6934},
    'NiII_1709': {'wave': 1709.6042, 'f': 0.0324, 'gamma': 4.350e8, 'mass': 58.6934},
    'NiII_1741': {'wave': 1741.5531, 'f': 0.0427, 'gamma': 5.000e8, 'mass': 58.6934},
    'NiII_1751': {'wave': 1751.9157, 'f': 0.0277, 'gamma': 3.700e8, 'mass': 58.6934},

    # --- Copper ---
    'CuII_1358': {'wave': 1358.773, 'f': 0.263, 'gamma': 7.200e8, 'mass': 63.546},
    
    # --- Zinc ---
    'ZnI_1589': {'wave': 1589.561, 'f': 0.1219, 'gamma': 1.073e8, 'mass': 65.39},
    'ZnI_2139': {'wave': 2139.2477, 'f': 1.47, 'gamma': 7.140e8, 'mass': 65.39},
    'ZnII_2026': {'wave': 2026.137, 'f': 0.501, 'gamma': 4.070e8, 'mass': 65.39}, 
    'ZnII_2062': {'wave': 2062.6604, 'f': 0.246, 'gamma': 3.860e8, 'mass': 65.39}, 

    # --- Gallium ---
    'GaII_1414': {'wave': 1414.402, 'f': 1.77, 'gamma': 1.970e9, 'mass': 60.723},

    # --- Tin ---
    'SnII_1400': {'wave': 1400.4, 'f': 0.7141, 'gamma': 0.0, 'mass': 118.71},
}

# --- 2. Standard Multiplet Definitions ---
# Gives a clean way to refer to groups without ambiguity.
METAL_MULTIPLETS = {
    # Doublets (strongest first for primary ID)
    'CIV': ['CIV_1548', 'CIV_1550'],
    'SiIV': ['SiIV_1393', 'SiIV_1402'],
    'NV': ['NV_1238', 'NV_1242'],
    'OVI': ['OVI_1031', 'OVI_1037'],
    'MgII': ['MgII_2796', 'MgII_2803'],
    'AlIII': ['AlIII_1854', 'AlIII_1862'],
    'CII': ['CII_1334', 'CII*_1335'],
    'ZnII': ['ZnII_2026', 'ZnII_2062'], 
    'NaI': ['NaI_5891', 'NaI_5897'],
    'CaII': ['CaII_3934', 'CaII_3969'],
    
    # Complex multiplets
    'NI': ['NI_963', 'NI_964', 'NI_965', 'NI_1134a', 'NI_1134b', 'NI_1134c', 'NI_1199', 'NI_1200', 'NI_1201'],
    'SiII': ['SiII_1260', 'SiII_1526', 'SiII_1193', 'SiII_1190', 'SiII*_1264', 'SiII_1304', 'SiII_1808'],
    'SII': ['SII_1250', 'SII_1253', 'SII_1259'],
    'CrII': ['CrII_2056', 'CrII_2062', 'CrII_2066'],
    'MnII': ['MnII_2576', 'MnII_2594', 'MnII_2606'],
    'FeII': ['FeII_1608', 'FeII_1611', 'FeII_2249', 'FeII_2260', 'FeII_2344', 'FeII_2367', 'FeII_2374', 'FeII_2382', 'FeII_2586', 'FeII_2600'],
    'NiII': ['NiII_1317', 'NiII_1370', 'NiII_1454', 'NiII_1709', 'NiII_1741', 'NiII_1751'],
}

HYDROGEN_MULTIPLETS = {
    # Series
    'Lyman': ['Ly_a', 'Ly_b', 'Ly_g', 'Ly_d', 'Ly_e', 'Ly_6', 'Ly_7', 'Ly_8'],
    'Ly_a': ['Ly_a'], 
    'Ly_ab': ['Ly_a', 'Ly_b'], 
}

STANDARD_MULTIPLETS = {**METAL_MULTIPLETS, **HYDROGEN_MULTIPLETS}

# Built from ATOM_DATA to be used by V2 functions.
# This replaces the dependency on legacy.vars.xem_d
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