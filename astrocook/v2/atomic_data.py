from astropy import units as au

# --- 1. Master Atomic Data Registry ---
# Key: Unique ID (Ion_Wavelength in Angstroms, usually)
# Data: {wave (Angstrom), f (oscillator strength), gamma (s^-1), mass (amu)}
# NOTE: Wavelengths here are in Angstroms as per standard convention, 
# but we will likely convert to nm for internal use in Astrocook V2.
ATOM_DATA = {
    # --- Hydrogen (HI Lyman Series) ---
    'Ly_a': {'wave': 1215.6700, 'f': 0.4164, 'gamma': 6.265e8, 'mass': 1.00794},
    'Ly_b': {'wave': 1025.7223, 'f': 0.07912, 'gamma': 1.897e8, 'mass': 1.00794},
    'Ly_g': {'wave': 972.5368, 'f': 0.02900, 'gamma': 8.127e7, 'mass': 1.00794},
    'Ly_d': {'wave': 949.7431, 'f': 0.01394, 'gamma': 4.204e7, 'mass': 1.00794},
    'Ly_e': {'wave': 937.8035, 'f': 0.007799, 'gamma': 2.450e7, 'mass': 1.00794},
    # (Add more Lyman series as needed)

    # --- Carbon ---
    'CII_1334': {'wave': 1334.5323, 'f': 0.1278, 'gamma': 2.880e8, 'mass': 12.0107},
    'CIV_1548': {'wave': 1548.2040, 'f': 0.1899, 'gamma': 2.643e8, 'mass': 12.0107},
    'CIV_1550': {'wave': 1550.7810, 'f': 0.09475, 'gamma': 2.628e8, 'mass': 12.0107},

    # --- Nitrogen ---
    'NV_1238': {'wave': 1238.821, 'f': 0.156, 'gamma': 3.391e8, 'mass': 14.0067},
    'NV_1242': {'wave': 1242.804, 'f': 0.0777, 'gamma': 3.356e8, 'mass': 14.0067},

    # --- Oxygen ---
    'OI_1302': {'wave': 1302.1685, 'f': 0.0480, 'gamma': 5.650e8, 'mass': 15.9994},
    'OVI_1031': {'wave': 1031.9261, 'f': 0.1325, 'gamma': 4.149e8, 'mass': 15.9994},
    'OVI_1037': {'wave': 1037.6167, 'f': 0.0658, 'gamma': 4.076e8, 'mass': 15.9994},

    # --- Magnesium ---
    'MgII_2796': {'wave': 2796.3543, 'f': 0.6155, 'gamma': 2.625e8, 'mass': 24.305},
    'MgII_2803': {'wave': 2803.5315, 'f': 0.3058, 'gamma': 2.595e8, 'mass': 24.305},

    # --- Silicon ---
    'SiII_1260': {'wave': 1260.4221, 'f': 1.18, 'gamma': 2.950e9, 'mass': 28.0855},
    'SiII_1526': {'wave': 1526.7066, 'f': 0.133, 'gamma': 1.130e9, 'mass': 28.0855},
    'SiIV_1393': {'wave': 1393.7602, 'f': 0.513, 'gamma': 8.800e8, 'mass': 28.0855},
    'SiIV_1402': {'wave': 1402.7729, 'f': 0.254, 'gamma': 8.620e8, 'mass': 28.0855},

    # --- Iron (Selected common ones) ---
    'FeII_2382': {'wave': 2382.765, 'f': 0.320, 'gamma': 3.13e8, 'mass': 55.845},
    'FeII_2600': {'wave': 2600.173, 'f': 0.239, 'gamma': 2.70e8, 'mass': 55.845},
}

# --- 2. Standard Multiplet Definitions ---
# Gives a clean way to refer to groups without ambiguity.
# Keys can be used in the GUI as "presets" for searching.
STANDARD_MULTIPLETS = {
    # Doublets
    'CIV': ['CIV_1548', 'CIV_1550'],
    'SiIV': ['SiIV_1393', 'SiIV_1402'],
    'NV': ['NV_1238', 'NV_1242'],
    'OVI': ['OVI_1031', 'OVI_1037'],
    'MgII': ['MgII_2796', 'MgII_2803'],
    
    # Series
    'Lyman': ['Ly_a', 'Ly_b', 'Ly_g', 'Ly_d', 'Ly_e']
}

# --- 3. V1 Backwards Compatibility Map ---
# Maps old awkward V1 tags to V2 lists of real lines.
# This solves your "CIV-1549" display issue: we know it maps to two real lines.
V1_TAG_MAP = {
    'CIV-1549': ['CIV_1548', 'CIV_1550'],
    'MgII-2799': ['MgII_2796', 'MgII_2803'],
    'SiIV-1398': ['SiIV_1393', 'SiIV_1402'],
    'NV-1240': ['NV_1238', 'NV_1242'],
    'OVI-1035': ['OVI_1031', 'OVI_1037'], # Approximate center
}