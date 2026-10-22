import csv
import functools
from pathlib import Path
from typing import Dict, List

import numpy as np
from astropy import constants as const
from astropy import units as au


# --- 1. Master Atomic Data Registry ---
#: The table of atomic data (#568): one row per transition, with the source
#: of each value (columns ``wave_ref``, ``f_ref``, ``gamma_ref``, keys of the
#: references listed in its header). It is an ECSV file (astropy), read here
#: with the csv module to keep the import light.
ATOMIC_DATA_PATH = Path(__file__).resolve().parent.parent / 'data' / 'atomic_data.ecsv'


def _read_atomic_data(path: Path) -> Dict[str, Dict[str, float]]:
    """Rows of the table, keyed by transition, in the order of the file."""
    lines = [line for line in path.read_text().splitlines()
             if line.strip() and not line.startswith('#')]
    reader = csv.reader(lines, delimiter=' ', quotechar='"',
                        skipinitialspace=True)
    header = next(reader)
    data: Dict[str, Dict[str, float]] = {}
    for values in reader:
        row = dict(zip(header, values))
        if row['name'] in data:
            raise ValueError(f"{path.name}: duplicate transition {row['name']}")
        data[row['name']] = {key: float(row[key])
                             for key in ('wave', 'f', 'gamma', 'mass')}
    return data


# Key: transition (ion and rounded wavelength)
# Data: {wave (Angstrom, vacuum), f (oscillator strength), gamma (s^-1),
#        mass (amu)}
ATOM_DATA = _read_atomic_data(ATOMIC_DATA_PATH)

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


# --- 3. Series names of V1 ---
# The series of the V1 syntax (``series_d`` of V1): the transitions of the
# table grouped by the part of the name before the underscore, in the order
# of the table, so that the first one is the reference of the series (e.g.
# 'Ly' goes from 'Ly_a' to 'Ly_lim'; in V1, which followed atom_par.dat, it
# went from 'Ly_lim' to 'Ly_a', F39), plus three Lyman aliases.
@functools.cache
def _v1_series() -> Dict[str, List[str]]:
    """Transitions of each V1 series, built once from the table."""
    series: Dict[str, List[str]] = {}
    for name in ATOM_DATA:
        series.setdefault(name.split('_')[0], []).append(name)
    series['Ly-a'] = ['Ly_a']
    series['Ly-ab'] = ['Ly_a', 'Ly_b']
    series['Ly-abg'] = ['Ly_a', 'Ly_b', 'Ly_g']
    return series

def trans_parse(series: str) -> List[str]:
    """
    Expand a series string into the list of its transitions (V1 syntax).

    The string is split at commas, semicolons and colons. Items with an
    underscore are taken as transitions (e.g. ``'CIV_1548'``), the others
    as V1 series names, expanded with all their transitions (e.g. ``'CIV'``
    gives ``['CIV_1548', 'CIV_1550']``).

    Parameters
    ----------
    series : str
        Series string, e.g. ``'CIV'``, ``'Ly_a'`` or ``'CIV_1548,SiIV'``.

    Returns
    -------
    list of str
        The transitions, in the order of the string and, within a series,
        of the table of atomic data (the first is the reference).

    Raises
    ------
    KeyError
        If an item without an underscore is not a V1 series name.
    """
    trans = []
    for s in series.replace(';', ',').replace(':', ',').split(','):
        if '_' in s:
            trans.append(s.strip())
        else:
            trans.extend(_v1_series()[s.strip()])
    return trans
