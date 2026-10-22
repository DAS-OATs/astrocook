"""The table of atomic data, ``trans_parse`` and the Ly-alpha velocity
conversion.

``trans_parse`` moved out of ``astrocook.legacy`` in phase 3, with the
series of V1. Since the single table (#568) its series are built from the
table: the same transitions as V1 within the table, but in the order of
the table, so that the first one, which the GUI takes as reference for the
redshift cursor, is the strongest line of H I and D I (F39).
``convert_axis_velocity`` must match V1 ``x_convert``, which it replaces in
the plot of the systems. The values of the table are checked against
standard line lists (#568, F38).
"""
import numpy as np
import astropy.units as au
import pytest
from astropy.table import Table

from astrocook.core.atomic_data import (
    ATOM_DATA, ATOMIC_DATA_PATH, _v1_series, trans_parse,
)
from astrocook.core.spectrum_operations import convert_axis_velocity
from reference.v1 import vars as vars_v1
from reference.v1.functions import trans_parse as trans_parse_v1
from reference.v1.functions import x_convert as x_convert_v1

_TABLE = Table.read(ATOMIC_DATA_PATH, format="ascii.ecsv")


#: Transitions of the table that V1 did not have (atom_par.dat named the
#: first two OIa_988 and OIb_988, series OIa and OIb)
NOT_IN_V1 = {"OI_988a", "OI_988b", "TiII_3067"}


def _v1_in_table(transitions):
    return [t for t in transitions if t in ATOM_DATA]


def test_series_are_those_of_v1_within_the_table():
    v1 = {k: _v1_in_table(v) for k, v in vars_v1.series_d.items()}
    series = _v1_series()
    # The V1 series with no transition in the table (H2, CO, Paschen,
    # Brackett, ...) are dropped: the GUI never drew them, since it keeps
    # only the transitions in xem_d
    assert set(series) == {k for k, v in v1.items() if v}
    for name, transitions in series.items():
        assert set(transitions) - NOT_IN_V1 == set(v1[name]), name


@pytest.mark.parametrize("series", [
    "CIV", "Ly", "Ly_a", "Ly-ab", "FeII", "OI", "CIV_1548,CIV_1550",
    "CIV, SiIV", "MgII;FeII", "Ly_b:Ly_a", "SiII*_1264",
])
def test_trans_parse_has_the_transitions_of_v1(series):
    new = [t for t in trans_parse(series) if t not in NOT_IN_V1]
    assert sorted(new) == sorted(_v1_in_table(trans_parse_v1(series)))


@pytest.mark.parametrize("series, expected", [
    ("CIV", ["CIV_1548", "CIV_1550"]),
    ("CIV_1548, SiIV", ["CIV_1548", "SiIV_1393", "SiIV_1402"]),
    ("Ly-ab", ["Ly_a", "Ly_b"]),
    ("Ly-abg", ["Ly_a", "Ly_b", "Ly_g"]),
    ("XYZ_123", ["XYZ_123"]),  # transitions are not checked
])
def test_trans_parse_examples(series, expected):
    assert trans_parse(series) == expected


@pytest.mark.parametrize("series, first, last, n", [
    ("Ly", "Ly_a", "Ly_lim", 16),
    ("D", "D_a", "D_12", 12),
    ("H", "H_a", "H_e", 5),
])
def test_trans_parse_hydrogen_order(series, first, last, n):
    # The redshift cursor of the GUI takes the first transition as
    # reference: Ly_a, not Ly_lim as in V1 (F39)
    trans = trans_parse(series)
    assert (trans[0], trans[-1], len(trans)) == (first, last, n)


@pytest.mark.parametrize("series", ["XYZ", "H2", "CIV-1549"])
def test_trans_parse_unknown_series(series):
    # H2 and CIV-1549 were V1 series with no transition in the table
    with pytest.raises(KeyError):
        trans_parse(series)


@pytest.mark.parametrize("zem", [0.0, 2.5])
@pytest.mark.parametrize("xunit", [au.nm, au.km / au.s, au.Angstrom])
def test_convert_axis_velocity_matches_v1_x_convert(zem, xunit):
    x = np.array([300.0, 425.4, 600.0]) * au.nm
    new = convert_axis_velocity(x, zem, xunit)
    old = x_convert_v1(x, zem=zem, xunit=xunit)
    assert new.unit == old.unit
    np.testing.assert_allclose(new.value, old.value, rtol=1e-12)


# --- Values against a standard line list (#568, F38) -------------------------

#: Morton (2003, ApJS 149, 205), vacuum wavelengths (Angstrom), oscillator
#: strengths and damping constants (s^-1), as distributed with linetools
#: (morton03_table2); gamma None where not checked. OIII 832 is below the
#: range of Morton (912 A): its value is the one of the V1 table, which the
#: V2 table had mistyped like the others.
REFERENCE = {
    "SiIV_1393": (1393.7602, 0.513, 8.80e8),
    "SiIV_1402": (1402.7729, 0.254, 8.62e8),
    "MgII_2796": (2796.3543, 0.6155, 2.625e8),
    "MgII_2803": (2803.5315, 0.3058, 2.595e8),
    "NV_1238": (1238.821, 0.156, 3.391e8),
    "OVI_1031": (1031.9261, 0.1325, 4.149e8),
    "AlII_1670": (1670.7886, 1.74, 1.39e9),
    # The four wavelengths mistyped in the conversion from nm to Angstrom
    # (91.39683 nm written 913.39683 A) and the damping constant off by ten
    "PIII_913": (913.9683, 0.203, 4.79e9),
    "PV_1117": (1117.9774, 0.472, 1.26e9),
    "FeI_2501": (2501.8858, 0.0493, 3.85e8),
    "ArI_1066": (1066.6598, 0.0675, 1.32e8),
    "OIII_832": (832.927, 0.107, None),
    # Also in the atom.dat of VPFIT (P. Noterdaeme), which has no line at
    # 3066.36 A: TiII_3066 is kept with f = 0
    "TiII_3067": (3067.2379, 0.0489, 2.50e8),
    # gamma was 1.17e8
    "CI_1560": (1560.3092, 0.0774, 1.27e8),
}

#: Newer values adopted instead of Morton (2003), as in the atom.dat of VPFIT
#: (Guido's decision, F38): name -> {column: (value, reference key)}. The
#: f of Zn II are 26% higher than in Morton (logN 0.10 dex lower), those of
#: Ni II up to 32% lower (1454: logN 0.17 dex higher), those of S II 5% lower
NEWER = {
    "CIV_1548": {"wave": (1548.2049, "PA04")},
    "CIV_1550": {"wave": (1550.7785, "PA04")},  # Morton: 1550.781
    "ZnII_2026": {"f": (0.630, "K15")},
    "ZnII_2062": {"f": (0.309, "K15")},
    "SII_1250": {"f": (0.00515, "K14")},
    "SII_1253": {"f": (0.01035, "K14")},
    "SII_1259": {"f": (0.015675, "K14")},
    "NiII_1317": {"f": (0.0596, "BB19")},
    "NiII_1370": {"f": (0.0616, "BB19")},
    "NiII_1454": {"f": (0.0220, "BB19")},
    "NiII_1709": {"f": (0.0318, "BB19")},
    "NiII_1741": {"f": (0.0432, "BB19")},
    "NiII_1751": {"f": (0.0280, "BB19")},
}

#: Forbidden emission lines, vacuum wavelengths (air: 3868.76, 3345.82)
VACUUM_EMISSION = {"NeIII_3869": 3869.86, "NeV_3346": 3346.79}


@pytest.mark.parametrize("name", sorted(REFERENCE))
def test_values_of_the_standard_line_list(name):
    wave, f, gamma = REFERENCE[name]
    atom = ATOM_DATA[name]
    assert atom["wave"] == pytest.approx(wave, abs=2e-3)  # 0.4 km/s at 1500 A
    assert atom["f"] == pytest.approx(f, rel=0.01)
    if gamma is not None:
        assert atom["gamma"] == pytest.approx(gamma, rel=0.02)


@pytest.mark.parametrize("name", sorted(VACUUM_EMISSION))
def test_emission_lines_are_in_vacuum(name):
    assert ATOM_DATA[name]["wave"] == pytest.approx(VACUUM_EMISSION[name], abs=0.01)


def test_no_invented_ti_ii_line():
    # The TiII_3066 entry had f and gamma of unknown origin for a line
    # that no reference lists
    assert ATOM_DATA["TiII_3066"]["f"] == 0.0


@pytest.mark.parametrize("name", sorted(NEWER))
def test_newer_values_with_their_source(name):
    row = _TABLE[_TABLE["name"] == name][0]
    for column, (value, ref) in NEWER[name].items():
        assert ATOM_DATA[name][column] == value
        assert row[column + "_ref"] == ref


# --- The table file -----------------------------------------------------------


def test_table_reads_the_same_with_astropy():
    # The package reads the ECSV file with the csv module: astropy must find
    # the same values, in the same order, with the declared units
    assert list(_TABLE["name"]) == list(ATOM_DATA)
    for column in ("wave", "f", "gamma", "mass"):
        assert list(_TABLE[column]) == [d[column] for d in ATOM_DATA.values()]
    assert _TABLE["wave"].unit == au.Angstrom
    assert _TABLE["gamma"].unit == 1 / au.s


def test_every_value_has_a_known_source():
    references = _TABLE.meta["references"]
    for column in ("wave", "f", "gamma"):
        refs = _TABLE[column + "_ref"]
        assert set(refs) <= set(references), column
        if column != "wave":
            # "-" (not applicable) exactly for the values that are 0
            assert list(refs == "-") == list(_TABLE[column] == 0), column
