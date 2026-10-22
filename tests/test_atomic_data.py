"""``trans_parse`` and the Ly-alpha velocity conversion moved out of
``astrocook.legacy`` (phase 3).

``trans_parse`` must give the same transitions, in the same order, as V1
(the GUI takes the first one as reference for the redshift cursor), and
``convert_axis_velocity`` must match V1 ``x_convert``, which it replaces
in the plot of the systems. The explicit cases document the behaviour
without V1, for when the V1 code leaves the package.
"""
import numpy as np
import astropy.units as au
import pytest

from astrocook.core.atomic_data import ATOM_DATA, _v1_series, trans_parse
from astrocook.core.spectrum_operations import convert_axis_velocity
from reference.v1 import vars as vars_v1
from reference.v1.functions import trans_parse as trans_parse_v1
from reference.v1.functions import x_convert as x_convert_v1


def test_series_match_v1():
    assert _v1_series() == {k: list(v) for k, v in vars_v1.series_d.items()}


@pytest.mark.parametrize("series", [
    "CIV", "Ly", "Ly_a", "Ly-ab", "FeII", "H2", "OI", "CIV-1549",
    "CIV_1548,CIV_1550", "CIV, SiIV", "MgII;FeII", "Ly_b:Ly_a", "SiII*_1264",
])
def test_trans_parse_matches_v1(series):
    assert trans_parse(series) == trans_parse_v1(series)


@pytest.mark.parametrize("series, expected", [
    ("CIV", ["CIV_1548", "CIV_1550"]),
    ("CIV_1548, SiIV", ["CIV_1548", "SiIV_1393", "SiIV_1402"]),
    ("Ly-abg", ["Ly_g", "Ly_b", "Ly_a"]),
    ("XYZ_123", ["XYZ_123"]),  # transitions are not checked
])
def test_trans_parse_examples(series, expected):
    assert trans_parse(series) == expected


def test_trans_parse_lyman_order():
    # Order of atom_par.dat: the redshift cursor of the GUI uses the first
    # transition known to xem_d as reference
    ly = trans_parse("Ly")
    assert ly[0] == "Ly_lim" and ly[-1] == "Ly_a" and len(ly) == 16


def test_trans_parse_unknown_series():
    with pytest.raises(KeyError):
        trans_parse("XYZ")


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
    "CIV_1548": (1548.204, 0.1899, 2.643e8),
    "CIV_1550": (1550.781, 0.09475, 2.628e8),
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
}


@pytest.mark.parametrize("name", sorted(REFERENCE))
def test_values_of_the_standard_line_list(name):
    wave, f, gamma = REFERENCE[name]
    atom = ATOM_DATA[name]
    assert atom["wave"] == pytest.approx(wave, abs=2e-3)  # 0.4 km/s at 1500 A
    assert atom["f"] == pytest.approx(f, rel=0.01)
    if gamma is not None:
        assert atom["gamma"] == pytest.approx(gamma, rel=0.02)
