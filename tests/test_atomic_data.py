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

from astrocook.core.atomic_data import _v1_series, trans_parse
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
