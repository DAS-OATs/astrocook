"""Series lookup of the GUI: redshift cursor and identification viewer.

Characterization tests for phase 3 (``trans_parse`` moved out of
``astrocook.legacy``): ``SpectrumPlotWidget.get_cursor_line_positions_at_z``
and ``calculate_z_from_x`` (the series typed in the cursor box) and
``IdentificationViewerDialog._get_ref_transition``. The series are looked
up in ``STANDARD_MULTIPLETS``, then with the V1 series syntax, then as a
single transition. The expected values were checked against the V1
``trans_parse`` before the move, then updated for the single table of
atomic data (#568: series in the order of the table, F39).
"""
from types import SimpleNamespace

import astropy.units as au
import pytest

from astrocook.core.atomic_data import xem_d
from test_recipes_edit_flux import _session

Z_CURSOR = 2.5
X_MOUSE = 500.0  # nm

# series typed by the user -> transitions drawn by the cursor (those in
# xem_d), the first one being the reference of calculate_z_from_x
CURSOR_CASES = {
    "CIV": ["CIV_1548", "CIV_1550"],  # STANDARD_MULTIPLETS
    "Ly_a": ["Ly_a"],  # single transition
    "Ly-ab": ["Ly_a", "Ly_b"],  # V1 alias
    # V1 series: order of the table, Ly_a first (it was Ly_lim, F39)
    "Ly": ["Ly_a", "Ly_b", "Ly_g", "Ly_d", "Ly_e", "Ly_6", "Ly_7", "Ly_8",
           "Ly_9", "Ly_10", "Ly_11", "Ly_12", "Ly_13", "Ly_14", "Ly_15",
           "Ly_lim"],
    # V1 series: with OI_988a/b since the single table (#568)
    "OI": ["OI_936", "OI_948", "OI_971", "OI_976", "OI_988a", "OI_988b",
           "OI_1025", "OI_1039", "OI_1302", "OI_6302", "OI_6365"],
    "CIV_1548,SiIV": ["CIV_1548", "SiIV_1393", "SiIV_1402"],
    # Only whole strings are looked up in STANDARD_MULTIPLETS: here MgII is
    # the V1 series (six lines), while "MgII" alone gives the doublet
    "MgII;Ly_b": ["MgII_1025", "MgII_1026", "MgII_1239", "MgII_1240",
                  "MgII_2796", "MgII_2803", "Ly_b"],
    "H2": [],  # V1 series with no transition in the table
    "XYZ": [],  # unknown
}


@pytest.fixture
def plot_widget_class(qtbot):
    # Import only once qtbot has created the QApplication: the GUI selects the
    # QtAgg matplotlib backend at import time, which fails without one (F7).
    from astrocook.gui.pyside_plot import SpectrumPlotWidget
    return SpectrumPlotWidget


def _widget(series):
    """Stand-in for the plot widget: the helpers only read the series box."""
    box = SimpleNamespace(text=lambda: series)
    return SimpleNamespace(
        main_window=SimpleNamespace(cursor_series_input=box))


@pytest.mark.parametrize("series, transitions", CURSOR_CASES.items())
def test_cursor_line_positions(plot_widget_class, series, transitions):
    positions = plot_widget_class.get_cursor_line_positions_at_z(
        _widget(series), _session(), Z_CURSOR)
    expected = [(1 + Z_CURSOR) * xem_d[t].to_value(au.nm) for t in transitions]
    assert positions == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("series, transitions", CURSOR_CASES.items())
def test_z_from_x(plot_widget_class, series, transitions):
    z = plot_widget_class.calculate_z_from_x(_widget(series), X_MOUSE)
    if not transitions:
        assert z is None
    else:
        ref = xem_d[transitions[0]].to_value(au.nm)
        assert z == pytest.approx(X_MOUSE / ref - 1, rel=1e-12)


@pytest.mark.parametrize("series, ref", [
    ("CIV", "CIV_1548"),  # STANDARD_MULTIPLETS
    ("Ly_a", "Ly_a"),  # in xem_d
    ("Ly-ab", "Ly_a"),  # V1 alias
    ("Ly", "Ly_a"),  # V1 series (F39: it was Ly_lim)
    ("OI", "OI_936"),
    ("H2", None),  # V1 series with no transition in the table: the caller
                   # skips the series with a warning (as when it was
                   # H2_1110, not in xem_d)
    ("CIV_1548,SiIV", "CIV_1548"),
    ("SiII*_1264", "SiII*_1264"),
    ("XYZ", None),
])
def test_identification_viewer_ref_transition(qtbot, series, ref):
    from astrocook.gui.identification_viewer_dialog import (
        IdentificationViewerDialog)
    assert IdentificationViewerDialog._get_ref_transition(None, series) == ref
