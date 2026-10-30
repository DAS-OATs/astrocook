"""Cues of the flags of the fit in the System Inspector (#517, PR 3):
the parameter at a bound in amber, the ID of a component without
uncertainties in red, a warning sign on the ID, "Only flagged" and
"Select Flagged"."""
from types import SimpleNamespace

from PySide6.QtCore import Qt

from astrocook.core.structures import ComponentDataV2, SystemListDataV2
from astrocook.core.system_list import SystemListV2
import astrocook.gui.system_inspector as gui
from astrocook.gui.inspector_table import FLAG_BOUND_COLOR, FLAG_NAN_COLOR, flag_cue


def _gui():
    return gui


def _systs():
    comps = [ComponentDataV2(id=1, z=3.0, dz=1e-5, logN=13.5, dlogN=0.05, b=50.0,
                             db=1.0, series="Ly_a", flags=("b_at_upper_bound",)),
             ComponentDataV2(id=2, z=3.01, dz=None, logN=13.0, dlogN=None, b=20.0,
                             db=None, series="Ly_a", flags=("nan_errors",)),
             ComponentDataV2(id=3, z=3.02, dz=1e-5, logN=13.2, dlogN=0.05, b=25.0,
                             db=1.0, series="Ly_a")]
    return SystemListV2(SystemListDataV2(components=comps))


def _inspector(with_spectrum=False):
    systs = _systs()
    inspector = _gui().SystemInspector(None)
    if with_spectrum:
        # Selecting rows draws the velocity plot
        from astrocook.core.session import SessionV2
        from astrocook.mocks import MockComponent, MockTruth, generate
        spec = generate(MockTruth(
            components=[MockComponent("Ly_a", 3.01, 13.0, 20.0)], lam_min_nm=484.0,
            lam_max_nm=490.0, dv_kms=2.5, resol=45000.0, snr=30.0, z_em=3.5, seed=1))
        inspector.current_session = SessionV2("s", None, spec=spec, systs=systs)
    else:
        inspector.current_session = SimpleNamespace(systs=systs)
    inspector.table_model.update_data(systs.components, {})
    return inspector


def _cell(inspector, row, column, role):
    model = inspector.table_model
    return model.data(model.index(row, _gui().COLUMNS.index(column)), role)


def test_flag_cue(qapp):
    flagged, nan, clean = _systs().components
    assert flag_cue(flagged, "b")[0] == FLAG_BOUND_COLOR
    assert flag_cue(flagged, "logN") == (None, [])
    assert flag_cue(nan)[0] == FLAG_NAN_COLOR
    assert flag_cue(clean) == (None, [])


def test_table_cues(qapp):
    inspector = _inspector()
    # The b at the bound in amber, not the other parameters
    assert _cell(inspector, 0, "b", Qt.BackgroundRole).name() == FLAG_BOUND_COLOR
    assert _cell(inspector, 0, "logN", Qt.BackgroundRole) is None
    assert "upper bound" in _cell(inspector, 0, "b", Qt.ToolTipRole)
    # No uncertainties: the ID in red
    assert _cell(inspector, 1, "ID", Qt.BackgroundRole).name() == FLAG_NAN_COLOR
    assert "could not be computed" in _cell(inspector, 1, "ID", Qt.ToolTipRole)
    # A warning sign on the ID of the flagged components only
    assert _cell(inspector, 0, "ID", Qt.DisplayRole).startswith("⚠")
    assert _cell(inspector, 2, "ID", Qt.DisplayRole) == "3"


def test_only_flagged_and_select_flagged(qapp):
    inspector = _inspector(with_spectrum=True)
    assert inspector.select_flagged() == 2
    selected = {inspector.proxy_model.mapToSource(i).row()
                for i in inspector.table_view.selectionModel().selectedRows()}
    assert selected == {0, 1}
    inspector.flagged_cb.setChecked(True)
    assert inspector.proxy_model.rowCount() == 2
    inspector.flagged_cb.setChecked(False)
    assert inspector.proxy_model.rowCount() == 3
