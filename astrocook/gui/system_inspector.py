import importlib.util
import json
import logging
import numpy as np
import astropy.units as au
from PySide6.QtCore import (QEvent, Qt, QLocale, QItemSelection, QItemSelectionModel, QSettings,
                            QTimer)
from PySide6.QtGui import QAction, QDoubleValidator
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QTableView, 
    QHeaderView, QAbstractItemView, QMenu, QInputDialog, QDialog, QDialogButtonBox,
    QMessageBox, QLabel, QLineEdit, QCheckBox, QPushButton, QFrame
)
from typing import List

from ..core.structures import ComponentDataV2
from ..core.atomic_data import STANDARD_MULTIPLETS, ATOM_DATA
from ..core.resolution import resolution_at
from ..fitting import bayesian_fitter

from .inspector_table import (COL_MAP, COLUMNS, LinkColorDelegate, SystemSortFilterProxyModel,
                              SystemTableModel)
from .velocity_plot import VelocityPlotWidget
from .theme import (FOCUS_BORDER, SELECTION, SELECTION_TEXT, TOOLTIP_BACKGROUND, TOOLTIP_BORDER,
                    TOOLTIP_TEXT)


#: The columns of the table shown or hidden by the user, in the settings
HIDDEN_COLUMNS_KEY = "inspector_hidden_columns"


def _corner_installed() -> bool:
    """Whether corner (optional dependency, for the corner plot) is installed."""
    return importlib.util.find_spec("corner") is not None


class SystemInspector(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.current_session = None
        self.active_transitions: List[str] = [] 
        self._manual_trans_edit = False 
        self._manual_base_text = ""  # Stores the persistent manual text
        self.setWindowTitle("System Inspector")
        self.resize(1200, 800) 

        # [CHANGE] Apply Rounded Tooltip Style
        self.setStyleSheet(f"""
            QTableView {{
                selection-background-color: {SELECTION};
                selection-color: {SELECTION_TEXT};
            }}
            QToolTip {{
                border: 1px solid {TOOLTIP_BORDER};
                border-radius: 5px;
                background-color: {TOOLTIP_BACKGROUND};
                color: {TOOLTIP_TEXT};
                padding: 4px;
                background-clip: border-box;
            }}
        """)

        self._setup_ui()
        
    def _setup_ui(self):
        # [Task 3] Zero margins for the main window to match Main Window style
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(0)
        
        # [Task 1] Reduce splitter size
        splitter = QSplitter(Qt.Horizontal)
        splitter.setHandleWidth(9)
        # Optional: Light gray line to make the thin split visible
        #splitter.setStyleSheet("QSplitter::handle { background-color: #D3D3D3; }")
        
        # --- Left Side: Table ---
        self.table_view = QTableView()
        self.table_model = SystemTableModel()
        
        # --- Apply the custom text color delegate ---
        self._link_delegate = LinkColorDelegate(self.table_view)
        self.table_view.setItemDelegate(self._link_delegate)

        self.proxy_model = SystemSortFilterProxyModel()
        self.proxy_model.setSourceModel(self.table_model)
        
        self.table_view.setModel(self.proxy_model)
        self.table_view.setSortingEnabled(True)

        # Default Sort: Redshift (z) Ascending
        # We dynamically find the index of "z" to be safe against column reordering
        if "z" in COLUMNS:
            self.table_view.sortByColumn(COLUMNS.index("z"), Qt.AscendingOrder)

        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_view.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.table_view.verticalHeader().setVisible(False)
        self.table_view.setAlternatingRowColors(True)
        self.table_view.setContextMenuPolicy(Qt.CustomContextMenu)
        # The columns fit their contents, refitted once when the rows shown
        # change. Not ResizeToContents: Qt then measured every row (all the
        # roles of every cell) at each layout pass, ~3 s on 764 components
        self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self._column_fit_pending = False
        for signal in (self.proxy_model.modelReset, self.proxy_model.layoutChanged,
                       self.proxy_model.rowsInserted, self.proxy_model.rowsRemoved):
            signal.connect(self._schedule_column_fit)
        self.table_view.selectionModel().selectionChanged.connect(self._on_selection_changed)
        self.table_view.customContextMenuRequested.connect(self._on_context_menu)
        # The columns shown: a menu on the header; btur, almost always a
        # frozen 0, is hidden unless a component uses it (visual review of
        # phase 5-bis), or as the user chose (saved in the settings)
        header = self.table_view.horizontalHeader()
        header.setContextMenuPolicy(Qt.CustomContextMenu)
        header.customContextMenuRequested.connect(self._on_header_menu)
        self.table_model.data_changed_request.connect(self._forward)
        
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(0)
        
        # Only the group of the selection: the components fitted with it
        # (get_connected_group). The old label, "Group View (Selected
        # Component Only)", did not say so (findings register F85)
        self.group_cb = QCheckBox("Only group")
        self.group_cb.setToolTip("Show only the group of the selected components: those "
                                 "linked to them, blended with them or close in velocity, "
                                 "which are fitted together")
        self.group_cb.toggled.connect(self._on_group_toggled)
        
        # Flags of the fit (#517): show only, or select, the components
        # with a parameter at a bound or without uncertainties
        self.flagged_cb = QCheckBox("Only flagged")
        self.flagged_cb.setToolTip("Show only the components with warnings of the fit "
                                   "(\u26a0: a parameter at a bound, or no uncertainties)")
        self.flagged_cb.toggled.connect(self.proxy_model.set_flagged_only)
        self.select_flagged_btn = QPushButton("Select Flagged")
        self.select_flagged_btn.setToolTip("Select the components with warnings of the fit, "
                                           "e.g. to refit them (table menu > Refit Selected)")
        self.select_flagged_btn.clicked.connect(self.select_flagged)

        l_toolbar = QHBoxLayout()
        l_toolbar.setContentsMargins(5, 5, 15, 15) # Exact match to Right Toolbar
        # Its own spacing: a nested layout takes that of its parent, 0 here,
        # and the boxes touched (F85)
        l_toolbar.setSpacing(12)
        l_toolbar.addWidget(self.group_cb)
        l_toolbar.addWidget(self.flagged_cb)
        l_toolbar.addWidget(self.select_flagged_btn)
        l_toolbar.addStretch() # Push checkbox to the left
        
        left_layout.addLayout(l_toolbar)
        left_layout.addWidget(self.table_view)
        
        splitter.addWidget(left_widget)
        
        # --- Right Side: Plot & Controls ---
        right_widget = QWidget() # Renamed locally for clarity
        r_layout = QVBoxLayout(right_widget)
        r_layout.setContentsMargins(0, 0, 0, 0)
        r_layout.setSpacing(0)
        
        # Right Toolbar
        c_layout = QHBoxLayout()
        c_layout.setContentsMargins(5, 5, 5, 5) # Exact match to Left Toolbar
        c_layout.setSpacing(5)

        # Three groups (visual review of phase 5-bis): what the plot shows,
        # the starting values and fit window of a new component, the
        # residuals; on one row, so that the plot starts where the table
        # does: short labels, explained by the tooltips
        trans_tip = "The transitions shown in the velocity plot (e.g. CIV, SiIV)"
        c_layout.addWidget(self._field_label("Trans.:", trans_tip))
        self.trans_in = QLineEdit()
        self.trans_in.setPlaceholderText("CIV, SiIV")
        self.trans_in.setToolTip(trans_tip)
        self.trans_in.returnPressed.connect(self._apply)
        self.trans_in.textEdited.connect(self._on_manual_edit)
        c_layout.addWidget(self.trans_in, 1) 
        
        z_tip = "The redshift at the centre of the velocity plot (Δv = 0)"
        c_layout.addWidget(self._field_label("z:", z_tip))
        self.z_in = QLineEdit()
        self.z_in.setToolTip(z_tip)
        self.z_in.setFixedWidth(70)
        validator = QDoubleValidator()
        validator.setLocale(QLocale.C) 
        self.z_in.setValidator(validator)
        self.z_in.returnPressed.connect(self._apply_z)
        c_layout.addWidget(self.z_in)

        dv_tip = "The velocity range of the plot (km/s), from and to"
        c_layout.addWidget(self._field_label("Δv:", dv_tip))
        self.vmin_in = QLineEdit("-300")
        self.vmin_in.setToolTip(dv_tip)
        self.vmin_in.setFixedWidth(50)
        self.vmin_in.setValidator(validator)
        self.vmin_in.returnPressed.connect(self._apply_limits)
        c_layout.addWidget(self.vmin_in)
        self.vmax_in = QLineEdit("300")
        self.vmax_in.setToolTip(dv_tip)
        self.vmax_in.setFixedWidth(50)
        self.vmax_in.setValidator(validator)
        self.vmax_in.returnPressed.connect(self._apply_limits)
        c_layout.addWidget(self.vmax_in)
        
        c_layout.addWidget(self._separator())
        logn_tip = ("The starting logN of a component added from the menu of the\n"
                    "velocity plot (right click); its profile follows the cursor")
        c_layout.addWidget(self._field_label("logN:", logn_tip))
        self.logn_in = QLineEdit("13.5")
        self.logn_in.setToolTip(logn_tip)
        self.logn_in.setFixedWidth(40)
        self.logn_in.setValidator(validator)
        self.logn_in.textChanged.connect(self._update_cursor)
        c_layout.addWidget(self.logn_in)
        
        b_tip = ("The starting b (km/s) of a component added from the menu of\n"
                 "the velocity plot (right click); its profile follows the cursor")
        c_layout.addWidget(self._field_label("b:", b_tip))
        self.b_in = QLineEdit("10.0")
        self.b_in.setToolTip(b_tip)
        self.b_in.setFixedWidth(40)
        self.b_in.setValidator(validator)
        self.b_in.textChanged.connect(self._update_cursor)
        c_layout.addWidget(self.b_in)

        # Velocity window of the manual fits (#518, F4), next to the
        # starting values of a new component: empty keeps the default of
        # each fit, so that nothing changes unless the user sets a value
        fit_window_tip = (
            "How far (±, km/s) the redshift of each component may move from\n"
            "its starting value in the fits run from the System Inspector:\n"
            "adding a component here, Refit Selected and Refit All (table menu).\n"
            "Empty: the default of each fit, 20 km/s when adding a component,\n"
            "100 km/s when refitting.")
        fit_window_label = QLabel("Fit ±:")
        fit_window_label.setToolTip(fit_window_tip)
        c_layout.addWidget(fit_window_label)
        self.fit_window_in = QLineEdit()
        self.fit_window_in.setPlaceholderText("def.")
        self.fit_window_in.setFixedWidth(40)  # as logN and b
        positive = QDoubleValidator(0.0, 10000.0, 3)
        positive.setLocale(QLocale.C)
        positive.setNotation(QDoubleValidator.StandardNotation)
        self.fit_window_in.setValidator(positive)
        self.fit_window_in.setToolTip(fit_window_tip)
        c_layout.addWidget(self.fit_window_in)
        
        c_layout.addWidget(self._separator())
        self.resid_cb = QCheckBox("Resid.")
        self.resid_cb.setToolTip("Show the residuals of the model below each panel")
        self.resid_cb.toggled.connect(self._apply) 
        c_layout.addWidget(self.resid_cb)
        
        self._apply_field_styles()
        
        r_layout.addLayout(c_layout)
        
        self.vel_plot = VelocityPlotWidget(self)
        r_layout.addWidget(self.vel_plot, 1) 
        
        splitter.addWidget(right_widget)
        splitter.setStretchFactor(0, 6); splitter.setStretchFactor(1, 4)
        layout.addWidget(splitter)

        # Add helper to set style based on state
        self._update_trans_box_visuals()
    
    def _schedule_column_fit(self, *args):
        """Fit the columns to their contents once the pending changes of the
        rows are done (they come as several signals)."""
        if not self._column_fit_pending:
            self._column_fit_pending = True
            QTimer.singleShot(0, self._fit_columns)

    def _fit_columns(self):
        self._column_fit_pending = False
        self.table_view.resizeColumnsToContents()

    def _on_manual_edit(self, text):
        """
        Called ONLY when the user types in the box. 
        Updates the 'Base' text which serves as the anchor for fluid prepending.
        """
        if not text.strip():
            self._manual_trans_edit = False
            self._manual_base_text = ""
        else:
            self._manual_trans_edit = True
            self._manual_base_text = text

        # Add helper to set style based on state
        self._update_trans_box_visuals()

    def set_session(self, session):
        """
        Updates the table with data from the session and handles row selection logic.
        Prioritizes the 'latest added' component if a new one is detected.
        """
        if session:
            self.setWindowTitle(f"System Inspector: {session.name}")
        else:
            self.setWindowTitle("System Inspector")

        # 1. Capture previous state (Selection & Max ID)
        sel = self.table_view.selectionModel().selectedRows()
        prev_uuid = None
        if sel:
            idx_source = self.proxy_model.mapToSource(sel[0])
            comp = self.table_model.get_component_at(idx_source.row())
            if comp: prev_uuid = comp.uuid
        
        # Detect existing max ID to identify if a new component is added later
        old_comps = self.table_model._components
        max_id_old = max([c.id for c in old_comps], default=0) if old_comps else 0

        # 2. Update Model with new data
        self.current_session = session
        new_comps = session.systs.components if session and session.systs else []
        
        # Extract Constraints Map from the Session
        constraints = {}
        if session and session.systs and session.systs.constraint_model:
            constraints = session.systs.constraint_model.v2_constraints_by_uuid
            
        self.table_model.update_data(new_comps, constraints)
        self._update_column_visibility()
        
        # 3. Determine Selection Target
        max_id_new = max([c.id for c in new_comps], default=0) if new_comps else 0
        target_uuid = prev_uuid
        
        # Rule: If the max ID increased, a component was likely added. Select it.
        # [FIX] Added 'old_comps and' check. 
        # This ensures we only auto-select the newest ID if we are UPDATING a session 
        # (e.g. adding a line), not when initially LOADING one.
        if old_comps and max_id_new > max_id_old:
            for c in new_comps:
                if c.id == max_id_new:
                    target_uuid = c.uuid
                    break
        
        # 4. Apply Selection (Find row by UUID)
        selection_applied = False
        if target_uuid:
            for r in range(self.table_model.rowCount()):
                comp = self.table_model.get_component_at(r)
                if comp and comp.uuid == target_uuid:
                    idx_src = self.table_model.index(r, 0)
                    idx_proxy = self.proxy_model.mapFromSource(idx_src)
                    
                    if idx_proxy.isValid():
                        self.table_view.selectRow(idx_proxy.row())
                        self.table_view.scrollTo(idx_proxy) # Ensure it's visible
                        selection_applied = True
                    break

        # 5. Fallback (Default to top row if no target found or table empty)
        if not selection_applied:
            if self.proxy_model.rowCount() > 0: self.table_view.selectRow(0)
            else: self.vel_plot.plot_spectrum(None)
                
        # --- Auto-Refresh Corner Plot if open ---
        if hasattr(self, 'corner_win') and self.corner_win.isVisible():
            # Grab the actively selected row
            sel = self.table_view.selectionModel().selectedRows()
            if sel:
                idx_source = self.proxy_model.mapToSource(sel[0])
                primary_comp = self.table_model.get_component_at(idx_source.row())
                
                if primary_comp:
                    meta = session.systs._data.meta if session and session.systs else {}
                    # Check if this newly selected component has Bayesian results
                    if 'bayesian_results' in meta and primary_comp.uuid in meta['bayesian_results']:
                        self._corner_uuid = primary_comp.uuid # Update our tracker
                        res = meta['bayesian_results'][primary_comp.uuid]
                        samples = res['samples']
                        labels = res.get('labels', [f"p_{i}" for i in range(samples.shape[1])])
                        title = f"Bayesian Corner Plot: {primary_comp.series} at z={primary_comp.z:.5f}"
                        
                        self.corner_win.update_data(samples, labels, title)

    def update_limit_boxes(self, v_min, v_max):
        self.vmin_in.blockSignals(True)
        self.vmax_in.blockSignals(True)
        self.vmin_in.setText(f"{v_min:.0f}")
        self.vmax_in.setText(f"{v_max:.0f}")
        self.vmin_in.blockSignals(False)
        self.vmax_in.blockSignals(False)

    def update_z_box(self, z_val):
        self.z_in.blockSignals(True)
        self.z_in.setText(f"{z_val:.5f}")
        self.z_in.blockSignals(False)

    def _update_cursor(self):
        try:
            n = float(self.logn_in.text())
            b = float(self.b_in.text())
            self.vel_plot.update_cursor_params(n, b)
        except ValueError: pass

    def _apply_limits(self):
        try:
            vmin = float(self.vmin_in.text())
            vmax = float(self.vmax_in.text())
            if vmin < vmax: self.vel_plot.set_velocity_limits(vmin, vmax)
        except ValueError: pass

    def _apply_z(self):
        try:
            z_target = float(self.z_in.text())
            self.vel_plot.set_center_redshift(z_target)
        except ValueError as e: logging.error(f"Error applying z: {e}")

    def _on_group_toggled(self, checked):
        self.proxy_model.set_group_filter_enabled(checked)

    def _on_selection_changed(self, sel, desel):
        indexes = self.table_view.selectionModel().selectedRows()
        if not indexes: 
            # Clear group if selection is cleared
            self._update_group_definition([])
            return
        
        idx_proxy = indexes[0]
        idx_source = self.proxy_model.mapToSource(idx_proxy)
        primary_row = idx_source.row()
        primary_comp = self.table_model.get_component_at(primary_row)
        
        selected_comps = []
        for idx in indexes:
            src = self.proxy_model.mapToSource(idx)
            c = self.table_model.get_component_at(src.row())
            if c: selected_comps.append(c)

        # Update the Group Logic (Highlighter + Filter List)
        self._update_group_definition(selected_comps)

        # Update Main Plot Highlights
        if self.main_window:
            self.main_window.update_session_highlights(selected_comps)

        if primary_comp and self.current_session:
            new_txt = ""
            
            # Use dict.fromkeys to preserve order while removing duplicates
            unique_series = list(dict.fromkeys([c.series for c in selected_comps if c.series]))
            series_str = ", ".join(unique_series)

            if not self._manual_trans_edit:
                # --- MODE A: Auto-Replace ---
                # The box reflects the FULL selection.
                new_txt = series_str
                # Update the BASE so if the user starts typing now, 
                # they append to the CURRENT selection.
                self._manual_base_text = new_txt 
                
            else:
                # --- MODE B: Smart Preserve or Reset ---
                # Check what is currently in the custom manual text
                base_tokens = [t.strip() for t in self._manual_base_text.split(',') if t.strip()]
                
                if primary_comp.series in base_tokens:
                    # The selected component is part of the custom list. Preserve the custom panels!
                    new_txt = self._manual_base_text
                else:
                    # The selected component is entirely new. Reset the view!
                    new_txt = series_str
                    self._manual_base_text = new_txt
                    self._manual_trans_edit = False # Turn off manual mode
            
            # Update the widget (Programmatic setText does NOT fire textEdited)
            self.trans_in.setText(new_txt)
            self.z_in.setText(f"{primary_comp.z:.5f}")
            
            # Add helper to set style based on state
            self._update_trans_box_visuals()

            # Update Internal State & Plot
            self._parse(new_txt)
            self.vel_plot.plot_system(self.current_session, primary_comp, selected_comps)
            
    def _apply(self):
        self._parse(self.trans_in.text())
        self._apply_limits()
        self.vel_plot.plot_spectrum()

    def _parse(self, text):
        self.active_transitions = [t.strip() for t in text.split(',') if t.strip()]

    def _forward(self, name, params):
        if self.main_window: self.main_window._on_recipe_requested("absorbers", name, params, {})

    def select_flagged(self):
        """Select the rows of the components with flags of the fit (#517),
        among those shown, e.g. to refit them with Refit Selected."""
        selection = QItemSelection()
        n = 0
        for row in range(self.proxy_model.rowCount()):
            src = self.proxy_model.mapToSource(self.proxy_model.index(row, 0))
            comp = self.table_model.get_component_at(src.row())
            if comp and comp.flags:
                selection.select(self.proxy_model.index(row, 0),
                                 self.proxy_model.index(row, self.proxy_model.columnCount() - 1))
                n += 1
        self.table_view.selectionModel().select(
            selection, QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)
        return n

    def focus_on_component(self, uuid: str, force_group_view: bool = False): # <--- Added Arg
        """
        Public API: Selects a component and scrolls to it.
        Args:
        
        force_group_view: If True, checks 'Only group' at the end (Main Window behavior). If False, restores previous state (Inspector Plot behavior).
        """
        # 1. Capture previous state
        was_checked = self.group_cb.isChecked()
        
        # 2. Temporarily disable grouping to ensure row is found/selectable
        if was_checked: self.group_cb.setChecked(False) 
        
        # 3. Find row in Source Model
        target_row = -1
        for r in range(self.table_model.rowCount()):
            c = self.table_model.get_component_at(r)
            if c.uuid == uuid:
                target_row = r
                break
        
        if target_row == -1: 
            if was_checked: self.group_cb.setChecked(True)
            return

        # 4. Map to Proxy Index
        source_idx = self.table_model.index(target_row, 0)
        proxy_idx = self.proxy_model.mapFromSource(source_idx)
        
        if not proxy_idx.isValid(): 
            if was_checked: self.group_cb.setChecked(True)
            return

        # 5. Select Row
        self.table_view.selectRow(proxy_idx.row())
        self.table_view.scrollTo(proxy_idx)
        
        # 6. [FIX] Logic for the 'Only group' State
        if force_group_view:
            # Case 1: Main Window -> Always Enforce
            self.group_cb.setChecked(True)
        else:
            # Case 2: Velocity Plot -> Restore User Preference
            if was_checked: self.group_cb.setChecked(True)
        
        # 7. Ensure window is active
        self.show()
        self.raise_()
        self.activateWindow()

    def _get_selected_components(self) -> List[ComponentDataV2]:
        rows = self.table_view.selectionModel().selectedRows()
        comps = []
        for idx in rows:
            src = self.proxy_model.mapToSource(idx)
            c = self.table_model.get_component_at(src.row())
            if c: comps.append(c)
        return comps

    def _on_context_menu(self, pos):
        index = self.table_view.indexAt(pos)
        if not index.isValid(): return
        
        # 1. Identify Context (Clicked Cell)
        idx_source = self.proxy_model.mapToSource(index)
        clicked_comp = self.table_model.get_component_at(idx_source.row())
        col_name = COLUMNS[idx_source.column()]
        param_attr = COL_MAP.get(col_name)

        if not clicked_comp: return
        
        selected_comps = self._get_selected_components()
        
        # Guard: If right-clicking outside the selection, reset selection to clicked row
        # (Standard GUI behavior: Right-click selects the target if not already selected)
        if clicked_comp not in selected_comps:
            selected_comps = [clicked_comp]
            
        count = len(selected_comps)
        is_multi = count > 1
        suffix = f" ({count} items)" if is_multi else ""
        
        m = QMenu(self)
        
        # 1. Linking Logic (Strictly 2 items)
        if count == 2 and col_name in ['z', 'logN', 'b', 'btur']:
            # ... (Existing Link Logic - only works for exactly 2) ...
            # Identify Source vs Target (Clicked is Target)
            source_comp = selected_comps[0] if selected_comps[0] != clicked_comp else selected_comps[1]
            
            m.addSection("Parameter Linking")
            
            expr_val = f"p['{source_comp.uuid}'].{param_attr}"
            act_link = QAction(f"Link {col_name} to {source_comp.series} (Value)", m)
            act_link.triggered.connect(
                lambda: self._toggle_constraint(
                    clicked_comp.uuid, param_attr, False, expr_val, source_comp.uuid
                )
            )
            m.addAction(act_link)
            
            if col_name == 'b':
                m_dep = self._get_mass(clicked_comp.series)
                m_src = self._get_mass(source_comp.series)
                if m_dep and m_src:
                    mass_ratio = np.sqrt(m_src / m_dep)
                    expr_therm = f"p['{source_comp.uuid}'].b * {mass_ratio:.4f}"
                    act_therm = QAction(f"Link {col_name} to {source_comp.series} (Thermal)", m)
                    act_therm.triggered.connect(
                        lambda: self._toggle_constraint(
                            clicked_comp.uuid, param_attr, False, expr_therm, source_comp.uuid
                        )
                    )
                    m.addAction(act_therm)
            m.addSeparator()

        # 2. Constraint Actions (Multi-aware)
        if col_name in ['z', 'logN', 'b', 'btur']:
            
            # Allows setting z, logN, b, btur for ALL selected rows to a single value
            act_set_val = QAction(f"Set Value of '{col_name}'...{suffix}", m)
            act_set_val.triggered.connect(
                lambda: self._batch_set_value(selected_comps, col_name)
            )
            m.addAction(act_set_val)
            m.addSeparator()

            if is_multi:
                # --- Batch Actions ---
                act_freeze_sel = QAction(f"Freeze '{col_name}'{suffix}", m)
                act_freeze_sel.triggered.connect(
                    lambda: self._toggle_constraint_list(selected_comps, param_attr, False)
                )
                m.addAction(act_freeze_sel)

                act_free_sel = QAction(f"Unfreeze '{col_name}'{suffix}", m)
                act_free_sel.triggered.connect(
                    lambda: self._toggle_constraint_list(selected_comps, param_attr, True)
                )
                m.addAction(act_free_sel)
                
            else:
                # --- Single Item Context (Smarter toggles) ---
                is_frozen = False
                is_linked = False
                if self.current_session and self.current_session.systs.constraint_model:
                    c_map = self.current_session.systs.constraint_model.v2_constraints_by_uuid
                    if clicked_comp.uuid in c_map and param_attr in c_map[clicked_comp.uuid]:
                        cons = c_map[clicked_comp.uuid][param_attr]
                        if not cons.is_free: is_frozen = True
                        if cons.expression: is_linked = True
                
                if is_linked:
                    act = QAction(f"Unlink '{col_name}' (Set Free)", m)
                    act.triggered.connect(lambda: self._toggle_constraint(clicked_comp.uuid, param_attr, True, None, None))
                    m.addAction(act)
                elif is_frozen:
                    act = QAction(f"Unfreeze '{col_name}'", m)
                    act.triggered.connect(lambda: self._toggle_constraint(clicked_comp.uuid, param_attr, True, None))
                    m.addAction(act)
                else:
                    act = QAction(f"Freeze '{col_name}'", m)
                    act.triggered.connect(lambda: self._toggle_constraint(clicked_comp.uuid, param_attr, False, None))
                    m.addAction(act)

            # Global Column Actions (Always available)
            m.addSeparator()
            act_freeze_col = QAction(f"Freeze All '{col_name}' (Column)", m)
            act_freeze_col.triggered.connect(lambda: self._bulk_toggle_constraint(param_attr, False))
            m.addAction(act_freeze_col)
            
            act_free_col = QAction(f"Unfreeze All '{col_name}' (Column)", m)
            act_free_col.triggered.connect(lambda: self._bulk_toggle_constraint(param_attr, True))
            m.addAction(act_free_col)
            
            m.addSeparator()

        # Set Resolution (Multi-aware via internal logic)
        act_resol = QAction(f"Set Resolving Power...{suffix}", m)
        act_resol.triggered.connect(self._set_component_resolution)
        m.addAction(act_resol)

        m.addSeparator()
        
        # 3. Operations (Multi-aware)
        act_refit = QAction(f"Refit Selected{suffix}", m)
        act_refit.triggered.connect(lambda: self._refit_list(selected_comps))
        m.addAction(act_refit)

        m.addAction("Refit All", self._refit_all)
        
        act_del = QAction(f"Delete{suffix}", m)
        act_del.triggered.connect(lambda: self._delete_list(selected_comps))
        m.addAction(act_del)

        m.addSeparator()
        m.addSection("Bayesian Fitting")
        # emcee and corner are optional: without them the entries are
        # disabled and their tooltip says what to install
        m.setToolTipsVisible(True)
        missing_tip = (f"Needs the optional dependencies: "
                       f"{bayesian_fitter.BAYES_INSTALL_HINT}")

        act_bayesian_mcmc = QAction(f"Refit Bayesian (MCMC){suffix}", m)
        act_bayesian_mcmc.triggered.connect(lambda: self._fit_bayesian(selected_comps, 'mcmc'))
        if bayesian_fitter.emcee is None:
            act_bayesian_mcmc.setEnabled(False)
            act_bayesian_mcmc.setToolTip(missing_tip)
        m.addAction(act_bayesian_mcmc)

        #act_bayesian_nested = QAction(f"Refit Bayesian (Nested){suffix}", m)
        #act_bayesian_nested.triggered.connect(lambda: self._fit_bayesian(selected_comps, 'nested'))
        #m.addAction(act_bayesian_nested)
        
        act_corner = QAction("Show Corner Plot", m)
        
        # Check if Bayesian results exist for this specific component
        has_results = False
        if self.current_session and self.current_session.systs:
            meta = self.current_session.systs._data.meta
            if 'bayesian_results' in meta and clicked_comp.uuid in meta['bayesian_results']:
                has_results = True
                
        act_corner.setEnabled(has_results)
        if not _corner_installed():
            act_corner.setEnabled(False)
            act_corner.setToolTip(missing_tip)
        act_corner.triggered.connect(lambda: self._show_corner_plot(clicked_comp))
        m.addAction(act_corner)
        
        m.exec(self.table_view.mapToGlobal(pos))

    # --- [NEW] Batch Helper Methods ---

    def _toggle_constraint_list(self, comps, param, is_free):
        """Iterates list and applies constraint update."""
        for c in comps:
            self._toggle_constraint(c.uuid, param, is_free, None, None)

    def fit_window_params(self) -> dict:
        """``z_window_kms`` for the manual fits, if the user set one (#518);
        empty, zero or not a number: nothing (the default of each fit)."""
        try:
            value = float(self.fit_window_in.text())
        except ValueError:
            return {}
        return {"z_window_kms": str(value)} if value > 0 else {}

    def _refit_list(self, comps):
        """Batch refit: Treats the entire selection as one fitting problem."""
        if self.main_window and comps:
            uuid_list = [c.uuid for c in comps]
            
            self.main_window._on_recipe_requested(
                "absorbers", "refit_all", 
                {
                    "selected_uuids": uuid_list,
                    "group_depth": "2",  # Standard FoF depth to catch blends
                    "max_nfev": "200",
                    **self.fit_window_params(),
                }, 
                {}
            )

    def _delete_list(self, comps):
        """Batch delete with a single confirmation dialog."""
        if not comps or not self.main_window: return
        
        count = len(comps)
        if count == 1:
            c = comps[0]
            txt = f"Series: {c.series}\nRedshift: {c.z:.5f}"
        else:
            txt = f"You are about to delete {count} components."

        confirm = self.main_window._show_custom_message(
            title="Delete Components",
            header=f"Delete {count} component(s)?",
            text=txt,
            buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            default_btn=QMessageBox.StandardButton.No,
            parent=self
        )
        
        if confirm == QMessageBox.Yes:
            # --- [FIX] Send single batch request with ALL uuids ---
            uuid_list = [c.uuid for c in comps]
            
            self.main_window._on_recipe_requested(
                "absorbers", "delete_component", {"uuids": uuid_list}, {}
            )

    def _fit_bayesian(self, comps, engine='mcmc'):
        """Batch Bayesian fit."""
        if self.main_window and comps:
            uuid_list = [c.uuid for c in comps]
            
            # Log the action
            logging.info(f"Triggering Bayesian {engine.upper()} fit for {len(uuid_list)} components...")
            
            self.main_window._on_recipe_requested(
                "absorbers", "fit_bayesian", 
                {
                    "uuids": uuid_list,
                    "engine": engine,
                    "nsteps": "500",
                    "nlive": "250",
                    "group_depth": "2"
                }, 
                {}
            )

    def _show_corner_plot(self, comp):
        """Show Corner Plot for a component if results exist."""
        if not self.current_session: return
        
        meta = self.current_session.systs._data.meta
        if 'bayesian_results' not in meta or comp.uuid not in meta['bayesian_results']:
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.information(self, "No Results", 
                                  f"No Bayesian sampling results found for {comp.series} (z={comp.z:.4f}).\n"
                                  "Please run 'Re-fit Bayesian' first.")
            return
            
        res = meta['bayesian_results'][comp.uuid]
        samples = res['samples']
        labels = res.get('labels', [])
        
        if not labels:
            # Fallback labels
            labels = [f"p_{i}" for i in range(samples.shape[1])]
            
        title = f"Bayesian Corner Plot: {comp.series} at z={comp.z:.5f}"
        # Track which component is currently being plotted
        self._corner_uuid = comp.uuid 
        
        # Reuse existing window if open, otherwise create a new one
        if hasattr(self, 'corner_win') and self.corner_win.isVisible():
            self.corner_win.update_data(samples, labels, title)
            self.corner_win.raise_()
            self.corner_win.activateWindow()
        else:
            # Imported here: corner is an optional dependency
            from .corner_plot import CornerPlotWindow
            self.corner_win = CornerPlotWindow(samples, labels, title=title, parent=self)
            self.corner_win.show()

    # Method to set resolution for specific components
    def _set_component_resolution(self):
        indexes = self.table_view.selectionModel().selectedRows()
        if not indexes: return

        # [NEW] Determine initial value to show
        initial_val = 0.0
        
        # 1. Try fetching from the first selected component
        first_src = self.proxy_model.mapToSource(indexes[0])
        first_comp = self.table_model.get_component_at(first_src.row())
        
        if first_comp and first_comp.resol is not None and first_comp.resol > 0:
            initial_val = first_comp.resol
        elif self.current_session:
            # 2. The resolution of the spectrum at the component (#526)
            r = self.current_session.spec.resolution()
            r = r[np.isfinite(r)]
            if r.size:
                initial_val = float(np.median(r))
                first_line = (STANDARD_MULTIPLETS.get(first_comp.series, [first_comp.series.split(',')[0].strip()])[0]
                              if first_comp else None)
                atom = ATOM_DATA.get(first_line) if first_line else None
                if atom is not None:
                    try:
                        x_obs = (atom['wave'] * (1 + first_comp.z) * au.Angstrom).to_value(
                            self.current_session.spec.x.unit)
                        r_at = resolution_at(self.current_session.spec, np.array([x_obs]))[0]
                        if np.isfinite(r_at):
                            initial_val = float(r_at)
                    except au.UnitConversionError:
                        pass

        # 1. Create Custom Dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Set Component Resolving Power")
        dialog.setMinimumWidth(300)
        
        # Apply strict styling to match other forms
        pal = QApplication.palette()
        text_col = pal.color(pal.ColorRole.Text).name()
        base_col = pal.color(pal.ColorRole.Base).name()
        dialog.setStyleSheet(f"""
            QLineEdit {{
                padding: 4px;
                border-radius: 4px;
                background-color: {base_col};
                color: {text_col};
            }}
            QLineEdit:focus {{
                {FOCUS_BORDER}
            }}
        """)
        layout = QVBoxLayout(dialog)
        layout.setSpacing(10)

        # 2. Form Layout (Label Left, Box Right)
        form_layout = QHBoxLayout()
        label = QLabel("Resolving power R (e.g. 50000):")
        
        self.res_input = QLineEdit()
        if initial_val > 0:
            self.res_input.setText(f"{initial_val:.0f}")
        else:
            #self.res_input.setPlaceholderText("e.g. 50000")
            pass
        validator = QDoubleValidator()
        validator.setBottom(0.0)
        validator.setNotation(QDoubleValidator.StandardNotation)
        self.res_input.setValidator(validator)
        
        form_layout.addWidget(label)
        form_layout.addWidget(self.res_input)
        layout.addLayout(form_layout)
        
        # 3. Standard Buttons
        btn_box = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        btn_box.accepted.connect(dialog.accept)
        btn_box.rejected.connect(dialog.reject)
        layout.addWidget(btn_box)

        # 4. Execute
        if dialog.exec() == QDialog.Accepted:
            try:
                val = float(self.res_input.text())
                
                # Update each selected component
                for idx in indexes:
                    src = self.proxy_model.mapToSource(idx)
                    comp = self.table_model.get_component_at(src.row())
                    if comp:
                         if self.main_window:
                             self.main_window._on_recipe_requested(
                                 "absorbers", "update_component", 
                                 {'uuid': comp.uuid, 'resol': val}, {}
                             )
            except ValueError:
                pass

    # Helper for Mass
    def _get_mass(self, series_name):
        # 1. Direct
        if series_name in ATOM_DATA: return ATOM_DATA[series_name]['mass']
        # 2. Multiplet
        if series_name in STANDARD_MULTIPLETS:
            prim = STANDARD_MULTIPLETS[series_name][0]
            if prim in ATOM_DATA: return ATOM_DATA[prim]['mass']
        return None

    def _toggle_constraint(self, uuid, param, set_free, expression=None, target_uuid=None): # <--- Added Arg
        if self.main_window:
            self.main_window._on_recipe_requested(
                "absorbers", "update_constraint", 
                {
                    "uuid": uuid, 
                    "param": param, 
                    "is_free": set_free,
                    "expression": expression,
                    "target_uuid": target_uuid # <--- Pass it
                }, 
                {}
            )
    
    def _bulk_toggle_constraint(self, param, is_free):
        if self.main_window:
            self.main_window._on_recipe_requested(
                "absorbers", "update_constraints_column",
                {"param": param, "is_free": is_free},
                {}
            )
    
    def _delete(self):
        idx = self.table_view.currentIndex()
        src = self.proxy_model.mapToSource(idx)
        comp = self.table_model.get_component_at(src.row())
        if comp and self.main_window:
            confirm = self.main_window._show_custom_message(
                title="Delete Component",
                header="Delete this component?",
                text=f"Series: {comp.series}\nRedshift: {comp.z:.5f}",
                buttons=QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                default_btn=QMessageBox.StandardButton.No,
                parent=self
            )
            if confirm == QMessageBox.Yes:
                self.main_window._on_recipe_requested(
                    "absorbers", "delete_component", {"uuid": comp.uuid}, {})

    def _refit(self):
        idx = self.table_view.currentIndex()
        src = self.proxy_model.mapToSource(idx)
        comp = self.table_model.get_component_at(src.row())
        if comp and self.main_window:
            self.main_window._on_recipe_requested(
                "absorbers", "fit_component",
                {"uuid": comp.uuid, **self.fit_window_params()}, {})
    
    def _refit_all(self):
        if self.main_window:
            # Trigger the newly created recipe
            self.main_window._on_recipe_requested(
                "absorbers", "refit_all", self.fit_window_params(), {}
            )
    
    def _update_group_definition(self, selected_comps: List[ComponentDataV2]):
        if not selected_comps or not self.current_session:
            # Clear Highlights
            self.table_model.set_highlighted_uuids(set(), set())
            self.proxy_model.set_allowed_uuids(set())
            return

        # 1. Identify Primary UUIDs (The ones explicitly selected)
        primary_uuids = {c.uuid for c in selected_comps}
        
        # 2. Identify Full Group
        group_list = self.current_session.systs.get_connected_group(
            list(primary_uuids), spec=getattr(self.current_session, 'spec', None))
        all_group_uuids = set(group_list)
        
        # 3. Identify Secondary UUIDs (Group members NOT selected)
        secondary_uuids = all_group_uuids - primary_uuids

        # 4. Update Model Colors
        self.table_model.set_highlighted_uuids(primary_uuids, secondary_uuids)
        
        # 5. Update Proxy Filter (Show entire group)
        self.proxy_model.set_allowed_uuids(all_group_uuids)

    @staticmethod
    def _field_label(text: str, tip: str) -> QLabel:
        label = QLabel(text)
        label.setToolTip(tip)
        return label

    @staticmethod
    def _separator() -> QFrame:
        """A thin vertical line between the groups of fields."""
        line = QFrame()
        line.setFrameShape(QFrame.Shape.VLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        return line

    # --- The columns shown ---

    def _hidden_columns(self) -> dict:
        """The columns shown or hidden by the user ({name: hidden}), from
        the settings; a column not in it follows its default."""
        raw = QSettings("Astrocook", "AstrocookV2").value(HIDDEN_COLUMNS_KEY, "")
        try:
            return {k: bool(v) for k, v in json.loads(raw).items()} if raw else {}
        except (TypeError, ValueError):
            return {}

    def _btur_in_use(self) -> bool:
        """Whether a component has btur other than 0, or linked. Not "free":
        the constraints of some sessions mark btur free at 0 everywhere,
        where it does nothing; a fit that moves it makes it non-zero."""
        model = self.table_model
        for comp in model._components:
            if comp.btur:
                return True
            constraint = model._constraints_map.get(comp.uuid, {}).get('btur')
            if constraint is not None and (constraint.target_uuid or constraint.expression):
                return True
        return False

    def _update_column_visibility(self):
        hidden = self._hidden_columns()
        for i, name in enumerate(COLUMNS):
            if name in hidden:
                hide = hidden[name]
            else:
                hide = name == 'btur' and not self._btur_in_use()
            self.table_view.setColumnHidden(i, hide)

    def _on_header_menu(self, pos):
        """Show or hide a column; the choice is saved in the settings."""
        menu = QMenu(self)
        for i, name in enumerate(COLUMNS):
            if name == 'ID':
                continue
            action = menu.addAction(name)
            action.setCheckable(True)
            action.setChecked(not self.table_view.isColumnHidden(i))
            action.toggled.connect(lambda shown, n=name: self._set_column_shown(n, shown))
        menu.addSeparator()
        menu.addAction("Default Columns", self._reset_columns)
        menu.exec(self.table_view.horizontalHeader().mapToGlobal(pos))

    def _set_column_shown(self, name: str, shown: bool):
        hidden = self._hidden_columns()
        hidden[name] = not shown
        QSettings("Astrocook", "AstrocookV2").setValue(HIDDEN_COLUMNS_KEY, json.dumps(hidden))
        self._update_column_visibility()

    def _reset_columns(self):
        QSettings("Astrocook", "AstrocookV2").remove(HIDDEN_COLUMNS_KEY)
        self._update_column_visibility()

    def _apply_field_styles(self):
        """The fields above the velocity plot, with the colours of the
        palette (again when it changes, dark mode)."""
        pal = QApplication.palette()
        style = f"""
            QLineEdit {{ 
                padding: 3px; 
                border-radius: 4px; 
                background: {pal.color(pal.ColorRole.Base).name()}; 
                color: {pal.color(pal.ColorRole.Text).name()}; }}
            QLineEdit:focus {{
                {FOCUS_BORDER}
            }}
        """
        for w in [self.z_in, self.vmin_in, self.vmax_in, self.logn_in, self.b_in,
                  self.fit_window_in]:
            w.setStyleSheet(style)
        # The transition box has its own style (automatic or manual)
        self._update_trans_box_visuals()

    def event(self, event):
        """The system switched between light and dark: the colours that
        come from the palette are applied again (Qt sends the change to the
        windows only, see MainWindowV2.event)."""
        if event.type() == QEvent.Type.ApplicationPaletteChange and hasattr(self, 'trans_in'):
            self._apply_field_styles()
            self.table_model.refresh_colours()
            self.vel_plot.toolbar.refresh_icons()
        return super().event(event)

    def _update_trans_box_visuals(self):
        """Updates the text color of the Trans box to indicate Auto vs Manual mode."""
        pal = QApplication.palette()
        base_col = pal.color(pal.ColorRole.Base).name()
        text_col = pal.color(pal.ColorRole.Text)
        
        # Calculate colors
        if self._manual_trans_edit:
            # Manual Mode: Standard Full Color
            col_str = text_col.name()
            font_weight = "bold"
        else:
            # Auto Mode: Dimmed (50% Alpha)
            # We can't use hex alpha easily in all Qt versions, so we use rgba
            col_str = f"rgba({text_col.red()},{text_col.green()},{text_col.blue()},0.6)"
            font_weight = "normal"

        self.trans_in.setStyleSheet(f"""
            QLineEdit {{
                padding: 3px;
                border-radius: 4px;
                background: {base_col};
                color: {col_str};
                font-weight: {font_weight};
            }}
            QLineEdit:focus {{
                {FOCUS_BORDER}
                color: {text_col.name()}; /* Always bright when focused for editing */
            }}
        """)

    def _batch_set_value(self, comps, col_name):
        """
        Prompts for a value (enforcing '.' decimal) and batch-updates components.
        """
        if not comps: return

        # 1. Determine default value
        first_val = getattr(comps[0], col_name, 0.0)
        
        # 2. Create Custom Dialog to enforce Locale
        dialog = QInputDialog(self)
        dialog.setWindowTitle(f"Batch Edit {col_name}")
        dialog.setLabelText(f"Enter new {col_name} for {len(comps)} components:")
        dialog.setInputMode(QInputDialog.DoubleInput)
        dialog.setDoubleValue(float(first_val))
        dialog.setDoubleDecimals(5)
        
        # [FIX] Enforce Standard 'C' Locale (Dot separator)
        dialog.setLocale(QLocale(QLocale.C)) 
        
        # Resize slightly for comfort
        dialog.resize(300, 150)

        if dialog.exec() == QDialog.Accepted:
            val = dialog.doubleValue()
            
            if self.main_window:
                # 3. Use the new BULK recipe
                # We collect all UUIDs and fire ONE recipe request.
                uuid_list = [c.uuid for c in comps]
                
                self.main_window._on_recipe_requested(
                    "absorbers", "update_components", 
                    {'uuids': uuid_list, col_name: val}, 
                    {}
                )