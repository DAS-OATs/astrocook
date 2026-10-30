"""The table of the System Inspector: its columns, the cues of the flags of
the fit, the model of the components, the delegate of the link colours and
the proxy that filters and sorts the rows.

Moved out of ``system_inspector.py`` in phase 5 (split of the GUI).
"""
import numpy as np
from PySide6.QtCore import Qt, QAbstractTableModel, QModelIndex, Signal, QSortFilterProxyModel
from PySide6.QtGui import QBrush, QColor, QFont, QPalette
from PySide6.QtWidgets import QStyledItemDelegate
from typing import List

from ..core.structures import ComponentDataV2
from .theme import FLAG_BOUND_COLOR, FLAG_NAN_COLOR, FROZEN_TEXT, LINK_COLOURS

# --- 1. Table Configuration (UPDATED) ---
# We removed the separate error columns (dz, dlogN, db)
COL_MAP = {
    "ID": "id", 
    "Trans.": "series", 
    "z": "z", 
    "logN": "logN", 
    "b": "b", 
    "btur": "btur"
}
COLUMNS = list(COL_MAP.keys())
EDITABLE_COLS = {"z", "logN", "b", "btur"}

# --- Flags of the fit (#517) ---
#: Parameter of each bound flag, and what it means for the user
_BOUND_FLAGS = {
    'z_at_window': ('z', "z at the edge of the fit window: the line may lie outside it; "
                         "refit with a wider window (Fit \u00b1)."),
    'logN_at_lower_bound': ('logN', "logN at the lower bound of the fit: the line is negligible."),
    'logN_at_upper_bound': ('logN', "logN at the upper bound of the fit (19 in the automatic "
                                    "pipelines): the line may be stronger, or stand for several "
                                    "lines; refit with Refit Selected (bound 23)."),
    'b_at_lower_bound': ('b', "b at the lower bound of the fit (1 km/s): the line is probably "
                              "unresolved or saturated."),
    'b_at_upper_bound': ('b', "b at the upper bound of the fit (50 km/s in the automatic "
                              "pipelines): the line may be broader; refit with Refit Selected "
                              "(bound 300 km/s)."),
}
_NAN_MESSAGE = "The uncertainties could not be computed: the value is not constrained."


def flag_cue(comp, attr: str = None):
    """Colour and message of the flags of ``comp`` (#517): for a parameter
    ``attr``, its bound flag; without, all its flags. ``(None, [])`` if
    there is none."""
    flags = getattr(comp, 'flags', ()) or ()
    messages, color = [], None
    for flag in flags:
        if flag in _BOUND_FLAGS:
            param, text = _BOUND_FLAGS[flag]
            if attr is None or attr == param:
                messages.append(text)
                color = FLAG_BOUND_COLOR
        elif flag == 'nan_errors' and attr is None:
            messages.append(_NAN_MESSAGE)
            color = color or FLAG_NAN_COLOR
    return color, messages


# Mapping to find the error attribute for a given value attribute
ERR_MAP = {
    "z": "dz",
    "logN": "dlogN",
    "b": "db",
    "btur": "dbtur"
}

# The Qt enums used by the table, read once: through the Qt namespace each
# access costs ~4 us with PySide6, and the model is called tens of thousands
# of times when the table is drawn, sorted or its columns fitted
_DISPLAY, _EDIT, _TOOLTIP, _FONT = (Qt.DisplayRole, Qt.EditRole, Qt.ToolTipRole,
                                    Qt.FontRole)
_FOREGROUND, _BACKGROUND, _ALIGNMENT, _USER = (Qt.ForegroundRole, Qt.BackgroundRole,
                                               Qt.TextAlignmentRole, Qt.UserRole)
_ALIGN_CENTER, _HORIZONTAL = Qt.AlignCenter, Qt.Horizontal
_ITEM_IS_EDITABLE, _NO_ITEM_FLAGS = Qt.ItemIsEditable, Qt.NoItemFlags
_HIGHLIGHTED_TEXT = QPalette.HighlightedText

#: The roles for which the table model has data; Qt asks for many others for
#: each cell (decoration, check state, size hint...), answered at once
_TABLE_ROLES = frozenset({_DISPLAY, _EDIT, _TOOLTIP, _FONT, _FOREGROUND,
                          _BACKGROUND, _ALIGNMENT, _USER})

class LinkColorDelegate(QStyledItemDelegate):
    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        # Fetch the custom text color we defined in the Model
        fg_brush = index.data(_FOREGROUND)
        if fg_brush:
            # Force Qt to use our color even when the row is natively selected!
            option.palette.setBrush(_HIGHLIGHTED_TEXT, fg_brush)

class SystemTableModel(QAbstractTableModel):
    data_changed_request = Signal(str, dict)

    def __init__(self, components: List[ComponentDataV2] = None):
        super().__init__()
        self._components = components if components else []
        self._primary_uuids = set()
        self._secondary_uuids = set()

        # Store constraints: {uuid: {param: ConstraintObj}}
        self._constraints_map = {} 
        # Link colours of the selected group: {(target_uuid, param): QColor}
        self._link_color_map = {}

    def update_data(self, components: List[ComponentDataV2], constraints_map: dict = None):
        """
        Updates the model with new components and constraint data.
        """
        self.beginResetModel()
        self._components = components
        self._constraints_map = constraints_map if constraints_map else {}
        # The links may have changed while the group stays selected, and
        # set_highlighted_uuids skips an unchanged selection (F13)
        self._update_link_colors()
        self.endResetModel()

    def set_highlighted_uuids(self, primary_uuids: set, secondary_uuids: set):
        if self._primary_uuids == primary_uuids and self._secondary_uuids == secondary_uuids:
            return
            
        self._primary_uuids = primary_uuids
        self._secondary_uuids = secondary_uuids
        self._update_link_colors()

        if self._components:
            top_left = self.index(0, 0)
            bottom_right = self.index(len(self._components) - 1, len(COLUMNS) - 1)
            self.dataChanged.emit(top_left, bottom_right, [_BACKGROUND])

    def _update_link_colors(self):
        """Map a colour to each link of the selected group.

        A link keeps its colour while it stays in the map, so that the colours
        do not change when another link is added or removed.
        """
        previous = self._link_color_map
        self._link_color_map = {}
        active_uuids = self._primary_uuids | self._secondary_uuids
        
        if active_uuids:
            link_clusters = set()
            for uid in active_uuids:
                if uid in self._constraints_map:
                    for attr, cons in self._constraints_map[uid].items():
                        # Track the target UUID and the parameter being linked
                        if cons.target_uuid:
                            link_clusters.add((cons.target_uuid, attr))

            # Color palette (Deep Blue, Magenta, Forest Green, Purple, Teal, Maroon)
            palette = LINK_COLOURS
            kept = {c: previous[c] for c in link_clusters if c in previous}
            used = [color.name() for color in kept.values()]
            free = [color for color in palette if color not in used] or palette
            # New links in a fixed order (the set order depends on the UUIDs)
            new = sorted(c for c in link_clusters if c not in kept)
            for i, cluster in enumerate(new):
                kept[cluster] = QColor(free[i % len(free)])
            self._link_color_map = kept

    def rowCount(self, parent=QModelIndex()): return len(self._components)
    def columnCount(self, parent=QModelIndex()): return len(COLUMNS)
    
    def flags(self, index):
        if not index.isValid(): return _NO_ITEM_FLAGS
        return super().flags(index) | _ITEM_IS_EDITABLE if COLUMNS[index.column()] in EDITABLE_COLS else super().flags(index)

    def data(self, index, role=Qt.DisplayRole):
        # Called for every cell and role when the table is drawn or sorted:
        # the cheap answers first
        if role not in _TABLE_ROLES or not index.isValid(): return None
        if role == _ALIGNMENT: return _ALIGN_CENTER
        comp = self._components[index.row()]
        col_name = COLUMNS[index.column()]
        attr = COL_MAP[col_name]

        # Sorting (raw float)
        if role == _USER:
            val = getattr(comp, attr, 0)
            if val is None: return -np.inf
            try: return float(val)
            except: return str(val)

        # 1. Background: the parameter at a bound of the fit (amber), the
        # ID of a component without uncertainties (red), #517; then the
        # group highlighting
        if role == _BACKGROUND:
            if attr in ('z', 'logN', 'b'):
                color, _ = flag_cue(comp, attr)
                if color:
                    c = QColor(color); c.setAlpha(110)
                    return c
            elif attr == 'id' and 'nan_errors' in (comp.flags or ()):
                c = QColor(FLAG_NAN_COLOR); c.setAlpha(90)
                return c
            if comp.uuid in self._secondary_uuids:
                return QColor(245, 161, 0, 50) 
            return None

        # --- CONSTRAINTS CHECK ---
        # (We calculate this upfront so we can use it for Font, Tooltip, AND Display)
        is_constrained_col = col_name in ['z', 'logN', 'b', 'btur']
        c_data = None
        is_free = True
        is_linked = False

        if is_constrained_col and comp.uuid in self._constraints_map:
            c_data = self._constraints_map[comp.uuid].get(attr)
            if c_data:
                is_free = c_data.is_free
                is_linked = (c_data.expression is not None) or (c_data.target_uuid is not None)

        if role == _FONT and is_constrained_col:
            font = QFont()
            if is_linked:
                font.setBold(True); return font
            elif not is_free:
                font.setItalic(True); return font
            return None

        # --- Foreground (Text Color) ---
        if role == _FOREGROUND and is_constrained_col:
            # 1. Global Behavior: Frozen parameters are always gray
            if not is_free and not is_linked:
                return QBrush(QColor(FROZEN_TEXT))

            # 2. Dynamic Behavior: Coloring Links within the selected group
            if comp.uuid in self._primary_uuids or comp.uuid in self._secondary_uuids:
                # If it's the dependent variable...
                if is_linked and c_data.target_uuid:
                    cluster = (c_data.target_uuid, attr)
                    if cluster in self._link_color_map:
                        return QBrush(self._link_color_map[cluster])
                
                # If it's the TARGET variable...
                elif is_free:
                    cluster = (comp.uuid, attr)
                    if cluster in self._link_color_map:
                        return QBrush(self._link_color_map[cluster])

            return None
        
        # 2. Tooltip (Keep existing logic)
        if role == _TOOLTIP:
             c_chi2 = f"{comp.chi2:.2f}" if comp.chi2 is not None else "N/A"
             c_resol = f"{comp.resol:.0f}" if comp.resol is not None else "spectrum"
             tooltip_html = (
                f"<div style='font-size: 13px'>"
                f"<b>ID:</b> {comp.id} | <b>{comp.series}</b><br>"
                f"<b>z:</b> {comp.z:.6f}<br>"
                f"<b>\u03c7\u00b2:</b> {c_chi2} | <b>R:</b> {c_resol}"
             )
             if is_constrained_col:
                tooltip_html += "<hr style='margin: 4px 0;'>"
                if is_linked:
                    target_name = "Unknown"; target_z = ""
                    if c_data.target_uuid:
                        for c_ref in self._components:
                            if c_ref.uuid == c_data.target_uuid:
                                target_name = c_ref.series; target_z = f" (z={c_ref.z:.5f})"; break
                    import re
                    expr_display = c_data.expression if c_data.expression else "Direct Link"
                    # Clean up expression for display
                    def replacer(match):
                        uid = match.group(1)
                        for c_ref in self._components:
                            if c_ref.uuid == uid: return c_ref.series
                        return "Unknown"
                    expr_pretty = re.sub(r"p\['([^']+)'\]", replacer, expr_display)
                    expr_pretty = re.sub(r'p\["([^"]+)"\]', replacer, expr_pretty)
                    tooltip_html += (f"<b>Status:</b> Linked<br>"
                                     f"<b>Target:</b> {target_name}{target_z}<br>"
                                     f"<b>Formula:</b> <code>{expr_pretty}</code>")
                elif not is_free:
                    tooltip_html += "<b>Status:</b> Frozen (Fixed)"
                else:
                    tooltip_html += "<b>Status:</b> Free"
             _, messages = flag_cue(comp, attr if is_constrained_col else None)
             if messages:
                 tooltip_html += "<hr style='margin: 4px 0;'>" + "<br>".join(
                     f"\u26a0 {m}" for m in messages)
             tooltip_html += "</div>"
             return tooltip_html

        # 3. DISPLAY & EDIT LOGIC (Refined)
        if role == _DISPLAY or role == _EDIT:
            val = getattr(comp, attr, None)
            
            val_fmt = "{}"
            err_fmt = "{}"
            
            if attr == 'z':
                val_fmt = "{:.6f}"; err_fmt = "{:.2e}"
            elif attr == 'logN':
                val_fmt = "{:.3f}"; err_fmt = "{:.3f}"
            elif attr == 'b':
                val_fmt = "{:.2f}"; err_fmt = "{:.2f}"
            elif attr == 'btur':
                val_fmt = "{:.2f}"; err_fmt = "{:.2f}"

            # --- EDIT ROLE: Clean Value Only ---
            if role == _EDIT:
                if isinstance(val, (int, float)):
                    return val_fmt.format(val)
                return str(val) if val is not None else ""
            
            # --- DISPLAY ROLE: Value ± Error ---
            if role == _DISPLAY:
                if val is None: return ""
                if attr == 'id' and comp.flags:
                    return f"\u26a0 {val}"
                
                if isinstance(val, (int, float)):
                    val_str = val_fmt.format(val)
                else:
                    val_str = str(val)

                # Append Error if it exists... AND if not linked!
                if attr in ERR_MAP and isinstance(val, (int, float)):
                    
                    # Suppress error if linked OR frozen
                    # We only show error if the parameter is truly free to vary
                    if not is_linked and is_free: 
                        err_attr = ERR_MAP[attr]
                        err_val = getattr(comp, err_attr, None)
                        
                        if err_val is not None and np.isfinite(err_val) and err_val > 0:
                            err_str = err_fmt.format(err_val)
                            return f"{val_str} ± {err_str}"
                
                return val_str

        return None

    def setData(self, index, value, role=Qt.EditRole):
        if not index.isValid() or role != _EDIT: return False
        attr = COL_MAP[COLUMNS[index.column()]]
        try:
            val = float(value) if attr in ['z', 'logN', 'b', 'btur'] else str(value)
            self.data_changed_request.emit("update_component", {attr: val, 'uuid': self._components[index.row()].uuid})
            return True
        except ValueError: return False

    def headerData(self, section, orientation, role):
        return COLUMNS[section] if role == _DISPLAY and orientation == _HORIZONTAL else None

    def get_component_at(self, row: int): 
        return self._components[row] if 0 <= row < len(self._components) else None
    
# --- Proxy Model ---
class SystemSortFilterProxyModel(QSortFilterProxyModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.group_filtering_enabled = False
        self.allowed_uuids = set()
        self.flagged_only = False

    def set_flagged_only(self, enabled: bool):
        """Show only the components with flags of the fit (#517)."""
        self.flagged_only = enabled
        # invalidate(), not invalidateFilter(): with PySide6 6.12 the latter
        # left rows of the sorted table out after the flags had changed
        self.invalidate()

    def set_group_filter_enabled(self, enabled: bool):
        self.group_filtering_enabled = enabled
        self.invalidateFilter()
        
    def set_allowed_uuids(self, uuids: set):
        """Updates the list of visible UUIDs (Fluid Group)."""
        # Break recursion: If the group hasn't changed, do nothing.
        if self.allowed_uuids == uuids: return
        self.allowed_uuids = uuids

        # Only re-filter if the checkbox is actually checked
        if self.group_filtering_enabled:
            self.invalidateFilter()

    def filterAcceptsRow(self, source_row, source_parent):
        comp = self.sourceModel().get_component_at(source_row)
        if not comp: return False
        if self.flagged_only and not comp.flags:
            return False
        # If group filtering is on, only show components in the Fluid Group
        if not self.group_filtering_enabled: return True
        return comp.uuid in self.allowed_uuids

    def lessThan(self, left, right):
        l_data = self.sourceModel().data(left, _USER)
        r_data = self.sourceModel().data(right, _USER)
        
        # Handle Nones/Types safely
        if l_data is None: l_data = -np.inf
        if r_data is None: r_data = -np.inf
        
        try: 
            return float(l_data) < float(r_data)
        except (ValueError, TypeError): 
            return str(l_data) < str(r_data)
