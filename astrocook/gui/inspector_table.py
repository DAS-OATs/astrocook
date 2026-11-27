"""The table of the System Inspector: its columns, the cues of the flags of
the fit, the model of the components, the delegate of the link colours and
the proxy that filters and sorts the rows.

Moved out of ``system_inspector.py`` in phase 5 (split of the GUI).
"""
import html
import re

import numpy as np
from PySide6.QtCore import Qt, QAbstractTableModel, QModelIndex, Signal, QSortFilterProxyModel
from PySide6.QtGui import QBrush, QColor, QFont, QPalette
from PySide6.QtWidgets import QStyledItemDelegate
from typing import Any, Dict, List, Optional, Set

from ..core.structures import ComponentData
from ..fitting.component_flags import large_error_params
from .plot_style import plot_colour
from .theme import FLAG_BOUND_COLOR, FLAG_NAN_COLOR, FROZEN_TEXT, link_colours  # noqa: F401 (re-exported)

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
# The invalid index, the default parent of a flat table (one instance, not
# a call in each default argument: ruff B008)
_NO_PARENT = QModelIndex()
EDITABLE_COLS = {"z", "logN", "b", "btur"}

#: The speed of light (km/s), for the error of z in velocity
C_KMS = 299792.458


def dz_kms(z: float, dz: float) -> float:
    """The error of a redshift as a velocity, c dz / (1 + z)."""
    return C_KMS * dz / (1.0 + z)


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
#: The values that are not measured (phase 7, PR 7d.2): the parameter of
#: each flag (None: the whole component) and one line for the tooltip. No
#: colour: the warning sign in the cell of the parameter
_MEASUREMENT_FLAGS = {
    'unresolved': ('b', "Unresolved: b below half the instrumental b, b is not measured."),
    'saturated': ('logN', "Saturated: the weakest line of the series is saturated, "
                          "logN is not measured."),
    'large_errors': (None, "Errors too large to have a meaning."),
    'telluric_overlap': (None, "On a telluric absorption: the errors are not to be trusted."),
}
_NAN_MESSAGE = "The uncertainties could not be computed: the value is not constrained."
_NOT_FITTED_MESSAGE = ("The fit of its group failed (Refit All): the values are those from "
                       "before the fit; refit it, e.g. with Refit Selected.")


def _measurement_params(comp, flag: str) -> tuple:
    """The parameters whose cells carry the warning sign of a flag of 7d.2
    (``large_errors``: those whose errors are too large)."""
    if flag == 'large_errors':
        return large_error_params(comp)
    param = _MEASUREMENT_FLAGS[flag][0]
    return (param,) if param else ()


def measurement_sign(comp, attr: str) -> bool:
    """Whether the cell of the parameter ``attr`` carries the warning sign
    of a value not measured (7d.2)."""
    return any(attr in _measurement_params(comp, f)
               for f in (getattr(comp, 'flags', ()) or ()) if f in _MEASUREMENT_FLAGS)


def flag_cue(comp, attr: Optional[str] = None):
    """Colour and message of the flags of ``comp`` (#517): for a parameter
    ``attr``, its flags (a bound; a value not measured, 7d.2, which have
    no colour); without, all its flags. ``(None, [])`` if there is none."""
    flags = getattr(comp, 'flags', ()) or ()
    messages, color = [], None
    for flag in flags:
        if flag in _BOUND_FLAGS:
            param, text = _BOUND_FLAGS[flag]
            if attr is None or attr == param:
                messages.append(text)
                color = plot_colour('flag_bound')
        elif flag in _MEASUREMENT_FLAGS:
            if attr is None or attr in _measurement_params(comp, flag):
                messages.append(_MEASUREMENT_FLAGS[flag][1])
        elif flag == 'nan_errors' and attr is None:
            messages.append(_NAN_MESSAGE)
            color = color or plot_colour('flag_nan')
        elif flag == 'not_fitted' and attr is None:
            # phase 6-bis, PR 4b
            messages.append(_NOT_FITTED_MESSAGE)
            color = color or plot_colour('flag_nan')
    return color, messages


_LINK_REF = re.compile(r"""p\[(['"])([^'"]+)\1\]""")


def pretty_expression(expression: str, components, name=None) -> str:
    """The expression of a link with ``p['uuid']`` written as the series
    of the component (or ``name(component)``), HTML-escaped."""
    by_uuid = {c.uuid: c for c in components}

    def replace(match):
        comp = by_uuid.get(match.group(2))
        if comp is None:
            return "Unknown"
        return name(comp) if name else html.escape(comp.series)
    parts, last = [], 0
    for m in _LINK_REF.finditer(expression):
        parts.append(html.escape(expression[last:m.start()]))
        parts.append(replace(m))
        last = m.end()
    parts.append(html.escape(expression[last:]))
    return "".join(parts)


def tooltip_line(comp, attr: str, constraint, components) -> str:
    """The one line of the tooltip of a cell (7d.9): the reason of its
    warning sign or colour, else the state of the parameter; the rest is in
    the Details panel."""
    _, messages = flag_cue(comp, attr if attr in ERR_MAP else None)
    if attr not in ERR_MAP and attr != 'id':
        messages = []
    if len(messages) == 1:
        # Rich text, so that Qt wraps a long message instead of running off
        # the screen
        return f"<qt>\u26a0 {html.escape(messages[0])}</qt>"
    if messages:
        names = ", ".join(f.replace('_', ' ') for f in (getattr(comp, 'flags', ()) or ()))
        return f"<qt>\u26a0 {len(messages)} warnings ({html.escape(names)}): see Details</qt>"
    if attr not in ERR_MAP:
        return ""
    target = getattr(constraint, 'target_uuid', None)
    expression = getattr(constraint, 'expression', None)
    if target or expression:
        if expression:
            return f"Linked: {pretty_expression(expression, components)}"
        ref = next((c for c in components if c.uuid == target), None)
        return f"Linked to {html.escape(ref.series) if ref else 'Unknown'}"
    if constraint is not None and not constraint.is_free:
        return "Frozen: not fitted"
    return "Free"


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
_DISPLAY, _EDIT, _TOOLTIP, _FONT = (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.EditRole, Qt.ItemDataRole.ToolTipRole,
                                    Qt.ItemDataRole.FontRole)
_FOREGROUND, _BACKGROUND, _ALIGNMENT, _USER = (Qt.ItemDataRole.ForegroundRole, Qt.ItemDataRole.BackgroundRole,
                                               Qt.ItemDataRole.TextAlignmentRole, Qt.ItemDataRole.UserRole)
_ALIGN_CENTER, _HORIZONTAL = Qt.AlignmentFlag.AlignCenter, Qt.Orientation.Horizontal
_ITEM_IS_EDITABLE, _NO_ITEM_FLAGS = Qt.ItemFlag.ItemIsEditable, Qt.ItemFlag.NoItemFlags
_HIGHLIGHTED_TEXT = QPalette.ColorRole.HighlightedText

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

    def __init__(self, components: Optional[List[ComponentData]] = None):
        super().__init__()
        self._components = components if components else []
        self._primary_uuids: Set[Any] = set()
        self._secondary_uuids: Set[Any] = set()

        # Store constraints: {uuid: {param: ConstraintObj}}
        self._constraints_map: Dict[Any, Any] = {} 
        # Link colours of the selected group: {(target_uuid, param): index
        # in the link colours}; the index, so that a link keeps its colour
        # when the palette switches between light and dark
        self._link_color_map: Dict[Any, Any] = {}

    def update_data(self, components: List[ComponentData], constraints_map: Optional[dict] = None):
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

            # Color palette (Deep Blue, Magenta, Forest Green, Purple, Teal,
            # Maroon; lighter on a dark palette, theme.link_colours)
            n_colours = len(link_colours())
            kept = {c: previous[c] for c in link_clusters if c in previous}
            used = set(kept.values())
            free = [i for i in range(n_colours) if i not in used] or list(range(n_colours))
            # New links in a fixed order (the set order depends on the UUIDs)
            new = sorted(c for c in link_clusters if c not in kept)
            for i, cluster in enumerate(new):
                kept[cluster] = free[i % len(free)]
            self._link_color_map = kept

    def refresh_colours(self):
        """Draw the table again with the colours of the palette (it switched
        between light and dark)."""
        if self._components:
            self.dataChanged.emit(self.index(0, 0),
                                  self.index(len(self._components) - 1, len(COLUMNS) - 1),
                                  [_FOREGROUND, _BACKGROUND])

    def rowCount(self, parent=_NO_PARENT): return len(self._components)
    def columnCount(self, parent=_NO_PARENT): return len(COLUMNS)
    
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
            except (TypeError, ValueError): return str(val)

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
                c = QColor(plot_colour('flag_nan')); c.setAlpha(90)
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
                        return QBrush(QColor(link_colours()[self._link_color_map[cluster]]))
                
                # If it's the TARGET variable...
                elif is_free:
                    cluster = (comp.uuid, attr)
                    if cluster in self._link_color_map:
                        return QBrush(QColor(link_colours()[self._link_color_map[cluster]]))

            return None
        
        # 2. Tooltip: one line (7d.9), the rest in the Details panel
        if role == _TOOLTIP:
            return tooltip_line(comp, attr, c_data, self._components) or None

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
                # The warning sign in the cell of a value not measured
                # (7d.2); the bounds have their colour
                sign = "\u26a0 " if comp.flags and measurement_sign(comp, attr) else ""
                
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
                            if attr == 'z':
                                # In velocity: dz in scientific notation was
                                # hard to compare between rows (visual review
                                # of phase 5-bis); dz is in the tooltip
                                err_str = f"{dz_kms(comp.z, err_val):.2f} km/s"
                            else:
                                err_str = err_fmt.format(err_val)
                            return f"{sign}{val_str} ± {err_str}"
                
                return f"{sign}{val_str}"

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
        # begin/endFilterChange: invalidateFilter() is deprecated since Qt
        # 6.10 (the minimum version of PySide6 in pyproject.toml)
        self.beginFilterChange()
        self.group_filtering_enabled = enabled
        self.endFilterChange()
        
    def set_allowed_uuids(self, uuids: set):
        """Updates the list of visible UUIDs (Fluid Group)."""
        # Break recursion: If the group hasn't changed, do nothing.
        if self.allowed_uuids == uuids: return
        # Only re-filter if the checkbox is actually checked
        if self.group_filtering_enabled:
            self.beginFilterChange()
            self.allowed_uuids = uuids
            self.endFilterChange()
        else:
            self.allowed_uuids = uuids

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
