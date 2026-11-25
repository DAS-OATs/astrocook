"""The *Settings...* dialog (phase 7-quater, PR U2).

Three tabs, *General*, *Appearance* and *Performance*. The settings are about
the interface only (``astrocook/gui/settings.py``): none changes a number,
none applies to scripts, the MCP server or the log. OK and Apply write them
and bring the window up to date; the theme and the colours of the plots go
through the actions of the View menu, which stay in sync. *Restore Defaults*
puts the defaults in the widgets of the dialog (nothing is saved until OK or
Apply).
"""
from typing import Any, Callable

from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                               QDoubleSpinBox, QFileDialog, QFormLayout, QHBoxLayout,
                               QLabel, QLineEdit, QPushButton, QSpinBox, QTabWidget,
                               QVBoxLayout, QWidget)

from . import settings as S
from .plot_style import PLOT_PALETTES
from .theme import APPEARANCES


class SettingsDialog(QDialog):
    def __init__(self, window):
        super().__init__(window)
        self.main_window = window
        self.setWindowTitle("Settings")
        self.setMinimumWidth(520)
        #: widget -> (getter of its value, setter, key): for load / save / reset
        self._fields: list[tuple[str, Callable[[], Any], Callable[[Any], None]]] = []

        layout = QVBoxLayout(self)
        self.tabs = QTabWidget()
        self.tabs.addTab(self._general_tab(), "General")
        self.tabs.addTab(self._appearance_tab(), "Appearance")
        self.tabs.addTab(self._performance_tab(), "Performance")
        layout.addWidget(self.tabs)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Apply
            | QDialogButtonBox.StandardButton.RestoreDefaults)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.StandardButton.Apply).clicked.connect(self.apply)
        buttons.button(QDialogButtonBox.StandardButton.RestoreDefaults).clicked.connect(
            self.restore_defaults)
        self.buttons = buttons
        layout.addWidget(buttons)
        self.load()

    # --- Fields ---
    def _add(self, key: str, getter: Callable[[], Any], setter: Callable[[Any], None]):
        self._fields.append((key, getter, setter))

    def _check(self, key: str, text: str, tip: str = "") -> QCheckBox:
        box = QCheckBox(text)
        box.setToolTip(tip)
        self._add(key, box.isChecked, box.setChecked)
        return box

    def _combo(self, key: str, choices: dict) -> QComboBox:
        """A combo of ``choices`` (value -> label); its value is the key."""
        combo = QComboBox()
        for value, label in choices.items():
            combo.addItem(label, value)
        self._add(key, combo.currentData,
                  lambda v: combo.setCurrentIndex(max(0, combo.findData(v))))
        return combo

    def _spin(self, key: str, limits: tuple, suffix: str = "", special: str = "") -> QSpinBox:
        spin = QSpinBox()
        spin.setRange(*limits)
        spin.setSuffix(suffix)
        if special:
            spin.setSpecialValueText(special)
        self._add(key, spin.value, spin.setValue)
        return spin

    # --- Tabs ---
    def _general_tab(self) -> QWidget:
        tab = QWidget()
        form = QFormLayout(tab)
        self.check_updates = self._check(
            S.CHECK_UPDATES_KEY, "Check for updates at start",
            "Asks GitHub for the latest release (only the version number is "
            "sent, nothing of you or your data). Off by default.")
        self.check_frequency = self._combo(S.CHECK_FREQUENCY_KEY, S.CHECK_FREQUENCIES)
        self.check_updates.toggled.connect(self.check_frequency.setEnabled)
        form.addRow(self.check_updates)
        form.addRow("How often:", self.check_frequency)
        self.max_recent = self._spin(S.MAX_RECENT_KEY, S.RECENT_RANGE)
        self.max_recent.setToolTip("The files shown in File > Open Recent and in the "
                                   "empty window (0: none is remembered).")
        form.addRow("Recent files:", self.max_recent)
        self.confirm_close = self._check(
            S.CONFIRM_CLOSE_KEY, "Ask before closing a session with unsaved changes",
            "If off, a session is closed without the question to save it.")
        form.addRow(self.confirm_close)
        self.dir_mode = self._combo(S.DIR_MODE_KEY, S.DIR_MODES)
        form.addRow("Starting folder of the file dialogs:", self.dir_mode)
        self.fixed_dir = QLineEdit()
        browse = QPushButton("Browse...")
        browse.clicked.connect(self._browse_folder)
        row = QHBoxLayout()
        row.addWidget(self.fixed_dir)
        row.addWidget(browse)
        self.fixed_row = QWidget()
        self.fixed_row.setLayout(row)
        row.setContentsMargins(0, 0, 0, 0)
        self._add(S.FIXED_DIR_KEY, self.fixed_dir.text, self.fixed_dir.setText)
        self.dir_mode.currentIndexChanged.connect(self._sync_enabled)
        form.addRow("Fixed folder:", self.fixed_row)
        return tab

    def _appearance_tab(self) -> QWidget:
        tab = QWidget()
        form = QFormLayout(tab)
        self.theme = self._combo(S.APPEARANCE_KEY, {m: m.capitalize() for m in APPEARANCES})
        form.addRow("Theme:", self.theme)
        self.plot_colours = self._combo(S.PLOT_COLOURS_KEY,
                                        {n: n.capitalize() for n in PLOT_PALETTES})
        form.addRow("Plot colours:", self.plot_colours)
        self.font_size = self._spin(S.FONT_SIZE_KEY, (0, S.FONT_RANGE[1]), " pt", "System")
        self.font_size.setToolTip("The font of tables, panels and menus (sizes written "
                                  "explicitly in some texts do not change).")
        form.addRow("Font size:", self.font_size)
        self.line_width = QDoubleSpinBox()
        self.line_width.setRange(*S.LINE_WIDTH_RANGE)
        self.line_width.setSingleStep(0.1)
        self.line_width.setPrefix("× ")
        self.line_width.setToolTip("A factor on the widths of the lines of the main plot "
                                   "(flux, error, continuum, model).")
        self._add(S.LINE_WIDTH_KEY, self.line_width.value, self.line_width.setValue)
        form.addRow("Line width:", self.line_width)
        self.x_unit = self._combo(S.INITIAL_X_UNIT_KEY, {u: u for u in S.X_UNITS})
        self.x_unit.setToolTip("The unit of the x axis when Astrocook starts "
                               "(the data stay in nm).")
        form.addRow("Initial x unit:", self.x_unit)
        self.decimals = self._spin(S.DECIMALS_KEY, S.DECIMALS_RANGE)
        self.decimals.setToolTip("Decimals of the wavelengths in the Data Inspector "
                                 "(display only: the values are not rounded).")
        form.addRow("Decimals of the tables:", self.decimals)
        return tab

    def _performance_tab(self) -> QWidget:
        tab = QWidget()
        form = QFormLayout(tab)
        self.downsample = self._spin(S.DOWNSAMPLE_KEY, S.DOWNSAMPLE_RANGE, " pixels")
        self.downsample.setSingleStep(100_000)
        self.downsample.setToolTip(
            "When more pixels than this are in view the plot draws 2 points out of "
            "every 10 (min-max), which is faster. The data are never changed.")
        form.addRow("Downsample the plot above:", self.downsample)
        note = QLabel("Downsampling is only about drawing: no result depends on it.")
        note.setWordWrap(True)
        form.addRow(note)
        return tab

    # --- Behaviour ---
    def _browse_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Fixed folder", self.fixed_dir.text())
        if folder:
            self.fixed_dir.setText(folder)

    def _sync_enabled(self):
        self.fixed_row.setEnabled(self.dir_mode.currentData() == "fixed")
        self.check_frequency.setEnabled(self.check_updates.isChecked())

    def load(self):
        """The widgets from the saved settings."""
        for key, _, setter in self._fields:
            setter(self._saved(key))
        self._sync_enabled()

    @staticmethod
    def _saved(key: str):
        default = S.DEFAULTS[key]
        if isinstance(default, bool):
            return S.get_bool(key)
        if key == S.CHECK_FREQUENCY_KEY:
            return S.get_choice(key, S.CHECK_FREQUENCIES)
        if key == S.DIR_MODE_KEY:
            return S.get_choice(key, S.DIR_MODES)
        if key == S.APPEARANCE_KEY:
            return S.get_choice(key, APPEARANCES)
        if key == S.PLOT_COLOURS_KEY:
            return S.get_choice(key, PLOT_PALETTES)
        if key == S.INITIAL_X_UNIT_KEY:
            return S.initial_x_unit()
        if key == S.FIXED_DIR_KEY:
            return str(S.app_settings().value(key, default) or "")
        limits = {S.MAX_RECENT_KEY: S.RECENT_RANGE, S.FONT_SIZE_KEY: (0, S.FONT_RANGE[1]),
                  S.LINE_WIDTH_KEY: S.LINE_WIDTH_RANGE, S.DECIMALS_KEY: S.DECIMALS_RANGE,
                  S.DOWNSAMPLE_KEY: S.DOWNSAMPLE_RANGE}[key]
        if key == S.DOWNSAMPLE_KEY and S.app_settings().value(key) is None:
            return _pyside_plot_threshold()  # the threshold of the plot, until chosen
        return S.get_number(key, limits, type(default))

    def restore_defaults(self):
        """The defaults in the widgets (saved only by OK or Apply)."""
        for key, _, setter in self._fields:
            setter(_pyside_plot_threshold() if key == S.DOWNSAMPLE_KEY else S.DEFAULTS[key])
        self._sync_enabled()

    def apply(self):
        """Saves the settings and updates the window."""
        settings = S.app_settings()
        values = {key: getter() for key, getter, _ in self._fields}
        for key, value in values.items():
            if key not in (S.APPEARANCE_KEY, S.PLOT_COLOURS_KEY):
                settings.setValue(key, value)
        settings.sync()
        # Theme and plot colours: through the View menu, which stays in sync
        self.main_window._set_appearance(values[S.APPEARANCE_KEY])
        self.main_window._set_plot_colours(values[S.PLOT_COLOURS_KEY])
        self.main_window._apply_settings()

    def accept(self):
        self.apply()
        super().accept()


def _pyside_plot_threshold() -> int:
    from .pyside_plot import DECIMATION_THRESHOLD
    return DECIMATION_THRESHOLD
