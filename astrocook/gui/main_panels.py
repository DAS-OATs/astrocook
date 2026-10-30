"""The panels of the main window: the floating cards (the session list on
the left, the plot controls on the right) with their tabs, opened and
closed with an animation; the plot controls, each redrawing the plot (plot
elements and their swatches, view, axis limits, redshift cursor, the
buttons of the continuum editor); the welcome view; the state of the
window with and without a session.

Moved out of ``main_window.py`` in phase 5 (split of the GUI), as a mixin
of ``MainWindowV2``.
"""
import logging
import os
from typing import Optional

import qtawesome as qta
from PySide6.QtCore import (QEasingCurve, QLocale, QParallelAnimationGroup, QPoint,
                            QPropertyAnimation, QRect, QSettings, QSize, Qt)
from PySide6.QtGui import QColor, QDoubleValidator, QIcon, QPainter, QPen, QPixmap
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QCheckBox, QComboBox,
                               QFormLayout, QGridLayout, QHBoxLayout, QLabel, QLineEdit,
                               QListView, QMessageBox, QPushButton, QSizePolicy, QSlider,
                               QStyle, QVBoxLayout, QWidget)

from astrocook.core.session import SessionV2
from astrocook.core.utils import resource_path

from .floating_cards import CardFader, CollapsibleSection, WrappedLabel
from .pyside_plot import SpectrumPlotWidget

# --- Constants for Sidebar Widths ---
LEFT_SIDEBAR_WIDTH = 250
RIGHT_SIDEBAR_WIDTH = 250
ANIMATION_DURATION = 150 # ** Speed up animation **
BUTTON_WIDTH = 28
BUTTON_HEIGHT = 28

# --- Floating cards (session list on the left, plot controls on the right) ---
# The panels float over the plot instead of shrinking it, so that the aspect
# of the spectrum does not change (Guido's choice). To cover as little as
# possible, they are only as tall as their contents and they fade out while
# the mouse is on the plot, so the effect of a change can be seen through.
# Cards and tabs start right below the menu bar, so that a closed tab sits in
# the corner of the window, above the axes, and does not cover the plot
CARD_MARGIN = 2        # Gap between the menu bar and the cards (px)
CARD_TAB_OFFSET = 0    # Collapse tab, from the top of its card (px)
SWATCH_SIZE = QSize(24, 12)


class MainPanelsMixin:
    """The panels of ``MainWindowV2``.

    Uses the window's ``active_history``, ``plot_viewer``,
    ``central_stack``, the animation groups of the sidebars, the menu
    actions (enabled by ``_update_ui_state``), the session list's model
    and slots (``_on_session_switched``, ``_on_session_list_context_menu``,
    ``_on_session_name_changed``) and the continuum editor's slots.
    """

    def resizeEvent(self, event):
        """Handle window resizing to reposition floating elements."""
        super().resizeEvent(event)
        self._reposition_floating_widgets()

    def _make_card_section(self, title: str, expanded: bool) -> CollapsibleSection:
        """A collapsible section of the plot controls; its state is remembered."""
        key = f"plot_controls/expanded/{title}"
        settings = QSettings("Astrocook", "AstrocookV2")
        section = CollapsibleSection(title, settings.value(key, expanded, type=bool))
        section.toggled.connect(
            lambda state: QSettings("Astrocook", "AstrocookV2").setValue(key, state))
        return section

    def _session_list_height(self) -> int:
        """Height of the session card: all the rows, at least one."""
        view = self.session_list_view
        rows = max(1, self.session_model.rowCount())
        row_height = view.sizeHintForRow(0) if self.session_model.rowCount() else 0
        if row_height <= 0:
            row_height = view.fontMetrics().height() + 20  # item padding 10 + 10
        return rows * row_height + 2 * view.frameWidth() + 8 + 8  # + layout margins

    def _reposition_floating_widgets(self):
        """Positions sidebars and buttons based on current window size and visibility."""
        window_height = self.height()
        window_width = self.width()
        menubar_height = self.menuBar().height() if self.menuBar() else 0
        content_y_start = menubar_height
        content_height = window_height - content_y_start

        # --- Left Sidebar (card) & Button ---
        left_width = self.left_sidebar_widget.width()
        # The tab is a handle just outside the card's edge (on the plot)
        button_x = left_width
        card_y = content_y_start + CARD_MARGIN
        left_height = min(self._session_list_height(), content_height - 2 * CARD_MARGIN)
        # The tab stays at the top of the card: it does not move when the
        # card grows or shrinks
        left_button_y = card_y + CARD_TAB_OFFSET
        self.session_collapse_button.setGeometry(button_x, left_button_y, BUTTON_WIDTH, BUTTON_HEIGHT)
        self.left_sidebar_widget.setGeometry(0, card_y, left_width, left_height)

        # --- Right Sidebar (card) & Button ---
        right_width = self.right_sidebar_widget.width()
        sidebar_x = window_width - right_width
        card_y = content_y_start + CARD_MARGIN
        layout = self.right_sidebar_widget.layout()
        # Labels that wrap (the hints of the continuum editor) need their
        # height at the width of the card, which sizeHint() does not give
        natural = (layout.heightForWidth(RIGHT_SIDEBAR_WIDTH) if layout.hasHeightForWidth()
                   else layout.sizeHint().height())
        card_height = min(natural, content_height - 2 * CARD_MARGIN)
        button_x = sidebar_x - BUTTON_WIDTH
        right_button_y = card_y + CARD_TAB_OFFSET
        self.plot_controls_collapse_button.setGeometry(button_x, right_button_y, BUTTON_WIDTH, BUTTON_HEIGHT)
        self.right_sidebar_widget.setGeometry(sidebar_x, card_y, right_width, card_height)
        

        self.central_stack.lower()

        # 1. Raise the sidebars (middle layer)
        self.left_sidebar_widget.raise_()
        self.right_sidebar_widget.raise_()
        
        # 2. Raise the buttons (top layer)
        self.session_collapse_button.raise_()
        self.plot_controls_collapse_button.raise_()

    def _setup_plot_view(self, session_for_plot: SessionV2):
        self.plot_viewer = SpectrumPlotWidget(session_for_plot, self)
        self.central_stack.addWidget(self.plot_viewer)

    def _setup_empty_view(self):
        empty_widget = QWidget()
        layout = QVBoxLayout(empty_widget)

        # 1. Add some stretching at the top so the content stays centered
        layout.addStretch()
        
        # 2. Add a welcome message
        label = QLabel("Welcome to")
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setStyleSheet("font-size: 36pt; color: #296bff;")
        layout.addWidget(label)

        layout.addSpacing(20)

        # 3. Add the Logo
        try:
            icon_path = resource_path(os.path.join("assets", "logo_3d_HR.png"))
            if os.path.exists(icon_path):
                logo_lbl = QLabel()
                # Scale the image, keeping it smooth and proportioned
                pixmap = QPixmap(icon_path).scaled(400, 400, Qt.KeepAspectRatio, Qt.SmoothTransformation)
                logo_lbl.setPixmap(pixmap)
                logo_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
                layout.addWidget(logo_lbl)
        except Exception as e:
            logging.warning(f"Could not load logo for empty view: {e}")

        # 4. Add instructions
        layout.addSpacing(40)

        label = QLabel("Use 'File > Open Session...'")
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setStyleSheet("font-size: 24pt; color: #ff7f0e;")

        layout.addWidget(label)

        # 5. Add stretching at the bottom
        layout.addStretch()

        self.central_stack.addWidget(empty_widget)

    def _setup_left_sidebar(self):
        """Creates the left sidebar widget as a child of the main window."""
        self.left_sidebar_widget = QWidget(self) # ** Parent is main window **
        sidebar_layout = QVBoxLayout(self.left_sidebar_widget)
        # Inset the list, so that the selected row keeps clear of the
        # rounded corners of the card
        sidebar_layout.setContentsMargins(4, 4, 4, 4); sidebar_layout.setSpacing(0)

        self.session_list_view = QListView()
        # ... (setup model, font, connection) ...
        self.session_list_view.setModel(self.session_model)
        font = self.session_list_view.font(); font.setPointSize(12); self.session_list_view.setFont(font)
        
        # ENABLE MULTI-SELECTION
        self.session_list_view.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        
        # 1. Set the context menu policy
        self.session_list_view.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        
        # 2. Connect the signal
        self.session_list_view.customContextMenuRequested.connect(self._on_session_list_context_menu)

        self.session_list_view.clicked.connect(self._on_session_switched)

        # 3a. Enable in-place editing
        self.session_list_view.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked | 
            QAbstractItemView.EditTrigger.EditKeyPressed
        )
        
        # 3b. Connect the signal that fires *after* editing is done
        self.session_model.dataChanged.connect(self._on_session_name_changed)

        # Long names: keep the end, which tells sessions apart (e.g. _2,
        # _3), and show the full name on hover
        self.session_list_view.setTextElideMode(Qt.TextElideMode.ElideMiddle)
        self.session_list_view.setMouseTracking(True)
        self.session_list_view.entered.connect(
            lambda index: self.session_list_view.setToolTip(index.data()))
        self.session_list_view.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.session_list_view.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        # The card is as tall as the list: resize it when sessions come and go
        self.session_model.rowsInserted.connect(self._reposition_floating_widgets)
        self.session_model.rowsRemoved.connect(self._reposition_floating_widgets)
        self.session_model.modelReset.connect(self._reposition_floating_widgets)

        sidebar_layout.addWidget(self.session_list_view)
        self.left_sidebar_widget.setObjectName("SessionContainer")
        self.session_list_view.setObjectName("SessionListView")
        self._sessions_fader = CardFader(self.left_sidebar_widget)

        # Initial geometry set later, start hidden/collapsed
        self.left_sidebar_widget.resize(0, self.height())
        self.left_sidebar_widget.setVisible(False)

    def _setup_right_sidebar(self):
        """Creates the right sidebar widget with plot toggles and cursor controls."""
        self.right_sidebar_widget = QWidget(self)
        sidebar_layout = QVBoxLayout(self.right_sidebar_widget)
        sidebar_layout.setContentsMargins(15, 10, 15, 10)
        sidebar_layout.setSpacing(10)

        # --- Helper: Expand Field ---
        def _expand_field(widget):
            widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            return widget

        # --- Plot Element Toggles ---
        # (This section remains unchanged)
        plot_toggles_layout = QVBoxLayout()
        plot_toggles_layout.setSpacing(4)

        # The spectrum is always drawn: a swatch with no checkbox, indented
        # to line up with the swatches of the checkboxes below
        spectrum_row = QHBoxLayout()
        indicator_w = self.style().pixelMetric(QStyle.PixelMetric.PM_IndicatorWidth)
        spectrum_row.setContentsMargins(15 + indicator_w + 5, 0, 0, 0)
        spectrum_row.setSpacing(4)
        self.spectrum_swatch = QLabel()
        self.spectrum_swatch.setPixmap(self._blank_swatch())
        spectrum_row.addWidget(self.spectrum_swatch)
        spectrum_row.addWidget(QLabel("Spectrum"))
        spectrum_row.addStretch(1)
        plot_toggles_layout.addLayout(spectrum_row)

        self.error_checkbox = QCheckBox("1-sigma error"); self.error_checkbox.setChecked(True)
        self.error_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.error_checkbox)

        # --- Error Style Toggle ---
        self.error_line_checkbox = QCheckBox("Show error as line")
        self.error_line_checkbox.setToolTip("If checked, plots the error value (dy) as a line instead of a band.")
        self.error_line_checkbox.setChecked(True) # Default to Band
        # Indent to show hierarchy
        self.error_line_checkbox.setStyleSheet("margin-left: 20px;")
        
        # Connect to main replot
        self.error_line_checkbox.stateChanged.connect(self._trigger_replot)
        # Also connect to System Inspector update if open
        self.error_line_checkbox.stateChanged.connect(
            lambda: self.system_inspector.vel_plot._update_plot() 
            if hasattr(self, 'system_inspector') and self.system_inspector and self.system_inspector.isVisible() 
            else None
        )
        
        plot_toggles_layout.addWidget(self.error_line_checkbox)

        self.continuum_checkbox = QCheckBox("Normalization"); self.continuum_checkbox.setChecked(True)
        self.continuum_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.continuum_checkbox)

        self.model_checkbox = QCheckBox("Absorption model"); self.model_checkbox.setChecked(True)
        self.model_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.model_checkbox)

        self.systems_checkbox = QCheckBox("Systems"); self.systems_checkbox.setChecked(True)
        self.systems_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.systems_checkbox)

        # Aux Column
        form_layout_aux = QFormLayout()
        form_layout_aux.setSpacing(5)
        self.aux_col_combo = QComboBox()
        self.aux_col_combo.setToolTip("Select an auxiliary column to plot")
        self.aux_col_combo.currentTextChanged.connect(self._trigger_replot)
        # Its swatch, as for the elements above (the plot has no legend, F86)
        self.aux_swatch = QLabel()
        self.aux_swatch.setPixmap(self._blank_swatch())
        aux_row = QHBoxLayout()
        aux_row.setSpacing(4)
        aux_row.addWidget(self.aux_swatch)
        aux_row.addWidget(self.aux_col_combo, 1)
        form_layout_aux.addRow("Aux. Column:", aux_row)
        plot_toggles_layout.addLayout(form_layout_aux)

        self.strong_lines_checkbox = QCheckBox("Strong emission lines only")
        self.strong_lines_checkbox.setToolTip("Show only major emission lines when z_em is active.")
        self.strong_lines_checkbox.setChecked(True)
        self.strong_lines_checkbox.stateChanged.connect(self._trigger_replot)
        plot_toggles_layout.addWidget(self.strong_lines_checkbox)

        # Its swatches show the styles of the plot, which has no legend (F86)
        self.plot_elements_section = self._make_card_section("Plot Elements", expanded=True)
        self.plot_elements_section.setContentLayout(plot_toggles_layout)
        sidebar_layout.addWidget(self.plot_elements_section)

        self.isomag_checkbox = QCheckBox("Show Iso-Mag Grid")
        self.isomag_checkbox.setToolTip("Show lines of constant AB magnitude")
        self.isomag_checkbox.toggled.connect(lambda b: self.plot_viewer.toggle_isomag_grid(b))
        self.isomag_checkbox.setObjectName("PlotControlCheckbox")
        plot_toggles_layout.addWidget(self.isomag_checkbox)
    
        # --- ** Axis & View Controls ** ---
        view_layout = QVBoxLayout()
        view_layout.setSpacing(4)

        # X-Axis Unit
        form_layout_x = QFormLayout()
        self.x_unit_combo = QComboBox()
        self.x_unit_options = ["nm", "Angstrom", "micron"]
        self.x_unit_combo.addItems(self.x_unit_options)
        self.x_unit_combo.setToolTip("Change λ display units (data remains nm).")
        form_layout_x.addRow("λ Unit:", self.x_unit_combo)
        view_layout.addLayout(form_layout_x)
        
        # X/Y Toggles
        self.norm_y_checkbox = QCheckBox("Normalize F")
        self.norm_y_checkbox.setToolTip("Plot F / Continuum and set F limits.")
        self.norm_y_checkbox.toggled.connect(self._on_view_toggle)
        view_layout.addWidget(self.norm_y_checkbox)

        self.snr_checkbox = QCheckBox("Show SNR (F / error)")
        self.snr_checkbox.setToolTip("Plot F / Error Column.")
        self.snr_checkbox.toggled.connect(self._on_view_toggle)
        view_layout.addWidget(self.snr_checkbox)
        
        form_layout_snr = QFormLayout()
        form_layout_snr.setSpacing(5)
        self.snr_col_combo = QComboBox()
        self.snr_col_combo.setToolTip("Select column to use as error for SNR")
        self.snr_col_combo.currentTextChanged.connect(self._trigger_replot)
        form_layout_snr.addRow("  Error Col.:", self.snr_col_combo)
        view_layout.addLayout(form_layout_snr)

        self.log_x_checkbox = QCheckBox("Logarithmic λ")
        view_layout.addWidget(self.log_x_checkbox)
        # Not in the panel before (findings register F78)
        self.log_y_checkbox = QCheckBox("Logarithmic F")
        self.log_y_checkbox.setToolTip(
            "Logarithmic scale of the flux. Pixels with F \u2264 0 are drawn at "
            "the bottom of the plot; the normalized flux starts from 0.01.")
        view_layout.addWidget(self.log_y_checkbox)
        
        view_section = self._make_card_section("View Controls", expanded=True)
        view_section.setContentLayout(view_layout)
        sidebar_layout.addWidget(view_section)

        self.x_unit_combo.currentTextChanged.connect(self._trigger_replot)
        self.log_x_checkbox.stateChanged.connect(self._trigger_replot)
        self.log_y_checkbox.stateChanged.connect(self._trigger_replot)

        # -- SECTION 1: Axis Limits --
        limits_layout = QGridLayout()
        limits_layout.setSpacing(5)

        validator = QDoubleValidator()
        validator.setLocale(QLocale.C)
        validator.setNotation(QDoubleValidator.ScientificNotation)

        self.xmin_input = _expand_field(QLineEdit("0.0")); self.xmin_input.setValidator(validator)
        self.xmax_input = _expand_field(QLineEdit("1.0")); self.xmax_input.setValidator(validator)
        self.ymin_input = _expand_field(QLineEdit("0.0")); self.ymin_input.setValidator(validator)
        self.ymax_input = _expand_field(QLineEdit("1.0")); self.ymax_input.setValidator(validator)

        self.xmin_input.editingFinished.connect(self._on_set_custom_limits)
        self.xmax_input.editingFinished.connect(self._on_set_custom_limits)
        self.ymin_input.editingFinished.connect(self._on_set_custom_limits)
        self.ymax_input.editingFinished.connect(self._on_set_custom_limits)
        
        # λ and F on same rows
        lbl_lambda = QLabel("λ:")
        lambda_layout = QHBoxLayout()
        lambda_layout.setSpacing(3)
        lambda_layout.addWidget(self.xmin_input)
        lambda_layout.addWidget(self.xmax_input)
        
        limits_layout.addWidget(lbl_lambda, 0, 0)
        limits_layout.addLayout(lambda_layout, 0, 1)
        
        self.lbl_flux = QLabel("F:")
        flux_layout = QHBoxLayout()
        flux_layout.setSpacing(3)
        flux_layout.addWidget(self.ymin_input)
        flux_layout.addWidget(self.ymax_input)
        
        limits_layout.addWidget(self.lbl_flux, 1, 0)
        limits_layout.addLayout(flux_layout, 1, 1)

        limits_section = self._make_card_section("Axis Limits", expanded=True)
        limits_section.setContentLayout(limits_layout)
        sidebar_layout.addWidget(limits_section)

        # -- SECTION 2: Redshift Cursor --
        cursor_group = self._make_card_section("Redshift Cursor", expanded=False)
        cursor_layout = QGridLayout()
        cursor_layout.setSpacing(5)

        self.cursor_series_input = _expand_field(QLineEdit("Ly_a"))
        self.cursor_z_input = _expand_field(QLineEdit("0.0"))
        z_validator = QDoubleValidator()
        z_validator.setLocale(QLocale.C) 
        z_validator.setNotation(QDoubleValidator.StandardNotation)
        self.cursor_z_input.setValidator(z_validator)

        cursor_layout.addWidget(QLabel("Trans.:"), 0, 0)
        cursor_layout.addWidget(self.cursor_series_input, 0, 1)
        cursor_layout.addWidget(QLabel("Redshift:"), 1, 0)
        cursor_layout.addWidget(self.cursor_z_input, 1, 1)

        self.cursor_show_checkbox = QCheckBox("Show Cursor Lines")
        self.cursor_show_checkbox.setChecked(False)
        cursor_layout.addWidget(self.cursor_show_checkbox, 2, 0, 1, 2)

        self.cursor_series_input.editingFinished.connect(self._update_cursor_and_replot)
        self.cursor_z_input.editingFinished.connect(self._update_cursor_and_replot)
        self.cursor_show_checkbox.stateChanged.connect(self._on_cursor_checkbox_toggled)
        self.cursor_series_input.returnPressed.connect(lambda: self.cursor_show_checkbox.setChecked(True))

        cursor_group.setContentLayout(cursor_layout)
        sidebar_layout.addWidget(cursor_group)

        # -- SECTION 4: Edit Continuum (Collapsible, hidden by default) --
        continuum_group = self._make_card_section("Edit Continuum", expanded=False)
        continuum_layout = QGridLayout()
        continuum_layout.setSpacing(5)

        # Buttons first
        button_layout = QHBoxLayout()
        button_layout.setSpacing(5)

        self.edit_cont_button = QPushButton("Start")
        self.edit_cont_button.clicked.connect(self._toggle_continuum_editor)
        self.edit_cont_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        button_layout.addWidget(self.edit_cont_button)

        self.cancel_cont_button = QPushButton("Cancel")
        self.cancel_cont_button.setEnabled(False)
        self.cancel_cont_button.clicked.connect(self._on_cancel_continuum_clicked)
        self.cancel_cont_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        button_layout.addWidget(self.cancel_cont_button)

        self.reset_cont_button = QPushButton("Reset")
        self.reset_cont_button.setEnabled(False)
        self.reset_cont_button.setToolTip("Discard changes and reset to default spacing.")
        self.reset_cont_button.clicked.connect(self._on_reset_continuum_clicked)
        self.reset_cont_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        button_layout.addWidget(self.reset_cont_button)

        continuum_layout.addLayout(button_layout, 0, 0, 1, 2)

        # Slider
        self.stride_slider = QSlider(Qt.Horizontal)
        self.stride_slider.setMinimum(0)
        self.stride_slider.setMaximum(100)
        self.stride_slider.setEnabled(False)
        self.stride_slider.valueChanged.connect(self._on_stride_change)
        self.stride_slider.sliderPressed.connect(self._on_stride_slider_pressed)
        self.stride_slider.sliderReleased.connect(self._on_stride_slider_released)
        self.stride_slider.setToolTip("Adjust density. Drag to smooth/roughen.")
        self._set_slider_from_stride(1500)
        
        continuum_layout.addWidget(self.stride_slider, 1, 0, 1, 2)

        # Spacing
        self.lbl_spacing = QLabel("Spacing:")
        continuum_layout.addWidget(self.lbl_spacing, 2, 0)
        
        self.stride_entry = QLineEdit()
        stride_validator = QDoubleValidator(0.0, 20000.0, 0)
        stride_validator.setLocale(QLocale.C)
        stride_validator.setNotation(QDoubleValidator.StandardNotation)
        self.stride_entry.setValidator(stride_validator)
        self.stride_entry.setAlignment(Qt.AlignLeft)
        self.stride_entry.setEnabled(False)
        self.stride_entry.setToolTip("Enter spacing in km/s")
        self.stride_entry.editingFinished.connect(self._on_stride_text_changed)
        
        continuum_layout.addWidget(self.stride_entry, 2, 1)

        # --- Edit Hints Label (Row 3, spanning 2 columns) ---
        # 1. Create a transparent container for the label
        self.edit_hints_container = QWidget()
        hints_layout = QVBoxLayout(self.edit_hints_container)
        hints_layout.setContentsMargins(0, 0, 0, 0)
        
        # 2. Create the label itself
        self.edit_hints_label = WrappedLabel(
            "<div style='font-size: 11pt;'>"
            "<p style='margin: 0px 0px 6px 0px;'><b>Controls:</b></p>"
            "<p style='margin: 0px 0px 2px 0px;'>• <b>Shift + Click:</b> Add knot</p>"
            "<p style='margin: 0px 0px 2px 0px;'>• <b>Right-Click:</b> Remove knot</p>"
            "<p style='margin: 0px 0px 0px 0px;'>• <b>Drag:</b> Move knot<br>"
            "&nbsp;&nbsp;<i>(Stay within neighboring knots)</i></p>"
            "</div>"
        )
        self.edit_hints_label.setWordWrap(True)
        self.edit_hints_label.setVisible(False) # Hidden by default
        
        # 3. Put the label in the container
        hints_layout.addWidget(self.edit_hints_label)

        # 4. Add the CONTAINER to the main layout
        continuum_layout.addWidget(self.edit_hints_container, 3, 0, 1, 2)

        continuum_group.setContentLayout(continuum_layout)
        sidebar_layout.addWidget(continuum_group)

        self.right_sidebar_widget.setObjectName("PlotControlsContainer")
        self._plot_controls_fader = CardFader(self.right_sidebar_widget,
                                              on_layout=self._reposition_floating_widgets)
        for checkbox in self._swatch_checkboxes():
            checkbox.setIcon(QIcon(self._blank_swatch()))
            checkbox.setIconSize(SWATCH_SIZE)
        self.error_checkbox.setObjectName("PlotControlCheckbox")
        self.error_line_checkbox.setObjectName("PlotControlCheckbox")
        self.continuum_checkbox.setObjectName("PlotControlCheckbox")
        self.model_checkbox.setObjectName("PlotControlCheckbox")
        self.systems_checkbox.setObjectName("PlotControlCheckbox")
        self.aux_col_combo.setObjectName("AuxColumnCombo")
        self.strong_lines_checkbox.setObjectName("PlotControlCheckbox")
        self.isomag_checkbox.setObjectName("PlotControlCheckbox")
        self.x_unit_combo.setObjectName("XUnitCombo")
        self.norm_y_checkbox.setObjectName("PlotControlCheckbox")
        self.snr_checkbox.setObjectName("PlotControlCheckbox")
        self.log_x_checkbox.setObjectName("PlotControlCheckbox")
        self.log_y_checkbox.setObjectName("PlotControlCheckbox")
        self.cursor_show_checkbox.setObjectName("PlotControlCheckbox")
        
        self._populate_snr_combo()

        self.right_sidebar_widget.resize(0, self.height())
        self.right_sidebar_widget.setVisible(False)

    def _setup_collapse_buttons(self):
        """Creates and positions the collapse buttons as children of main window."""
        # Left Button
        self.session_collapse_button = QPushButton(self) # ** Parent is main window **
        self.session_collapse_button.setToolTip("Collapse/Expand Session Panel")
        self.session_collapse_button.setFixedSize(QSize(BUTTON_WIDTH, BUTTON_HEIGHT))
        self.session_collapse_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed) # Fill height
        self.session_collapse_button.clicked.connect(self._toggle_left_sidebar)
        self.session_collapse_button.setObjectName("CollapseButton")
        self.session_collapse_button.setProperty("side", "left")
        self.session_collapse_button.setVisible(False) # Start hidden

        # Right Button
        self.plot_controls_collapse_button = QPushButton(self) # ** Parent is main window **
        self.plot_controls_collapse_button.setToolTip("Collapse/Expand Plot Controls")
        self.plot_controls_collapse_button.setFixedSize(QSize(BUTTON_WIDTH, BUTTON_HEIGHT))
        self.plot_controls_collapse_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed) # Fill height
        self.plot_controls_collapse_button.clicked.connect(self._toggle_right_sidebar)
        self.plot_controls_collapse_button.setObjectName("CollapseButton")
        self.plot_controls_collapse_button.setProperty("side", "right")
        self.plot_controls_collapse_button.setVisible(False) # Start hidden

    def _trigger_replot(self):
        """Slot called when plot toggles or cursor changes."""
        if self.plot_viewer:
            # Check if we are currently in manual edit mode
            if getattr(self.plot_viewer, '_edit_mode_active', False):
                # Tell the plot viewer to refresh using the live buffer
                # instead of the session's static 'cont' column.
                self.plot_viewer.refresh_edit_layer()
            else:
                # Standard behavior: draw from the session state
                session_state = self.active_history.current_state if self.active_history else None
                self.plot_viewer.plot_spectrum(session_state=session_state)

    def _update_cursor_and_replot(self):
        """Slot called when cursor series or redshift changes."""
        # Optional: Add validation for series input here if needed
        # (e.g., check if series exists in xem_d)
        if self.cursor_show_checkbox.isChecked(): # Only replot if cursor is visible
            self._trigger_replot()

    def _on_cursor_checkbox_toggled(self, checked):
        """Intercepts the redshift cursor toggle to prevent conflicts with continuum editing."""
        if checked and hasattr(self, 'plot_viewer') and getattr(self.plot_viewer, '_edit_mode_active', False):
            # Prompt the user to save or discard continuum edits
            ret = self._show_custom_message(
                title="Continuum Edit Active",
                header="Save Continuum Edits?",
                text="You must close the Continuum Editor to use the Redshift Cursor.<br><br>Do you want to save your current continuum edits?",
                buttons=QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
                default_btn=QMessageBox.StandardButton.Save
            )
            
            if ret == QMessageBox.StandardButton.Save:
                self._toggle_continuum_editor() # This saves and turns off edit mode
                self._trigger_replot()
            elif ret == QMessageBox.StandardButton.Discard:
                # Cancel edit without saving
                self.plot_viewer.stop_continuum_edit(save=False)
                self._reset_continuum_ui()
                self._trigger_replot()
            else: # Cancel
                # Revert the checkbox visually without triggering signals
                self.cursor_show_checkbox.blockSignals(True)
                self.cursor_show_checkbox.setChecked(False)
                self.cursor_show_checkbox.blockSignals(False)
                return
        else:
            # Standard behavior
            self._trigger_replot()

    def _on_set_custom_limits(self):
        """ Reads values from limit boxes and applies them to the plot. """
        if not self.plot_viewer:
            return
            
        try:
            # Read values as-typed (they are in display units)
            display_xmin = float(self.xmin_input.text())
            display_xmax = float(self.xmax_input.text())
            ymin = float(self.ymin_input.text())
            ymax = float(self.ymax_input.text())
        except ValueError:
            logging.warning("Invalid custom limit values entered.")
            return # Ignore if values are not valid floats
            
        selected_unit = self.x_unit_combo.currentText()
        
        if selected_unit == 'Angstrom':
            xmin_nm = display_xmin / 10.0
            xmax_nm = display_xmax / 10.0
        elif selected_unit == 'micron':
            xmin_nm = display_xmin * 1000.0
            xmax_nm = display_xmax * 1000.0
        else: # 'nm'
            xmin_nm = display_xmin
            xmax_nm = display_xmax

        if xmax_nm <= xmin_nm or ymax <= ymin:
            logging.warning("Invalid custom limit range (max <= min).")
            return
        if self.log_y_checkbox.isChecked() and ymin <= 0:
            # A log scale needs a positive limit: Matplotlib would ignore it
            # with a warning (F78)
            ymin = self.plot_viewer.canvas.axes.get_ylim()[0]
            logging.warning("With the logarithmic F scale the lower limit of F must "
                            f"be positive: kept at {ymin:.4g}.")
            if ymax <= ymin:
                return

        # Set limits on the Matplotlib axes
        # This will automatically trigger the 'on_lim_changed' callback
        # in pyside_plot.py, which schedules a new plot_spectrum draw.
        self.plot_viewer.canvas.axes.set_xlim(xmin_nm, xmax_nm)
        self.plot_viewer.canvas.axes.set_ylim(ymin, ymax)

    def _update_limit_boxes_from_plot(self):
        """ Updates the limit text boxes with the plot's current limits. """
        if not self.plot_viewer:
            return
            
        try:
            xlim_nm = self.plot_viewer.canvas.axes.get_xlim() # These are always in nm
            ylim_plot = self.plot_viewer.canvas.axes.get_ylim()
            
            # --- NEW: Convert nm limits to selected display unit ---
            selected_unit = self.x_unit_combo.currentText()

            if selected_unit == 'Angstrom':
                display_xlim = (xlim_nm[0] * 10.0, xlim_nm[1] * 10.0)
            elif selected_unit == 'micron':
                display_xlim = (xlim_nm[0] / 1000.0, xlim_nm[1] / 1000.0)
            else: # 'nm'
                display_xlim = xlim_nm
            # --- END NEW ---
            
            # Block signals to prevent a feedback loop
            self.xmin_input.blockSignals(True)
            self.xmax_input.blockSignals(True)
            self.ymin_input.blockSignals(True)
            self.ymax_input.blockSignals(True)
            
            # Set text using the converted display values
            self.xmin_input.setText(f"{display_xlim[0]:.4g}")
            self.xmax_input.setText(f"{display_xlim[1]:.4g}")
            self.ymin_input.setText(f"{ylim_plot[0]:.4g}")
            self.ymax_input.setText(f"{ylim_plot[1]:.4g}")
            
        except Exception as e:
            logging.warning(f"Failed to update limit boxes: {e}")
        finally:
            # Always unblock signals
            self.xmin_input.blockSignals(False)
            self.xmax_input.blockSignals(False)
            self.ymin_input.blockSignals(False)
            self.ymax_input.blockSignals(False)

    def _force_restack_floating_widgets(self):
        """
        Slot to be called by QTimer.
        This is the programmatic equivalent of the user manually
        resizing the window, which we know fixes the bug.
        """
        logging.debug("QTimer: Forcing native resize event to fix widget stacking.")
        try:
            # 1. Get the current size
            current_size = self.size()
            
            # 2. "Jiggle" the size: resize by +1 pixel
            # This is invisible but forces the native resizeEvent to run
            self.resize(current_size.width() + 1, current_size.height())
            
            # 3. Immediately resize back to the original size
            self.resize(current_size)
            
            # The resizeEvent itself will call _reposition_floating_widgets,
            # which will handle the final Z-stacking.
            
        except Exception as e:
            logging.warning(f"Failed to force resize: {e}")

    def _toggle_left_sidebar(self):
        """Animates the left sidebar overlay and its button together."""
        is_currently_visible = self.left_sidebar_widget.width() > 0 # Use width check
        start_sidebar_geom = QRect(self.left_sidebar_widget.geometry()) # Get current geometry
        end_sidebar_geom = QRect(start_sidebar_geom) # Copy for end state
        start_button_x = self.session_collapse_button.x()
        end_button_x = 0 # Default button X when closed

        if is_currently_visible: # Closing
            end_sidebar_geom.setWidth(0) # Animate width to 0
            # Button X already calculated (end_button_x = 0)
            target_visible = False
        else: # Opening
            # Start geom should reflect the hidden state (width 0)
            start_sidebar_geom.setWidth(0)
            end_sidebar_geom.setWidth(LEFT_SIDEBAR_WIDTH) # Animate width open
            end_button_x = LEFT_SIDEBAR_WIDTH # Button moves out with the card's edge
            target_visible = True
            self.left_sidebar_widget.setVisible(True) # Show before animating

        # --- Animation Group ---
        # Stop previous animation if running
        if self.left_animation_group and self.left_animation_group.state() == QParallelAnimationGroup.Running:
            self.left_animation_group.stop()
        self.left_animation_group = QParallelAnimationGroup(self)

        # ** Animate sidebar GEOMETRY **
        sidebar_anim = QPropertyAnimation(self.left_sidebar_widget, b"geometry")
        sidebar_anim.setDuration(ANIMATION_DURATION)
        sidebar_anim.setStartValue(start_sidebar_geom) # Start geometry
        sidebar_anim.setEndValue(end_sidebar_geom)     # End geometry
        sidebar_anim.setEasingCurve(QEasingCurve.InOutQuad)
        self.left_animation_group.addAnimation(sidebar_anim)

        # Animate button X position
        button_anim = QPropertyAnimation(self.session_collapse_button, b"pos")
        button_anim.setDuration(ANIMATION_DURATION)
        # Use current Y, animate X
        button_anim.setStartValue(QPoint(start_button_x, self.session_collapse_button.y()))
        button_anim.setEndValue(QPoint(end_button_x, self.session_collapse_button.y()))
        button_anim.setEasingCurve(QEasingCurve.InOutQuad)
        self.left_animation_group.addAnimation(button_anim)

        # Hide sidebar widget AFTER closing animation finishes
        if not target_visible:
            # Ensure lambda captures the correct widget state if needed
            self.left_animation_group.finished.connect(lambda: self.left_sidebar_widget.setVisible(False))
        # Opaque when it opens, and its tab opaque when it closes
        self._sessions_fader.timer.stop()
        self._set_sessions_faded(False, animate=False)

        self.left_animation_group.start()
        self._update_sidebar_button_icon(self.session_collapse_button, target_visible, is_left=True)

    def _toggle_right_sidebar(self):
        # ... (get start/end width, window_width, start/end sidebar X, start/end button X) ...
        start_width = self.right_sidebar_widget.width(); end_width = 0 if start_width > 0 else RIGHT_SIDEBAR_WIDTH
        window_width = self.width()
        start_sidebar_x = self.right_sidebar_widget.x(); end_sidebar_x = window_width if end_width == 0 else window_width - end_width
        start_button_x = self.plot_controls_collapse_button.x(); end_button_x = window_width - end_width - BUTTON_WIDTH
        target_visible = end_width > 0

        # ** Calculate the constant button Y **
        button_y = self.plot_controls_collapse_button.y() # Keep current Y

        if target_visible: self.right_sidebar_widget.setVisible(True)
        self.right_animation_group = QParallelAnimationGroup(self)
        # Sidebar Animation (geometry)
        sidebar_geom_anim = QPropertyAnimation(self.right_sidebar_widget, b"geometry"); sidebar_geom_anim.setDuration(ANIMATION_DURATION)
        sidebar_geom_anim.setStartValue(QRect(start_sidebar_x, self.right_sidebar_widget.y(), start_width, self.right_sidebar_widget.height()))
        sidebar_geom_anim.setEndValue(QRect(end_sidebar_x, self.right_sidebar_widget.y(), end_width, self.right_sidebar_widget.height()))
        sidebar_geom_anim.setEasingCurve(QEasingCurve.InOutQuad)
        self.right_animation_group.addAnimation(sidebar_geom_anim)
        # Button Animation (pos)
        button_anim = QPropertyAnimation(self.plot_controls_collapse_button, b"pos"); button_anim.setDuration(ANIMATION_DURATION)
        button_anim.setStartValue(QPoint(start_button_x, button_y)) # Start at current Y
        button_anim.setEndValue(QPoint(end_button_x, button_y))     # End at same Y
        button_anim.setEasingCurve(QEasingCurve.InOutQuad)
        self.right_animation_group.addAnimation(button_anim)

        if not target_visible: self.right_animation_group.finished.connect(lambda: self.right_sidebar_widget.setVisible(False))
        # Opaque when it opens, and its tab opaque when it closes
        self._plot_controls_fader.timer.stop()
        self._set_plot_controls_faded(False, animate=False)
        # No redraw: it moved the entries of the legend out of the card, and
        # the plot has no legend now (F86)
        self.right_animation_group.start()
        self._update_sidebar_button_icon(self.plot_controls_collapse_button, target_visible, is_left=False)

    # --- Floating cards: fading and legend swatches ---

    def _set_plot_controls_faded(self, faded: bool, animate: bool = True):
        self._plot_controls_fader.set_faded(faded, animate)

    def _set_sessions_faded(self, faded: bool, animate: bool = True):
        self._sessions_fader.set_faded(faded, animate)

    def _swatch_checkboxes(self):
        return [self.error_checkbox, self.continuum_checkbox,
                self.model_checkbox, self.systems_checkbox]

    def _blank_swatch(self) -> QPixmap:
        ratio = self.devicePixelRatioF()
        pix = QPixmap(SWATCH_SIZE * ratio)
        pix.setDevicePixelRatio(ratio)
        pix.fill(Qt.transparent)
        return pix

    def _swatch(self, artists) -> QPixmap:
        """Draws a legend swatch for one or more matplotlib artists."""
        from matplotlib.colors import to_rgba
        pix = self._blank_swatch()
        painter = QPainter(pix)
        painter.setRenderHint(QPainter.Antialiasing)
        w, h = SWATCH_SIZE.width(), SWATCH_SIZE.height()
        n = len(artists)
        for i, artist in enumerate(artists):
            if hasattr(artist, 'lines'):  # errorbar container
                coll = artist.lines[2][0] if artist.lines[2] else None
                rgba = coll.get_colors()[0] if coll is not None else (0.5, 0.5, 0.5, 1)
                linestyle, marker, lw = '-', 'None', 1.0
            elif not hasattr(artist, 'get_linestyle') or hasattr(artist, 'get_paths'):
                # A shaded area (regions, masks): a filled rectangle, more
                # opaque than in the plot to read at this size
                rgba = to_rgba(artist.get_facecolor()[0])
                color = QColor.fromRgbF(rgba[0], rgba[1], rgba[2], max(rgba[3], 0.5))
                painter.setPen(Qt.NoPen)
                painter.setBrush(color)
                painter.drawRect(int(i * w / n) + 1, 2, int(w / n) - 2, h - 4)
                continue
            else:
                rgba = to_rgba(artist.get_color(), artist.get_alpha())
                linestyle, marker = artist.get_linestyle(), artist.get_marker()
                lw = artist.get_linewidth()
            color = QColor.fromRgbF(*rgba)
            if marker == '|':
                # Ticks of the systems: one group per artist, side by side
                pen = QPen(color, 1.5)
                painter.setPen(pen)
                x0, x1 = i * w / n, (i + 1) * w / n
                for k in range(3):
                    x = x0 + (x1 - x0) * (k + 0.5) / 3
                    painter.drawLine(int(x), 1, int(x), h - 1)
            elif linestyle not in ('None', '', ' '):
                # Thicker than in the plot, to read as a colour at this size
                pen = QPen(color, max(1.5, lw * 1.5))
                if linestyle in ('--', 'dashed'):
                    pen.setStyle(Qt.DashLine)
                elif linestyle in (':', 'dotted'):
                    pen.setStyle(Qt.DotLine)
                painter.setPen(pen)
                painter.drawLine(1, h // 2, w - 1, h // 2)
            else:
                painter.setPen(Qt.NoPen)
                painter.setBrush(color)
                painter.drawEllipse(w // 2 - 3, h // 2 - 3, 6, 6)
        painter.end()
        return pix

    def update_plot_swatches(self, handles_by_label: dict, ghosts: Optional[dict] = None):
        """Shows the plotted styles next to the plot controls, and those of
        the other sessions shown (``ghosts``: row of the list -> line) next
        to their names in the list of the sessions. The plot has no legend:
        these swatches replace it (F86).

        Called by the plot after each redraw. An element that is not drawn
        although its box is checked (e.g. no continuum) gets a blank swatch;
        an unchecked one keeps its last swatch, to show what it would add.
        """
        def pick(*labels):
            return [handles_by_label[lab] for lab in labels if lab in handles_by_label]

        spectrum = pick('Spectrum', 'Spectrum (decimated)')
        self.spectrum_swatch.setPixmap(self._swatch(spectrum[:1]) if spectrum else self._blank_swatch())
        for checkbox, labels in (
                (self.error_checkbox, ('1-sigma error',)),
                (self.continuum_checkbox, ('Normalization',)),
                (self.model_checkbox, ('Absorption model',)),
                (self.systems_checkbox, ('HI Systems', 'Metal Systems'))):
            artists = pick(*labels)
            if artists:
                checkbox.setIcon(QIcon(self._swatch(artists)))
            elif checkbox.isChecked():
                checkbox.setIcon(QIcon(self._blank_swatch()))
        aux = self.aux_col_combo.currentText()
        artists = pick(aux, 'Detected Regions', f'Mask: {aux}')
        self.aux_swatch.setPixmap(self._swatch(artists[:1]) if artists
                                  else self._blank_swatch())
        # The colour of each session shown, opaque and thick: its line is
        # thin and translucent, and the swatch is drawn on the selection
        from matplotlib.lines import Line2D
        swatches = {}
        for row, line in (ghosts or {}).items():
            pix = self._swatch([Line2D([], [], color=line.get_color(), lw=3)])
            icon = QIcon(pix)
            # The same colours on a selected row (Qt tints the icon there)
            icon.addPixmap(pix, QIcon.Mode.Selected)
            swatches[row] = icon
        self.session_model.set_swatches(swatches, QIcon(self._blank_swatch()))

    def _update_sidebar_button_icon(self, button, sidebar_is_visible, is_left):
        """Updates collapse button icon based on sidebar visibility."""
        # Icon logic adjusted: Left button shows '>' when sidebar is hidden (points right to expand)
        # Right button shows '<' when sidebar is hidden (points left to expand)
        # The tabs show what they open (a list of sessions, the settings of
        # the plot), not a direction: the triangles are left to the sections
        color = QApplication.palette().color(QApplication.palette().ColorRole.ButtonText)
        if is_left:
            icon = qta.icon('mdi.menu', color=color)
            tooltip = "Collapse Session Panel" if sidebar_is_visible else "Expand Session Panel"
        else: # Right sidebar button
            icon = qta.icon('mdi.cog', color=color)
            tooltip = "Collapse Plot Controls" if sidebar_is_visible else "Expand Plot Controls"

        button.setIcon(icon); button.setToolTip(tooltip)

    def _on_view_toggle(self, checked: bool):
        """ Handles toggling the main view (Norm vs SNR vs Raw). """
        sender = self.sender()
        
        # 1. Enforce mutual exclusion
        if checked:
            if sender is self.norm_y_checkbox:
                self.snr_checkbox.blockSignals(True)
                self.snr_checkbox.setChecked(False)
                self.snr_checkbox.blockSignals(False)
            elif sender is self.snr_checkbox:
                self.norm_y_checkbox.blockSignals(True)
                self.norm_y_checkbox.setChecked(False)
                self.norm_y_checkbox.blockSignals(False)
        
        # 2. Show/hide the SNR error column selector
        #self.snr_col_combo.parentWidget().setVisible(self.snr_checkbox.isChecked())
        
        # Determine symbol
        is_norm = self.norm_y_checkbox.isChecked()
        sym = "ƒ" if is_norm else "F"
        
        # Update Labels
        if hasattr(self, 'lbl_flux'):
            self.lbl_flux.setText(f"{sym}:")

        # 3. Trigger the redraw
        self._trigger_replot()

    # Helper to populate the SNR combo with "dF" instead of "dy"
    # Call this whenever you load a session or refresh columns
    def _populate_snr_combo(self):
        self.snr_col_combo.blockSignals(True)
        current_text = self.snr_col_combo.currentText() 
        self.snr_col_combo.clear()
        
        # 1. Always add the standard error as "dF"
        self.snr_col_combo.addItem("dF") 
        
        # 2. Add filtered aux columns
        if self.active_history and self.active_history.current_state.spec:
            spec = self.active_history.current_state.spec
            
            if hasattr(spec, '_data') and hasattr(spec._data, 'aux_cols'):
                for col in spec._data.aux_cols:
                    # [CHANGE] Strict filtering: only 'running_std' or '*err*'
                    is_error_like = (col == 'running_std') or ('err' in col.lower())
                    
                    if is_error_like and col != 'dy': 
                        self.snr_col_combo.addItem(col)

        # Restore selection if it still exists
        index = self.snr_col_combo.findText(current_text)
        if index >= 0:
            self.snr_col_combo.setCurrentIndex(index)
        else:
            self.snr_col_combo.setCurrentIndex(0) # Default to dF

        self.snr_col_combo.blockSignals(False)

    def _update_ui_state(self, is_valid_session, is_startup=False):
        is_valid_session = bool(is_valid_session)

        # ... (Enable/Disable File menu actions) ...
        if hasattr(self, 'save_action'): self.save_action.setEnabled(is_valid_session)
        if hasattr(self, 'export_action'): self.export_action.setEnabled(is_valid_session)
        if hasattr(self, 'import_action'): self.import_action.setEnabled(is_valid_session)
        if hasattr(self, 'close_session_action'): self.close_session_action.setEnabled(is_valid_session)

        # ... (Enable/Disable View menu actions) ...
        if hasattr(self, 'toggle_left_action'): self.toggle_left_action.setEnabled(is_valid_session)
        if hasattr(self, 'toggle_right_action'): self.toggle_right_action.setEnabled(is_valid_session)
        if hasattr(self, 'view_system_inspector_action'): self.view_system_inspector_action.setEnabled(is_valid_session)
        if hasattr(self, 'view_log_action'): self.view_log_action.setEnabled(is_valid_session)
        if hasattr(self, 'view_script_console_action'): self.view_script_console_action.setEnabled(is_valid_session)
        has_ids = False
        if is_valid_session and self.active_history.current_state.spec:
            # Check if the key exists AND is not None
            has_ids = self.active_history.current_state.spec.meta.get('region_identifications') is not None
        if hasattr(self, 'view_identifications_action'): 
            self.view_identifications_action.setEnabled(has_ids)

        # ** Enable/Disable Recipe Actions **
        enable_recipes = is_valid_session
        # Check if actions exist before enabling/disabling
        if hasattr(self, 'set_properties_action'): self.set_properties_action.setEnabled(enable_recipes)
        if hasattr(self, 'apply_expression_action'): self.apply_expression_action.setEnabled(enable_recipes)
        if hasattr(self, 'mask_expression_action'): self.mask_expression_action.setEnabled(enable_recipes)
        if hasattr(self, 'smooth_column_action'): self.smooth_column_action.setEnabled(enable_recipes)
        if hasattr(self, 'import_systems_action'): self.import_systems_action.setEnabled(enable_recipes)
        if hasattr(self, 'import_systems_ascii_action'): self.import_systems_ascii_action.setEnabled(enable_recipes)
        if hasattr(self, 'delete_action'): self.delete_action.setEnabled(enable_recipes)
        if hasattr(self, 'split_action'): self.split_action.setEnabled(enable_recipes)
        if hasattr(self, 'extract_preset_action'): self.extract_preset_action.setEnabled(enable_recipes)
        if hasattr(self, 'coadd_action'): self.coadd_action.setEnabled(enable_recipes)
        
        if hasattr(self, 'calculate_running_std_action'): self.calculate_running_std_action.setEnabled(enable_recipes)
        if hasattr(self, 'smooth_action'): self.smooth_action.setEnabled(enable_recipes)
        if hasattr(self, 'rebin_action'): self.rebin_action.setEnabled(enable_recipes)
        if hasattr(self, 'equalize_action'): self.equalize_action.setEnabled(enable_recipes)
        if hasattr(self, 'resample_action'): self.resample_action.setEnabled(enable_recipes)
        if hasattr(self, 'calibrate_action'): self.calibrate_action.setEnabled(enable_recipes)

        if hasattr(self, 'auto_cont_action'): self.auto_cont_action.setEnabled(enable_recipes)
        if hasattr(self, 'find_abs_action'): self.find_abs_action.setEnabled(enable_recipes)
        if hasattr(self, 'fit_cont_action'): self.fit_cont_action.setEnabled(enable_recipes)
        if hasattr(self, 'fit_powerlaw_action'): self.fit_powerlaw_action.setEnabled(enable_recipes)

        if hasattr(self, 'doublets_auto_action'): self.doublets_auto_action.setEnabled(enable_recipes)
        if hasattr(self, 'lya_auto_action'): self.lya_auto_action.setEnabled(enable_recipes)
        if hasattr(self, 'forest_auto_action'): self.forest_auto_action.setEnabled(enable_recipes)
        if hasattr(self, 'metals_auto_action'): self.metals_auto_action.setEnabled(enable_recipes)
        if hasattr(self, 'identify_lines_action'): self.identify_lines_action.setEnabled(enable_recipes)
        if hasattr(self, 'refit_all_action'): self.refit_all_action.setEnabled(enable_recipes)
        if hasattr(self, 'clean_negligible_action'): self.clean_negligible_action.setEnabled(enable_recipes)
        # --- *** END NEW RECIPES *** ---
        
        # ... enable/disable other recipe actions ...

        # Enable Save action only if valid session
        if hasattr(self, 'save_action'): self.save_action.setEnabled(is_valid_session)

        # --- Enable/Disable Normalize F based on continuum ---
        if is_valid_session and self.active_history and self.active_history.current_state.spec:
            has_cont = self.active_history.current_state.spec.cont is not None
            if hasattr(self, 'norm_y_checkbox'):
                self.norm_y_checkbox.setEnabled(has_cont)
                if has_cont:
                    self.norm_y_checkbox.setToolTip("Plot F / Continuum and set F limits.")
                else:
                    self.norm_y_checkbox.setToolTip("Requires a continuum to normalize F.")
        else:
            if hasattr(self, 'norm_y_checkbox'):
                self.norm_y_checkbox.setEnabled(False)

        if is_valid_session:
            # --- State when a valid session IS loaded ---
            if self.central_stack.currentIndex() != 0: # If switching from empty
                self.central_stack.setCurrentIndex(0)
            
            """
            # 2. Resize logic (DE-NESTED from the stack index check)
            # Check if this is the very first session being loaded
            was_previously_empty = len(self.session_histories) <= 1
            if is_startup or was_previously_empty:
                logging.debug("Resizing and centering window for first valid session.")
                self.resize(1400, 900)
                # Recenter after resize
                screen_geometry = QApplication.primaryScreen().geometry()
                x = (screen_geometry.width() - self.width()) // 2
                y = (screen_geometry.height() - self.height()) // 2
                self.move(x, y)
            """

            # Show buttons
            self.session_collapse_button.setVisible(True)
            self.plot_controls_collapse_button.setVisible(True)

            # Update icons based on *current width* (more reliable during animation)
            self._update_sidebar_button_icon(self.session_collapse_button, self.left_sidebar_widget.width() > 0, is_left=True)
            self._update_sidebar_button_icon(self.plot_controls_collapse_button, self.right_sidebar_widget.width() > 0, is_left=False)

            # Ensure sidebars are visible if they have width > 0
            # (Animation handles the actual hiding/showing based on width)
            if self.left_sidebar_widget.width() > 0: self.left_sidebar_widget.setVisible(True)
            if self.right_sidebar_widget.width() > 0: self.right_sidebar_widget.setVisible(True)         

        else:
            if self.central_stack.currentIndex() != 1: self.central_stack.setCurrentIndex(1)

            # Hide buttons
            self.session_collapse_button.setVisible(False)
            self.plot_controls_collapse_button.setVisible(False)

            # Hide sidebars immediately (stop animations if running)
            if self.left_animation_group and self.left_animation_group.state() == QParallelAnimationGroup.Running: self.left_animation_group.stop()
            if self.right_animation_group and self.right_animation_group.state() == QParallelAnimationGroup.Running: self.right_animation_group.stop()
            self.left_sidebar_widget.setVisible(False); self.left_sidebar_widget.setGeometry(0, self.left_sidebar_widget.y(), 0, self.left_sidebar_widget.height()) # Reset geom
            self.right_sidebar_widget.setVisible(False); self.right_sidebar_widget.setGeometry(self.width(), self.right_sidebar_widget.y(), 0, self.right_sidebar_widget.height()) # Reset geom

        # Reposition elements after visibility/state changes might affect layout needs
        # QTimer.singleShot(0, self._reposition_floating_widgets) # Schedule reposition slightly later
        self._reposition_floating_widgets() # Try immediate reposition first

        self._update_undo_redo_actions()
