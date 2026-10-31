"""The menu bar of the main window: its menus and entries (each recipe
entry opens its dialog), undo and redo, the windows of the View menu,
export and import of ASCII files, the documentation and the About dialog.

Moved out of ``main_window.py`` in phase 5 (split of the GUI), as a mixin
of ``MainWindowV2``.
"""
import logging

from PySide6.QtCore import QSettings, QTimer, QUrl, Qt
from PySide6.QtGui import QAction, QActionGroup, QDesktopServices, QKeySequence
from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox

from .about_dialog import AboutDialog
from .script_dialog import ScriptDialog
from .system_inspector import SystemInspector
from .theme import APPEARANCES, apply_appearance


class MainMenusMixin:
    """The menus of ``MainWindowV2``.

    Uses the window's ``active_history``, the recipe machinery
    (``_launch_recipe_dialog``, ``_on_recipe_requested``), the sidebars,
    the dialogs it opens (``system_inspector``, ``script_dialog``, the log
    scripter, the identification viewer) and the view update
    (``_update_view_for_session``); the actions it creates are enabled by
    ``_update_ui_state``.
    """

    def _set_appearance(self, mode: str, save: bool = True):
        """View > Appearance: light, dark, or as the system. The windows
        apply the new palette by themselves (ApplicationPaletteChange)."""
        if mode not in APPEARANCES:
            mode = "system"
        apply_appearance(mode)
        self.appearance_actions[mode].setChecked(True)
        if save:
            QSettings("Astrocook", "AstrocookV2").setValue("appearance", mode)

    def _set_toolbar_visible(self, visible: bool, save: bool = True):
        """View > Show Toolbar: the toolbar of the main plot (with the
        coordinates under the mouse) shown or hidden."""
        if getattr(self, 'plot_viewer', None) is not None:
            self.plot_viewer.toolbar.set_shown(visible)
        if self.show_toolbar_action.isChecked() != visible:
            self.show_toolbar_action.setChecked(visible)
        if save:
            QSettings("Astrocook", "AstrocookV2").setValue("show_toolbar", visible)

    def _restore_view_settings(self):
        """The appearance and the toolbar as the user left them."""
        settings = QSettings("Astrocook", "AstrocookV2")
        self._set_appearance(settings.value("appearance", "system"), save=False)
        visible = settings.value("show_toolbar", True)
        self._set_toolbar_visible(visible not in (False, "false"), save=False)

    def _create_menubar(self):
        menu_bar = self.menuBar()
        
        # --- Define all 6 Scientific Menus ---
        
        # 1. Standard Application Menus
        file_menu = menu_bar.addMenu("&File")
        view_menu = menu_bar.addMenu("&View")

        # 2. Scientific Interaction Menus
        edit_menu = menu_bar.addMenu("&Edit")
        flux_menu = menu_bar.addMenu("&Flux")
        continuum_menu = menu_bar.addMenu("&Continuum")
        absorbers_menu = menu_bar.addMenu("&Absorbers")

        # 3. Help Menu
        help_menu = menu_bar.addMenu("&Help")
        
        # --- File and View operations ---

        # --- File Menu Actions
        open_action = QAction("&Open Session...", self); open_action.setShortcut("Ctrl+O"); open_action.triggered.connect(self._on_open_spectrum)
        file_menu.addAction(open_action)
        save_action = QAction("&Save Session...", self); save_action.setShortcut("Ctrl+S"); save_action.triggered.connect(self._on_save_session); 
        save_action.setEnabled(False)
        file_menu.addAction(save_action)
        self.save_action = save_action # Assign to self

        export_action = QAction("&Export to ASCII...", self)
        export_action.triggered.connect(self._on_export_requested)
        # Insert it after Save
        file_menu.insertAction(self.save_action, export_action) 
        self.export_action = export_action
        self.export_action.setEnabled(False) # Default disabled

        import_action = QAction("&Import from ASCII...", self)
        import_action.triggered.connect(self._on_import_requested)
        # Insert it after Save
        file_menu.insertAction(self.save_action, import_action) 
        self.import_action = import_action
        self.import_action.setEnabled(False) # Default disabled

        file_menu.addSeparator()
        close_action = QAction("&Close Session", self)
        close_action.setShortcut("Ctrl+W")
        close_action.triggered.connect(self._on_close_session_requested)
        close_action.setEnabled(False)
        file_menu.addAction(close_action)
        self.close_session_action = close_action

        # 'VIEW' MENU ACTIONS
        # Checkable, as Show Toolbar: the check follows the panel, also when it
        # is opened or closed by its tab (_update_sidebar_button_icon)
        toggle_left_action = QAction("Show Session Panel", self); toggle_left_action.setCheckable(True)
        toggle_left_action.triggered.connect(self._toggle_left_sidebar)
        view_menu.addAction(toggle_left_action)
        toggle_right_action = QAction("Show Plot Controls", self); toggle_right_action.setCheckable(True)
        toggle_right_action.triggered.connect(self._toggle_right_sidebar)
        view_menu.addAction(toggle_right_action)
        self.toggle_left_action = toggle_left_action
        self.toggle_right_action = toggle_right_action

        # The toolbar of the plot (its tools have keys: it can be hidden)
        self.show_toolbar_action = QAction("Show &Toolbar", self)
        self.show_toolbar_action.setCheckable(True)
        self.show_toolbar_action.setChecked(True)
        self.show_toolbar_action.toggled.connect(self._set_toolbar_visible)
        # The tab of the floating toolbar does the same
        if getattr(self, 'plot_viewer', None) is not None:
            self.plot_viewer.toolbar.on_toggled = self._set_toolbar_visible
        view_menu.addAction(self.show_toolbar_action)

        # Light, dark, or as the system (saved in the settings)
        appearance_menu = view_menu.addMenu("&Appearance")
        self.appearance_actions = {}
        group = QActionGroup(self)
        group.setExclusive(True)
        for mode in APPEARANCES:
            action = QAction(mode.capitalize(), self)
            action.setCheckable(True)
            action.triggered.connect(lambda checked=False, m=mode: self._set_appearance(m))
            group.addAction(action)
            appearance_menu.addAction(action)
            self.appearance_actions[mode] = action

        view_menu.addSeparator()
        self.view_system_inspector_action = QAction("View System &Inspector", self)
        self.view_system_inspector_action.triggered.connect(self._on_view_system_inspector)
        self.view_system_inspector_action.setEnabled(False) # Disabled by default
        view_menu.addAction(self.view_system_inspector_action)

        view_menu.addSeparator()
        self.view_log_action = QAction("View Session &Log", self)
        self.view_log_action.triggered.connect(self._on_view_log)
        self.view_log_action.setEnabled(False) # Disabled by default
        view_menu.addAction(self.view_log_action)

        self.view_script_console_action = QAction("View Script &Console", self)
        self.view_script_console_action.triggered.connect(self._on_view_script_console)
        self.view_script_console_action.setEnabled(False) # Disabled until session loaded
        view_menu.addAction(self.view_script_console_action)

        self.view_identifications_action = QAction("View &Identifications", self)
        self.view_identifications_action.triggered.connect(self._on_view_identifications)
        self.view_identifications_action.setEnabled(False) # Disabled by default
        view_menu.addAction(self.view_identifications_action)

        # 'EDIT' MENU ACTIONS

        # ** Undo Action **
        self.undo_action = QAction("&Undo", self)
        self.undo_action.setShortcut(QKeySequence.Undo) # Standard shortcut (Ctrl+Z / Cmd+Z)
        self.undo_action.setShortcutContext(Qt.ShortcutContext.ApplicationShortcut)
        self.undo_action.triggered.connect(self._undo_last_action)
        edit_menu.addAction(self.undo_action)

        # ** Redo Action **
        self.redo_action = QAction("&Redo", self)
        self.redo_action.setShortcut(QKeySequence.Redo) # Standard shortcut (Ctrl+Y / Cmd+Shift+Z)
        self.redo_action.setShortcutContext(Qt.ShortcutContext.ApplicationShortcut)
        self.redo_action.triggered.connect(self._redo_last_action)
        edit_menu.addAction(self.redo_action)

        edit_menu.addSeparator()

        # RECIPES FOR 'EDIT' MENU

        # set_properties Action
        set_props_action = QAction("Set &Properties...", self)
        set_props_action.setToolTip("Set core session properties (e.g., z_em)")
        set_props_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "set_properties"))
        edit_menu.addAction(set_props_action)
        self.set_properties_action = set_props_action; self.set_properties_action.setEnabled(False)
        edit_menu.addSeparator()
        
        # apply_expression Action
        apply_expr_action = QAction("Apply &Expression...", self)
        apply_expr_action.setToolTip("Apply a numerical expression to columns (e.g., 'y / 2.0')")
        apply_expr_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "apply_expression"))
        edit_menu.addAction(apply_expr_action)
        self.apply_expression_action = apply_expr_action; self.apply_expression_action.setEnabled(False)
        
        # mask_expression Action
        mask_expr_action = QAction("&Mask by Expression...", self)
        mask_expr_action.setToolTip("Mask a column using a boolean expression (e.g., 'x < 300')")
        mask_expr_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "mask_expression"))
        edit_menu.addAction(mask_expr_action)
        self.mask_expression_action = mask_expr_action; self.mask_expression_action.setEnabled(False)

        # smooth_column Action
        smooth_col_action = QAction("Smooth &Column...", self)
        smooth_col_action.setToolTip("Apply Gaussian smoothing to a single column")
        smooth_col_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "smooth_column"))
        edit_menu.addAction(smooth_col_action)
        self.smooth_column_action = smooth_col_action; self.smooth_column_action.setEnabled(False)

        # import_systems Action
        import_sys_action = QAction("&Import Systems...", self)
        import_sys_action.setToolTip("Copy absorption systems from another open session.")
        # Note: We call the 'absorbers' category because that's where the recipe logic lives
        import_sys_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "import_systems"))
        edit_menu.addAction(import_sys_action)
        self.import_systems_action = import_sys_action; self.import_systems_action.setEnabled(False)

        # import_systems_ascii Action (issue #480)
        import_sys_ascii_action = QAction("Import Systems from &ASCII...", self)
        import_sys_ascii_action.setToolTip("Read absorption systems from an ASCII/CSV table (series, z, logN, b, ...).")
        import_sys_ascii_action.triggered.connect(self._import_systems_ascii)
        edit_menu.addAction(import_sys_ascii_action)
        self.import_systems_ascii_action = import_sys_ascii_action
        self.import_systems_ascii_action.setEnabled(False)

        # delete Action
        delete_action = QAction("&Delete Elements...", self)
        delete_action.setToolTip("Delete columns or clear the line list")
        delete_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "delete"))
        edit_menu.addAction(delete_action)
        self.delete_action = delete_action; self.delete_action.setEnabled(False)

        edit_menu.addSeparator()

        # split Action
        split_action = QAction("S&plit Out Region...", self)
        split_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "split"))
        edit_menu.addAction(split_action)
        self.split_action = split_action; split_action.setEnabled(False)

        # extract_preset Action
        extract_preset_action = QAction("Extract &Preset Region...", self)
        extract_preset_action.setToolTip("Extract standard regions (Forest, etc.) based on z_em")
        extract_preset_action.triggered.connect(lambda: self._launch_recipe_dialog("edit", "extract_preset"))
        edit_menu.addAction(extract_preset_action)
        self.extract_preset_action = extract_preset_action; self.extract_preset_action.setEnabled(False)

        edit_menu.addSeparator()

        # coadd Action
        coadd_action = QAction("Co-&add Spectra...", self)
        coadd_action.setToolTip("Combine multiple sessions and rebin to a single grid")
        coadd_action.triggered.connect(lambda: self._launch_recipe_dialog(
            "edit", "coadd",
            initial_params={'session_names': self.active_history.display_name if self.active_history else ""}
        ))
        edit_menu.addAction(coadd_action)
        self.coadd_action = coadd_action; self.coadd_action.setEnabled(False)

        # RECIPES FOR 'FLUX' MENU
        
        std_action = QAction("Calculate Running &StdDev...", self)
        std_action.setToolTip("Calculate running standard deviation on a column (e.g., y)")
        std_action.triggered.connect(
            lambda: self._launch_recipe_dialog("flux", "calculate_running_std")
        )
        flux_menu.addAction(std_action)
        self.calculate_running_std_action = std_action; self.calculate_running_std_action.setEnabled(False)

        flux_menu.addSeparator()

        # smooth Action
        smooth_action = QAction("&Smooth Spectrum...", self)
        smooth_action.setToolTip("Apply Gaussian smoothing to the flux")
        smooth_action.triggered.connect(lambda: self._launch_recipe_dialog("flux", "smooth"))
        flux_menu.addAction(smooth_action)
        self.smooth_action = smooth_action; self.smooth_action.setEnabled(False)

        # rebin Action
        rebin_action = QAction("&Rebin Spectrum...", self)
        rebin_action.triggered.connect(lambda: self._launch_recipe_dialog("flux", "rebin"))
        flux_menu.addAction(rebin_action)
        self.rebin_action = rebin_action; rebin_action.setEnabled(False)

        # equalize Action
        equalize_action = QAction("&Equalize to Reference...", self)
        equalize_action.setToolTip("Scale this spectrum to match another session's flux level.")
        equalize_action.triggered.connect(lambda: self._launch_recipe_dialog("flux", "equalize"))
        flux_menu.addAction(equalize_action)
        self.equalize_action = equalize_action; self.equalize_action.setEnabled(False)

        # resample Action
        resample_action = QAction("Re&sample on Grid...", self)
        resample_action.setToolTip("Resample this spectrum onto another session's grid")
        resample_action.triggered.connect(lambda: self._launch_recipe_dialog("flux", "resample"))
        flux_menu.addAction(resample_action)
        self.resample_action = resample_action; self.resample_action.setEnabled(False)

        flux_menu.addSeparator()

        # flux calibrate Action
        calib_action = QAction("Flux &Calibrate...", self)
        calib_action.setToolTip("Rescale spectrum to match a photometric magnitude")
        calib_action.triggered.connect(lambda: self._launch_recipe_dialog("flux", "calibrate_from_magnitudes"))
        flux_menu.addAction(calib_action)
        self.calibrate_action = calib_action; self.calibrate_action.setEnabled(False)

        # RECIPES FOR 'CONTINUUM' MENU

        auto_cont_action = QAction("&Auto-estimate Continuum...", self)
        auto_cont_action.setToolTip("Find Absorbed regions and fit a continuum")
        auto_cont_action.triggered.connect(lambda: self._launch_recipe_dialog("continuum", "estimate_auto"))
        continuum_menu.addAction(auto_cont_action)
        self.auto_cont_action = auto_cont_action; self.auto_cont_action.setEnabled(False)

        continuum_menu.addSeparator()

        find_abs_action = QAction("&Find Absorbed Regions...", self)
        find_abs_action.setToolTip("Create a 'abs_mask' column using V1 'clip_flux' logic")
        find_abs_action.triggered.connect(lambda: self._launch_recipe_dialog("continuum", "find_absorbed"))
        continuum_menu.addAction(find_abs_action)
        self.find_abs_action = find_abs_action; self.find_abs_action.setEnabled(False)

        fit_cont_action = QAction("Fit &Continuum to Mask...", self)
        fit_cont_action.setToolTip("Fit a continuum to the 'abs_mask' column (V1 logic)")
        fit_cont_action.triggered.connect(lambda: self._launch_recipe_dialog("continuum", "fit_continuum"))
        continuum_menu.addAction(fit_cont_action)
        self.fit_cont_action = fit_cont_action; self.fit_cont_action.setEnabled(False)

        continuum_menu.addSeparator()

        fit_pl_action = QAction("Fit &Power-Law...", self)
        fit_pl_action.setToolTip("Fit a power-law to specified rest-frame regions")
        fit_pl_action.triggered.connect(lambda: self._launch_recipe_dialog("continuum", "fit_powerlaw"))
        continuum_menu.addAction(fit_pl_action)
        self.fit_powerlaw_action = fit_pl_action; self.fit_powerlaw_action.setEnabled(False)

        # RECIPES FOR 'ABSORBERS' MENU

        # --- AUTOMATED PIPELINES (Prominent at top) ---
        doublets_auto_action = QAction("Auto-detect &Doublets...", self)
        doublets_auto_action.setToolTip("Automated pipeline: Identify, populate, and fit metal doublets.")
        doublets_auto_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "doublets_auto"))
        absorbers_menu.addAction(doublets_auto_action)
        self.doublets_auto_action = doublets_auto_action; self.doublets_auto_action.setEnabled(False)
        
        lya_auto_action = QAction("Auto-fit &Ly-alpha Forest...", self)
        lya_auto_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "lya_auto"))
        lya_auto_action.setToolTip("Automated pipeline: Identify, populate, and fit Ly-alpha lines.")
        absorbers_menu.addAction(lya_auto_action)
        self.lya_auto_action = lya_auto_action; self.lya_auto_action.setEnabled(False)

        forest_auto_action = QAction("Auto-fit &Full Lyman Series...", self)
        forest_auto_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "forest_auto"))
        forest_auto_action.setToolTip("Automated pipeline: Identify, populate, and fit the Lyman series.")
        #absorbers_menu.addAction(forest_auto_action)
        self.forest_auto_action = forest_auto_action; self.forest_auto_action.setEnabled(False)

        metals_auto_action = QAction("Check &Associated Metals...", self)
        metals_auto_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "metals_auto"))
        metals_auto_action.setToolTip("Scans all existing systems and attempts to add other ions at the same redshift.")
        #absorbers_menu.addAction(metals_auto_action)
        self.metals_auto_action = metals_auto_action; self.metals_auto_action.setEnabled(False)

        absorbers_menu.addSeparator()

        identify_action = QAction("&Identify Absorption Lines...", self)
        identify_action.setToolTip("Automatically identify absorption regions using correlation signals")
        identify_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "identify_lines"))
        absorbers_menu.addAction(identify_action)
        self.identify_lines_action = identify_action; self.identify_lines_action.setEnabled(False)
        
        refit_all_action = QAction("&Refit All Systems...", self)
        refit_all_action.setToolTip("Optimizes all components in the spectrum by fitting disjoint groups.")
        refit_all_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "refit_all"))
        absorbers_menu.addAction(refit_all_action)
        self.refit_all_action = refit_all_action; self.refit_all_action.setEnabled(False)

        clean_negligible_action = QAction("&Clean Negligible Components...", self)
        clean_negligible_action.setToolTip("Removes components with negligible column densities or broadening.")
        clean_negligible_action.triggered.connect(lambda: self._launch_recipe_dialog("absorbers", "clean_negligible"))
        absorbers_menu.addAction(clean_negligible_action)
        self.clean_negligible_action = clean_negligible_action; self.clean_negligible_action.setEnabled(False)
        self._update_undo_redo_actions()

        # --- HELP MENU ACTIONS ---
        
        # 1. Documentation Action (Keeps the menu alive on macOS)
        doc_action = QAction("&Online Documentation", self)
        doc_action.setShortcut(QKeySequence.HelpContents) # Standard F1 / Cmd+? shortcut
        doc_action.triggered.connect(self._open_docs)
        help_menu.addAction(doc_action)

        help_menu.addSeparator()

        # 2. About Action (Will move to App Menu on macOS)
        about_action = QAction("&About Astrocook...", self)
        about_action.triggered.connect(self._on_about_requested)
        help_menu.addAction(about_action)

    def _undo_last_action(self):
        """Switches the view to the previous state in the active history."""
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            if self.active_history:
                previous_state = self.active_history.undo()
                if previous_state:
                    logging.debug(f"Undo: Switched to state index {self.active_history.current_index}")
                    self._update_view_for_session(previous_state, set_current_list_item=False) 
                    self._update_undo_redo_actions()

                    # --- MODIFIED CHECK ---
                    if self.log_scripter_dialog and self.log_scripter_dialog.isVisible():
                        self.log_scripter_dialog.refresh()
                    # --- END MODIFIED CHECK ---
                    QTimer.singleShot(0, self._force_restack_floating_widgets)
                else:
                    logging.debug("Undo: Already at oldest state for this session.")
            else:
                 logging.debug("Undo: No active session.")
        finally:
            QApplication.restoreOverrideCursor()

    def _redo_last_action(self):
        """Switches the view to the next state in the active history."""
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            if self.active_history:
                next_state = self.active_history.redo()
                if next_state:
                    logging.debug(f"Redo: Switched to state index {self.active_history.current_index}")
                    self._update_view_for_session(next_state, set_current_list_item=False)
                    self._update_undo_redo_actions()

                    # --- MODIFIED CHECK ---
                    if self.log_scripter_dialog and self.log_scripter_dialog.isVisible():
                        self.log_scripter_dialog.refresh()
                    # --- END MODIFIED CHECK ---
                    QTimer.singleShot(0, self._force_restack_floating_widgets)
                else:
                    logging.debug("Redo: Already at newest state for this session.")
            else:
                 logging.debug("Redo: No active session.")
        finally:
            QApplication.restoreOverrideCursor()

    def _update_undo_redo_actions(self):
        """Enables/disables Undo/Redo based on the active history manager."""
        can_undo = self.active_history.can_undo() if self.active_history else False
        can_redo = self.active_history.can_redo() if self.active_history else False

        if hasattr(self, 'undo_action'): self.undo_action.setEnabled(can_undo)
        if hasattr(self, 'redo_action'): self.redo_action.setEnabled(can_redo)

    def _on_about_requested(self):
        """Launches the About dialog."""
        dlg = AboutDialog(self)
        dlg.exec()

    def _open_docs(self):
        """Opens the online documentation in the default browser."""
        # Replace with your actual documentation URL (e.g., ReadTheDocs or GitHub Wiki)
        url = "https://das-oats.github.io/astrocook/" 
        QDesktopServices.openUrl(QUrl(url))

    def _on_view_system_inspector(self):
        if not hasattr(self, 'system_inspector') or self.system_inspector is None:
            self.system_inspector = SystemInspector(self)
        
        # Initialize with current session
        if self.active_history:
            self.system_inspector.set_session(self.active_history.current_state)
            
        self.system_inspector.show()
        self.system_inspector.raise_()

    def _on_view_log(self):
        """
        Handles the "View > View Session Log" menu action.
        Shows/raises the *single, persistent* log scripter.
        """
        if self.active_history:
            self._launch_log_scripter(self.active_history)
        else:
            logging.warning("View Log called with no active history.")

    def _on_view_script_console(self):
        """Launches the Script Console dialog."""
        if not self.script_dialog:
            self.script_dialog = ScriptDialog(self)
        
        self.script_dialog.show()
        self.script_dialog.raise_()
        self.script_dialog.activateWindow()

    def _on_view_identifications(self):
        """
        Handles the "View > View Identifications" menu action.
        """
        if self.active_history:
            self._launch_identification_viewer(self.active_history)
        else:
            logging.warning("View Identifications called with no active history.")

    def _on_export_requested(self):
        """Launches the selection dialog first."""
        if not self.active_history:
            return
        # No QFileDialog here anymore! Just the dialog.
        self._launch_recipe_dialog("file", "export_ascii")

    def _on_import_requested(self):
        if not self.active_history: return

        # 1. The Disclaimer
        header = "Import and Overwrite Data?"
        text = ("You are about to import data from an ASCII file. "
                "Columns with matching names will be <b>overwritten</b>.<br><br>"
                "<i>Note: The number of rows in the file must exactly match the current spectrum.</i>")

        confirm = self._show_custom_message(
            title="Import ASCII",
            header=header,
            text=text,
            buttons=QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel,
            default_btn=QMessageBox.StandardButton.Cancel
        )

        if confirm != QMessageBox.StandardButton.Ok: return

        # 2. File Browser
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select ASCII/CSV File to Import", "", 
            "CSV Files (*.csv);;Text Files (*.txt *.dat);;All Files (*)"
        )

        if file_path:
            # 3. Execute via standard recipe pipeline to ensure it's in the Undo history
            self._on_recipe_requested("file", "import_ascii", {"file_path": file_path}, {})

    def _import_systems_ascii(self):
        """Choose an ASCII/CSV system table and import it (appended, not refitted)."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select ASCII/CSV System Table", "",
            "CSV Files (*.csv);;Text Files (*.txt *.dat);;All Files (*)"
        )
        if file_path:
            self._on_recipe_requested("edit", "import_systems_ascii",
                                      {"file_path": file_path, "append": "True",
                                       "refit": "False"}, {})
