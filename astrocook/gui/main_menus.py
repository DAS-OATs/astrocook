"""The menu bar of the main window: its menus and entries (each recipe
entry opens its dialog), undo and redo, the windows of the View menu,
export and import of ASCII files, the documentation and the About dialog.

Moved out of ``main_window.py`` in phase 5 (split of the GUI), as a mixin
of ``MainWindow``.
"""
import logging
import os
from typing import Any, Optional

from PySide6.QtCore import QSettings, QTimer, QUrl, Qt
from PySide6.QtGui import QAction, QActionGroup, QDesktopServices, QKeySequence
from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox

from .about_dialog import AboutDialog
from .script_dialog import ScriptDialog
from .system_inspector import SystemInspector
from .plot_style import PLOT_PALETTE_NAMES, plot_palette_name, set_plot_palette
from .theme import APPEARANCES, apply_appearance
from astrocook.gui._window_host import HostBase
from astrocook.recipes.registry import DOCS_URL, find_recipe

logger = logging.getLogger(__name__)

#: The files opened last, in the settings (File > Open Recent, the empty window)
RECENT_FILES_KEY = "recent_files"
MAX_RECENT_FILES = 10
#: The size and position of the main window, in the settings
GEOMETRY_KEY = "window_geometry"
#: View > Plot Colours: the palette of the plots ("standard", "colorblind")
PLOT_COLOURS_KEY = "plot_colours"


#: The tooltips of the entries of the recipes (a longer text than the label)
RECIPE_TOOLTIPS = {
    'set_properties': 'Set core session properties (e.g., z_em)',
    'apply_expression': "Apply a numerical expression to columns (e.g., 'y / 2.0')",
    'mask_expression': "Mask a column using a boolean expression (e.g., 'x < 300')",
    'smooth_column': 'Smooth one column of the spectrum (e.g. the continuum or an auxiliary column); Flux > Smooth Spectrum smooths the flux and the columns that go with it',
    'delete': 'Delete columns or clear the line list',
    'extract_preset': 'Extract standard regions (Forest, etc.) based on z_em',
    'coadd': 'Combine multiple sessions and rebin to a single grid',
    'calculate_running_std': 'Calculate running standard deviation on a column (e.g., y)',
    'smooth': 'Smooth the flux and the columns that go with it (continuum, model, transmission); Edit > Smooth Column smooths one column',
    'equalize': "Scale this spectrum to match another session's flux level.",
    'resample': "Resample this spectrum onto another session's grid",
    'calibrate_from_magnitudes': 'Rescale spectrum to match a photometric magnitude',
    'estimate_auto': 'Find Absorbed regions and fit a continuum',
    'find_absorbed': "Create a 'abs_mask' column using V1 'clip_flux' logic",
    'fit_continuum': "Fit a continuum to the 'abs_mask' column (V1 logic)",
    'fit_powerlaw': 'Fit a power-law to specified rest-frame regions',
    'doublets_auto': 'Automated pipeline: Identify, populate, and fit metal doublets.',
    'lya_auto': 'Automated pipeline: Identify, populate, and fit Ly-alpha lines.',
    'forest_auto': 'Automated pipeline: Identify, populate, and fit the Lyman series.',
    'metals_auto': 'Scans all existing systems and attempts to add other ions at the same redshift.',
    'identify_lines': 'Automatically identify absorption regions using correlation signals',
    'refit_all': 'Optimizes all components in the spectrum by fitting disjoint groups.',
    'clean_negligible': 'Removes components with negligible column densities or broadening.',
    'import_systems': 'Copy absorption systems from another open session.',
    'import_systems_ascii': 'Read absorption systems from an ASCII/CSV table (series, z, logN, b, ...).',
}

#: The recipe entries of each menu, in order (``None``: a separator); the
#: label of each is the ``menu`` key of the recipe in the registry. The File
#: menu is arranged at the end of ``_create_menubar``.
RECIPE_LAYOUT = {
    "Edit": ("set_properties", None, "apply_expression", "mask_expression",
             "smooth_column", "delete", None, "split", "extract_preset", None, "coadd"),
    "Flux": ("calculate_running_std", None, "smooth", "rebin", "equalize", "resample",
             None, "calibrate_from_magnitudes"),
    "Continuum": ("estimate_auto", None, "find_absorbed", "fit_continuum", None,
                  "fit_powerlaw"),
    "Absorbers": ("doublets_auto", "lya_auto", None, "identify_lines", "refit_all",
                  "clean_negligible"),
}
#: Recipes with an action but no entry in the menu, and no ``menu`` key in
#: the registry (kept for when they come back)
RECIPES_OFF_MENU = ("forest_auto", "metals_auto")


class MainMenusMixin(HostBase):
    """The menus of ``MainWindow``.

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

    def _set_plot_colours(self, name: str, save: bool = True):
        """View > Plot Colours: the standard palette of the plots or the
        colour-blind safe one (``plot_style.PLOT_PALETTES``). The main plot
        and the System Inspector are drawn again."""
        name = set_plot_palette(name)
        self.plot_colours_actions[name].setChecked(True)
        if save:
            QSettings("Astrocook", "AstrocookV2").setValue(PLOT_COLOURS_KEY, name)
        if getattr(self, 'plot_viewer', None) is not None and self.active_history:
            self._trigger_replot()
        inspector = getattr(self, 'system_inspector', None)
        if inspector is not None:
            inspector.refresh_plot_colours()

    def _set_toolbar_visible(self, visible: bool, save: bool = True):
        """View > Show Toolbar: the toolbar of the main plot (with the
        coordinates under the mouse) shown or hidden."""
        if getattr(self, 'plot_viewer', None) is not None:
            self.plot_viewer.toolbar.set_shown(visible)
            self._reposition_floating_widgets()  # the session card stops above it
        if self.show_toolbar_action.isChecked() != visible:
            self.show_toolbar_action.setChecked(visible)
        if save:
            QSettings("Astrocook", "AstrocookV2").setValue("show_toolbar", visible)

    # --- Recent files ---

    def _recent_files(self) -> list:
        """The files opened last, most recent first, those that still exist."""
        value: Any = QSettings("Astrocook", "AstrocookV2").value(RECENT_FILES_KEY, [])
        if isinstance(value, str):  # one entry comes back as a string
            value = [value]
        return [p for p in (value or []) if os.path.exists(p)]

    def _add_recent_file(self, path: str):
        path = os.path.abspath(path)
        recent = [p for p in self._recent_files() if p != path]
        QSettings("Astrocook", "AstrocookV2").setValue(
            RECENT_FILES_KEY, ([path] + recent)[:MAX_RECENT_FILES])
        if hasattr(self, '_refresh_recent_list'):
            self._refresh_recent_list()

    def _clear_recent_files(self):
        QSettings("Astrocook", "AstrocookV2").remove(RECENT_FILES_KEY)
        if hasattr(self, '_refresh_recent_list'):
            self._refresh_recent_list()

    def _rebuild_recent_menu(self):
        """File > Open Recent: the recent files, and Clear Menu."""
        self.recent_menu.clear()
        recent = self._recent_files()
        for path in recent:
            action = self.recent_menu.addAction(os.path.basename(path))
            action.setToolTip(path)
            action.triggered.connect(lambda checked=False, p=path: self.open_files([p]))
        if not recent:
            self.recent_menu.addAction("No Recent Files").setEnabled(False)
        self.recent_menu.addSeparator()
        clear = self.recent_menu.addAction("Clear Menu", self._clear_recent_files)
        clear.setEnabled(bool(recent))
        self.recent_menu.setToolTipsVisible(True)

    def _restore_view_settings(self):
        """The appearance and the toolbar as the user left them."""
        settings = QSettings("Astrocook", "AstrocookV2")
        self._set_appearance(settings.value("appearance", "system"), save=False)
        self._set_plot_colours(settings.value(PLOT_COLOURS_KEY, "standard"), save=False)
        visible = settings.value("show_toolbar", True)
        self._set_toolbar_visible(visible not in (False, "false"), save=False)
        # The size and position of the window, as the user left them (they
        # were not saved)
        geometry = settings.value(GEOMETRY_KEY)
        if geometry is not None:
            self.restoreGeometry(geometry)

    def _save_window_geometry(self):
        QSettings("Astrocook", "AstrocookV2").setValue(GEOMETRY_KEY, self.saveGeometry())

    def _recipe_action(self, name: str, title: Optional[str], handler=None) -> QAction:
        """The (disabled) action of the entry of recipe ``name`` in the menu
        ``title``: its text, ``Set Properties...``, comes from the ``menu``
        key of the registry (the only place where it is written), and its
        action is ``self.<name>_action``. It opens the dialog of the recipe
        unless ``handler`` says otherwise. ``ValueError`` if the registry
        puts the recipe in another menu. With ``title=None`` (the recipes kept
        out of the menus, which have no ``menu``) the text is the brief."""
        spec = find_recipe(name)
        if spec is None:
            raise ValueError(f"No recipe {name!r}")
        if spec.menu_title != (title or ""):
            raise ValueError(f"{name} is in the menu {spec.menu_title!r}, not in {title!r}")
        action = QAction((spec.menu_label or spec.brief.rstrip('.')) + "...", self)
        if name in RECIPE_TOOLTIPS:
            action.setToolTip(RECIPE_TOOLTIPS[name])
        if handler is None:
            category = spec.category

            def launch():
                initial = self._recipe_initial_params(name)
                if initial:
                    self._launch_recipe_dialog(category, name, initial_params=initial)
                else:
                    self._launch_recipe_dialog(category, name)
            action.triggered.connect(launch)
        else:
            action.triggered.connect(handler)
        action.setEnabled(False)
        setattr(self, f"{name}_action", action)
        return action

    def _recipe_initial_params(self, name: str) -> dict:
        """What the dialog of the recipe starts from, besides its defaults:
        Co-add Spectra starts from the active session."""
        if name == "coadd" and self.active_history:
            return {'session_names': self.active_history.display_name}
        return {}

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
        # The tooltips of the entries (e.g. the two Smooth entries) shown
        for menu in (file_menu, view_menu, edit_menu, flux_menu, continuum_menu,
                     absorbers_menu, help_menu):
            menu.setToolTipsVisible(True)
        
        # --- File and View operations ---

        # --- File Menu Actions
        open_action = QAction("&Open Session...", self); open_action.setShortcut("Ctrl+O"); open_action.triggered.connect(self._on_open_spectrum)
        file_menu.addAction(open_action)
        # The files opened last (saved in the settings), rebuilt each time
        # the menu opens
        self.recent_menu = file_menu.addMenu("Open &Recent")
        self.recent_menu.aboutToShow.connect(self._rebuild_recent_menu)
        if hasattr(self, '_refresh_recent_list'):
            self._refresh_recent_list()  # the button of the welcome screen
        save_action = QAction("&Save Session...", self); save_action.setShortcut("Ctrl+S"); save_action.triggered.connect(self._on_save_session); 
        save_action.setEnabled(False)
        file_menu.addAction(save_action)
        self.save_action = save_action # Assign to self

        export_action = self._recipe_action("export_ascii", "File", self._on_export_requested)
        # Insert it after Save
        file_menu.insertAction(self.save_action, export_action) 
        self.export_action = export_action

        # The log as the text of a methods section (phase 6-ter, PR 5c)
        methods_action = QAction("Export &Log as Text...", self)
        methods_action.setToolTip("The session log as the text of a methods section, "
                                  "to copy or save.")
        methods_action.triggered.connect(self._on_export_log_text)
        self.export_log_text_action = methods_action
        self.export_log_text_action.setEnabled(False)

        import_action = self._recipe_action("import_ascii", "File", self._on_import_requested)
        # Insert it after Save
        file_menu.insertAction(self.save_action, import_action) 
        self.import_action = import_action

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

        # The palette of the plots: standard or colour-blind safe (saved)
        colours_menu = view_menu.addMenu("Plot &Colours")
        self.plot_colours_actions = {}
        group = QActionGroup(self)
        group.setExclusive(True)
        for name, label in PLOT_PALETTE_NAMES.items():
            action = QAction(label, self)
            action.setCheckable(True)
            action.setChecked(name == plot_palette_name())
            action.triggered.connect(lambda checked=False, n=name: self._set_plot_colours(n))
            group.addAction(action)
            colours_menu.addAction(action)
            self.plot_colours_actions[name] = action

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

        # The entries of the recipes: label from the registry (its ``menu``
        # key), in the order of RECIPE_LAYOUT
        for title, menu in (("Edit", edit_menu), ("Flux", flux_menu),
                            ("Continuum", continuum_menu), ("Absorbers", absorbers_menu)):
            for name in RECIPE_LAYOUT[title]:
                if name is None:
                    menu.addSeparator()
                else:
                    menu.addAction(self._recipe_action(name, title))
        # Two pipelines kept out of the menu (their actions stay, enabled with the others)
        for name in RECIPES_OFF_MENU:
            self._recipe_action(name, None)

        # In File with the other imports (they were in Edit; phase 5-bis,
        # Guido), before Save
        import_sys_action = self._recipe_action("import_systems", "File")
        file_menu.insertAction(self.save_action, import_sys_action)
        import_sys_ascii_action = self._recipe_action("import_systems_ascii", "File",
                                                      self._import_systems_ascii)
        file_menu.insertAction(import_sys_action, import_sys_ascii_action)
        self._update_undo_redo_actions()

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

        # The File menu in groups (Guido): open, save, import, export, close
        # (removed one by one: clear() deletes the actions the menu owns)
        for action in file_menu.actions():
            file_menu.removeAction(action)
        for group in ((open_action, self.recent_menu.menuAction()), (save_action,),
                      (import_action, self.import_systems_ascii_action,
                       self.import_systems_action),
                      (export_action, methods_action), (close_action,)):
            if file_menu.actions():
                file_menu.addSeparator()
            for action in group:
                file_menu.addAction(action)

    def _undo_last_action(self):
        """Switches the view to the previous state in the active history."""
        if self.is_followed() and self._refuse_if_followed("Undo"):
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            if self.active_history:
                previous_state = self.active_history.undo()
                if previous_state:
                    logger.debug(f"Undo: Switched to state index {self.active_history.current_index}")
                    self._update_view_for_session(previous_state, set_current_list_item=False) 
                    self._update_undo_redo_actions()

                    # --- MODIFIED CHECK ---
                    if self.log_scripter_dialog and self.log_scripter_dialog.isVisible():
                        self.log_scripter_dialog.refresh()
                    # --- END MODIFIED CHECK ---
                    QTimer.singleShot(0, self._force_restack_floating_widgets)
                else:
                    logger.debug("Undo: Already at oldest state for this session.")
            else:
                 logger.debug("Undo: No active session.")
        finally:
            QApplication.restoreOverrideCursor()

    def _redo_last_action(self):
        """Switches the view to the next state in the active history."""
        if self.is_followed() and self._refuse_if_followed("Redo"):
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            if self.active_history:
                next_state = self.active_history.redo()
                if next_state:
                    logger.debug(f"Redo: Switched to state index {self.active_history.current_index}")
                    self._update_view_for_session(next_state, set_current_list_item=False)
                    self._update_undo_redo_actions()

                    # --- MODIFIED CHECK ---
                    if self.log_scripter_dialog and self.log_scripter_dialog.isVisible():
                        self.log_scripter_dialog.refresh()
                    # --- END MODIFIED CHECK ---
                    QTimer.singleShot(0, self._force_restack_floating_widgets)
                else:
                    logger.debug("Redo: Already at newest state for this session.")
            else:
                 logger.debug("Redo: No active session.")
        finally:
            QApplication.restoreOverrideCursor()

    def _update_undo_redo_actions(self):
        """Enables/disables Undo/Redo based on the active history manager."""
        can_undo = self.active_history.can_undo() if self.active_history else False
        can_redo = self.active_history.can_redo() if self.active_history else False
        if self.is_followed():
            can_undo = can_redo = False  # read only (phase 6-ter, PR 5-bis)

        if hasattr(self, 'undo_action'): self.undo_action.setEnabled(can_undo)
        if hasattr(self, 'redo_action'): self.redo_action.setEnabled(can_redo)

    def _on_about_requested(self):
        """Launches the About dialog."""
        dlg = AboutDialog(self)
        dlg.exec()

    def _open_docs(self):
        """Opens the online documentation in the default browser."""
        QDesktopServices.openUrl(QUrl(DOCS_URL))

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
            logger.debug("View Log without an active session.")

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
            logger.debug("View Identifications without an active session.")

    def _on_export_requested(self):
        """Launches the selection dialog first."""
        if not self.active_history:
            return
        # No QFileDialog here anymore! Just the dialog.
        self._launch_recipe_dialog("file", "export_ascii")

    def _on_export_log_text(self):
        """File > Export Log as Text: the log of the active session as the
        text of a methods section, in a dialog to copy or save it."""
        if not self.active_history:
            return
        from astrocook.gui.methods_text_dialog import MethodsTextDialog
        dialog = MethodsTextDialog(self.active_history, self)
        self.methods_text_dialog = dialog
        dialog.show()
        dialog.raise_()

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
