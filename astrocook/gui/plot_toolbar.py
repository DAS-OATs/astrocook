"""The toolbar of the plots (main plot, velocity plot, corner plot): the
matplotlib toolbar without the buttons Astrocook does not use, and a Home
that replots the current state.

Moved out of ``pyside_plot.py`` in phase 5 (split of the GUI).
"""
import logging
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from PySide6.QtGui import QAction
import qtawesome as qta


class AstrocookToolbar(NavigationToolbar):
    """
    Custom Matplotlib toolbar that inherits from NavigationToolbar2QT 
    to remove unwanted buttons and simplify the interface.
    """
    # Define the order of toolbar items to keep (based on V1 functionality)
    toolitems = [
        # Label, Tooltip, Icon Filename, Tool ID
        ('Home', 'Reset original view', 'home', 'home'),
        ('Back', 'Back to previous view', 'back', 'back'),
        ('Forward', 'Forward to next view', 'forward', 'forward'),
        (None, None, None, None), # Separator placeholder
        
        # FIX: Use Matplotlib's internal symbolic names for icons (third element).
        # We rely on the base class finding 'hand' and 'magnifying_glass' equivalent icons.
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'), 
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
        
        (None, None, None, None), # Separator placeholder
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
    ]

    def __init__(self, canvas, parent, select: bool = False):
        """``select``: add the Select tool (selection of a region, main plot
        only, F73)."""
        # 1. Base Class Initialization: This CREATES and stores the C++ QAction objects.
        # This MUST run first.
        super().__init__(canvas, parent, coordinates=True) 
        self.viewer_parent = parent
        self.select_action = None

        # 2. Final cleanup (Remove unwanted actions created by the base class)
        actions_to_remove = ['Customize', 'Subplots', 'Edit axis, curve and image parameters']
        
        for action in self.actions():
            if action.text() in actions_to_remove or action.iconText() in actions_to_remove:
                self.removeAction(action)

        # 3. The Select tool, after Zoom, as a third mode exclusive with
        # Pan and Zoom
        if select:
            # FontAwesome, the family of the icons of Matplotlib (the
            # dashed rectangle of Material Design looked thinner than them)
            self.select_action = QAction(qta.icon('fa5s.arrows-alt-h', color='black'),
                                         'Select', self)
            self.select_action.setCheckable(True)
            self.select_action.setToolTip(
                "Select a region: drag along the wavelength axis, then choose "
                "what to do with it")
            self.select_action.triggered.connect(self.toggle_select_mode)
            actions = self.actions()
            after_zoom = actions[actions.index(self._actions['zoom']) + 1]
            self.insertAction(after_zoom, self.select_action)


    def deactivate_tools(self):
        """Switch off the active tool (Pan or Zoom), if any.

        Matplotlib's ``pan`` and ``zoom`` toggle their tool and ignore their
        argument: ``pan(False)`` with no tool on switches pan on (findings
        register F81).
        """
        if self.mode == 'pan/zoom':
            self.pan()
        elif self.mode == 'zoom rect':
            self.zoom()

    def toggle_select_mode(self):
        """Toggles the region selection mode on/off (the Select button)."""
        if not self.viewer_parent: return
        self.viewer_parent.toggle_region_selector()

    def sync_select_action(self):
        """Check the Select button as the selection mode is."""
        if self.select_action is not None:
            self.select_action.setChecked(bool(self.viewer_parent._is_selecting_region))

    def set_select_enabled(self, enabled: bool):
        """The Select tool is offered only on a wavelength axis (the region
        becomes an expression on λ): switched off otherwise."""
        if self.select_action is None:
            return
        if not enabled and self.viewer_parent._is_selecting_region:
            self.viewer_parent.set_region_selection(False)
        self.select_action.setEnabled(enabled)

    def _leave_select(self):
        if self.select_action is not None and self.viewer_parent._is_selecting_region:
            self.viewer_parent.set_region_selection(False)

    def pan(self, *args):
        """Pan, switching the Select tool off (the three tools exclude each
        other)."""
        self._leave_select()
        super().pan(*args)

    def zoom(self, *args):
        """Zoom, switching the Select tool off."""
        self._leave_select()
        super().zoom(*args)

    def home(self):
        """ Overrides the default 'home' to reset to the *current*
            data state, respecting normalization.
        """
        logging.debug("Custom 'Home' button pressed.")
        if self.viewer_parent and self.viewer_parent.main_window:
            session_state = None
            if self.viewer_parent.main_window.active_history:
                session_state = self.viewer_parent.main_window.active_history.current_state
            
            # Call our plot function with the new flag
            self.viewer_parent.plot_spectrum(
                session_state=session_state,
                force_autoscale=True
            )
        else:
            # Fallback to original behavior
            super().home()
