"""The toolbar of the plots (main plot, velocity plot, corner plot): the
matplotlib toolbar without the buttons Astrocook does not use, and a Home
that replots the current state.

Moved out of ``pyside_plot.py`` in phase 5 (split of the GUI).
"""
import logging
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
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

    def __init__(self, canvas, parent):
        # 1. Base Class Initialization: This CREATES and stores the C++ QAction objects.
        # This MUST run first.
        super().__init__(canvas, parent, coordinates=True) 
        self.viewer_parent = parent

        # 2. Final cleanup (Remove unwanted actions created by the base class)
        actions_to_remove = ['Customize', 'Subplots', 'Edit axis, curve and image parameters']
        
        for action in self.actions():
            if action.text() in actions_to_remove or action.iconText() in actions_to_remove:
                self.removeAction(action)

        for action in self.actions():
            if action.text() == 'Select':
                # Use a standard Qt icon. 
                # SP_FileDialogDetailedView often looks like a list/selection.
                # You can try others like SP_ArrowRight or SP_ToolBarHorizontalExtensionButton
                icon = qta.icon('mdi.selection-drag', color='black')
                action.setIcon(icon)
                break
            
    def toggle_select_mode(self):
        """Toggles the region selection mode on/off."""
        if not self.viewer_parent: return
        
        # Toggle the mode in the parent widget
        self.viewer_parent.toggle_region_selector()
        
        # Visual feedback: "press" or "unpress" the button based on new state
        # (We need to find our action in the toolbar list)
        for action in self.actions():
             if action.text() == 'Select':
                 action.setChecked(self.viewer_parent._is_selecting_region)
                 break

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
