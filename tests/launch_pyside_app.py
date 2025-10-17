import sys
import os
import logging
from PySide6.QtWidgets import QApplication

# Assuming your session and GUI modules are correctly set up on the path
# (This usually means the project root needs to be added for module resolution)

from astrocook.v2.session import SessionV2
from astrocook.v2.gui.main_window import MainWindowV2 # Import your new main window

# --- Minimal V1/V2 Context Mocks (Still needed for initialization) ---
class MockGUI:
    """A minimal mock for the V1 GUI object to satisfy V2 session initialization."""
    def __init__(self):
        self._flags = []
    def _flags_cond(self, flag):
        return False
    def _flags_extr(self, flag):
        return None
# ----------------------------------------------------------------------

def main():
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(module)s: %(message)s')
    
    # 1. Initialize PySide Application (Must be first)
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)

    # 1. Create a TRULY empty V2 session 
    mock_gui = MockGUI()
    empty_session = SessionV2(name="Astrocook V2", gui=mock_gui, current_spectrum=None) 
    
    # 2. Initialize the Main Window (Start with the empty state)
    # The MainWindowV2.__init__ must be adapted to handle the truly empty list.
    main_window = MainWindowV2(empty_session) 
    main_window.show()

    # --- DELETE/COMMENT OUT ALL DATA LOADING LOGIC ---
    # Delete the logic that attempts to load the J0103-1305 file implicitly.
    # The initial state must be empty.
    
    # 3. Start the PySide Event Loop (CRITICAL)
    sys.exit(app.exec())


if __name__ == '__main__':
    # Ensure this is runnable
    main()