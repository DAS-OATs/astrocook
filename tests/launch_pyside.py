import sys
import os
import logging

# Ensure the root Astrocook package is on the path (standard practice for external test scripts)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from astrocook.v2.session import SessionV2
from astrocook.v2.gui.pyside_plot import launch_viewer

# --- Minimal V1/V2 Context Mocks ---
# These mocks satisfy the minimal requirements of the V2 Session constructor
class MockGUI:
    """A minimal mock for the V1 GUI object."""
    def __init__(self):
        self._flags = []
    def _flags_cond(self, flag):
        return False
    def _flags_extr(self, flag):
        return None

def main():
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(module)s: %(message)s')
    
    # 1. Define Test Parameters
    TEST_FILE_PATH = os.path.join(os.path.dirname(__file__), 'data', 'J0103-1305_metallines_lya_new_spec.fits')
    if not os.path.exists(TEST_FILE_PATH):
        logging.error(f"Test file not found at: {TEST_FILE_PATH}")
        sys.exit(1)
        
    FORMAT_NAME = 'generic_spectrum' # Name of the V1 format loader to use

    # 2. Create Initial Session and Context
    mock_gui = MockGUI()
    initial_sess = SessionV2(name="PySide_Test", gui=mock_gui) 

    # 3. Execute V2 Immutable Loading (Bypassing V1 GUI dialogs)
    try:
        # Load the final, fully-populated V2 session object
        # NOTE: We rely on the SessionV2.open_new method to handle all I/O
        sess_loaded = initial_sess.open_new(
            path=TEST_FILE_PATH, 
            format_name=FORMAT_NAME,
            #gui=mock_gui # Pass context to the I/O adapter chain
        )
        logging.info("V2 Session loaded successfully. Launching PySide viewer...")

    except Exception as e:
        logging.error(f"Fatal error during V2 session loading: {e}")
        # Display the full error for immediate debugging
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # 4. Launch the PySide Viewer
    launch_viewer(sess_loaded)

if __name__ == '__main__':
    main()