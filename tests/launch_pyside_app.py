import argparse
import os
import logging
import sys
from PySide6.QtWidgets import QApplication

# Assuming your session and GUI modules are correctly set up on the path
# (This usually means the project root needs to be added for module resolution)

from astrocook.v2.session import SessionV2, load_session_from_file
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
    
    # 1. Set up Argument Parser
    parser = argparse.ArgumentParser(description="Launch Astrocook V2 PySide GUI.")
    parser.add_argument('session_file', nargs='?', default=None,
                        help="Optional path to a .acs or .acs2 session file to load.")
    args = parser.parse_args()

    # 2. Initialize PySide Application (Must be first)
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)

    # 3. Create a TRULY empty V2 session
    mock_gui = MockGUI()
    empty_session = SessionV2(name="Astrocook V2", gui=mock_gui, spec=None)

    # 4. Load Session or Create Empty
    session_to_load = None
    if args.session_file:
        file_path = os.path.realpath(args.session_file)
        if os.path.exists(file_path):
            logging.info(f"Attempting to load session from: {file_path}")
            try:
                session_name = os.path.splitext(os.path.basename(file_path))[0]
                session, log_str = load_session_from_file(
                    archive_path=file_path,
                    name=session_name,
                    format_name='generic_spectrum', # V1 loader default
                    gui_context=mock_gui
                )
                if session != 0:
                    initial_session = session
                    initial_log_string = log_str
                else:
                    logging.error(f"Failed to load initial session from {file_path}.")
            except Exception as e:
                logging.error(f"Failed to load session file {file_path}: {e}", exc_info=True)
                initial_session = None # Fallback to empty
        else:
            logging.warning(f"File not found: {args.session_file}. Starting empty.")

    if initial_session is None:
        # No file provided, or loading failed.
        logging.info("Starting with an empty session.")
        initial_session = SessionV2(name="Astrocook V2", gui=mock_gui, spec=None)

    # 5. Initialize the Main Window (Start with the empty state)
    # The MainWindowV2.__init__ must be adapted to handle the truly empty list.
    main_window = MainWindowV2(initial_session, initial_log_string)
    main_window.show()

    # 6. Start the PySide Event Loop (CRITICAL)
    sys.exit(app.exec())


if __name__ == '__main__':
    # Ensure this is runnable
    main()