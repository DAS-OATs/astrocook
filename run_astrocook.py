import argparse
import glob
import logging
import os
from pathlib import Path
import sys
import time

# PySide6 Imports
from PySide6.QtWidgets import QApplication, QSplashScreen
from PySide6.QtGui import QPixmap, QColor, QPainter, QIcon
from PySide6.QtCore import Qt, QEvent, Signal

def resource_path(relative_path):
    if getattr(sys, 'frozen', False):
        base_path = Path(sys._MEIPASS)
    else:
        base_path = Path(__file__).resolve().parent
    return str(base_path / relative_path)

# Minimal V1/V2 Context Mocks
class MockGUI:
    """A minimal mock for the V1 GUI object to satisfy V2 session initialization."""
    def __init__(self):
        self._flags = []
    def _flags_cond(self, flag):
        return False
    def _flags_extr(self, flag):
        return None

# Class to handle macOS file open events
class AstroCookApp(QApplication):
    # Signal emitted when macOS asks to open a file at runtime
    file_open_requested = Signal(str)

    def __init__(self, argv):
        super().__init__(argv)
        self.startup_files = [] # Change to list to handle multiple files

    def event(self, event):
        # Intercept the specific macOS event
        if event.type() == QEvent.FileOpen:
            file_path = event.file()
            logging.info(f"MacOS FileOpen Event received: {file_path}")
            
            # If the main window isn't ready, queue it
            if not self.activeWindow():
                self.startup_files.append(file_path)
            
            # Emit the signal anyway (useful if the app is already open)
            self.file_open_requested.emit(file_path)
            return True
            
        return super().event(event)
    
def main():
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(module)s: %(message)s')
    
    # 1. Set up Argument Parser
    parser = argparse.ArgumentParser(description="Launch Astrocook V2 PySide GUI.")
    parser.add_argument('session_files', nargs='*', default=[],
                        help="Optional paths to .acs, .acs2, or FITS files.")
    args = parser.parse_args()

    # 2. Initialize PySide Application (Must be first for Splash Screen)
    app = QApplication.instance()
    if not app:
        app = AstroCookApp(sys.argv)

    is_frozen = getattr(sys, 'frozen', False)

    if is_frozen:
        # App mode
        # Use the official (square) icon for windows
        icon_name = "icon_3d_HR.icns"
    else:
        # Script mode (Development)
        # Use the round icon to distinguish it immediately
        icon_name = "icon_3d_HR_round.png"

    try:
        app_icon_path = resource_path(os.path.join("assets", icon_name))
        app_icon = QIcon(app_icon_path)
        app.setWindowIcon(app_icon)
        
        if sys.platform == 'darwin' and not is_frozen:
            try:
                from AppKit import NSApplication, NSImage
                
                # Calcoliamo il percorso assoluto per Cocoa
                abs_path = os.path.abspath(app_icon_path)
                image = NSImage.alloc().initWithContentsOfFile_(abs_path)
                
                if image:
                    NSApplication.sharedApplication().setApplicationIconImage_(image)
                    logging.info(f"Dev Dock Icon set to: {icon_name}")
            except ImportError:
                pass
            except Exception as e:
                logging.warning(f"Failed to set Dock icon: {e}")

    except Exception as e:
        logging.warning(f"Could not set app icon: {e}")

    import time
    min_splash_duration = 1.0  # Durata minima in secondi
    start_time = time.time()

    # 3. Setup Splash Screen
    splash = None
    try:
        logo_path = resource_path(os.path.join("assets", "logo_3d_LR.png"))
        logo_pixmap = QPixmap(logo_path)

        if not logo_pixmap.isNull():
            hpadding = 80
            vpadding = 40
            
            new_width = logo_pixmap.width() + (hpadding * 2)
            new_height = logo_pixmap.height() + (vpadding * 2)
            
            splash_pixmap = QPixmap(new_width, new_height)
            splash_pixmap.fill(Qt.transparent) 

            painter = QPainter(splash_pixmap)
            painter.setRenderHint(QPainter.Antialiasing)

            bg_color = QColor(255, 255, 255, 150) 
            
            painter.setBrush(bg_color)
            painter.setPen(Qt.NoPen) 

            painter.drawRoundedRect(splash_pixmap.rect(), 15, 15)

            painter.drawPixmap(hpadding, vpadding, logo_pixmap)
            
            painter.end() 

            splash = QSplashScreen(splash_pixmap)
            splash.setAttribute(Qt.WA_TranslucentBackground)
            splash.setWindowFlags(Qt.WindowStaysOnTopHint | Qt.FramelessWindowHint)
            splash.show()
            app.processEvents()
            
            time.sleep(1.5)
        else:
            logging.warning("Logo pixmap is null.")
            
    except Exception as e:
        logging.warning(f"Could not load splash screen: {e}")

    logging.info("Importing heavy modules...")
    
    from astrocook.core.session import SessionV2, load_session_from_file
    from astrocook.gui.main_window import MainWindowV2 
    
    elapsed = time.time() - start_time
    if elapsed < min_splash_duration:
        # If it was too fast, let's wait
        time.sleep(min_splash_duration - elapsed)

    # 5. Load Session Logic
    mock_gui = MockGUI()
    loaded_sessions = [] # List to hold all loaded sessions

    # A. Gather file paths from CLI (handling wildcards explicitly for Windows)
    files_to_load = []

    if args.session_files:
        for pattern in args.session_files:
            # Glob expands wildcards (e.g. *.fits) if the shell didn't already
            matched = glob.glob(os.path.expanduser(pattern))
            if not matched:
                # If no match (maybe it's a specific file that doesn't exist yet?), keep literal
                files_to_load.append(pattern)
            else:
                files_to_load.extend(matched)
    
    # B. Check for macOS startup event
    elif hasattr(app, 'startup_files') and app.startup_files:
        files_to_load.extend(app.startup_files)

    # C. Load loop
    for file_path in files_to_load:
        real_path = os.path.realpath(file_path)
        if os.path.exists(real_path):
            if splash:
                splash.showMessage(f"Loading {os.path.basename(real_path)}...", 
                                   Qt.AlignBottom | Qt.AlignCenter, Qt.white)
                app.processEvents()
            try:
                name = os.path.splitext(os.path.basename(real_path))[0]
                
                # Note: We now use 'auto' as the standard thanks to your previous fix
                new_sess = load_session_from_file(real_path, name, mock_gui, format_name='auto')
                
                if new_sess and new_sess != 0:
                    loaded_sessions.append(new_sess)
                    logging.info(f"Successfully loaded: {name}")
                else:
                    logging.error(f"Failed to load: {real_path}")

            except Exception as e:
                logging.error(f"Error loading {real_path}: {e}")

    # 6. GUI Startup
    if splash:
        splash.showMessage("Initializing Interface...", Qt.AlignBottom | Qt.AlignCenter, Qt.white)
        app.processEvents()

    # Pass the LIST of sessions to the main window
    main_window = MainWindowV2(loaded_sessions, None)

    # Connecting the signal to open files when requested by macOS
    app.file_open_requested.connect(lambda path: main_window.open_session_from_path(path))

    main_window.show()

    # 7. Close Splash Screen
    if splash:
        # finish() closes the splash only when main_window is visible
        splash.finish(main_window)

    # 8. Start the PySide Event Loop
    sys.exit(app.exec())


if __name__ == '__main__':
    main()