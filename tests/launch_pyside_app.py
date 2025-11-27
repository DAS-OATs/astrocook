import argparse
import os
import logging
import sys
import time

# PySide6 Imports
from PySide6.QtWidgets import QApplication, QSplashScreen
from PySide6.QtGui import QPixmap, QColor, QPainter
from PySide6.QtCore import Qt

# Astrocook Imports
# Assuming your session and GUI modules are correctly set up on the path
from astrocook.v2.session import SessionV2, load_session_from_file
from astrocook.v2.gui.main_window import MainWindowV2 
from astrocook.v2.utils import resource_path


# --- Minimal V1/V2 Context Mocks ---
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

    # 2. Initialize PySide Application (Must be first for Splash Screen)
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)

    # --- 3. SPLASH SCREEN SETUP ---
    try:
        # Percorso del logo originale
        logo_path = resource_path(os.path.join("assets", "logo_3d_LR.png"))
        logo_pixmap = QPixmap(logo_path)

        if not logo_pixmap.isNull():
            # CONFIGURAZIONE MARGINI
            padding = 30  # 30px di spazio extra su ogni lato
            
            # Calcoliamo le nuove dimensioni totali
            new_width = logo_pixmap.width() + (padding * 2)
            new_height = logo_pixmap.height() + (padding * 2)
            
            # A. Creiamo una "tela" più grande
            splash_pixmap = QPixmap(new_width, new_height)
            splash_pixmap.fill(Qt.transparent) 

            # B. Iniziamo a dipingere
            painter = QPainter(splash_pixmap)
            painter.setRenderHint(QPainter.Antialiasing)

            # Sfondo semitrasparente (Alpha 100 come richiesto)
            bg_color = QColor(255, 255, 255, 100) 
            
            painter.setBrush(bg_color)
            painter.setPen(Qt.NoPen) 

            # Disegniamo lo sfondo su tutta la nuova area estesa
            # rect() restituisce il rettangolo 0,0,new_width,new_height
            painter.drawRoundedRect(splash_pixmap.rect(), 15, 15)

            # --- DISEGNO LOGO CENTRATO ---
            # Disegniamo il logo spostato di (padding, padding)
            painter.drawPixmap(padding, padding, logo_pixmap)
            
            painter.end() 

            # C. Creiamo lo Splash Screen
            splash = QSplashScreen(splash_pixmap)
            splash.setAttribute(Qt.WA_TranslucentBackground)
            splash.setWindowFlags(Qt.WindowStaysOnTopHint | Qt.FramelessWindowHint)
            splash.show()
            app.processEvents()
            
            time.sleep(1.5)
        else:
            logging.warning("Logo pixmap is null.")
            splash = None
            
    except Exception as e:
        logging.warning(f"Could not load splash screen: {e}")
        splash = None

    # --- 4. Load Session Logic ---
    mock_gui = MockGUI()
    initial_session = None # Default is None (triggers Empty View in MainWindow)

    if args.session_file:
        file_path = os.path.realpath(args.session_file)
        if os.path.exists(file_path):
            logging.info(f"Attempting to load session from: {file_path}")
            
            # Aggiorna messaggio sullo splash
            if splash:
                splash.showMessage(f"Loading {os.path.basename(file_path)}...", 
                                   Qt.AlignBottom | Qt.AlignCenter, Qt.white)
                app.processEvents()

            try:
                session_name = os.path.splitext(os.path.basename(file_path))[0]
                session = load_session_from_file(
                    archive_path=file_path,
                    name=session_name,
                    format_name='generic_spectrum', 
                    gui_context=mock_gui
                )
                if session != 0:
                    initial_session = session
                else:
                    logging.error(f"Failed to load initial session from {file_path}.")
            except Exception as e:
                logging.error(f"Failed to load session file {file_path}: {e}", exc_info=True)
                initial_session = None 
        else:
            logging.warning(f"File not found: {args.session_file}. Starting empty.")

    # NOTA: Se initial_session è None, passiamo None alla MainWindow.
    # La MainWindowV2 è programmata per mostrare la schermata di benvenuto in questo caso.

    # 5. Initialize the Main Window
    # initial_log_object is passed as None for now based on previous context
    main_window = MainWindowV2(initial_session, None)
    main_window.show()

    # 6. Close Splash Screen
    if splash:
        # finish() chiude lo splash solo quando main_window è visibile
        splash.finish(main_window)

    # 7. Start the PySide Event Loop
    sys.exit(app.exec())


if __name__ == '__main__':
    main()