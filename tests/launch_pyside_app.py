import argparse
import os
import logging
import sys
import time

# PySide6 Imports
from PySide6.QtWidgets import QApplication, QSplashScreen
from PySide6.QtGui import QPixmap, QColor, QPainter
from PySide6.QtCore import Qt, QEvent, Signal

# Astrocook Imports
# Assuming your session and GUI modules are correctly set up on the path
from astrocook.core.session import SessionV2, load_session_from_file
from astrocook.gui.main_window import MainWindowV2 
from astrocook.core.utils import resource_path


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
    # Segnale emesso quando macOS chiede di aprire un file a runtime
    file_open_requested = Signal(str)

    def __init__(self, argv):
        super().__init__(argv)
        self.startup_file = None # Per salvare il file se arriva durante l'avvio

    def event(self, event):
        # Intercettiamo l'evento specifico di macOS
        if event.type() == QEvent.FileOpen:
            file_path = event.file()
            logging.info(f"MacOS FileOpen Event received: {file_path}")
            
            # Se la main window non è ancora pronta, salviamo il path per dopo
            self.startup_file = file_path
            
            # Emettiamo comunque il segnale (utile se l'app è già aperta)
            self.file_open_requested.emit(file_path)
            return True
            
        return super().event(event)
    
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
        app = AstroCookApp(sys.argv)

    # --- 3. SPLASH SCREEN SETUP ---
    try:
        # Percorso del logo originale
        logo_path = resource_path(os.path.join("assets", "logo_3d_LR.png"))
        logo_pixmap = QPixmap(logo_path)

        if not logo_pixmap.isNull():
            # CONFIGURAZIONE MARGINI
            hpadding = 80
            vpadding = 40
            
            # Calcoliamo le nuove dimensioni totali
            new_width = logo_pixmap.width() + (hpadding * 2)
            new_height = logo_pixmap.height() + (vpadding * 2)
            
            # A. Creiamo una "tela" più grande
            splash_pixmap = QPixmap(new_width, new_height)
            splash_pixmap.fill(Qt.transparent) 

            # B. Iniziamo a dipingere
            painter = QPainter(splash_pixmap)
            painter.setRenderHint(QPainter.Antialiasing)

            # Sfondo semitrasparente (Alpha 100 come richiesto)
            bg_color = QColor(255, 255, 255, 150) 
            
            painter.setBrush(bg_color)
            painter.setPen(Qt.NoPen) 

            # Disegniamo lo sfondo su tutta la nuova area estesa
            # rect() restituisce il rettangolo 0,0,new_width,new_height
            painter.drawRoundedRect(splash_pixmap.rect(), 15, 15)

            # --- DISEGNO LOGO CENTRATO ---
            # Disegniamo il logo spostato di (padding, padding)
            painter.drawPixmap(hpadding, vpadding, logo_pixmap)
            
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
    file_to_load = None

    # A. Priorità: Argomento da riga di comando (Windows/Linux/Debug)
    if args.session_file:
        file_to_load = os.path.realpath(args.session_file)
    
    # B. Se vuoto, controlliamo se macOS ha inviato un evento durante lo splash
    elif hasattr(app, 'startup_file') and app.startup_file:
        file_to_load = app.startup_file

    # Caricamento effettivo
    if file_to_load and os.path.exists(file_to_load):
        if splash:
            splash.showMessage(f"Loading {os.path.basename(file_to_load)}...", 
                               Qt.AlignBottom | Qt.AlignCenter, Qt.white)
            app.processEvents()
        try:
            name = os.path.splitext(os.path.basename(file_to_load))[0]
            
            # [FIX] Auto-detect format based on extension for CLI loading
            if file_to_load.lower().endswith(('.txt', '.dat')):
                format_name = 'ascii_resvel_header'
            else:
                format_name = 'generic_spectrum'

            # Pass the detected format_name instead of hardcoded 'generic_spectrum'
            initial_session = load_session_from_file(file_to_load, name, mock_gui, format_name)
        except Exception as e:
            logging.error(f"Error loading {file_to_load}: {e}")

    # 3. AVVIO GUI
    main_window = MainWindowV2(initial_session, None)

    # Connecting the signal to open files when requested by macOS
    app.file_open_requested.connect(lambda path: main_window.open_session_from_path(path))

    main_window.show()

    # 6. Close Splash Screen
    if splash:
        # finish() chiude lo splash solo quando main_window è visibile
        splash.finish(main_window)

    # 7. Start the PySide Event Loop
    sys.exit(app.exec())


if __name__ == '__main__':
    main()