import sys
from PySide6.QtWidgets import QApplication
from astrocook.gui.main_window import MainWindowV2

def main():
    app = QApplication(sys.argv)
    window = MainWindowV2(initial_session=None, initial_log_object=None)
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()