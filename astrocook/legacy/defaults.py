import json
import numpy as np
import os
import pathlib
import wx
import sys  # <--- AGGIUNTO: Necessario per PyInstaller

class Defaults(object):

    def __init__(self,
                 gui):
        self._gui = gui
        self._indent = 2

        self._extend = {"voigt": {"z": 0.0, "logN": 13, "b": 10.0, "btur": 0.0,
                        "resol": 35000}}
        self.open()

    def open(self, path=None, file='defaults.json'):

        if path is None:
            # --- FIX START ---
            # Rileva se stiamo correndo dentro PyInstaller o in sviluppo
            if getattr(sys, 'frozen', False):
                # Se "congelato" (App), la root e' la cartella temporanea _MEIPASS
                base_dir = pathlib.Path(sys._MEIPASS)
            else:
                # Se in sviluppo (astrocook/v1/defaults.py), risaliamo di 3 livelli alla root
                base_dir = pathlib.Path(__file__).resolve().parent.parent.parent
            
            # Costruiamo il percorso assoluto: root -> astrocook -> data -> file
            path = base_dir / 'astrocook' / 'data' / file
            # --- FIX END ---

        # Convertiamo path in stringa per sicurezza
        with open(str(path)) as json_file:
            self.str = json_file.read()
            self.str = self.str.replace('“', '"')
            self.str = self.str.replace('”', '"')
            self.str = self.str.replace('—', '--')
        self.dict = json.loads(self.str)
        self._dict_extend()


    def update(self, str):
        self.str = str
        self.dict = json.loads(self.str)
        self._dict_extend()

    def _dict_extend(self):
        for k1 in self._extend:
            for k2 in self._extend[k1]:
                self.dict[k1][k2] = self._extend[k1][k2]