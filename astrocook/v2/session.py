# astrocook/v2/session.py

from .spectrum import SpectrumV2
from ..v1.format import Format # <-- Importiamo la classe Format V1 per l'I/O reale
from .io_adapter import load_spectrum_to_v2_format # Importiamo l'adapter V2

import astropy.units as au
import numpy as np
from typing import Optional

# ... (altre importazioni, es. LineListV2, SystListV2, ecc.)

class SessionV2:
    """
    Sessione Astrocook V2: Contenitore di stato immutabile.
    """

    def __init__(self, name: str, current_spectrum: Optional[SpectrumV2] = None, 
                 lines=None, systs=None, history: list = None, **kwargs):
        self.name = name
        # I dati principali sono gestiti tramite composizione e immutabilità
        self._current_spectrum = current_spectrum
        self._lines = lines
        self._systs = systs
        self.history = history if history is not None else []
        self._gui = kwargs.get('gui') # Manteniamo il link alla GUI

        # --- Adapter Attributo per V1 (Temporaneo) ---
        # Inizializza l'attributo legacy che il codice V1 si aspetta.
        # Lo usiamo come attributo diretto per la massima compatibilità
        # con il modo in cui il codice legacy accede agli attributi interni.
        self._shade = False # Il valore predefinito sicuro per 'non ombreggiato'
        self._open_twin = kwargs.get('twin', False)
        self._clicks = []  # Used by graph._on_click to track points
        self._shade = False # Used by graph._refresh and graph._on_click
        #self._ztrans = None # Used by graph._refresh (Redshift axis)
        self._z_sel = 0.0   # Used by graph._on_syst_new and graph._on_move
        self._series_sel = 'Ly_a' # Used by graph._on_syst_new for cursor
        # ----------------------------------------------

    @property
    def spec(self) -> Optional[SpectrumV2]: # Mantiene il nome V1 'spec' per l'adapter
        return self._current_spectrum

    @property
    def lines(self):
        return self._lines # Sarà LineListV2

    # Metodo cruciale: Crea una NUOVA sessione quando si carica un file
    def open_new(self, path: str, format_name: str, **kwargs) -> 'SessionV2':
        """
        Carica un nuovo spettro utilizzando l'adapter V2 e restituisce una NUOVA SessionV2.
        (Sostituisce il vecchio 'sess.open()' che modificava in-place)
        """
        # 1. Carica il nuovo spettro V2
        new_spec_v2 = load_spectrum_to_v2_format(
            path, 
            format_name, 
            gui=self._gui,
            **kwargs
        ) 

        # 2. Aggiorna lo stato. (Gli altri oggetti sono None per ora)
        new_history = self.history + [f"Opened file {path} with {format_name}"]
        
        # 3. Restituisce una NUOVA istanza (immutabilità)
        return SessionV2(name=self.name, 
                         current_spectrum=new_spec_v2, 
                         lines=self._lines, 
                         systs=self._systs,
                         history=new_history,
                         gui=self._gui)

    # ... (altri metodi, come apply_operation, che restituiscono NUOVE sessioni)