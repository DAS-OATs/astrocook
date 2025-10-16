# astrocook/v2/session.py

from ..v1.defaults import Defaults
from ..v1.format import Format # Import the Format V1 class for I/O 
from ..v1.gui_log import GUILog # Import the V1 logger for GUI compatibility
from .io_adapter import load_spectrum_to_v2_format # Import V2 adapter for loading
from .recipes.edit import RecipeEditV2
from .spectrum import SpectrumV2
from .utils import guarded_deepcopy_v1_state

import astropy.units as au
from copy import deepcopy
import numpy as np
from typing import Optional


class SessionV2:
    """
    Sessione Astrocook V2: Contenitore di stato immutabile.
    """

    def __init__(self, name: str, current_spectrum: Optional[SpectrumV2] = None, 
                 lines=None, systs=None, log=None, history: list = None, **kwargs):
        self.name = name
        # I dati principali sono gestiti tramite composizione e immutabilità
        self._current_spectrum = current_spectrum
        self._lines = lines
        self._systs = systs
        self.history = history if history is not None else []
        self._gui = kwargs.get('gui') # Manteniamo il link alla GUI
        
        # The log needs to be attached to the session for V1 compatibility.
        # We deepcopy the GUILog object if passed during copying.
        log_instance = kwargs.get('log')
        if log_instance is not None:
            self.log = deepcopy(log_instance)
        else:
            # Initialize a new GUILog instance (requires the GUI object)
            self.log = GUILog(self._gui) 

        defs_instance = kwargs.get('defs')
        if defs_instance is not None:
            # When creating an immutable copy, use the deepcopied defs passed
            self.defs = defs_instance 
        else:
            # When creating a brand new session, initialize a new Defaults instance
            self.defs = Defaults(self._gui)

        # The V1 GUI expects the Cookbook to be named 'cb'
        self._cb = RecipeEditV2(self)

        # --- Adapter Attributo per V1 (Temporaneo) ---
        # Inizializza l'attributo legacy che il codice V1 si aspetta.
        # Lo usiamo come attributo diretto per la massima compatibilità
        # con il modo in cui il codice legacy accede agli attributi interni.
        self._shade = False # Il valore predefinito sicuro per 'non ombreggiato'
        self._open_twin = kwargs.get('twin', False)
        self._clicks = []  # Used by graph._on_click to track points
        self._shade = False # Used by graph._refresh and graph._on_click
        self._z_sel = 0.0   # Used by graph._on_syst_new and graph._on_move
        self._series_sel = 'Ly_a' # Used by graph._on_syst_new for cursor
        # ----------------------------------------------

    @property
    def spec(self) -> Optional[SpectrumV2]: # Mantiene il nome V1 'spec' per l'adapter
        return self._current_spectrum

    @property
    def lines(self):
        return self._lines # Sarà LineListV2
    
    @property
    def cb(self) -> RecipeEditV2:
        return self._cb

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

    def with_new_spectrum(self, new_spec_v2: SpectrumV2) -> 'SessionV2':
        """
        Helper method to create a NEW SessionV2 instance with the updated spectrum,
        preserving all other state (immutability pattern).
        """

        # 1. Use the new utility for both log and defs
        copied_log = guarded_deepcopy_v1_state(self.log)
        copied_defs = guarded_deepcopy_v1_state(self.defs)

        # 1. Mutate the V2 constructor parameters
        kwargs = {
            'name': self.name,
            'current_spectrum': new_spec_v2,
            
            # Preserve all V1 adapter attributes and essential state
            'lines': deepcopy(self._lines),
            'systs': deepcopy(self._systs),
            'history': new_spec_v2.history, # Use the updated history from the spectrum
            'gui': self._gui,
            'log': copied_log, 
            'defs': copied_defs, 

            # --- V1 ADAPTER ATTRIBUTES (Must be copied to maintain state) ---
            'twin': self._open_twin, 
            '_shade': self._shade,
            '_z_sel': self._z_sel,
            '_series_sel': self._series_sel,
            # ... and any other internal attributes the V1 code might rely on!
            # The simplest and safest approach is to manually pass the adapters like this.
        }
        
        # 2. Return the new, copied instance
        return SessionV2(**kwargs)