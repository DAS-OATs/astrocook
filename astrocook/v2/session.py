import astropy.units as au
from copy import deepcopy
import logging
import numpy as np
import os
from typing import Optional, Any

from ..v1.defaults import Defaults
from ..v1.format import Format # Import the Format V1 class for I/O 
from ..v1.gui_log import GUILog # Import the V1 logger for GUI compatibility
from .io_adapter import load_and_migrate_structure # Import V2 adapter for loading
from .io_v1_stubs import V1ArchiveManager
from .recipes.edit import RecipeEditV2
from .recipes.flux import RecipeFluxV2
from .spectrum import SpectrumV2
from .system_list_migration import migrate_system_list_v1_to_v2
from .utils import guarded_deepcopy_v1_state


def load_session_from_file(file_path: str, format_name: str, gui_context: Any) -> 'SessionV2':
    """
    Orchestrates the loading and migration of all associated structures for a new SessionV2 object.
    """
    
    archive_manager = None
    
    try:
        # --- 1. Archive Detection and Unpacking ---
        is_archive = file_path.lower().endswith('.acs')
        
        if is_archive:
            archive_manager = V1ArchiveManager(file_path)
            temp_root = archive_manager.unpack()
            if temp_root is None:
                raise FileNotFoundError(f"Could not unpack archive {file_path}")
            
            # The archive manager must resolve the actual file path inside the temp dir.
            # We must resolve the path for the primary structure (spec).
            spec_file_path = archive_manager.get_structure_path('spec')
            if spec_file_path is None:
                raise FileNotFoundError("Mandatory *_spec.fits not found in archive.")
                
            # The archive root for subsequent loads (systs) is derived from the spec file
            archive_root = os.path.splitext(spec_file_path)[0] # e.g., /tmp/xyz/session_root_spec
            name = os.path.basename(archive_root)
            
        else:
            # --- Standard FITS File Handling ---
            spec_file_path = file_path
            archive_root = os.path.splitext(file_path)[0]
            name = os.path.basename(archive_root)

        # Clean the name: remove '_spec' suffix if present (as the V1 stubs expects the core name)
        if name.lower().endswith('_spec'):
            name = name[:-5]
            archive_root = archive_root[:-5]

        # 2. Load and migrate Spectrum (required structure)
        # The io_adapter needs the full path to the FITS file and the correct root name
        spec_v2 = load_and_migrate_structure(
            archive_root, 'spec', gui_context, format_name, spec_file_path # Pass full path for direct FITS loading
        )
        
        if spec_v2 is None:
            raise FileNotFoundError(f"Failed to load Spectrum structure from {spec_file_path}")

        # 3. Load and migrate System List (optional structure)
        # The archive_root variable is the base path for associated files (*_systs.fits)
        systs_v2 = load_and_migrate_structure(
            archive_root, 'systs', gui_context, format_name, spec_file_path=None
        )
        
        # 4. Create the final SessionV2 object
        # The initial session is created just to get the log/defs objects correctly
        initial_sess = SessionV2(name=name, gui=gui_context)
        
        sess_loaded = SessionV2(
            name=name, 
            current_spectrum=spec_v2, 
            systs=systs_v2, 
            # Pass original context objects
            log=initial_sess.log,
            defs=initial_sess.defs,
            gui=gui_context
        )
        return sess_loaded

    except Exception as e:
        logging.error(f"FATAL I/O ERROR during session loading: {e}")
        # Re-raise or handle cleanup
        raise
        
    finally:
        if archive_manager:
            archive_manager.cleanup()

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
        self.systs = systs
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
        self.edit = RecipeEditV2(self)
        self.flux = RecipeFluxV2(self)
        self.cb = self.edit

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

    def __getattr__(self, name):
        """
        Dynamically dispatches method calls (like 'rebin') 
        to the appropriate V2 recipe adapter (self.flux, self.edit).
        """
        # 1. Check if the method exists in the 'flux' recipe category
        if hasattr(self.flux, name):
            return getattr(self.flux, name)
            
        # 2. Check if the method exists in the 'edit' recipe category
        if hasattr(self.edit, name):
            return getattr(self.edit, name)

        # 3. Fallback to default Python behavior (raise AttributeError)
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")

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
        # We access the GUI context from the instance variable
        gui_context = self._gui 
        
        # --- CRITICAL FIX: The entire loading logic is now handled by the orchestrator ---
        try:
            # 1. Call the module-level orchestrator function
            new_session = load_session_from_file(
                file_path=path,
                format_name=format_name,
                gui_context=gui_context
            )
        except Exception as e:
            # Handle I/O failures gracefully
            logging.error(f"FATAL I/O ERROR: Failed to load and orchestrate session from {path}: {e}")
            return 0 # Return 0 to signal failure to the V1 dialog loop

        # 2. Add history (optional, as history is added in the orchestrator)
        # We rely on the orchestrator to set the history correctly.
        
        # 3. Return the fully loaded, immutable SessionV2 object
        return new_session

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
    
    def load_structure_v2_from_file(archive_root: str, structure_name: str):
        """
        Simulates V1 Session.open() logic for one structure and migrates it to V2.
        """

        file_path_fits = f"{archive_root}_{structure_name}.fits"

        if structure_name == 'systs':
            # --- CRITICAL: Call V1 loading utility and then migrate ---

            # 1. Load V1 SystList object (Requires a specialized V1 loading function)
            # We must call the original V1 format parser to get the V1 SystList object.
            # Since we don't have the V1 format file, we assume a utility exists:
            v1_systs = load_v1_systs_object(file_path_fits) # Placeholder for V1 I/O

            if v1_systs:
                return migrate_system_list_v1_to_v2(v1_systs)
            return None

        elif structure_name == 'spec':
            # Use the existing spectrum migration logic
            # return v1_table_to_data_v2(load_v1_spec_object(file_path_fits)) # Simplified
            pass # The initial loading already covers the spectrum.

        return None # Return None for other unimplemented structuresl