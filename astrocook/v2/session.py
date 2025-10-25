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
from .io_v1_stubs import V1ArchiveManager, save_archive_v1
from .recipes.edit import RecipeEditV2
from .recipes.flux import RecipeFluxV2
from .spectrum import SpectrumV2
from .structures import SystemListDataV2
from .system_list import SystemListV2
from .system_list_migration import migrate_system_list_v1_to_v2
from .utils import guarded_deepcopy_v1_state

def load_session_from_file(archive_path: str, name: str, gui_context: Any, format_name: str) -> 'SessionV2':
    """
    Orchestrates the loading of a session from an archive (.acs or .acs2)
    and handles the V1-to-V2 migration if necessary.
    """
    archive_manager = V1ArchiveManager(archive_path)
    temp_dir = archive_manager.unpack()
    archive_root = ""
    v2_metadata = None
    
    if temp_dir is None:
        # Not an archive, assume it's a single FITS file
        archive_root = os.path.splitext(archive_path)[0]
        spec_file_path = archive_path # The path *is* the spec file
    else:
        # It's an archive, find the components
        spec_file_path = archive_manager.get_structure_path('spec')
        if not spec_file_path:
             raise FileNotFoundError("Could not find a _spec.fits file in the archive.")
        
        archive_root = os.path.splitext(spec_file_path)[0].replace('_spec', '')
        
        # --- V2 METADATA LOADING ---
        # V2 archive name is predictable: {base_name}_meta.json
        meta_fname = f"{os.path.basename(archive_root)}_meta.json"
        meta_file_path = os.path.join(temp_dir, meta_fname)

        if os.path.exists(meta_file_path):
            try:
                with open(meta_file_path, 'r') as f:
                    v2_metadata = json.load(f)
                logging.info("Loaded V2 metadata from _meta.json.")
            except Exception as e:
                logging.error(f"Failed to load _meta.json: {e}")
        else:
            logging.debug("No _meta.json found, assuming V1 archive.")
        # -------------------------

    # 1. Load Spectrum (This logic works for both V1 and V2 FITS files)
    spectrum_v2 = load_and_migrate_structure(
        archive_root, 'spec', gui_context, format_name, spec_file_path=spec_file_path
    )

    # 2. Load System List (Now handles V1 or V2)
    system_list_v2 = load_and_migrate_structure(
        archive_root, 'systs', gui_context, format_name,
        v2_metadata=v2_metadata  # <<< PASS THE LOADED V2 METADATA
    )

    # 3. Create the SessionV2 wrapper (V1 compatibility attributes)
    
    # Initialize V1 GUI logger (with GUI context) and Defaults
    v1_log = GUILog(gui_context)
    v1_defs = Defaults(gui_context)
    
    # This is the 'try' block you were referring to, which loads the log.
    try:
        if v2_metadata and 'log_history_json' in v2_metadata:
            # V2 path: Load log from the metadata dict
            v1_log.str = v2_metadata['log_history_json']
            logging.info("V2 log history restored.")
        elif temp_dir:
            # V1 path: Find and load the _log.json file from the archive
            log_path = archive_manager.get_structure_path('log.json')
            if not log_path:
                 # V1 log files might not follow the _suffix pattern, try basename
                 log_path_alt = os.path.join(temp_dir, f"{os.path.basename(archive_root)}_log.json")
                 if os.path.exists(log_path_alt):
                     log_path = log_path_alt

            if log_path and os.path.exists(log_path):
                with open(log_path, 'r') as f:
                    v1_log.str = f.read()
                logging.info(f"V1 log history restored from {log_path}.")
            else:
                 logging.debug(f"No V1 _log.json found (archive root: {archive_root}).")
        else:
            logging.debug("No archive or V2 metadata, starting with empty log.")
            
    except Exception as e:
        logging.error(f"Failed to restore log history: {e}")
            
    # Clean up the temporary directory
    archive_manager.cleanup()

    # Instantiate and return the new SessionV2 object
    # (We assume SessionV2 class is defined in this file)
    # Instantiate and return the new SessionV2 object
    new_session = SessionV2(
        name=name, 
        gui=gui_context, 
        log=v1_log, 
        defs=v1_defs,
        spec=spectrum_v2, 
        systs=system_list_v2
    )
    return new_session

class SessionV2:
    """
    Sessione Astrocook V2: Contenitore di stato immutabile.
    """

    def __init__(self,
                 name: str,
                 gui: Any, # Pass GUI context directly for loading
                 spec: Optional[SpectrumV2] = None, # Use 'spec' argument name matching the loader
                 systs: Optional[SystemListV2] = None,
                 log: Optional[GUILog] = None, # Optional log/defs for loading/copying
                 defs: Optional[Defaults] = None):

        self._gui = gui
        self.name = name

        # --- V1 Compatibility Adapters ---
        # Initialize log: Use passed one (e.g., from deepcopy) or create new.
        # load_session_from_file passes None here and sets .str later.
        self.log = log if log is not None else GUILog(self._gui)

        # Initialize defs: Use passed one or create new.
        self.defs = defs if defs is not None else Defaults(self._gui)

        # --- V2 Core Components ---
        # Store spectrum internally to avoid property conflict
        self._current_spectrum = spec
        # Ensure systs is always a valid SystemListV2 object
        self.systs = systs if systs is not None else SystemListV2(data=SystemListDataV2())

        # --- V2 Recipe Initialization ---
        self.edit = RecipeEditV2(self)
        self.flux = RecipeFluxV2(self)
        self.cb = self.edit

        # --- Adapter Attributo per V1 (Temporaneo) ---
        # Inizializza l'attributo legacy che il codice V1 si aspetta.
        # Lo usiamo come attributo diretto per la massima compatibilità
        # con il modo in cui il codice legacy accede agli attributi interni.
        self._shade = False # Il valore predefinito sicuro per 'non ombreggiato'
        self._open_twin = False
        self._clicks = []  # Used by graph._on_click to track points
        self._shade = False # Used by graph._refresh and graph._on_click
        self._z_sel = 0.0   # Used by graph._on_syst_new and graph._on_move
        self._series_sel = 'Ly_a' # Used by graph._on_syst_new for cursor
        # ----------------------------------------------

        if systs is None:
            # If the loading fails (systs=None), initialize with an empty SystemListDataV2 core
            empty_data_core = SystemListDataV2()
            self.systs = SystemListV2(data=empty_data_core)
        else:
            # Assume 'systs' is a valid SystemListV2 object passed from open_new
            self.systs = systs
        
        logging.debug(f"SessionV2 '{self.name}' initialized. "
                      f"Internal spec type: {type(self._current_spectrum)}, "
                      f"Property spec type: {type(self.spec)}") # Check both

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
    @classmethod
    def open_new(cls, file_path: str, name: str, gui_context: Any, format_name: str) -> 'SessionV2':        
        """
        Carica un nuovo spettro utilizzando l'adapter V2 e restituisce una NUOVA SessionV2.
        (Sostituisce il vecchio 'sess.open()' che modificava in-place)
        """
        try:
            # 1. Call the module-level orchestrator function
            session = load_session_from_file(
                archive_path=file_path,  # <<< ARGUMENT RENAMED
                name=name, 
                gui_context=gui_context, 
                format_name=format_name
            )
            logging.info(f"SessionV2 '{name}' loaded successfully from {file_path}.")
            return session
        except Exception as e:
            # Handle I/O failures gracefully
            logging.error(f"FATAL I/O ERROR: Failed to load and orchestrate session from {file_path}: {e}")
            return 0 # Return 0 to signal failure to the V1 dialog loop

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
    
    def save(self, file_path: str, models: bool = False):
        """
        Saves the current session state.
        Calls V2 saver (.acs2) or V1 saver (.acs) based on file extension.
        """
        
        # Determine which archive format to use
        if file_path.lower().endswith('.acs2'):
            # --- V2 NATIVE SAVE ---
            try:
                # Call the V2 archive writer, passing the SessionV2 instance
                save_archive_v2(self, file_path)
                logging.info(f"V2 Session state saved successfully to {file_path}.")
                
            except Exception as e:
                logging.error(f"FATAL: V2 archive saving failed: {e}")
                return 0 # V1 failure code
                
        else:
            # --- V1 LEGACY SAVE (.acs) ---
            
            # 1. Convert V2 immutable structures to V1 mutable/saveable format
            try:
                v1_spec_for_save = self.spec.to_v1_spectrum() 
                v1_systs_for_save = self.systs.to_v1_systlist() 

            except Exception as e:
                logging.error(f"FATAL: V2-to-V1 conversion failed during saving: {e}")
                return 0
            
            # 2. Get the V1-compatible JSON log string
            # save_archive_v1 expects a string, not a dict
            json_log_str = self.log.str 
            
            # 3. Call V1 Archive Writer Utility
            try:
                save_archive_v1(
                    v1_spec_for_save, 
                    v1_systs_for_save, 
                    json_log_str,  # CRITICAL FIX: Pass the correct JSON string
                    file_path
                )
                logging.info(f"V1 compatibility session saved successfully to {file_path}.")
            
            except Exception as e:
                logging.error(f"FATAL: V1 archive writing failed: {e}")
                return 0
        
        return 0 # V1 success code