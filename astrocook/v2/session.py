import astropy.units as au
from copy import deepcopy
import json
import logging
import numpy as np
import os
import tarfile  # <<< Import
import tempfile # <<< Import
from typing import Any, Optional, Tuple, Union

from ..v1.defaults import Defaults
from ..v1.format import Format # Import the Format V1 class for I/O 
from ..v1.gui_log import GUILog # Import the V1 logger for GUI compatibility
from .io.adapter import load_and_migrate_structure, save_archive_v2
from .io.loaders import get_loader
from .io.v1_stubs import V1ArchiveManager, save_archive_v1
from .recipes.continuum import RecipeContinuumV2
from .recipes.absorbers import RecipeAbsorbersV2
from .recipes.edit import RecipeEditV2
from .recipes.flux import RecipeFluxV2
from .spectrum import SpectrumV2
from .structures import (
    SystemListDataV2, HistoryLogV2, LogEntryV2, V1LogArtifact
)
from .system_list import SystemListV2
from .system_list_migration import migrate_system_list_v1_to_v2
from .utils import guarded_deepcopy_v1_state

# Define the type for the log manager
LogManager = Union[HistoryLogV2, V1LogArtifact, GUILog]

def load_session_from_file(archive_path: str, name: str, gui_context: Any, format_name: str) -> Optional['SessionV2']:
    """
    Orchestrates the loading of a session from an archive (.acs or .acs2)
    and handles the V1-to-V2 migration if necessary.
    
    Per Point 3, this function now *ignores* any saved log.
    """
    archive_manager = None 
    temp_dir = None
    archive_root = ""
    v2_metadata = None
    # log_object_to_return has been removed
    
    try:
        archive_path_lower = archive_path.lower()
        
        if archive_path_lower.endswith('.acs2') or archive_path_lower.endswith('.tar.gz'):
            # ... (unpacking logic) ...
            logging.debug(f"Unpacking V2 archive: {archive_path}")
            temp_dir = tempfile.mkdtemp()
            with tarfile.open(archive_path, 'r:gz') as tar:
                tar.extractall(path=temp_dir)
            
            spec_file_path = None
            for f in os.listdir(temp_dir):
                if f.endswith('_spec.fits'):
                    spec_file_path = os.path.join(temp_dir, f)
                    break
            if not spec_file_path:
                 raise FileNotFoundError("Could not find a _spec.fits file in the .acs2 archive.")
            
            archive_root = os.path.splitext(spec_file_path)[0].replace('_spec', '')

        elif archive_path_lower.endswith('.acs'):
            # ... (unpacking logic) ...
            logging.debug(f"Unpacking V1 archive: {archive_path}")
            archive_manager = V1ArchiveManager(archive_path)
            temp_dir = archive_manager.unpack()
            if not temp_dir:
                 raise RuntimeError("V1ArchiveManager failed to unpack .acs file.")

            spec_file_path = archive_manager.get_structure_path('spec')
            if not spec_file_path:
                 raise FileNotFoundError("Could not find a _spec.fits file in the .acs archive.")
            archive_root = os.path.splitext(spec_file_path)[0].replace('_spec', '')

        else:
            # ... (single FITS logic) ...
            logging.debug("Loading single FITS file.")
            temp_dir = None
            archive_root = os.path.splitext(archive_path)[0]
            spec_file_path = archive_path
        
        # --- V2 METADATA LOADING (Still needed for constraints) ---
        if temp_dir:
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
                logging.debug("No _meta.json found (assuming V1 archive or single FITS).")

        # 1. Check if a V2 Native Loader exists for this format
        v2_loader = get_loader(format_name)
    
        if v2_loader:
            logging.info(f"Using V2 Native Loader for format '{format_name}'")
            try:
                # Call the pure loader
                spec_data = v2_loader(spec_file_path)

                # Wrap in API object
                spectrum_v2 = SpectrumV2(data=spec_data)

                # Create session (SystemList will be empty initially for raw loads)
                return SessionV2(name=name, gui=gui_context, spec=spectrum_v2)

            except Exception as e:
                logging.error(f"V2 Loader failed: {e}")
                return 0

        # 1. Load Spectrum
        spectrum_v2 = load_and_migrate_structure(
            archive_root, 'spec', gui_context, format_name, 
            spec_file_path=spec_file_path,
            v2_metadata=v2_metadata 
        )

        # 2. Load System List
        system_list_v2 = load_and_migrate_structure(
            archive_root, 'systs', gui_context, format_name,
            v2_metadata=v2_metadata
        )

        if not spectrum_v2: 
            raise RuntimeError("Spectrum loading failed, aborting session creation.")

        # --- *** (Point 3) *** ---
        # 3. Load V2/V1 Log Object
        #    This entire block is now REMOVED. We no longer parse the log.
        logging.debug("Ignoring saved log on load, as per new design.")
            
        
        # Call the simplified constructor
        new_session = SessionV2(
            name=name, 
            gui=gui_context, 
            spec=spectrum_v2, 
            systs=system_list_v2
        )
        
        return new_session # <<< *** MODIFIED: Only return session

    except Exception as e: 
        logging.error(f"FATAL: load_session_from_file failed: {e}", exc_info=True)
        return 0 # <<< *** MODIFIED: Return 0 on failure (V1 style)

    finally: 
        # ... (cleanup logic) ...
        if archive_manager:
            archive_manager.cleanup()
        elif temp_dir and os.path.exists(temp_dir):
            import shutil
            shutil.rmtree(temp_dir)
            
# --- SessionV2 Class ---
# (The class definition from the previous step is correct)
class SessionV2:
    def __init__(self,
                 name: str,
                 gui: Any, 
                 spec: Optional[SpectrumV2] = None, 
                 systs: Optional[SystemListV2] = None):

        self._gui = gui
        self.name = name
        self.log: Optional[GUILog] = GUILog(self._gui)
        self.defs: Optional[Defaults] = Defaults(self._gui)
        self.log_manager: Optional[LogManager] = None
        self._current_spectrum = spec
        self.systs = systs if systs is not None else SystemListV2(data=SystemListDataV2())
        
        self.edit = RecipeEditV2(self)
        self.flux = RecipeFluxV2(self)
        self.continuum = RecipeContinuumV2(self)
        self.absorbers = RecipeAbsorbersV2(self)

        self.cb = self.edit
        self._shade = False 
        self._open_twin = False
        self._clicks = []  
        self._z_sel = 0.0   
        self._series_sel = 'Ly_a' 
        logging.debug(f"SessionV2 '{self.name}' initialized.")

    def __getattr__(self, name):
        if hasattr(self.flux, name):
            return getattr(self.flux, name)
        if hasattr(self.edit, name):
            return getattr(self.edit, name)
        if hasattr(self.continuum, name):
            return getattr(self.continuum, name)
        if hasattr(self.absorbers, name):
            return getattr(self.absorbers, name)
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")

    @property
    def spec(self) -> Optional[SpectrumV2]: 
        return self._current_spectrum
    @property
    def lines(self):
        return self._lines 

    @classmethod
    def open_new(cls, file_path: str, name: str, gui_context: Any, format_name: str) -> 'SessionV2':        
        try:
            session = load_session_from_file(
                archive_path=file_path,
                name=name, 
                gui_context=gui_context, 
                format_name=format_name
            )
            return session
        except Exception as e:
            logging.error(f"FATAL I/O ERROR: Failed to load and orchestrate session from {file_path}: {e}")
            return 0

    def with_new_spectrum(self, new_spectrum: SpectrumV2) -> 'SessionV2':
        constructor_args = {
            'name': self.name,
            'gui': self._gui,
            'spec': new_spectrum,
            'systs': self.systs,
        }
        new_session = SessionV2(**constructor_args)
        new_session.log = self.log 
        new_session.defs = guarded_deepcopy_v1_state(self.defs)
        new_session.log_manager = self.log_manager 
        return new_session

    def with_new_system_list(self, new_systs: SystemListV2) -> 'SessionV2':
        constructor_args = {
            'name': self.name,
            'gui': self._gui,
            'spec': self.spec,
            'systs': new_systs,
        }
        new_session = SessionV2(**constructor_args)
        new_session.log = self.log 
        new_session.defs = guarded_deepcopy_v1_state(self.defs)
        new_session.log_manager = self.log_manager
        return new_session
    
    def save(self, file_path: str, 
             initial_session: Optional['SessionV2'] = None, 
             models: bool = False):
        """
        Saves the session.
        
        Args:
            file_path (str): The destination path.
            initial_session (Optional[SessionV2]): The *initial* state of the
                session, to be saved as "original data".
            models (bool): V1 compatibility flag.
        """
        if file_path.lower().endswith('.acs2'):
            # Pass *both* the initial and final (self) session states
            save_archive_v2(
                session_v2_final=self, 
                session_v2_initial=initial_session,
                file_path=file_path
            )
            logging.info(f"V2 Session state saved successfully to {file_path}.")
        else:
            # ... (V1 save logic is unchanged, it only saves final state) ...
            try:
                v1_spec_for_save = self.spec.to_v1_spectrum() 
                v1_systs_for_save = self.systs.to_v1_systlist() 
            except Exception as e:
                logging.error(f"FATAL: V2-to-V1 conversion failed during saving: {e}")
                raise 
            
            json_log_str = ""
            if isinstance(self.log_manager, GUILog):
                json_log_str = self.log_manager.str
            elif isinstance(self.log_manager, V1LogArtifact):
                json_log_str = json.dumps(self.log_manager.v1_json)
            elif isinstance(self.log_manager, HistoryLogV2):
                # We can serialize the V2 log for documentation in V1
                try:
                    json_log_str = json.dumps(dataclasses.asdict(self.log_manager))
                except Exception:
                    json_log_str = json.dumps({"set_menu": []})
            else:
                logging.warning("Cannot save V2 log to V1 .acs archive. Log will be empty.")
                json_log_str = json.dumps({"set_menu": []})
            
            try:
                save_archive_v1(
                    v1_spec_for_save, 
                    v1_systs_for_save, 
                    json_log_str,
                    file_path
                )
                logging.info(f"V1 compatibility session saved successfully to {file_path}.")
            except Exception as e:
                logging.error(f"FATAL: V1 archive writing failed: {e}")
                raise 
        return 0