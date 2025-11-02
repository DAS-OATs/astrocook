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
from .io_adapter import load_and_migrate_structure, save_archive_v2
from .io_v1_stubs import V1ArchiveManager, save_archive_v1
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

def load_session_from_file(archive_path: str, name: str, gui_context: Any, format_name: str) -> Tuple[Optional['SessionV2'], Optional[LogManager]]:
    """
    Orchestrates the loading of a session from an archive (.acs or .acs2)
    and handles the V1-to-V2 migration if necessary.
    """
    archive_manager = None 
    temp_dir = None
    archive_root = ""
    v2_metadata = None
    log_object_to_return: Optional[LogManager] = None # <<< *** BUG FIX IS HERE ***
    
    try:
        archive_path_lower = archive_path.lower()
        
        if archive_path_lower.endswith('.acs2') or archive_path_lower.endswith('.tar.gz'):
            # This is a V2 archive (tar.gz)
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
            # This is a V1 archive
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
            # This is a single FITS file
            logging.debug("Loading single FITS file.")
            temp_dir = None
            archive_root = os.path.splitext(archive_path)[0]
            spec_file_path = archive_path
        
        # --- V2 METADATA LOADING ---
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

        # 1. Load Spectrum (This will now use the V2 path if v2_metadata is not None)
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

        # 3. Load V2/V1 Log Object
        try:
            if v2_metadata:
                log_data = v2_metadata.get('log_history_json')
                log_type = v2_metadata.get('log_type', 'v1_legacy')
                if log_type == 'v2' and isinstance(log_data, dict):
                    log_obj = HistoryLogV2()
                    log_obj.current_index = log_data.get('current_index', -1)
                    log_obj.entries = [LogEntryV2(**entry) for entry in log_data.get('entries', [])]
                    log_object_to_return = log_obj
                    logging.info(f"Restored V2 log with {len(log_obj.entries)} entries.")
                elif log_type in ('v1_artifact', 'v1_legacy'):
                    if isinstance(log_data, str): log_data = json.loads(log_data)
                    log_object_to_return = V1LogArtifact(log_data or {"set_menu": []})
                    logging.info("Restored V1 log artifact.")
            elif temp_dir and archive_manager: # <<< Make sure this is V1 path
                log_path = archive_manager.get_structure_path('log.json')
                if not log_path:
                    log_path_alt = os.path.join(temp_dir, f"{os.path.basename(archive_root)}_log.json")
                    if os.path.exists(log_path_alt): log_path = log_path_alt
                if log_path and os.path.exists(log_path):
                    with open(log_path, 'r') as f:
                        log_string = f.read()
                    log_object_to_return = V1LogArtifact(json.loads(log_string))
                    logging.info(f"V1 log artifact restored from {log_path}.")
            
            if log_object_to_return is None:
                 logging.debug("No log found, creating new HistoryLogV2.")
                 log_object_to_return = HistoryLogV2()
        except Exception as e:
            logging.error(f"Failed to restore log history: {e}. Creating new V2 log.")
            log_object_to_return = HistoryLogV2()
            
        
        # Call the simplified constructor
        new_session = SessionV2(
            name=name, 
            gui=gui_context, 
            spec=spectrum_v2, 
            systs=system_list_v2
        )
        
        return new_session, log_object_to_return 

    except Exception as e: 
        logging.error(f"FATAL: load_session_from_file failed: {e}", exc_info=True)
        return 0, None

    finally: 
        # Cleanup
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
            session, log_object = load_session_from_file(
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
    
    def save(self, file_path: str, models: bool = False):
        if file_path.lower().endswith('.acs2'):
            save_archive_v2(self, file_path)
            logging.info(f"V2 Session state saved successfully to {file_path}.")
        else:
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