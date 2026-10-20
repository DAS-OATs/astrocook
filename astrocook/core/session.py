import json
import logging
import os
import shutil
import tarfile 
import tempfile
from typing import Any, Optional, Union

from astrocook.io.adapter import load_and_migrate_structure, save_archive_v2
from astrocook.io.v1_stubs import V1ArchiveManager
from astrocook.recipes.continuum import RecipeContinuumV2
from astrocook.recipes.absorbers import RecipeAbsorbersV2
from astrocook.recipes.edit import RecipeEditV2
from astrocook.recipes.file import RecipeFileV2
from astrocook.recipes.flux import RecipeFluxV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import (
    SystemListDataV2, HistoryLogV2, V1LogArtifact
)
from astrocook.core.system_list import SystemListV2

# Define the type for the log manager
LogManager = Union[HistoryLogV2, V1LogArtifact]

def load_session_from_file(archive_path: str, name: str, gui_context: Any, format_name: str) -> Optional['SessionV2']:
    """
    Orchestrate the loading of a session from a file or archive.

    This function handles various input formats, including legacy Astrocook V1 archives (``.acs``),
    modern V2 archives (``.acs2``), and raw FITS files. It automatically handles
    V1-to-V2 migration of structures.

    Parameters
    ----------
    archive_path : str
        The full path to the file to load.
    name : str
        The display name to assign to the new session.
    gui_context : Any
        A reference to the main GUI window (for logging and dialog parenting),
        or a mock object if running in a script.
    format_name : str
        A format identifier string (e.g., ``'generic_spectrum'``, ``'eso_midas'``)
        used to select the appropriate loader for raw files. Use ``'auto'`` to attempt detection.

    Returns
    -------
    SessionV2 or int
        The loaded :class:`~astrocook.core.session.SessionV2` object if successful.
        Returns ``0`` (integer) if the loading process failed.

    Raises
    ------
    FileNotFoundError
        If the archive path is invalid or internal files are missing.
    RuntimeError
        If the V1 archive unpacking fails.
    """
    archive_manager = None 
    temp_dir = None
    archive_root = ""
    v2_metadata = None
    
    archive_path_lower = archive_path.lower()
    is_v2_archive = archive_path_lower.endswith('.acs2') or archive_path_lower.endswith('.tar.gz')

    # --- 1. V2 NATIVE LOADERS (Only for RAW files, NOT archives) ---
    from astrocook.io.loaders import get_loader, detect_file_format

    if not is_v2_archive and not archive_path_lower.endswith('.acs'):
        if format_name == 'auto' or format_name is None:
            try:
                format_name = detect_file_format(archive_path)
            except Exception:
                format_name = 'generic_spectrum'

        v2_loader = get_loader(format_name)
        if v2_loader:
            try:
                spec_data = v2_loader(archive_path)
                return SessionV2(name=name, gui=gui_context, spec=SpectrumV2(data=spec_data))
            except Exception as e:
                logging.error(f"V2 Loader failed: {e}")
                return 0

    # --- 2. ARCHIVE & LEGACY HANDLING ---
    try:
        if is_v2_archive:
            logging.debug(f"Unpacking V2 archive: {archive_path}")
            temp_dir = tempfile.mkdtemp()
            with tarfile.open(archive_path, 'r:gz') as tar:
                # "data" filter: refuse members that would land outside
                # temp_dir (absolute paths, "..", links), F35
                tar.extractall(path=temp_dir, filter="data")
            
            # 1. FIND THE ACTUAL FITS FILE
            spec_file_path = None
            for f in os.listdir(temp_dir):
                if f.endswith('_spec.fits'):
                    spec_file_path = os.path.join(temp_dir, f)
                    break
            
            if not spec_file_path:
                 raise FileNotFoundError("Could not find _spec.fits in .acs2 archive.")
            
            # 2. SET THE ROOT TO THE UNPACKED FILE (Crucial for the Adapter)
            archive_root = os.path.join(temp_dir, os.path.basename(spec_file_path).replace('_spec.fits', ''))
            
            # 3. LOAD METADATA
            meta_file_path = f"{archive_root}_meta.json"
            if os.path.exists(meta_file_path):
                with open(meta_file_path, 'r') as f:
                    v2_metadata = json.load(f)

        elif archive_path_lower.endswith('.acs'):
            # ... (V1 logic remains same) ...
            archive_manager = V1ArchiveManager(archive_path)
            temp_dir = archive_manager.unpack()
            spec_file_path = archive_manager.get_structure_path('spec')
            archive_root = os.path.splitext(spec_file_path)[0].replace('_spec', '')

            if format_name == 'auto' or format_name == 'acs':
                format_name = 'archive'

        else:
            # Single FITS file load
            temp_dir = None
            archive_root = os.path.splitext(archive_path)[0]
            spec_file_path = archive_path

        # --- 3. EXECUTE MIGRATION ---
        # archive_root MUST point to the directory/prefix where the .fits live
        spectrum_v2 = load_and_migrate_structure(
            archive_root, 'spec', gui_context, format_name, 
            spec_file_path=spec_file_path,
            v2_metadata=v2_metadata 
        )

        system_list_v2 = load_and_migrate_structure(
            archive_root, 'systs', gui_context, format_name,
            v2_metadata=v2_metadata
        )

        session = SessionV2(name=name, gui=gui_context, spec=spectrum_v2, systs=system_list_v2)
        if v2_metadata and v2_metadata.get('log_type') == 'v2':
            try:
                session.log_manager = HistoryLogV2.from_dict(
                    v2_metadata.get('log_history_json'))
            except (KeyError, TypeError, AttributeError) as e:
                logging.warning(f"Could not restore the session log: {e}")
        return session

    except Exception as e: 
        logging.error(f"Load failed: {e}", exc_info=True)
        return 0 
    finally:
        # Cleanup temp files
        if archive_manager: archive_manager.cleanup()
        if temp_dir and os.path.exists(temp_dir): shutil.rmtree(temp_dir)
            

class SessionV2:
    """
    The central coordinator for an Astrocook analysis session.

    This class encapsulates the state of a single spectrum and its associated
    models (systems, components). It serves as the primary entry point for
    high-level analysis workflows ("Recipes").

    .. note::
       **Immutability Pattern**: In Astrocook V2, the ``SessionV2`` object is designed
       to be treated as immutable within the undo/redo history stack. Operations that
       modify data (like smoothing flux or fitting a continuum) do not modify the object
       in-place; instead, they return a *new* ``SessionV2`` instance containing the
       updated data structures.

    Attributes
    ----------
    name : str
        The display name of the session (e.g., filename without extension).
    systs : astrocook.core.system_list.SystemListV2, optional
        The container for identified absorption systems and Voigt profile components.
    file : astrocook.recipes.file.RecipeFileV2
        Access to file import/export recipes (e.g., ``import_ascii``, ``export_ascii``).
    edit : astrocook.recipes.edit.RecipeEditV2
        Access to editing recipes (e.g., ``set_properties``, ``convert_x_axis``).
    flux : astrocook.recipes.flux.RecipeFluxV2
        Access to flux manipulation recipes (e.g., ``smooth``, ``rebin``).
    continuum : astrocook.recipes.continuum.RecipeContinuumV2
        Access to continuum estimation and fitting recipes.
    absorbers : astrocook.recipes.absorbers.RecipeAbsorbersV2
        Access to line identification and fitting recipes.
    """
    def __init__(self,
                 name: str,
                 gui: Any, 
                 spec: Optional[SpectrumV2] = None, 
                 systs: Optional[SystemListV2] = None):

        self._gui = gui
        self.name = name
        self.log_manager: Optional[LogManager] = None
        self._current_spectrum = spec
        self.systs = systs if systs is not None else SystemListV2(data=SystemListDataV2())
        
        self.file = RecipeFileV2(self)
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
        """
        The main container for spectral data (flux, wavelength, error, etc.).
        """
        return self._current_spectrum
    
    @property
    def lines(self):
        return self._lines 

    @classmethod
    def open_new(cls, file_path: str, name: str, gui_context: Any, format_name: str) -> 'SessionV2':        
        """
        Factory method: Create a new SessionV2 by loading data from a file.

        Parameters
        ----------
        file_path : str
            The system path to the file (FITS or ACS/ACS2 archive).
        name : str
            The name to give the new session.
        gui_context : Any
            The GUI parent object (used for logging/dialogs).
        format_name : str
            Hint for the file format loader (e.g. ``'generic_spectrum'``).

        Returns
        -------
        SessionV2 or int
            A new SessionV2 instance, or 0 on failure.
        """
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
        """
        Return a new SessionV2 instance with an updated spectrum.

        This is used by recipes that modify spectral data (e.g., smoothing,
        continuum fitting) to ensure history preservation.

        Parameters
        ----------
        new_spectrum : SpectrumV2
            The new spectral data object.

        Returns
        -------
        SessionV2
            A shallow copy of the session pointing to the new spectrum.
        """
        constructor_args = {
            'name': self.name,
            'gui': self._gui,
            'spec': new_spectrum,
            'systs': self.systs,
        }
        new_session = SessionV2(**constructor_args)
        new_session.log_manager = self.log_manager 
        return new_session

    def with_new_system_list(self, new_systs: SystemListV2) -> 'SessionV2':
        """
        Return a new SessionV2 instance with an updated system list.

        Parameters
        ----------
        new_systs : SystemListV2
            The new system list object.

        Returns
        -------
        SessionV2
            A shallow copy of the session pointing to the new system list.
        """
        constructor_args = {
            'name': self.name,
            'gui': self._gui,
            'spec': self.spec,
            'systs': new_systs,
        }
        new_session = SessionV2(**constructor_args)
        new_session.log_manager = self.log_manager
        return new_session
    
    def save(self, file_path: str, 
             initial_session: Optional['SessionV2'] = None, 
             models: bool = False):
        """
        Save the current session state to an ``.acs2`` archive.

        The archive contains the full session state, metadata, the session
        log and (if ``initial_session`` is provided) the original data.
        Saving in the V1 ``.acs`` format is no longer supported: V1
        archives can still be opened, and Astrocook 1.9.0 can still write
        them.

        Parameters
        ----------
        file_path : str
            The destination path. ``.acs2`` is appended if missing.
        initial_session : SessionV2, optional
            The state of the session *before* any modifications, saved in the
            archive so that it contains the original raw data for
            reproducibility.
        models : bool, optional
            Unused; kept for backward compatibility of the signature.

        Returns
        -------
        int
            Returns 0 on success.

        Raises
        ------
        ValueError
            If ``file_path`` has the V1 ``.acs`` extension.
        """
        if file_path.lower().endswith('.acs'):
            raise ValueError(
                "Saving in the V1 .acs format is no longer supported: save as "
                ".acs2, or use Astrocook 1.9.0 to write .acs archives.")
        if not file_path.lower().endswith('.acs2'):
            file_path += '.acs2'
        # Pass *both* the initial and final (self) session states
        save_archive_v2(
            session_v2_final=self, 
            session_v2_initial=initial_session,
            file_path=file_path
        )
        logging.info(f"V2 Session state saved successfully to {file_path}.")
        return 0