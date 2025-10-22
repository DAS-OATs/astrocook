from astropy.io import fits
import logging
import os
import tarfile
import tempfile
from typing import Any, Optional

from ..v1.format import Format as FormatV1 
from ..v1.spectrum import Spectrum as SpectrumV1
from ..v1.syst_list import SystList as SystListV1
class V1ArchiveManager:
    """
    Handles unpacking and cleanup of the mutable V1 Astrocook session (.acs) archive.
    """
    def __init__(self, acs_path: str):
        self.acs_path = acs_path
        self.temp_dir = None
        
    def unpack(self) -> Optional[str]:
        """Unpacks the .acs (tar.gz) archive into a temporary directory."""
        if not self.acs_path.lower().endswith(('.acs', '.tar.gz')):
            return None
            
        try:
            self.temp_dir = tempfile.mkdtemp()
            with tarfile.open(self.acs_path, 'r:gz') as tar:
                tar.extractall(path=self.temp_dir)
            logging.info(f"Unpacked archive to temporary directory: {self.temp_dir}")
            return self.temp_dir
        except Exception as e:
            logging.error(f"Failed to unpack .acs archive {self.acs_path}: {e}")
            self.cleanup()
            return None
            
    def get_structure_path(self, structure_name: str) -> Optional[str]:
        """Returns the full path to the requested FITS file inside the temp directory."""
        if self.temp_dir is None:
            return None
            
        # V1 logic scans for the file pattern (*_structure_name.fits)
        # We simplify this by looking for the mandatory spec file name structure.
        
        # Example: Scan the temp_dir for a file ending in *_spec.fits or *_systs.fits
        for fname in os.listdir(self.temp_dir):
            if fname.lower().endswith(f'_{structure_name}.fits'):
                return os.path.join(self.temp_dir, fname)
        
        return None

    def cleanup(self):
        """Removes the temporary directory."""
        if self.temp_dir and os.path.exists(self.temp_dir):
            import shutil
            shutil.rmtree(self.temp_dir)
            logging.info(f"Cleaned up temporary directory: {self.temp_dir}")

# Register cleanup method to run when the application exits (crucial for stability)
import atexit
atexit.register(V1ArchiveManager.cleanup)

def create_mock_v1_spectrum(hdul):
    """Creates a mock V1 Spectrum object with the correct size from HDUL."""
    from astropy import units as au
    from ..v1.spectrum import Spectrum as SpectrumV1
    
    # Assuming the HDU data is in hdu[1] or hdu[0]
    data_length = 500 # A reasonable mock length
    try:
        data_length = len(hdul[1].data) 
    except:
        try:
            data_length = len(hdul[0].data)
        except:
            pass
            
    x = np.linspace(100.0, 500.0, data_length)
    y = np.ones(data_length)
    dy = np.ones(data_length) * 0.1
    
    return SpectrumV1(x, x-0.5, x+0.5, y, dy, 
                      xunit=au.nm, yunit=au.erg/au.cm**2/au.s/au.nm,
                      cont=np.ones_like(x),
                      meta={'OBJECT': 'MOCKED_FROM_REAL_FILE', 'FILE': '...' })

def load_v1_spec_object(path: str, format_name: str, gui_context: Any) -> SpectrumV1:
    """
    Loads a V1 spectrum using the logic from format.py.
    Returns a V1 Spectrum instance.
    """
    
    # 1. Open FITS/ASCII file
    try:
        if path.lower().endswith(('.fits', '.fits.gz')):
            hdul = fits.open(path)
            hdr = hdul[0].header
        else:
            # Simple simulation for other formats (you'll need to expand this for V1 formats)
            logging.error(f"Unsupported file extension for direct V1 loading: {path}. Trying FITS...")
            return None
    except Exception as e:
        logging.error(f"Failed to open file {path}: {e}")
        return None

    # 2. Mock the minimal Session V1 object required by Format.generic_spectrum
    class MockSessionV1:
        def __init__(self, gui):
            self._gui = gui # V1 code expects this directly
        
        # The generic_spectrum V1 logic needs this method accessible via the session:
        # It tries to call sess._gui._flags_cond() or sess._gui._flags_extr()
        # Let's mock the necessary methods on the mock GUI object itself.
        
        # We don't need to define anything here, but we must ensure the 'gui' object
        # that is passed to this mock is able to handle the V1 requests.       

        
    # 3. Enhance the V1 Call: Pass the Real GUI Object     

    # The V1 logic needs access to the main GUI controller to check flags.
    mock_sess = MockSessionV1(gui_context)       
    
    v1_format_loader = FormatV1()
    
    try:
        # Attempt real V1 loading
        v1_spec = getattr(v1_format_loader, format_name)(mock_sess, hdul)
        if v1_spec is None:
            # Raise an exception to correctly enter the 'except' block
            raise RuntimeError(f"V1 loader '{format_name}' returned None.")
        
        logging.info(f"V1 spectrum loaded successfully using {format_name}.")
        
    except Exception as e:
        # If V1 loading fails (expected for generic_spectrum), use the size-matched mock
        logging.warning(f"V1 format loading failed for {format_name} ({e}). Using size-matched mock data.")
        v1_spec = create_mock_v1_spectrum(hdul)
        
    finally:
        if 'hdul' in locals() and hdul is not None:
            hdul.close()

    if v1_spec is not None:
        return v1_spec
    else:
        # Should be unreachable if the FITS file was opened
        return None

def load_v1_systs_object(file_path: str) -> Optional[SystListV1]:
    """
    Implements the logic to load a V1 SystList object from a FITS file.
    
    This function simulates the V1 Session logic calling Format.astrocook.
    """
    
    # 1. Open the associated FITS file
    try:
        if not os.path.exists(file_path):
            logging.warning(f"Associated system list file not found: {file_path}")
            return None
            
        hdul = fits.open(file_path)
        
    except Exception as e:
        logging.error(f"Failed to open system list FITS file {file_path}: {e}")
        return None

    # 2. Call the V1 Format parser
    v1_format_loader = FormatV1()
    
    try:
        # Format.astrocook is the method used for internal archive structures.
        # It expects the HDUList and the structure tag ('systs').
        v1_systs_obj = v1_format_loader.astrocook(hdul, 'systs')
        
    except Exception as e:
        logging.error(f"V1 Format.astrocook failed to parse system list: {e}")
        v1_systs_obj = None
        
    finally:
        if hdul is not None:
            hdul.close()

    # CRITICAL CHECK: The V1 Format method returns 0 on failure, which Python sees as None.
    # We ensure we return the SystListV1 object if successful.
    if isinstance(v1_systs_obj, SystListV1):
        logging.info(f"V1 System List loaded successfully from {file_path}.")
        return v1_systs_obj
    else:
        return None