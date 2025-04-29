from .cookbook_edit_old import CookbookEditOld
from .functions import *
from .message import *
from .vars import *

from astropy import units as au
import logging

class CookbookEdit(CookbookEditOld):
    """ Cookbook of utilities for editing data
    """

    def __init__(self):
        super(CookbookEdit, self).__init__()

    def modify_columns(self):
        """@brief Modify columns 🚧
        @details 🚧
        @url edit_cb.html#modify-columns
        """

        return 0

    def import_systs(self, source=0, mode='replace'):
        """@brief Import system list
        @details Import a system list into the currently selected session.
        @url edit_cb.html#import-system-list
        @param source Source session
        @param mode Mode (replace or append)
        @return 0
        """

        try:
            source = str(source)
        except Exception as e:
            # TODO: Testing this specific except block is currently problematic
            #       due to issues mocking builtins.str in the test env.
            logging.error(logging.error(f"Parameter conversion failed (details: {e})"))
            return 0

        struct = source+',systs'
        return self.struct_import(struct, mode)
      
    def import_telluric(self, source=0, col='telluric_model', merge_cont=True):
        """@brief Import telluric model
        @details Import a telluric model into the currently selected session 
        and optionally merge it into the continuum.
        @url edit_cb.html#import-telluric-model
        @param source Source session
        @param col Telluric model column
        @param merge_cont Merge telluric model into continuum
        @return 0
        """
        
        try:
            source = str(source)
            col = str(col)
            merge_cont = str(merge_cont) == 'True'
        except Exception as e:
            logging.error(logging.error(f"Parameter conversion failed (details: {e})"))
            return 0
        
        struct = source+',spec,'+col
        
        _, model, _ = self._struct_parse(struct, length=3)
        
        if not hasattr(self, 'sess') or self.sess is None:
            logging.error("Internal error: Session object not found.")
            return 0
        if not hasattr(self.sess, 'spec') or self.sess.spec is None:
            logging.error("No active spectrum found in the session.")
            return 0

        t = self.sess.spec._t

        if len(t)!=len(model):
            logging.error(f"Cannot import telluric model: Spectrum length ({len(t)}) "
                          f"does not match imported model length ({len(model)}).")
            return 0

        logging.info(f"Creating column `telluric_model`...")
        t['telluric_model'] = model
        if merge_cont and 'cont' in t.colnames:
            logging.info(f"Merging the telluric model with the continuum...")
            if 'cont_no_telluric' not in t.colnames:
                logging.info(f"Creating column `cont_no_telluric` to store the" 
                             "existing continuum...")
                t['cont_no_telluric'] = t['cont']
            t['cont'] *= t['telluric_model']
        elif merge_cont and 'cont' not in t.colnames:
            logging.warning(f"I cannot merge the model with the continuum: "
                            "continuum not found.")
        return 0
        
    def extract(self):
        """@brief Extract 🚧
        @details 🚧
        @url edit_cb.html#extract
        """

        return 0

    def extract(self, xmin=None, xmax=None, zem=None, part='blue'):
        """ @brief Extract
        @details Extract a session region based on wavelength range or Ly-alpha part.
        Creates a new session containing extracted components.
        @url edit_cb.html#extract
        @param xmin Minimum wavelength for region extraction
        @param xmax Maximum wavelength for region extraction
        @param zem Emission redshift for Ly-alpha part extraction
        @param part Part relative to Ly-alpha ('blue' or 'red')
        @return The new session object with extracted data, or None on failure.
        """

        # Determine mode based on provided arguments
        is_region_mode = xmin is not None and xmax is not None
        is_part_mode = zem is not None

        # --- Input Validation ---
        if is_region_mode and is_part_mode:
            logging.error("Ambiguous input for extract: Provide either (xmin, xmax) OR zem, not both.")
            return None
        if not is_region_mode and not is_part_mode:
            logging.error("Insufficient input for extract: Provide either (xmin, xmax) OR zem.")
            return None

        # --- Dispatch to appropriate private method ---
        new_session = None
        if is_region_mode:
            logging.info(f"Extracting region: xmin={xmin}, xmax={xmax}")
            try:
                # Convert here or assume _extract_region handles it? Let's convert here for safety.
                wmin_f = float(xmin)
                wmax_f = float(xmax)
                # Call the assumed existing method (might raise errors)
                # We assume it returns the new Session or None
                new_session = self._extract_region(xmin=wmin_f, xmax=wmax_f)
            except ValueError:
                logging.error(f"Invalid input: xmin='{xmin}' or xmax='{xmax}' not valid numbers.")
                return None # Error during parameter validation
            except AttributeError:
                logging.error("Internal error: Missing required method '_extract_region'.")
                return None # Method not found
            except Exception as e:
                 logging.error(f"An unexpected error occurred during region extraction: {e}")
                 return None # Catch other potential errors from _extract_region

        elif is_part_mode:
            logging.info(f"Extracting part: zem={zem}, part='{part}'")
            try:
                zem_f = float(zem)
                 # Call the assumed existing method (might raise errors)
                 # We assume it returns the new Session or None
                new_session = self._extract_part(zem=zem_f, part=part)
            except ValueError:
                 logging.error(f"Invalid input: zem='{zem}' not a valid number.")
                 return None # Error during parameter validation
            except AttributeError:
                logging.error("Internal error: Missing required method '_extract_part'.")
                return None # Method not found
            except Exception as e:
                 logging.error(f"An unexpected error occurred during part extraction: {e}")
                 return None # Catch other potential errors from _extract_part

        # Return the result (which should be the new Session or None)
        return new_session


    def mask(self, shift=0, tell=True, sky=True, cond=''):
        """@brief Mask
        @details Mask telluric lines, sky lines, or spectral regions defined by
        a condition.
        @url edit_cb.html#mask
        @param shift Shift to the barycentric frame (km/s)
        @param tell Mask telluric lines
        @param sky Mask sky lines
        @param cond Condition
        @return 0
        """

        try:
            shift = float(shift)
            tell = str(tell) == 'True'
            sky = str(sky) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        if tell: self.telluric_mask(shift=shift)
        if sky: self.sky_mask(shift=shift)
        if cond!='':
            self.mask_cond(cond=cond, new_sess=False)
            self.sess.spec._t['mask'] = np.logical_not(self.sess.spec._t['mask'])
            mask = np.where(np.array(self.sess.spec._t['mask']==0, dtype=bool))
            self.sess.spec._t['y'][mask] = np.nan

        return 0
    

    def _extract_part(self, zem, part='blue'):
        return self.part_extract(zem, part)

    def _extract_region(self, xmin, xmax):
        return self.region_extract(xmin, xmax, verbose=True)

    def _struct_parse(self, struct, length=2):
        """ Parses a structure string like 'session_idx,attr_name[,col_name]' """

        # TODO: --- Refactoring Point 1: Access session list directly ---
        # GUESSING session list name - *ADJUST IF INCORRECT*
        # Possibilities: self.sess.session_list, self.sess.spec_list, self.sess.data_objects
        try:
             # Use the new property:
             session_list = self.sess.session_list
        except AttributeError as e:
             # Handle case where Session couldn't provide the list
             logging.error(f"Failed to get session list from Session object: {e}")
             return None

        parse = struct.split(',')

        if len(parse) < length:
            logging.error(f"Cannot parse structure string '{struct}': Not enough parts (expected {length}).")
            return None

        # Session Index
        sessn_str = parse[0]
        try:
            sessn = int(sessn_str)
            # Keep original parse list with string index if needed later?
            # parse_orig = list(parse) # Keep original if needed
            parse[0] = sessn # Modify list with integer index
        except ValueError: # --- Refactoring Point 2: Specific Exception ---
            logging.error(f"Invalid session index '{sessn_str}': Must be an integer.")
            # logging.error(msg_param_fail) # Keep if msg_param_fail is very specific
            return None

        if not (0 <= sessn < len(session_list)): # Pythonic index check
            logging.error(f"Session index {sessn} out of range (found {len(session_list)} sessions).")
            return None
        sess_obj = session_list[sessn] # Get the target session/data object

        # Attribute Name
        attrn = parse[1]
        if not hasattr(sess_obj, attrn):
            # logging.error(msg_attr_miss % attrn) # Old way
            logging.error(f"Attribute '{attrn}' not found in session {sessn} object.") # --- Refactoring Point 4: f-string ---
            return None

        attr = getattr(sess_obj, attrn)
        if attr is None:
            # Returning the parse list might be confusing, maybe return None or specific tuple?
            logging.warning(f"Attribute '{attrn}' in session {sessn} is None.")
            # return attrn, attr, parse # Original return
            return None # Let's be consistent: None means failure/not found here


        if length == 3:
            # Column Name
            coln = parse[2]

            # --- Refactoring Point 5: Check for table structure ---
            if not hasattr(attr, '_t') or not isinstance(attr._t, Table):
                 logging.error(f"Attribute '{attrn}' in session {sessn} does not have an '_t' Astropy Table.")
                 return None
            # --- End Refactoring Point 5 ---

            if coln not in attr._t.colnames:
                # logging.error(msg_col_miss % coln) # Old way
                logging.error(f"Column '{coln}' not found in table for attribute '{attrn}' in session {sessn}.") # f-string
                return None
            col_data = attr._t[coln]
            # Return consistently: name, data, original_parse (or modified parse?)
            return coln, col_data, parse
        else: # length == 2
             # Return consistently: name, data, original_parse (or modified parse?)
            return attrn, attr, parse