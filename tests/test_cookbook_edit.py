import numpy as np
import pytest
from astrocook.cookbook import Cookbook
from astrocook.spectrum import Spectrum
from astropy.table import Table
# Assuming logging and msg_param_fail are accessible
# If msg_param_fail is specific to the module, we might need to import it
# or mock logging more intricately. For now, let's focus on calls.
import logging

try:
    from astrocook.functions import create_xmin_xmax
except ImportError:
    logging.warning("Could not import create_xmin_xmax from astrocook.functions. Using fallback for tests.")
    # Basic fallback implementation (adjust if needed based on your real one)
    def create_xmin_xmax(x):
        x = np.asarray(x)
        if len(x) < 2:
             # Handle edge case with single point or empty array
             if len(x) == 1:
                  return np.array([x[0]-0.5]), np.array([x[0]+0.5])
             else:
                  return np.array([]), np.array([])
        x_c = (x[:-1] + x[1:]) * 0.5
        xmin = np.concatenate(([x[0] - (x_c[0] - x[0])], x_c))
        xmax = np.concatenate((x_c, [x[-1] + (x[-1] - x_c[-1])]))
        # Ensure same dtype as x if important (e.g., float)
        return xmin.astype(x.dtype, copy=False), xmax.astype(x.dtype, copy=False)


# --- Test Setup ---

def create_test_spectrum(length=5):
    """Creates a minimal Spectrum object with y and mask columns."""
    x = np.arange(length)
    xmin, xmax = create_xmin_xmax(x)
    y = np.linspace(1.0, 5.0, length) # Using 'y' as flux here
    dy = np.ones(length) * 0.1
    t = Table({'x': x, 'xmin': xmin, 'xmax': xmax, 'y': y, 'dy': dy})
    spec = Spectrum(
        x=t['x'],
        xmin=t['xmin'],
        xmax=t['xmax'],
        y=t['y'],
        dy=t['dy']
    )

    # Access the internal table (assuming it's spec._t)
    # This relies on internal knowledge of Spectrum class structure
    if hasattr(spec, '_t') and isinstance(spec._t, Table):
        if 'mask' not in spec._t.colnames:
            # Add a default mask column (all False) if it doesn't exist
            mask_column = np.zeros(len(spec.x), dtype=bool)
            spec._t['mask'] = mask_column
            # logging.debug("Added default 'mask' column to spec._t in test setup")
        # else: mask column already exists, do nothing
    else:
        # This would indicate our assumption about spec._t is wrong
        # or Spectrum init failed unexpectedly in the test context.
        raise AttributeError("Spectrum object in test setup doesn't have "
                           "an Astropy Table attribute named '_t'. "
                           "Cannot add 'mask' column.")

    # spec._meta = {'expr': 'Test Spectrum'} # Example if needed
    return spec


@pytest.fixture
def cookbook_instance(mocker):
    """ Fixture to provide a Cookbook instance with mocked methods. """
    # Mock the methods that mask calls internally within Cookbook
    mocker.patch.object(Cookbook, 'telluric_mask')
    mocker.patch.object(Cookbook, 'sky_mask')
    mocker.patch.object(Cookbook, 'mask_cond')
    return Cookbook()


# --- Fake Session Object Wrapper ---
class FakeSessionObject:
    """ A dummy object to hold data within the fake session list.
        It mimics having a '.spec' attribute.
    """
    def __init__(self, spectrum_object):
        self.spec = spectrum_object # REQUIRED: Has a .spec attribute
        # Add other potential attributes if needed by other parsing structures
        # self.systs = None
        # self.metadata = {}

# --- Fake GUI & Fake Session (Revised for Wrapper) ---
class FakeGui:
    def __init__(self):
        self._sess_list = [] # Now holds FakeSessionObject instances

class FakeSession:
    def __init__(self):
        """ Initializes an empty fake session environment. """
        self._gui = FakeGui()
        # This 'active' spec should probably point to the .spec of the active FakeSessionObject
        self.spec = None # Represents the active Spectrum object

    def add_session_object(self, session_obj_wrapper):
        """ Adds a mock session object (wrapper) to the list. """
        self._gui._sess_list.append(session_obj_wrapper)
        # Update the 'active' spec pointer if this is the first one
        if self.spec is None and self._gui._sess_list:
             # Point to the .spec attribute of the first session object wrapper
             self.spec = self._gui._sess_list[0].spec

    # Helper to add a Spectrum directly by wrapping it
    def add_spectrum(self, spectrum_object):
        """ Convenience method to wrap a Spectrum in a FakeSessionObject and add it. """
        wrapper = FakeSessionObject(spectrum_object=spectrum_object)
        self.add_session_object(wrapper)

    @property
    def session_list(self):
        """ Provides access to the list in the fake GUI object. """
        # Mirrors the structure of the real Session property
        if hasattr(self, '_gui') and self._gui is not None and hasattr(self._gui, '_sess_list'):
             return self._gui._sess_list
        else:
             # Match the error/fallback behavior of the real property
             raise AttributeError("FakeSession object cannot provide session_list")
             # Or return [] if the real one does


# --- Fixture for mask/systs/general tests (Revised) ---
@pytest.fixture
def fake_session(mocker):
    """ Fixture providing a FakeSession with one default Spectrum wrapped in a FakeSessionObject. """
    # mocker.patch('logging.error') # Mock logging if needed by mask/systs
    session = FakeSession()
    # Use the helper to wrap and add the default spectrum
    session.add_spectrum(create_test_spectrum())
    return session

# --- Fixture for telluric tests (Revised) ---
@pytest.fixture
def session_for_telluric(mocker):
    """ Creates a FakeSession with a spectrum (in a wrapper) prepared for telluric tests. """
    mocker.patch('logging.error')
    mocker.patch('logging.info')
    mocker.patch('logging.warning')

    # Create the specific spectrum needed
    length = 10
    x = np.arange(length, dtype=float)
    xmin, xmax = create_xmin_xmax(x)
    y = np.ones(length) * 5.0
    dy = np.ones(length) * 0.1
    cont = np.ones(length) * 4.0
    spec = Spectrum(x=x, xmin=xmin, xmax=xmax, y=y, dy=dy)
    if hasattr(spec, '_t') and isinstance(spec._t, Table):
         spec._t['cont'] = cont
         if 'mask' not in spec._t.colnames:
              spec._t['mask'] = np.zeros(len(spec.x), dtype=bool)
    else:
         raise AttributeError("Spectrum setup failed in fixture")

    # Create empty session and add the wrapped spectrum
    session = FakeSession()
    session.add_spectrum(spec) # Use helper
    return session

# --- Fixture for parsing tests (Revised) ---
@pytest.fixture
def session_for_parsing(mocker):
    """ Creates an enhanced FakeSession for parsing tests, mocks logging """
    mocker.patch('logging.error')
    mocker.patch('logging.warning')
    mocker.patch('logging.info')
    session = FakeSession() # Create empty session
    return session # Tests will add FakeSessionObjects as needed


# --- Tests for _struct_parse ---

def test_struct_parse_success_len2(cookbook_instance, session_for_parsing):
    """ Test successful parsing for length 2 (attribute '.spec'). """
    mock_spec = create_test_spectrum(length=5)
    # Use the add_spectrum helper which creates the wrapper
    session_for_parsing.add_spectrum(mock_spec)
    cookbook_instance.sess = session_for_parsing

    struct = "0,spec" # Ask for the .spec attribute of the object at index 0
    result = cookbook_instance._struct_parse(struct, length=2)

    assert result is not None
    attr_name, attr_obj, parse_list = result
    assert attr_name == 'spec'
    # Check it returned the *actual Spectrum object* held within the wrapper
    assert attr_obj is mock_spec
    assert parse_list == [0, 'spec']

def test_struct_parse_success_len3(cookbook_instance, session_for_parsing):
    """ Test successful parsing for length 3 (column within '.spec'). """
    mock_spec = create_test_spectrum(length=5)
    session_for_parsing.add_spectrum(mock_spec)
    cookbook_instance.sess = session_for_parsing

    struct = "0,spec,y" # Ask for column 'y' from the .spec attribute of object 0
    result = cookbook_instance._struct_parse(struct, length=3)

    assert result is not None
    col_name, col_data, parse_list = result
    assert col_name == 'y'
    np.testing.assert_allclose(col_data, np.linspace(1.0, 5.0, 5))
    assert parse_list == [0, 'spec', 'y']

def test_struct_parse_fail_short_string(cookbook_instance, session_for_parsing):
    """ Test failure when struct string is too short. """
    cookbook_instance.sess = session_for_parsing
    result = cookbook_instance._struct_parse("0", length=2)
    assert result is None
    logging.error.assert_called_once()

def test_struct_parse_fail_non_integer_session(cookbook_instance, session_for_parsing):
    """ Test failure with non-integer session index. """
    cookbook_instance.sess = session_for_parsing
    result = cookbook_instance._struct_parse("abc,spec", length=2)
    assert result is None
    logging.error.assert_called_once() # Checks if the ValueError except block logs

def test_struct_parse_fail_session_out_of_range(cookbook_instance, session_for_parsing):
    """ Test failure when session index is out of range. """
    # Setup: Add only one object (index 0)
    session_for_parsing.add_spectrum(create_test_spectrum())
    cookbook_instance.sess = session_for_parsing

    result = cookbook_instance._struct_parse("1,spec", length=2) # Ask for index 1
    assert result is None
    logging.error.assert_called_once()

def test_struct_parse_fail_attribute_not_found(cookbook_instance, session_for_parsing):
    """ Test failure when the attribute doesn't exist on the session object wrapper. """
    mock_spec = create_test_spectrum()
    # Add the spectrum normally (it gets wrapped with a .spec attribute)
    session_for_parsing.add_spectrum(mock_spec)
    cookbook_instance.sess = session_for_parsing

    # Ask for an attribute that doesn't exist on the FakeSessionObject wrapper
    result = cookbook_instance._struct_parse("0,nonexistent_attr", length=2)
    assert result is None
    logging.error.assert_called_once()


def test_struct_parse_fail_column_not_found(cookbook_instance, session_for_parsing):
    """ Test failure for length 3 when column doesn't exist in the .spec attribute. """
    mock_spec = create_test_spectrum(length=5)
    session_for_parsing.add_spectrum(mock_spec)
    cookbook_instance.sess = session_for_parsing

    # Ask for column in .spec that doesn't exist
    result = cookbook_instance._struct_parse("0,spec,nonexistent_col", length=3)
    assert result is None
    logging.error.assert_called_once()

# --- Mock Object for testing missing _t ---
class MockSessionObjectWithSpecNoT:
     """ Represents an object in the session list.
         It has a .spec attribute, but that attribute lacks '_t'.
     """
     def __init__(self):
         self.spec = "This is the spec attribute, but it's just a string"


# --- Mock Object for testing wrong type of _t ---
class MockSessionObjectWithSpecWrongT:
     """ Represents an object in the session list.
         It has a .spec attribute, and spec._t exists, but it's not a Table.
     """
     def __init__(self):
         # Create a dummy object with a _t that is a list
         dummy_spec_attr = type('DummySpecAttr', (object,), {'_t': [10, 20, 30]})()
         self.spec = dummy_spec_attr

def test_struct_parse_fail_attr_no_table(cookbook_instance, session_for_parsing):
     """ Test failure for len 3 when attribute doesn't have _t table. """
     # Create an object that doesn't have _t for its relevant attribute
     mock_obj_wrapper = MockSessionObjectWithSpecNoT() # This IS the session object
     session_for_parsing.add_session_object(mock_obj_wrapper) # Add it directly
     cookbook_instance.sess = session_for_parsing

     # Ask for a column within the .some_attr attribute (which doesn't have _t)
     result = cookbook_instance._struct_parse("0,spec,any_col", length=3)
     assert result is None
     logging.error.assert_called_once()
     logging.warning.assert_not_called()

def test_struct_parse_fail_attr_t_not_table(cookbook_instance, session_for_parsing):
     """ Test failure for len 3 when attribute._t is not an Astropy Table. """
     mock_obj_wrapper = MockSessionObjectWithSpecWrongT() # This wrapper has .spec, but spec._t is wrong
     session_for_parsing.add_session_object(mock_obj_wrapper) # Add it directly
     cookbook_instance.sess = session_for_parsing

     # Ask for a column within the .spec attribute (whose _t is not a Table)
     result = cookbook_instance._struct_parse("0,spec,any_col", length=3)
     assert result is None
     logging.error.assert_called_once()
     logging.warning.assert_not_called()


# --- Tests for import_systs ---

def test_import_systs_calls_struct_import(cookbook_instance, fake_session, mocker):
    """ Test that import_systs calls struct_import with correct args. """
    cookbook_instance.sess = fake_session # Attach session
    # Mock the method that import_systs delegates to
    mock_struct_import = mocker.patch.object(cookbook_instance, 'struct_import')
    mock_struct_import.return_value = 0
    source_session = 'other_session_name'
    mode = 'append'
    result = cookbook_instance.import_systs(source=source_session, mode=mode)

    assert result == 0 # Assuming struct_import returns 0 on success
    expected_struct = f"{source_session},systs"
    mock_struct_import.assert_called_once_with(expected_struct, mode)

def test_import_systs_param_conversion(cookbook_instance, fake_session, mocker):
    """ Test parameter conversion for source. """
    cookbook_instance.sess = fake_session
    mock_struct_import = mocker.patch.object(cookbook_instance, 'struct_import')
    mock_struct_import.return_value = 0
    result = cookbook_instance.import_systs(source=123, mode='replace') # Use int source

    assert result == 0
    expected_struct = "123,systs" # Source should be converted to string
    mock_struct_import.assert_called_once_with(expected_struct, 'replace')

#def test_import_systs_invalid_param(cookbook_instance, fake_session, mocker):
#    """ Test handling of invalid source parameter during conversion. """
#    cookbook_instance.sess = fake_session
#    # Mock struct_import first (so we can check it's not called)
#    mock_struct_import = mocker.patch.object(cookbook_instance, 'struct_import')
#    # Mock logging error
#    mock_log_error = mocker.patch('logging.error')
#    # Mock str to raise an error *only for this test*
#    str_patch = mocker.patch('builtins.str', side_effect=TypeError("Fake conversion error"))
#
#    # Use a value that *would* normally convert fine, to ensure the mock is hit
#    source_value = 123
#
#    # --- Run the code ---
#    result = cookbook_instance.import_systs(source=source_value)
#
#    # --- Assertions ---
#    assert result == 0 # Should return 0 on error
#    mock_log_error.assert_called_once() # Check error was logged
#    mock_struct_import.assert_not_called() # Check delegate method not called
#
#    # --- Optional: Explicitly check the mock call count for str ---
#    # This ensures our mock was actually triggered by the code under test
#    str_patch.assert_called_once_with(source_value)
#
#    # pytest-mock should automatically stop the patch here


# --- Tests for import_telluric ---

def test_import_telluric_success_no_merge(cookbook_instance, session_for_telluric, mocker):
    """ Test successful import without merging continuum. """
    cookbook_instance.sess = session_for_telluric # Attach session
    spec_table = session_for_telluric.spec._t
    source_session = 'source_sess'
    source_col = 'sky_model'
    expected_length = len(spec_table)
    # Mock the return value of _struct_parse: _, model, _
    mock_model_data = np.linspace(0.9, 1.0, expected_length) # Example model
    mock_struct_parse = mocker.patch.object(cookbook_instance, '_struct_parse',
                                           return_value=(None, mock_model_data, None))

    result = cookbook_instance.import_telluric(source=source_session, col=source_col, merge_cont=False)

    assert result == 0
    expected_struct = f"{source_session},spec,{source_col}"
    mock_struct_parse.assert_called_once_with(expected_struct, length=3)
    logging.error.assert_not_called()
    logging.warning.assert_not_called()
    # Check 'telluric_model' column was added and has correct data
    assert 'telluric_model' in spec_table.colnames
    np.testing.assert_allclose(spec_table['telluric_model'], mock_model_data)
    # Check 'cont' column was NOT modified
    assert 'cont' in spec_table.colnames
    np.testing.assert_allclose(spec_table['cont'], np.ones(expected_length) * 4.0)
    assert 'cont_no_telluric' not in spec_table.colnames # Backup shouldn't exist

def test_import_telluric_success_with_merge(cookbook_instance, session_for_telluric, mocker):
    """ Test successful import with merging continuum. """
    cookbook_instance.sess = session_for_telluric
    spec_table = session_for_telluric.spec._t
    original_cont = spec_table['cont'].copy() # Save original continuum
    source_session = 'source_sess'
    source_col = 'sky_model'
    expected_length = len(spec_table)
    # Mock the return value of _struct_parse
    mock_model_data = np.linspace(0.9, 1.0, expected_length) # Example model (0.9 to 1.0)
    mock_struct_parse = mocker.patch.object(cookbook_instance, '_struct_parse',
                                           return_value=(None, mock_model_data, None))

    result = cookbook_instance.import_telluric(source=source_session, col=source_col, merge_cont=True)

    assert result == 0
    expected_struct = f"{source_session},spec,{source_col}"
    mock_struct_parse.assert_called_once_with(expected_struct, length=3)
    logging.error.assert_not_called()
    logging.warning.assert_not_called()
    # Check 'telluric_model' column exists
    assert 'telluric_model' in spec_table.colnames
    np.testing.assert_allclose(spec_table['telluric_model'], mock_model_data)
    # Check backup continuum exists and is correct
    assert 'cont_no_telluric' in spec_table.colnames
    np.testing.assert_allclose(spec_table['cont_no_telluric'], original_cont)
    # Check 'cont' column WAS modified (original * model)
    expected_merged_cont = original_cont * mock_model_data
    np.testing.assert_allclose(spec_table['cont'], expected_merged_cont)

def test_import_telluric_merge_cont_already_merged(cookbook_instance, session_for_telluric, mocker):
    """ Test merging when cont_no_telluric already exists. """
    cookbook_instance.sess = session_for_telluric
    spec_table = session_for_telluric.spec._t
    # Add a dummy backup column beforehand
    spec_table['cont_no_telluric'] = spec_table['cont'] * 0.5 # Make it different
    original_cont_no_telluric = spec_table['cont_no_telluric'].copy()
    original_cont = spec_table['cont'].copy()

    expected_length = len(spec_table)
    mock_model_data = np.linspace(0.9, 1.0, expected_length)
    mocker.patch.object(cookbook_instance, '_struct_parse', return_value=(None, mock_model_data, None))

    result = cookbook_instance.import_telluric(source='s', col='c', merge_cont=True)

    assert result == 0
    logging.error.assert_not_called()
    # Check the backup column was NOT overwritten
    np.testing.assert_allclose(spec_table['cont_no_telluric'], original_cont_no_telluric)
    # Check cont was still updated
    expected_merged_cont = original_cont * mock_model_data
    np.testing.assert_allclose(spec_table['cont'], expected_merged_cont)

def test_import_telluric_merge_no_cont_column(cookbook_instance, session_for_telluric, mocker):
    """ Test merging when the original spectrum has no 'cont' column. """
    cookbook_instance.sess = session_for_telluric
    spec_table = session_for_telluric.spec._t
    # Remove the 'cont' column for this test
    del spec_table['cont']

    expected_length = len(spec_table)
    mock_model_data = np.linspace(0.9, 1.0, expected_length)
    mocker.patch.object(cookbook_instance, '_struct_parse', return_value=(None, mock_model_data, None))

    result = cookbook_instance.import_telluric(source='s', col='c', merge_cont=True)

    assert result == 0
    logging.error.assert_not_called()
    # Check warning was issued
    logging.warning.assert_called_once()
    # Check 'telluric_model' still added
    assert 'telluric_model' in spec_table.colnames
    np.testing.assert_allclose(spec_table['telluric_model'], mock_model_data)
    # Check no continuum columns exist
    assert 'cont' not in spec_table.colnames
    assert 'cont_no_telluric' not in spec_table.colnames

def test_import_telluric_length_mismatch(cookbook_instance, session_for_telluric, mocker):
    """ Test error handling when imported model length doesn't match spectrum. """
    cookbook_instance.sess = session_for_telluric
    spec_table = session_for_telluric.spec._t
    wrong_length = len(spec_table) + 5
    # Mock _struct_parse to return model with wrong length
    mock_model_data = np.linspace(0.9, 1.0, wrong_length)
    mock_struct_parse = mocker.patch.object(cookbook_instance, '_struct_parse',
                                           return_value=(None, mock_model_data, None))

    result = cookbook_instance.import_telluric(source='s', col='c', merge_cont=False)

    assert result == 0 # Returns 0 even on handled error
    mock_struct_parse.assert_called_once()
    # Check error was logged
    logging.error.assert_called_once()
    # Check 'telluric_model' was NOT added
    assert 'telluric_model' not in spec_table.colnames


# --- Tests for mask ---

def test_mask_param_valid(cookbook_instance, fake_session):
    """ Test mask with default valid parameters, ensure no exceptions. """
    cookbook_instance.sess = fake_session # Attach session to cookbook
    # Default: tell=True, sky=True, cond=''
    result = cookbook_instance.mask()
    assert result == 0
    # Check default calls happened
    Cookbook.telluric_mask.assert_called_once_with(shift=0.0)
    Cookbook.sky_mask.assert_called_once_with(shift=0.0)
    Cookbook.mask_cond.assert_not_called()

def test_mask_param_conversion(cookbook_instance, fake_session, mocker):
    """ Test correct conversion of parameters. """
    cookbook_instance.sess = fake_session
    result = cookbook_instance.mask(shift='-10.5', tell='False', sky=False, cond='')
    assert result == 0
    Cookbook.telluric_mask.assert_not_called() # tell='False'
    Cookbook.sky_mask.assert_not_called()    # sky=False
    Cookbook.mask_cond.assert_not_called()

def test_mask_param_invalid_shift(cookbook_instance, fake_session, mocker):
    """ Test mask with invalid shift parameter. """
    # Mock logging.error to check if it's called
    mock_log_error = mocker.patch('logging.error')
    cookbook_instance.sess = fake_session

    result = cookbook_instance.mask(shift='not_a_number')
    assert result == 0 # Should return 0 on error
    # Check that logging.error was called (ideally check the message too)
    mock_log_error.assert_called_once()
    # Check that masking methods were NOT called because of early exit
    Cookbook.telluric_mask.assert_not_called()
    Cookbook.sky_mask.assert_not_called()
    Cookbook.mask_cond.assert_not_called()

def test_mask_telluric_only(cookbook_instance, fake_session):
    """ Test mask with only tell=True. """
    cookbook_instance.sess = fake_session
    result = cookbook_instance.mask(tell=True, sky=False, cond='')
    assert result == 0
    Cookbook.telluric_mask.assert_called_once_with(shift=0.0)
    Cookbook.sky_mask.assert_not_called()
    Cookbook.mask_cond.assert_not_called()

def test_mask_sky_only(cookbook_instance, fake_session):
    """ Test mask with only sky=True. """
    cookbook_instance.sess = fake_session
    result = cookbook_instance.mask(tell=False, sky=True, cond='')
    assert result == 0
    Cookbook.telluric_mask.assert_not_called()
    Cookbook.sky_mask.assert_called_once_with(shift=0.0)
    Cookbook.mask_cond.assert_not_called()

def test_mask_condition_only(cookbook_instance, fake_session):
    """ Test mask with only cond provided. """
    cookbook_instance.sess = fake_session
    condition_string = "wave > 5000"
    # We need mask_cond to *do* something to the mask for the rest to work
    def mock_mask_cond_action(cond, new_sess):
        # Simulate mask_cond setting indices 2 and 3 to True
        fake_session.spec._t['mask'][2] = True
        fake_session.spec._t['mask'][3] = True
    Cookbook.mask_cond.side_effect = mock_mask_cond_action

    result = cookbook_instance.mask(tell=False, sky=False, cond=condition_string)
    assert result == 0
    Cookbook.telluric_mask.assert_not_called()
    Cookbook.sky_mask.assert_not_called()
    Cookbook.mask_cond.assert_called_once_with(cond=condition_string, new_sess=False)

    # Now check the critical logic: inversion and NaN setting
    # Initial mask: [F, F, F, F, F]
    # After mock_mask_cond: [F, F, T, T, F]
    # After logical_not: [T, T, F, F, T] -> This is the final mask state
    expected_final_mask = np.array([True, True, False, False, True])
    np.testing.assert_array_equal(fake_session.spec._t['mask'], expected_final_mask)

    # Find indices where final mask is False (indices 2, 3)
    # Check that y at these indices is NaN
    expected_y = np.array([1.0, 2.0, np.nan, np.nan, 5.0]) # Original y: [1, 2, 3, 4, 5]
    np.testing.assert_allclose(fake_session.spec._t['y'], expected_y, equal_nan=True)


def test_mask_telluric_and_condition(cookbook_instance, fake_session):
    """ Test mask with tell=True and cond provided. """
    cookbook_instance.sess = fake_session
    condition_string = "y < 2.5" # Should affect indices 0, 1

    # Simulate telluric_mask masking index 4
    def mock_telluric_action(shift):
        fake_session.spec._t['mask'][4] = True
    Cookbook.telluric_mask.side_effect = mock_telluric_action

    # Simulate mask_cond masking indices 0, 1
    def mock_mask_cond_action(cond, new_sess):
        fake_session.spec._t['mask'][0] = True
        fake_session.spec._t['mask'][1] = True
    Cookbook.mask_cond.side_effect = mock_mask_cond_action

    result = cookbook_instance.mask(tell=True, sky=False, cond=condition_string, shift=10)
    assert result == 0
    Cookbook.telluric_mask.assert_called_once_with(shift=10.0)
    Cookbook.sky_mask.assert_not_called()
    Cookbook.mask_cond.assert_called_once_with(cond=condition_string, new_sess=False)

    # Trace the mask state:
    # Initial mask:      [F, F, F, F, F]
    # After telluric:    [F, F, F, F, T]
    # After mask_cond:   [T, T, F, F, T] (mask is ORed)
    # After logical_not: [F, F, T, T, F] -> This is the final mask state
    expected_final_mask = np.array([False, False, True, True, False])
    np.testing.assert_array_equal(fake_session.spec._t['mask'], expected_final_mask)

    # Find indices where final mask is False (indices 0, 1, 4)
    # Check that y at these indices is NaN
    expected_y = np.array([np.nan, np.nan, 3.0, 4.0, np.nan]) # Original y: [1, 2, 3, 4, 5]
    np.testing.assert_allclose(fake_session.spec._t['y'], expected_y, equal_nan=True)