import numpy as np
import pytest
from astrocook.cookbook_edit import CookbookEdit
from astrocook.spectrum import Spectrum
from astropy.table import Table
# Assuming logging and msg_param_fail are accessible
# If msg_param_fail is specific to the module, we might need to import it
# or mock logging more intricately. For now, let's focus on calls.
import logging

# --- Test Setup ---

# Re-using the helper from before, ensuring 'y' column exists
def create_test_spectrum(length=5):
    """Creates a minimal Spectrum object with y and mask columns."""
    wave = np.arange(length)
    flux = np.linspace(1.0, 5.0, length) # Using 'y' as flux here
    err = np.ones(length) * 0.1
    mask = np.zeros(length, dtype=bool) # Start with no mask
    t = Table({'wave': wave, 'y': flux, 'err': err, 'mask': mask})
    spec = Spectrum(t=t)
    # spec._meta = {'expr': 'Test Spectrum'} # Example if needed
    return spec

# Fake Session, slightly simpler as cookbook methods access sess.spec directly
class FakeSession:
    def __init__(self, active_spectrum):
        self.spec = active_spectrum # Direct access

# --- Tests for CookbookEdit.mask ---

@pytest.fixture
def cookbook_instance(mocker):
    """ Fixture to provide a CookbookEdit instance with mocked methods. """
    # Mock the methods that mask calls internally within CookbookEdit
    mocker.patch.object(CookbookEdit, 'telluric_mask')
    mocker.patch.object(CookbookEdit, 'sky_mask')
    mocker.patch.object(CookbookEdit, 'mask_cond')
    return CookbookEdit()

@pytest.fixture
def fake_session():
    """ Fixture to provide a basic FakeSession with a Spectrum. """
    return FakeSession(create_test_spectrum())

# --- Parameter Handling Tests ---

def test_mask_param_valid(cookbook_instance, fake_session):
    """ Test mask with default valid parameters, ensure no exceptions. """
    cookbook_instance.sess = fake_session # Attach session to cookbook
    # Default: tell=True, sky=True, cond=''
    result = cookbook_instance.mask()
    assert result == 0
    # Check default calls happened
    CookbookEdit.telluric_mask.assert_called_once_with(shift=0.0)
    CookbookEdit.sky_mask.assert_called_once_with(shift=0.0)
    CookbookEdit.mask_cond.assert_not_called()

def test_mask_param_conversion(cookbook_instance, fake_session, mocker):
    """ Test correct conversion of parameters. """
    cookbook_instance.sess = fake_session
    result = cookbook_instance.mask(shift='-10.5', tell='False', sky=False, cond='')
    assert result == 0
    CookbookEdit.telluric_mask.assert_not_called() # tell='False'
    CookbookEdit.sky_mask.assert_not_called()    # sky=False
    CookbookEdit.mask_cond.assert_not_called()

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
    CookbookEdit.telluric_mask.assert_not_called()
    CookbookEdit.sky_mask.assert_not_called()
    CookbookEdit.mask_cond.assert_not_called()

# --- Conditional Call Tests ---

def test_mask_telluric_only(cookbook_instance, fake_session):
    """ Test mask with only tell=True. """
    cookbook_instance.sess = fake_session
    result = cookbook_instance.mask(tell=True, sky=False, cond='')
    assert result == 0
    CookbookEdit.telluric_mask.assert_called_once_with(shift=0.0)
    CookbookEdit.sky_mask.assert_not_called()
    CookbookEdit.mask_cond.assert_not_called()

def test_mask_sky_only(cookbook_instance, fake_session):
    """ Test mask with only sky=True. """
    cookbook_instance.sess = fake_session
    result = cookbook_instance.mask(tell=False, sky=True, cond='')
    assert result == 0
    CookbookEdit.telluric_mask.assert_not_called()
    CookbookEdit.sky_mask.assert_called_once_with(shift=0.0)
    CookbookEdit.mask_cond.assert_not_called()

def test_mask_condition_only(cookbook_instance, fake_session):
    """ Test mask with only cond provided. """
    cookbook_instance.sess = fake_session
    condition_string = "wave > 5000"
    # We need mask_cond to *do* something to the mask for the rest to work
    def mock_mask_cond_action(cond, new_sess):
        # Simulate mask_cond setting indices 2 and 3 to True
        fake_session.spec._t['mask'][2] = True
        fake_session.spec._t['mask'][3] = True
    CookbookEdit.mask_cond.side_effect = mock_mask_cond_action

    result = cookbook_instance.mask(tell=False, sky=False, cond=condition_string)
    assert result == 0
    CookbookEdit.telluric_mask.assert_not_called()
    CookbookEdit.sky_mask.assert_not_called()
    CookbookEdit.mask_cond.assert_called_once_with(cond=condition_string, new_sess=False)

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
    CookbookEdit.telluric_mask.side_effect = mock_telluric_action

    # Simulate mask_cond masking indices 0, 1
    def mock_mask_cond_action(cond, new_sess):
        fake_session.spec._t['mask'][0] = True
        fake_session.spec._t['mask'][1] = True
    CookbookEdit.mask_cond.side_effect = mock_mask_cond_action

    result = cookbook_instance.mask(tell=True, sky=False, cond=condition_string, shift=10)
    assert result == 0
    CookbookEdit.telluric_mask.assert_called_once_with(shift=10.0)
    CookbookEdit.sky_mask.assert_not_called()
    CookbookEdit.mask_cond.assert_called_once_with(cond=condition_string, new_sess=False)

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