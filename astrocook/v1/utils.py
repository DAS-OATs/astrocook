from contextlib import contextmanager
import logging
import numpy as np
import sys


def parse_range(string, spec=None):
    if string=='all':
        if spec is None:
            xmin = 0
            xmax = np.inf
        else:
            xmin, xmax = np.min(spec._t['x']), np.max(spec._t['x'])
    elif string[0]=='[' and string[-1]==']':
        xmin, xmax = tuple([float(s) for s in string[1:-1].split(',')])
    else:
        logging.error("I cannot parse the range. The syntax is `[•,•]`.")
        return 0

    return xmin, xmax


@contextmanager
def recursion_limit(limit: int):
    """
    A context manager to temporarily set a new recursion limit.
    
    Example:
        with recursion_limit(2000):
            # Your deep recursive call here
            some_library.recursive_function(data)
    """
    original_limit = sys.getrecursionlimit()
    if limit <= original_limit:
        # If the new limit is not higher, do nothing.
        # This avoids a potential SystemError on some platforms.
        yield
        return
        
    try:
        sys.setrecursionlimit(limit)
        yield
    finally:
        # Always restore the original limit
        sys.setrecursionlimit(original_limit)


def run_with_recursion_retry(func, *args, new_limit=10000, **kwargs):
    """
    Executes a function, retrying once with an increased recursion limit
    if a RecursionError is caught.

    Args:
        func: The function to execute.
        *args: Positional arguments to pass to the function.
        new_limit (int): The recursion limit for the retry attempt.
        **kwargs: Keyword arguments to pass to the function.

    Returns:
        The result of the function call.

    Raises:
        RecursionError: If the function fails even with the increased limit.
        Any other exception raised by the function will propagate normally.
    """
    try:
        # First attempt with the default limit
        #logging.info(
        #    f"Attempting '{func.__name__}' with default limit ({sys.getrecursionlimit()})."
        #)
        return func(*args, **kwargs)
    except RecursionError:
        # First attempt failed, log it and prepare to retry
        #logging.warning(
        #    f"Caught RecursionError. Retrying '{func.__name__}' with limit={new_limit}."
        #)
        try:
            # Second attempt within the context manager
            with recursion_limit(new_limit):
                return func(*args, **kwargs)
        except RecursionError:
            # The second attempt also failed. This is a real problem.
            logging.error(
                f"Retry failed. '{func.__name__}' exceeded even the "
                f"increased recursion limit of {new_limit}."
            )
            # Re-raise the exception to signal a definitive failure
            raise