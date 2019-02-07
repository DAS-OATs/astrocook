from . import Spectrum

class Session(object):
    """ Class for sessions.

    A Session is a self-sufficient set of analysis operations."""

    def __init__(self):
        self.spec = Spectrum()
