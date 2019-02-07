from . import Frame

class Spectrum(Frame):
    """Class for spectra

    A Spectrum is a Frame with methods for handling spectral operations."""

    def __init__(self, **kwargs):
        super(Spectrum, self).__init__(**kwargs)
