from . import Frame

class Spec(Frame):
    """Class for spectra

    A spectrum is a Frame with methods for handling spectral operations."""

    def __init__(self, **kwargs):
        super(Spec, self).__init__(**kwargs)
