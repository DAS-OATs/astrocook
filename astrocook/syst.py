import copy
import numpy as np

class Syst():

    def __init__(self, lines):

        self._lines = copy.deepcopy(lines)
