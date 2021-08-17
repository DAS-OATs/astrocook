import json
import numpy as np
import os
import pathlib
import wx

class Defaults(object):

    def __init__(self,
                 gui):
        self._gui = gui
        self._indent = 2

        self._extend = {"voigt": {"z": 0.0, "logN": 13, "b": 10.0, "btur": 0.0,
                        "resol": 35000}}
        self.open()

    def open(self, file='defaults.json'):

        p = '/'.join(pathlib.PurePath(os.path.realpath(__file__)).parts[0:-1]) \
            + '/../' + file
        with open(p) as json_file:
            self.str = json_file.read()
            self.str = self.str.replace('“', '"')
            self.str = self.str.replace('”', '"')
            self.str = self.str.replace('—', '--')
        self.dict = json.loads(self.str)
        self._dict_extend()


    def update(self, str):
        self.str = str
        self.dict = json.loads(self.str)
        self._dict_extend()

    def _dict_extend(self):
        for k1 in self._extend:
            for k2 in self._extend[k1]:
                self.dict[k1][k2] = self._extend[k1][k2]
