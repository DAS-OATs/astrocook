import json
import numpy as np
import os
import pathlib
import wx


def _dict_merge(dict1, dict2):
    for k1 in dict1:
        if k1 not in dict2:
            dict2[k1] = {}
        for k2 in dict1[k1]:
            dict2[k1][k2] = dict1[k1][k2]
    return dict2


class Defaults(object):

    def __init__(self,
                 gui):
        self._gui = gui
        self._indent = 2

        self._extend = {"voigt": {"z": 0.0, "logN": 13, "b": 10.0, "btur": 0.0,
                        "resol": 35000}}
        self.open()

    def open(self, path=None, file='defaults.json'):

        if path is None:
            pathfile = '/'.join(pathlib.PurePath(os.path.realpath(__file__)).parts[0:-1]) \
                       + '/../' + file
        else:
            pathfile = path+'/'+file
        with open(pathfile) as json_file:
            self.str = json_file.read()
            self.str = self.str.replace('“', '"')
            self.str = self.str.replace('”', '"')
            self.str = self.str.replace('—', '--')
        self.dict = json.loads(self.str)
        #self._dict_extend()
        self.dict = _dict_merge(self._extend, self.dict)


    def update(self, str):
        self.str = str
        self.dict = json.loads(self.str)
        #self._dict_extend()
        self.dict = _dict_merge(self._extend, self.dict)
