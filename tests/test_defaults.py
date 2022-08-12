import pytest
import wx

app = wx.App(False)
from astrocook.gui import GUI
gui = GUI()

from astrocook.defaults import Defaults, _dict_merge
defs = Defaults(gui)
defs.open(path='tests', file='test_defaults.json')

def test__dict_merge():
    dict1 = {'a':{'a0':0, 'a1':1}, 'b':{'b0':0, 'b1':1}}
    dict2 = {'b':{'b2':2, 'b3':3}, 'c':{'c0':0, 'c1':1}}
    assert _dict_merge(dict1, dict2) == {'a':{'a0':0, 'a1':1},
                                         'b':{'b0':0, 'b1':1,'b2':2, 'b3':3},
                                         'c':{'c0':0, 'c1':1}}


class TestDefaults:

    def test_open(self):
        dict = {'key': 'value--'}
        dict = _dict_merge(defs._extend, dict)
        assert defs.dict == dict


    def test_update(self):
        str = '{"key": "value--"}'
        defs.update(str)
        dict = {'key': 'value--'}
        dict = _dict_merge(defs._extend, dict)
        assert defs.dict == dict
