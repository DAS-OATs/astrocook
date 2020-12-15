from .vars import log_seed
from copy import deepcopy as dc
import json
import numpy as np
import pathlib
import wx

class GUILog(object):

    def __init__(self,
                 gui):
        self._gui = gui
        self._indent = 2
        self._init()

    def _init(self):
        self.json = dc(log_seed)
        self.str = json.dumps(self.json, indent=self._indent)


    def _sort(self, sess_list):
        json = {'set_menu': []}
        for s in sess_list:
            for m in self.json['set_menu']:
                if m['cookbook']=='_panel_sess' and m['recipe']=='_on_open':
                    if m['params']['path'] == s.path:
                        json['set_menu'].append(m)
        for m in self.json['set_menu']:
            if m['cookbook']!='_panel_sess' or m['recipe']!='_on_open':
                json['set_menu'].append(m)
        return json


    def _trim(self, menu):
        menu_new = []
        for m in menu:
            if m['cookbook']!='' or m['recipe']!='_refresh' or m!=menu[-1]:
                menu_new.append(m)
        return menu_new


    def append(self, cb, rec, params):
        #if not isinstance(params, list):
        #    params = [params]
        self.trim()
        method_dict = {'cookbook': cb, 'recipe': rec, 'params': dict(params)}
        self.json['set_menu'].append(method_dict)
        self.str = json.dumps(self.json, indent=self._indent)


    def append_full(self, cb, rec, params):
        self.append(cb, rec, params)
        self.close()


    def clear(self):
        self._init()


    def close(self):
        self.append('', '_refresh', {'autosort': False})


    def merge(self, sess_list, sess_sel):
        """
        self.trim()
        for s in sess_list:
            if s != self._gui._sess_sel:
                menu = dc(s.log.json['set_menu'])
                trim = self._trim(menu)
                for m in trim:
                    self.json['set_menu'].append(m)
        self.json = self._sort(sess_list)
        """
        json_new = {'set_menu': []}
        for s in sess_list:
            if s != sess_sel:
                menu = dc(s.log.json['set_menu'])
                trim = self._trim(menu)
                for m in trim:
                    #print(m)
                    if m not in json_new['set_menu']:
                        json_new['set_menu'].append(m)
        #print(json_new)
        self.trim()
        for m in self.json['set_menu']:
            if m not in json_new['set_menu']:
                json_new['set_menu'].append(m)
        self.json = json_new
        #print(self.json)
        #self.json = self._sort(sess_list)
        #print(self.json)
        self.str = json.dumps(self.json, indent=self._indent)

    def merge_full(self, cb, rec, params, sess_list, sess_sel):
        self.append(cb, rec, params)
        self.merge(sess_list, sess_sel)
        self.close()


    def save(self, path):

        root = path[:-5]
        stem = pathlib.PurePath(path[:-5]).parts[-1]

        file = open(root+'.json', "w")
        n = file.write(self.str)
        file.close()


    def update(self, sess_list, sess_sel, sess_item_sel, cb, rec, params):
        sess_sel.log.trim()
        #print(self._gui._sess_sel.log.json)
        if cb=='_panel_sess' and rec=='combine':
            sess_sel.log.append(cb, rec, params)
            sl = [sess_list[s] for s in np.sort(sess_item_sel)]
            sess_sel.log.merge(sl, sess_sel)
        elif cb=='_panel_sess' and rec=='equalize':
            sess_sel.log.append(cb, rec, params)
            sl = [sess_list[s] for s in np.sort(sess_item_sel)]
            sess_sel.log.merge(sl, sess_sel)
        elif cb=='cb' and rec=='rebin':
            sess_sel.log.append(cb, rec, params)
            sl = [sess_list[s] for s in np.sort(sess_item_sel)]
            sess_sel.log.merge(sl, sess_sel)
        else:
            sess_sel.log.append(cb, rec, params)
        sess_sel.log.close()


    def trim(self):
        self.json['set_menu'] = self._trim(self.json['set_menu'])
