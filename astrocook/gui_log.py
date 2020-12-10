from .vars import log_seed
from copy import deepcopy as dc
import json

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
            if m['cookbook']!='' or m['recipe']!='_refresh' or m['params']!={}\
                or m!=menu[-1]:
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


    def close(self):
        self.append('', '_refresh', {})


    def merge(self, sess_list):
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
            if s != self._gui._sess_sel:
                menu = dc(s.log.json['set_menu'])
                trim = self._trim(menu)
                for m in trim:
                    if m not in json_new['set_menu']:
                        json_new['set_menu'].append(m)
        self.trim()
        for m in self.json['set_menu']:
            if m not in json_new['set_menu']:
                json_new['set_menu'].append(m)
        self.json = json_new
        self.json = self._sort(sess_list)
        self.str = json.dumps(self.json, indent=self._indent)


    def trim(self):
        self.json['set_menu'] = self._trim(self.json['set_menu'])
