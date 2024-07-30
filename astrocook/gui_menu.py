from . import *
from .vars import *
from .gui_dialog import *
#from .session import Session
#from .model import Model
from .model_list import ModelList
from .graph import GraphCursorZSeries
from astropy.io import ascii
import datetime
import logging
import os
import sys
import wx

class GUIMenu(object):

    def __init__(self,
                 gui):
        self._gui = gui
        self._gui._menu = self
        self._params_last = None
        self._togg_set()
        """
        self._menus_togg = {'attr': [menus_assoc[m] \
                                     for m in self._gui._defs.dict['menus']],
                            'title': self._gui._defs.dict['menus']}
        """

    def bar(self):
        bar = wx.MenuBar()

        self._file = GUIMenuFile(self._gui)
        bar.Append(self._file._menu, "File")
        self._edit = GUIMenuEdit(self._gui)
        bar.Append(self._edit._menu, "Edit")
        self._view = GUIMenuView(self._gui)
        bar.Append(self._view._menu, "View")

        for a, t in zip(self._menus_togg['attr'], self._menus_togg['title']):
            setattr(self, a,
                getattr(sys.modules[__name__], 'GUIMenu%s' % t)(self._gui))
            bar.Append(getattr(self, a)._menu, t)
        self._key_list = ['spec', 'lines', 'systs', 'meta', 'defs', 'log',
                          'graph', 'legend', 'norm']

        return bar

    def _create(self, menu, rec, cb, start_id):
        subtitle = ''
        for i, r in enumerate(rec):
            id = start_id+i
            if isinstance(r, str):
                if r == '--':
                    menu.AppendSeparator()
                elif r[0] == '>':
                    submenu = wx.Menu()
                    subtitle = r[2:]
                elif r[0] == '<':
                    menu.AppendSubMenu(submenu, subtitle)
                    subtitle = ''
            if isinstance(r, dict):
                m = menu if subtitle == '' else submenu
                if 'func' in r:
                    setattr(m, r['targ'], {'start_id': id, 'func': r['func'],
                                           'value': r['value']})
                    r['enable'] = self._enable(r['func'], r['value'])
                if 'type' not in r: r['type'] = '_item_method'
                if 'append' not in r: r['append'] = None
                if 'title' not in r: r['title'] = self._get_doc(getattr(cb, r['targ']))

                if r['type'] == '_item':
                    if 'event' not in r: r['event'] = None
                    if 'key' not in r: r['key'] = None
                    if 'enable' not in r: r['enable'] = True
                    self._item(m, id, r['append'], r['title'], r['event'],
                               r['key'], r['enable'])

                if r['type'] == '_item_graph':
                    if 'key' not in r: r['key'] = None
                    if 'enable' not in r: r['enable'] = False
                    if 'dlg_mini' not in r: r['dlg_mini'] = None
                    if 'targ' not in r: r['targ'] = None
                    if 'alt_title' not in r: r['alt_title'] = r['title']
                    self._item_graph(m, id, r['append'], r['title'], r['key'],
                                     r['enable'], r['dlg_mini'], r['targ'],
                                     r['alt_title'])

                if r['type'] == '_item_method':
                    if 'enable' not in r: r['enable'] = False
                    if 'obj' not in r: r['obj'] = None
                    self._item_method(m, id, r['append'], r['title'], r['targ'],
                                      r['enable'], r['obj'])

        return 0

    def _enable(self, func, value):
        return getattr(len(self._gui._sess_item_sel), func)(value)


    def _get_doc(self, method):
        full = inspect.getdoc(method)
        split = full.split('@')
        return [s[6:-1] for s in split if s[0:5]=='brief'][0].replace('\n', ' ')


    def _item(self, menu, id, append, title, event, key=None, enable=True):
        if key is not None:
            item = wx.MenuItem(menu, id, title, kind=wx.ITEM_CHECK)
            #item.Check(False)
            item.key = key
        else:
            item = wx.MenuItem(menu, id, title)

        self._gui._panel_sess.Bind(wx.EVT_MENU, event, item)
        menu.Append(item)
        if append is not None:
            if isinstance(append, list):
                for a in append:
                    getattr(self._gui, '_menu_'+a+'_id').append(id)
            else:
                getattr(self._gui, '_menu_'+append+'_id').append(id)
            item.Enable(False)
        else:
            item.Enable(enable)


    def _item_graph(self, menu, id, append, title, key=None, enable=False,
                    dlg_mini=None, targ=None, alt_title=None):
        if alt_title == None: alt_title = title
        item = wx.MenuItem(menu, id, title, kind=wx.ITEM_CHECK)
        item.key = key
        if targ == GraphCursorZSeries:
            self._gui._cursor = item
        if dlg_mini == "graph":
            self._gui._graph_elem = item
        self._gui._panel_sess.Bind(
            wx.EVT_MENU,
            lambda e: self._on_graph(e, alt_title, key, item, dlg_mini, targ),
            item)
        menu.Append(item)
        if append is not None:
            if isinstance(append, list):
                for a in append:
                    getattr(self._gui, '_menu_'+a+'_id').append(id)
            else:
                getattr(self._gui, '_menu_'+append+'_id').append(id)
            item.Enable(False)
        else:
            item.Enable(enable)

    def _item_method(self, menu, id, append, title, targ, enable=False,
                     obj=None):
        item = wx.MenuItem(menu, id, title+'...')
        self._gui._panel_sess.Bind(
            wx.EVT_MENU,
            lambda e: self._on_dialog(e, title, targ, obj), item)
        menu.Append(item)
        if append is not None:
            if isinstance(append, list):
                for a in append:
                    getattr(self._gui, '_menu_'+a+'_id').append(id)
            else:
                getattr(self._gui, '_menu_'+append+'_id').append(id)
            item.Enable(False)
        else:
            item.Enable(enable)

    def _on_dialog(self, event, title, attr, obj=None):
        if isinstance(attr, list):
            dlg = GUIDialogMethods(self._gui, title, attr, obj,
                                   params_last=self._params_last)
        else:
            dlg = GUIDialogMethod(self._gui, title, attr, obj,
                                  params_last=self._params_last)
        self._params_last = dlg._params

    def _on_dialog_mini_defs(self, event, title, targ, log=True):
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_menu', '_on_dialog_mini_defs',
                                 {'event': None, 'title': title, 'targ': targ})
        if hasattr(self._gui, '_dlg_mini_graph'):
            self._gui._dlg_mini_defs._refresh()
        else:
            dlg = GUIDialogMiniDefaults(self._gui, title)


    def _on_dialog_mini_graph(self, event, title, targ, log=True):
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_menu', '_on_dialog_mini_graph',
                                 {'event': None, 'title': title, 'targ': targ})
        if hasattr(self._gui, '_dlg_mini_graph'):
            try:
                self._gui._dlg_mini_graph._refresh()
            except:
                dlg = GUIDialogMiniGraph(self._gui, title)
        else:
            dlg = GUIDialogMiniGraph(self._gui, title)

    def _on_dialog_mini_log(self, event, title, targ):
        dlg = GUIDialogMiniLog(self._gui, title)


    def _on_dialog_mini_meta(self, event, title, targ, log=True):
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_menu', '_on_dialog_mini_meta',
                                 {'event': None, 'title': title, 'targ': targ})
        dlg = GUIDialogMiniMeta(self._gui, title)


    def _on_dialog_mini_systems(self, event, title, targ):
        dlg = GUIDialogMiniSystems(self._gui, title, targ)


    def _on_graph(self, event, title, key, item, dlg_mini, targ):
        sel = self._gui._graph_main._sel
        if key is not None:
            if key in sel:
                sel.remove(key)
            else:
                sel.append(key)
        #item.IsChecked() == False
        self._gui._refresh(init_tab=False)
        if dlg_mini is not None:
            #self._gui._cursor = item
            if item.IsChecked() or key is None:
                if not hasattr(self._gui, '_dlg_mini_'+dlg_mini) \
                    or getattr(self._gui, '_dlg_mini_'+dlg_mini) == None \
                    or not getattr(getattr(self._gui, '_dlg_mini_'+dlg_mini), '_shown'):
                    getattr(self, '_on_dialog_mini_'+dlg_mini)\
                        (event, title, targ)
                gui_dlg_mini = getattr(self._gui, '_dlg_mini_'+dlg_mini)
                gui_dlg_mini._shown = True
                if dlg_mini == 'systems':
                    gui_dlg_mini._on_apply(event, refresh=True)
                    gui_dlg_mini._cursor_button.SetLabel("Hide cursor")
                else:
                    gui_dlg_mini._on_apply(event, refresh=False)
            else:
                if hasattr(self._gui, '_dlg_mini_'+dlg_mini):
                    gui_dlg_mini = getattr(self._gui, '_dlg_mini_'+dlg_mini)
                    gui_dlg_mini._shown = False
                    gui_dlg_mini._on_cancel(event)
                    if dlg_mini == 'systems':
                        gui_dlg_mini._cursor_button.SetLabel("Show cursor")


    def _on_open(self, event, path=None, wildcard=None,
                 action='_on_open_session'):
        """ Behaviour for Session > Open """

        if path is None:
            if hasattr(self._gui, '_path'):
                path=self._gui._path
            else:
                path='.'
        if wildcard is None:
            wildcard = "Astrocook sessions (*.acs)|*.acs|" \
                       "FITS files (*.fits)|*.fits|" \
                       "JSON files (*.json)|*.json|" \
                       "CSV files (*.csv)|*.csv|" \
                       "Data files (*.dat)|*.dat|" \
                       "Text files (*.txt)|*.txt"
        with wx.FileDialog(self._gui._panel_sess, "Open file", path,
                           wildcard=wildcard,
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return
            self._gui._path = fileDialog.GetPath()
        self._gui._panel_sess._open_path = self._gui._path
        if self._gui._path[-4:] == 'json':
            self._gui._panel_sess._open_rec = 'json_load'
            self._gui._panel_sess.json_load(os.path.realpath(self._gui._path))
        else:
            self._gui._panel_sess._open_rec = '_on_open'
            self._gui._panel_sess._on_open(os.path.realpath(self._gui._path))

    def _refresh(self, init_bar=False):
        # Nested loops! WOOOO!
        sess = self._gui._sess_sel
        sel = self._gui._graph_main._sel
        self._togg_set()
        if init_bar:
            self._gui._panel_sess.SetMenuBar(self.bar())
            for k in self._key_list:
                try:
                    getattr(self._gui, '_tab_'+k)._on_close(None)
                except:
                    pass

        for a in seq_menu:  # from .vars
            for i in getattr(self._gui, '_menu_'+a+'_id'):
                for m in ['_file', '_edit', '_view']+self._menus_togg['attr']:
                    try:
                        item = getattr(self, m)._menu.FindItemById(i)
                        if m == '_view' and item.IsCheckable() \
                            and item.key not in self._key_list:
                            item.Check(False)

                        if hasattr(sess, a):
                            cond = getattr(sess, a) != None
                        else:
                            if hasattr(sess, 'systs') and sess.systs != None:
                                cond = a in sess.systs.t.colnames \
                                           or a in sess.spec.t.colnames
                            else:
                                cond = a in sess.spec.t.colnames
                        if cond:
                            item.Enable(True)
                            if m == '_view' and item.IsCheckable():
                                if item.key not in self._key_list:
                                    item.Check(item.key in sel)
                        else:
                            item.Enable(False)

                    except:
                        pass


    def _sel_graph_cols(self, cols=graph_cols_sel):
        """ @brief Select graph columns
        @details Select columns to be displayed on graph
        @param cols Columns to be displayed
        @return 0
        """
        self._gui._graph_main._cols_sel = cols

        return 0

    def _togg_set(self):
        self._menus_togg = {'attr': [menus_assoc[m] \
                             for m in self._gui._defs.dict['menus']],
                            'title': self._gui._defs.dict['menus']}


class GUIMenuCook(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=11000,
                 **kwargs):
        super(GUIMenuCook, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to Cook menu here
        self._item(self._menu, start_id+1, None,
                   "Fit CIV forest in all QSOs...",
                   self._on_civ_full)
        #self._item(self._menu, start_id+2, 'spec', "Test...", self._on_test)

    def _on_civ_full(self, event):
        from .cookbook import Cookbook
        from .workflow import Workflow
        wf = Workflow(self._gui, Cookbook())
        wf.civ_full()


class GUIMenuAbsorbers(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=7000,
                 **kwargs):
        super(GUIMenuAbsorbers, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        self._rec = [{'targ': 'systs_new_from_like', 'append': 'cont'},
                     {'targ': 'systs_new_from_lines', 'append': 'lines'},
                     {'targ': 'systs_complete', 'append': ['z0']},
                     {'targ': 'systs_complete_from_z', 'append': ['z0']},
                     {'targ': 'lya_fit', 'append': 'cont'},
                     '> Other',
                     {'targ': 'cands_find', 'append': 'z0'},
                     {'targ': 'systs_improve', 'append': 'z0'},
                     '<',
                     '--',
                     {'targ': 'systs_fit', 'append': 'z0'},
                     {'targ': 'systs_clean', 'append': 'z0'},
                     {'targ': 'mods_recreate', 'append': 'z0'},
                     '> Other',
                     {'targ': 'systs_snr', 'append': 'z0'},
                     {'targ': 'systs_select', 'append': 'z0'},
                     {'targ': 'comp_extract', 'append': 'z0'},
                     {'targ': 'systs_merge', 'append': 'z0'},
                     '<',
                     '--',
                     {'targ': 'feats', 'append': 'z0'},
                     '--',
                     {'targ': 'mods_ccf_max', 'append': 'z0'},
                     {'targ': 'systs_sigmav', 'append': 'z0'},
                    ]

        from .cookbook_absorbers import CookbookAbsorbers as cbc
        self._cb = cbc()

        self._create(self._menu, self._rec, self._cb, start_id)


class GUIMenuContinuum(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=6000,
                 **kwargs):
        super(GUIMenuContinuum, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        self._rec = [{'targ': 'flux_clip', 'append': 'spec'},
                     '--',
                     {'targ': 'lines_find', 'append': 'spec'},
                     {'targ': 'nodes_cont', 'append': 'spec'},
                     {'targ': 'lines_update', 'append': 'z0'},
                     '> Other',
                     {'targ': 'peaks_find', 'append': 'spec'},
                     {'targ': 'nodes_extract', 'append': 'cont'},
                     {'targ': 'nodes_clean', 'append': 'lines'},
                     {'targ': 'nodes_interp', 'append': 'nodes'},
                     '<',
                     '--',
                     {'targ': 'lya_corr', 'append': 'spec'},
                     {'targ': 'abs_cont', 'append': 'spec'},
                     ]

        from .cookbook_continuum import CookbookContinuum as cbc
        self._cb = cbc()

        self._create(self._menu, self._rec, self._cb, start_id)


class GUIMenuEdit(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=2000,
                 **kwargs):
        super(GUIMenuEdit, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        self._rec = [{'targ': 'x_convert', 'append': 'spec'},
                     {'targ': 'y_convert', 'append': 'spec'},
                     '--',
                     {'targ': 'shift_bary', 'append': 'spec'},
                     {'targ': 'shift_to_rf', 'append': 'spec'},
                     {'targ': 'shift_from_rf', 'append': 'spec'},
                     '--',
                     {'targ': 'struct_import', 'append': 'spec'},#'func': '__gt__', 'value': 0},
                     {'targ': 'struct_modify', 'append': 'spec'},#'func': '__gt__', 'value': 0},
                     ]


        from .cookbook_edit import CookbookEdit as cbc
        self._cb = cbc()

        self._create(self._menu, self._rec, self._cb, start_id)


class GUIMenuFile(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=1000,
                 **kwargs):
        super(GUIMenuFile, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()
        self._gui._menu_file = self
        self._start_id = start_id

        # Add items to File menu here
        self._item(self._menu, start_id, None, "Open...\tCtrl+O",
                   lambda e: self._on_open(e, **kwargs))
        self._menu.AppendSeparator()
        self._item(self._menu, start_id+101, None, "Save session...\tCtrl+S",
                   lambda e: self._on_save(e, **kwargs))
        self._item(self._menu, start_id+102, 'systs', "Save session with models...\tCtrl+S",
                   lambda e: self._on_save(e, models=True, **kwargs))
        self._item(self._menu, start_id+103, None, "Save spectrum as PDF...\tCtrl+S",
                   lambda e: self._on_save_pdf(e, **kwargs))
        self._menu.AppendSeparator()
        self._item(self._menu, start_id+400, None, "Quit\tCtrl+Q",
                   self._gui._panel_sess._on_close)

    def _on_combine(self, event):
        self._gui._panel_sess._combine()


    def _on_save(self, event, path=None, models=False):
        """ Behaviour for Session > Save """

        if path is None:
            if hasattr(self._gui, '_path'):
                path=os.path.basename(self._gui._path)
            else:
                path='.'
        name = self._gui._sess_sel.name
        with wx.FileDialog(self._gui._panel_sess, "Save session", path, name,
                           wildcard="Astrocook session (*.acs)|*.acs",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return

            path = fileDialog.GetPath()
            dir = fileDialog.GetDirectory()
            logging.info("I'm saving session %s..." % path)
            self._gui._sess_sel.save(path, models)

    def _on_save_pdf(self, event, path=None):
        if path is None:
            if hasattr(self._gui, '_path'):
                path=os.path.basename(self._gui._path)
            else:
                path='.'
        name = self._gui._sess_sel.name
        with wx.FileDialog(self._gui._panel_sess, "Save spectrum as PDF", path, name,
                           wildcard="PDF (*.pdf)|*.pdf",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return

            path = fileDialog.GetPath()
            dir = fileDialog.GetDirectory()
            self._gui._sess_sel.save_pdf(path)





class GUIMenuFlux(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=5100,
                 **kwargs):
        super(GUIMenuFlux, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        self._rec = [{'targ': 'y_scale', 'append': 'spec'},
                     {'targ': 'y_scale_med', 'append': 'spec'},
                     {'targ': 'y_scale_x', 'append': 'spec'},
                     '--',
                     {'targ': 'deredden', 'append': 'spec'},
                     '--',
                     {'targ': 'mags_adjust', 'append': 'spec'},
                     ]

        from .cookbook_flux import CookbookFlux as cbc
        self._cb = cbc()

        self._create(self._menu, self._rec, self._cb, start_id)


class GUIMenuGeneral(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=4000,
                 **kwargs):
        super(GUIMenuGeneral, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()


        self._rec = [{'targ': 'region_extract', 'append': 'spec'},
                     {'targ': 'part_extract', 'append': 'spec'},
                     '--',
                     {'targ': 'equalize', 'func': '__eq__', 'value': 2},
                     {'targ': 'combine', 'func': '__gt__', 'value': 1},
                     '--',
                     {'targ': 'outliers_clean', 'append': 'spec'},
                     '--',
                     {'targ': 'x_mask', 'append': 'spec'},
                     {'targ': 'mask_cond', 'append': 'spec'},
                     {'targ': 'sky_mask', 'append': 'spec'},
                     {'targ': 'telluric_mask', 'append': 'spec'},
                     '--',
                     {'targ': 'snr_est', 'append': 'spec'},
                     {'targ': 'resol_est', 'append': 'spec'},
                     {'targ': 'rms_est', 'append': 'spec'},
                     {'targ': 'dx_est', 'append': 'spec'},
                     '--',
                     {'targ': 'rebin', 'append': 'spec'},
                     {'targ': 'gauss_convolve', 'append': 'spec'},
                     {'targ': 'flux_ccf', 'append': 'spec'},
                     {'targ': 'flux_ccf_stats', 'append': 'spec'},
                     ]

        from .cookbook_general import CookbookGeneral as cbc
        self._cb = cbc()

        self._create(self._menu, self._rec, self._cb, start_id)


class GUIMenuSynthetic(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=8000,
                 **kwargs):
        super(GUIMenuSynthetic, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        self._rec = [{'targ': 'spec_from_struct', 'append': 'spec'},
                     {'targ': 'spec_from_systs', 'append': 'systs'},
                     {'targ': 'spec_from_systs_random', 'append': 'systs'},
                     {'targ': 'systs_random', 'append': 'systs'},
                     ]

        from .cookbook_synthetic import CookbookSynthetic as cbc
        self._cb = cbc()

        self._create(self._menu, self._rec, self._cb, start_id)


class GUIMenuTemplates(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=9000,
                 **kwargs):
        super(GUIMenuTemplates, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        self._rec = [{'targ': 'bb', 'append': 'spec'},
                     {'targ': 'pl', 'append': 'spec'},
                     ]

        from .cookbook_templates import CookbookTemplates as cbc
        self._cb = cbc()

        self._create(self._menu, self._rec, self._cb, start_id)


class GUIMenuView(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=3000,
                 **kwargs):
        super(GUIMenuView, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()
        self._menu_view = self
        self._gui._menu_view = self
        #tab_id = [start_id+1, start_id+2, start_id+3, start_id+4, start_id+5]
        tab_id = [start_id+0, start_id+1, start_id+2]
        dlg_id = [start_id+5, start_id+6, start_id+7, start_id+8, start_id+9]
        self._gui._menu_tab_id = tab_id
        self._gui._menu_dlg_id = dlg_id

        self._rec = [{'type': '_item',
                      'event': lambda e: self._on_tab(e, 'spec'),
                      'title': "Spectrum table", 'append': 'spec', 'key': 'spec'},
                     {'type': '_item',
                      'event': lambda e: self._on_tab(e, 'lines'),
                      'title': "Line table", 'append': 'lines', 'key': 'lines'},
                     {'type': '_item',
                      'event': lambda e: self._on_tab(e, 'systs'),
                      'title': "System table", 'append': 'systs', 'key': 'systs'},
                     {'type': '_item', 'event': self._on_compress,
                      'title': "Compress system table", 'append': 'systs'},
                     '--',
                     {'type': '_item', 'title': "Session metadata",
                      'event': lambda e: self._on_dlg_mini(e, 'meta'),
                      'key': 'meta', 'append': 'spec', 'dlg_mini': 'meta'},
                     {'type': '_item', 'title': "Session defaults",
                      'event': lambda e: self._on_dlg_mini(e, 'defs'),
                      'key': 'defs', 'append': 'spec', 'dlg_mini': 'defs'},
                     {'type': '_item', 'title': "Session log",
                      'event': lambda e: self._on_dlg_mini(e, 'log'),
                      'key': 'log', 'append': 'spec', 'dlg_mini': 'log'},
                     {'type': '_item', 'title': "Graph elements",
                      'event': lambda e: self._on_dlg_mini(e, 'graph'),
                      'key': 'graph', 'append': 'spec', 'dlg_mini': 'graph'},
                     {'type': '_item_graph', 'title': "Redshift cursor",
                      'key': 'cursor_z_series', 'append': 'spec',
                      'dlg_mini': 'systems', 'targ': GraphCursorZSeries,
                      'alt_title': "System controls"},
                     '--',
                     {'type': '_item', 'event': self._on_logx,
                      'title': "Toggle log x axis", 'append': 'spec'},
                     {'type': '_item', 'event': self._on_logy,
                      'title': "Toggle log y axis", 'append': 'spec'},
                     {'type': '_item', 'event': self._on_norm,
                      'title': "Toggle normalization", 'append': 'cont'},
                     '> Toggle graph add-ons',
                     {'type': '_item_graph', 'title': "Saturated H2O regions",
                      'key': 'spec_h2o_reg', 'append': 'spec'},
                     {'type': '_item', 'event': self._on_legend,
                      'key': 'legend', 'title': "Legend", 'append': 'spec'},
                     '<',
                     '--',
                     {'targ': 'z_ax', 'append': 'spec'},
                     {'type': '_item', 'event': self._on_z_ax_remove,
                      'title': "Hide redshift axis", 'append': 'spec'},
                     '--',
                     ]

        from .cookbook_view import CookbookView as cbc
        self._cb = cbc()

        self._create(self._menu, self._rec, self._cb, start_id)

        return


    def _on_compress(self, event, log=True):
        if self._menu.GetLabel(self._start_id+5) == "Compress system table":
            self._menu.SetLabel(self._start_id+5, "Uncompress system table")
        else:
            self._menu.SetLabel(self._start_id+5, "Compress system table")

        self._gui._sess_sel.systs._compress()
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_menu_view', '_on_compress',
                                 {'event': None, 'log': False})
        self._gui._refresh()

    def _on_ima(self, event, obj):
        method = '_ima_'+obj
        getattr(self._gui, method)._on_view(event)

    def _on_legend(self, event):
        self._gui._graph_main._legend = ~self._gui._graph_main._legend
        self._gui._refresh()

    def _on_logx(self, event, log=False):
        self._gui._graph_main._logx = ~self._gui._graph_main._logx
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_menu_view', '_on_logx',
                                 {'event': None, 'log': False})
        self._gui._refresh()

    def _on_logy(self, event, log=False):
        self._gui._graph_main._logy = ~self._gui._graph_main._logy
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_menu_view', '_on_logy',
                                 {'event': None, 'log': False})
        self._gui._refresh()

    def _on_norm(self, event, log=True):
        self._gui._graph_main._norm = ~self._gui._graph_main._norm
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_menu_view', '_on_norm',
                                 {'event': None, 'log': False})
        self._gui._refresh()

    def _on_dlg_mini(self, event, obj, check=None, log=True):
        index = ['meta', 'defs', 'log', 'graph'].index(obj)
        title = ['Session metadata', 'Session defaults', 'Session log',
                 'Graph elements'][index]
        item = self._menu.FindItemById(self._gui._menu_dlg_id[index])

        if check is not None:
            view = check
            item.Check(view)
        else:
            view = item.IsChecked()

        sess = self._gui._sess_sel
        if view:
            getattr(self, '_on_dialog_mini_'+obj)(event, title, None)
            setattr(getattr(self._gui, '_dlg_mini_'+obj), '_shown', True)
        else:
            setattr(getattr(self._gui, '_dlg_mini_'+obj), '_shown', False)
            getattr(self._gui, '_dlg_mini_'+obj)._on_cancel(event)


    def _on_tab(self, event, obj, check=None, log=True):
        method = '_tab_'+obj
        index = ['spec', 'lines', 'systs'].index(obj)
        item = self._menu.FindItemById(self._gui._menu_tab_id[index])

        if check is not None:
            view = check
            item.Check(view)
        else:
            view = item.IsChecked()

        sess = self._gui._sess_sel
        if view:
            if log:
                sess.log.append_full('_menu_view', '_on_tab',
                    {'event': None, 'obj': obj, 'check': True, 'log': False})
            setattr(self._gui, '_tab_'+obj+'_shown', True)
            getattr(self._gui, method)._on_view(event)
        else:
            if log:
                sess.log.append_full('_menu_view', '_on_tab',
                    {'event': None, 'obj': obj, 'check': False, 'log': False})
            setattr(self._gui, '_tab_'+obj+'_shown', False)
            getattr(self._gui, method)._on_close(event)
        if hasattr(self._gui, '_dlg_mini_log') \
            and self._gui._dlg_mini_log._shown:
            self._gui._dlg_mini_log._refresh()


    def _on_z_ax_remove(self, event, log=False):
        try:
            delattr(self._gui._sess_sel, '_ztrans')
            self._gui._refresh()
        except:
            pass
