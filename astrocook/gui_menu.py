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
import wx

class GUIMenu(object):

    def __init__(self,
                 gui):
        self._gui = gui
        self._gui._menu = self
        self._params_last = None

    def bar(self):
        bar = wx.MenuBar()
        self._file = GUIMenuFile(self._gui)
        self._edit = GUIMenuEdit(self._gui)
        self._view = GUIMenuView(self._gui)
        self._recipes = GUIMenuRecipes(self._gui)
        self._courses = GUIMenuCourses(self._gui)
        self._cook = GUIMenuCook(self._gui)
        self._key_list = ['spec', 'lines', 'systs', 'legend', 'norm']
        bar.Append(self._file._menu, "File")
        bar.Append(self._edit._menu, "Edit")
        bar.Append(self._view._menu, "View")
        bar.Append(self._recipes._menu, "Recipes")
        #bar.Append(self._courses._menu, "Courses")
        bar.Append(self._courses._menu, "Set menus")
        #bar.Append(self._cook._menu, "Cook")
        return bar

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
            getattr(self._gui, '_menu_'+append+'_id').append(id)
            item.Enable(False)
        else:
            item.Enable(enable)
        return item

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

    def _on_dialog_mini_graph(self, event, title, targ, log=True):
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_menu', '_on_dialog_mini_graph',
                                 {'event': None, 'title': title, 'targ': targ})
        if hasattr(self._gui, '_dlg_mini_graph'):
            self._gui._dlg_mini_graph._refresh()
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
                gui_dlg_mini._on_apply(event, refresh=False)
                if dlg_mini == 'systems':
                    gui_dlg_mini._cursor_button.SetLabel("Hide cursor")
            else:
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
                       "Text files (*.txt)|*.txt"
        with wx.FileDialog(self._gui._panel_sess, "Open file", path,
                           wildcard=wildcard,
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return
            self._gui._path = fileDialog.GetPath()
            """
            try:
                getattr(self, action)(self._gui._path)
            except:
                getattr(self._gui._panel_sess, action)(self._gui._path)
            """
        """
        try:
            getattr(self, action)(self._gui._path)
        except:
            getattr(self._gui._panel_sess, action)(self._gui._path)
        """
        self._gui._panel_sess._open_path = self._gui._path
        if self._gui._path[-4:] == 'json':
            self._gui._panel_sess._open_rec = 'json_load'
            self._gui._panel_sess.json_load(os.path.realpath(self._gui._path))
        else:
            self._gui._panel_sess._open_rec = '_on_open'
            self._gui._panel_sess._on_open(os.path.realpath(self._gui._path))

    """
    def _on_open_session(self, path):
        name = path.split('/')[-1].split('.')[0]
        #logging.info("I'm loading session %s..." % path)
        sess = Session(gui=self._gui, path=path, name=name)
        self._gui._panel_sess._on_add(sess, open=True)
        if sess._open_twin:
            #logging.info("I'm loading twin session %s..." % path)
            sess = Session(gui=self._gui, path=path, name=name, twin=True)
            self._gui._panel_sess._on_add(sess, open=True)
        #self._gui._path = path
    """

    def _refresh(self):
        # Nested loops! WOOOO!
        sel = self._gui._graph_main._sel

        for a in seq_menu:  # from .vars
            for i in getattr(self._gui, '_menu_'+a+'_id'):
                for m in ['_edit', '_view', '_recipes', '_courses', 'cook']:
                    try:
                        item = getattr(self, m)._menu.FindItemById(i)
                        if m == '_view' and item.IsCheckable() \
                            and item.key not in self._key_list:
                            item.Check(False)
                        if hasattr(self._gui._sess_sel, a):
                            cond = getattr(self._gui._sess_sel, a) != None
                        else:
                            try:
                                cond = a in self._gui._sess_sel.systs.t.colnames
                            except:
                                cond = a in self._gui._sess_sel.spec.t.colnames
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


class GUIMenuCook(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=6000,
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

    """
    def _on_civ_full_old(self, event):

        targ_list = ascii.read('/data/cupani/CIV/targets_3.csv')

        for l in targ_list:
            time_start = datetime.datetime.now()
            t = l['name']
            zem = l['z']
            xmin = l['lambdamin']
            xmax = l['lambdamax']
            self._gui._panel_sess._on_open('/data/cupani/CIV/reduced/'+t\
                                           +'.fits')


            sess_start = self._gui._sess_sel
            if sess_start.spec.meta['object'] == 'J2123-0050':
                sess = sess_start.region_extract(xmin=xmin, xmax=xmax)
            else:
                sess = sess_start

            sess.gauss_convolve(std=10)
            sess.peaks_find(kappa=3.0)
            sess.lines._t.remove_rows(sess.lines.y == 0)
            if np.mean(sess.spec._t['y'])<1 and np.std(sess.spec._t['y'])<1:
                sess.spec._t['cont'] = [1] * len(sess.spec._t)*sess.spec.y.unit
            if 'cont' not in sess.spec._t.colnames:
                sess.nodes_extract(delta_x=1000)
                sess.nodes_interp()

            sess_reg = sess.region_extract(xmin=xmin, xmax=xmax)
            self._gui._panel_sess._on_add(sess_reg, open=False)
            sess_center = dc(sess_reg)
            sess_center.add_syst_from_lines(z_end=20, maxfev=10)#series='unknown')

            sess_reg.lines.t['x'] = (1+sess_center.systs.t['z'])\
                                    *xem_d['Ly_a'].to(sess_reg.spec.x.unit)
            sess_reg.lines.t['logN'] = sess_center.systs.t['logN']

            #sess_reg.add_syst_from_lines(series='SiII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='SiIV', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            sess_reg.add_syst_from_lines(series='CIV', logN=None, b=20.0,
                                         dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='FeII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='MgII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)

            sess_reg.add_syst_from_resids(chi2r_thres=2.0, logN=13.0, b=10.0,
                                          maxfev=10)

            sess_reg.compl_syst(n=10)#, z_start=2.128, z_end=2.1372)
            sess_reg.add_syst_slide(col='deabs')#, z_start=1.6, z_end=1.61)
            sess_reg.syst_merge()
            self._gui._refresh()
            sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_'
                          +datetime.date.today().isoformat()+'.xxx')
            sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_latest.xxx')
            time_end = datetime.datetime.now()
            print("%s; computation time: %s" \
                  % (datetime.datetime.now(), time_end-time_start))

    def _on_test(self, event):
        targ_list = ascii.read('/data/cupani/CIV/targets_prova.csv')

        for l in targ_list:
            time_start = datetime.datetime.now()
            t = l['name']
            zem = l['z']
            xmin = l['lambdamin']
            xmax = l['lambdamax']
            #xmin = 381
            #xmax = 382.5
            #xmin = 592
            #xmax = 598
            #xmin = 625
            #xmax = 630
            #xmin = 480
            #xmax = 490
            self._gui._panel_sess._on_open('/data/cupani/CIV/reduced/'+t\
                                           +'.fits')
            sess_start = self._gui._sess_sel
            if sess_start.spec.meta['object'] == 'J2123-0050':
                sess = sess_start.region_extract(xmin=xmin, xmax=xmax)
            else:
                sess = sess_start

            sess.gauss_convolve(std=10)
            sess.peaks_find(kappa=3.0)
            sess.lines._t.remove_rows(sess.lines.y == 0)
            if np.mean(sess.spec._t['y'])<1 and np.std(sess.spec._t['y'])<1:
                sess.spec._t['cont'] = [1] * len(sess.spec._t)*sess.spec.y.unit
            if 'cont' not in sess.spec._t.colnames:
                sess.nodes_extract(delta_x=1000)
                sess.nodes_interp()

            sess_reg = sess.region_extract(xmin=xmin, xmax=xmax)
            self._gui._panel_sess._on_add(sess_reg, open=False)
            #sess_center = dc(sess_reg)
            #sess_center.add_syst_from_lines(z_end=20, maxfev=10)#series='unknown')
            #sess_reg.lines.t['x'] = (1+sess_center.systs.t['z'])\
            #                        *xem_d['Ly_a'].to(sess_reg.spec.x.unit)
            #sess_reg.lines.t['logN'] = sess_center.systs.t['logN']
            #sess_reg.add_syst_from_lines(series='SiII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='SiIV', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            sess_reg.add_syst_from_lines(series='CIV', logN=None, b=20.0,
                                         dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='FeII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='MgII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)

            sess_reg.add_syst_from_resids(chi2r_thres=1.0, logN=13.0, b=2.0,
                                          maxfev=10)
            #sess_reg.compl_syst(n=10)#, z_start=2.128, z_end=2.1372)
            #sess_reg.add_syst_slide(col='deabs')#, z_start=1.6, z_end=1.61)
            #sess_reg.syst_merge()
            self._gui._refresh()
            #sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_'
            #              +datetime.date.today().isoformat()+'.xxx')
            #sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_latest.xxx')
            time_end = datetime.datetime.now()
            print("%s; computation time: %s" \
                  % (datetime.datetime.now(), time_end-time_start))
    """

class GUIMenuEdit(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=2000,
                 **kwargs):
        super(GUIMenuEdit, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()
        self._start_id = start_id

        # Add items to Edit menu here
        #print(len(self._gui._sess_list), len(self._gui._sess_item_sel))
        self._item_method(self._menu, start_id+300, None,
                          "Equalize sessions", 'equalize',
                          enable=len(self._gui._sess_item_sel)==2,
                          obj=self._gui._panel_sess)
        self._item_method(self._menu, start_id+301, None,
                          "Combine sessions", 'combine',
                          enable=len(self._gui._sess_item_sel)>1,
                          obj=self._gui._panel_sess)
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+310, None,
                          "Import structure", 'struct_import',
                          enable=len(self._gui._sess_list)>0,
                          obj=self._gui._panel_sess)
        self._item_method(self._menu, start_id+311, None,
                          "Modify structures", 'struct_modify2',
                          enable=len(self._gui._sess_list)>0,
                          obj=self._gui._panel_sess)
        self._item_method(self._menu, start_id+312, None,
                          "Synthetic spectrum from structure",
                          'spec_from_struct',
                          enable=len(self._gui._sess_list)>0,
                          obj=self._gui._panel_sess)
        submenu = wx.Menu()
        self._item_method(submenu, start_id+312, 'spec',
                          "Blackbody", 'bb')
        self._item_method(submenu, start_id+313, 'spec',
                          "Power-law", 'pl')
        self._menu.AppendSubMenu(submenu, "Apply template")
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+320, 'spec',
                          "Extract region", 'region_extract')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+330, 'spec',
                          "Convert x axis", 'x_convert')
        self._item_method(self._menu, start_id+331, 'spec',
                          "Convert y axis", 'y_convert')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+340, 'spec',
                          "Scale y axis by median", 'y_scale_med')
        self._item_method(self._menu, start_id+341, 'spec',
                          "Scale y axis", 'y_scale')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+350, 'spec',
                          "Shift to barycentric frame", 'shift_bary')
        self._item_method(self._menu, start_id+351, 'spec',
                          "Shift to rest frame", 'shift_to_rf')
        self._item_method(self._menu, start_id+352, 'spec',
                          "Shift from rest frame", 'shift_from_rf')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+360, 'spec',
                          "Deredden", 'deredden')


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
        self._item(self._menu, start_id+101, None, "Save...\tCtrl+S",
                   lambda e: self._on_save(e, **kwargs))
        self._menu.AppendSeparator()
        self._item(self._menu, start_id+400, None, "Quit\tCtrl+Q",
                   self._gui._panel_sess._on_close)

    def _on_combine(self, event):
        self._gui._panel_sess._combine()


    def _on_save(self, event, path=None):
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
            self._gui._sess_sel.save(path)
            """
            try:
                acs = self
                self.IO.acs_write(self.acs, name, dir)

            except IOError:
                wx.LogError("Cannot save session '%s'." % name)
            """


class GUIMenuRecipes(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=4000,
                 **kwargs):
        super(GUIMenuRecipes, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()
        #sess = self._gui._sess_sel


        # Add items to Recipes menu here
        self._item_method(self._menu, start_id+100, 'spec',
                          "Create spectral mask", 'mask')
        self._item_method(self._menu, start_id+101, 'spec',
                          "Mask telluric absorption", 'telluric_mask')
        self._item_method(self._menu, start_id+102, 'spec',
                          "Rebin spectrum", 'rebin')
        self._item_method(self._menu, start_id+103, 'spec',
                          "Convolve with gaussian", 'gauss_convolve')
        self._item_method(self._menu, start_id+104, 'spec',
                          "Estimate resolution", 'resol_est')
        self._item_method(self._menu, start_id+105, 'spec',
                          "Estimate SNR", 'snr_est')
        submenu = wx.Menu()
        self._item_method(submenu,start_id+110, 'spec', "Compute CCF", 'flux_ccf')
        self._menu.AppendSubMenu(submenu, "Other general recipes")

        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+150, 'spec', "Clip flux",
                          'flux_clip')

        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+200, 'spec', "Find lines",
                          'lines_find')
        self._item_method(self._menu, start_id+201, 'spec',
                          "Continuum from nodes", 'nodes_cont')
        submenu = wx.Menu()
        self._item_method(submenu, start_id+210, 'spec', "Find peaks",
                          'peaks_find')
        self._item_method(submenu, start_id+211, 'lines', "Extract nodes",
                          'nodes_extract')
        self._item_method(submenu, start_id+212, 'lines', "Clean nodes",
                          'nodes_clean')
        self._item_method(submenu, start_id+213, 'nodes',
                          "Interpolate nodes", 'nodes_interp')
        self._menu.AppendSubMenu(submenu, "Other recipes for continuum")
        self._menu.AppendSeparator()

        #self._item_method(self._menu, start_id+301, 'lines',
        #                  "Add and fit a system", 'add_syst')
        self._item_method(self._menu, start_id+300, 'z0',
                          "New system", 'syst_new')
        #self._item_method(self._menu, start_id+302, 'cont',
        #                  "Add and fit systems from line list",
        #                  'add_syst_from_lines')
        self._item_method(self._menu, start_id+301, 'lines',
                          "New systems from lines",
                          'systs_new_from_lines')
        self._item_method(self._menu, start_id+302, 'lines',
                          "Find candidate systems", 'cands_find')
        self._item_method(self._menu, start_id+303, 'z0',
                          "Improve systems", 'systs_improve')
        self._item_method(self._menu, start_id+304, 'z0',
                          "Complete systems", 'systs_complete')
        self._item_method(self._menu, start_id+305, 'z0',
                          "Fit systems", 'systs_fit')
        submenu = wx.Menu()
        #self._item_method(submenu, start_id+310, 'systs',
        #                  "Fit system", 'syst_fit')
        self._item_method(submenu, start_id+311, 'z0',
                          "Recreate models", 'mods_recreate')
        self._item_method(submenu, start_id+312, 'z0',
                          "Estimate SNR of systems", 'systs_snr')
        self._item_method(submenu, start_id+313, 'z0', "Update lines",
                          'lines_update')
        self._item_method(submenu, start_id+314, 'z0',
                          "Estimate position uncertainty", 'systs_sigmav')
        submenu.AppendSeparator()
        self._item_method(submenu, start_id+322, 'z0',
                          "Clean systems", 'systs_clean')
        self._item_method(submenu, start_id+323, 'z0',
                          "Extract systems", 'comp_extract')
        self._item_method(submenu, start_id+324, 'systs',
                          "Select systems", 'systs_select')
        submenu.AppendSeparator()
        self._item_method(submenu,start_id+331, 'z0', "Compute CCF",
                          'mods_ccf_max')
        self._menu.AppendSubMenu(submenu, "Other recipes for absorbers")
        #self._item_method(self._menu, start_id+303, 'systs',
        #                  "Add and fit systems from residuals",
        #                  'add_syst_from_resids')
        #self._item_method(self._menu, start_id+303, 'systs',
        #                  "New systems from residuals",
        #                  'systs_new_from_resids')
        #self._item_method(self._menu, start_id+304, 'systs',
        #                  "New systems from sliding technique",
        #                  'systs_new_from_slide')
        submenu = wx.Menu()

        """
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+401, 'lines',
                          "Simulate a system", 'syst_simul')
        self._item_method(self._menu, start_id+402, 'systs',
                          "Estimate completeness with simulated systems",
                          'systs_compl')
        """

class GUIMenuCourses(GUIMenu):
    def __init__(self,
                 gui,
                 start_id=5000,
                 **kwargs):
        super(GUIMenuCourses, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to Courses menu here
        #self._item_method(self._menu, start_id, 'spec', "Find lines",
        #                  ['gauss_convolve', 'peaks_find'])
        self._item_method(self._menu, start_id, 'spec', "Guess continuum",
                          ['lines_find', 'nodes_cont'])
        self._item_method(self._menu, start_id+1, 'spec', "Fit Ly-a forest",
                          ['lines_find', 'nodes_cont', 'systs_new_from_lines'])
        #self._item_method(self._menu, start_id+2, 'lines', "Fit systems",
        #                  ['add_syst_from_lines', 'add_syst_from_resids',
        #                   'add_syst_slide', 'compl_syst'])
        self._menu.AppendSeparator()
        #self._item_method(self._menu, start_id+101, 'spec', "From file",
        #                  'from_file', obj=self._gui._panel_sess)
        self._item(self._menu, start_id+101, None, "From JSON...\tCtrl+J",
                   lambda e: \
                   self._on_open(e, wildcard="JSON file (*.json)|*.json",
                                 action='json_load'))

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
        self._start_id = start_id

        # Add items to View menu here
        tab_id = [start_id+1, start_id+2, start_id+3, start_id+4, start_id+5]
        self._gui._menu_tab_id = tab_id
        self._item(self._menu, tab_id[0], 'spec', "Spectrum table",
                   lambda e: self._on_tab(e, 'spec'), key='spec')
        self._item(self._menu, tab_id[1], 'lines', "Line table",
                   lambda e: self._on_tab(e, 'lines'), key='lines')
        self._item(self._menu, tab_id[2], 'systs', "System table",
                   lambda e: self._on_tab(e, 'systs'), key='systs')
        self._item_graph(self._menu, tab_id[3], 'spec', "Metadata",
                         dlg_mini='meta', alt_title="Metadata")
        self._item_graph(self._menu, tab_id[4], 'spec', "Session log",
                         dlg_mini='log', alt_title="Session log")
        self._menu.AppendSeparator()
        """
        self._item(self._menu, start_id+101, 'systs',
                   "System detection correctness",
                   lambda e: self._on_ima(e, 'corr'))
        self._item(self._menu, start_id+102, 'systs',
                   "System detection completeness",
                   lambda e: self._on_ima(e, 'compl'))
        """
        self._item(self._menu, start_id+101, 'systs', "Compress system table",
                   self._on_compress)
        self._menu.AppendSeparator()
        self._item(self._menu, start_id+201, 'spec',
                   "Toggle log x axis", self._on_logx)
        self._item(self._menu, start_id+202, 'spec',
                   "Toggle log y axis", self._on_logy)
        self._norm = self._item(self._menu, start_id+203, 'spec', "Toggle normalization",
                                self._on_norm, key='norm')
        self._menu.AppendSeparator()
        self._submenu = wx.Menu()
        self._item_graph(self._menu, start_id+402, 'spec', "Edit graph elements",
                         dlg_mini='graph', alt_title="Graph elements")
        self._item_graph(self._submenu, start_id+314, 'spec', "Saturated H2O regions",
                         'spec_h2o_reg')
        self._item_graph(self._submenu, start_id+313, 'spec', "Redshift cursor",
                         'cursor_z_series', dlg_mini='systems',
                         targ=GraphCursorZSeries, alt_title="System controls")
        self._legend = self._item(self._submenu, start_id+315, 'spec', "Legend",
                                  self._on_legend, key='legend')
        self._menu.AppendSubMenu(self._submenu, "Toggle graph add-ons")

        #self._item_method(self._menu, start_id+401, 'spec',
        #                  "Edit graph details", '_sel_graph_cols', obj=self)


    def _on_compress(self, event):
        if self._menu.GetLabel(self._start_id+101) == "Compress system table":
            self._menu.SetLabel(self._start_id+101, "Uncompress system table")
        else:
            self._menu.SetLabel(self._start_id+101, "Compress system table")

        self._gui._sess_sel.systs._compress()
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
