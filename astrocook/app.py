from . import * #Cont, IO, Line, Plot, Recipe, Spec1DReader, System
from .utils import *
from .recipe import *
from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import Column, vstack
from collections import OrderedDict as od
from copy import deepcopy as dc
from datetime import datetime
import inspect
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg
from matplotlib.figure import Figure
import numpy as np
import os
import random
import re
import sys
import wx
import wx.grid as gridlib
import wx.lib.mixins.listctrl as listmix
import yaml

proc_descr = {'convolve': "Convolve with a custom profile",
              'extract_forest': "Extract forest",
              'extract_reg': "Extract spectral region",
              'find_extrema': "Find flux extrema (maxima and minima)",
              'select_extrema': "Select the most prominent extrema",
}


class MainFrame(wx.Frame):
    def __init__(self, parent=None, title="Astrocook", **kwargs):
        """ Constructor for the Frame class """ 

        size = (wx.DisplaySize()[0]*0.708, wx.DisplaySize()[1]*0.9)
        self.pad = 10
        super(MainFrame, self).__init__(parent, title=title, size=size)
        self.targ = None
        self.spec = None
        self.line = None
        self.cont = None
        self.syst = None

        self.rec = None
        self.plot = None
        
        self.targ_list = []
        self.spec_dict = {}
        self.z_dict = {}
        self.part_dict = {}
        self.line_dict = {}
        self.cont_dict = {}
        self.syst_dict = {}
        self.proc_dict = {}
        self.rec_dict = {}

        self.count = 0
        
        self.init_UI(**kwargs)
        self.IO = IO()

        try:
            bck = open("astrocook_app.bck", "r")
            lines = bck.readlines()
            self.path_chosen = lines[0][:-1]
            targ = lines[1]
            self.on_file_open(None, targ=targ, **kwargs)
        except:
            pass

        
    def dialog_proc(self, obj, proc):
        """ Run a procedure through a dialog window """

        self.obj = obj
        self.procs = [proc]
        self.descr = proc_descr[proc]
        dialog = ParamDialog(self)
        dialog.ShowModal()
        dialog.Destroy()
        self.params = dialog.params
        self.proc_dict[proc] = dialog.execute
        if dialog.execute == True:
            out = getattr(obj, proc)(**self.params)
            self.log[proc] = self.params
        else:
            out = None
        return out

    def dialog_rec(self, obj, rec):
        """ Run a recipe through a dialog window """

        self.obj = obj
        self.rec = Recipe(obj, rec)
        self.procs = self.rec.procs
        self.descr = self.rec.descr
        dialog = ParamDialog(self)
        dialog.ShowModal()
        dialog.Destroy()
        self.params = dialog.params
        self.rec_dict[rec] = dialog.execute
        if dialog.execute == True:
            #out = getattr(obj, rec)(**self.params)
            out = self.rec.line_find(**self.params)
            self.log[rec] = self.params
        else:
            out = None
        return out
    
    def init_line(self, panel):
        """ Create the Lines panel """

        self.line_gr = gridlib.Grid(panel)
        self.line_gr.CreateGrid(0, 5)
        self.line_gr.SetColLabelValue(0, "X")
        self.line_gr.SetColLabelValue(1, "XMIN")
        self.line_gr.SetColLabelValue(2, "XMAX")
        self.line_gr.SetColLabelValue(3, "Y")
        self.line_gr.SetColLabelValue(4, "DY")
        self.line_gr.Bind(gridlib.EVT_GRID_RANGE_SELECT, self.on_line_select)
        self.line_gr.Bind(gridlib.EVT_GRID_CELL_CHANGED, self.on_line_edit)
       
    def init_plot(self, panel):
        """ Create the Plot panel """
        
        self.fig = Figure()#figsize=(20,20))
        self.ax = self.fig.add_subplot(111)
        self.fig.tight_layout(rect=[-0.03, 0.02, 1.03, 1])
        self.plot = Plot(self.ax)
        self.plot_fig = FigureCanvasWxAgg(panel, -1, self.fig)
        self.plot_tb = NavigationToolbar2WxAgg(self.plot_fig)
        self.plot_tb.Realize()
        self.plot_pb = wx.Button(panel, label="Plot", size=(100,38))
        self.plot_cb = wx.Button(panel, label="Clear", size=(100,38))
        #self.plot_pb.Bind(wx.EVT_BUTTON,
        #                  lambda e: self.on_plot_draw(e, self.spec))
        self.plot_pb.Bind(wx.EVT_BUTTON, self.on_plot_draw)
        self.plot_cb.Bind(wx.EVT_BUTTON, self.on_plot_clear)
        
    def init_spec(self, panel):
        """ Create the Spectra panel """

        self.spec_lc = EditableListCtrl(panel, -1, style=wx.LC_REPORT)
        self.spec_lc.Bind(wx.EVT_LIST_BEGIN_LABEL_EDIT, self.on_spec_begin_edit)
        self.spec_lc.Bind(wx.EVT_LIST_END_LABEL_EDIT, self.on_spec_end_edit)
        self.spec_lc.Bind(wx.EVT_LIST_ITEM_SELECTED, self.on_spec_select)
        self.spec_lc.InsertColumn(0, 'target', width=150)
        self.spec_lc.InsertColumn(1, 'object', width=150)
        self.spec_lc.InsertColumn(2, 'redshift', width=150)
        self.spec_lc.InsertColumn(3, 'active range [nm]', width=150)
        self.spec_lc.InsertColumn(4, '# lines', width=150)
        self.spec_lc.InsertColumn(5, '# systems', width=150)
        
    def init_syst(self, panel):
        """ Create the Systems panel """

        self.syst_gr = gridlib.Grid(panel)
        self.syst_gr.CreateGrid(0, 9)
        self.syst_gr.SetColLabelValue(0, "SERIES")
        self.syst_gr.SetColLabelValue(1, "Z")
        self.syst_gr.SetColLabelValue(2, "N")
        self.syst_gr.SetColLabelValue(3, "B")
        self.syst_gr.SetColLabelValue(4, "BTUR")
        self.syst_gr.SetColLabelValue(5, "DZ")
        self.syst_gr.SetColLabelValue(6, "DN")
        self.syst_gr.SetColLabelValue(7, "DB")
        self.syst_gr.SetColLabelValue(8, "DBTUR")
        self.syst_gr.Bind(gridlib.EVT_GRID_RANGE_SELECT, self.on_syst_select)
        self.syst_gr.Bind(gridlib.EVT_GRID_CELL_CHANGED, self.on_syst_edit)
        #self.syst_gr.Bind(gridlib.EVT_BUTTON, self.on_syst_menu)
        
    def init_UI(self, **kwargs):
        """ Initialize the main frame """
        
        self.menu(**kwargs)

        panel = wx.Panel(self)
        self.init_spec(panel)
        self.init_line(panel)
        self.init_syst(panel)
        self.init_plot(panel)

        self.spec_lc.SetMaxSize((3000,120))
        self.line_gr.SetMaxSize((502,3000))
        self.syst_gr.SetMinSize((822,3000))

        
        # Spectra panel
        box_spec = wx.BoxSizer(wx.VERTICAL)
        box_spec.Add(wx.StaticText(panel, label="Spectra"))
        box_spec.Add(self.spec_lc, 1, wx.EXPAND)

        # Lines panel
        box_line = wx.BoxSizer(wx.VERTICAL)
        box_line.Add(wx.StaticText(panel, label="Lines"))
        box_line.Add(self.line_gr, 1, wx.EXPAND|wx.RIGHT, self.pad)

        # Systems panel
        box_syst = wx.BoxSizer(wx.VERTICAL)
        box_syst.Add(wx.StaticText(panel, label="Systems"))
        box_syst.Add(self.syst_gr, 1, wx.EXPAND)

        # Plot control panel
        box_ctrl = wx.BoxSizer(wx.HORIZONTAL)
        box_ctrl.Add(self.plot_tb, 1, wx.RIGHT, border=5)
        box_ctrl.Add(self.plot_pb, 0, wx.ALIGN_RIGHT|wx.RIGHT, border=5)
        box_ctrl.Add(self.plot_cb, 0, wx.ALIGN_RIGHT|wx.RIGHT, border=5)

        # Plot panel (includinc controls)
        box_plot = wx.BoxSizer(wx.VERTICAL)
        box_plot.Add(self.plot_fig, 1, wx.EXPAND)
        box_plot.Add(box_ctrl, 0, wx.TOP, self.pad)

        # Tables panel (Lines + Systems)
        box_table = wx.BoxSizer(wx.HORIZONTAL)
        box_table.Add(box_line, 1, wx.EXPAND)
        box_table.Add(box_syst, 1, wx.EXPAND|wx.ALIGN_LEFT)

        # Lists panel (Spectra + Tables)
        box_list = wx.BoxSizer(wx.VERTICAL)
        box_list.Add(box_spec, 0, wx.EXPAND|wx.BOTTOM, self.pad)
        box_list.Add(box_table, 1, wx.EXPAND)

        # Main panel (Lists + Plot)
        box_main = wx.GridSizer(2, 1, 0, 0)
        box_main.Add(box_list, 1, wx.EXPAND|wx.ALL, self.pad)
        box_main.Add(box_plot, 1, wx.EXPAND|wx.BOTTOM|wx.LEFT|wx.RIGHT,
                     self.pad)
        panel.SetSizer(box_main)
        
        self.Centre()
        self.Show()

    def menu(self, **kwargs):
        """ Create a menu in the frame """

        # Menu item IDs
        self.id_spec = 100
        self.id_line = 200
        self.id_cont = 300
        self.id_syst = 400
        self.id_syst_sel = 500
        
        # File menu
        self.file_menu = wx.Menu()
        
        file_open = wx.MenuItem(self.file_menu, wx.ID_OPEN, "&Open\tCtrl+O")
        file_save = wx.MenuItem(self.file_menu, self.id_spec, "&Save\tCtrl+S")
        file_quit = wx.MenuItem(self.file_menu, wx.ID_EXIT, "&Quit\tCtrl+Q")
        self.Bind(wx.EVT_MENU, lambda e: self.on_file_open(e, **kwargs),
                  file_open)
        self.Bind(wx.EVT_MENU, lambda e: self.on_file_save(e, **kwargs),
                  file_save)
        self.Bind(wx.EVT_MENU, self.on_quit, file_quit)

        self.file_menu.Append(file_open)
        self.file_menu.AppendSeparator()
        self.file_menu.Append(file_save)
        self.file_menu.AppendSeparator()
        self.file_menu.Append(file_quit)

        if (hasattr(self, 'spec') == False):
            self.menu_disable(self.file_menu, self.id_spec)

        # Procedures menu
        # N.B. Procedures are simple methods applied to objects (spectra, lines,
        # ecc.)
        self.proc_menu = wx.Menu()
        proc_spec_convolve = wx.MenuItem(
            self.proc_menu, 11, proc_descr['convolve']+"...")
        proc_spec_extract_forest = wx.MenuItem(
            self.proc_menu, 12, proc_descr['extract_forest']+"...")
        proc_spec_extract_reg = wx.MenuItem(
            self.proc_menu, 13, proc_descr['extract_reg']+"...")
        proc_spec_find_extrema = wx.MenuItem(
            self.proc_menu, 14, proc_descr['find_extrema'])
        proc_spec_select_extrema = wx.MenuItem(
            self.proc_menu, 15, proc_descr['select_extrema']+"...")
        
        self.Bind(wx.EVT_MENU, self.on_proc_spec_convolve,
                  proc_spec_convolve)
        self.Bind(wx.EVT_MENU, self.on_proc_spec_extract_forest,
                  proc_spec_extract_forest)
        self.Bind(wx.EVT_MENU, self.on_proc_spec_extract_reg,
                  proc_spec_extract_reg)
        self.Bind(wx.EVT_MENU, self.on_proc_spec_find_extrema,
                  proc_spec_find_extrema)
        self.Bind(wx.EVT_MENU, self.on_proc_spec_select_extrema,
                  proc_spec_select_extrema)
        
        self.proc_menu.Append(proc_spec_convolve)
        self.proc_menu.Append(proc_spec_extract_forest)
        self.proc_menu.Append(proc_spec_extract_reg)
        self.proc_menu.Append(proc_spec_find_extrema)
        self.proc_menu.Append(proc_spec_select_extrema)
        
        # Recipes menu
        self.rec_menu = wx.Menu()

        rec_line_find = wx.MenuItem(self.rec_menu, 21,
                                    rec_descr['line_find']+"...")

        self.Bind(wx.EVT_MENU, self.on_rec_line_find, rec_line_find)
        
        self.rec_menu.Append(rec_line_find)
        
        """
        rec_spec_extract = wx.MenuItem(self.rec_menu, self.id_spec+1,
                                       "E&xtract Spectral Region...")
        rec_cont_line_rem = wx.MenuItem(self.rec_menu, self.id_line,
                                        rec_descr['cont_line_rem']+"...")
        rec_cont_max_smooth = wx.MenuItem(self.rec_menu, self.id_line+1,
                                    "Find Continuum by Smoothing the "
                                    "Flux Maxima...")
        rec_syst_find = wx.MenuItem(self.rec_menu, self.id_cont,
                                    "Find &Systems...")
        rec_syst_def = wx.MenuItem(self.rec_menu, self.id_cont+1,
                                    "&Define System...")
        rec_syst_fit = wx.MenuItem(self.rec_menu, self.id_syst_sel,
                                   "&Fit Selected System...")
        #self.Bind(wx.EVT_MENU, self.on_spec_extract, rec_spec_extract)
        self.Bind(wx.EVT_MENU, self.on_line_find, rec_line_find)
        self.Bind(wx.EVT_MENU, self.on_cont_line_rem, rec_cont_line_rem)
        self.Bind(wx.EVT_MENU, self.on_cont_max_smooth, rec_cont_max_smooth)
        self.Bind(wx.EVT_MENU, self.on_syst_find, rec_syst_find)
        self.Bind(wx.EVT_MENU, self.on_syst_def, rec_syst_def)
        self.Bind(wx.EVT_MENU, self.on_syst_fit, rec_syst_fit)
        
        self.rec_menu.Append(rec_spec_extract)
        self.rec_menu.AppendSeparator()
        self.rec_menu.Append(rec_line_find)
        self.rec_menu.AppendSeparator()
        self.rec_menu.Append(rec_cont_line_rem)
        self.rec_menu.Append(rec_cont_max_smooth)
        self.rec_menu.AppendSeparator()
        self.rec_menu.Append(rec_syst_find)
        self.rec_menu.Append(rec_syst_def)
        self.rec_menu.Append(rec_syst_fit)
        """

        # Utilities menu
        self.util_menu = wx.Menu()

        util_log_view = wx.MenuItem(self.proc_menu, 41, "View log")

        self.Bind(wx.EVT_MENU, self.on_util_log_view, util_log_view)

        self.util_menu.Append(util_log_view)
        
        
        # Menu bar
        self.update_menu()
        menu_bar = wx.MenuBar()
        menu_bar.Append(self.file_menu, '&File')
        menu_bar.Append(self.proc_menu, '&Procedures')
        menu_bar.Append(self.rec_menu, '&Recipes')
        menu_bar.Append(self.util_menu, '&Utilities')
        self.SetMenuBar(menu_bar)        

        
    def menu_disable(self, menu, id):
        for i in range(10):
            try:
                menu.Enable(id+i, False)
            except:
                pass

    def menu_enable(self, menu, id):
        for i in range(100):
            try:
                menu.Enable(id+i, True)
            except:
                pass

    def on_proc_spec_convolve(self, event):
        proc = 'convolve'
        out = self.dialog_proc(self.spec, proc)
        if self.proc_dict[proc]:
            self.spec = out
            self.targ = self.targ + '_' + proc
            self.row = self.spec_lc.GetItemCount()
            self.spec_lc.insert_string_item(self.row, self.targ)
            self.spec_dict[self.targ] = self.spec
            self.update_all()

    def on_proc_spec_extract_forest(self, event):
        proc = 'extract_forest'
        out = self.dialog_proc(self.spec, proc)
        if self.proc_dict[proc]:
            self.spec = out
            self.targ = self.targ + '_' + self.params['ion']
            self.row = self.spec_lc.GetItemCount()
            self.spec_lc.insert_string_item(self.row, self.targ)
            self.spec_dict[self.targ] = self.spec
            self.update_all()
        
    def on_proc_spec_extract_reg(self, event):
        proc = 'extract_reg'
        out = self.dialog_proc(self.spec, proc)
        if self.proc_dict[proc]:
            self.spec = out
            self.targ = self.targ + '_%3.0f-%3.0f' \
                        % (self.params['xmin'], self.params['xmax'])
            self.row = self.spec_lc.GetItemCount()
            self.spec_lc.insert_string_item(self.row, self.targ)
            self.spec_dict[self.targ] = self.spec
            self.update_all()

    def on_proc_spec_find_extrema(self, event):
        proc = 'find_extrema'
        getattr(self.spec, proc)()
        self.plot.line(self.spec._mins)
        self.plot.line(self.spec._maxs, c='b')
        self.plot_fig.draw()
        self.log[proc] = None

    def on_proc_spec_select_extrema(self, event):
        proc = 'select_extrema'
        self.dialog_proc(self.spec, proc)
        if self.proc_dict[proc]:
            self.plot.line(self.spec._exts_sel, c='g', marker='o')
            self.plot_fig.draw()
            
    def on_rec_line_find(self, event):
        rec = 'line_find'
        run = self.dialog_rec(self.spec, rec)
        self.line = self.rec.line
        if self.rec_dict[rec]:            
            self.line_num = len(self.line.t)
            self.plot.line(self.line.t, c='g')
            self.plot_fig.draw()

    def on_util_log_view(self, event):
        dialog = wx.MessageDialog(None, yaml.safe_dump(self.log), 'Log', wx.OK)
        dialog.ShowModal()
        
    def on_cont_line_rem(self, event):
        """ Behaviour for Recipes > Find Continuum by Removing Lines """
        
        self.cont = Cont(self.spec, self.line)
        run = self.recipe_dialog(self.cont, "cont_line_rem")
        if run:
            self.cont_dict[self.targ] = self.cont
            self.update_plot()
            self.menu_enable(self.rec_menu, self.id_cont)

    def on_cont_max_smooth(self, event):
        """ Behaviour for Recipes > Find Continuum by Smoothing the Flux 
        Maxima """

        self.cont = Cont(self.spec, self.line)
        run = self.recipe_dialog(self.cont, 'cont_max_smooth')
        if run:
            self.cont_dict[self.targ] = self.cont
            self.update_plot()
            self.menu_enable(self.rec_menu, self.id_cont)
            
    def on_file_open(self, event, path='.', targ=None):
        """ Behaviour for File > Open """


        # otherwise ask the user what new file to open
        if (targ == None):
            wildcard = "Astrocook sessions (*.acs)|*.acs|" \
                       "FITS files (*.fits)|*.fits"
            with wx.FileDialog(self, "Open file", path,
                               wildcard=wildcard,
                               style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) \
                               as fileDialog:

                if fileDialog.ShowModal() == wx.ID_CANCEL:
                    return
                name = fileDialog.GetPath()
                self.path_chosen = fileDialog.GetDirectory()
                if (name[-4:] == '.acs'):
                    self.targ = fileDialog.GetFilename()
                    r = re.compile('_.{4}-.{2}-.{2}_.{2}-.{2}-.{2}.acs')
                    tail = self.targ[-24:]
                    if r.match(tail):
                        self.targ = self.targ[:-24]
                    else:
                        self.targ = self.targ[:-4]
                    try:
                        acs = self.IO.acs_read(name, self.path_chosen)
                        self.spec = acs.spec
                        self.spec_name = acs.spec_name
                    except IOError:
                        wx.LogError("Cannot open archive '%s'." % name)
                else:
                    self.targ = fileDialog.GetFilename()[:-5] 
                    try:
                        self.spec = self.IO.spec_read(name)
                        self.spec_name = name
                    except IOError:
                        wx.LogError("Cannot open file '%s'." % name)
        else:
            self.targ = targ
            name = self.path_chosen+'/'+targ+'.acs'
            acs = self.IO.acs_read(name, path=self.path_chosen)
            self.spec = acs.spec
            self.spec_name = acs.spec_name
            

        if self.spec != None:
            if self.targ in self.targ_list:
                self.count = self.count + 1
                self.targ = self.targ + '_%i' % self.count
            self.targ_list.append(self.targ)
                    
            self.spec_dict[self.targ] = self.spec
            self.row = self.spec_lc.GetItemCount()
            self.spec_lc.insert_string_item(self.row, self.targ)
            self.menu_enable(self.file_menu, self.id_spec)
            self.menu_enable(self.rec_menu, self.id_spec)
            try:
                self.line = acs.line
                self.line_name = acs.line_name
                self.line_dict[self.targ] = self.line
                self.update_line()
                self.menu_enable(self.file_menu, self.id_line)
                self.menu_enable(self.rec_menu, self.id_line)
            except:
                self.line = None
            
            try:
                self.cont = acs.cont
                self.cont_name = acs.cont_name
                self.cont_dict[self.targ] = self.cont
                self.menu_enable(self.file_menu, self.id_cont)
                self.menu_enable(self.rec_menu, self.id_cont)
            except:
                self.cont = None
                
            try:
                self.syst = acs.syst
                self.syst_name = acs.syst_name
                self.syst_dict[self.targ] = self.syst
                self.update_syst()
                self.menu_enable(self.file_menu, self.id_syst)
                self.menu_enable(self.rec_menu, self.id_syst)
            except:
                self.syst = None

            try:
                self.log = acs.log
                self.log_name = acs.log_name
            except:
                self.log = dict()
                self.log["targ"] = self.targ
                
            self.update_syst()
            self.update_spec()
            self.update_plot()

    def on_file_save(self, event, path='.'):
        """ Behaviour for File > Save """

        #timestamp = \
        #    '_'+str(datetime.now()).replace(" ", "_").replace(":", "-")[:-7]
        snapshot = self.targ #+ timestamp
        root = path + snapshot
        with wx.FileDialog(self, "Save session", path, snapshot,
                           wildcard="Astrocook session (*.acs)|*.acs",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return

            name = fileDialog.GetPath()
            self.path_chosen = fileDialog.GetDirectory()
            try:
                acs = self
                self.IO.acs_write(acs, name, self.path_chosen)

            except IOError:
                wx.LogError("Cannot save session '%s'." % name)

    def on_line_edit(self, event):
        row = event.GetRow()
        col = event.GetCol()
        label = self.line_gr.GetColLabelValue(col)
        data = self.line_gr.GetCellValue(row, col)
        self.line._t[label][row] = data
        
    def on_line_find(self, event):
        """ Behaviour for Recipes > Find Lines """

        self.line = Line(self.spec)
        run = self.method_dialog(self.line, "find_special")
        if run:            
            self.line_num = len(self.line.t)
            self.update_spec()
            self.update_line()
            self.update_plot()
            self.menu_enable(self.rec_menu, self.id_line)
        
    def on_line_select(self, event):
        """ Behaviour for line selection """
        
        if event.GetTopRow() == event.GetBottomRow():            
            sel = event.GetTopRow()
            try:
                self.line_focus.remove()
            except:
                pass

            # Move to line.py
            x = self.line.x[sel]
            y = self.line.y[sel]
            xmin = self.line.xmin[sel]
            xmax = self.line.xmax[sel]
            #self.line_focus, = self.ax.plot(x, y, c='C0', marker='o', ms=20,
            #                                alpha=0.2)
            self.update_plot(xmin.value-0.3, xmax.value+0.3)
            self.ax.axvline(x=x.value, color='C3', alpha=0.5)
            self.ax.axvline(x=xmin.value, color='C3', alpha=0.5, linestyle=':')
            self.ax.axvline(x=xmax.value, color='C3', alpha=0.5, linestyle=':')
            self.plot_fig.draw()
            
    def on_plot_clear(self, event):
        self.ax.clear()
        self.plot_fig.draw()

    #def on_plot_draw(self, event, obj):
    def on_plot_draw(self, event, xmin=None, xmax=None):
        if xmin == None:
            xmin = self.xmin
        if xmax == None:
            xmax = self.xmax
        """
        self.spec = self.spec_dict[self.targ]
        self.spec.plot(ax=self.ax)
        try:
            self.line.plot_new(ax=self.ax)
        except:
            pass
        try:
            self.ax.plot(self.cont.t['X'], self.cont.t['Y'])
        except:
            pass
        try:
            self.syst = self.syst_dict[self.targ]
        except:
            pass
	"""
        self.plot.spec(self.spec.t)
        self.ax.set_xlim(xmin, xmax)
        self.plot_fig.draw()
        
    def on_quit(self, event):
        """ Behaviour for File > Quit """
        
        bck = open("astrocook_app.bck", "wb")
        bck.write(self.path_chosen+'\n')
        bck.write(self.targ)
        bck.close()
        self.Close()

    def on_spec_begin_edit(self, event):
        """ Veto the editing of some columns of the spectrum list """
        if event.GetColumn() in [0,3,4,5]:
            event.Veto()
        else:
            event.Skip()

    def on_spec_end_edit(self, event):
        """ Behaviour when spectrum is edited on list """
        
        index = self.spec_lc.GetFocusedItem()
        row = event.GetIndex()
        col = event.GetColumn()
        data = event.GetLabel()
        self.spec_lc.SetItem(row, col, data)
        try:
            self.z_dict[self.targ] = self.spec_lc.GetItem(index, 2).GetText()
        except:
            pass
        
    def on_spec_extract(self, event):
        """ Behavior for Recipes > Extract Spectral Region """

        #"""
        try:
            z = self.z_dict[self.targ]
        except:
            z = 0.0
        #"""
        self.rec = Recipe(self.spec, 'spec_extract')
        dialog = ParamDialog(self, title="Extract Spectral Region")
        dialog.ShowModal()
        dialog.Destroy()
        if dialog.execute == True:
            if self.rec.params['forest'] == True:
                params = {k: self.rec.params[k] for k in ['ion', 'zem',
                                                        'prox_vel']}
                forest = self.spec.extract_forest(**params)
                self.targ = self.targ + '_' + params['ion']
                self.z_dict[self.targ] = params['zem']
            else:
                params = {k: self.rec.params[k] for k in ['xmin', 'xmax']}
                self.targ = self.targ + '_%3.0f-%3.0f' \
                            % (params['xmin'], params['xmax'])
                if float(params['zem']) != 0.0:
                    self.z_dict[self.targ] = float(params['zem'])
                forest = self.spec.extract_reg(**params)

                
            self.row = self.spec_lc.GetItemCount()
            self.spec_lc.insert_string_item(self.row, self.targ)

            self.spec = forest
            self.spec_dict[self.targ] = self.spec
            self.update_all()
            self.log["spec_extract"] = dict(self.params)
            #self.update_spec()
            #self.update_line()
            #self.update_syst()
            #self.update_plot()
        #"""
        
    def on_spec_select(self, event):
        """ Behaviour when spectrum is selected from list """

        item = self.spec_lc.GetItem(self.spec_lc.GetFirstSelected(), 0)
        self.targ = item.GetText()
        self.row = event.GetIndex()
        self.spec = self.spec_dict[self.targ]
        self.update_all()

    def on_syst_def(self, event):
        """ Behaviour for Systems > Define System """
        self.params = od([('series', 'Ly'), ('z', '0.0'), ('N', '1e14'),
                          ('b', '20')])
        dialog = ParamDialog(self, title="Define System")
        dialog.ShowModal()
        dialog.Destroy()
        if dialog.execute == True:
            #vary = [bool(v) for v in self.params['vary']]
            syst = System(self.spec, self.line, self.cont,
                          series=self.params['series'], 
                          z=self.params['z'], N=self.params['N'],
                          b=self.params['b'])
            dx = 0.5
            syst.create_line(dx)

            if self.syst != None:
                self.syst.merge(syst)
            else:
                self.syst = syst
            self.syst._line._t = vstack([self.line._t, syst._line._t])

            self.line = self.syst._line#vstack([self.line._t, syst._line._t])
            self.line_dict[self.targ] = self.line
            self.syst_dict[self.targ] = self.syst
            self.log["syst_def"] = dict(self.params)
            self.update_spec()
            self.update_line()
            self.update_syst()

    def on_syst_edit(self, event):
        row = event.GetRow()
        col = event.GetCol()
        label = self.syst_gr.GetColLabelValue(col)
        data = self.syst_gr.GetCellValue(row, col)
        self.syst._t[label][row] = data
        
    def on_syst_find(self, event):
        """ Behaviour for Recipes > Find Systems """

        syst = System(self.spec, self.line, self.cont)
        run = self.recipe_dialog(syst, "syst_find")
        if run:
            if self.syst != None:
                self.syst.merge(syst)
            else:
                self.syst = syst
            self.syst_dict[self.targ] = self.syst
            self.update_spec()
            self.update_syst()
            
    def on_syst_fit(self, event):
        self.syst = self.syst_dict[self.targ]
        
        dialog = SystDialog(self, title="Fit selected system")
        dialog.ShowModal()
        dialog.Close()
        while dialog.execute == True:
            self.syst.fit(self.z_sel, norm=False)
            new_z = (np.abs(self.syst._t['Z']-self.z_sel)).argmin()
            self.z_sel = self.syst._t['Z'][new_z]
            self.ax.plot(self.syst._chunk['X'],
                         self.syst._chunk['MODEL'])
            self.plot_fig.draw()

            self.update_syst()
            dialog = SystDialog(self, title="Fit selected system")
            dialog.ShowModal()
            dialog.Close()
        dialog.Destroy()
        self.update_syst()
            
    def on_syst_select(self, event):
        """ Behaviour for system selection """        
        if event.GetTopRow() == event.GetBottomRow():  
            sel = event.GetTopRow()
            try:
                self.syst_focus.remove()
            except:
                pass
            z = self.syst.t['Z'][sel]
            x = [(1 + z) * dict_wave[i].value \
                 for i in dict_series[self.syst.t['SERIES'][sel]]]
            dx = 0.5
            h = np.max(self.spec.y)
            if self.ax.get_xlim() != (self.xmin, self.xmax):
                self.update_plot()
            self.syst.model(z, norm=False, dx=dx)
            self.syst_focus = self.ax.bar(x, h, dx, 0, color='C1', alpha=0.2)
            self.plot_fig.draw()

            """
            self.syst_sel, self.syst_idx = self.syst.extract(z)
            self.syst_sel._group = self.syst._group
            self.syst_sel._group_t = self.syst._group_t
            self.syst_sel._group_map = self.syst._group_map
            self.syst_sel._map = self.syst._map
            self.syst_sel.model(z, dx)            
            """
            self.menu_enable(self.rec_menu, self.id_syst_sel)
            self.z_sel = z

        
    def update_all(self):
        """ Update all panels """
        
        self.spec = self.spec_dict[self.targ]
        try:
            self.line = self.line_dict[self.targ]
        except:
            self.line = None
        try:
            self.cont = self.cont_dict[self.targ]
        except:
            self.cont = None

        try:
            self.syst = self.syst_dict[self.targ]
        except:
            self.syst = None

        self.update_spec()
        self.update_line()
        self.update_syst()
        self.update_plot()
        self.update_menu()
        
        
    def update_line(self):
        """ Update the line table """
        
        try:
            self.line_gr.DeleteRows(pos=0, numRows=self.line_gr.GetNumberRows())
        except:
            pass
        try:
            self.line_gr.AppendRows(len(self.line.t))
            for i, l in enumerate(self.line.t):
                self.line_gr.SetCellValue(i, 0, "%3.3f" % l['X'])
                self.line_gr.SetCellValue(i, 1, "%3.3f" % l['XMIN'])
                self.line_gr.SetCellValue(i, 2, "%3.3f" % l['XMAX'])
                self.line_gr.SetCellValue(i, 3, "%3.3f" % l['Y'])
                self.line_gr.SetCellValue(i, 4, "%3.3f" % l['DY'])
            self.line_dict[self.targ] = self.line
        except:
            pass

    def update_menu(self):
        """ Update the menus """
        
        #if (hasattr(self, 'spec') == False):
        if (self.spec is None):
            self.menu_disable(self.rec_menu, self.id_spec)
        else:
            self.menu_enable(self.rec_menu, self.id_spec)
        if (self.line is None):
            self.menu_disable(self.rec_menu, self.id_line)
        else:
            self.menu_enable(self.rec_menu, self.id_line)
        if (self.cont is None):
            self.menu_disable(self.rec_menu, self.id_cont)
        else:
            self.menu_enable(self.rec_menu, self.id_cont)
        if (self.syst is None):
            self.menu_disable(self.rec_menu, self.id_syst)            
        else:
            self.menu_enable(self.rec_menu, self.id_syst)

        self.menu_disable(self.rec_menu, self.id_syst_sel)
                        
    
            
    def update_plot(self, *kwargs):
        """ Update the plot panel """

        self.on_plot_clear(None)
        self.on_plot_draw(None, *kwargs)

    def update_spec(self):
        """ Update the spec list """

        self.spec = self.spec_dict[self.targ]

        try:
            self.spec_lc.SetItem(self.row, 2, str(self.z_dict[self.targ]))
        except:
            pass

        self.xmin = self.spec.t['X'][0]
        self.xmax = self.spec.t['X'][-1]
        self.spec_lc.SetItem(self.row, 3, "[%3.2f, %3.2f]"
                             % (self.xmin, self.xmax))

        try:
            self.spec_lc.SetItem(self.row, 4, str(len(self.line.t)))
        except:
            pass
            
        try:
            self.ax.plot(self.cont.t['X'], self.cont.t['Y'])
        except:
            pass
        try:
            self.spec_lc.SetItem(self.row, 5, str(len(self.syst.t)))
        except:
            pass
        
    def update_syst(self):
        """ Update the system table """

        try:
            self.syst_gr.DeleteRows(pos=0,
                                    numRows=self.syst_gr.GetNumberRows())
        except:
            pass

        try:
            self.syst_gr.AppendRows(len(self.syst.t))
            for i, s in enumerate(self.syst.t):
                self.syst_gr.SetCellValue(i, 0, str(s['SERIES']))
                self.syst_gr.SetCellValue(i, 1, "%3.5f" % s['Z'])
                self.syst_gr.SetCellValue(i, 2, "%3.3e" % s['N'])
                self.syst_gr.SetCellValue(i, 3, "%3.3f" % s['B'])
                self.syst_gr.SetCellValue(i, 4, "%3.3f" % s['BTUR'])
                self.syst_gr.SetCellValue(i, 5, "%3.5f" % s['DZ'])
                self.syst_gr.SetCellValue(i, 6, "%3.3e" % s['DN'])
                self.syst_gr.SetCellValue(i, 7, "%3.3f" % s['DB'])
                self.syst_gr.SetCellValue(i, 8, "%3.3f" % s['DBTUR'])
            self.syst_dict[self.targ] = self.syst
        except:
            pass


class EditableListCtrl(wx.ListCtrl, listmix.TextEditMixin):
    def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        """ Constructor for the EditableListCtrl class """
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.TextEditMixin.__init__(self)


    def insert_string_item(self, *args):
        self.InsertItem(*args)
        listmix.TextEditMixin.__init__(self)
        
class ParamDialog(wx.Dialog):

    def __init__(self, parent=None, size=(250,500), **kwargs):
        """ Constructor for the ParamDialog class """

        self.p = parent

        self.params = {}
        for p in self.p.procs:
            method = getattr(self.p.obj, p)
            try:
                self.params.update({k: v for (k,v) in
                                    zip(inspect.getargspec(method)[0][1:],
                                        inspect.getargspec(method)[3])})
            # When the procedure has no parameters
            except:
                pass  
        self.dialog = self.params
        """
        except:
            method = getattr(self.p.obj, self.p.proc)
            title = self.p.proc
            self.params = {k: v for (k,v) in
                           zip(inspect.getargspec(method)[0][1:],
                               inspect.getargspec(method)[3])}
            self.dialog = self.params
        """
        super(ParamDialog, self).__init__(parent=self.p, title=self.p.descr)
        self.init_UI()
        
    def init_UI(self):
        """ Initialize the main frame """
        
        panel = wx.Panel(self)
        box_main = wx.BoxSizer(wx.VERTICAL)

        box_params = wx.BoxSizer(wx.VERTICAL)
        self.par = []
        self.dial = []
        self.ctrl = []
        self.type = []
        #for d, v in self.p.rec.dialog.iteritems():
        #for d in self.dialog:
        for p in self.params:
            box_param = wx.BoxSizer(wx.HORIZONTAL)
            #p = self.dialog[d]
            v = self.params[p]
            self.par.append(p)
            self.dial.append(p)
            self.type.append(type(v))
            if type(v) == bool and 1==0:
                rb = wx.RadioButton(panel, -1, label=p)
                box_param.Add(rb, 1, 0)
                self.ctrl.append(rb)
            else:
                st = wx.StaticText(panel, -1, label=p)
                tc = wx.TextCtrl(panel, -1, value=str(v), size=(150,25))
                box_param.Add(st, 1, 0)
                box_param.Add(tc, 1, 0)
                self.ctrl.append(tc)


            box_params.Add(box_param, 1, 0, 0)
        panel.SetSizer(box_params)
        
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        cancel_button = wx.Button(self, label='Cancel')
        run_button = wx.Button(self, label='Run')
        run_button.SetDefault()
        buttons.Add(cancel_button, 0, wx.RIGHT, border=5)
        buttons.Add(run_button, 0)
 
        box_main.Add(panel, 0, wx.EXPAND|wx.ALL, border=10)
        box_main.Add(buttons, 0, wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT|wx.BOTTOM,
                     border=10)
        box_main.SetSizeHints(self)

        self.SetSizer(box_main)
        
        cancel_button.Bind(wx.EVT_BUTTON, self.on_cancel)
        run_button.Bind(wx.EVT_BUTTON, self.on_run)

        self.Centre()
        self.Show()
        
    def on_cancel(self, e):
        self.execute = False
        self.Close()

    def on_run(self, e):
        for p, t, ctrl in zip(self.par, self.type, self.ctrl):
            try:
                self.params[p] = t(ctrl.GetValue())
            except:
                self.params[p] = None
        self.execute = True
        self.Close()

                
class SystDialog(wx.Dialog):

    def __init__(self, parent=None, title="Parameters", **kwargs):
        """ Constructor for the ParamDialog class """

        size = (wx.DisplaySize()[0]*0.5, wx.DisplaySize()[1]*0.7)
        super(SystDialog, self).__init__(parent, title=title, size=size) 

        self.p = parent
        self.syst = self.p.syst#_sel
        self.group = self.p.syst._group
        self.z = self.p.z_sel
        cond = np.logical_and(self.syst._t['Z'] > self.z-0.002,
                                   self.syst._t['Z'] < self.z+0.002)
        cond = self.syst._t['Z'] == self.z
        self.sel = np.where(cond)[0]
        self.zs = self.syst._t['Z'][self.sel]
        self.ions = np.unique([dict_series[i] \
                                 for i in self.group['SERIES']])#[self.sel]])
        self.ions = np.unique(np.ravel(self.group['ION']))
        self.ions = self.ions[np.where(self.ions != 'unknown')]
        self.init_UI()

    def init_buttons(self, panel):
        self.syst_b = wx.Button(panel, label="Add system", size=(100,38))
        self.line_b = wx.Button(panel, label="Add line", size=(100,38))
        self.syst_b.Bind(wx.EVT_BUTTON, self.on_tab_add)
        self.line_b.Bind(wx.EVT_BUTTON, lambda e: self.on_tab_add(e, 'unknown'))

        
    def init_plot(self, panel):
        """ Create the spectrum panel """
        self.fig = Figure()#figsize=(20,20))


        
        #sel = np.where(self.p.syst._t['Z'] == self.p.z_sel)[0]
        
        #ions = self.p.syst._t['SERIES'][sel][0] #self.p.syst._group['ION']
        #waves = [dict_wave[i].value for i in ions]
        #ions = ions[np.argsort(waves)]
        rown = 5.
        self.pn = len(self.ions)
        row = int(min(self.pn,rown))
        col = int(np.ceil(self.pn/rown))
        #fig = plt.figure(figsize=(col*6, n*3.5))
        grid = gs(row,col)
        self.ax = []
        for p in range(self.pn):
            self.ax.append(self.fig.add_subplot(grid[int(p%rown),
                                                int(np.floor(p/rown))]))
        try:
            self.fig.suptitle(r"$\chi_r^2$ = %3.1f" % self._fit.redchi)
        except:
            pass

        #grid.update(wspace=0.2, hspace=0.0)
        #grid.tight_layout(self.fig, rect=[0.01, 0.01, 1, 0.9], h_pad=0.0)
        self.plot_fig = FigureCanvasWxAgg(panel, -1, self.fig)
        self.plot_tb = NavigationToolbar2WxAgg(self.plot_fig)
        self.plot_tb.Realize()
        
    def init_tab(self, panel):
        """ Create the system list panel """

        gr = gridlib.Grid(panel)
        gr.CreateGrid(0, 11)
        gr.SetColLabelValue(0, "X")
        gr.SetColLabelValue(1, "XMIN")
        gr.SetColLabelValue(2, "XMAX")
        gr.SetColLabelValue(3, "ION")
        gr.SetColLabelValue(4, "Z")
        gr.SetColLabelValue(5, "N")
        gr.SetColLabelValue(6, "B")
        gr.SetColLabelValue(7, "BTUR")
        gr.SetColLabelValue(8, "VARY")
        gr.SetColLabelValue(9, "EXPR")
        gr.SetColLabelValue(10, "#")
        gr.Bind(gridlib.EVT_GRID_CELL_CHANGED,
                lambda e: self.on_tab_edit(e, gr))
        #self.focus_gr.Bind(gridlib.EVT_GRID_RANGE_SELECT, self.on_line_select)

        return gr

    def init_UI(self):
        """ Initialize the main frame """

        panel = wx.Panel(self)

        #self.focus_gr = self.init_tab(panel)
        self.add_gr = self.init_tab(panel)
        self.init_buttons(panel)
        #self.add_gr.SetColLabelSize(0)
        self.init_plot(panel)
        #self.update_group()

        self.box_main = wx.BoxSizer(wx.VERTICAL)

        box_focus = wx.BoxSizer(wx.VERTICAL)
        #box_focus.Add(wx.StaticText(panel, label="Focus"))
        #box_focus.Add(self.focus_gr, 1, wx.BOTTOM, border=5)
        box_focus.Add(self.add_gr, 1, wx.BOTTOM, border=10)

        box_button = wx.BoxSizer(wx.HORIZONTAL)
        box_button.Add(self.syst_b, 0, wx.RIGHT, border=5)
        box_button.Add(self.line_b, 0, wx.RIGHT, border=5)

        box_ctrl = wx.BoxSizer(wx.HORIZONTAL)
        box_ctrl.Add(self.plot_tb, 0, wx.RIGHT)        

        box_add = wx.BoxSizer(wx.HORIZONTAL)
        #box_add.Add(wx.StaticText(panel, label="Additional lines"))
        #box_add.Add(self.add_gr, 1, wx.BOTTOM, border=10)
        box_add.Add(box_button)
        box_add.Add(box_ctrl)

        box_plot = wx.BoxSizer(wx.VERTICAL)
        box_plot.Add(self.plot_fig, 1, wx.EXPAND)

        box_disp = wx.BoxSizer(wx.VERTICAL)
        box_disp.Add(box_focus)
        box_disp.Add(box_add)
        box_disp.Add(box_plot, 1, wx.EXPAND|wx.TOP|wx.BOTTOM, border=10)
        #box_disp.Add(box_ctrl)


        panel.SetSizer(box_disp)
        
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        cancel_button = wx.Button(self, label='Cancel')
        run_button = wx.Button(self, label='Run')
        run_button.SetDefault()
        buttons.Add(cancel_button, 0, wx.RIGHT, border=5)
        buttons.Add(run_button, 0)
 
        self.box_main.Add(panel, 0, wx.EXPAND|wx.ALL, border=10)
        self.box_main.Add(
            buttons, 0, wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT|wx.BOTTOM, border=10)
        self.box_main.SetSizeHints(self)

        self.SetSizer(self.box_main)
        
        cancel_button.Bind(wx.EVT_BUTTON, self.on_cancel)
        run_button.Bind(wx.EVT_BUTTON, self.on_run)

        self.update_tab()
        self.update_plot()
        
        self.Centre()
        self.Show()

    def on_cancel(self, event):
        self.execute = False
        self.Close()

        
    def on_run(self, event):
        self.execute = True
        self.Close()

    def on_tab_add(self, event, series=None):


        dx = 0.5
        
        chunk = self.syst._chunk
        group = self.syst._group

        # Minimum X
        x_min = chunk['X'][(chunk['Y']-chunk['MODEL']).argmin()]
    
        # Minimum positions in system list and in group
        group_min = np.abs(x_min-self.group['X']).argmin()
        syst_min = np.abs(self.group['Z'][group_min]\
                          -self.syst._t['Z']).argmin()
        z_min = self.syst._t['Z'][syst_min]

        if series == None:
            ion_min = self.group['ION'][group_min]
            #ions = dict_series[self.group['SERIES'][group_min]]
            ions = np.array(self.group['ION'])
        elif series == 'unknown':
            ion_min = dict_series[series][0]
            ions = [ion_min]
        z_add = x_min/dict_wave[ion_min].value - 1     
        
        # Create a duplicate of the system
        dupl = dc(self.syst._t[syst_min])
        dupl['Z'] = z_add
        if series is not None:
            dupl['SERIES'] = series
        
        # Update line and map tables
        match_z = np.where(self.syst._map['Z'] == z_min)[0]
        x_add = [(1+z_add)*dict_wave[i].value for i in ions]
        if series is not None:
            if series is 'unknown':
                match_z = [match_z[0]]
        match_x = np.where(np.in1d(self.syst._line.t['X'],
                                   self.syst._map['X'][match_z]))[0]
        for i, r in enumerate(self.syst._map[match_z]):
            self.syst._map.add_row(r)
            self.syst._map['X'][-1] = x_add[i]
            self.syst._map['Z'][-1] = z_add
        for i, r in enumerate(self.syst._line.t[match_x]):
            self.syst._line._t.add_row(r)
            self.syst._line._t['X'][-1] = x_add[i]
        
        self.syst._map.sort('Z')
        self.syst._line._t.sort('X')

        self.syst._t = vstack([self.syst._t, dupl])#merge(dupl)
        self.syst.group(self.z, 0.5)
        self.syst.model(self.z, dx)
        self.update_tab()
        self.update_plot()
        
    def on_tab_edit(self, event, tab):
        """ Behaviour when table is edited """
        
        row = event.GetRow()
        col = event.GetCol()
        label = tab.GetColLabelValue(col)
        data = tab.GetCellValue(row, col)
        idx = int(tab.GetCellValue(row, 10))
        group = self.syst._group
        t = self.syst._t
        line = self.syst._line._t
        map = self.syst._map
        
        x = group['X'][idx]
        z = group['Z'][idx]
        dx = 0.5
        group_z = np.where(group['Z'] == z)[0]
        group[label][group_z] = data
        if label in t.colnames:
            t[label][np.where(t['Z'] == z)[0]] = data
        if label in line.colnames:
            line[label][np.where(line['X'] == x)[0]] = data
        if label == 'Z':
            map[label][np.where(map['Z'] == z)[0]] = data
        #"""
        if label == 'XMIN':
            group['XMIN'][group_z] = data
            group['ZMIN'][group_z] = [g['XMIN']/dict_wave[g['ION']].value - 1 \
                                      for g in group[group_z]]
        if label == 'XMAX':
            group['XMAX'][group_z] = data
            group['ZMAX'][group_z] = [g['XMAX']/dict_wave[g['ION']].value - 1 \
                                      for g in group[group_z]]
        #"""
        self.syst.model(self.z, dx)
        self.update_plot()
        
    def update_plot(self):
        for p in range(self.pn):
            self.ax[p].clear()
        self.syst.plot(z=self.z, ax=self.ax, ions=self.ions)
        self.plot_fig.draw()
        
        
    def update_tab(self):
        """ Update the system table """

        try:
            self.focus_gr.DeleteRows(pos=0,
                                     numRows=self.focus_gr.GetNumberRows())
        except:
            pass
        try:
            self.add_gr.DeleteRows(pos=0,
                                   numRows=self.add_gr.GetNumberRows())
        except:
            pass

        group = self.syst._group
        for i, g in enumerate(group):
            #if g['Z'] == self.z:
            #if g['Z'] in self.zs:
            #    tab = self.focus_gr
            #else:
            #    tab = self.add_gr
            tab = self.add_gr
            r = tab.GetNumberRows()
            tab.AppendRows(1)
            tab.SetCellValue(r, 0, "%3.3f" % g['X'])
            tab.SetCellValue(r, 1, "%3.3f" % g['XMIN'])
            tab.SetCellValue(r, 2, "%3.3f" % g['XMAX'])
            tab.SetCellValue(r, 3, str(g['ION']))
            tab.SetCellValue(r, 4, "%3.5f" % g['Z'])
            tab.SetCellValue(r, 5, "%3.3e" % g['N'])
            tab.SetCellValue(r, 6, "%3.3f" % g['B'])
            tab.SetCellValue(r, 7, "%3.3f" % g['BTUR'])
            tab.SetCellValue(r, 8, str(g['VARY']))
            tab.SetCellValue(r, 9, str(g['EXPR']))
            tab.SetCellValue(r, 10, str(i))
        
        self.box_main.Fit(self)
