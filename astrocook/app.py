from . import Cont, IO, Line, Spec1DReader, System
from .utils import *
from astropy.io import fits
from astropy.table import vstack
from collections import OrderedDict as od
from datetime import datetime
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg
from matplotlib.figure import Figure
import numpy as np
import os
import random
import sys
import wx
import wx.grid as gridlib
import wx.lib.mixins.listctrl as listmix

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
        
        self.targ_list = []
        self.spec_dict = {}
        self.z_dict = {}
        self.part_dict = {}
        self.line_dict = {}
        self.cont_dict = {}
        self.syst_dict = {} 

        self.count = 0
        
        self.init_UI(**kwargs)
        self.IO = IO()

        
    def init_line(self, panel):
        """ Create the line list panel """

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
        """ Create the spectrum panel """
        self.fig = Figure()#figsize=(20,20))
        self.ax = self.fig.add_subplot(111)
        self.fig.tight_layout(rect=[-0.03, 0.02, 1.03, 1])
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
        """ Create the spectrum panel """

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
        """ Create the system list panel """

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

        #box_main = wx.BoxSizer(wx.VERTICAL)
        box_main = wx.GridSizer(2, 1, 0, 0)
        
        box_list = wx.BoxSizer(wx.VERTICAL)
        #box_list = wx.GridSizer(3, 1, 0, 0)
        box_displ = wx.BoxSizer(wx.VERTICAL)

        box_spec = wx.BoxSizer(wx.VERTICAL)
        #box_table = wx.GridSizer(2, self.pad, self.pad)
        box_table = wx.BoxSizer(wx.HORIZONTAL)
        box_line = wx.BoxSizer(wx.VERTICAL)
        box_syst = wx.BoxSizer(wx.VERTICAL)
        
        box_plot = wx.BoxSizer(wx.VERTICAL)
        box_ctrl = wx.BoxSizer(wx.HORIZONTAL)

        box_spec.Add(wx.StaticText(panel, label="Spectra"))
        box_spec.Add(self.spec_lc, 1, wx.EXPAND)
        box_line.Add(wx.StaticText(panel, label="Lines"))
        #box_line.Add(self.line_lc, 1, wx.EXPAND)
        box_line.Add(self.line_gr, 1, wx.EXPAND|wx.RIGHT, self.pad)
        box_syst.Add(wx.StaticText(panel, label="Systems"))
        box_syst.Add(self.syst_gr, 1, wx.EXPAND)
        box_ctrl.Add(self.plot_tb, 1, wx.RIGHT, border=5)
        box_ctrl.Add(self.plot_pb, 0, wx.ALIGN_RIGHT|wx.RIGHT, border=5)
        box_ctrl.Add(self.plot_cb, 0, wx.ALIGN_RIGHT|wx.RIGHT, border=5)
        box_plot.Add(self.plot_fig, 1, wx.EXPAND)
        box_plot.Add(box_ctrl, 0, wx.TOP, self.pad)

        box_table.Add(box_line, 1, wx.EXPAND)
        box_table.Add(box_syst, 1, wx.EXPAND|wx.ALIGN_LEFT)
        
        box_list.Add(box_spec, 0, wx.EXPAND|wx.BOTTOM, self.pad)
        #box_list.Add(box_line, 1, wx.EXPAND|wx.RIGHT, self.pad)
        #box_list.Add(box_syst, 1, wx.EXPAND)
        box_list.Add(box_table, 1, wx.EXPAND)
        
        #box_displ.Add(box_plot, 1, wx.EXPAND)
        #box_displ.Add(box_ctrl, 1)

        box_main.Add(box_list, 1, wx.EXPAND|wx.ALL, self.pad)
        box_main.Add(box_plot, 1, wx.EXPAND|wx.BOTTOM|wx.LEFT|wx.RIGHT,
                     self.pad)
        #box_main.Add(box_ctrl, 0, wx.BOTTOM|wx.LEFT|wx.RIGHT,
        #             self.pad)
        
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

        # Recipes menu
        self.rec_menu = wx.Menu()

        rec_spec_extract = wx.MenuItem(self.rec_menu, self.id_spec+1,
                                       "E&xtract Spectral Region...")
        rec_line_find = wx.MenuItem(self.rec_menu, self.id_spec+2,
                                    "Find &Lines...")
        rec_line_cont = wx.MenuItem(self.rec_menu, self.id_line,
                                    "Find &Continuum by Removing Lines...")
        rec_syst_find = wx.MenuItem(self.rec_menu, self.id_cont,
                                    "Find &Systems...")
        rec_syst_fit = wx.MenuItem(self.rec_menu, self.id_syst_sel,
                                   "&Fit Selected System...")

        self.Bind(wx.EVT_MENU, self.on_spec_extract, rec_spec_extract)
        self.Bind(wx.EVT_MENU, self.on_line_find, rec_line_find)
        self.Bind(wx.EVT_MENU, self.on_line_cont, rec_line_cont)
        self.Bind(wx.EVT_MENU, self.on_syst_find, rec_syst_find)
        self.Bind(wx.EVT_MENU, self.on_syst_fit, rec_syst_fit)
        
        self.rec_menu.Append(rec_spec_extract)
        self.rec_menu.AppendSeparator()
        self.rec_menu.Append(rec_line_find)
        self.rec_menu.Append(rec_line_cont)
        self.rec_menu.AppendSeparator()
        self.rec_menu.Append(rec_syst_find)
        self.rec_menu.Append(rec_syst_fit)

        """
        if (hasattr(self, 'spec') == False):
            self.menu_disable(self.rec_menu, self.id_spec)
        if (hasattr(self, 'line') == False):
            self.menu_disable(self.rec_menu, self.id_line)
        if (hasattr(self, 'cont') == False):
            self.menu_disable(self.rec_menu, self.id_cont)
        if (hasattr(self, 'syst') == False):
            self.menu_disable(self.rec_menu, self.id_syst)            
            self.menu_disable(self.rec_menu, self.id_syst_sel)            
        """
        
        # Menu bar
        self.update_menu()
        menu_bar = wx.MenuBar()
        menu_bar.Append(self.file_menu, '&File')
        menu_bar.Append(self.rec_menu, '&Recipes')
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
            
    def on_file_open(self, event, path='.'):
        """ Behaviour for File > Open """

        wildcard = "Astrocook sessions (*.acs)|*.acs|" \
                   "FITS files (*.fits)|*.fits"
        
        # otherwise ask the user what new file to open
        with wx.FileDialog(self, "Open file", path,
                           wildcard=wildcard,
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) \
                           as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return

            name = fileDialog.GetPath()
            if (name[-4:] == '.acs'):
                self.targ = fileDialog.GetFilename()[:-24]
                try:
                    acs = self.IO.acs_read(name, path)
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

            self.update_syst()
            self.update_spec()
            self.update_plot()
            

    def on_file_save(self, event, path='.'):
        """ Behaviour for File > Save """

        timestamp = \
            '_'+str(datetime.now()).replace(" ", "_").replace(":", "-")[:-7]
        snapshot = self.targ + timestamp
        root = path + snapshot
        with wx.FileDialog(self, "Save session", path, snapshot,
                           wildcard="Astrocook session (*.acs)|*.acs",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return

            name = fileDialog.GetPath()
            try:
                acs = self
                self.IO.acs_write(acs, name, path)

            except IOError:
                wx.LogError("Cannot save session '%s'." % newfile)

    def on_line_cont(self, event):
        self.cont = Cont(self.spec, self.line)
        
        #self.line.cont()
        self.cont.smooth_maxima()

        # Only until the Cont class is rewritten
        #self.spec._cont = self.cont._y
        
        #self.ax.plot(self.spec.x, self.line._cont.y)
        self.ax.plot(self.spec.x, self.cont._t['Y'])
        self.plot_fig.draw()

        self.menu_enable(self.rec_menu, self.id_cont)
        
    def on_line_edit(self, event):
        row = event.GetRow()
        col = event.GetCol()
        label = self.line_gr.GetColLabelValue(col)
        data = self.line_gr.GetCellValue(row, col)
        self.line._t[label][row] = data
        
    def on_line_find(self, event):
        """ Behaviour for Recipes > Find Lines """
        self.params = od([('Mode:', 'abs'), ('Difference:', 'max'),
                          ('Threshold (sigma):', 5.0),
                          ('Smoothing (km/s):', 40.0)])
        dialog = ParamDialog(self, title="Find Lines")
        dialog.ShowModal()
        dialog.Destroy()
        if dialog.execute == True:
            val = self.params.values()
            self.line = Line(self.spec)
            self.line_dict[self.targ] = self.line
            self.line.find(mode=val[0], diff=val[1], kappa=float(val[2]),
                           sigma=float(val[3]))

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
            x = self.line.x[sel]
            y = self.line.y[sel]
            self.line_focus, = self.ax.plot(x, y, c='C0', marker='o', ms=20,
                                            alpha=0.2)
            self.plot_fig.draw()
            
    def on_plot_clear(self, event):
        self.ax.clear()
        self.plot_fig.draw()

    #def on_plot_draw(self, event, obj):
    def on_plot_draw(self, event):
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
        self.plot_fig.draw()
        
    def on_quit(self, event):
        """ Behaviour for File > Quit """
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
        try:
            z = self.z_dict[self.targ]
        except:
            z = 0.0
        self.params = od(
            [('Use Forest', True), ('Ion', 'Ly'), ('Emission redshift', z),
             ('Use Range', False), ('Min. wavelength', 0.0),
             ('Max. wavelength', 0.0), ('Prox. velocity', 0.0)])
        dialog = ParamDialog(self, title="Extract Spectral Region")
        dialog.ShowModal()
        dialog.Destroy()
        if dialog.execute == True:
            val = self.params.values()
            if val[0] == 'True':
                forest = self.spec.extract(forest=val[1], zem=float(val[2]))
                self.targ = self.targ + '_' + val[1]
                self.z_dict[self.targ] = float(val[2])
            else:
                forest = self.spec.extract(xmin=float(val[4])*u.nm,
                                           xmax=float(val[5])*u.nm)
                self.targ = self.targ + '_%3.0f-%3.0f' \
                            % (float(val[4]), float(val[5]))
                if float(val[2]) != 0.0:
                    self.z_dict[self.targ] = float(val[2])

            self.row = self.spec_lc.GetItemCount()
            self.spec_lc.insert_string_item(self.row, self.targ)

            self.spec = forest
            self.spec_dict[self.targ] = self.spec
            self.update_all()
            #self.update_spec()
            #self.update_line()
            #self.update_syst()
            #self.update_plot()

    def on_spec_select(self, event):
        """ Behaviour when spectrum is selected from list """

        item = self.spec_lc.GetItem(self.spec_lc.GetFirstSelected(), 0)
        self.targ = item.GetText()
        self.row = event.GetIndex()
        self.spec = self.spec_dict[self.targ]
        self.update_all()

    def on_syst_edit(self, event):
        row = event.GetRow()
        col = event.GetCol()
        label = self.syst_gr.GetColLabelValue(col)
        data = self.syst_gr.GetCellValue(row, col)
        self.syst._t[label][row] = data
        
    def on_syst_find(self, event):
        """ Behaviour for Recipes > Find Lines """
        self.params = od([('series', 'CIV')])
        dialog = ParamDialog(self, title="Find Lines")
        dialog.ShowModal()
        dialog.Destroy()
        if dialog.execute == True:
            syst = System(self.spec, self.line, self.cont)
            syst.find(tag=self.params['series'])

            if self.syst != None:
                self.syst.merge(syst)
            else:
                self.syst = syst

            self.syst_dict[self.targ] = self.syst
            self.update_spec()
            self.update_syst()

        # Only until the Cont class is rewritten
        #self.syst._cont = self.line._cont
        #self.syst._cont = self.cont._y
            
    def on_syst_fit(self, event):
        self.syst = self.syst_dict[self.targ]
        
        dialog = SystDialog(self, title="Fit selected system")
        dialog.ShowModal()
        dialog.Destroy()
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
                 for i in self.syst.t['SERIES'][sel]]
            dx = 0.5
            h = np.max(self.spec.y)
            #self.syst.group(z, dx)
            self.syst.model(z, dx)
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
                        
    
            
    def update_plot(self):
        """ Update the plot panel """

        self.on_plot_clear(None)
        self.on_plot_draw(None)

    def update_spec(self):
        """ Update the spec list """

        self.spec = self.spec_dict[self.targ]

        try:
            self.spec_lc.SetItem(self.row, 2, str(self.z_dict[self.targ]))
        except:
            pass

        xmin = self.spec.t['X'][0]
        xmax = self.spec.t['X'][-1]
        self.spec_lc.SetItem(self.row, 3, "[%3.2f, %3.2f]" % (xmin, xmax))

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

    def __init__(self, parent=None, title="Parameters", size=(250,500),
                 **kwargs):
        """ Constructor for the ParamDialog class """
        super(ParamDialog, self).__init__(parent, title=title)#, size=size) 

        self.p = parent
        self.init_UI()
        #self.SetSize((250, 200))
        #self.SetTitle("Change Color Depth")

        
    def init_UI(self):
        """ Initialize the main frame """
        
        panel = wx.Panel(self)
        box_main = wx.BoxSizer(wx.VERTICAL)

        box_params = wx.BoxSizer(wx.VERTICAL)
        self.par = []
        self.ctrl = []
        for p, v in self.p.params.iteritems():
            box_param = wx.BoxSizer(wx.HORIZONTAL)
            self.par.append(p)
            if type(v) == bool:
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
        self.Destroy()

    def on_run(self, e):
        for p, ctrl in zip(self.par, self.ctrl):
            self.p.params[p] = str(ctrl.GetValue())
        self.execute = True
        self.Destroy()
                      
class SystDialog(wx.Dialog):

    def __init__(self, parent=None, title="Parameters", **kwargs):
        """ Constructor for the ParamDialog class """

        size = (wx.DisplaySize()[0]*0.5, wx.DisplaySize()[1]*0.7)
        super(SystDialog, self).__init__(parent, title=title, size=size) 

        self.p = parent
        self.syst = self.p.syst#_sel
        self.group = self.p.syst._group
        self.z = self.p.z_sel
        self.sel = np.where(self.syst._t['Z'] == self.z)[0]
        self.series = self.syst._t['SERIES'][self.sel][0]
        self.init_UI()

    def init_buttons(self, panel):
        self.syst_b = wx.Button(panel, label="Add system", size=(100,38))
        self.line_b = wx.Button(panel, label="Add line", size=(100,38))
        self.syst_b.Bind(wx.EVT_BUTTON, self.on_syst_add)

        
    def init_plot(self, panel):
        """ Create the spectrum panel """
        self.fig = Figure()#figsize=(20,20))


        
        #sel = np.where(self.p.syst._t['Z'] == self.p.z_sel)[0]
        
        #ions = self.p.syst._t['SERIES'][sel][0] #self.p.syst._group['ION']
        #waves = [dict_wave[i].value for i in ions]
        #ions = ions[np.argsort(waves)]
        rown = 5.
        self.pn = len(self.series)
        row = min(self.pn,rown)
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
        grid.tight_layout(self.fig, rect=[0.01, 0.01, 1, 0.9], h_pad=0.0)
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

        self.focus_gr = self.init_tab(panel)
        self.add_gr = self.init_tab(panel)
        self.init_buttons(panel)
        self.add_gr.SetColLabelSize(0)
        self.init_plot(panel)
        #self.update_group()

        self.box_main = wx.BoxSizer(wx.VERTICAL)

        box_focus = wx.BoxSizer(wx.VERTICAL)
        #box_focus.Add(wx.StaticText(panel, label="Focus"))
        box_focus.Add(self.focus_gr, 1, wx.BOTTOM, border=5)

        box_button = wx.BoxSizer(wx.HORIZONTAL)
        box_button.Add(self.syst_b, 0, wx.RIGHT, border=5)
        box_button.Add(self.line_b, 0, wx.RIGHT, border=5)

        box_add = wx.BoxSizer(wx.VERTICAL)
        #box_add.Add(wx.StaticText(panel, label="Additional lines"))
        box_add.Add(self.add_gr, 1, wx.BOTTOM, border=10)
        box_add.Add(box_button)

        box_plot = wx.BoxSizer(wx.VERTICAL)
        box_plot.Add(self.plot_fig, 1, wx.EXPAND)

        box_ctrl = wx.BoxSizer(wx.HORIZONTAL)
        box_ctrl.Add(self.plot_tb, 1)

        box_disp = wx.BoxSizer(wx.VERTICAL)
        box_disp.Add(box_focus)
        box_disp.Add(box_add)
        box_disp.Add(box_plot, 1, wx.EXPAND|wx.TOP|wx.BOTTOM, border=10)
        box_disp.Add(box_ctrl)


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
        self.Destroy()

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
        #print z
        dx = 0.5
        group_z = np.where(group['Z'] == z)[0]
        print group_z
        group[label][group_z] = data
        if label in t.colnames:
            t[label][np.where(t['Z'] == z)[0]] = data
        if label in line.colnames:
            line[label][np.where(map['Z'] == z)[0]] = data
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
        print group
        self.syst.model(self.z, dx)
        self.update_plot()
        
    def on_run(self, event):
        self.execute = True
        self.Destroy()

    def on_syst_add(self, event):
        z_shift = 0.001*random.random()

        chunk = self.syst._chunk
        group = self.syst._group
        x_min = chunk['X'][(chunk['Y']-chunk['MODEL']).argmin()]
        z_min_arr = np.array([])
        for r in range(self.focus_gr.GetNumberRows()):
            ion = self.focus_gr.GetCellValue(r, 3)
            z_min_arr = np.append(z_min_arr, x_min/dict_wave[ion].value - 1)
        z_min = z_min_arr[np.abs(z_min_arr-self.z).argmin()]
        print z_min
            
        dupl = self.syst._t[self.sel]
        #dupl['Z'] = dupl['Z']+z_shift
        dupl['Z'] = z_min
        print dupl

        match = np.where(self.syst._map['Z'] == self.z)[0]
        for r in self.syst._map[match]:
            self.syst._map.add_row(r)
            self.syst._map['Z'][-1] = z_min#self.syst._map['Z'][-1]+z_shift

        self.syst._map.sort('X')
        
        self.syst._t = vstack([self.syst._t, dupl])#merge(dupl)
        self.syst.group(self.z, 0.5)
        dx = 0.5
        self.syst.model(self.z, dx)
        self.update_tab()
        self.update_plot()
        
    def update_plot(self):
        for p in range(self.pn):
            self.ax[p].clear()
        self.syst.plot(z=self.z, ax=self.ax)
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
            if g['Z'] == self.z:
                tab = self.focus_gr
            else:
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
