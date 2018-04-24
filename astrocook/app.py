from . import Line, Spec1DReader, System
from .utils import *
from astropy.io import fits
from collections import OrderedDict as od
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg
from matplotlib.figure import Figure
import numpy as np
import os
import sys
import tarfile
import wx
import wx.grid as gridlib
import wx.lib.mixins.listctrl as listmix

class MainFrame(wx.Frame):
    def __init__(self, parent=None, title="Astrocook", **kwargs):
        """ Constructor for the Frame class """ 

        size = (wx.DisplaySize()[0]*0.708, wx.DisplaySize()[1]*0.9)
        self.pad = 10
        super(MainFrame, self).__init__(parent, title=title, size=size)
        self.init_UI(**kwargs)

        self.spec_dict = {}
        self.line_dict = {}
        self.syst_dict = {}
        
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
        box_ctrl.Add(self.plot_tb, 1)
        box_ctrl.Add(self.plot_pb, 0, wx.ALIGN_RIGHT)
        box_ctrl.Add(self.plot_cb, 0, wx.ALIGN_RIGHT)
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
        self.spec_lc.Bind(wx.EVT_LIST_END_LABEL_EDIT, self.on_spec_edit)
        self.spec_lc.Bind(wx.EVT_LIST_ITEM_SELECTED, self.on_spec_select)
        self.spec_lc.InsertColumn(0, 'path', width=200)
        self.spec_lc.InsertColumn(1, 'object', width=100)
        self.spec_lc.InsertColumn(2, 'redshift')
        #self.spec_lc.Select(0)
        
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
        
    def menu(self, **kwargs):
        """ Create a menu in the frame """

        # File menu
        file_menu = wx.Menu()
        
        file_open = wx.MenuItem(file_menu, wx.ID_OPEN, "&Open\tCtrl+O")
        file_save = wx.MenuItem(file_menu, wx.ID_SAVE, "&Save\tCtrl+S")
        file_quit = wx.MenuItem(file_menu, wx.ID_EXIT, "&Quit\tCtrl+Q")
        self.Bind(wx.EVT_MENU, lambda e: self.on_file_open(e, **kwargs),
                  file_open)
        self.Bind(wx.EVT_MENU, lambda e: self.on_file_save(e, **kwargs),
                  file_save)
        self.Bind(wx.EVT_MENU, self.on_quit, file_quit)

        file_menu.Append(file_open)
        file_menu.AppendSeparator()
        file_menu.Append(file_save)
        file_menu.AppendSeparator()
        file_menu.Append(file_quit)

        # Recipes menu
        rec_menu = wx.Menu()

        rec_spec_extract = wx.MenuItem(rec_menu, wx.ID_ANY,
                                       "E&xtract Spectral Region...")
        rec_line_find = wx.MenuItem(rec_menu, wx.ID_ANY, "Find &Lines...")
        rec_line_cont = wx.MenuItem(rec_menu, wx.ID_ANY,
                                    "Find &Continuum by Removing Lines...")
        rec_syst_find = wx.MenuItem(rec_menu, wx.ID_ANY, "Find &Systems...")
        rec_syst_fit = wx.MenuItem(rec_menu, wx.ID_ANY,
                                   "&Fit Selected System...")
        self.Bind(wx.EVT_MENU, self.on_spec_extract, rec_spec_extract)
        self.Bind(wx.EVT_MENU, self.on_line_find, rec_line_find)
        self.Bind(wx.EVT_MENU, self.on_line_cont, rec_line_cont)
        self.Bind(wx.EVT_MENU, self.on_syst_find, rec_syst_find)
        self.Bind(wx.EVT_MENU, self.on_syst_fit, rec_syst_fit)
        
        rec_menu.Append(rec_spec_extract)
        rec_menu.AppendSeparator()
        rec_menu.Append(rec_line_find)
        rec_menu.Append(rec_line_cont)
        rec_menu.AppendSeparator()
        rec_menu.Append(rec_syst_find)
        rec_menu.Append(rec_syst_fit)
        
        # Menu bar
        menu_bar = wx.MenuBar()
        menu_bar.Append(file_menu, '&File')
        menu_bar.Append(rec_menu, '&Recipes')
        self.SetMenuBar(menu_bar)        

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
                self.session = fileDialog.GetFilename()[:-24]
                try:
                    with tarfile.open(name) as arch:
                        arch.extractall(path=path)

                        # Open spectrum
                        self.spec_name = name[:-4]+'_spec.fits'
                        self.spec = Spec1DReader().ac(self.spec_name)
                        self.z = 0.0
                        self.spec_dict[self.spec_name] = self.spec 
                        self.spec_lc.insert_string_item(
                            self.spec_lc.GetItemCount(), self.spec_name)
                        self.spec_lc.SetStringItem(
                            self.spec_lc.GetItemCount()-1, 2, str(self.z))
                        #self.on_plot_draw(None, self.spec)
                        os.remove(self.spec_name)

                        # Open line list
                        try:
                            line_name = name[:-4]+'_line.fits'

                            # Put in a dedicated class
                            hdulist =  fits.open(line_name)
                            data = hdulist[1].data
                            names = np.array(hdulist[1].columns.names)
                            units = np.array(hdulist[1].columns.units)
                            xunit = units[np.where(names == 'X')][0]
                            yunit = units[np.where(names == 'Y')][0]
                            x = data['X']
                            xmin = data['XMIN']
                            xmax = data['XMAX']
                            y = data['Y']            
                            dy = data['DY']
                            self.line = Line(x=x, xmin=xmin, xmax=xmax, y=y,
                                             dy=dy, xunit=xunit, yunit=yunit)
                            self.line_dict[self.spec_name] = self.line
                            self.refresh_line()
                            os.remove(line_name)
                        except:
                            pass
                        
                        # Open system list
                        try:
                            syst_name = name[:-4]+'_syst.fits'

                            # Put in a dedicated class
                            hdulist =  fits.open(syst_name)
                            data = hdulist[1].data
                            names = np.array(hdulist[1].columns.names)
                            units = np.array(hdulist[1].columns.units)
                            Nunit = units[np.where(names == 'N')][0]
                            bunit = units[np.where(names == 'B')][0]
                            series = data['SERIES']
                            z = data['Z']
                            N = data['N']
                            b = data['B']
                            btur = data['BTUR']
                            dz = data['DZ']
                            dN = data['DN']
                            db = data['DB']
                            dbtur = data['DBTUR']
                            vary = data['VARY']
                            expr = data['EXPR']
                            self.syst = System(
                                self.spec, series=series,
                                z=z, N=N, b=b, btur=btur,
                                dz=dz, dN=dN, db=db, dbtur=dbtur,
                                vary=vary, expr=expr, Nunit=Nunit, bunit=bunit)
                            self.syst_dict[self.spec_name] = self.syst
                            self.refresh_syst()
                            os.remove(syst_name)
                        except:
                            pass

                        self.on_plot_draw(None)
                    
                except IOError:
                    wx.LogError("Cannot open file '%s'." % arch)
                
            else:
                self.session = fileDialog.GetFilename()[:-5] 
                try:
                    with open(name, 'r') as file:
                        try:
                            self.spec = Spec1DReader().ac(name)
                        except:
                            self.spec = Spec1DReader().xq(name)
                        self.z = 0.0
                        self.spec_dict[name] = self.spec 
                        self.spec_lc.insert_string_item(
                            self.spec_lc.GetItemCount(), name)
                        self.spec_lc.SetStringItem(
                            self.spec_lc.GetItemCount()-1, 2, str(self.z))
                        #self.on_plot_draw(None, self.spec)
                        self.on_plot_draw(None)
                    
                except IOError:
                    wx.LogError("Cannot open file '%s'." % file)

    def on_file_save(self, event, path='.'):
        """ Behaviour for File > Save """

        timestamp = \
            '_'+str(datetime.now()).replace(" ", "_").replace(":", "-")[:-7]
        snapshot = self.session + timestamp
        root = path + snapshot
        with wx.FileDialog(self, "Save session", path, snapshot,
                           wildcard="Astrocook session (*.acs)|*.acs",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return

            try:
                with tarfile.open(root+'.acs', 'w:gz') as arch:
                    spec_name = root + '_spec.fits'
                    spec_arcname = snapshot + '_spec.fits'
                    self.spec.save(spec_name)
                    arch.add(spec_name, arcname=spec_arcname)
                    os.remove(spec_name)
                    if hasattr(self, 'line'):
                        line_name = root + '_line.fits'
                        line_arcname = snapshot + '_line.fits'
                        self.line.save(line_name)
                        arch.add(line_name, arcname=line_arcname)
                        os.remove(line_name)
                    if hasattr(self, 'syst'):
                        syst_name = root + '_syst.fits'
                        syst_arcname = snapshot + '_syst.fits'
                        self.syst.save(syst_name)
                        arch.add(syst_name, arcname=syst_arcname)
                        os.remove(syst_name)

            except IOError:
                wx.LogError("Cannot save session '%s'." % newfile)

    def on_line_cont(self, event):
        self.line.cont()

        # Only until the Cont class is rewritten
        self.spec._cont = self.line._cont
        
        self.ax.plot(self.spec.x, self.line._cont.y)
        self.plot_fig.draw()

    def on_line_find(self, event):
        """ Behaviour for Recipes > Find Lines """
        self.params = od([('kappa', 5.0), ('sigma', 40.0)])
        param = ParamDialog(self, title="Find Lines")
        param.ShowModal()
        param.Destroy()
        if param.execute == True:
            self.line = Line(self.spec)
            self.line_dict[self.spec_name] = self.line
            self.line.find(kappa=float(self.params['kappa']),
                           sigma=float(self.params['sigma']))

        self.refresh_line()
        self.refresh_plot()
        
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
        self.spec.plot(ax=self.ax)
        try:
            self.line.plot_new(ax=self.ax)
        except:
            pass
        self.plot_fig.draw()
        
    def on_quit(self, event):
        """ Behaviour for File > Quit """
        self.Close()

    def on_spec_edit(self, event):
        """ Behaviour when spectrum is selected from list """

        index = self.spec_lc.GetFocusedItem()
        row = event.GetIndex()
        col = event.GetColumn()
        data = event.GetLabel()
        self.spec_lc.SetStringItem(row, col, data)
        self.z = self.spec_lc.GetItem(index, 2).GetText()
        
    def on_spec_extract(self, event):
        self.params = od([('zem', self.z), ('forest', 'Ly')])#('Lyman', True), ('CIV', False)])
#        xmin=None, xmax=None, prox=False, forest=[], zem=[],
#                prox_vel=[]
        param = ParamDialog(self, title="Extract Spectral Region")
        param.ShowModal()
        param.Destroy()
        if param.execute == True:
            forest = self.spec.extract(forest=str(self.params['forest']),
                                       zem=float(self.params['zem']))
            self.spec = forest
            self.spec.plot(ax=self.ax)
            self.plot_fig.draw()

    def on_spec_select(self, event):
        """ Behaviour when spectrum is selected from list """

        item = self.spec_lc.GetItem(self.spec_lc.GetFirstSelected(), 0)
        self.spec_name = item.GetText()
        self.spec = self.spec_dict[self.spec_name]
        try:
            self.line = self.line_dict[self.spec_name]
            self.refresh_line()
        except:
            pass
        
    def on_syst_find(self, event):
        """ Behaviour for Recipes > Find Lines """
        self.params = od([('series', 'CIV')])
        param = ParamDialog(self, title="Find Lines")
        param.ShowModal()
        param.Destroy()
        if param.execute == True:
            self.syst = System(self.spec, self.line)#, doubl=self.params['doubl'])
            self.syst.find(tag=self.params['series'])
            self.refresh_syst()

        # Only until the Cont class is rewritten
        self.syst._cont = self.line._cont
            
    def on_syst_fit(self, event):
        syst_sel = self.syst.extract(self.syst_sel)

        # Only until the Cont class is rewritten
        syst_sel._cont = self.syst._cont

        z_sel = syst_sel.t['Z']
        syst_sel.fit(z_sel, norm=False)
        chunk = syst_sel._chunk
        self.ax.plot(chunk['X'], chunk['MODEL'])
        self.plot_fig.draw()

        for c in ['Z', 'N', 'B', 'BTUR', 'DZ', 'DN', 'DB', 'DBTUR', 'VARY',
                  'EXPR']:
            self.syst._t[self.syst_sel][c] = syst_sel.t[c]
        self.refresh_syst()
        
    def on_syst_select(self, event):
        """ Behaviour for line selection """        
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
            self.syst_focus = self.ax.bar(x, h, dx, 0, color='C1', alpha=0.2)
            self.plot_fig.draw()
            self.syst_sel = sel
            
        
    def refresh_line(self):
        try:
            self.line_gr.DeleteRows(pos=0, numRows=self.line_gr.GetNumberRows())
        except:
            pass
        self.line_gr.AppendRows(len(self.line.t))
        for i, l in enumerate(self.line.t):
            self.line_gr.SetCellValue(i, 0, "%3.3f" % l['X'])
            self.line_gr.SetCellValue(i, 1, "%3.3f" % l['XMIN'])
            self.line_gr.SetCellValue(i, 2, "%3.3f" % l['XMAX'])
            self.line_gr.SetCellValue(i, 3, "%3.3f" % l['Y'])
            self.line_gr.SetCellValue(i, 4, "%3.3f" % l['DY'])
        #self.line.plot_new(ax=self.ax)
        #self.plot_fig.draw()

    def refresh_plot(self):
        self.on_plot_clear(None)
        self.on_plot_draw(None)
        
    def refresh_syst(self):
        try:
            self.syst_gr.DeleteRows(pos=0, numRows=self.syst_gr.GetNumberRows())
        except:
            pass
        self.syst_gr.AppendRows(len(self.syst.t))
        for i, l in enumerate(self.syst.t):
            self.syst_gr.SetCellValue(i, 0, str(l['SERIES']))
            self.syst_gr.SetCellValue(i, 1, "%3.3f" % l['Z'])
            self.syst_gr.SetCellValue(i, 2, "%3.3e" % l['N'])
            self.syst_gr.SetCellValue(i, 3, "%3.3f" % l['B'])
            self.syst_gr.SetCellValue(i, 4, "%3.3f" % l['BTUR'])
            self.syst_gr.SetCellValue(i, 5, "%3.3e" % l['DZ'])
            self.syst_gr.SetCellValue(i, 6, "%3.3f" % l['DN'])
            self.syst_gr.SetCellValue(i, 7, "%3.3f" % l['DB'])
            self.syst_gr.SetCellValue(i, 8, "%3.3f" % l['DBTUR'])
        

class EditableListCtrl(wx.ListCtrl, listmix.TextEditMixin):
    def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        """ Constructor for the EditableListCtrl class """
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.TextEditMixin.__init__(self)


    def insert_string_item(self, *args):
        self.InsertStringItem(*args)
        listmix.TextEditMixin.__init__(self)
        
class ParamDialog(wx.Dialog):

    def __init__(self, parent=None, title="Parameters", size=(250,500),
                 **kwargs):
        """ Constructor for the ParamDialog class """
        super(ParamDialog, self).__init__(parent, title=title)#, size=size) 

        self.parent = parent
        self.init_UI()
        #self.SetSize((250, 200))
        #self.SetTitle("Change Color Depth")

        
    def init_UI(self):
        
        panel = wx.Panel(self)
        box_main = wx.BoxSizer(wx.VERTICAL)

        box_params = wx.BoxSizer(wx.VERTICAL)
        self.p = []
        self.ctrl = []
        for p, v in self.parent.params.iteritems():
            box_param = wx.BoxSizer(wx.HORIZONTAL)
            self.p.append(p)
            if type(v) == bool:
                rb = wx.RadioButton(panel, label=p)
                box_param.Add(rb, 1, wx.LEFT, border=5)
                self.ctrl.append(rb)
            else:
                st = wx.StaticText(panel, -1, label=p)
                tc = wx.TextCtrl(panel, -1, value=str(v))
                box_param.Add(st, 1, 0)
                box_param.Add(tc, 1, 0)
                self.ctrl.append(tc)


            box_params.Add(box_param, 1, 0, 0)

        panel.SetSizer(box_params)
        #for tc, p in zip(tc_list, p_list):
        #    print(tc, p)
        #    self.Bind(wx.EVT_TEXT, lambda e, p=p: self.on_edit(e, p), tc)
        
        #sb = wx.StaticBox(panel, label='Colors')
        #sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)        
        #sbs.Add(wx.RadioButton(panel, label='256 Colors', 
        #    style=wx.RB_GROUP))
        #sbs.Add(wx.RadioButton(panel, label='16 Colors'))
        #sbs.Add(wx.RadioButton(panel, label='2 Colors'))
        
        #box_main.Add(wx.TextCtrl(panel), flag=wx.LEFT, border=5)
        
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
        for p, ctrl in zip(self.p, self.ctrl):
            self.parent.params[p] = str(ctrl.GetValue())
        self.execute = True
        self.Destroy()
                      
        
