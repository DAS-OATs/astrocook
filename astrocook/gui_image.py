from astropy import units as au
import wx
import matplotlib.pyplot as plt

prefix = "GUI:"

class GUIImage(wx.Frame):
    """ Class for the GUI image frame """

    def __init__(self,
                 gui,
                 attr,
                 title="Image",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.5):

        super(GUIImage, self).__init__(parent=None, title=title,
                                       size=(size_x, size_y))

        self._gui = gui
        self._attr = attr
        self._title = title

    def _on_view(self, event):
        data = getattr(self._gui._sess_sel, self._attr)
        extent = getattr(self._gui._sess_sel, self._attr+'_e')
        fig, ax = plt.subplots()
        fig.canvas.set_window_title(self._title)
        im = ax.imshow(data, aspect='auto', extent=extent)
        fig.colorbar(im)
        ax.set_xlabel('b', str(au.km/au.s))
        ax.set_ylabel('log N/'+str(au.cm**2))
        plt.show()


class GUIImageCompleteness(GUIImage):
    """ Class for the GUI system completeness image """

    def __init__(self,
                 gui,
                 title="System completeness",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.5):

        super(GUIImageCompleteness, self).__init__(gui, 'compl', title, size_x,
                                                  size_y)

        self._gui = gui
        self._title = title
        self._gui._ima_compl = self

class GUIImageCorrectness(GUIImage):
    """ Class for the GUI system correctness image """

    def __init__(self,
                 gui,
                 title="System correctness",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.5):

        super(GUIImageCorrectness, self).__init__(gui, 'corr', title, size_x,
                                                  size_y)

        self._gui = gui
        self._title = title
        self._gui._ima_corr = self
