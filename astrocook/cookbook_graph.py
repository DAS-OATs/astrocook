from .graph import GraphCursorZSeries

class CookbookGraph(object):
    """ Cookbook of graph utilities
    """

    def __init__(self):
        super(CookbookGraph, self).__init__()


    def menu_view(self, attr):
        getattr(self.sess._gui._panel_sess._menu._view, attr)(event=None)


    def menu_view_cursor_z_series(self, z=2.0, series="Ly_a"):
        gui = self.sess._gui
        menu = gui._panel_sess._menu
        menu._on_graph(None, "Redshift cursor", 'cursor_z_series', gui._cursor, True, GraphCursorZSeries)
        gui._dlg_mini._ctrl_z.SetValue(str(z))
        gui._dlg_mini._ctrl_series.SetValue(str(series))
        gui._dlg_mini._on_apply(None)
        """
        if not hasattr(gui, '_dlg_mini') \
            or gui._dlg_mini == None:
            gui._panel_sess._menu._on_dialog_mini(None, "System controls", GraphCursorZSeries)
        gui._dlg_mini._shown = True
        gui._dlg_mini._ctrl_z.SetValue(str(z))
        gui._dlg_mini._ctrl_series.SetValue(str(series))
        print(gui._dlg_mini._targ)
        gui._dlg_mini._on_apply(None)
        gui._dlg_mini._cursor_button.SetLabel("Hide cursor")
        """
        """
        self.sess._gui._sess_sel._series_sel = series
        GraphCursorZSeries(self.sess._gui._sess_sel)
        self.sess._gui._refresh(init_cursor=True, init_tab=False)
        """
