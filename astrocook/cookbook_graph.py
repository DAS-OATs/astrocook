from .graph import GraphCursorZSeries

class CookbookGraph(object):
    """ Cookbook of graph utilities

    This utilities are not directly called from the GUI, but are meant as
    wrapper to make GUI functionalities available also through JSON
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
