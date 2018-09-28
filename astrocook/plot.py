from .utils import dict_wave, xunit_def, yunit_def
import numpy as np

class Plot():

    def __init__(self, ax, xlabel=None, ylabel=None):
        self.ax = ax

        # Axis labels
        if xlabel == None:
            try:
                self.ax.set_xlabel("Wavelength [" + str(xunit_def) + "]")
            except:
                pass
        else:
            self.ax.set_xlabel(xlabel)
        if ylabel == None:
            try:
                self.ax.set_ylabel("Flux [" + str(yunit_def) + "]")
            except:
                pass
        else:
            self.ax.set_ylabel(ylabel)
        
        (self.xmin, self.xmax) = self.ax.get_xlim()
        
        # Lists of all previously plotted data, to clear them selectively
        self.cont_p = []
        self.model_p = []
        self.line_p = []
        self.spec_p = []
        self.sel_p = []

    def clean(self, attr):
        """ Clear all previously plotted data from a list, if requested """

        attr_p = getattr(self, attr)
        try:
            for p in attr_p:
                p.remove()                
            setattr(self, attr, [])
        except:
            pass

    def clear(self):
        xlabel = self.ax.get_xlabel()
        ylabel = self.ax.get_ylabel()
        self.ax.clear()
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        
    def cont(self, tab, replace=True, cont=None, ion=None, c='C2', lw=1.0,
             **kwargs):
        if replace:
            self.clean('cont_p')
        
        (p, p_other) = self.spec(tab, replace=False, cont=cont, ion=ion,
                                 dy=False, c=c, lw=lw, **kwargs)
        self.cont_p.append(p)

        return p

    def model(self, tab, replace=True, cont=None, ion=None, c='C4', lw=1.0,
            **kwargs):
        if replace:
            self.clean('model_p')
        
        (p, p_other) = self.spec(tab, replace=False, cont=cont, ion=ion,
                                 dy=False, c=c, lw=lw, **kwargs)
        self.model_p.append(p)

        return p
        
    def line(self, tab, replace=True, cont=None, ion=None, s=100, c='C3',
             marker='+', fill=False, **kwargs):
        if replace:
            self.clean('line_p')

        tab_x = self.x_mode(tab['X'], ion)
        (tab_y, tab_dy) = self.y_mode(tab, cont)
        if marker == 'o':
            if fill:
                p = self.ax.scatter(tab_x, tab_y, s=s, marker=marker,
                                    facecolors=c, edgecolors=c, **kwargs)
            else:
                p = self.ax.scatter(tab_x, tab_y, s=s, marker=marker,
                                    facecolors='none', edgecolors=c,
                                    **kwargs)
        else:
            p = self.ax.scatter(tab_x, tab_y, s=s, c=c, marker=marker, **kwargs)

        self.line_p.append(p)

        if ion:
            self.ax.text(0.05, 0.5, ion, transform=self.ax.transAxes,
                         fontsize=13)
        
        return p

    def sel(self, obj, rows, replace=True, ion=None, extra_width=1.0, c='C3',
            **kwargs):
        if replace:
            self.clean('sel_p')

        xmins = []
        xmaxs = []
        for r in rows:
            x = self.x_mode(obj.x[r], ion)
            xmin = self.x_mode(obj.xmin[r], ion)
            xmax = self.x_mode(obj.xmax[r], ion)
            xmins.append(xmin.value)
            xmaxs.append(xmax.value)
        
            p = self.ax.axvline(x=x.value, color=c, **kwargs)
            pmin = self.ax.axvline(x=xmin.value, color=c, linestyle=':',
                                   **kwargs)
            pmax = self.ax.axvline(x=xmax.value, color=c, linestyle=':',
                                   **kwargs)

            self.sel_p.append(p)
            self.sel_p.append(pmin)
            self.sel_p.append(pmax)

        if extra_width > 0:
            xmin_a = min(xmins)
            xmax_a = max(xmaxs)
            w = extra_width*(xmax_a-xmin_a)
            self.ax.set_xlim(xmin_a-w, xmax_a+w)
        else:
            self.ax.set_xlim(self.xmin, self.xmax)
        
    def spec(self, tab, replace=True, dy=True, cont=None, ion=None, xmin=None,
             xmax=None, c='C0', c_other='C1', lw=1.0, **kwargs):
        if replace:
            self.clean('spec_p')

        if xmin != None and xmax != None:
            self.xmin = xmin
            self.xmax = xmax
            self.ax.set_xlim(xmin, xmax)
            
        tab_x = self.x_mode(tab['X'], ion)
        (tab_y, tab_dy) = self.y_mode(tab, cont)
        p, = self.ax.plot(tab_x, tab_y, c=c, lw=lw, **kwargs)
        self.spec_p.append(p)
        if dy:
            p_other, = self.ax.plot(tab_x, tab_dy, c=c_other, lw=lw, **kwargs)
            self.spec_p.append(p_other)
        else:
            p_other = None

        return (p, p_other)
            
    def x_mode(self, x, ion):
        if ion != None:
            x = x/dict_wave[ion]-1
        return x

    def y_mode(self, tab, cont):
        tab_y = tab['Y']
        tab_dy = tab['DY']
        if cont != None:
            if (len(cont) == len(tab_y)):
                cont_y = cont['Y']
            else:
                cont_y = np.interp(tab['X'], cont['X'], cont['Y'])
            tab_y = tab_y/cont_y
            tab_dy = tab_dy/cont_y
        return (tab_y, tab_dy)

