class Plot():

    def __init__(self, ax):
        self.ax = ax

        (self.xmin, self.xmax) = self.ax.get_xlim()
        
        # Lists of all previously plotted data, to clear them selectively
        self.cont_p = []
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
        
    def cont(self, tab, replace=True, c='C2', lw=2.0, **kwargs):
        if replace:
            self.clean('cont_p')
        
        (p, p_other) = self.spec(tab, replace=False, dy=False, c=c, lw=lw,
                                 **kwargs)
        self.cont_p.append(p)

        return p
        
    def line(self, tab, replace=True, s=100, c='C3', marker='+', fill=False):
        if replace:
            self.clean('line_p')
        
        if marker == 'o':
            if fill:
                p = self.ax.scatter(tab['X'], tab['Y'], s=s, marker=marker,
                                    facecolors=c, edgecolors=c)
            else:
                p = self.ax.scatter(tab['X'], tab['Y'], s=s, marker=marker,
                                    facecolors='none', edgecolors=c)
        else:
            p = self.ax.scatter(tab['X'], tab['Y'], s=s, c=c, marker=marker)

        self.line_p.append(p)
        return p

    def sel(self, obj, rows, replace=True, extra_width=1.0, c='C3', **kwargs):
        if replace:
            self.clean('sel_p')

        xmins = []
        xmaxs = []
        for r in rows:
            x = obj.x[r]
            y = obj.y[r]
            xmin = obj.xmin[r]
            xmax = obj.xmax[r]
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
        
    def spec(self, tab, replace=True, dy=True, xmin=None, xmax=None, c='C0',
             c_other='C1',  lw=1.0, **kwargs):
        if replace:
            self.clean('spec_p')

        if xmin != None and xmax != None:
            self.xmin = xmin
            self.xmax = xmax
            self.ax.set_xlim(xmin, xmax)
            
        # Plot new spectra
        try:
            self.ax.set_xlabel("Wavelength [" + str(tab['X'].unit) + "]")
        except:
            pass
        try:
            self.ax.set_ylabel("Flux [" + str(tab['Y'].unit) + "]")
        except:
            pass
        p, = self.ax.plot(tab['X'], tab['Y'], c=c, lw=lw, **kwargs)
        self.spec_p.append(p)
        if dy:
            p_other, = self.ax.plot(tab['X'], tab['DY'], c=c_other, lw=lw, **kwargs)
            self.spec_p.append(p_other)
        else:
            p_other = None

        return (p, p_other)

            
