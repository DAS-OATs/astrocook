class Plot():

    def __init__(self, ax):
        self.ax = ax

        # Lists of all previously plotted data, to clear them selectively
        self.cont_p = []
        self.line_p = []
        self.spec_p = []

    def cont(self, tab, replace=False, c='C2', lw=2.0, **kwargs):
        try:
            for p in self.cont_p:
                p.remove()                
            self.cont_p = []
        except:
            pass

        (p, p_other) = self.spec(tab, replace, dy=False, c=c, lw=lw, **kwargs)
        self.cont_p.append(p)

        return p
        
    def line(self, tab, replace=True, s=100, c='r', marker='+', fill=False):
        
        # Clear all previously plotted lines, if requested
        if replace:
            try:
                for p in self.line_p:
                    p.remove()
                self.line_p = []
            except:
                pass

        # Plot new lines
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
            
    def spec(self, tab, replace=True, dy=True, xmin=None, xmax=None, c='C0',
             c_other='C1',  lw=1.0, **kwargs):

        # Clear all previously plotted spectra, if requested
        if replace:
            try:
                for p in self.spec_p:
                    p.remove()
                self.spec_p = []
            except:
                pass

        if xmin != None and xmax != None:
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

            
