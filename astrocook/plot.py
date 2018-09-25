class Plot():

    def __init__(self, ax, xmin=None, xmax=None):
        self.ax = ax
        self.xmin = xmin
        self.xmax = xmax
    
    def spec(self, tab, dy=True):
        try:
            self.ax.set_xlabel("Wavelength [" + str(spec.x.unit) + "]")
        except:
            pass
        try:
            self.ax.set_ylabel("Flux [" + str(spec.y.unit) + "]")
        except:
            pass
        self.ax.plot(tab['X'], tab['Y'], lw=1.0)
        if dy:
            self.ax.plot(tab['X'], tab['DY'], lw=1.0)

    def line(self, tab, s=100, c='r', marker='+', fill=False):
        if marker == 'o':
            if fill == True:
                self.ax.scatter(tab['X'], tab['Y'], s=s, marker=marker,
                                facecolors=c, edgecolors=c)
            else:
                self.ax.scatter(tab['X'], tab['Y'], s=s, marker=marker,
                                facecolors='none', edgecolors=c)
        else:
            self.ax.scatter(tab['X'], tab['Y'], s=s, c=c, marker=marker)
