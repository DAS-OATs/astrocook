from astrocook.cookbook import Cookbook
from astrocook.line_list import LineList
from astrocook.session import Session
from astrocook.spectrum import Spectrum
from astrocook.syst_model import SystModel
from astrocook.vars import *
import astropy.units as au
import matplotlib.pyplot as plt
import numpy as np

xstep = 0.003
xrange = np.arange(300, 302, xstep)
n_lines = 10
texp = 2
path = 'outfile.acs'

x_lines = np.random.random(n_lines)*(xrange[-1]-xrange[0])+xrange[0]
logN_lines = -np.random.power(3, n_lines)*5+17
b_lines = np.random.poisson(10, n_lines)+np.random.random()-0.5

s = Session()
s.spec = Spectrum(x=xrange,
                     xmin=xrange-xstep*0.5,
                     xmax=xrange+xstep*0.5,
                     y=np.ones(xrange.shape),
                     dy=np.ones(xrange.shape)/(50*np.sqrt(texp)))
s.spec.t['cont'] = np.ones(xrange.shape)*s.spec.y.unit
s.lines = LineList(x=x_lines,
                   xmin=x_lines-xstep*5,
                   xmax=x_lines+xstep*5,
                   y=np.zeros(x_lines.shape),
                   dy=np.zeros(x_lines.shape))

for x, logN, b in zip(x_lines, logN_lines, b_lines):
    z = x*au.nm/xem_d['Ly_a']-1
    s.cb._append_syst()
    s.cb._mod_syst('Ly_a', z, logN, b)
s.cb._update_spec()

s.spec.y = s.spec.t['model']+np.random.normal(scale=s.spec.dy)

s.save(path)

plt.plot(s.spec.x, s.spec.y)
plt.plot(s.spec.x, s.spec.dy)
plt.plot(s.spec.x, s.spec.t['model'])
plt.show()
