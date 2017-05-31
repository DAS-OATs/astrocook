import sys
if (len(sys.argv) < 2):
    print("Usage: python testCont.py <UVES FITS file path>")
    sys.exit()

#Import libraries
from astrocook import spec_1d_reader
from astrocook import spec1dContinuum
import matplotlib.pyplot as plt
import time

f = sys.argv[1]

r = Spec1dReader()
s = r.uves(f)

c = spec1dContinuum()
start = time.time()
cont = c.findStretchable(s)
#cont = c.findStretchable(s, stiff=0.8)
end = time.time()
print("Time ", end - start, " s")

plt.title('File: ' + f)
plt.plot(s.x, s.y, marker='.', linestyle='None', markersize=1)
plt.plot(s.x, cont, color='red', linewidth=5)
plt.show()
