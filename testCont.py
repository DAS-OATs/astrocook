import sys
if (len(sys.argv) < 2):
    print("Usage: python testCont.py <UVES FITS file path>")
    sys.exit()

#Import libraries
<<<<<<< HEAD
from astrocook import spec_1d_reader
from astrocook import spec1dContinuum
=======
from astrocook import Spec1D
from astrocook import Spec1DReader
from astrocook import Spec1DCont
>>>>>>> f78c5b36c1d239c084fb35467ce4af5adead4cc3
import matplotlib.pyplot as plt
import time

f = sys.argv[1]

<<<<<<< HEAD
r = Spec1dReader()
=======
r = Spec1DReader()
>>>>>>> f78c5b36c1d239c084fb35467ce4af5adead4cc3
s = r.uves(f)

c = Spec1DCont(s)
start = time.time()
cont = c.findStretchable(s)
#cont = c.findStretchable(s, stiff=0.8)
end = time.time()
print("Time ", end - start, " s")

plt.title('File: ' + f)
plt.plot(s.x, s.y, marker='.', linestyle='None', markersize=1)
plt.plot(s.x, cont, color='red', linewidth=5)
plt.show()
