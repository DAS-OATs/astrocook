from astrocook import Frame
from copy import deepcopy as dc
import numpy as np
import unittest

x = [1, 2, 3]
xmin = [0.9, 1.8, 2.7]
xmax = [1.1, 2.2, 3.3]
y = [5, 6, 7]
dy = [0.3, 0.4, 0.5]
mod1 = [2, 3, 4]
mod2 = dc(mod1)
mod2[0] = 3
obj_empty = Frame()
obj_full = Frame(x, xmin, xmax, y, dy)
check_empty = []


class FrameTest(unittest.TestCase):

    def equal(self, attr, check_full):
        self.assertListEqual(list(getattr(obj_empty, attr)), check_empty)
        self.assertListEqual(list(getattr(obj_full, attr)), check_full)

    def test_read(self):
        self.equal('x', x)
        self.equal('xmin', xmin)
        self.equal('xmax', xmax)
        self.equal('y', y)
        self.equal('dy', dy)

    def test_write(self):
        obj_full.x = mod1
        self.equal('x', mod1)
        obj_full.x = mod2
        self.equal('x', mod2)

        obj_full.xmin = mod1
        self.equal('xmin', mod1)
        obj_full.xmin = mod2
        self.equal('xmin', mod2)

        obj_full.xmax = mod1
        self.equal('xmax', mod1)
        obj_full.xmax = mod2
        self.equal('xmax', mod2)

        obj_full.y = mod1
        self.equal('y', mod1)
        obj_full.y = mod2
        self.equal('y', mod2)

        obj_full.dy = mod1
        self.equal('dy', mod1)
        obj_full.dy = mod2
        self.equal('dy', mod2)

if __name__ == '__main__':
    unittest.main()
