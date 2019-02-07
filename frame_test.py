from astrocook import Frame
import numpy as np
from test_global import *
import unittest

obj_empty = Frame()
obj_full = Frame(x, xmin, xmax, y, dy)

class FrameTest(unittest.TestCase):

    def equal(self, attr, check_full):
        self.assertListEqual(list(getattr(obj_empty, attr)), check_empty)
        self.assertListEqual(list(getattr(obj_full, attr)), check_full)
    


class MainTest():
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

class MainTest(unittest.TestCase):
    ft = FrameTest()
    ft.test_read()

if __name__ == '__main__':
    unittest.main()
