from astrocook import Frame, Spec
from frame_test import FrameTest
from test_global import *
import unittest

obj_empty = Spec()
obj_full = Spec(x=[7,9,11], xmin=xmin, xmax=xmax, y=y, dy=dy)

class SpecTest(FrameTest):

    def test_read(self):
        self.equal('x', x)
        self.equal('xmin', xmin)
        self.equal('xmax', xmax)
        self.equal('y', y)
        self.equal('dy', dy)


if __name__ == '__main__':
    unittest.main()
