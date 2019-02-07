from copy import deepcopy as dc
import unittest

x = [1, 2, 3]
xmin = [0.9, 1.8, 2.7]
xmax = [1.1, 2.2, 3.3]
y = [5, 6, 7]
dy = [0.3, 0.4, 0.5]
mod1 = [2, 3, 4]
mod2 = dc(mod1)
mod2[0] = 3
check_empty = []
