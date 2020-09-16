
# OE surface in form of conic equation:
#      ccc[0]*X^2 + ccc[1]*Y^2 + ccc[2]*Z^2 +
#      ccc[3]*X*Y + ccc[4]*Y*Z + ccc[5]*X*Z  +
#      ccc[6]*X   + ccc[7]*Y   + ccc[8]*Z + ccc[9] = 0


from shadow4.optical_surfaces.conic import Conic

import numpy
import unittest
from numpy.testing import assert_equal, assert_almost_equal

class TestSurfaceConic(unittest.TestCase):



    def testsInitializers(self):
        #
        # initializers
        #
        a = Conic()
        print(a.info())

        a = Conic.initialize_from_coefficients([1.0,1,1,1,1,1,1,1,1,1])
        assert_equal(a.get_coefficients(),numpy.array([1.0,1,1,1,1,1,1,1,1,1]))
        print(a.info())

        a = Conic.initialize_as_plane()
        assert_equal(a.get_coefficients(),numpy.array([0,0,0,0,0,0,0,0,-1.,0]))
        print(a.info())


