__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "17/01/2017"


import numpy
from numpy.testing import assert_equal, assert_almost_equal
import unittest
from shadow4.beam.beam import Beam

class TestBeam(unittest.TestCase):

    def tests_initializars(self):
        print("# ")
        print("# initializers ")
        print("# ")
        a = Beam(N=100)
        print(a.info())
        self.assertEqual(100,a.get_number_of_rays())

        a = Beam(array=numpy.zeros( (1000,18) ))
        print(a.info())
        print(a.info())
        self.assertEqual(1000,a.get_number_of_rays())

        a = Beam.initialize_from_array(numpy.zeros( (500,18) ))
        print(a.info())
        print(a.info())
        self.assertEqual(500,a.get_number_of_rays())

        a = Beam.initialize_as_pencil(200)
        print(a.info())
        self.assertEqual(200,a.get_number_of_rays())

    def test_setters_and_getters(self):

        print("# ")
        print("# setters and getters ")
        print("# ")
        a = Beam.initialize_as_pencil(200)
        b= a.duplicate()
        assert_equal (a.get_number_of_rays() - b.get_number_of_rays(), 0)
        assert_equal (a.get_column(1).mean() - b.get_column(1).mean(), 0)

        a.set_photon_energy_eV(1.0)
        assert_equal(a.get_photon_energy_eV(),1.0)

        a.set_photon_wavelength(1.51e-10)
        assert_equal(a.get_photon_wavelength(),1.51e-10)

        for i in range(18):
            a.set_column(i,numpy.pi)
            assert_equal (a.get_column(i).mean(),numpy.pi)

        a = Beam.initialize_as_pencil(200)
        assert_equal (a.get_intensity(nolost=1),200)
        assert_equal (a.get_intensity(nolost=2),0)
        flag = a.get_column(10)
        flag[50:100] = -1 # remember flag[100] is NOT changed!!
        a.set_column(10,flag)
        assert_equal (a.get_intensity(nolost=1),150)

    def test_rotations_and_translations(self):


        print("# ")
        print("# rotations and translations ")
        print("# ")
        a = Beam.initialize_as_pencil(200)
        a.translation([10,100.0,20])
        assert_equal (a.get_column(1).mean(),10)
        assert_equal (a.get_column(2).mean(),100)
        assert_equal (a.get_column(3).mean(),20)

        a = Beam.initialize_as_pencil(200)
        a.rotate(-45.*numpy.pi/180,axis=1)
        assert_equal(a.get_column(4).mean(),0)
        assert_almost_equal(a.get_column(5).mean(),numpy.sqrt(2)/2)
        assert_almost_equal(a.get_column(6).mean(),numpy.sqrt(2)/2)

        a = Beam.initialize_as_pencil(200)
        a.rotate(-45.*numpy.pi/180,axis=2)
        assert_equal(a.get_column(4).mean(),0.0)
        assert_equal(a.get_column(5).mean(),1.0)
        assert_equal(a.get_column(6).mean(),0.0)

        a = Beam.initialize_as_pencil(200)
        a.rotate(45.*numpy.pi/180,axis=3)
        # print(a.get_column(4).mean(),a.get_column(5).mean(),a.get_column(6).mean(),)
        assert_almost_equal(a.get_column(4).mean(),numpy.sqrt(2)/2)
        assert_almost_equal(a.get_column(5).mean(),numpy.sqrt(2)/2)
        assert_equal(a.get_column(6).mean(),0)

        a = Beam.initialize_as_pencil(200)
        a.rotate(-45.*numpy.pi/180,axis=1)
        a.retrace(5.0)
        assert_equal(a.get_column(1).mean(),0)
        assert_almost_equal(a.get_column(2).mean(),5.0)
        assert_almost_equal(a.get_column(3).mean(),5.0)
        #
        a = Beam.initialize_as_pencil(200)
        a.rotate(-45.*numpy.pi/180,axis=1)
        a.retrace(15.0,resetY=True)
        assert_equal(a.get_column(1).mean(),0)
        assert_almost_equal(a.get_column(2).mean(),0.0)
        assert_almost_equal(a.get_column(3).mean(),15.0)

    def test_write_and_load(self):
        print("# ")
        print("# file write/load ")
        print("# ")
        a = Beam.initialize_as_pencil(200)
        a.rays[:, 0] = numpy.random.rand(200)
        a.rays[:, 2] = numpy.random.rand(200)

        a.write("tmp.h5", simulation_name="run1", beam_name="begin", overwrite=True)
        a.write("tmp.h5", simulation_name="run1", beam_name="star01", overwrite=False)
        a.write("tmp.h5", simulation_name="run2", beam_name="begin", overwrite=False)

        b = Beam.load("tmp.h5",simulation_name="run2", beam_name="begin")

        print("a is equal to b ? ", a.identical(b))
        assert(a.identical(b))

        a.difference(b)



