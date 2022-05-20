
import numpy

from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian
from shadow4.beamline.optical_elements.additive_elements.s4_additive import S4Additive, S4AdditiveElement
from shadow4.beam.beam import Beam

from shadow4.beamline.optical_elements.ideal_elements.s4_ideal_lens import S4IdealLens, S4IdealLensElement


from shadow4.tools.graphics import plotxy

from shadow4.syned.element_coordinates import ElementCoordinates

from shadow4.beamline.optical_elements.mirrors.s4_toroidal_mirror import S4ToroidalMirror, S4ToroidalMirrorElement
from shadow4.beamline.s4_optical_element import SurfaceCalculation
from shadow4.beamline.optical_elements.mirrors.s4_surface_data_mirror import S4SurfaceDataMirror, S4SurfaceDataMirrorElement

from shadow4.syned.shape import Rectangle, Direction, Side  # TODO from syned.beamline.shape


def get_sigmas_radiation(photon_energy,undulator_length):
    import scipy.constants as codata
    lambdan = 1e-10 * codata.h*codata.c/codata.e*1e10 / photon_energy # in m
    print("wavelength in m",lambdan)
    return 1e6*2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),1e6*0.69*numpy.sqrt(lambdan/undulator_length)


# def lens_with_collimated_beam(do_plot=True):
#
#     #
#     # collimated source
#     #
#     src = SourceGaussian.initialize_collimated_source(number_of_rays=10000,sigmaX=1e-6,sigmaZ=1e-6)
#
#
#
#
#     beam = src.get_beam()
#
#     print(beam.info())
#     SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))
#
#     if do_plot:
#         plotxy(beam,1,3,nbins=100,title="SOURCE")
#         # histo1(beam, 1, nbins=100)
#
#     #
#     # lens definition
#     #
#
#
#
#     lens1e = S4AdditiveElement(optical_element=S4Additive(
#                                 name="Undefined",
#                                 optical_elements_list=[S4IdealLens(name="Undefined",focal_x=10.0, focal_y=10.0)]
#                                 ),
#                                 coordinates=ElementCoordinates(p=100.0, q=10.0))
#
#
#     print(lens1e.info())
#
#     #
#     # trace
#     #
#
#     beam2, tmp = lens1e.trace_beam(beam)
#
#     #
#     if do_plot:
#         plotxy(beam2,1,3,nbins=100,title="FOCAL PLANE")
#
#     FX, FZ = (1e6*beam2.get_standard_deviation(1),1e6*beam2.get_standard_deviation(3))
#     print("Source dimensions: %f %f um"%(SX,SZ))
#     print("Focal dimensions: %f %f um"%(FX,FZ))
#     print("Demagnification: %g %g"%(SX/FX,SX/FZ))

def example_branch_2(do_plot=True):
    source = SourceGaussian.initialize_from_keywords(number_of_rays=10000,
                 sigmaX=0.0,
                 sigmaY=0.0,
                 sigmaZ=0.0,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,)
    beam0 = Beam()
    beam0.genSource(source)
    print(beam0.info())

    if do_plot:
        plotxy(beam0, 4, 6, title="Image 0", nbins=201)

    #
    # syned definitopns
    #

    # boundaries
    boundary_shape = None #Rectangle(x_left=-rwidx2,x_right=rwidx1,y_bottom=-rlen2,y_top=rlen1)

    #
    # shadow definitions
    #
    # boundaries
    rlen1 = 0.6
    rlen2 = 0.6
    rwidx1 = 0.05
    rwidx2 = 0.05




    error_element = S4SurfaceDataMirror(name="M1",
                                          surface_data_file="../../../../oasys_workspaces/test_shadow4.hdf5",
                                          boundary_shape=Rectangle(x_left=-rwidx2, x_right=rwidx1, y_bottom=-rlen2,
                                                                   y_top=rlen1))

    added_element = S4ToroidalMirror(name="M1",
                                                        surface_calculation=SurfaceCalculation.EXTERNAL,
                                                        min_radius=0.157068,
                                                        maj_radius=358.124803 - 0.157068,
                                                        boundary_shape=boundary_shape)

    mirror1 = S4AdditiveElement(optical_element=S4Additive(name="",optical_elements_list=
                                      [error_element, added_element]
                                      ),
                                      coordinates=ElementCoordinates(p = 10.0,
                                                                     q = 6.0,
                                                                     angle_radial = numpy.radians(88.8)))


    print(mirror1.info())

    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam(beam0)
    print(mirr1.info())

    #
    # check
    #

    if do_plot:
        plotxy(beam1, 1, 3, title="Image 1", nbins=101, nolost=1)
        plotxy(mirr1, 2, 1, title="Footprint 1", nbins=101, nolost=1)

if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()
    do_plot = True

    # lens_with_collimated_beam(do_plot=do_plot)
    example_branch_2(do_plot=1)
