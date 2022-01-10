
import numpy

from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian
from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
from shadow4.beamline.optical_elements.ideal_elements.s4_ideal_lens import S4IdealLens, S4IdealLensElement
from shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterface, S4ConicInterfaceElement


from shadow4.tools.graphics import plotxy

from shadow4.syned.element_coordinates import ElementCoordinates


def get_sigmas_radiation(photon_energy,undulator_length):
    import scipy.constants as codata
    lambdan = 1e-10 * codata.h*codata.c/codata.e*1e10 / photon_energy # in m
    print("wavelength in m",lambdan)
    return 1e6*2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),1e6*0.69*numpy.sqrt(lambdan/undulator_length)


def refractive_interface_with_collimated_beam(do_plot=True):

    #
    # collimated source
    #
    src = SourceGaussian.initialize_collimated_source(number_of_rays=10000,sigmaX=1e-6,sigmaZ=1e-6)

    src = SourceGeometrical()
    src.set_energy_distribution_singleline(value=5000, unit='A')
    src.set_spatial_type_rectangle(width=1e-3, height=1e-3)
    src.set_angular_distribution_uniform(0,0,0,0)

    beam = src.get_beam()

    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))

    # if do_plot:
    #     plotxy(beam,1,3,nbins=100,title="SOURCE")
        # histo1(beam, 1, nbins=100)

    #
    # lens definition
    #

    # interface1 = S4IdealLensElement(optical_element=S4IdealLens(name="Undefined", focal_x=5.0, focal_y=5.0),
    #                             coordinates=ElementCoordinates(p=0.0, q=5.0))
    #
    #
    # print(interface1.info())


    interface1 = S4ConicInterfaceElement(
                                optical_element=S4ConicInterface(
                                    name="Conic Refractive Interface",
                                    boundary_shape=None,
                                    material_object=1.0,
                                    material_image=1.5,
                                    conic_coefficients=[1.0, 1.0, 1.0, 0.0, -0.0, -0.0, 0.0, 0.0, 3350.0e-3, 0.0],
                                    ),
                                coordinates=ElementCoordinates(p=0.0, q=5000.0e-3,
                                            angle_radial=0.0, angle_azimuthal=0.0, angle_radial_out=numpy.pi))

    print(interface1.info())
    print(interface1.get_optical_element().get_surface_shape().get_conic_coefficients())
    #
    # trace
    #

    beam2, mirr2 = interface1.trace_beam(beam)

    #
    if do_plot:
        plotxy(beam2, 1, 3, nbins=100, title="FOCAL PLANE")
        plotxy(mirr2, 1, 3, nbins=100, title="FOOTPRINT")
        plotxy(mirr2, 4, 5, nbins=100, title="FOOT DIV")

    FX, FZ = (1e6*beam2.get_standard_deviation(1),1e6*beam2.get_standard_deviation(3))
    print("Source dimensions: %f %f um"%(SX,SZ))
    print("Focal dimensions: %f %f um"%(FX,FZ))
    print("Demagnification: %g %g"%(SX/FX,SX/FZ))



# def refractive_interface_with_divergent_beam(do_plot=True):
#
#
#
#     src = SourceGaussian.initialize_from_keywords(number_of_rays=100000,
#                                                   sigmaX=     20e-6/2.35,sigmaZ=     10e-6/2.35,
#                                                   sigmaXprime=50e-6/2.35,sigmaZprime=10e-6/2.35)
#
#     beam = src.get_beam()
#
#     SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))
#
#     if do_plot:
#         plotxy(beam,4,6,nbins=100,title="SOURCE DIVERGENCES")
#
#     beam_tmp = beam.duplicate()
#     beam_tmp.retrace(100.0)
#     if do_plot:
#         plotxy(beam_tmp,1,3,nbins=100,title="SOURCE AFTER 10m")
#
#
#     #
#     # lens definition
#     #
#
#     p = 10.0
#     q = 10.0
#     F = 1.0 / (1/p + 1/q)
#
#     lens1e = S4IdealLensElement(optical_element=S4IdealLens(name="Undefined", focal_x=F, focal_y=F),
#                                 coordinates=ElementCoordinates(p=p, q=q))
#
#     print(lens1e.info())
#     beam2, tmp = lens1e.trace_beam(beam)
#
#
#     X = beam2.get_column(1)
#     Y = beam2.get_column(3)
#     print("X: ",X.min(),X.max(),X.std())
#     print("Y: ",Y.min(),Y.max(),Y.std())
#
#
#     if do_plot:
#         plotxy(beam2,1,3,nbins=100,title="FOCAL PLANE",xrange=[-5e-5,5e-5],yrange=[-5e-5,5e-5])
#
#     FX, FZ = (1e6*beam2.get_standard_deviation(1),1e6*beam2.get_standard_deviation(3))
#     print("Source dimensions (rms): %f %f um"%(SX,SZ))
#     print("Focal dimensions (rms): %f %f um"%(FX,FZ))
#     print("Demagnification: H:%g V:%g (theoretical: %g) "%(SX/FX,SZ/FZ,p/q))


if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()
    do_plot = True

    refractive_interface_with_collimated_beam(do_plot=do_plot)
    # refractive_interface_with_divergent_beam(do_plot=do_plot)

