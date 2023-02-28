
import numpy

from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
from shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterface, S4ConicInterfaceElement


from shadow4.tools.graphics import plotxy

from syned.beamline.element_coordinates import ElementCoordinates


def get_sigmas_radiation(photon_energy,undulator_length):
    import scipy.constants as codata
    lambdan = 1e-10 * codata.h*codata.c/codata.e*1e10 / photon_energy # in m
    print("wavelength in m",lambdan)
    return 1e6*2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),1e6*0.69*numpy.sqrt(lambdan/undulator_length)


def refractive_interface_with_collimated_beam(do_plot=True):

    #
    # collimated source
    #

    src = SourceGeometrical()
    src.set_energy_distribution_singleline(value=5000, unit='A')
    src.set_spatial_type_rectangle(width=1e-3, height=1e-3)
    src.set_angular_distribution_uniform(0,0,0,0)

    beam = src.get_beam()

    print(beam.info())
    SX, SZ = (1e6*beam.get_standard_deviation(1),1e6*beam.get_standard_deviation(3))


    #
    # lens definition
    #


    interface1 = S4ConicInterfaceElement(
                                optical_element=S4ConicInterface(
                                    name="Conic Refractive Interface",
                                    boundary_shape=None,
                                    material_object="vacuum",
                                    material_image="glass",
                                    f_r_ind = 0,
                                    r_ind_obj = 1.0,
                                    r_ind_ima = 1.5,
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
        plotxy(mirr2, 1, 3, nbins=100, title="LENS HEIGHT")
        # plotxy(mirr2, 4, 5, nbins=100, title="FOOT DIV")

    FX, FZ = (1e6*beam2.get_standard_deviation(1),1e6*beam2.get_standard_deviation(3))
    print("Source dimensions: %f %f um"%(SX,SZ))
    print("Focal dimensions: %f %f um"%(FX,FZ))
    print("Demagnification: %g %g"%(SX/FX,SX/FZ))




if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()

    refractive_interface_with_collimated_beam(do_plot=True)


