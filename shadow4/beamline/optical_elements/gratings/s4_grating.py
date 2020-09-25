import numpy

from syned.beamline.optical_elements.gratings.grating import GratingVLS
from syned.beamline.shape import Plane

from shadow4.syned.element_coordinates import ElementCoordinates # TODO from shadow4.syned.element_coordinates
from shadow4.beamline.s4_optical_element import S4OpticalElement
from shadow4.beamline.s4_beamline_element import S4BeamlineElement

class S4Grating(GratingVLS, S4OpticalElement):

    def __init__(self,
                 name="Undefined",
                 surface_shape=None,
                 boundary_shape=None,
                 ruling=800e3,
                 ruling_coeff_linear=0.0,
                 ruling_coeff_quadratic=0.0,
                 ruling_coeff_cubic=0.0,
                 ruling_coeff_quartic=0.0,
                 coating=None,
                 coating_thickness=None,):

        super.__init__(name=name,
                 surface_shape=surface_shape,
                 boundary_shape=boundary_shape,
                 ruling=ruling,
                 ruling_coeff_linear=ruling_coeff_linear,
                 ruling_coeff_quadratic=ruling_coeff_quadratic,
                 ruling_coeff_cubic=ruling_coeff_cubic,
                 ruling_coeff_quartic=ruling_coeff_quartic,
                 coating=coating,
                 coating_thickness=coating_thickness,)

class S4GratingElement(S4BeamlineElement):

    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4Grating(),
                         coordinates if coordinates is not None else ElementCoordinates())


    def trace_beam(self,beam1):
        p, q, theta1, theta2, alpha = self.get_coordinates().get_positions()

        beam = beam1.duplicate()

        print("# to be implemented...")

        return beam


def test_plane_vls():


    from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian
    from shadow4.beam.beam import Beam

    #
    # source
    #
    src = SourceGaussian.initialize_from_keywords(
                 number_of_rays=10000,
                 sigmaX=1.0e-6,
                 sigmaY=0.0,
                 sigmaZ=1.0e-6,
                 sigmaXprime=0.0002,
                 sigmaZprime=0.0002,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0]
                                 )
    beam = Beam()

    beam.genSource(src)
    beam.set_photon_energy_eV(1000.0)

    print(beam.info())

    # plotxy(Beam3.initialize_from_shadow4_beam(beam),1,3,nbins=100,title="SOURCE")

    #
    # grating
    #
    grating_syned = GratingVLS(
        name = "my_grating",
        surface_shape =  Plane(), # SurfaceShape(),
        boundary_shape = None, # BoundaryShape(),
        ruling = 600000.0,
        ruling_coeff_linear = 260818.35944225,
        ruling_coeff_quadratic = 260818.35944225,
        ruling_coeff_cubic = 13648.21037618,
        ruling_coeff_quartic = 0.0,
        coating = None,
        coating_thickness = None,
    )

    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 6.0,
                                           angle_radial = 88.840655 * numpy.pi / 180,
                                           angle_radial_out= 87.588577 * numpy.pi / 180,
                                           angle_azimuthal = 0.0)



    grating_s4 = S4GratingElement(optical_element=grating_syned, coordinates=coordinates_syned)

    print(grating_s4.info())

    beam_out = grating_s4.trace_beam(beam)



if __name__ == "__main__":

    test_plane_vls()
