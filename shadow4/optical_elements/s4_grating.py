
import numpy
from collections import OrderedDict

from syned.beamline.optical_elements.gratings.grating import GratingVLS
from syned.beamline.shape import SurfaceShape, Plane
from shadow4.syned.element_coordinates import ElementCoordinates # from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.beamline import BeamlineElement

class Grating(object):

    def __init__(self, beamline_element_syned=None):
        if beamline_element_syned is None:
            self._beamline_element_syned = BeamlineElement(SurfaceShape(),ElementCoordinates())
        else:
            self._beamline_element_syned = beamline_element_syned

    def info(self):
        if self._beamline_element_syned is not None:
            return (self._beamline_element_syned.info())
    # def to_dictionary(self):
    #     #returns a dictionary with the variable names as keys, and a tuple with value, unit and doc string
    #     mytuple = [ ("focal_x"   ,( self._focal_x ,"m",  "Ideal lens focal length (horizontal)" ) ),
    #                 ("focal_y"   ,( self._focal_y ,"m",  "Ideal lens focal length (vertical)"  ) )]
    #     return(OrderedDict(mytuple))

    def trace_beam(self,beam1):
        p = self._beamline_element_syned.get_coordinates().p()
        q = self._beamline_element_syned.get_coordinates().q()
        theta1 = self._beamline_element_syned.get_coordinates().angle_radial()
        theta1 = self._beamline_element_syned.get_coordinates().angle_radial_2()
        alpha = self._beamline_element_syned.get_coordinates().angle_azimuthal()

        beam = beam1.duplicate()

        return beam


def test_plane_vls():


    from shadow4.sources.source_geometrical.gaussian import SourceGaussian
    from shadow4.beam.beam import Beam
    from shadow4.compatibility.beam3 import Beam3
    from Shadow.ShadowTools import plotxy

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

    beamline_element_syned = BeamlineElement(optical_element=grating_syned, coordinates=coordinates_syned)

    grating_s4 = Grating(beamline_element_syned=beamline_element_syned)

    print(grating_s4.info())

    beam_out = grating_s4.trace_beam(beam)

    #
    # lens1 = LensIdeal("test",focal_x=10.0,focal_z=10.0,p=100.0,q=10.0)
    #
    # method = 2 # 0:direct, 1:interface with overwrite, 2: no overwrite
    # if method == 0:
    #     beam2 = lens1.trace_beam(beam)
    # elif method == 1:
    #     beam.traceOE(lens1,1,overwrite=True)
    #     beam2 = beam
    # elif method == 2:
    #     beam2 = beam.traceOE(lens1,1,overwrite=True)
    # else:
    #     raise Exception("Undefined method")
    #
    # #
    #
    # plotxy(beam2.get_shadow3_beam(),1,3,nbins=100,title="FOCAL PLANE")
    # FX, FZ = (1e6*beam2.get_standard_deviation(1),1e6*beam2.get_standard_deviation(3))
    # print("Source dimensions: %f %f um"%(SX,SZ))
    # print("Focal dimensions: %f %f um"%(FX,FZ))
    # print("Demagnification: %g %g"%(SX/FX,SX/FZ))


if __name__ == "__main__":

    test_plane_vls()