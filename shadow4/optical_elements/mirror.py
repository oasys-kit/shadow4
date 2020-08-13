
import numpy
from collections import OrderedDict


from syned.beamline.optical_elements.mirrors.mirror import Mirror as SMirror
from syned.beamline.shape import SurfaceShape, BoundaryShape
from syned.beamline.shape import Plane, Sphere, SphericalCylinder, Ellipsoid, EllipticalCylinder, Paraboloid
from syned.beamline.shape import ParabolicCylinder, Hyperboloid, HyperbolicCylinder, Toroidal, Conic

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.beamline import BeamlineElement

from shadow4.optical_surfaces.conic import Conic as S4Conic
from shadow4.optical_surfaces.toroid import Toroid as S4Toroid

class Mirror(object):

    def __init__(self, beamline_element_syned=None):
        if beamline_element_syned is None:
            self._beamline_element_syned = BeamlineElement(SurfaceShape(),ElementCoordinates())
        else:
            self._beamline_element_syned = beamline_element_syned

    def info(self):
        if self._beamline_element_syned is not None:
            return (self._beamline_element_syned.info())

    def trace_beam(self,beam1,undo_shadow_orientation_angle_rotation=False):

        p = self._beamline_element_syned.get_coordinates().p()
        q = self._beamline_element_syned.get_coordinates().q()
        theta_grazing1 = numpy.pi - self._beamline_element_syned.get_coordinates().angle_radial()
        alpha1 = self._beamline_element_syned.get_coordinates().angle_azimuthal()

        #
        beam = beam1.duplicate()

        #
        # put beam in mirror reference system
        #
        beam.rotate(alpha1, axis=2)
        beam.rotate(theta_grazing1, axis=1)
        beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        #
        # reflect beam in the mirror surface
        #

        if isinstance(self._beamline_element_syned.get_optical_element().get_surface_shape(),Plane):
            print(">>>>> Plane mirror")
            ccc = S4Conic.initialize_as_plane()
            mirr = ccc.apply_specular_reflection_on_beam(beam)
        elif isinstance(self._beamline_element_syned.get_optical_element().get_boundary_shape(),Conic):
            print(">>>>> Conic (no plane) mirror")
            conic = self._beamline_element_syned.get_optical_element().get_boundary_shape()
            ccc = S4Conic.initialize_from_coefficients(conic._conic_coefficients)
            mirr = ccc.apply_specular_reflection_on_beam(beam)
        elif isinstance(self._beamline_element_syned.get_optical_element().get_boundary_shape(), Toroidal):
            print(">>>>> Toroidal mirror")
            #...........
        else:
            raise Exception("cannot trace this surface shape")

        #
        # TODO: apply mirror boundaries...
        #

        #
        # TODO: apply mirror reflectivity...
        #


        #
        # from mirror reference system to image plane
        #

        beam_out = mirr.duplicate()
        beam_out.rotate(theta_grazing1, axis=1)
        # do not undo alpha rotation: newbeam.rotate(-alpha, axis=2)
        if undo_shadow_orientation_angle_rotation:
            beam_out.rotate(-alpha1, axis=2)
            beam_out.retrace(q, resetY=True)

        return beam_out, mirr

    #
    # i/o utilities
    #
    def set_positions(self, p, q, theta_grazing, theta_azimuthal=None):
        self._beamline_element_syned.get_coordinates()._p = p
        self._beamline_element_syned.get_coordinates()._q = q
        self._beamline_element_syned.get_coordinates()._angle_radial = numpy.pi - theta_grazing
        if theta_azimuthal is not None:
            self._beamline_element_syned.get_coordinates()._angle_azimuthal = numpy.pi - theta_azimuthal

    def set_surface(self, surface='plane'):
        if surface == 'plane':
            self._beamline_element_syned._optical_element = Plane()
        elif surface == 'conic':
            self._beamline_element_syned._optical_element = Conic()
        elif surface == 'toroid':
            self._beamline_element_syned._optical_element = Toroidal()
        else:
            raise("Not implemented")

if __name__ == "__main__":
    from shadow4.beam.beam import Beam
    from shadow4.sources.source_geometrical.gaussian import SourceGaussian

    # beam0 = Beam.initialize_as_pencil(N=500)
    source = SourceGaussian.initialize_from_keywords(500, 0, 0, 0, 1e-6, 1e-6)

    beam0 = Beam()
    beam0.genSource(source)
    print(beam0.info())

    surface_shape = Plane() # SurfaceShape()
    boundary_shape = None   # BoundaryShape()

    smirror1 = SMirror(
                name="Undefined",
                surface_shape=surface_shape,
                boundary_shape=boundary_shape,
                coating=None,
                coating_thickness=None)

    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 6.0,
                                           angle_radial = 88.840655 * numpy.pi / 180,)

    beamline_element_syned = BeamlineElement(optical_element=smirror1, coordinates=coordinates_syned)
    #

    mirror1 = Mirror(beamline_element_syned=beamline_element_syned)

    # mirror2 = Mirror()
    # mirror2.set_positions(10,100,3e-3)
    # mirror2.set_surface('conic')

    # print(mirror1.info())
    print(mirror1.info())

    beam1, mirr1 = mirror1.trace_beam(beam0)
    print(mirr1.info())
    print(beam1.info())





