
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

from syned.beamline.optical_element_with_surface_shape import OpticalElementsWithSurfaceShape
from syned.beamline.optical_elements.mirrors.mirror import Mirror as SyMirror

from syned.beamline.shape import Rectangle, Ellipse


class Mirror(object):

    def __init__(self, beamline_element_syned=None):
        if beamline_element_syned is None:
            self._beamline_element_syned = BeamlineElement(
                SyMirror(name="Undefined",
                 surface_shape=Plane(),
                 boundary_shape=None,
                 coating=None,
                 coating_thickness=None),
                ElementCoordinates())
        else:
            self._beamline_element_syned = beamline_element_syned

    def info(self):
        if self._beamline_element_syned is not None:
            return (self._beamline_element_syned.info())

    def trace_beam(self,beam_in,undo_shadow_orientation_angle_rotation=False):

        p = self._beamline_element_syned.get_coordinates().p()
        q = self._beamline_element_syned.get_coordinates().q()
        theta_grazing1 = numpy.pi / 2 - self._beamline_element_syned.get_coordinates().angle_radial()
        alpha1 = self._beamline_element_syned.get_coordinates().angle_azimuthal()

        #
        beam = beam_in.duplicate()

        #
        # put beam in mirror reference system
        #
        beam.rotate(alpha1, axis=2)
        beam.rotate(theta_grazing1, axis=1)
        beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        #
        # reflect beam in the mirror surface
        #
        soe = self._beamline_element_syned.get_optical_element()


        if not isinstance(soe, OpticalElementsWithSurfaceShape): # undefined
            print("Undefined mirror, set as plane")
            ccc = S4Conic.initialize_as_plane()
            mirr = ccc.apply_specular_reflection_on_beam(beam)
        else:
            if isinstance(self._beamline_element_syned.get_optical_element().get_surface_shape(),Plane):
                print(">>>>> Plane mirror")
                ccc = S4Conic.initialize_as_plane()
                mirr = ccc.apply_specular_reflection_on_beam(beam)
                mirr.apply_boundaries_syned(self._beamline_element_syned.get_optical_element().get_boundary_shape())
                # mirr.apply_boundaries_shadow(fhit_c=1, fshape=1, rlen1=5e-05, rlen2=5e-05, rwidx1=2e-05, rwidx2=2e-05, flag_lost_value=-1)
            elif isinstance(self._beamline_element_syned.get_optical_element().get_surface_shape(),Conic):
                print(">>>>> Conic mirror")
                conic = self._beamline_element_syned.get_optical_element().get_surface_shape()
                ccc = S4Conic.initialize_from_coefficients(conic._conic_coefficients)
                mirr = ccc.apply_specular_reflection_on_beam(beam)
                mirr.apply_boundaries_syned(self._beamline_element_syned.get_optical_element().get_boundary_shape())
            elif isinstance(self._beamline_element_syned.get_optical_element().get_surface_shape(), Toroidal):
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
        self._beamline_element_syned.get_coordinates()._angle_radial = numpy.pi / 2 - theta_grazing
        if theta_azimuthal is not None:
            self._beamline_element_syned.get_coordinates()._angle_azimuthal = numpy.pi - theta_azimuthal

    def set_surface_plane(self):
        self._beamline_element_syned._optical_element._surface_shape = Plane()

    def set_surface_conic(self, conic_coefficients=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0]):
        self._beamline_element_syned._optical_element._surface_shape = Conic(conic_coefficients=conic_coefficients)

    def set_surface_toroid(self, min_radius=0.0, maj_radius=0.0):
        self._beamline_element_syned._optical_element = Toroidal(min_radius=min_radius, maj_radius=maj_radius)

    def set_boundaries_rectangle(self, x_left=-1e3, x_right=1e3, y_bottom=-1e3, y_top=1e3):
        self._beamline_element_syned._optical_element._boundary_shape = \
            Rectangle(x_left=x_left, x_right=x_right, y_bottom=y_bottom, y_top=y_top)



if __name__ == "__main__":
    pass