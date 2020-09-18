import numpy

from shadow4.syned.shape import Ellipsoid, EllipticalCylinder, Convexity, Direction

from shadow4.optical_elements.s4_optical_element import SurfaceCalculation, S4EllipsoidOpticalElement
from shadow4.optical_elements.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates
from shadow4.optical_surfaces.s4_conic import S4Conic

class S4EllispoidMirror(S4Mirror, S4EllipsoidOpticalElement):
    def __init__(self,
                 name="Ellipsoid Mirror",
                 boundary_shape=None,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 min_axis=0.0,
                 maj_axis=0.0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 # inputs related to mirror reflectivity
                 f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                 f_refl=0,  # 0=prerefl file
                 # 1=electric susceptibility
                 # 2=user defined file (1D reflectivity vs angle)
                 # 3=user defined file (1D reflectivity vs energy)
                 # 4=user defined file (2D reflectivity vs energy and angle)
                 file_refl="",  # preprocessor file fir f_refl=0,2,3,4
                 refraction_index=1.0  # refraction index (complex) for f_refl=1
                 ):
        S4EllipsoidOpticalElement.__init__(surface_calculation, is_cylinder, cylinder_direction, convexity,
                                           min_axis, maj_axis, p_focus, q_focus, grazing_angle)
        S4Mirror.__init__(name, boundary_shape, self._curved_surface_shape,
                          f_reflec, f_refl, file_refl, refraction_index)

class S4EllipsoidMirrorElement(S4MirrorElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4EllispoidMirror(),
                         coordinates if coordinates is not None else ElementCoordinates())

    def analyze_surface_shape(self, beam):
        surshape = self.get_optical_element().get_surface_shape()

        switch_convexity = 0 if surshape.get_convexity() == Convexity.UPWARD else 1

        if isinstance(surshape, EllipticalCylinder):
            print(">>>>> EllipticalCylinder mirror", surshape)
            cylindrical = 1
            cylangle = 0.0 if surshape.get_cylinder_direction() == Direction.TANGENTIAL else (0.5 * numpy.pi)
        elif isinstance(surshape, Ellipsoid):
            print(">>>>> EllipticalCylinder mirror", surshape)
            cylindrical = 0
            cylangle    = 0.0
        else:
            raise ValueError("Surface shape is not Ellipsoid or Elliptical Cylinder")

        ccc = S4Conic.initialize_as_ellipsoid_from_focal_distances(surshape.get_p(), surshape.get_q(), surshape.get_grazing_angle(),
                                                                   cylindrical=cylindrical, cylangle=cylangle, switch_convexity=switch_convexity)

        mirr, normal = ccc.apply_specular_reflection_on_beam(beam)

        return mirr, normal
