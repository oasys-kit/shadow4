import numpy

from shadow4.syned.shape import Sphere, SphericalCylinder, Convexity, Direction

from shadow4.beamline.s4_optical_element import SurfaceCalculation, S4SphereOpticalElement
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates
from shadow4.optical_surfaces.s4_conic import S4Conic

class S4SphereMirror(S4Mirror, S4SphereOpticalElement):
    def __init__(self,
                 name="Sphere Mirror",
                 boundary_shape=None,
                 surface_calculation=SurfaceCalculation.EXTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 radius=1.0,
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
        S4SphereOpticalElement.__init__(self, surface_calculation, is_cylinder, cylinder_direction, convexity,
                                        radius, p_focus, q_focus, grazing_angle)
        S4Mirror.__init__(self, name, boundary_shape, self.get_surface_shape_instance(),
                          f_reflec, f_refl, file_refl, refraction_index)

    def get_optical_surface_instance(self):
        surface_shape = self.get_surface_shape()

        switch_convexity = 0 if surface_shape.get_convexity() == Convexity.DOWNWARD else 1

        if isinstance(surface_shape, SphericalCylinder):
            print(">>>>> SphericalCylinder mirror", surface_shape)
            cylindrical = 1
            cylangle = 0.0 if surface_shape.get_cylinder_direction() == Direction.TANGENTIAL else (0.5 * numpy.pi)
        elif isinstance(surface_shape, Sphere):
            print(">>>>> Sphere mirror", surface_shape)
            cylindrical = 0
            cylangle    = 0.0

        return S4Conic.initialize_as_sphere_from_curvature_radius(surface_shape.get_radius(),
                                                                 cylindrical=cylindrical, cylangle=cylangle, switch_convexity=switch_convexity)


    def apply_geometrical_model(self, beam):
        ccc = self.get_optical_surface_instance()
        print(">>>>>>> ccc: ", ccc.info())
        mirr, normal = ccc.apply_specular_reflection_on_beam(beam)
        return mirr, normal

class S4SphereMirrorElement(S4MirrorElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4SphereMirror(),
                         coordinates if coordinates is not None else ElementCoordinates())
        if not (isinstance(self.get_optical_element().get_surface_shape(), SphericalCylinder) or
                isinstance(self.get_optical_element().get_surface_shape(), Sphere)):
            raise ValueError("Wrong Optical Element: only Sphere or Spherical Cylinder shape is accepted")

    def apply_local_reflection(self, beam):
        return self.get_optical_element().apply_geometrical_model(beam)

