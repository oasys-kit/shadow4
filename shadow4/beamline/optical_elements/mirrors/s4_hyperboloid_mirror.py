import numpy
from shadow4.syned.shape import Hyperboloid, HyperbolicCylinder, Convexity, Direction
from shadow4.beamline.s4_optical_element import SurfaceCalculation, S4HyperboloidOpticalElement
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates

class S4HyperboloidMirror(S4Mirror, S4HyperboloidOpticalElement):
    def __init__(self,
                 name="Hyperboloid Mirror",
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
        S4HyperboloidOpticalElement.__init__(self, surface_calculation, is_cylinder, cylinder_direction, convexity,
                                             min_axis, maj_axis, p_focus, q_focus, grazing_angle)
        S4Mirror.__init__(self, name, boundary_shape, self.get_surface_shape_instance(),
                          f_reflec, f_refl, file_refl, refraction_index)

    def apply_geometrical_model(self, beam):
        ccc = self.get_optical_surface_instance()
        mirr, normal = ccc.apply_specular_reflection_on_beam(beam)
        return mirr, normal

class S4HyperboloidMirrorElement(S4MirrorElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4HyperboloidMirror(),
                         coordinates if coordinates is not None else ElementCoordinates())
        if not (isinstance(self.get_optical_element().get_surface_shape(), HyperbolicCylinder) or
                isinstance(self.get_optical_element().get_surface_shape(), Hyperboloid)):
            raise ValueError("Wrong Optical Element: only Hyperboloid or Hyperbolic Cylinder shape is accepted")


if __name__=="__main__":
    from shadow4.syned.shape import Rectangle

    angle_radial = 88.0

    el = S4HyperboloidMirrorElement(optical_element=S4HyperboloidMirror(boundary_shape=Rectangle(),
                                                                        surface_calculation=SurfaceCalculation.INTERNAL,
                                                                        is_cylinder=True,
                                                                        cylinder_direction=Direction.TANGENTIAL,
                                                                        convexity=Convexity.UPWARD,
                                                                        p_focus=20000,
                                                                        q_focus=1000,
                                                                        grazing_angle=numpy.radians(90-angle_radial)),
                                    coordinates=ElementCoordinates(p=20000, q=1000, angle_radial=88.0, angle_azimuthal=0.0))

    print(el.get_optical_element().get_surface_shape())
