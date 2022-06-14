from shadow4.syned.shape import Paraboloid, ParabolicCylinder, Convexity, Direction, Side
from shadow4.beamline.s4_optical_element import SurfaceCalculation, S4ParaboloidOpticalElement
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates

class S4ParaboloidMirror(S4Mirror, S4ParaboloidOpticalElement):
    def __init__(self,
                 name="Paraboloid Mirror",
                 boundary_shape=None,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 parabola_parameter=0.0,
                 at_infinity=Side.SOURCE,
                 pole_to_focus=None,
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
        S4ParaboloidOpticalElement.__init__(self, surface_calculation, is_cylinder, cylinder_direction, convexity,
                                            parabola_parameter, at_infinity, pole_to_focus, p_focus, q_focus, grazing_angle)
        S4Mirror.__init__(self, name, boundary_shape, self.get_surface_shape_instance(),
                          f_reflec, f_refl, file_refl, refraction_index)

    def apply_geometrical_model(self, beam):
        ccc = self.get_optical_surface_instance()
        mirr, normal = ccc.apply_specular_reflection_on_beam(beam)
        return mirr, normal

class S4ParaboloidMirrorElement(S4MirrorElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4ParaboloidMirror(),
                         coordinates if coordinates is not None else ElementCoordinates())
        if not (isinstance(self.get_optical_element().get_surface_shape(), ParabolicCylinder) or
                isinstance(self.get_optical_element().get_surface_shape(), Paraboloid)):
            raise ValueError("Wrong Optical Element: only Paraboloid or Parabolic Cylinder shape is accepted")

