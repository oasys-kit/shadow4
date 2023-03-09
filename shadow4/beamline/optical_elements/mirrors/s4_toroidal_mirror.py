from syned.beamline.shape import Toroid
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates
from shadow4.beamline.s4_optical_element import SurfaceCalculation, S4ToroidalOpticalElementDecorator

class S4ToroidalMirror(S4Mirror, S4ToroidalOpticalElementDecorator):
    def __init__(self,
                 name="Toroidal Mirror",
                 boundary_shape=None,
                 surface_calculation=SurfaceCalculation.EXTERNAL,
                 min_radius=0.1,
                 maj_radius=1.0,
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
        S4ToroidalOpticalElementDecorator.__init__(self, surface_calculation,
                                                   min_radius, maj_radius, p_focus, q_focus, grazing_angle)
        S4Mirror.__init__(self, name, boundary_shape, self.get_surface_shape_instance(),
                          f_reflec, f_refl, file_refl, refraction_index)


    def apply_geometrical_model(self, beam):
        toroid = self.get_optical_surface_instance()
        footprint, normal = toroid.apply_specular_reflection_on_beam(beam)
        return footprint, normal


class S4ToroidalMirrorElement(S4MirrorElement):
    def __init__(self,
                 optical_element : S4ToroidalMirror = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4ToroidalMirror(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)
        if not isinstance(self.get_optical_element().get_surface_shape(), Toroid):
            raise ValueError("Wrong Optical Element: only Toroid shape is accepted")

