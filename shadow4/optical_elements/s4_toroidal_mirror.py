from shadow4.syned.shape import Toroidal
from shadow4.optical_elements.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates
from shadow4.optical_surfaces.s4_toroid import S4Toroid

class S4ToroidalMirror(S4Mirror):
    def __init__(self,
                 name="Toroidal Mirror",
                 boundary_shape=None,
                 min_radius=0.0,
                 maj_radius=0.0,
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

        S4Mirror.__init__(name, boundary_shape, Toroidal(min_radius, maj_radius),
                          f_reflec, f_refl, file_refl, refraction_index)

class S4ToroidalMirrorElement(S4MirrorElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4ToroidalMirror(),
                         coordinates if coordinates is not None else ElementCoordinates())

    def analyze_surface_shape(self, beam):
        surshape = self.get_optical_element().get_surface_shape()

        if isinstance(surshape, Toroidal):
            print(">>>>> Toroidal mirror", surshape._min_radius, surshape._maj_radius)
        else:
            raise ValueError("Surface shape is not Toroidal")

        toroid = S4Toroid()
        toroid.set_toroid_radii(surshape._maj_radius, surshape._min_radius)

        mirr, normal = toroid.apply_specular_reflection_on_beam(beam)

        return mirr, normal
