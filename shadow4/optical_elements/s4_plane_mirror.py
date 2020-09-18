import numpy

from shadow4.syned.shape import Plane

from shadow4.optical_elements.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates
from shadow4.optical_surfaces.s4_conic import S4Conic

class S4PlaneMirror(S4Mirror):
    def __init__(self,
                 name="Plane Mirror",
                 boundary_shape=None,
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

        S4Mirror.__init__(name, boundary_shape, Plane(), f_reflec, f_refl, file_refl, refraction_index)

class S4PlaneMirrorElement(S4MirrorElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4PlaneMirror(),
                         coordinates if coordinates is not None else ElementCoordinates())

    def analyze_surface_shape(self, beam):
        surshape = self.get_optical_element().get_surface_shape()

        if isinstance(surshape, Plane):
            print(">>>>> Plane mirror")
        else:
            raise ValueError("Surface shape is not Plane")

        ccc = S4Conic.initialize_as_plane()

        mirr, normal = ccc.apply_specular_reflection_on_beam(beam)

        return mirr, normal
