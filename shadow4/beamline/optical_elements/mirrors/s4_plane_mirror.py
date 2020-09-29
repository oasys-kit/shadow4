from shadow4.syned.shape import Plane

from shadow4.beamline.s4_optical_element import S4PlaneOpticalElement
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates
from shadow4.optical_surfaces.s4_conic import S4Conic

class S4PlaneMirror(S4Mirror, S4PlaneOpticalElement):
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
        S4PlaneOpticalElement.__init__(self)
        S4Mirror.__init__(self, name, boundary_shape, self._plane_surface_shape, f_reflec, f_refl, file_refl, refraction_index)

class S4PlaneMirrorElement(S4MirrorElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4PlaneMirror(),
                         coordinates if coordinates is not None else ElementCoordinates())
        if not isinstance(self.get_optical_element().get_surface_shape(), Plane):
            raise ValueError("Wrong Optical Element: only Plane shape is accepted")

    def apply_local_reflection(self, beam):
        print(">>>>> Plane mirror")

        ccc = S4Conic.initialize_as_plane()

        mirr, normal = ccc.apply_specular_reflection_on_beam(beam)

        return mirr, normal
