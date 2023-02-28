from shadow4.syned.shape import Plane
from shadow4.beamline.s4_optical_element import S4PlaneOpticalElement
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates


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
        S4Mirror.__init__(self, name, boundary_shape, self.get_surface_shape_instance(), f_reflec, f_refl, file_refl, refraction_index)


    def apply_geometrical_model(self, beam):
        ccc = self.get_optical_surface_instance()
        mirr, normal = ccc.apply_specular_reflection_on_beam(beam)
        return mirr, normal

    def to_python_code(self, data=None):
        txt = "\nfrom shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror"
        txt_pre = "\noptical_element = S4PlaneMirror(name='{name:s}',boundary_shape=None,f_reflec={f_reflec:d},f_refl={f_refl:d},file_refl='{file_refl:s}',refraction_index={refraction_index:g})"
        txt += txt_pre.format(name=self.get_name(),
                              f_reflec=self._f_reflec,
                              f_refl=self._f_refl,
                              file_refl=self._file_refl,
                              refraction_index=self._refraction_index)
        return txt

class S4PlaneMirrorElement(S4MirrorElement):
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4PlaneMirror(),
                         coordinates if coordinates is not None else ElementCoordinates())
        if not isinstance(self.get_optical_element().get_surface_shape(), Plane):
            raise ValueError("Wrong Optical Element: only Plane shape is accepted")

    def to_python_code(self, data=None):
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        coordinates = self.get_coordinates()
        txt += "\nfrom shadow4.syned.element_coordinates import ElementCoordinates"
        txt += "\ncoordinates=ElementCoordinates(p=%g,q=%g,angle_radial=%g)" % \
               (coordinates.p(), coordinates.q(), coordinates.angle_radial())
        txt += "\nfrom shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement"
        txt += "\nbeamline_element = S4PlaneMirrorElement(optical_element=optical_element,coordinates=coordinates)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam(beam)"
        return txt

if __name__ == "__main__":
    m = S4PlaneMirror(refraction_index=1+1e-7j)
    me = S4PlaneMirrorElement(optical_element=m, coordinates=ElementCoordinates(p=10))
    print(me.info())
    print(me.to_python_code())
