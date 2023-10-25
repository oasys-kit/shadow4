import numpy
from syned.beamline.shape import Conic
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates

from shadow4.beamline.s4_optical_element_decorators import S4ConicOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4ConicMirror(S4Mirror, S4ConicOpticalElementDecorator):
    def __init__(self,
                 name="Conic Mirror",
                 boundary_shape=None,
                 conic_coefficients=[0.0]*10,
                 # inputs related to mirror reflectivity
                 f_reflec=0, # reflectivity of surface: 0=no reflectivity, 1=full polarization
                 f_refl=0,   # 0=prerefl file
                             # 1=electric susceptibility
                             # 2=user defined file (1D reflectivity vs angle)
                             # 3=user defined file (1D reflectivity vs energy)
                             # 4=user defined file (2D reflectivity vs energy and angle)
                             # 5=direct calculation using xraylib
                             # 6=direct calculation using dabax
                 file_refl="",  # preprocessor file fir f_refl=0,2,3,4
                 refraction_index=1.0,  # refraction index (complex) for f_refl=1
                 coating_material="",   # string with coating material formula for f_refl=5,6
                 coating_density=1.0,   # coating material density for f_refl=5,6
                 coating_roughness=0.0, # coating material roughness in A for f_refl=5,6
                 ):
        S4ConicOpticalElementDecorator.__init__(self, conic_coefficients)
        S4Mirror.__init__(self, name, boundary_shape, self.get_surface_shape_instance(),
                          f_reflec, f_refl, file_refl, refraction_index, coating_material, coating_density, coating_roughness)

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "conic_coefficients": repr(conic_coefficients),
            "f_reflec": f_reflec,
            "f_refl": f_refl,
            "file_refl": file_refl,
            "refraction_index": refraction_index,
            "coating_material": coating_material,
            "coating_density": coating_density,
            "coating_roughness": coating_roughness,
        }

    def to_python_code(self, **kwargs):
        txt = self.to_python_code_boundary_shape()
        txt_pre = """
        
from shadow4.beamline.optical_elements.mirrors.s4_conic_mirror import S4ConicMirror
optical_element = S4ConicMirror(name='{name:s}',boundary_shape=boundary_shape,
    conic_coefficients={conic_coefficients:s},
    f_reflec={f_reflec:d},f_refl={f_refl:d},file_refl='{file_refl:s}',refraction_index={refraction_index:g},
    coating_material='{coating_material:s}',coating_density={coating_density:g},coating_roughness={coating_roughness:g})
"""
        txt += txt_pre.format(**self.__inputs)
        return txt


class S4ConicMirrorElement(S4MirrorElement):
    def __init__(self,
                 optical_element: S4ConicMirror = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4ConicMirror(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)
        if not isinstance(self.get_optical_element().get_surface_shape(), Conic):
            raise ValueError("Wrong Optical Element: only Conic shape is accepted")

    def to_python_code(self, **kwargs):
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        coordinates = self.get_coordinates()
        txt += "\nfrom syned.beamline.element_coordinates import ElementCoordinates"
        txt += "\ncoordinates = ElementCoordinates(p=%.12g, q=%.12g, angle_radial=%.12g, angle_azimuthal=%.12g, angle_radial_out=%.12g)" % \
               (coordinates.p(), coordinates.q(), coordinates.angle_radial(), coordinates.angle_azimuthal(), coordinates.angle_radial_out())

        txt += self.to_python_code_movements()

        txt += "\nfrom shadow4.beamline.optical_elements.mirrors.s4_conic_mirror import S4ConicMirrorElement"
        txt += "\nbeamline_element = S4ConicMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    print(S4ConicMirrorElement().to_python_code())
