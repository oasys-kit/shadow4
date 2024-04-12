import numpy
from syned.beamline.shape import Conic
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates

from shadow4.beamline.s4_optical_element_decorators import S4ConicOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4ConicMirror(S4Mirror, S4ConicOpticalElementDecorator):
    """
    Constructor.

    Parameters
    ----------
    name : str, optional
        The name of the mirror.
    boundary_shape : instance of BoundaryShape, optional
        The boundary shape of the mirror.
    conic_coefficients : list, ndarray, optional
        The list of the 10 conic coefficients.
    f_reflec : int, optional
         the reflectivity of surface:
            - 0=no reflectivity,
            - 1=full polarization.
    f_refl : int, optional
        A flag to indicate the source of reflectivities:
            - 0=prerefl file
            - 1=electric susceptibility
            - 2=user defined file (1D angle in mrad, reflectivity)
            - 3=user defined file (1D energy in eV, reflectivity)
            - 4=user defined file (2D energy in eV, angle in mrad, reflectivity)
    file_refl : str, optional
            name of user defined file (for f_refl=0).
    refraction_index : complex, optional
            complex scalar with refraction index n (for f_refl=1).
    material : str, optional
            string with material formula (for f_refl=5,6)
    density : float, optional
            material density in g/cm^3 (for f_refl=5,6)

    Returns
    -------
    instance of S4ConicMirror.
    """
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
        """
        Creates the python code for defining the element.

        Parameters
        ----------
        **kwargs

        Returns
        -------
        str
            Python code.
        """
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
    """
    Constructor.

    Parameters
    ----------
    optical_element : instance of OpticalElement, optional
        The syned optical element.
    coordinates : instance of ElementCoordinates, optional
        The syned element coordinates.
    movements : instance of S4BeamlineElementMovements, optional
        The S4 element movements.
    input_beam : instance of S4Beam, optional
        The S4 incident beam.

    Returns
    -------
    instance of S4ConicMirrorElement
    """
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
        """
        Creates the python code for defining the element.

        Parameters
        ----------
        **kwargs

        Returns
        -------
        str
            Python code.
        """
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        txt += self.to_python_code_coordinates()
        txt += self.to_python_code_movements()
        txt += "\nfrom shadow4.beamline.optical_elements.mirrors.s4_conic_mirror import S4ConicMirrorElement"
        txt += "\nbeamline_element = S4ConicMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    print(S4ConicMirrorElement().to_python_code())
