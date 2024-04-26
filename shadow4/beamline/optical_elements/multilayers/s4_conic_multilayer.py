import numpy
from syned.beamline.shape import Conic
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.multilayers.s4_multilayer import S4MultilayerElement, S4Multilayer, ElementCoordinates

from shadow4.beamline.s4_optical_element_decorators import S4ConicOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4ConicMultilayer(S4Multilayer, S4ConicOpticalElementDecorator):
    """
    Constructor.

    Parameters
    ----------
    name : str, optional
        The name of the multilayer.
    boundary_shape : instance of BoundaryShape, optional
        The boundary shape of the multilayer.
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
    instance of S4ConicMultilayer.
    """
    def __init__(self,
                 name="Conic Multilayer",
                 boundary_shape=None,
                 conic_coefficients=[0.0]*10,
                 # inputs related to multilayer reflectivity
                 f_refl=0,  # 0=pre_mlayer file
                 # 1=user defined file (1D reflectivity vs angle)
                 # 2=user defined file (1D reflectivity vs energy)
                 # 3=user defined file (2D reflectivity vs energy and angle)
                 # 4=direct calculation using xraylib
                 # 5=direct calculation using dabax
                 file_refl="",  # preprocessor file fir f_refl=0
                 structure='[B/W]x50+Si',
                 period=25.0,
                 Gamma=0.5,
                 ):
        S4ConicOpticalElementDecorator.__init__(self, conic_coefficients)
        S4Multilayer.__init__(self,
                              name=name,
                              boundary_shape=boundary_shape,
                              surface_shape=self.get_surface_shape_instance(),
                              f_refl=f_refl,
                              file_refl=file_refl,
                              structure=structure,
                              period=period,
                              Gamma=Gamma,)

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "conic_coefficients": repr(conic_coefficients),
            "f_refl": f_refl,
            "file_refl": file_refl,
            "structure": structure,
            "period": period,
            "Gamma": Gamma,
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
        
from shadow4.beamline.optical_elements.multilayers.s4_conic_multilayer import S4ConicMultilayer
optical_element = S4ConicMultilayer(name='{name:s}',boundary_shape=boundary_shape,
    conic_coefficients={conic_coefficients:s},
    f_refl={f_refl:d},file_refl='{file_refl:s}', structure='{structure:s}', period={period:f}, Gamma={Gamma:f})
"""
        txt += txt_pre.format(**self.__inputs)
        return txt


class S4ConicMultilayerElement(S4MultilayerElement):
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
    instance of S4ConicMultilayerElement
    """
    def __init__(self,
                 optical_element: S4ConicMultilayer = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4ConicMultilayer(),
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
        txt += "\nfrom shadow4.beamline.optical_elements.multilayers.s4_conic_multilayer import S4ConicMultilayerElement"
        txt += "\nbeamline_element = S4ConicMultilayerElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    print(S4ConicMultilayerElement().to_python_code())
