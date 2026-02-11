"""
The s4 plane multilayer (optical element and beamline element).
"""
from syned.beamline.shape import Plane, Rectangle, Ellipse
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beam.s4_beam import S4Beam

from shadow4.beamline.s4_optical_element_decorators import S4PlaneOpticalElementDecorator
from shadow4.beamline.optical_elements.multilayers.s4_multilayer import S4MultilayerElement, S4Multilayer
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4PlaneMultilayer(S4Multilayer, S4PlaneOpticalElementDecorator):
    """
    Constructor.

    Parameters
    ----------
    name : str, optional
        The name of the mirror.
    boundary_shape : instance of BoundaryShape, optional
        The boundary shape of the mirror.
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
    dabax : None or instance of DabaxXraylib,
        The pointer to the dabax library  (used for f_refl=6).

    Returns
    -------
    instance of S4PlaneMirror.
    """
    def __init__(self,
                 name="Plane Multilayer",
                 boundary_shape=None,
                 # inputs related to multilayer reflectivity
                 f_refl=0,   # 0=pre_mlayer file
                             # 1=user defined file (1D reflectivity vs angle)
                             # 2=user defined file (1D reflectivity vs energy)
                             # 3=user defined file (2D reflectivity vs energy and angle)
                             # 4=direct calculation using xraylib
                             # 5=direct calculation using dabax
                 file_refl="",  # preprocessor file fir f_refl=0
                 structure='[B/W]x50+Si',
                 period=25.0,
                 Gamma=0.5,
                 dabax=None,
                 ):
        S4PlaneOpticalElementDecorator.__init__(self)
        S4Multilayer.__init__(self,
                              name=name,
                              boundary_shape=boundary_shape,
                              surface_shape=self.get_surface_shape_instance(),
                              f_refl=f_refl,
                              file_refl=file_refl,
                              structure=structure,
                              period=period,
                              Gamma=Gamma,
                              dabax=dabax,
                              )

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "f_refl": f_refl,
            "file_refl": file_refl,
            "structure": structure,
            "period": period,
            "Gamma": Gamma,
            "dabax": self._get_dabax_txt(),
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
   
from shadow4.beamline.optical_elements.multilayers.s4_plane_multilayer import S4PlaneMultilayer
optical_element = S4PlaneMultilayer(name='{name:s}', boundary_shape=boundary_shape,
    f_refl={f_refl:d}, # 0=pre_mlayer, 1=(mrad, refl), 2=(eV, refl), 3=(eV, mrad, refl); 4=xraylib, 5=dabax
    file_refl='{file_refl:s}', # for f_refl=0,1,2,3
    structure='{structure:s}', period={period:f}, Gamma={Gamma:f}, # for f_refl=4,5
    dabax={dabax:s}, # if using dabax (f_refl=5), instance of DabaxXraylib() (use None for default)
    )
"""
        txt += txt_pre.format(**self.__inputs)
        return txt


class S4PlaneMultilayerElement(S4MultilayerElement):
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
    instance of S4PlaneMirrorElement
    """
    def __init__(self,
                 optical_element: S4PlaneMultilayer = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4PlaneMultilayer(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)
        if not isinstance(self.get_optical_element().get_surface_shape(), Plane):
            raise ValueError("Wrong Optical Element: only Plane shape is accepted")

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
        txt += "\nfrom shadow4.beamline.optical_elements.multilayers.s4_plane_multilayer import S4PlaneMultilayerElement"
        txt += "\nbeamline_element = S4PlaneMultilayerElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    m = S4PlaneMultilayer(boundary_shape=Ellipse())
    me = S4PlaneMultilayerElement(optical_element=m, coordinates=ElementCoordinates(p=10, q=20, angle_radial=30, angle_azimuthal=40))
    print(me.info())
    print(me.to_python_code())
    print(me.duplicate().to_python_code())

