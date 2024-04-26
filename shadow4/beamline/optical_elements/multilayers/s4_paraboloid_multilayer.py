from syned.beamline.shape import Paraboloid, ParabolicCylinder, Convexity, Direction, Side
from shadow4.beamline.s4_optical_element_decorators import SurfaceCalculation, S4ParaboloidOpticalElementDecorator
from shadow4.beamline.optical_elements.multilayers.s4_multilayer import S4MultilayerElement, S4Multilayer, ElementCoordinates
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4ParaboloidMultilayer(S4Multilayer, S4ParaboloidOpticalElementDecorator):
    """
    Constructor.

    Parameters
    ----------
    name : str, optional
        The name of the multilayer.
    boundary_shape : instance of BoundaryShape, optional
        The boundary shape of the multilayer.
    surface_calculation : int, optional
        flag:
            0 = SurfaceCalculation.INTERNAL,
            1 = SurfaceCalculation.EXTERNAL.
    is_cylinder : int, optional
        flag:
            0=No (there is revolution symmetry along Y)
            1=Yes (flat surface along X or Y).
    cylinder_direction : int (as defined by Direction), optional
        NONE = -1, UPWARD = 0, DOWNWARD = 1.
    convexity : int (as defined by Convexity), optional
        NONE = -1, UPWARD = 0, DOWNWARD = 1.
    parabola_parameter : float, optional
        For surface_calculation=0, The parabola parameter PARAM (y^2 = 2 PARAM z)
    at_infinity : int, optional
        For surface_calculation=0, flag to indicate that
            0: the source is at infinity (focusing paraboloid) = Side.SOURCE.
            1: the image is at infinity (collimating paraboloid) = Side.IMAGE.
    pole_to_focus : float, optional
        For surface_calculation=0, the p or q distance (from focus to center of the optical element).
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
    instance of S4ParaboloidMultilayer.
    """
    def __init__(self,
                 name="Paraboloid Multilayer",
                 boundary_shape=None,
                 at_infinity=Side.SOURCE, # used by both internal/external
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 p_focus=0.0,             # for surface_calculation=SurfaceCalculation.INTERNAL
                 q_focus=0.0,             # for surface_calculation=SurfaceCalculation.INTERNAL
                 grazing_angle=0.0,       # for surface_calculation=SurfaceCalculation.INTERNAL
                 parabola_parameter=0.0,  # for surface_calculation=SurfaceCalculation.EXTERNAL
                 pole_to_focus=0.0,       # for surface_calculation=SurfaceCalculation.EXTERNAL
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
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
                 ):
        S4ParaboloidOpticalElementDecorator.__init__(self, surface_calculation, is_cylinder, cylinder_direction, convexity,
                                                     parabola_parameter, at_infinity, pole_to_focus, p_focus, q_focus, grazing_angle)
        S4Multilayer.__init__(self,
                              name=name,
                              boundary_shape=boundary_shape,
                              surface_shape=self.get_surface_shape_instance(),
                              f_refl=f_refl,
                              file_refl=file_refl,
                              structure=structure,
                              period=period,
                              Gamma=Gamma,
                              )

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "surface_calculation": surface_calculation,
            "is_cylinder": is_cylinder,
            "cylinder_direction": cylinder_direction,
            "convexity": convexity,
            "parabola_parameter": parabola_parameter,
            "at_infinity": at_infinity,
            "pole_to_focus": pole_to_focus,
            "p_focus": p_focus,
            "q_focus": q_focus,
            "grazing_angle": grazing_angle,
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
from shadow4.beamline.optical_elements.multilayers.s4_paraboloid_multilayer import S4ParaboloidMultilayer
optical_element = S4ParaboloidMultilayer(name='{name:s}', boundary_shape=boundary_shape,
    at_infinity={at_infinity:d}, 
    surface_calculation={surface_calculation:d},
    p_focus={p_focus:f}, q_focus={q_focus:f}, grazing_angle={grazing_angle:f}, # for internal
    parabola_parameter={parabola_parameter:f}, pole_to_focus={pole_to_focus:f}, # for external
    is_cylinder={is_cylinder:d}, cylinder_direction={cylinder_direction:d}, convexity={convexity:d},
    f_refl={f_refl:d},file_refl='{file_refl:s}', structure='{structure:s}', period={period:f}, Gamma={Gamma:f})
"""
        txt += txt_pre.format(**self.__inputs)
        return txt


class S4ParaboloidMultilayerElement(S4MultilayerElement):
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
    instance of S4ParaboloidMultilayerElement
    """
    def __init__(self,
                 optical_element : S4ParaboloidMultilayer = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4ParaboloidMultilayer(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)
        if not (isinstance(self.get_optical_element().get_surface_shape(), ParabolicCylinder) or
                isinstance(self.get_optical_element().get_surface_shape(), Paraboloid)):
            raise ValueError("Wrong Optical Element: only Paraboloid or Parabolic Cylinder shape is accepted")

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
        txt += "\nfrom shadow4.beamline.optical_elements.multilayers.s4_paraboloid_multilayer import S4ParaboloidMultilayerElement"
        txt += "\nbeamline_element = S4ParaboloidMultilayerElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        return txt

if __name__=="__main__":
    from syned.beamline.shape import Rectangle
    import numpy

    angle_radial = 88.0

    el = S4ParaboloidMultilayerElement(optical_element=S4ParaboloidMultilayer(),
                                    coordinates=ElementCoordinates(p=20000, q=1000, angle_radial=88.0, angle_azimuthal=0.0))
    print(el.to_python_code())