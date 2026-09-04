"""
The s4 sphere multilayer (optical element and beamline element).
"""
from syned.beamline.shape import Sphere, SphericalCylinder, Convexity, Direction
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_optical_element_decorators import SurfaceCalculation, S4SphereOpticalElementDecorator
from shadow4.beamline.optical_elements.multilayers.s4_multilayer import S4MultilayerElement, S4Multilayer, ElementCoordinates
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4SphereMultilayer(S4Multilayer, S4SphereOpticalElementDecorator):
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
    radius : float, optional
        The radius of the sphere.
    p_focus : float, optional
        For surface_calculation=1, the p distance (from source plane to center of the optical element).
    q_focus : float, optional
        For surface_calculation=1, the q distance (from the center of the optical element to the image plane).
    grazing_angle : float, optional
        For surface_calculation=1, the grazing angle in rad.
    f_refl : int, optional
        A flag to indicate the source of reflectivities:
            - 0=prerefl (pre_mlayer) file
            - 1=user defined file (1D reflectivity vs angle, in mrad)
            - 2=user defined file (1D reflectivity vs energy, in eV)
            - 3=user defined file (2D reflectivity vs energy in eV and angle in mrad)
            - 4=direct calculation using xraylib
            - 5=direct calculation using dabax
    file_refl : str, optional
            name of user defined file (for f_refl=0,1,2,3).
    structure : str, optional
            A compact string defining the odd material, even material, number of bilayers and
            substrate material, in the form "[Odd,Even]xNpairs+Substrate" (e.g. "[B,W]x50+Si");
            used for f_refl=4,5.
    period : float, optional
            The bilayer thickness in Angstroms; used for f_refl=4,5.
    Gamma : float, optional
            The gamma factor thickness(even) / (thickness(odd) + thickness(even)); used for f_refl=4,5.
    density_O : None or float, optional
            The density in g/cm3 for the ODD layers; used for f_refl=4,5. If None, it is calculated
            from the material formula (required if the ODD material in `structure` is a compound).
    roughness_O : float, optional
            The roughness RMS in Angstroms for the ODD layers; used for f_refl=4,5.
    density_E : None or float, optional
            The density in g/cm3 for the EVEN layers; used for f_refl=4,5. If None, it is calculated
            from the material formula (required if the EVEN material in `structure` is a compound).
    roughness_E : float, optional
            The roughness RMS in Angstroms for the EVEN layers; used for f_refl=4,5.
    density_S : None or float, optional
            The density in g/cm3 for the SUBSTRATE; used for f_refl=4,5. If None, it is calculated
            from the material formula (required if the SUBSTRATE material in `structure` is a compound).
    roughness_S : float, optional
            The roughness RMS in Angstroms for the SUBSTRATE; used for f_refl=4,5.
    dabax : None or instance of DabaxXraylib,
        The pointer to the dabax library  (used for f_refl=5).

    Returns
    -------
    instance of S4SphereMultilayer.
    """
    def __init__(self,
                 name="Sphere Multilayer",
                 boundary_shape=None,
                 surface_calculation=SurfaceCalculation.EXTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 radius=1.0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 # inputs related to multilayer reflectivity
                 f_refl=0,   # 0=pre_mlayer file
                             # 1=user defined file (1D reflectivity vs angle)
                             # 2=user defined file (1D reflectivity vs energy)
                             # 3=user defined file (2D reflectivity vs energy and angle)
                             # 4=direct calculation using xraylib
                             # 5=direct calculation using dabax
                 file_refl="",  # preprocessor file fir f_refl=0
                 structure='[B,W]x50+Si',
                 period=25.0,
                 Gamma=0.5,
                 density_O=None,
                 roughness_O=0.0,
                 density_E=None,
                 roughness_E=0.0,
                 density_S=None,
                 roughness_S=0.0,
                 dabax=None,
                 ):
        S4SphereOpticalElementDecorator.__init__(self, surface_calculation, is_cylinder, cylinder_direction, convexity,
                                                 radius, p_focus, q_focus, grazing_angle)
        S4Multilayer.__init__(self,
                              name=name,
                              boundary_shape=boundary_shape,
                              surface_shape=self.get_surface_shape_instance(),
                              f_refl=f_refl,
                              file_refl=file_refl,
                              structure=structure,
                              period=period,
                              Gamma=Gamma,
                              density_O=density_O,
                              roughness_O=roughness_O,
                              density_E=density_E,
                              roughness_E=roughness_E,
                              density_S=density_S,
                              roughness_S=roughness_S,
                              dabax=dabax,
                              )

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "surface_calculation": surface_calculation,
            "is_cylinder": is_cylinder,
            "cylinder_direction": cylinder_direction,
            "convexity": convexity,
            "radius": radius,
            "p_focus": p_focus,
            "q_focus": q_focus,
            "grazing_angle": grazing_angle,
            "f_refl": f_refl,
            "file_refl": file_refl,
            "structure": structure,
            "period": period,
            "Gamma": Gamma,
            "density_O": density_O,
            "roughness_O": roughness_O,
            "density_E": density_E,
            "roughness_E": roughness_E,
            "density_S": density_S,
            "roughness_S": roughness_S,
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
        
from shadow4.beamline.optical_elements.multilayers.s4_sphere_multilayer import S4SphereMultilayer
optical_element = S4SphereMultilayer(name='{name:s}', boundary_shape=boundary_shape,
    surface_calculation={surface_calculation:d}, # 0=Internal calculation, 1=External 
    is_cylinder={is_cylinder:d},
    cylinder_direction={cylinder_direction:d},
    convexity={convexity:d},
    radius={radius:f}, # for surface_calculation=1
    p_focus={p_focus:f}, q_focus={q_focus:f}, grazing_angle={grazing_angle:.10f}, # for surface_calculation=0
    f_refl={f_refl:d}, # 0=pre_mlayer, 1=(mrad, refl), 2=(eV, refl), 3=(eV, mrad, refl); 4=xraylib, 5=dabax
    file_refl='{file_refl:s}', # for f_refl=0,1,2,3
    structure='{structure:s}', period={period:f}, Gamma={Gamma:f}, # for f_refl=4,5
    density_O={density_O}, roughness_O={roughness_O:f}, density_E={density_E}, roughness_E={roughness_E:f},
    density_S={density_S}, roughness_S={roughness_S:f}, # for f_refl=4,5
    dabax={dabax:s}, # if using dabax (f_refl=5), instance of DabaxXraylib() (use None for default)
    )
"""
        txt += txt_pre.format(**self.__inputs)
        return txt

class S4SphereMultilayerElement(S4MultilayerElement):
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
    instance of S4SphereMultilayerElement
    """
    def __init__(self,
                 optical_element: S4SphereMultilayer = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4SphereMultilayer(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)
        if not (isinstance(self.get_optical_element().get_surface_shape(), SphericalCylinder) or
                isinstance(self.get_optical_element().get_surface_shape(), Sphere)):
            raise ValueError("Wrong Optical Element: only Sphere or Spherical Cylinder shape is accepted")

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
        txt += "\nfrom shadow4.beamline.optical_elements.multilayers.s4_sphere_multilayer import S4SphereMultilayerElement"
        txt += "\nbeamline_element = S4SphereMultilayerElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        return txt
