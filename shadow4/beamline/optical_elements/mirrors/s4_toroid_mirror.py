from syned.beamline.shape import Toroid
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.mirrors.s4_mirror import S4MirrorElement, S4Mirror, ElementCoordinates
from shadow4.beamline.s4_optical_element_decorators import SurfaceCalculation, S4ToroidOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4ToroidMirror(S4Mirror, S4ToroidOpticalElementDecorator):
    """
    Constructor.

    Parameters
    ----------
    name : str, optional
        The name of the mirror.
    boundary_shape : instance of BoundaryShape, optional
        The boundary shape of the mirror.
    surface_calculation : int, optional
        flag:
            0 = SurfaceCalculation.INTERNAL,
            1 = SurfaceCalculation.EXTERNAL.
    min_radius : float, optional
        The minor axis of the toroid in m. This corresponds to the sagittal optical surface.
    maj_radius : float, optional
        The optical surface major radius in m. This corresponds to the tangential optical surface
        (it is **not** the radius of the toroid).
    f_torus : int, optional
        Flag to indicate which optical surface is in use (where the mirror pole is located). Values are:
            * 0=lower/outer (tangential: concave / sagittal:concave),
            * 1=lower/inner (tangential: concave / sagittal:convex) ,
            * 2=upper/inner (tangential: convex  / sagittal:concave),
            * 3=upper/outer (tangential: convex  / sagittal:convex).
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
    instance of S4ToroidMirror.
    """
    def __init__(self,
                 name="Toroid Mirror",
                 boundary_shape=None,
                 surface_calculation=SurfaceCalculation.EXTERNAL,
                 min_radius=0.1,
                 maj_radius=1.0,
                 f_torus=0, # mirror pole location:
                            #  lower/outer (concave/concave) (0),
                            # lower/inner (concave/convex) (1),
                            # upper/inner (convex/concave) (2),
                            # upper/outer (convex/convex) (3).
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 # inputs related to mirror reflectivity
                 f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                 f_refl=0,  # 0=prerefl file
                 # 1=electric susceptibility
                 # 2=user defined file (1D reflectivity vs angle)
                 # 3=user defined file (1D reflectivity vs energy)
                 # 4=user defined file (2D reflectivity vs energy and angle)
                 file_refl="",  # preprocessor file fir f_refl=0,2,3,4
                 refraction_index=1.0,  # refraction index (complex) for f_refl=1
                 coating_material="",   # string with coating material formula for f_refl=5,6
                 coating_density=1.0,   # coating material density for f_refl=5,6
                 coating_roughness=0.0, # coating material roughness in A for f_refl=5,6
                 ):
        S4ToroidOpticalElementDecorator.__init__(self, surface_calculation,
                                                 min_radius, maj_radius, f_torus, p_focus, q_focus, grazing_angle)
        S4Mirror.__init__(self, name, boundary_shape, self.get_surface_shape_instance(),
                          f_reflec, f_refl, file_refl, refraction_index, coating_material, coating_density, coating_roughness)

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "surface_calculation": surface_calculation,
            "min_radius" : min_radius,
            "maj_radius" : maj_radius,
            "f_torus": f_torus,
            "p_focus": p_focus,
            "q_focus": q_focus,
            "grazing_angle": grazing_angle,
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
        
from shadow4.beamline.optical_elements.mirrors.s4_toroid_mirror import S4ToroidMirror
optical_element = S4ToroidMirror(name='{name:s}',boundary_shape=boundary_shape,
    surface_calculation={surface_calculation:d},
    min_radius={min_radius:.6g},  # min_radius = sagittal
    maj_radius={maj_radius:.6g},  # maj_radius = tangential
    f_torus={f_torus},
    p_focus={p_focus:.6g},q_focus={q_focus:.6g},grazing_angle={grazing_angle:g},
    f_reflec={f_reflec:d},f_refl={f_refl:d},file_refl='{file_refl:s}',refraction_index={refraction_index:g},
    coating_material='{coating_material:s}',coating_density={coating_density:g},coating_roughness={coating_roughness:g})
"""
        txt += txt_pre.format(**self.__inputs)
        return txt


class S4ToroidMirrorElement(S4MirrorElement):
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
    instance of S4ToroidMirrorElement
    """
    def __init__(self,
                 optical_element : S4ToroidMirror = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4ToroidMirror(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)
        if not isinstance(self.get_optical_element().get_surface_shape(), Toroid):
            raise ValueError("Wrong Optical Element: only Toroid shape is accepted")

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
        txt += "\nfrom shadow4.beamline.optical_elements.mirrors.s4_toroid_mirror import S4ToroidMirrorElement"
        txt += "\nbeamline_element = S4ToroidMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    a = S4ToroidMirror(refraction_index=6j)
    b= S4ToroidMirrorElement(optical_element=a, coordinates=ElementCoordinates(1.11111111111111111111111111111111,2,3,4))
    print(b.to_python_code())
