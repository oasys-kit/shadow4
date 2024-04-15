from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.crystals.s4_crystal import S4CrystalElement, S4Crystal
from shadow4.beamline.s4_optical_element_decorators import SurfaceCalculation, S4ToroidOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

from syned.beamline.shape import Toroid, Convexity, Direction

class S4ToroidCrystal(S4Crystal, S4ToroidOpticalElementDecorator):
    """
    Shadow4 Toroid Crystal Class
    This is a toroid curved perfect crystal in reflection geometry (Bragg), using the diffracted beam.

    Constructor.

    Parameters
    ----------
    name :  str, optional
        A name for the crystal
    boundary_shape : instance of BoundaryShape, optional
        The information on the crystal boundaries.
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
    material : str, optional
        The crystal material name (a name accepted by crystalpy).
    miller_index_h : int, optional
        The Miller index H.
    miller_index_k : int, optional
        The Miller index K.
    miller_index_l : int, optional
        The Miller index L.
    f_bragg_a : int, optional
        Asymmetric crystal 0:No, 1:Yes.
    asymmetry_angle : float, optional
        For f_bragg_a=1, the asymmetry angle (angle between crystal planes and surface) in rads.
    is_thick : int, optional
        Use thick crystal approximation.
    thickness : float, optional
        For is_thick=0, the crystal thickness in m.
    f_central : int, optional
        Flag for autosetting the crystal to the corrected Bragg angle.
    f_phot_cent : int, optional
        0: setting photon energy in eV, 1:setting photon wavelength in m.
    phot_cent : float, optional
        for f_central=1, the value of the photon energy (f_phot_cent=0) or photon wavelength (f_phot_cent=1).
    f_ext : inf, optional
        Flag for autosetting the crystal surface parameters.
        0: internal/calculated parameters, 1:external/user defined parameters. TODO: delete?
    material_constants_library_flag : int, optional
        Flag for indicating the origin of the crystal data:
        0: xraylib, 1: dabax, 2: preprocessor file v1, 3: preprocessor file v2.
    file_refl : str, optional
        for material_constants_library_flag=2,3, the name of the file containing the crystal parameters.

    Returns
    -------
    instance of S4ToroidCrystal.
    """
    def __init__(self,
                 name="Toroid crystal",
                 boundary_shape=None,
                 material=None,
                 # diffraction_geometry=DiffractionGeometry.BRAGG,  # ?? not supposed to be in syned...
                 miller_index_h=1,
                 miller_index_k=1,
                 miller_index_l=1,
                 asymmetry_angle=0.0,
                 is_thick=0,
                 thickness=0.010,
                 f_central=False,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 file_refl="",
                 f_bragg_a=False,
                 f_ext=0,
                 material_constants_library_flag=0,  # 0=xraylib, 1=dabax
                                                     # 2=shadow preprocessor file v1
                                                     # 3=shadow preprocessor file v2
                 min_radius=0.1,
                 maj_radius=1.0,
                 f_torus=0,
                 ):
        p_focus, q_focus, grazing_angle = 1.0, 1.0, 1e-3
        S4ToroidOpticalElementDecorator.__init__(self, SurfaceCalculation.EXTERNAL,
                                                 min_radius, maj_radius, f_torus, p_focus, q_focus, grazing_angle)

        S4Crystal.__init__(self,
                           name=name,
                           boundary_shape=boundary_shape,
                           surface_shape=self.get_surface_shape_instance(),
                           material=material,
                           # diffraction_geometry=diffraction_geometry,  # ?? not supposed to be in syned...
                           miller_index_h=miller_index_h,
                           miller_index_k=miller_index_k,
                           miller_index_l=miller_index_l,
                           asymmetry_angle=asymmetry_angle,
                           is_thick=is_thick,
                           thickness=thickness,
                           f_central=f_central,
                           f_phot_cent=f_phot_cent,
                           phot_cent=phot_cent,
                           file_refl=file_refl,
                           f_bragg_a=f_bragg_a,
                           f_ext=f_ext,
                           material_constants_library_flag=material_constants_library_flag,
                           )

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "material": material,
            # "diffraction_geometry": diffraction_geometry,
            "miller_index_h": miller_index_h,
            "miller_index_k": miller_index_k,
            "miller_index_l": miller_index_l,
            "asymmetry_angle": asymmetry_angle,
            "is_thick": is_thick,
            "thickness": thickness,
            "f_central": f_central,
            "f_phot_cent": f_phot_cent,
            "phot_cent": phot_cent,
            "file_refl": file_refl,
            "f_bragg_a": f_bragg_a,
            "f_ext": f_ext,
            "material_constants_library_flag": material_constants_library_flag,
            "min_radius": min_radius,
            "maj_radius": maj_radius,
            "f_torus": f_torus,
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
        txt = "\nfrom shadow4.beamline.optical_elements.crystals.s4_toroid_crystal import S4ToroidCrystal"

        txt_pre = """\noptical_element = S4ToroidCrystal(name='{name}',
    boundary_shape=None, material='{material}',
    miller_index_h={miller_index_h}, miller_index_k={miller_index_k}, miller_index_l={miller_index_l},
    f_bragg_a={f_bragg_a}, asymmetry_angle={asymmetry_angle},
    is_thick={is_thick}, thickness={thickness},
    f_central={f_central}, f_phot_cent={f_phot_cent}, phot_cent={phot_cent},
    file_refl='{file_refl}',
    f_ext={f_ext},
    material_constants_library_flag={material_constants_library_flag}, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
    min_radius={min_radius:g}, # min_radius = sagittal
    maj_radius={maj_radius:g}, # maj_radius = tangential
    f_torus={f_torus},
    )"""
        txt += txt_pre.format(**self.__inputs)

        return txt

class S4ToroidCrystalElement(S4CrystalElement):
    """
    The Shadow4 toroid crystal element.
    It is made of a S4ToroidCrystal and an ElementCoordinates instance. It also includes the input beam.

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

    """
    def __init__(self,
                 optical_element : S4ToroidCrystal = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4ToroidCrystal(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)

        if not (isinstance(self.get_optical_element().get_surface_shape(), Toroid)):
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
        txt += "\nfrom shadow4.beamline.optical_elements.crystals.s4_toroid_crystal import S4ToroidCrystalElement"
        txt += "\nbeamline_element = S4ToroidCrystalElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    c = S4ToroidCrystal(
            name="Undefined",
            boundary_shape=None,
            material="Si",
            # diffraction_geometry=DiffractionGeometry.BRAGG, #?? not supposed to be in syned...
            miller_index_h=1,
            miller_index_k=1,
            miller_index_l=1,
            asymmetry_angle=0.0,
            thickness=0.010,
            f_central=False,
            f_phot_cent=0,
            phot_cent=8000.0,
            file_refl="",
            f_bragg_a=False,
            f_ext=0,
            min_radius=1.0,
            maj_radius=10.0)
    # print(c.info())
    # print(c.to_python_code())


    ce = S4ToroidCrystalElement(optical_element=c)
    print(ce.info())
    print(ce.to_python_code())

    cc = S4ToroidCrystalElement()
