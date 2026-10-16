from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.phase_plate_crystals.s4_phase_plate_crystal import S4PhasePlateCrystalElement, S4PhasePlateCrystal
from shadow4.beamline.s4_optical_element_decorators import S4PlaneOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements


class S4PlanePhasePlateCrystal(S4PhasePlateCrystal, S4PlaneOpticalElementDecorator):
    """
    Shadow4 Plane Phase Plate Crystal Class
    This is a plane perfect crystal plate in Bragg geometry, using the transmitted (forward-diffracted) beam.

    Constructor.

    Parameters
    ----------
    name :  str, optional
        A name for the crystal
    boundary_shape : instance of BoundaryShape, optional
        The information on the crystal boundaries.
    material : str, optional
        The crystal material name (a name accepted by crystalpy).
    miller_index_h : int, optional
        The Miller index H.
    miller_index_k : int, optional
        The Miller index K.
    miller_index_l : int, optional
        The Miller index L.
    asymmetry_angle : float, optional
        For f_bragg_a=1, the asymmetry angle (angle between crystal planes and surface) in rads.
    is_thick : int, optional
        Use thick crystal approximation. Must be 0.
    thickness : float, optional
        The crystal thickness in m.
    f_central : int, optional
        Flag for autosetting the crystal to the corrected Bragg angle (plus the detuning).
    f_phot_cent : int, optional
        0: setting photon energy in eV, 1:setting photon wavelength in A.
    phot_cent : float, optional
        for f_central=1, the value of the photon energy (f_phot_cent=0) or photon wavelength (f_phot_cent=1).
    f_detune : int, optional
        For f_central=1, detuning: 0: no detuning, 1: detune_angle, 2: detune_half_darwin_widths Darwin half widths.
    detune_angle : float, optional
        For f_detune=1, the detuning angle in rad.
    detune_half_darwin_widths : float, optional
        For f_detune=2, the detuning in units of the Darwin half width (sigma polarization).
    file_refl : str, optional
        for material_constants_library_flag=2,3, the name of the file containing the crystal parameters.
    f_bragg_a : int, optional
        Asymmetric crystal 0:No, 1:Yes.
    f_ext : inf, optional
        Flag for autosetting the crystal surface parameters.
    material_constants_library_flag : int, optional
        Flag for indicating the origin of the crystal data:
        0: xraylib, 1: dabax, 2: preprocessor file v1, 3: preprocessor file v2.
    dabax : None or instance of DabaxXraylib,
        The pointer to the dabax library  (used for material_constants_library_flag=1).

    Returns
    -------
    instance of S4PlanePhasePlateCrystal.
    """
    def __init__(self,
                 name="Plane phase plate crystal",
                 boundary_shape=None,
                 material=None,
                 miller_index_h=1,
                 miller_index_k=1,
                 miller_index_l=1,
                 asymmetry_angle=0.0,
                 is_thick=0,
                 thickness=100e-6,
                 f_central=False,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 f_detune=0,
                 detune_angle=0.0,
                 detune_half_darwin_widths=0.0,
                 file_refl="",
                 f_bragg_a=False,
                 f_ext=0,
                 material_constants_library_flag=0,  # 0=xraylib, 1=dabax
                                                     # 2=shadow preprocessor file v1
                                                     # 3=shadow preprocessor file v2
                 dabax=None,
                 ):
        S4PlaneOpticalElementDecorator.__init__(self)
        S4PhasePlateCrystal.__init__(self,
                           name=name,
                           surface_shape=self.get_surface_shape_instance(),
                           boundary_shape=boundary_shape,
                           material=material,
                           miller_index_h=miller_index_h,
                           miller_index_k=miller_index_k,
                           miller_index_l=miller_index_l,
                           asymmetry_angle=asymmetry_angle,
                           is_thick=is_thick,
                           thickness=thickness,
                           f_central=f_central,
                           f_phot_cent=f_phot_cent,
                           phot_cent=phot_cent,
                           f_detune=f_detune,
                           detune_angle=detune_angle,
                           detune_half_darwin_widths=detune_half_darwin_widths,
                           file_refl=file_refl,
                           f_bragg_a=f_bragg_a,
                           f_ext=f_ext,
                           material_constants_library_flag=material_constants_library_flag,
                           dabax=dabax,
                           )

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "material": material,
            "miller_index_h": miller_index_h,
            "miller_index_k": miller_index_k,
            "miller_index_l": miller_index_l,
            "asymmetry_angle": repr(float(asymmetry_angle)),
            "is_thick": is_thick,
            "thickness": repr(float(thickness)),
            "f_central": f_central,
            "f_phot_cent": f_phot_cent,
            "phot_cent": repr(float(phot_cent)),
            "f_detune": f_detune,
            "detune_angle": repr(float(detune_angle)),
            "detune_half_darwin_widths": repr(float(detune_half_darwin_widths)),
            "file_refl": file_refl,
            "f_bragg_a": f_bragg_a,
            "f_ext": f_ext,
            "material_constants_library_flag": material_constants_library_flag,
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

from shadow4.beamline.optical_elements.phase_plate_crystals.s4_plane_phase_plate_crystal import S4PlanePhasePlateCrystal
optical_element = S4PlanePhasePlateCrystal(name='{name}',
    boundary_shape=boundary_shape, material='{material}',
    miller_index_h={miller_index_h}, miller_index_k={miller_index_k}, miller_index_l={miller_index_l},
    f_bragg_a={f_bragg_a}, asymmetry_angle={asymmetry_angle},
    is_thick={is_thick}, thickness={thickness},
    f_central={f_central}, f_phot_cent={f_phot_cent}, phot_cent={phot_cent},
    f_detune={f_detune}, detune_angle={detune_angle}, detune_half_darwin_widths={detune_half_darwin_widths}, # 0=No, 1=angle, 2=Darwin half widths
    file_refl='{file_refl}',
    f_ext={f_ext},
    material_constants_library_flag={material_constants_library_flag}, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
    dabax={dabax}, # used when material_constants_library_flag=1,
    )"""
        txt += txt_pre.format(**self.__inputs)

        return txt


class S4PlanePhasePlateCrystalElement(S4PhasePlateCrystalElement):
    """
    The Shadow4 plane phase plate crystal element.
    It is made of a S4PlanePhasePlateCrystal and an ElementCoordinates instance. It also includes the input beam.

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
                 optical_element : S4PlanePhasePlateCrystal = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4PlanePhasePlateCrystal(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)

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
        txt += "\nfrom shadow4.beamline.optical_elements.phase_plate_crystals.s4_plane_phase_plate_crystal import S4PlanePhasePlateCrystalElement"
        txt += "\nbeamline_element = S4PlanePhasePlateCrystalElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        return txt


if __name__ == "__main__":
    c = S4PlanePhasePlateCrystal(material="Si", miller_index_h=2, miller_index_k=2, miller_index_l=0,
                                 thickness=100e-6, f_central=1, phot_cent=10000.0, f_detune=2,
                                 detune_half_darwin_widths=-5.0)
    ce = S4PlanePhasePlateCrystalElement(optical_element=c)
    print(ce.info())
    print(ce.to_python_code())
