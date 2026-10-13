"""
The s4 hyperboloid mosaic crystal (optical element and beamline element).
"""
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.mosaic_crystals.s4_mosaic_crystal import S4MosaicCrystalElement, S4MosaicCrystal
from shadow4.beamline.s4_optical_element_decorators import SurfaceCalculation, S4HyperboloidOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

from syned.beamline.shape import Hyperboloid, HyperbolicCylinder, Convexity, Direction

class S4HyperboloidMosaicCrystal(S4MosaicCrystal, S4HyperboloidOpticalElementDecorator):
    """
    Shadow4 Hyperboloid Mosaic Crystal Class
    This is a hyperboloid mosaic crystal in reflection geometry (Bragg), using the diffracted beam.

    Constructor.

    Parameters
    ----------
    name :  str, optional
        A name for the crystal
    boundary_shape : instance of BoundaryShape, optional
        The information on the crystal boundaries.
    is_cylinder : int, optional
        Flag to indicate that the surface has cylindrical symmetry (it is flat in one direction).
    cylinder_direction : int, optional
       For is_cylinder=1, the direction where the surface is flat.
       Use synedDirection.TANGENTIAL (0) or Direction.SAGITTAL (1).
    convexity : int, optional
        The surface is concave (0) or convex (1).
        Use syned Convexity.UPWARD (0) for concave or Convexity.DOWNWARD (1).
    min_axis : float, optional
        The hyperbola/hyperboloid minor axis.
    maj_axis : float, optional
        The hyperbola/hyperboloid major axis.
    pole_to_focus : float, optional
        The distance from focus 1 (locus of the source) to the crystal pole.
    material : str, optional
        The crystal material name (a name accepted by crystalpy).
    miller_index_h : int, optional
        The Miller index H.
    miller_index_k : int, optional
        The Miller index K.
    miller_index_l : int, optional
        The Miller index L.
    thickness : float, optional
        The diffracting crystal thickness in m.
    f_central : int, optional
        Flag for autosetting the crystal to the corrected Bragg angle.
    f_phot_cent : int, optional
        0: setting photon energy in eV, 1:setting photon wavelength in Angstrom.
    phot_cent : float, optional
        for f_central=1, the value of the photon energy (f_phot_cent=0) or photon wavelength (f_phot_cent=1).
    mosaicity_fwhm_deg : float, optional
        The crystal mosaicity (FWHM) in degrees.
    mosaicity_profile_flag : int, optional
        The crystallite distribution flag: 0=Gaussian, 1=External.
    material_constants_library_flag : int, optional
        Flag for indicating the origin of the crystal data:
        0: xraylib, 1: dabax, 2: preprocessor file v1, 3: preprocessor file v2.
    file_refl : str, optional
        for material_constants_library_flag=2,3, the name of the file containing the crystal parameters.
    dabax : None or instance of DabaxXraylib,
        The pointer to the dabax library  (used for material_constants_library_flag=1).
    calculation_method : int, optional
        The model used for the diffraction in the mosaic crystal:
        0: macroscopic model (1992 paper), 1: Monte Carlo model (crystallite by crystallite).
        See S4MosaicCrystal.
    mc_crystallite_thickness_flag : int, optional
        For calculation_method=1, how the crystallite thickness t0 is set:
        0: automatic (mc_crystallite_factor times the primary extinction depth),
        1: user-defined (mc_crystallite_thickness).
    mc_crystallite_factor : float, optional
        For calculation_method=1 and mc_crystallite_thickness_flag=0, the ratio of the
        crystallite thickness to the primary extinction depth (amplitude, normal to the
        surface, sigma polarization).
    mc_crystallite_thickness : float, optional
        For calculation_method=1 and mc_crystallite_thickness_flag=1, the crystallite
        thickness in m.
    mc_max_energies : int, optional
        For calculation_method=1, the maximum number of photon energies where the crystallite
        rocking curve is computed (otherwise interpolated). See S4MosaicCrystal.

    Returns
    -------
    instance of S4HyperboloidMosaicCrystal.
    """

    def __init__(self,
                 name="Hyperboloid mosaic crystal",
                 boundary_shape=None,
                 material="graphite",
                 miller_index_h=1,
                 miller_index_k=1,
                 miller_index_l=1,
                 thickness=0.010,
                 f_central=False,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 file_refl="",
                 material_constants_library_flag=0,  # 0=xraylib, 1=dabax
                                                     # 2=shadow preprocessor file v1
                                                     # 3=shadow preprocessor file v2
                 min_axis=0.0,
                 maj_axis=0.0,
                 pole_to_focus=0.0,  # for external calculation
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 dabax=None,
                 mosaicity_fwhm_deg=0.4,
                 mosaicity_profile_flag=0,  # 0=Gaussian, 1=External
                 calculation_method=0,            # 0=macroscopic (1992 paper), 1=Monte Carlo
                 mc_crystallite_thickness_flag=0, # for calculation_method=1: 0=automatic, 1=user-defined
                 mc_crystallite_factor=0.3,       # for mc_crystallite_thickness_flag=0: t0 / extinction depth
                 mc_crystallite_thickness=1e-6,   # for mc_crystallite_thickness_flag=1: t0 in m
                 mc_max_energies=21,              # for calculation_method=1: max number of crystallite curves
                 ):
        p_focus, q_focus, grazing_angle = 1.0, 1.0, 1e-3
        S4HyperboloidOpticalElementDecorator.__init__(self, SurfaceCalculation.EXTERNAL, is_cylinder, cylinder_direction, convexity,
                                                 min_axis, maj_axis, pole_to_focus, p_focus, q_focus, grazing_angle)

        S4MosaicCrystal.__init__(self,
                           name=name,
                           boundary_shape=boundary_shape,
                           surface_shape=self.get_surface_shape_instance(),
                           material=material,
                           miller_index_h=miller_index_h,
                           miller_index_k=miller_index_k,
                           miller_index_l=miller_index_l,
                           thickness=thickness,
                           f_central=f_central,
                           f_phot_cent=f_phot_cent,
                           phot_cent=phot_cent,
                           file_refl=file_refl,
                           material_constants_library_flag=material_constants_library_flag,
                           dabax=dabax,
                           mosaicity_fwhm_deg=mosaicity_fwhm_deg,
                           mosaicity_profile_flag=mosaicity_profile_flag,
                           calculation_method=calculation_method,
                           mc_crystallite_thickness_flag=mc_crystallite_thickness_flag,
                           mc_crystallite_factor=mc_crystallite_factor,
                           mc_crystallite_thickness=mc_crystallite_thickness,
                           mc_max_energies=mc_max_energies,
                           )

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "material": material,
            "miller_index_h": miller_index_h,
            "miller_index_k": miller_index_k,
            "miller_index_l": miller_index_l,
            "thickness": thickness,
            "f_central": f_central,
            "f_phot_cent": f_phot_cent,
            "phot_cent": phot_cent,
            "file_refl": file_refl,
            "material_constants_library_flag": material_constants_library_flag,
            "min_axis": min_axis,
            "maj_axis": maj_axis,
            "pole_to_focus": pole_to_focus,
            "is_cylinder": is_cylinder,
            "cylinder_direction": cylinder_direction,
            "convexity": convexity,
            "dabax": self._get_dabax_txt(),
            "mosaicity_fwhm_deg": mosaicity_fwhm_deg,
            "mosaicity_profile_flag": mosaicity_profile_flag,
            "calculation_method": calculation_method,
            "mc_crystallite_thickness_flag": mc_crystallite_thickness_flag,
            "mc_crystallite_factor": mc_crystallite_factor,
            "mc_crystallite_thickness": mc_crystallite_thickness,
            "mc_max_energies": mc_max_energies,
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
        if self._material_constants_library_flag == 1:
            txt += "\nfrom dabax.dabax_xraylib import DabaxXraylib"
        txt_pre = """

from shadow4.beamline.optical_elements.mosaic_crystals.s4_hyperboloid_mosaic_crystal import S4HyperboloidMosaicCrystal
optical_element = S4HyperboloidMosaicCrystal(name={name!r},
    boundary_shape=boundary_shape, material={material!r},
    miller_index_h={miller_index_h}, miller_index_k={miller_index_k}, miller_index_l={miller_index_l},
    thickness={thickness},
    f_central={f_central}, f_phot_cent={f_phot_cent}, phot_cent={phot_cent},
    file_refl={file_refl!r},
    material_constants_library_flag={material_constants_library_flag}, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
    min_axis={min_axis}, maj_axis={maj_axis}, pole_to_focus={pole_to_focus}, is_cylinder={is_cylinder}, cylinder_direction={cylinder_direction}, convexity={convexity},
    dabax={dabax}, # used when material_constants_library_flag=1
    mosaicity_fwhm_deg={mosaicity_fwhm_deg},
    mosaicity_profile_flag={mosaicity_profile_flag},  # 0=Gaussian, 1=External
    calculation_method={calculation_method},  # 0=macroscopic (1992 paper), 1=Monte Carlo
    mc_crystallite_thickness_flag={mc_crystallite_thickness_flag},  # for Monte Carlo: 0=automatic, 1=user-defined
    mc_crystallite_factor={mc_crystallite_factor},  # for automatic: crystallite thickness / extinction depth
    mc_crystallite_thickness={mc_crystallite_thickness},  # for user-defined: crystallite thickness in m
    mc_max_energies={mc_max_energies},  # for Monte Carlo: max number of energies for the crystallite curves
    )"""
        txt += txt_pre.format(**self.__inputs)

        return txt

class S4HyperboloidMosaicCrystalElement(S4MosaicCrystalElement):
    """
    The Shadow4 hyperboloid crystal element.
    It is made of a S4HyperboloidMosaicCrystal and an ElementCoordinates instance. It also includes the input beam.

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
                 optical_element : S4HyperboloidMosaicCrystal = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4HyperboloidMosaicCrystal(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)

        if not (isinstance(self.get_optical_element().get_surface_shape(), HyperbolicCylinder) or
                isinstance(self.get_optical_element().get_surface_shape(), Hyperboloid)):
            raise ValueError("Wrong Optical Element: only Hyperboloid or Hyperbolic Cylinder shape is accepted")

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
        txt += "\nfrom shadow4.beamline.optical_elements.mosaic_crystals.s4_hyperboloid_mosaic_crystal import S4HyperboloidMosaicCrystalElement"
        txt += "\nbeamline_element = S4HyperboloidMosaicCrystalElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        if self.get_optical_element()._calculation_method == 1: txt += self.to_python_code_mc_crystallite_plot()
        return txt
