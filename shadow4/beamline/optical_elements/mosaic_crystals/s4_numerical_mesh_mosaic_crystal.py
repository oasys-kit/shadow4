"""
The s4 numerical mesh mosaic crystal (optical element and beamline element).
"""
import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.shape import NumericalMesh
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.mosaic_crystals.s4_mosaic_crystal import S4MosaicCrystalElement, S4MosaicCrystal

from shadow4.beamline.s4_optical_element_decorators import S4NumericalMeshOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4NumericalMeshMosaicCrystal(S4MosaicCrystal, S4NumericalMeshOpticalElementDecorator):
    """
    Shadow4 Mesh Mosaic Crystal Class
    This is a curved mosaic crystal in reflection geometry (Bragg), using the diffracted beam.
    The surface shape is defined as a numerical mesh, either in an hdf5 file or in arrays.

    Constructor.

    Parameters
    ----------
    name :  str, optional
        A name for the crystal
    boundary_shape : instance of BoundaryShape, optional
        The information on the crystal boundaries.
    xx : numpy array, optional
        The array with the X (width) coordinates.
    yy : numpy array, optional
        The array with the Y (length) coordinates.
    zz : numpy array, optional
        The heights with shape (yy.size, xx.size), following the OASYS convention.
    surface_data_file : str, optional
        The h5 file name with the mesh following the Oasys convention.
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
    instance of S4NumericalMeshMosaicCrystal.
    """
    def __init__(self,
                 name="Numerical Mesh mosaic crystal",
                 boundary_shape=None,
                 xx=None,
                 yy=None,
                 zz=None,
                 surface_data_file="",
                 # inputs related to crystal
                 material="graphite",
                 miller_index_h=1,
                 miller_index_k=1,
                 miller_index_l=1,
                 thickness=0.010,
                 f_central=False,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 material_constants_library_flag=0,  # 0=xraylib, 1=dabax
                                                     # 2=shadow preprocessor file v1
                                                     # 3=shadow preprocessor file v2
                 file_refl="",
                 dabax=None,
                 mosaicity_fwhm_deg=0.4,
                 mosaicity_profile_flag=0,  # 0=Gaussian, 1=External
                 calculation_method=0,            # 0=macroscopic (1992 paper), 1=Monte Carlo
                 mc_crystallite_thickness_flag=0, # for calculation_method=1: 0=automatic, 1=user-defined
                 mc_crystallite_factor=0.3,       # for mc_crystallite_thickness_flag=0: t0 / extinction depth
                 mc_crystallite_thickness=1e-6,   # for mc_crystallite_thickness_flag=1: t0 in m
                 mc_max_energies=21,              # for calculation_method=1: max number of crystallite curves
                 ):

        S4NumericalMeshOpticalElementDecorator.__init__(self, xx, yy, zz, surface_data_file)
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
            "xx": xx,
            "yy": yy,
            "zz": zz,
            "surface_data_file": surface_data_file,
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
        txt += "\nimport numpy"
        for key in ("xx", "yy", "zz"):
            value = self.__inputs[key]
            code = "None" if value is None else "numpy.array(%r)" % numpy.asarray(value).tolist()
            txt += "\n%s = %s" % (key, code)
        txt_pre = """

from shadow4.beamline.optical_elements.mosaic_crystals.s4_numerical_mesh_mosaic_crystal import S4NumericalMeshMosaicCrystal
optical_element = S4NumericalMeshMosaicCrystal(name={name!r},boundary_shape=boundary_shape,
    xx=xx, yy=yy, zz=zz, surface_data_file={surface_data_file!r},
    material={material!r},
    miller_index_h={miller_index_h}, miller_index_k={miller_index_k}, miller_index_l={miller_index_l},
    thickness={thickness},
    f_central={f_central}, f_phot_cent={f_phot_cent}, phot_cent={phot_cent},
    file_refl={file_refl!r},
    material_constants_library_flag={material_constants_library_flag}, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
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

class S4NumericalMeshMosaicCrystalElement(S4MosaicCrystalElement):
    """
    The Shadow4 mesh crystal element.
    It is made of a S4NumericalMeshMosaicCrystal and an ElementCoordinates instance. It also includes the input beam.

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
                 optical_element: S4NumericalMeshMosaicCrystal = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4NumericalMeshMosaicCrystal(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)
        if not isinstance(self.get_optical_element().get_surface_shape(), NumericalMesh):
            raise ValueError("Wrong Optical Element: only Surface Data shape is accepted")

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
        txt += "\nfrom shadow4.beamline.optical_elements.mosaic_crystals.s4_numerical_mesh_mosaic_crystal import S4NumericalMeshMosaicCrystalElement"
        txt += "\nbeamline_element = S4NumericalMeshMosaicCrystalElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        if self.get_optical_element()._calculation_method == 1: txt += self.to_python_code_mc_crystallite_plot()
        return txt
