"""
The s4 additional numerical mesh mosaic crystal (optical element and beamline element).
"""
import numpy
from syned.beamline.shape import NumericalMesh
from syned.beamline.element_coordinates import ElementCoordinates
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.mosaic_crystals.s4_mosaic_crystal import S4MosaicCrystalElement, S4MosaicCrystal
from shadow4.beamline.optical_elements.mosaic_crystals.s4_numerical_mesh_mosaic_crystal import S4NumericalMeshMosaicCrystal
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements


class S4AdditionalNumericalMeshMosaicCrystal(S4NumericalMeshMosaicCrystal):
    """
    A mosaic crystal with surface errors added to its ideal shape.
    Diffraction settings and boundaries are inherited from the ideal crystal.

    Parameters
    ----------
    ideal_crystal : instance of S4MosaicCrystal
        The crystal baseline.
    numerical_mesh_crystal : instance of S4NumericalMeshMosaicCrystal
        The numerical mesh to be added to the ideal crystal.
    name : str, optional
        The name of the crystal.

    Returns
    -------
    instance of S4AdditionalNumericalMeshMosaicCrystal.
    """
    def __init__(self,
                 ideal_crystal : S4MosaicCrystal = None,
                 numerical_mesh_crystal : S4NumericalMeshMosaicCrystal = None,
                 name="Mosaic crystal with Additional Numerical Mesh"):
        if ideal_crystal is not None and not isinstance(ideal_crystal, S4MosaicCrystal):
            raise ValueError("ideal_crystal must be a mosaic crystal")
        if numerical_mesh_crystal is not None and not isinstance(numerical_mesh_crystal, S4NumericalMeshMosaicCrystal):
            raise ValueError("numerical_mesh_crystal must be a numerical mesh mosaic crystal")
        S4NumericalMeshMosaicCrystal.__init__(self, name=name,
                 boundary_shape=None if ideal_crystal is None else ideal_crystal.get_boundary_shape(),
                 xx=None if numerical_mesh_crystal is None else numerical_mesh_crystal._curved_surface_shape._xx,
                 yy=None if numerical_mesh_crystal is None else numerical_mesh_crystal._curved_surface_shape._yy,
                 zz=None if numerical_mesh_crystal is None else numerical_mesh_crystal._curved_surface_shape._zz,
                 surface_data_file="" if numerical_mesh_crystal is None else numerical_mesh_crystal._curved_surface_shape._surface_data_file,
                 # inputs related to crystal
                 material                       ="graphite" if ideal_crystal is None else ideal_crystal._material,
                 miller_index_h                 =1 if ideal_crystal is None else ideal_crystal._miller_index_h,
                 miller_index_k                 =1 if ideal_crystal is None else ideal_crystal._miller_index_k,
                 miller_index_l                 =1 if ideal_crystal is None else ideal_crystal._miller_index_l,
                 thickness                      =0.010 if ideal_crystal is None else ideal_crystal._thickness,
                 f_central                      =0 if ideal_crystal is None else ideal_crystal._f_central,
                 f_phot_cent                    =0 if ideal_crystal is None else ideal_crystal._f_phot_cent,
                 phot_cent                      =8000.0 if ideal_crystal is None else ideal_crystal._phot_cent,
                 material_constants_library_flag=0 if ideal_crystal is None else ideal_crystal._material_constants_library_flag,
                 file_refl                      ="" if ideal_crystal is None else ideal_crystal._file_refl,
                 dabax                          =None if ideal_crystal is None else ideal_crystal._dabax,
                 mosaicity_fwhm_deg              =0.4 if ideal_crystal is None else ideal_crystal._mosaicity_fwhm_deg,
                 mosaicity_profile_flag          =0 if ideal_crystal is None else ideal_crystal._mosaicity_profile_flag,
                 calculation_method              =0 if ideal_crystal is None else ideal_crystal._calculation_method,
                 mc_crystallite_thickness_flag   =0 if ideal_crystal is None else ideal_crystal._mc_crystallite_thickness_flag,
                 mc_crystallite_factor           =0.3 if ideal_crystal is None else ideal_crystal._mc_crystallite_factor,
                 mc_crystallite_thickness        =1e-6 if ideal_crystal is None else ideal_crystal._mc_crystallite_thickness,
                 mc_max_energies                 =21 if ideal_crystal is None else ideal_crystal._mc_max_energies,
                 )

        self.__ideal_crystal         = ideal_crystal
        self.__numerical_mesh_crystal = numerical_mesh_crystal

        self.__inputs = {
            "name": name,
            "ideal_crystal": ideal_crystal,
            "numerical_mesh_crystal": numerical_mesh_crystal,
        }

    def ideal_crystal(self):
        """
        get the ideal optical element.

        Returns
        -------
        instance of S4MosaicCrystal
        """
        return self.__ideal_crystal

    def get_ideal(self):
        """
        get the ideal optical element.

        Returns
        -------
        instance of S4MosaicCrystal
        """
        return self.__ideal_crystal

    def get_optical_surface_instance(self):
        """Return a fresh mesh containing the ideal heights plus surface errors."""
        numerical_mesh = super().get_optical_surface_instance()
        if self.__ideal_crystal is not None:
            x, y = numerical_mesh.get_mesh_x_y()
            if x is None or y is None:
                raise ValueError("A numerical mesh is required to add errors to the ideal crystal")
            X, Y = numpy.meshgrid(x, y, indexing="ij")
            ideal = self.__ideal_crystal.get_optical_surface_instance()
            numerical_mesh.add_to_mesh(ideal.surface_height(X, Y))
        return numerical_mesh

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
        txt = "\nfrom dabax.dabax_xraylib import DabaxXraylib"
        if self.__ideal_crystal is None:
            txt += "\nideal_crystal = None"
        else:
            txt += self.__ideal_crystal.to_python_code()
            txt += "\nideal_crystal = optical_element"
        if self.__numerical_mesh_crystal is None:
            txt += "\nnumerical_mesh_crystal = None"
        else:
            txt += self.__numerical_mesh_crystal.to_python_code()
            txt += "\nnumerical_mesh_crystal = optical_element"

        txt += self.to_python_code_boundary_shape()
        txt_pre = """

from shadow4.beamline.optical_elements.mosaic_crystals.s4_additional_numerical_mesh_mosaic_crystal import S4AdditionalNumericalMeshMosaicCrystal
optical_element = S4AdditionalNumericalMeshMosaicCrystal(name={name!r}, ideal_crystal=ideal_crystal, numerical_mesh_crystal=numerical_mesh_crystal)
    """
        txt += txt_pre.format(**self.__inputs)
        return txt

class S4AdditionalNumericalMeshMosaicCrystalElement(S4MosaicCrystalElement):
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
    """
    def __init__(self,
                 optical_element: S4AdditionalNumericalMeshMosaicCrystal = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4AdditionalNumericalMeshMosaicCrystal(),
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

        txt += "\nfrom shadow4.beamline.optical_elements.mosaic_crystals.s4_additional_numerical_mesh_mosaic_crystal import S4AdditionalNumericalMeshMosaicCrystalElement"
        txt += "\nbeamline_element = S4AdditionalNumericalMeshMosaicCrystalElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        if self.get_optical_element()._calculation_method == 1: txt += self.to_python_code_mc_crystallite_plot()
        return txt
