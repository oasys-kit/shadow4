from syned.beamline.shape import NumericalMesh
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.crystals.s4_crystal import S4CrystalElement, S4Crystal, ElementCoordinates

from shadow4.beamline.s4_optical_element_decorators import S4NumericalMeshOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

class S4NumericalMeshCrystal(S4Crystal, S4NumericalMeshOpticalElementDecorator):
    """
    Shadow4 Mesh Crystal Class
    This is a curved perfect crystal in reflection geometry (Bragg), using the diffracted beam.
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
        The array with the Z(X,Y) coordinates.
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
    instance of S4NumericalMeshCrystal.
    """
    def __init__(self,
                 name="Numerical Mesh Crystal",
                 boundary_shape=None,
                 xx=None,
                 yy=None,
                 zz=None,
                 surface_data_file="",
                 # inputs related to crystal
                 material=None,
                 miller_index_h=1,
                 miller_index_k=1,
                 miller_index_l=1,
                 f_bragg_a=False,
                 asymmetry_angle=0.0,
                 is_thick=0,  # 1=Use thick crystal approximation
                 thickness=0.010,
                 f_central=False,
                 f_phot_cent=0,
                 phot_cent=8000.0,
                 f_ext=0,
                 material_constants_library_flag=0,  # 0=xraylib, 1=dabax
                                                     # 2=shadow preprocessor file v1
                                                     # 3=shadow preprocessor file v1
                 file_refl="",
                 ):

        S4NumericalMeshOpticalElementDecorator.__init__(self, xx, yy, zz, surface_data_file)
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
            "xx": xx,
            "yy": yy,
            "zz": zz,
            "surface_data_file": surface_data_file,
            "material": material,
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
        
from shadow4.beamline.optical_elements.crystals.s4_numerical_mesh_crystal import S4NumericalMeshCrystal
optical_element = S4NumericalMeshCrystal(name='{name:s}',boundary_shape=boundary_shape,
    xx=None,yy=None,zz=None,surface_data_file='{surface_data_file:s}',
    material='{material}',
    miller_index_h={miller_index_h}, miller_index_k={miller_index_k}, miller_index_l={miller_index_l},
    f_bragg_a={f_bragg_a}, asymmetry_angle={asymmetry_angle},
    is_thick={is_thick}, thickness={thickness},
    f_central={f_central}, f_phot_cent={f_phot_cent}, phot_cent={phot_cent},
    file_refl='{file_refl}',
    f_ext={f_ext},
    material_constants_library_flag={material_constants_library_flag}, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
    )"""
        txt += txt_pre.format(**self.__inputs)
        return txt

class S4NumericalMeshCrystalElement(S4CrystalElement):
    """
    The Shadow4 mesh crystal element.
    It is made of a S4meshCrystal and an ElementCoordinates instance. It also includes the input beam.

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
                 optical_element: S4NumericalMeshCrystal = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4NumericalMeshCrystal(),
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
        txt += "\nfrom shadow4.beamline.optical_elements.crystals.s4_numerical_mesh_crystal import S4NumericalMeshCrystalElement"
        txt += "\nbeamline_element = S4NumericalMeshCrystalElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":
    a = S4NumericalMeshCrystal(name="")
    b = S4NumericalMeshCrystalElement(optical_element=a)
    print(b.info())
    print(b.to_python_code())

    #
    #
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

    light_source = SourceGeometrical(name='SourceGeometrical', nrays=25000, seed=5676561)
    light_source.set_spatial_type_gaussian(sigma_h=5.70000011e-07, sigma_v=0.000007)
    light_source.set_depth_distribution_off()
    light_source.set_angular_distribution_gaussian(sigdix=0.000088, sigdiz=0.000007)
    light_source.set_energy_distribution_uniform(value_min=9990.000000, value_max=10010.000000, unit='eV')
    light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
    beam = light_source.get_beam()

    # optical element number XX

    #write file
    import numpy

    # calculate a Gaussian bump
    # create an array of 2 cm length
    npoints = 51
    length = 2.0e-2

    x = numpy.linspace(-0.5 * length, 0.5 * length, npoints)
    y = numpy.linspace(-0.5 * length * 2, 0.5 * length * 2, npoints * 2)
    z = numpy.zeros((x.size, y.size))
    # write file for h5
    from oasys.util.oasys_util import write_surface_file
    write_surface_file(z.T, x, y, 'bump.h5', overwrite=True)

    # run
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal

    optical_element = S4NumericalMeshCrystal(name='Plane Crystal',
                                     boundary_shape=None, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     f_bragg_a=False, asymmetry_angle=0.0,
                                     is_thick=1, thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=10000.0,
                                     file_refl='/users/srio/Oasys/Si5_55.111',
                                     f_ext=0,
                                     material_constants_library_flag=2,
                                     surface_data_file="bump.h5",
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=30, q=10, angle_radial=1.37174, angle_azimuthal=0, angle_radial_out=1.37174)
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement

    beamline_element = S4NumericalMeshCrystalElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    # test plot
    if True:
        from srxraylib.plot.gol import plot_scatter

        plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
                     title='(Intensity,Photon Energy)', plot_histograms=0)
        # plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns')



