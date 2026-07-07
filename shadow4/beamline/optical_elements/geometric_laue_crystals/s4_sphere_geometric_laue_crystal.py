import numpy
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.geometric_laue_crystals.s4_geometric_laue_crystal import \
    S4GeometricLaueCrystalElement, S4GeometricLaueCrystal
from shadow4.beamline.s4_optical_element_decorators import SurfaceCalculation, S4SphereOpticalElementDecorator
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

from syned.beamline.shape import Sphere, SphericalCylinder, Convexity, Direction

from srxraylib.plot.gol import plot_scatter, plot, plot_image, plot_show

class S4SphereGeometricLaueCrystal(S4GeometricLaueCrystal, S4SphereOpticalElementDecorator):
    """
    Shadow4 Sphere Geometric Laue Crystal Class
    This is a spherically curved perfect crystal in Laue (transmission) diffraction geometry, treated with a
    purely geometric (ray-tracing) formalism.

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
    radius : float, optional
        The surface spherical radius.
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
    dabax : None or instance of DabaxXraylib,
        The pointer to the dabax library  (used for material_constants_library_flag=1).

    Returns
    -------
    instance of S4SphereGeometricLaueCrystal.
    """
    def __init__(self,
                 name="Sphere geometric Laue crystal",
                 boundary_shape=None,
                 material=None,
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
                 radius=1.0,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 poisson_ratio=0.22,
                 dabax=None,
                 ):
        p_focus, q_focus, grazing_angle = 1.0, 1.0, 1e-3
        S4SphereOpticalElementDecorator.__init__(self, SurfaceCalculation.EXTERNAL, is_cylinder, cylinder_direction, convexity,
                                                 radius, p_focus, q_focus, grazing_angle)

        S4GeometricLaueCrystal.__init__(self,
                           name=name,
                           boundary_shape=boundary_shape,
                           surface_shape=self.get_surface_shape_instance(),
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
                           file_refl=file_refl,
                           f_bragg_a=f_bragg_a,
                           f_ext=f_ext,
                           material_constants_library_flag=material_constants_library_flag,
                           poisson_ratio=poisson_ratio,
                           dabax=dabax,
                           )

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
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
            "radius": radius,
            "is_cylinder": is_cylinder,
            "cylinder_direction": cylinder_direction,
            "convexity": convexity,
            "poisson_ratio": poisson_ratio,
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

from shadow4.beamline.optical_elements.geometric_laue_crystals.s4_sphere_geometric_laue_crystal import S4SphereGeometricLaueCrystal
optical_element = S4SphereGeometricLaueCrystal(name='{name}',
    boundary_shape=boundary_shape, material='{material}',
    miller_index_h={miller_index_h}, miller_index_k={miller_index_k}, miller_index_l={miller_index_l},
    f_bragg_a={f_bragg_a}, asymmetry_angle={asymmetry_angle},
    is_thick={is_thick}, thickness={thickness},
    f_central={f_central}, f_phot_cent={f_phot_cent}, phot_cent={phot_cent},
    file_refl='{file_refl}',
    f_ext={f_ext},
    material_constants_library_flag={material_constants_library_flag}, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
    radius={radius:f}, is_cylinder={is_cylinder:d}, cylinder_direction={cylinder_direction:d}, convexity={convexity:d},
    poisson_ratio={poisson_ratio:f},
    dabax={dabax}, # used when material_constants_library_flag=1,
    )"""
        txt += txt_pre.format(**self.__inputs)

        return txt

class S4SphereGeometricLaueCrystalElement(S4GeometricLaueCrystalElement):
    """
    The Shadow4 sphere geometric Laue crystal element.
    It is made of a S4SphereGeometricLaueCrystal and an ElementCoordinates instance. It also includes the input beam.

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
                 optical_element : S4SphereGeometricLaueCrystal = None,
                 coordinates : ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4SphereGeometricLaueCrystal(),
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
        txt += "\nfrom shadow4.beamline.optical_elements.geometric_laue_crystals.s4_sphere_geometric_laue_crystal import S4SphereGeometricLaueCrystalElement"
        txt += "\nbeamline_element = S4SphereGeometricLaueCrystalElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        return txt

def caustic(beam):
    #
    # caustic
    #
    print(">>> Intensity: ", beam.get_intensity(nolost=1))
    beam_to_analyze = beam.duplicate()
    y_min, y_max, npositions = 0.01, 10.0, 100
    x_min, x_max, npoints_x = -0.02, 0.02, 300
    nolost = 1

    positions = numpy.linspace(y_min, y_max, npositions)
    out_x = numpy.zeros((npoints_x, npositions))
    fwhm = numpy.zeros(npositions)
    center = numpy.zeros(npositions)
    col = 3
    ref = 23
    col_title = "Z (col 3)"

    for i in range(npositions):
        beami = beam_to_analyze.duplicate()
        beami.retrace(positions[i], resetY=True)
        tkt_x = beami.histo1(col, xrange=[x_min, x_max], nbins=npoints_x, nolost=nolost, ref=ref)
        out_x[:, i] = tkt_x['histogram']
        fwhm[i] = tkt_x['fwhm']
        if ref == 23:
            center[i] = numpy.average(beami.get_column(col, nolost=nolost),
                                      weights=beami.get_column(23, nolost=nolost))
        else:
            center[i] = numpy.average(beami.get_column(col, nolost=nolost))
    #
    # plots
    #
    print('Result arrays X,Y (shapes): ', out_x.shape, tkt_x['bin_center'].shape, positions.shape)
    x = tkt_x['bin_center']
    y = positions

    # 2D
    plot_image(out_x.T, y, 1e6 * x,
               title="", ytitle="%s [um] (%d pixels)" % (col_title, x.size),
               xtitle="Y [m] (%d pixels)" % (y.size), aspect="auto")
    # FWHM
    fwhm[fwhm == 0] = "nan"
    plot(y, 1e6 * fwhm, title="FWHM",
         xtitle="y [m]", ytitle="FHWH [um]", marker=".")
    # I0
    nx, ny = out_x.shape
    I0 = out_x.T[:, nx // 2]
    plot(y, I0, title="I at central profile", xtitle="y [m]", ytitle="I0", marker=".")
    # center
    plot(y, 1e6 * center, title="CENTER",
         xtitle="y [m]", ytitle="CENTER [um]", marker=".", yrange=[1e6 * x_min, 1e6 * x_max])



if __name__ == "__main__":
    from srxraylib.plot.gol import plot_scatter, plot, plot_image, plot_show

    if 1: # example 1


        from shadow4.tools.logger import set_verbose, printlog
        set_verbose(1)

        import numpy as np
        from dabax.dabax_xraylib import DabaxXraylib
        from shadow4.beamline.s4_beamline import S4Beamline

        beamline = S4Beamline()

        #
        #
        #
        from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

        light_source = SourceGeometrical(name='Geometrical Source', nrays=5000, seed=5676561)
        light_source.set_spatial_type_point()
        light_source.set_depth_distribution_off()
        light_source.set_angular_distribution_gaussian(sigdix=0.000885, sigdiz=7.2e-05)
        light_source.set_energy_distribution_uniform(value_min=9900, value_max=10100, unit='eV')
        light_source.set_polarization(polarization_degree=1, phase_diff=0, coherent_beam=0)
        beam = light_source.get_beam()

        beamline.set_light_source(light_source)

        # optical element number XX
        boundary_shape = None

        from shadow4.beamline.optical_elements.crystals.s4_sphere_crystal import S4SphereCrystal

        # optical_element = S4SphereCrystal(name='Generic Crystal 1:1 focusing (Rowland)',
        optical_element = S4SphereGeometricLaueCrystal(name='Generic Crystal 1:1 focusing (Rowland)',
                                          boundary_shape=boundary_shape, material='Si',
                                          miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                          f_bragg_a=True, asymmetry_angle=numpy.radians(90),
                                          is_thick=0, thickness=0.01,
                                          f_central=1, f_phot_cent=0, phot_cent=10000.0,
                                          file_refl='bragg.dat',
                                          f_ext=0,
                                          material_constants_library_flag=0,
                                          # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                          radius=1.0, is_cylinder=0, cylinder_direction=0, convexity=1,
                                          dabax=None,  # used when material_constants_library_flag=1,
                                          )
        from syned.beamline.element_coordinates import ElementCoordinates

        coordinates = ElementCoordinates(p=3, q=0, angle_radial=1.371743969, angle_azimuthal=0,
                                         angle_radial_out=1.371743969)
        movements = None
        from shadow4.beamline.optical_elements.crystals.s4_sphere_crystal import S4SphereCrystalElement

        beamline_element = S4SphereGeometricLaueCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                                  movements=movements, input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        plot_scatter(footprint.get_column(2, nolost=1), footprint.get_column(1, nolost=1), title='footprint (Y,X) in m', show=False)

        from shadow4.tools.beamline_tools import focnew

        t = focnew(beamline, beam)
        print(t['text'])

        ###
        ticket = beam.histo1(26, nbins=100, xrange=[np.float64(9950.010761521635), np.float64(10049.981772997284)],
                             nolost=1, ref=23)
        title = "I: %.1f " % ticket['intensity']
        if ticket['fwhm'] is not None: title += "FWHM: %f " % ticket['fwhm']
        plot(ticket['bin_path'], ticket['histogram_path'],
             title=title, xtitle="column 26", show=0)

        ###
        ticket = footprint.histo1(3, nbins=100, xrange=[-0.002, 0.002],
                             nolost=1, ref=23)
        title = "I: %.1f " % ticket['intensity']
        if ticket['fwhm'] is not None: title += "FWHM: %f " % ticket['fwhm']
        plot(ticket['bin_path'], ticket['histogram_path'],
             title=title, xtitle="footprint column 3", show=0)

        plot_show()
        # caustic(beam)

    if 0: # Fig 11 in Qi 2021
        from srxraylib.plot.gol import plot_scatter, plot, plot_show

        from shadow4.tools.logger import set_verbose, printlog

        set_verbose(1)

        import numpy as np
        from dabax.dabax_xraylib import DabaxXraylib
        from shadow4.beamline.s4_beamline import S4Beamline

        beamline = S4Beamline()

        #
        #
        #
        from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

        light_source = SourceGeometrical(name='Geometrical Source', nrays=5000, seed=5676561)
        light_source.set_spatial_type_point()
        light_source.set_depth_distribution_off()
        light_source.set_angular_distribution_gaussian(sigdix=0.000885, sigdiz=7.2e-05)
        # light_source.set_angular_distribution_flat(hdiv1=0, hdiv2=0, vdiv1=0, vdiv2=0)
        # light_source.set_energy_distribution_singleline(34561, unit='eV')
        light_source.set_energy_distribution_uniform(value_min=34561.0+100, value_max=34561.0-100, unit='eV')
        light_source.set_polarization(polarization_degree=1, phase_diff=0, coherent_beam=0)
        beam = light_source.get_beam()

        beamline.set_light_source(light_source)

        # optical element number XX
        boundary_shape = None

        from shadow4.beamline.optical_elements.crystals.s4_sphere_crystal import S4SphereCrystal

        # optical_element = S4SphereCrystal(name='Generic Crystal 1:1 focusing (Rowland)',
        optical_element = S4SphereGeometricLaueCrystal(name='Generic Crystal 1:1 focusing (Rowland)',
                                                       boundary_shape=boundary_shape, material='Si',
                                                       miller_index_h=3, miller_index_k=1, miller_index_l=1,
                                                       f_bragg_a=True, asymmetry_angle=numpy.radians(90),
                                                       is_thick=0, thickness=0.0005,
                                                       f_central=1, f_phot_cent=0, phot_cent=34561.0,
                                                       file_refl='bragg.dat',
                                                       f_ext=0,
                                                       material_constants_library_flag=0,
                                                       # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                                       radius=0.5, is_cylinder=1, cylinder_direction=0, convexity=1,
                                                       dabax=None,  # used when material_constants_library_flag=1,
                                                       )
        from syned.beamline.element_coordinates import ElementCoordinates

        coordinates = ElementCoordinates(p=3, q=0, angle_radial=1.371743969, angle_azimuthal=0,
                                         angle_radial_out=1.371743969)
        movements = None
        from shadow4.beamline.optical_elements.crystals.s4_sphere_crystal import S4SphereCrystalElement

        beamline_element = S4SphereGeometricLaueCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                                               movements=movements, input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # plot_scatter(footprint.get_column(2, nolost=1), footprint.get_column(1, nolost=1), title='footprint (Y,X) in m')

        from shadow4.tools.beamline_tools import focnew

        t = focnew(beamline, beam)
        print(t['text'])


    plot_show()

