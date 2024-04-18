import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.refractors.crl import CRL
from syned.beamline.shape import Rectangle, Ellipse, Circle

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
from shadow4.beamline.optical_elements.refractors.s4_conic_interface import S4ConicInterfaceElement
from shadow4.beamline.s4_optical_element_decorators import S4RefractiveLensOpticalElementDecorator
from shadow4.beamline.optical_elements.refractors.s4_lens import _get_lens_interfaces

class S4CRL(CRL, S4RefractiveLensOpticalElementDecorator):
    """
    Constructor.

    Parameters
    ----------
    name : str, optional
        The name of the mirror.
    boundary_shape : instance of BoundaryShape, optional
        The boundary shape of the mirror.
    n_lens : int, optional
        The number of (identical) lenses.
    piling_thickness : float, optional
        The distance from one lens to the next one in m.
    material : str, optional
        A string with the material element symbol or compound formula.
    density : float, optional
        The density of the material in the lens in g/cm^3.
    thickness : float, optional
        The thickness of a single lens in m. (Distance between the two interfaces at the center of the lenses.)
    surface_shape : int, optional
        A flag to indicate the shape of the optical surfaces: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients.
    convex_to_the_beam : int, optional
        A flag to indicate the convexity of the first optical surface. Used for surface_shape > 0.
        The first interface exposed to the beam is convex: 0=No, 1=Yes.
        The second interface has opposite convexity.
    cylinder_angle : int, optional
        A flag to indicate is the CRL is 2D0fucusing, aor 1D focusing and in which direction:
        Used for surface_shape > 0. Values are:
            0=CRL is focusing in 2D (not cylindrical),
            1=CRL is focusing in 1D (meridional focusing),
            2=CRL is focusing in 2D (sagittal focusing).
    ri_calculation_mode : int, optional
        A flag to indicate the source of the refraction index. Values are:
            * 0=User,
            * 1=prerefl file,
            * 2=direct calculation using xraylib,
            * 3=direct calculation using dabax.
    prerefl_file : str, optional
        For ri_calculation_mode=1, the prerefl preprocessor file name.
    refraction_index : float, optional
        For ri_calculation_mode=0, the real part of the refraction index.
    attenuation_coefficient : float, optional
        For ri_calculation_mode=0, the attenuation coefficient in m^-1 !!!.
    dabax : None or instance of DabaxXraylib,
        The pointer to the dabax library  (used for f_r_ind > 6).
    radius : float, optional
        For surface_shape=(1,2), the lens radius in m. (For parabolic lenses, it is the radius at the tip for paraboloid.)
    conic_coefficients1 : None or list, optional
        For surface_shape=3, A list with the 10 conic coefficients of interface 1. None is considered as Plane.
    conic_coefficients2 : None or list, optional
        For surface_shape=3, A list with the 10 conic coefficients of interface 2. None is considered as Plane.

    Returns
    -------
    instance of S4CRL.
    """
    def __init__(self,
                 name="Undefined",
                 n_lens=1,
                 piling_thickness=0.0,  # syned stuff,
                 boundary_shape=None,  # syned stuff, replaces "diameter" in the shadow3 append_lens
                 material="",
                 density=1.0,
                 thickness=0.0,        # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
                 surface_shape=1,      # now: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                                       # (in shadow3: 1=sphere 4=paraboloid, 5=plane)
                 convex_to_the_beam=1,  # for surface_shape: convexity of the first interface exposed to the beam 0=No, 1=Yes
                                        # the second interface has opposite convexity
                 cylinder_angle=0,      # for surface_shape: 0=not cylindricaL, 1=meridional 2=sagittal
                 ri_calculation_mode=0,   # source of refraction indices and absorption coefficients
                                 # 0=User
                                 # 1=prerefl file
                                 # 2=direct calculation using xraylib
                                 # 3=direct calculation using dabax
                 prerefl_file=None,    # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
                 refraction_index=1.0, # for ri_calculation_mode=1: n (real)
                 attenuation_coefficient=0.0, # for ri_calculation_mode=1: mu in cm^-1 (real)
                 dabax=None,
                 radius=500e-6,        # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
                 conic_coefficients1=None,   # for surface_shape = 3: the conic coefficients of the first interface
                 conic_coefficients2=None,  # for surface_shape = 3: the conic coefficients of the second interface
                 ):
        """

        Parameters
        ----------

        """
        S4RefractiveLensOpticalElementDecorator.__init__(self,
                                                         surface_shape,
                                                         convex_to_the_beam,
                                                         cylinder_angle,
                                                         ri_calculation_mode,
                                                         prerefl_file,
                                                         refraction_index,
                                                         attenuation_coefficient,
                                                         density,
                                                         dabax,
                                                         radius,
                                                         conic_coefficients1,
                                                         conic_coefficients2)
        CRL.__init__(self,
                     name=name,
                     n_lens=n_lens,
                     surface_shape1=self.get_surface_shape_instance()[0],
                     surface_shape2=self.get_surface_shape_instance()[1],
                     boundary_shape=boundary_shape,
                     material=material,
                     thickness=thickness,
                     piling_thickness=piling_thickness)

        self.__inputs = {
            "name": name,
            "n_lens" : n_lens,
            "boundary_shape": boundary_shape,
            "material": material,
            "thickness": thickness,
            "piling_thickness" : piling_thickness,
            "surface_shape": surface_shape,
            "convex_to_the_beam": convex_to_the_beam,
            "cylinder_angle": cylinder_angle,
            "ri_calculation_mode": ri_calculation_mode,
            "prerefl_file": prerefl_file,
            "refraction_index": refraction_index,
            "attenuation_coefficient": attenuation_coefficient,
            "density": density,
            "dabax": repr(dabax),
            "radius": radius,
            "conic_coefficients1": repr(conic_coefficients1),
            "conic_coefficients2": repr(conic_coefficients2),
        }

    def interthickness(self):
        """
        Returns the interthickness of the beamline element, which is the distance covered by the element along the
        optical axis.
        Elements with a single optical surface (mirrors, crystals, etc.) have interthickness zero.
        Elements like lenses, CRL, transfocators, etc. have interthickness > 0. It is redefined in this method.
        Note that the interthickness is the projection along the (image) optical axis.

        Returns
        -------
        float
        """
        return self._piling_thickness * self._n_lens

    def get_info(self):
        """
        Returns the specific information of the S4 CRL optical element.

        Returns
        -------
        str
        """
        txt = "\n\n"
        txt += "CRL (COMPOUND REFRACTIVE LENS)\n"

        txt += "\n"

        txt += "Source of material refraction index for all lenses:\n"
        if self._ri_calculation_mode == 0:
            txt += "   constant value\n"
            txt += "   index of refraction (real part): %f \n"  % self._refraction_index
            txt += "   attenuation coeff: %f \n" % self._attenuation_coefficient
        elif self._ri_calculation_mode == 1:
            txt += "   file generated by prerefl preprocessor: %s \n" % self._prerefl_file
        elif self._ri_calculation_mode == 2:
            txt += "   calculated using xraylib\n"
            txt += "       material: %s\n" % self.__inputs["material"]
            txt += "       density: %f g/cm3\n" % self.__inputs["density"]
        elif self._ri_calculation_mode == 3:
            txt += "   calculated using dabax\n"
            txt += "       material: %s\n" % self.__inputs["material"]
            txt += "       density: %f g/cm3\n" % self.__inputs["density"]

        txt += "Number of lenses: %d \n" % self._n_lens
        txt += "Radius: %f m\n" % self.__inputs["radius"]
        txt += "Piling thickness: %f m\n" % self._piling_thickness
        txt += "thickness: %f m\n" % self._thickness
        txt += "Total thickness [interthickness]: %f m\n" % self.interthickness()
        txt += "convex_to_the_beam (0:No, 1=Yes): %d \n" % self.__inputs["convex_to_the_beam"]

        if self.__inputs["cylinder_angle"] == 0:
            txt += "2D focusing lens\n"
        elif self.__inputs["cylinder_angle"] == 1:
            txt += "1D focusing lens (meridional)\n"
        elif self.__inputs["cylinder_angle"] == 2:
            txt += "1D focusing lens (sagittal)\n"

        if self._surface_shapes == 0:
            txt = "Lenses are PLANE"
        elif self._surface_shapes == 1:
            txt = "Lenses are SPHERE"
        elif self._surface_shapes == 2:
            txt = "Lenses are PARABOLA"
        elif self._surface_shapes == 3:
            txt = "Lenses are CONIC COEFFICIENTS"



                                       # (in shadow3: 1=sphere 4=paraboloid, 5=plane
        txt += "\n"
        ss = self.get_surface_shape_instance()[0]
        if ss is None:
            txt += "Surface shape 1 is: Plane (** UNDEFINED?? **)\n"
        else:
            txt += "Surface shape 1 is: %s\n" % ss.__class__.__name__
        txt += "Parameters:\n %s\n" % ss.info()

        ss = self.get_surface_shape_instance()[1]
        if ss is None:
            txt += "Surface shape 2 is: Plane (** UNDEFINED?? **)\n"
        else:
            txt += "Surface shape 2 is: %s\n" % ss.__class__.__name__
        txt += "Parameters:\n %s\n" % ss.info()


        #
        # txt += self.get_optical_surface_instance().info() + "\n"
        #
        boundary = self.get_boundary_shape()
        if boundary is None:
            txt += "Surface boundaries not considered (infinite)"
        else:
            txt += "Surface boundaries are: %s\n" % boundary.__class__.__name__
            txt += "    Limits: " + repr( boundary.get_boundaries()) + "\n"
            txt += boundary.info()

        return txt

    def to_python_code_boundary_shape(self):
        """
        Creates a code block with information of boundary shape.

        Returns
        -------
        str
            The text with the code.
        """
        txt = ""
        bs = self._boundary_shape
        if bs is None:
            txt += "\nboundary_shape = None"
        elif isinstance(bs, Rectangle):
            txt += "\nfrom syned.beamline.shape import Rectangle"
            txt += "\nboundary_shape = Rectangle(x_left=%g, x_right=%g, y_bottom=%g, y_top=%g)" % bs.get_boundaries()
        elif isinstance(bs, Circle):
            txt += "\nfrom syned.beamline.shape import Circle"
            txt += "\nboundary_shape = Circle(radius=%g, x_center=%g,y_center=%g)" % bs.get_boundaries()
        elif isinstance(bs, Ellipse):
            txt += "\nfrom syned.beamline.shape import Ellipse"
            txt += "\nboundary_shape = Ellipse(a_axis_min=%g, a_axis_max=%g, b_axis_min=%g, b_axis_max=%g)" % bs.get_boundaries()
        return txt

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
from shadow4.beamline.optical_elements.refractors.s4_crl import S4CRL

optical_element = S4CRL(name='{name:s}',
     n_lens={n_lens:d},
     piling_thickness={piling_thickness:g}, # syned stuff
     boundary_shape=boundary_shape,         # syned stuff, replaces "diameter" in the shadow3 append_lens
     material='{material:s}',  # the material for ri_calculation_mode > 1
     density={density:g}, # the density for ri_calculation_mode > 1
     thickness={thickness}, # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
     surface_shape={surface_shape}, # now: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                                    # (in shadow3: 1=sphere 4=paraboloid, 5=plane)
     convex_to_the_beam={convex_to_the_beam}, # for surface_shape: convexity of the first interface exposed to the beam 0=No, 1=Yes
     cylinder_angle={cylinder_angle}, # for surface_shape: 0=not cylindricaL, 1=meridional 2=sagittal
     ri_calculation_mode={ri_calculation_mode},   # source of refraction indices and absorption coefficients
                                     # 0=User, 1=prerefl file, 2=xraylib, 3=dabax
     prerefl_file='{prerefl_file:s}', # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
     refraction_index={refraction_index:g}, # for ri_calculation_mode=1: n (real)
     attenuation_coefficient={attenuation_coefficient:g}, # for ri_calculation_mode=1: mu in cm^-1 (real)
     dabax={dabax:s}, # the pointer to dabax library
     radius={radius:g}, # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
     conic_coefficients1={conic_coefficients1}, # for surface_shape = 3: the conic coefficients of the single lens interface 1
     conic_coefficients2={conic_coefficients2}, # for surface_shape = 3: the conic coefficients of the single lens interface 2
     )
    """
        txt += txt_pre.format(**self.__inputs)
        return txt

    def get_lens_interfaces(self):
        single_lens_optical_surfaces = self.get_optical_surface_instance()
        lens_interfaces              = numpy.full((self._n_lens, 2), None)
        boundary_shape               = self.get_boundary_shape()

        for lens_index in range(self._n_lens):
            lens_interfaces[lens_index, 0], \
            lens_interfaces[lens_index, 1] = \
                _get_lens_interfaces(lens_optical_surfaces=single_lens_optical_surfaces,
                                     boundary_shape=boundary_shape,
                                     ri_calculation_mode=self._ri_calculation_mode,
                                     refraction_index=self._refraction_index,
                                     attenuation_coefficient=self._attenuation_coefficient,
                                     prerefl_file=self._prerefl_file,
                                     material=self.get_material(),
                                     density=self._density,
                                     dabax=self._dabax,
                                     )

        return lens_interfaces


class S4CRLElement(S4BeamlineElement):
    """
    Constructor.

    Parameters
    ----------
    optical_element : instance of OpticalElement, optional
        The syned optical element.
    coordinates : instance of ElementCoordinates, optional
        The syned element coordinates.
    movements : instance of S4BeamlineElementMovements, optional
        The S4 element movements. (The same movements are applied to the two interfaces. Therefore, each rotation is
        applied around the local axes of each interface, which are different.)
    input_beam : instance of S4Beam, optional
        The S4 incident beam.

    Returns
    -------
    instance of S4CRLElement.
    """
    def __init__(self,
                 optical_element: S4CRL = None,
                 coordinates: ElementCoordinates = None,
                 movements: S4BeamlineElementMovements = None,
                 input_beam: S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4CRL(),
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
        txt += "\nimport numpy"
        txt += self.to_python_code_coordinates()
        txt += self.to_python_code_movements()
        txt += "\nfrom shadow4.beamline.optical_elements.refractors.s4_crl import S4CRLElement"
        txt += "\nbeamline_element = S4CRLElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, mirr = beamline_element.trace_beam()"
        return txt

    def trace_beam(self, **params):
        """
        Runs (ray tracing) the input beam through the element.

        Parameters
        ----------
        **params

        Returns
        -------
        tuple
            (output_beam, footprint) instances of S4Beam.
        """
        input_beam = self.get_input_beam().duplicate()
        movements  = self.get_movements()
        oe         = self.get_optical_element()

        optical_surfaces = oe.get_lens_interfaces()
        p, q, angle_radial, angle_radial_out, angle_azimuthal = self.get_coordinates().get_positions()
        n_lens = oe.get_n_lens()
        for lens_index in range(n_lens):
            if lens_index==0: source_plane = p
            else:             source_plane = oe.get_piling_thickness()
            if lens_index==n_lens-1: image_plane = q
            else:                    image_plane = 0.0

            coordinates_1 = ElementCoordinates(p=source_plane, q=oe.get_thickness()*0.5, angle_radial=angle_radial, angle_radial_out=numpy.pi,          angle_azimuthal=angle_azimuthal)
            coordinates_2 = ElementCoordinates(p=oe.get_thickness()*0.5, q=image_plane,  angle_radial=0.0,          angle_radial_out=angle_radial_out, angle_azimuthal=0.0)

            beamline_element_1 = S4ConicInterfaceElement(optical_element=optical_surfaces[lens_index, 0], coordinates=coordinates_1, movements=movements, input_beam=input_beam)
            if lens_index==0:
                beam1, footprint1 = beamline_element_1.trace_beam()
                n1, mu1, n2, mu2 = beamline_element_1.get_stored_optical_constants()
            else:
                beam1, _          = beamline_element_1.trace_beam(reused_stored_optical_constants=(n1, mu1, n2, mu2))

            beamline_element_2 = S4ConicInterfaceElement(optical_element=optical_surfaces[lens_index, 1], coordinates=coordinates_2, movements=movements, input_beam=beam1)
            if lens_index==n_lens-1:
                beam2, footprint2 = beamline_element_2.trace_beam(reused_stored_optical_constants=(n2, mu2, n1, mu1))
            else:
                beam2, _          = beamline_element_2.trace_beam(reused_stored_optical_constants=(n2, mu2, n1, mu1))

            if lens_index < n_lens-1: input_beam = beam2.duplicate()

        return beam2, [footprint1, footprint2]


if __name__ == "__main__":
    if False:
        from shadow4.physical_models.prerefl.prerefl import PreRefl

        #
        # collimated source
        #
        from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

        src = SourceGeometrical()
        src.set_energy_distribution_singleline(value=5000, unit='eV')
        src.set_spatial_type_rectangle(width=1e-3, height=1e-3)
        src.set_angular_distribution_uniform(0, 0, 0, 0)

        beam = src.get_beam()

        filename = "Be.dat"

        PreRefl.prerefl(interactive=False, SYMBOL="Be", FILE=filename, DENSITY=1.848, E_MIN=4500, E_MAX=5500, E_STEP=1)

        #
        # lens
        #
        lens = S4CRL(n_lens=3, piling_thickness=2.5e-3, prerefl_file=filename)
        e = S4CRLElement(optical_element=lens,
                          coordinates=ElementCoordinates(p=10, q=20, angle_radial=0, angle_azimuthal=0, angle_radial_out=numpy.pi),
                          input_beam=beam)

        e.trace_beam()

        print(e.to_python_code())


    if True:
        from shadow4.tools.logger import set_verbose
        set_verbose()
        from shadow4.beamline.s4_beamline import S4Beamline

        beamline = S4Beamline()

        #
        #
        #
        from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

        light_source = SourceGeometrical(name='SourceGeometrical', nrays=5000, seed=5676561)
        light_source.set_spatial_type_gaussian(sigma_h=6.38000015e-05, sigma_v=0.000064)
        light_source.set_depth_distribution_off()
        light_source.set_angular_distribution_cone(cone_max=0.000010, cone_min=0.000000)
        light_source.set_energy_distribution_singleline(14000.000000, unit='eV')
        light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
        beam = light_source.get_beam()

        beamline.set_light_source(light_source)

        # optical element number XX
        boundary_shape = None
        from shadow4.beamline.optical_elements.refractors.s4_crl import S4CRL

        optical_element = S4CRL(name='Compound Refractive Lens',
                                n_lens=30,
                                piling_thickness=0.000625,  # syned stuff
                                boundary_shape=boundary_shape,
                                # syned stuff, replaces "diameter" in the shadow3 append_lens
                                material='Al',  # the material for ri_calculation_mode > 1
                                density=2.6989,  # the density for ri_calculation_mode > 1
                                thickness=2.4999999999999998e-05,
                                # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
                                surface_shape=1,  # now: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                                # (in shadow3: 1=sphere 4=paraboloid, 5=plane)
                                convex_to_the_beam=0,
                                # for surface_shape: convexity of the first interface exposed to the beam 0=No, 1=Yes
                                cylinder_angle=1,  # for surface_shape: 0=not cylindricaL, 1=meridional 2=sagittal
                                ri_calculation_mode=2,  # source of refraction indices and absorption coefficients
                                # 0=User, 1=prerefl file, 2=xraylib, 3=dabax
                                prerefl_file='Al5_55.dat',
                                # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
                                refraction_index=1,  # for ri_calculation_mode=1: n (real)
                                attenuation_coefficient=0,  # for ri_calculation_mode=1: mu in cm^-1 (real)
                                dabax=None,  # the pointer to dabax library
                                radius=0.0003,
                                # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
                                conic_coefficients1=None,
                                # for surface_shape = 3: the conic coefficients of the single lens interface 1
                                conic_coefficients2=None,
                                # for surface_shape = 3: the conic coefficients of the single lens interface 2
                                )

        import numpy
        from syned.beamline.element_coordinates import ElementCoordinates

        coordinates = ElementCoordinates(p=30, q=1.8987, angle_radial=0, angle_azimuthal=0,
                                         angle_radial_out=3.141592654)
        movements = None
        from shadow4.beamline.optical_elements.refractors.s4_crl import S4CRLElement

        beamline_element = S4CRLElement(optical_element=optical_element, coordinates=coordinates, movements=movements,
                                        input_beam=beam)

        beam, mirr = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # test plot
        if True:
            from srxraylib.plot.gol import plot_scatter

            # plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
            #              title='(Intensity,Photon Energy)', plot_histograms=0)
            plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1),
                         title='(X,Z) in microns')

        print(beam.intensity())
        print(optical_element.get_info())