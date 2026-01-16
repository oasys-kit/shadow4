import numpy

from syned.beamline.optical_element import OpticalElement

from shadow4.beamline.s4_optical_element_decorators import S4OpticalElementDecorator

from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

from shadow4.tools.logger import is_verbose

from shadow4.optical_surfaces.s4_conic import S4Conic
from shadow4.beamline.optical_elements.crystals.s4_conic_crystal import S4ConicCrystal, S4ConicCrystalElement
from shadow4.beamline.optical_elements.mirrors.s4_conic_mirror import S4ConicMirror, S4ConicMirrorElement


class S4Compound(OpticalElement, S4OpticalElementDecorator):
    """
    Shadow4 Compound Element Class
    This is a base class for compound elements, that is a list of elements.

    It is intended to group some elements into a single one, with the idea of
    doing the inter-element ray tracing without changing reference system.

    An examples may be the channel-cut monochromator.

    Note that not all S4 optical elements can be used.

    Constructor.

    Parameters
    ----------
    name : str, optional
        The name of the mirror.
    oe_list : list, optional
        A list of S4OpticalElement instances.

    Returns
    -------
    instance of S4Mirror.
    """
    def __init__(self,
                 name="Undefined",
                 oe_list=None,
                 ):

        OpticalElement.__init__(self,
                        name=name,
                        )

        if oe_list is None:
            self.oe_list = []
        else:
            self._oe_list = oe_list

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._add_support_text([
            ("oe_list",            "S4: list of optical elements",                 ""),
        ] )

    def get_info(self):
        """
        Returns the specific information of the S4 mirror optical element.

        Returns
        -------
        str
        """
        txt = "\n\n"
        txt += "COMPOUND OPTICAL ELEMENT\n"

        for i, element in enumerate(self._oe_list):
            txt += "   ELEMENT number %d\n" % (i + 1)
            txt += "\n" + element.get_info()

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
        txt_pre = ""
        for i, oe in enumerate(self._oe_list):
            script = oe.to_python_code()
            indented_script = '\n'.join('    ' + line for line in script.splitlines())
            txt_pre += "\n\ndef get_oe_%d():" % (i + 1)
            txt_pre += indented_script
            txt_pre += "\n    return optical_element"

        txt_pre += "\noe_list=["
        for i, oe in enumerate(self._oe_list):
            txt_pre += "get_oe_%d()," % (i+1)
        txt_pre += "]\n"


        txt = "" # self.to_python_code_boundary_shape()
        txt += "\nfrom shadow4.beamline.optical_elements.compound.s4_compound import S4Compound"
        txt += "\noptical_element = S4Compound(name='%s', oe_list=oe_list)" % self.get_name()
        return txt_pre + txt

class S4CompoundElement(S4BeamlineElement):
    """
    The base class for Shadow4 compound element.
    It is made of a S4Compound and an ElementCoordinates instance. It also includes the input beam.

    Use derived classes for plane or other curved crystal surfaces.

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
    instance of S4MirrorElement.
    """
    def __init__(self,
                 optical_element : S4Compound = None,
                 coordinates : ElementCoordinates = None,
                 movements : S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4Compound(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)

    def trace_beam(self, **params):
        """
        Runs (ray tracing) the input beam through the compound element.

        Parameters
        ----------
        **params

        Returns
        -------
        tuple
            (output_beam, footprint) instances of S4Beam.
        """
        flag_lost_value = params.get("flag_lost_value", -1)
        change_reference_system_in = params.get("change_reference_system_in", True)
        change_reference_system_out = params.get("change_reference_system_out", True)

        if is_verbose():
            if not change_reference_system_in:
                print("change_reference_system_in = False: skipping reference change to o.e.")
            if not change_reference_system_out:
                print("change_reference_system_out = False: skipping reference change from o.e. to image")

        p = self.get_coordinates().p()
        q = self.get_coordinates().q()
        theta_grazing1 = numpy.pi / 2 - self.get_coordinates().angle_radial()
        theta_grazing2 = numpy.pi / 2 - self.get_coordinates().angle_radial_out()
        alpha1 = self.get_coordinates().angle_azimuthal()

        #
        input_beam = self.get_input_beam().duplicate()

        #
        # put beam in mirror reference system
        #
        if change_reference_system_in:
            input_beam.rotate(alpha1, axis=2)
            input_beam.rotate(theta_grazing1, axis=1)
            input_beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])
            print(">>> compound changes to in")

        # mirror movement:
        movements = self.get_movements()
        if movements is not None:
            if movements.f_move:
                input_beam.rot_for(OFFX=movements.offset_x,
                                   OFFY=movements.offset_y,
                                   OFFZ=movements.offset_z,
                                   X_ROT=movements.rotation_x,
                                   Y_ROT=movements.rotation_y,
                                   Z_ROT=movements.rotation_z)


        #
        # reflect beam in the mirror surface
        #
        elements = self.get_optical_element()._oe_list

        print(">>> elements: ", type(elements), elements)

        beam = input_beam.duplicate()
        coordinates = ElementCoordinates(p=0, q=0,
                                         angle_radial=0,
                                         angle_azimuthal=0,
                                         angle_radial_out=0)

        footprints = []
        for i, element in enumerate(elements):
            print(">>> element number %d : " % (i + 1), element)

            print("    >>>> vin: ", beam.get_columns([4, 5, 6])[:, 0])
            if isinstance(element, S4ConicMirror):
                print("    >>> Conic mirror")
                be = S4ConicMirrorElement(optical_element=element,
                                           coordinates=coordinates,
                                           movements=None,
                                           input_beam=beam,
                                           )
                beam, footprint = be.trace_beam(flag_lost_value=flag_lost_value,
                                                change_reference_system_in=False,
                                                change_reference_system_out=False)
            elif isinstance(element, S4ConicCrystal):
                print("    >>> Conic crystal")
                be = S4ConicCrystalElement(optical_element=element,
                                           coordinates=coordinates,
                                           movements=None,
                                           input_beam=beam,
                                           )
                beam, footprint = be.trace_beam(flag_lost_value=flag_lost_value,
                                                change_reference_system_in=False,
                                                change_reference_system_out=False)

            else:
                raise Exception("Not implemented %s in CompoundElement", type(element))

            footprints.append(footprint)
            print("    >>>> intercept: ", footprint.get_columns([1, 2, 3])[:, 0])
            print("    >>>> vout: ", footprint.get_columns([4, 5, 6])[:, 0])

        #
        # from mirror reference system to image plane
        #
        output_beam = footprint.duplicate()
        if change_reference_system_out:
            output_beam.change_to_image_reference_system(theta_grazing2, q)
            print(">>> compound changes to out")

        return output_beam, footprints

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
        txt += "\nfrom shadow4.beamline.optical_elements.compound.s4_compound import S4CompoundElement"
        txt += "\nbeamline_element = S4CompoundElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam, movements=movements)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        return txt



if __name__ == "__main__":

    def get_optical_element_instance_channel_cut(
            crystal_separation=0,
            roll=0, pitch=0, yaw=0,  # applied in this order
            T=[0, 0, 0],
            use_mirrors=False,
    ):

        try:
            name = self.getNode().title
        except:
            name = "Channel Cut Crystal Monochromator"

        boundary_shape = None

        from shadow4.beamline.optical_elements.crystals.s4_conic_crystal import S4ConicCrystal

        ccc1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -0.5 * crystal_separation]
        ccc2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.5 * crystal_separation]

        # Rx, beta
        R_pitch = [[1, 0, 0],
                   [0, numpy.cos(pitch), -numpy.sin(pitch)],
                   [0, numpy.sin(pitch), numpy.cos(pitch)]]

        # Ry, gamma
        R_roll = [[numpy.cos(roll), 0, numpy.sin(roll)],
                  [0, 1, 0],
                  [-numpy.sin(roll), 0, numpy.cos(roll)]]

        # Rz, alpha
        R_yaw = [[numpy.cos(yaw), -numpy.sin(yaw), 0],
                 [numpy.sin(yaw), numpy.cos(yaw), 0],
                 [0, 0, 1]]

        # R = Rz Rx Ry
        R = numpy.array(R_yaw) @ numpy.array(R_pitch) @ numpy.array(R_roll)

        print(">>> R: ", R)

        conic_coefficients1 = S4Conic.rotate_and_translate_coefficients(ccc1, R, T)

        conic_coefficients2 = S4Conic.rotate_and_translate_coefficients(ccc2, R, T)

        if use_mirrors:

            optical_element1 = S4ConicMirror(name='Generic Mirror', boundary_shape=boundary_shape,
                                             conic_coefficients=conic_coefficients1,
                                             f_reflec=0, f_refl=6, file_refl='<none>',
                                             refraction_index=0.99999 + 0.001j,
                                             coating_material='Ni', coating_density=8.902, coating_roughness=0)
            optical_element2 = S4ConicMirror(name='Generic Mirror', boundary_shape=boundary_shape,
                                             conic_coefficients=conic_coefficients2,
                                             f_reflec=0, f_refl=6, file_refl='<none>',
                                             refraction_index=0.99999 + 0.001j,
                                             coating_material='Ni', coating_density=8.902, coating_roughness=0)
        else:
            optical_element1 = S4ConicCrystal(name='Generic Crystal',
                                              boundary_shape=boundary_shape,
                                              conic_coefficients=conic_coefficients1,
                                              material='Si', miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                              f_bragg_a=False, asymmetry_angle=0.0,
                                              is_thick=1, thickness=0.001,
                                              f_central=1, f_phot_cent=0, phot_cent=5000.0,
                                              file_refl='bragg.dat',
                                              f_ext=0,
                                              material_constants_library_flag=1,
                                              # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                              )

            optical_element2 = S4ConicCrystal(name='Generic Crystal',
                                              boundary_shape=boundary_shape,
                                              conic_coefficients=conic_coefficients2,
                                              material='Si', miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                              f_bragg_a=False, asymmetry_angle=0.0,
                                              is_thick=1, thickness=0.001,
                                              f_central=0, f_phot_cent=0, phot_cent=5000.0,
                                              file_refl='bragg.dat',
                                              f_ext=0,
                                              material_constants_library_flag=1,
                                              # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                              )
        #
        #
        #
        return S4Compound(name=name, oe_list=[optical_element1, optical_element2])

    def plot_2d(beam, footprints, irange=[-0.0001, 0.0001], frange=[-0.009,0.009]):
        from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms, plot_show

        # image
        ticket = beam.histo2(1, 3, nbins_h=100, nbins_v=100,
                             xrange=irange, yrange=irange, nolost=1, ref=23)

        title = "BEAM I: %.1f " % ticket['intensity']
        if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
        if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

        plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                                   title=title, xtitle="column 1", ytitle="column 3",
                                   cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=0)
        # footprint 1
        ticket = footprints[0].histo2(2, 1, nbins_h=100, nbins_v=100,
                                      xrange=frange, yrange=frange, nolost=1, ref=23)

        title = "FOOTPRINT FIRST CRYSTAL: %.1f " % ticket['intensity']
        if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
        if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

        plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                                   title=title, xtitle="column 2", ytitle="column 1",
                                   cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=0)

        # footprint 2
        ticket = footprints[1].histo2(2, 1, nbins_h=100, nbins_v=100,
                                      xrange=frange, yrange=frange, nolost=1, ref=23)

        title = "FOOTPRINT SECOND CRYSTAL: %.1f " % ticket['intensity']
        if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
        if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

        plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                                   title=title, xtitle="column 2", ytitle="column 1",
                                   cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=1)

    def print_centroids(beam1, title='', factor=1.0):
        x = beam1.get_column(1, nolost=1)
        y = beam1.get_column(2, nolost=1)
        z = beam1.get_column(3, nolost=1)
        xp = beam1.get_column(4, nolost=1)
        yp = beam1.get_column(5, nolost=1)
        zp = beam1.get_column(6, nolost=1)
        w = beam1.get_column(23, nolost=1)

        arrays = [x, y, z, xp, yp, zp]
        titles = ['x', 'y', 'z', 'xp', 'yp', 'zp']

        arrays = [x, y, z]
        titles = ['x', 'y', 'z']

        print("\n-------------------------", title, "factor= %f" % factor)
        for j, array in enumerate(arrays):
            average = numpy.average(array, weights=w)
            variance = numpy.average((array - average) ** 2, weights=w)
            print(titles[j], factor * average, "+/-", factor * numpy.sqrt(variance))

    # channel-cut
    if 1:
        from shadow4.beamline.s4_beamline import S4Beamline

        beamline = S4Beamline()

        if 0: # undulator
            # electron beam
            from shadow4.sources.s4_electron_beam import S4ElectronBeam

            electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
            electron_beam.set_sigmas_all(sigma_x=3.01836e-05, sigma_y=3.63641e-06, sigma_xp=4.36821e-06,
                                         sigma_yp=1.37498e-06)
            electron_beam.set_dispersion_all(0, 0, 0, 0)

            # magnetic structure
            from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian

            source = S4UndulatorGaussian(
                period_length=0.042,  # syned Undulator parameter (length in m)
                number_of_periods=38.571,  # syned Undulator parameter
                photon_energy=5000.0,  # Photon energy (in eV)
                delta_e=2.0,  # Photon energy width (in eV)
                ng_e=100,  # Photon energy scan number of points
                flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
                flag_energy_spread=0,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
                harmonic_number=1,  # harmonic number
                flag_autoset_flux_central_cone=0,  # value to set the flux peak
                flux_central_cone=10000000000.0,  # value to set the flux peak
            )

            # light source
            from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource

            light_source = S4UndulatorGaussianLightSource(name='Undulator Gaussian', electron_beam=electron_beam,
                                                          magnetic_structure=source, nrays=50000, seed=5676561)
        else: # pencil

            from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

            light_source = SourceGeometrical(name='Geometrical Source', nrays=500, seed=5676561)
            light_source.set_spatial_type_point()
            light_source.set_depth_distribution_off()
            light_source.set_angular_distribution_flat(hdiv1=0.000000, hdiv2=0.000000, vdiv1=0.000000, vdiv2=0.000000)
            light_source.set_energy_distribution_uniform(value_min=4999.000000, value_max=5001.000000, unit='eV')
            light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)


        beam = light_source.get_beam()

        beamline.set_light_source(light_source)

        optical_element = get_optical_element_instance_channel_cut(
            crystal_separation = 0.0005,
            roll=0, pitch=0, yaw= 0, # applied in this order
            T=[0,0,0],
            use_mirrors=1,
        )

        from syned.beamline.element_coordinates import ElementCoordinates

        coordinates = ElementCoordinates(p=38.1, q=1,
                                         angle_radial=numpy.radians(66.703995),
                                         angle_azimuthal=0,
                                         angle_radial_out=numpy.radians(180 - 66.703995))

        from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

        movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=5e-3, offset_z=0,
                                               rotation_x=0, rotation_y=0, rotation_z=0)

        beamline_element = S4CompoundElement(
            optical_element=optical_element,
            coordinates=coordinates,
            movements=movements,
            input_beam=beam)

        beam, footprints = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        print_centroids(footprints[0], title='CRYSTAL 1', factor=1e3)
        print_centroids(footprints[1], title='CRYSTAL 2', factor=1e3)
        print_centroids(beam, title='IMAGE', factor=1e6)
        # print(beamline.to_python_code())
        plot_2d(beam, footprints, irange=[-5e-3, 5e-3], frange=[-20e-3, 20e-3])

