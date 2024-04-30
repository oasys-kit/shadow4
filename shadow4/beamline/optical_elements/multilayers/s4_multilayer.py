import numpy

from syned.beamline.shape import Rectangle, Ellipse
from syned.beamline.element_coordinates import ElementCoordinates

# TODO: change from syned.beamline.optical_elements.multilayers.multilayer import Multilayer
from shadow4.beamline.optical_elements.multilayers.syned_multilayer import Multilayer

from shadow4.physical_models.mlayer.mlayer import MLayer
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

from shadow4.optical_surfaces.s4_conic import S4Conic
from shadow4.optical_surfaces.s4_toroid import S4Toroid
from shadow4.optical_surfaces.s4_mesh import S4Mesh

class S4Multilayer(Multilayer):
    """
    Shadow4 Mirror Class
    This is a base class for mirrors.
    Use derived classes for plane or other curved mirror surfaces.

    Constructor.

    Parameters
    ----------
    name : str, optional
        The name of the mirror.
    boundary_shape : instance of BoundaryShape, optional
        The boundary shape of the mirror.
    surface_shape : instance of SurfaceShape, optional
        The surface shape of the mirror.
    f_reflec : int, optional
         the reflectivity of surface:
            - 0=no reflectivity,
            - 1=full polarization.
    f_refl : int, optional
        A flag to indicate the source of reflectivities:
            * 0=prerefl file,
            * 1=electric susceptibility,
            * 2=user defined file (1D angle in mrad, reflectivity),
            * 3=user defined file (1D energy in eV, reflectivity),
            * 4=user defined file (2D energy in eV, angle in mrad, reflectivity),
            * 5=direct calculation using xraylib,
            * 6=direct calculation using dabax.
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
    instance of S4Mirror.
    """
    def __init__(self,
                 name="Undefined",
                 boundary_shape=None,
                 surface_shape=None,
                 # inputs related to mirror reflectivity
                 f_refl=0,      # 0=prerefl file
                 file_refl="<not defined>",  # pre_mlayer preprocessor file for f_refl=0
                 structure='[B/W]x50+Si',
                 period=25.0,
                 Gamma=0.5,
                 ):

        Multilayer.__init__(self,
                        name=name,
                        surface_shape=surface_shape,
                        boundary_shape=boundary_shape,
                        structure=structure,
                        period=period,
                        Gamma=Gamma,
                        )

        # reflectivity
        self._f_refl = f_refl
        self._file_refl = file_refl

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._add_support_text([
                    ("f_refl",             "S4: refl. source for f_reflec=1: 0=prerefl, 5=xraylib, 6=dabax", ""),
                    ("file_refl",          "S4: for f_refl=0: file name",               ""),
            ] )

    def get_info(self):
        """
        Returns the specific information of the S4 mirror optical element.

        Returns
        -------
        str
        """
        txt = "\n\n"
        txt += "MULTILAYER\n"

        txt += "Reflectivity calculation:\n"
        if self._f_refl == 0:
            txt += "   Calculated reflectivity from preprocessor (pre_mlayer) file: %s\n" % self._file_refl
        elif self._f_refl == 1:
            txt += "   Calculated reflectivity from user defined file (1D reflectivity vs angle): %s\n" % self._file_refl
        elif self._f_refl == 2:
            txt += "   Calculated reflectivity from user defined file (1D reflectivity vs energy): %s\n" % self._file_refl
        elif self._f_refl == 3:
            txt += "   Calculated reflectivity from user defined file (2D reflectivity vs energy and angle): %s\n" % self._file_refl
        elif self._f_refl == 4:
            txt += "   Calculated reflectivity using xraylib\n"
        elif self._f_refl == 5:
            txt += "   Calculated reflectivity using dabax\n"

        txt += "\n"
        ss = self.get_surface_shape()
        if ss is None:
            txt += "Surface shape is: Plane (** UNDEFINED?? **)\n"
        else:
            txt += "Surface shape is: %s\n" % ss.__class__.__name__

        try:
            sc = self._surface_calculation
            if sc == 0:
                txt += "Multilayer optical surface parameters COMPUTED\n"
                txt += "    Objective focus at p: %f m\n" % ss.get_p_focus()
                txt += "    Image focus at p: %f m\n" % ss.get_q_focus()
                txt += "    Incidence angle (grazing): %f mrad\n" % (1e3 * ss.get_grazing_angle())
            else:
                txt += "Multilayer optical surface parameters EXTERNAL\n"
        except:
            pass

        txt += "\nParameters:\n %s\n" % ss.info()

        txt += self.get_optical_surface_instance().info() + "\n"

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
        txt = "" # "\nfrom shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror"
        bs = self._boundary_shape
        if bs is None:
            txt += "\nboundary_shape = None"
        elif isinstance(bs, Rectangle):
            txt += "\nfrom syned.beamline.shape import Rectangle"
            txt += "\nboundary_shape = Rectangle(x_left=%g, x_right=%g, y_bottom=%g, y_top=%g)" % bs.get_boundaries()
        elif isinstance(bs, Ellipse):
            txt += "\nfrom syned.beamline.shape import Ellipse"
            txt += "\nboundary_shape = Ellipse(a_axis_min=%g, a_axis_max=%g, b_axis_min=%g, b_axis_max=%g)" % bs.get_boundaries()
        return txt

    def _apply_multilayer_reflection(self, beam):
        sur = self.get_optical_surface_instance()
        footprint, normal, _, _, _, _, _ = sur.apply_specular_reflection_on_beam(beam)
        return footprint, normal


class S4MultilayerElement(S4BeamlineElement):
    """
    The base class for Shadow4 mirror element.
    It is made of a S4Mirror and an ElementCoordinates instance. It also includes the input beam.

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
                 optical_element : S4Multilayer = None,
                 coordinates : ElementCoordinates = None,
                 movements : S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4Multilayer(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
                         input_beam=input_beam)
    
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
        flag_lost_value = params.get("flag_lost_value", -1)

        p = self.get_coordinates().p()
        q = self.get_coordinates().q()
        theta_grazing1 = numpy.pi / 2 - self.get_coordinates().angle_radial()
        alpha1 = self.get_coordinates().angle_azimuthal()

        #
        input_beam = self.get_input_beam().duplicate()

        #
        # put beam in mirror reference system
        #
        input_beam.rotate(alpha1, axis=2)
        input_beam.rotate(theta_grazing1, axis=1)
        input_beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

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
        soe = self.get_optical_element() #._optical_element_syned

        v_in = input_beam.get_columns([4,5,6])

        # footprint, normal = self.apply_local_reflection(input_beam)
        footprint, normal = self.get_optical_element()._apply_multilayer_reflection(input_beam)

        if movements is not None:
            if movements.f_move:
                footprint.rot_back(OFFX=movements.offset_x,
                                   OFFY=movements.offset_y,
                                   OFFZ=movements.offset_z,
                                   X_ROT=movements.rotation_x,
                                   Y_ROT=movements.rotation_y,
                                   Z_ROT=movements.rotation_z)
        #
        # apply mirror boundaries
        #
        footprint.apply_boundaries_syned(soe.get_boundary_shape(), flag_lost_value=flag_lost_value)

        #
        # apply mirror reflectivity
        # TODO: add phase
        #

        v_out = input_beam.get_columns([4, 5, 6])
        angle_in = numpy.arccos(v_in[0,:] * normal[0,:] +
                                v_in[1,:] * normal[1,:] +
                                v_in[2,:] * normal[2,:])

        angle_out = numpy.arccos(v_out[0,:] * normal[0,:] +
                                 v_out[1,:] * normal[1,:] +
                                 v_out[2,:] * normal[2,:])

        grazing_angle_mrad = 1e3 * (numpy.pi / 2 - angle_in)

        # TODO: it should be checked why s4_conic gives a downwards normal and s4_mesh an upwards normal
        # This causes negative angles with s4_mesh. Therefore abs() is used
        grazing_angle_mrad = numpy.abs(grazing_angle_mrad)

        if soe._f_refl == 0: # prerefl
            preprocessor_file = soe._file_refl
            pr = MLayer()
            print(">>>>>>>>>>>>>>>>>>>>>>> ", preprocessor_file )
            pr.read_preprocessor_file(preprocessor_file)

            print(">>>> grazing angle mrad: ", grazing_angle_mrad)
            print(">>>> grazing angle deg: ", numpy.degrees(grazing_angle_mrad * 1e-3) )
            print(">>>> energy eV: ", input_beam.get_column(26))
            # grazing_angle_deg, photon_energy_ev
            Rs, Rp, phase_s, phase_p = pr.reflectivity(numpy.degrees(grazing_angle_mrad*1e-3),
                                                       input_beam.get_column(26),
                                                       Y=footprint.get_column(2))
            # from srxraylib.plot.gol import plot
            # plot(input_beam.get_column(26), Rs**2)
            footprint.apply_reflectivities(Rs, Rp)
            # todo: apply phases

        elif soe._f_refl == 1:  # user angle, mrad ref

            values = numpy.loadtxt(soe._file_refl)

            mirror_grazing_angles = values[:, 0]
            mirror_reflectivities = values[:, 1]

            if mirror_grazing_angles[-1] < mirror_grazing_angles[0]: # XOPPY MLayer gives angles in descendent order
                mirror_grazing_angles = values[:, 0][::-1]
                mirror_reflectivities = values[:, 1][::-1]

            Rs = numpy.interp(grazing_angle_mrad,
                              mirror_grazing_angles,
                              mirror_reflectivities,
                              left=mirror_reflectivities[0],
                              right=mirror_reflectivities[-1])
            Rp = Rs
            footprint.apply_reflectivities(numpy.sqrt(Rs), numpy.sqrt(Rp))

        elif soe._f_refl == 2:  # user energy

            beam_energies = input_beam.get_photon_energy_eV()

            values = numpy.loadtxt(soe._file_refl)

            mirror_energies = values[:, 0]
            mirror_reflectivities = values[:, 1]

            Rs = numpy.interp(beam_energies,
                              mirror_energies,
                              mirror_reflectivities,
                              left=mirror_reflectivities[0],
                              right=mirror_reflectivities[-1])
            Rp = Rs
            footprint.apply_reflectivities(numpy.sqrt(Rs), numpy.sqrt(Rp))

        elif soe._f_refl == 3:  # user 2D
            values = numpy.loadtxt(soe._file_refl)

            beam_energies = input_beam.get_photon_energy_eV()

            mirror_energies       = values[:, 0]
            mirror_grazing_angles = values[:, 1]
            mirror_energies         = numpy.unique(mirror_energies)
            mirror_grazing_angles   = numpy.unique(mirror_grazing_angles)
            # if self.user_defined_angle_units  == 0: mirror_grazing_angles = numpy.degrees(1e-3*mirror_grazing_angles)
            # if self.user_defined_energy_units == 1: mirror_energies *= 1e3 # KeV to eV

            def get_interpolator_weight_2D(mirror_energies, mirror_grazing_angles, mirror_reflectivities):
                mirror_reflectivities = numpy.reshape(mirror_reflectivities, (mirror_energies.shape[0], mirror_grazing_angles.shape[0]))
                from scipy.interpolate import  RectBivariateSpline
                interpolator = RectBivariateSpline(mirror_energies, mirror_grazing_angles, mirror_reflectivities, kx=2, ky=2)

                interpolated_weight = numpy.zeros(beam_energies.shape[0])
                for energy, angle, i in zip(beam_energies, grazing_angle_mrad, range(interpolated_weight.shape[0])):
                    interpolated_weight[i] = interpolator(energy, angle)
                interpolated_weight[numpy.where(numpy.isnan(interpolated_weight))] = 0.0

                return interpolated_weight

            if values.shape[1] == 3:
                mirror_reflectivities = values[:, 2]

                Rs = get_interpolator_weight_2D(mirror_energies, mirror_grazing_angles, mirror_reflectivities)
                Rp = Rs
                footprint.apply_reflectivities(numpy.sqrt(Rs), numpy.sqrt(Rp))
                footprint.apply_reflectivities(numpy.sqrt(Rs), numpy.sqrt(Rp))

            elif values.shape[1] == 4:
                mirror_reflectivities_s = values[:, 2]
                mirror_reflectivities_p = values[:, 3]

                Rs = get_interpolator_weight_2D(mirror_energies, mirror_grazing_angles, mirror_reflectivities_s)
                Rp = get_interpolator_weight_2D(mirror_energies, mirror_grazing_angles, mirror_reflectivities_p)

            footprint.apply_reflectivities(numpy.sqrt(Rs), numpy.sqrt(Rp))

        elif soe._f_refl == 4: # xraylib
            print(">>>>>>>>>>>>>>>>>>>>>>> xraylib")
            pr = MLayer.initialize_from_bilayer_stack_in_compressed_format(
                                    structure=soe._structure,
                                    density_O=None,  roughness_O=0.0,
                                    density_E=None,  roughness_E=0.0,
                                    density_S=None,  roughness_S=0.0,
                                    bilayer_thickness=soe._period,
                                    bilayer_gamma=soe._Gamma,
                                    use_xraylib_or_dabax=0,
                                    )

            print(">>>> grazing angle mrad: ", grazing_angle_mrad)
            print(">>>> grazing angle deg: ", numpy.degrees(grazing_angle_mrad * 1e-3) )
            print(">>>> energy eV: ", input_beam.get_column(26))
            # grazing_angle_deg, photon_energy_ev
            Rs, Rp, phase_s, phase_p = pr.reflectivity(numpy.degrees(grazing_angle_mrad*1e-3),
                                                       input_beam.get_column(26),
                                                       Y=footprint.get_column(2))
            # from srxraylib.plot.gol import plot
            # plot(input_beam.get_column(26), Rs**2)
            footprint.apply_reflectivities(Rs, Rp)
            # todo: apply phases

        elif soe._f_refl == 5: # dabax
            print(">>>>>>>>>>>>>>>>>>>>>>> dabax")
            pr = MLayer.initialize_from_bilayer_stack_in_compressed_format(
                                    structure=soe._structure,
                                    density_O=None,  roughness_O=0.0,
                                    density_E=None,  roughness_E=0.0,
                                    density_S=None,  roughness_S=0.0,
                                    bilayer_thickness=soe._period,
                                    bilayer_gamma=soe._Gamma,
                                    use_xraylib_or_dabax=1,
                                    dabax=None,
                                    )

            print(">>>> grazing angle mrad: ", grazing_angle_mrad)
            print(">>>> grazing angle deg: ", numpy.degrees(grazing_angle_mrad * 1e-3) )
            print(">>>> energy eV: ", input_beam.get_column(26))
            # grazing_angle_deg, photon_energy_ev
            Rs, Rp, phase_s, phase_p = pr.reflectivity(numpy.degrees(grazing_angle_mrad*1e-3),
                                                       input_beam.get_column(26),
                                                       Y=footprint.get_column(2))
            # from srxraylib.plot.gol import plot
            # plot(input_beam.get_column(26), Rs**2)
            footprint.apply_reflectivities(Rs, Rp)
            # todo: apply phases

        else:
            raise Exception("Not implemented source of multilayer reflectivity")


        #
        # from multilayer reference system to image plane
        #

        output_beam = footprint.duplicate()
        output_beam.change_to_image_reference_system(theta_grazing1, q)

        return output_beam, footprint

    #
    # i/o utilities
    #
    def set_grazing_angle(self, theta_grazing, theta_azimuthal=None): # TODO: remove if not used
        """
        Sets the grazing angle.

        Parameters
        ----------
        theta_grazing : float, optional
            The grazing angle in rad.
        theta_azimuthal : float, optional
            The azimuthal angle in rad.
        """
        self.get_coordinates()._angle_radial     = numpy.pi / 2 - theta_grazing
        self.get_coordinates()._angle_radial_out = numpy.pi / 2 - theta_grazing
        if theta_azimuthal is not None: self.get_coordinates()._angle_azimuthal = theta_azimuthal

if __name__ == "__main__":

    #
    # from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
    #
    # light_source = SourceGeometrical(name='SourceGeometrical', nrays=10000, seed=5676561)
    # light_source.set_spatial_type_gaussian(sigma_h=5e-06, sigma_v=0.000001)
    # light_source.set_angular_distribution_gaussian(sigdix=0.000001, sigdiz=0.000001)
    # light_source.set_energy_distribution_singleline(1.000000, unit='A')
    # light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
    # beam = light_source.get_beam()

    # optical element number XX
    boundary_shape = None

    # from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror

