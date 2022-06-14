import numpy

from shadow4.syned.shape import Rectangle

from shadow4.syned.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.mirrors.mirror import Mirror

from shadow4.physical_models.prerefl.prerefl import PreRefl

from shadow4.beamline.s4_beamline_element import S4BeamlineElement


class S4Mirror(Mirror):

    def __init__(self,
                 name="Undefined",
                 boundary_shape=None,
                 surface_shape=None,
                 # inputs related to mirror reflectivity
                 f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                 f_refl=0,  # 0=prerefl file
                 # 1=electric susceptibility
                 # 2=user defined file (1D reflectivity vs angle)
                 # 3=user defined file (1D reflectivity vs energy)
                 # 4=user defined file (2D reflectivity vs energy and angle)
                 file_refl="",  # preprocessor file fir f_refl=0,2,3,4
                 refraction_index=1.0,  # refraction index (complex) for f_refl=1
                 ):
        """

        Parameters
        ----------
        name
        boundary_shape
        element_coordinates_syned
        surface_shape
        f_reflec
                reflectivity of surface: 0=no reflectivity, 1=full polarization
        f_refl
                0=prerefl file
                1=electric susceptibility
                2=user defined file (1D angle in mrad, reflectivity)
                3=user defined file (1D energy in eV, reflectivity)
                4=user defined file (2D energy in eV, angle in mrad, reflectivity)
        file_refl
                name if user defined file
        refraction_index
                complex scalar with refraction index n (for f_refl=1)
        """
        Mirror.__init__(self,
                        name=name,
                        surface_shape=surface_shape,
                        boundary_shape=boundary_shape,
                        coating=None,  # not used
                        coating_thickness=None,  #not used
                        )


        # reflectivity
        self._f_reflec = f_reflec
        self._f_refl = f_refl
        self._file_refl = file_refl
        self._refraction_index = refraction_index

    def set_boundaries_rectangle(self, x_left=-1e3, x_right=1e3, y_bottom=-1e3, y_top=1e3):
        self._boundary_shape = Rectangle(x_left=x_left, x_right=x_right, y_bottom=y_bottom, y_top=y_top)

    def apply_geometrical_model(self, beam):
        raise Exception("To be implemented in the children class")

class S4MirrorElement(S4BeamlineElement):
    
    def __init__(self, optical_element=None, coordinates=None):
        super().__init__(optical_element if optical_element is not None else S4Mirror(),
                         coordinates if coordinates is not None else ElementCoordinates())
    
    def trace_beam(self, beam_in, flag_lost_value=-1):

        p = self.get_coordinates().p()
        q = self.get_coordinates().q()
        theta_grazing1 = numpy.pi / 2 - self.get_coordinates().angle_radial()
        alpha1 = self.get_coordinates().angle_azimuthal()

        #
        beam = beam_in.duplicate()

        #
        # put beam in mirror reference system
        #
        beam.rotate(alpha1, axis=2)
        beam.rotate(theta_grazing1, axis=1)
        beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        #
        # reflect beam in the mirror surface
        #
        soe = self.get_optical_element() #._optical_element_syned

        v_in = beam.get_columns([4,5,6])
        if not isinstance(soe, Mirror): # undefined
            raise Exception("Undefined mirror")
        else:
            mirr, normal = self.apply_local_reflection(beam)

        #
        # apply mirror boundaries
        #
        mirr.apply_boundaries_syned(soe.get_boundary_shape(), flag_lost_value=flag_lost_value)

        #
        # apply mirror reflectivity
        # TODO: add phase
        #

        if soe._f_reflec == 0:
            pass
        elif soe._f_reflec == 1: # full polarization
            v_out = beam.get_columns([4, 5, 6])
            angle_in = numpy.arccos( v_in[0,:] * normal[0,:] +
                                     v_in[1,:] * normal[1,:] +
                                     v_in[2,:] * normal[2,:])

            angle_out = numpy.arccos( v_out[0,:] * normal[0,:] +
                                     v_out[1,:] * normal[1,:] +
                                     v_out[2,:] * normal[2,:])

            grazing_angle_mrad = 1e3 * (numpy.pi / 2 - angle_in)

            if soe._f_refl == 0: # prerefl
                prerefl_file = soe._file_refl
                pr = PreRefl()
                pr.read_preprocessor_file(prerefl_file)
                print(pr.info())

                Rs, Rp, Ru = pr.reflectivity_fresnel(grazing_angle_mrad=grazing_angle_mrad,
                                                     photon_energy_ev=beam.get_column(-11),
                                                     roughness_rms_A=0.0)

            elif soe._f_refl == 1:  # alpha, gamma, electric susceptibilities
                Rs, Rp, Ru = self.reflectivity_fresnel(soe._refraction_index , grazing_angle_mrad=grazing_angle_mrad)

            elif soe._f_refl == 2:  # user angle, mrad ref
                # raise Exception("Not implemented f_refl == 2")

                # values = numpy.loadtxt(self._file_refl)
                #
                # beam_incident_angles = 90.0 - values[:, 1]

                values = numpy.loadtxt(soe._file_refl)

                mirror_grazing_angles = values[:, 0]
                mirror_reflectivities = values[:, 1]

                if mirror_grazing_angles[-1] < mirror_grazing_angles[0]: # XOPPY MLayer gives angles in descendent order
                    mirror_grazing_angles = values[:, 0][::-1]
                    mirror_reflectivities = values[:, 1][::-1]

                # mirror_grazing_angles = numpy.degrees(1e-3*mirror_grazing_angles) # mrad to deg

                Rs = numpy.interp(grazing_angle_mrad,
                                  mirror_grazing_angles,
                                  mirror_reflectivities,
                                  left=mirror_reflectivities[0],
                                  right=mirror_reflectivities[-1])
                Rp = Rs


            elif soe._f_refl == 3:  # user energy

                beam_energies = beam.get_photon_energy_eV()

                values = numpy.loadtxt(soe._file_refl)

                mirror_energies = values[:, 0]
                mirror_reflectivities = values[:, 1]

                Rs = numpy.interp(beam_energies,
                                  mirror_energies,
                                  mirror_reflectivities,
                                  left=mirror_reflectivities[0],
                                  right=mirror_reflectivities[-1])
                Rp = Rs

            elif soe._f_refl == 4:  # user 2D
                values = numpy.loadtxt(soe._file_refl)

                beam_energies = beam.get_photon_energy_eV()

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
                elif values.shape[1] == 4:
                    mirror_reflectivities_s = values[:, 2]
                    mirror_reflectivities_p = values[:, 3]

                    Rs = get_interpolator_weight_2D(mirror_energies, mirror_grazing_angles, mirror_reflectivities_s)
                    Rp = get_interpolator_weight_2D(mirror_energies, mirror_grazing_angles, mirror_reflectivities_p)

            else:
                raise Exception("Not implemented source of mirror reflectivity")

            beam.apply_reflectivities(numpy.sqrt(Rs), numpy.sqrt(Rp))


        #
        # TODO: write angle.xx for comparison
        #


        #
        # from mirror reference system to image plane
        #

        beam_out = mirr.duplicate()
        beam_out.change_to_image_reference_system(theta_grazing1, q)

        return beam_out, mirr

    def apply_local_reflection(self, beam):
        return self.get_optical_element().apply_geometrical_model(beam)

    #
    # i/o utilities
    #
    def set_grazing_angle(self, theta_grazing, theta_azimuthal=None):
        self.get_coordinates()._angle_radial     = numpy.pi / 2 - theta_grazing
        self.get_coordinates()._angle_radial_out = numpy.pi / 2 - theta_grazing
        if theta_azimuthal is not None: self.get_coordinates()._angle_azimuthal = theta_azimuthal

    @classmethod
    def reflectivity_fresnel(cls, refraction_index, grazing_angle_mrad=3.0, debyewaller=1.0):
        theta1 = grazing_angle_mrad * 1e-3     # in rad

        alpha = 2 * (1.0 - refraction_index.real)
        gamma = 2 * refraction_index.imag

        rho = (numpy.sin(theta1))**2 - alpha
        rho += numpy.sqrt((numpy.sin(theta1)**2 - alpha)**2 + gamma**2)
        rho = numpy.sqrt(rho / 2)

        rs1 = 4 * (rho**2) * (numpy.sin(theta1) - rho)**2 + gamma**2
        rs2 = 4 * (rho**2) * (numpy.sin(theta1) + rho)**2 + gamma**2
        rs = rs1 / rs2

        ratio1 = 4 * rho**2 * (rho * numpy.sin(theta1) - numpy.cos(theta1)**2)**2 + gamma**2 * numpy.sin(theta1)**2
        ratio2 = 4 * rho**2 * (rho * numpy.sin(theta1) + numpy.cos(theta1)**2)**2 + gamma**2 * numpy.sin(theta1)**2
        ratio = ratio1 / ratio2

        rp = rs * ratio
        runp = 0.5 * (rs + rp)

        return rs*debyewaller, rp*debyewaller, runp*debyewaller

if __name__ == "__main__":
    pass
