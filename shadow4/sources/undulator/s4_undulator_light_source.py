"""
Undulator light source.

Computes undulator radiation distributions and samples rays according to them.

The radiation divergences (far field) are computed in polar coordinates for a more efiicient sampling.
"""
import numpy

from srxraylib.util.inverse_method_sampler import Sampler1D, Sampler2D, Sampler3D
import scipy.constants as codata
from scipy import interpolate
import scipy.constants as codata

from syned.storage_ring.magnetic_structures.undulator import Undulator
from syned.storage_ring.light_source import LightSource

from shadow4.sources.s4_electron_beam import S4ElectronBeam
from shadow4.beam.s4_beam import S4Beam
from shadow4.sources.undulator.calculate_undulator_emission import calculate_undulator_emission # SourceUndulatorFactory
from shadow4.sources.undulator.calculate_undulator_emission_srw import calculate_undulator_emission_srw # SourceUndulatorFactorySrw
from shadow4.sources.undulator.calculate_undulator_emission_pysru import calculate_undulator_emission_pysru # SourceUndulatorFactoryPysru
from shadow4.sources.undulator.s4_undulator import S4Undulator
from shadow4.sources.s4_light_source import S4LightSource
from shadow4.tools.arrayofvectors import vector_cross, vector_norm

from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource # to get q_a and q_s

from shadow4.tools.logger import is_verbose, is_debug

INTEGRATION_METHOD = 1 # 0=sum, 1=trapz

class S4UndulatorLightSource(S4LightSource):
    """
    Defines an undulator light source and implements the mechanism of sampling rays.

    Parameters
    ----------
    name : str, optional
        The name of the light source.
    electron_beam : instance of ElectronBeam
        The electron beam parameters.
    magnetic_structure : instance of S4BendingMagnet
        The shadow4 bending magnet magnetic structure.
    nrays : int, optional
        The number of rays.
    seed : int, optional
        The Monte Carlo seed.
    """
    def __init__(self,
                 name="Undefined",
                 electron_beam=None,
                 magnetic_structure=None,
                 nrays=5000,
                 seed=12345,
                 ):

        super().__init__(name,
                         electron_beam=electron_beam if not electron_beam is None else S4ElectronBeam(),
                         magnetic_structure=magnetic_structure if not magnetic_structure is None else S4Undulator(),
                         nrays=nrays,
                         seed=seed,
                         )

        # results of calculations
        self.__result_radiation = None
        self.__result_photon_size_distribution = None
        self.__result_photon_size_sigma = None
        # self.__result_photon_size_farfield = None



    def get_beam(self, F_COHER=0):
        """
        Creates the beam as emitted by the undulator.

        Parameters
        ----------
        F_COHER : int, optional
            A flag to indicate that the phase for the s-component is set to zero (coherent_beam=1) or is random for incoherent.

        Returns
        -------
        instance of S4Beam
        """

        return S4Beam.initialize_from_array(self.__calculate_rays(
            user_unit_to_m=1.0, F_COHER=F_COHER))


    def get_resonance_ring(self, harmonic_number=1, ring_order=1):
        """
        Return the angle at which a ring is found.

        Parameters
        ----------
        harmonic_number : int, optional
            The harmonic number.
        ring_order : int, optional
            The ring order (first ring, second ring, etc.).

        Returns
        -------
        float
        """
        return 1.0 / self.get_electron_beam().gamma() * numpy.sqrt(
            ring_order / harmonic_number * (1 + 0.5 * self.get_magnetic_structure().K_vertical()**2))

    def set_energy_monochromatic_at_resonance(self, harmonic_number=1):
        """
        Sets the undulator to a resonance photon energy.

        It changes the K_vertical and MAXANGLE.

        Parameters
        ----------
        harmonic_number : int, optinal
            Harmonic number.
        """
        e = self.get_electron_beam()
        u = self.get_magnetic_structure()

        e_res = u.resonance_energy(e.gamma(), harmonic=harmonic_number)
        maxangle = 3 * 0.69 * u.gaussian_central_cone_aperture(e.gamma(), harmonic_number) # take 3*sigma - _MAXANGLE is in RAD  TODO: why 0.69??

        if is_debug(): print(">>> Setting monochromatic: n=%d, e0=%f eV, maxangle=%f urad" % (harmonic_number, e_res, 1e6 * maxangle))

        u.set_energy_monochromatic(e_res)
        u.set_maxangle(maxangle)

    def set_energy_at_resonance(self, harmonic_number=1, delta_e=0.0):
        """
        Sets the undulator to a resonance photon energy.

        It changes the K_vertical and MAXANGLE.

        Parameters
        ----------
        harmonic_number : int, optinal
            Harmonic number.
        delta_e : float, optinal
            Harmonic number.
        """
        e = self.get_electron_beam()
        u = self.get_magnetic_structure()


        e0 = u.resonance_energy(e.gamma(), harmonic=harmonic_number)
        maxangle = 3 * 0.69 * u.gaussian_central_cone_aperture(e.gamma(), harmonic_number) # take 3*sigma - _MAXANGLE is in RAD  TODO: why 0.69??

        if is_debug(): print(">>> Setting monochromatic: n=%d, e0=%f eV, maxangle=%f urad" % (harmonic_number, e0, 1e6 * maxangle))

        u.set_energy_box(e0 - delta_e / 2, e0 + delta_e / 2)
        u.set_maxangle(3 * 0.69 * u.gaussian_central_cone_aperture(e.gamma(), harmonic_number)) # take 3*sigma - _MAXANGLE is in RAD  TODO: why 0.69??



    #
    # get from results
    #
    def __calculate_photon_size_sigma(self, photon_energy=None):
        u = self.get_magnetic_structure()
        if photon_energy is None: photon_energy = u._EMIN
        lambda1 = codata.h * codata.c / codata.e / photon_energy
        s_phot = 2.740 / (4e0 * numpy.pi) * numpy.sqrt(u.length() * lambda1)
        self.__result_photon_size_sigma = s_phot

    def get_result_dictionary(self):
        """
        Returns a dictionary with plenty of parameters.

        Returns
        -------
        dict
        """
        if self.__result_radiation is None: self.__calculate_radiation()
        return self.__result_radiation

    def get_result_trajectory(self):
        """
        Returns the trajectory arrays y, x, beta_x [in SHADOW frame].

        Returns
        -------
        tuple
            (y, x, beta_x) arrays with coordinates y and x and beta (velocity in c units) along x.
        """
        traj = self.get_result_dictionary()["trajectory"]
        y = traj[3].copy() * codata.c
        x = traj[1].copy() * codata.c
        beta_x = traj[4].copy()
        return y, x, beta_x

    def get_power_density(self):
        """
        Returns the power density of the radiation stack (integrated oner photon energies).

        Returns
        -------
        tuple
            (power_density, theta, phi).
        """

        dict_results = self.get_result_dictionary()
        radiation = dict_results['radiation'].copy()
        photon_energy = dict_results['photon_energy']
        theta = dict_results['theta']
        phi = dict_results['phi']

        if self.get_magnetic_structure().is_monochromatic():
            step_e = 1.0
            power_density = radiation[0] * step_e * codata.e * 1e3  # W/rad2/eV
        else:
            step_e = photon_energy[1] - photon_energy[0]

            for i in range(radiation.shape[0]):
                radiation[i] *= 1e-3 * photon_energy[i] # photons/eV/rad2 -> photons/0.1%bw/rad2

            if INTEGRATION_METHOD == 0:
                power_density = radiation.sum(axis=0) * step_e * codata.e * 1e3 # W/rad2
            else:
                power_density = numpy.trapz(radiation, photon_energy, axis=0) * codata.e * 1e3 # W/rad2

        return power_density, theta, phi


    def get_power_density_interpolated_cartesian(self, npointsx=100, npointsz=100, thetamax=None):
        """
        Returns the power density interpolated to cartesian coordinates.

        Parameters
        ----------
        npointsx : int, optional
            The number of points in X.
        npointsz : int, optional
            The number of points in Z.
        thetamax : None or float, optional
            Maximum value of theta. By default (None) it uses the maximum theta.

        Returns
        -------
        tuple
            (power density, array_x, array_y) in units W/rad2.

        """

        power_density_polar, theta, phi = self.get_power_density()

        if thetamax is None:
            thetamax = theta.max()

        vx = numpy.linspace(-1.1 * thetamax, 1.1 * thetamax, npointsx)
        vz = numpy.linspace(-1.1 * thetamax, 1.1 * thetamax, npointsz)
        VX = numpy.outer(vx, numpy.ones_like(vz))
        VZ = numpy.outer(numpy.ones_like(vx), vz)
        VY = numpy.sqrt(1 - VX**2 - VZ**2)

        THETA = numpy.abs(numpy.arctan( numpy.sqrt(VX**2 + VZ**2) / VY))
        PHI = numpy.arctan2(numpy.abs(VZ), numpy.abs(VX))

        interpolator_value = interpolate.RectBivariateSpline(theta, phi, power_density_polar)
        power_density_cartesian = interpolator_value.ev(THETA, PHI)

        return power_density_cartesian, vx, vz

    def get_flux_and_spectral_power(self):
        """
        Return the integrated flux (photons/s/0.1%bw) and spectral power (W/eV) versus photon energy.

        Returns
        -------
        tuple
            (flux, spectral_power, photon_energy).
        """
        dict_results = self.get_result_dictionary()
        radiation2 = dict_results['radiation'].copy()
        theta = dict_results['theta']
        phi = dict_results['phi']
        photon_energy = dict_results['photon_energy']

        THETA = numpy.outer(theta, numpy.ones_like(phi))

        for i in range(radiation2.shape[0]):
            radiation2[i] *= THETA

        if INTEGRATION_METHOD == 0:
            flux = radiation2.sum(axis=2).sum(axis=1) * (1e-3 * photon_energy) # photons/eV -> photons/0.1%bw
            flux *= 4 * (theta[1] - theta[0]) * (phi[1] - phi[0]) # adding the four quadrants!
        else:
            flux = 4 * numpy.trapz(numpy.trapz(radiation2, phi, axis=2), theta, axis=1) * (1e-3 * photon_energy) # photons/eV -> photons/0.1%bw


        spectral_power = flux * codata.e * 1e3

        return flux, spectral_power, photon_energy

    def calculate_spectrum(self):
        """
        Calculates the spectrum.

        Returns
        -------
        tuple
            (e, f, w) numpy arrays with photon energy (e), photon flux (f) and spectral power (w).
        """
        flux, spectral_power, photon_energy = self.get_flux_and_spectral_power()
        return photon_energy, flux, spectral_power

    def get_flux(self):
        """
        Return the integrated flux (photons/s/0.1%bw) versus photon energy.

        Returns
        -------
        tuple
            (flux, photon_energy).
        """
        photon_energy, flux, spectral_power  = self.calculate_spectrum()
        return flux, photon_energy

    def get_spectral_power(self):
        """
        Return the integrated spectral power (W/eV) versus photon energy.

        Returns
        -------
        tuple
            (spectral_power, photon_energy).
        """
        photon_energy, flux, spectral_power = self.calculate_spectrum()
        return spectral_power, photon_energy

    def get_info(self, debug=False): # todo: clean and merge with S4Undulator.get_info()
        """
        Returns the specific information for the wiggler light source.

        Returns
        -------
        str
        """
        syned_electron_beam = self.get_electron_beam()
        undulator = self.get_magnetic_structure()

        txt = ""
        txt += "-----------------------------------------------------\n"

        txt += "Input Electron parameters: \n"
        txt += "        Electron energy: %f geV\n"%syned_electron_beam._energy_in_GeV
        txt += "        Electron current: %f A\n"%syned_electron_beam._current
        if undulator.get_flag_emittance():
            sigmas = syned_electron_beam.get_sigmas_all()
            txt += "        Electron sigmaX: %g um\n"%(1e6*sigmas[0])
            txt += "        Electron sigmaZ: %g um\n"%(1e6*sigmas[2])
            txt += "        Electron sigmaX': %f urad\n"%(1e6*sigmas[1])
            txt += "        Electron sigmaZ': %f urad\n"%(1e6*sigmas[3])

        txt += "Lorentz factor (gamma): %f\n"%syned_electron_beam.gamma()
        txt += "Electron velocity: %.12f c units\n"%(numpy.sqrt(1.0 - 1.0 / syned_electron_beam.gamma() ** 2))

        txt += "\n" + undulator.get_info()

        txt += "\nResonances: \n"
        txt += "        harmonic number [n]                   %10d %10d %10d \n"%(1,3,5)
        txt += "        wavelength [A]:                       %10.6f %10.6f %10.6f   \n"%(\
                                                                1e10*undulator.resonance_wavelength(syned_electron_beam.gamma(),harmonic=1),
                                                                1e10*undulator.resonance_wavelength(syned_electron_beam.gamma(),harmonic=3),
                                                                1e10*undulator.resonance_wavelength(syned_electron_beam.gamma(),harmonic=5))
        txt += "        energy [eV]   :                       %10.3f %10.3f %10.3f   \n"%(\
                                                                undulator.resonance_energy(syned_electron_beam.gamma(),harmonic=1),
                                                                undulator.resonance_energy(syned_electron_beam.gamma(),harmonic=3),
                                                                undulator.resonance_energy(syned_electron_beam.gamma(),harmonic=5))
        txt += "        frequency [Hz]:                       %10.3g %10.3g %10.3g   \n"%(\
                                                                1e10*undulator.resonance_frequency(syned_electron_beam.gamma(),harmonic=1),
                                                                1e10*undulator.resonance_frequency(syned_electron_beam.gamma(),harmonic=3),
                                                                1e10*undulator.resonance_frequency(syned_electron_beam.gamma(),harmonic=5))
        txt += "        central cone 'half' width [urad]:     %10.6f %10.6f %10.6f   \n"%(\
                                                                1e6*undulator.gaussian_central_cone_aperture(syned_electron_beam.gamma(),1),
                                                                1e6*undulator.gaussian_central_cone_aperture(syned_electron_beam.gamma(),3),
                                                                1e6*undulator.gaussian_central_cone_aperture(syned_electron_beam.gamma(),5))
        txt += "        first ring at [urad]:                 %10.6f %10.6f %10.6f   \n"%(\
                                                                1e6*self.get_resonance_ring(1,1),
                                                                1e6*self.get_resonance_ring(3,1),
                                                                1e6*self.get_resonance_ring(5,1))

        if self.__result_radiation is None:
            txt += "\nradiation: NOT YET CALCULATED\n"
        else:
            txt += "\nradiation: CALCULATED\n"

        return txt

    def to_python_code(self, **kwargs):
        """
        Returns the python code for calculating the wiggler source.

        Returns
        -------
        str
            The python code.
        """
        script = ''
        try:
            script += self.get_electron_beam().to_python_code()
        except:
            script += "\n\n#Error retrieving electron_beam code"

        try:
            script += self.get_magnetic_structure().to_python_code()
        except:
            script += "\n\n#Error retrieving magnetic structure code"


        script += "\n\n\n# light source\nfrom shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource"
        script += "\nlight_source = S4UndulatorLightSource(name='%s', electron_beam=electron_beam, magnetic_structure=source,nrays=%s,seed=%s)" % \
                                                          (self.get_name(),self.get_nrays(),self.get_seed())

        script += "\nbeam = light_source.get_beam()"

        return script

    #
    # internal functions
    #
    def __calculate_radiation(self):
        """
        Calculates the radiation (emission) as a function of theta (elevation angle) and phi (azimuthal angle)
        This radiation will be sampled to create the source

        It calls calculate_undulator_emission*

        Sets the results in self.__result_radiation dictionary.
        """
        syned_electron_beam = self.get_electron_beam()
        undulator = self.get_magnetic_structure()

        self.__result_radiation = None
        if undulator.code_undul_phot == 'internal':
            undul_phot_dict = calculate_undulator_emission(
                electron_energy                  = syned_electron_beam.energy(),
                electron_current                 = syned_electron_beam.current(),
                undulator_period                 = undulator.period_length(),
                undulator_nperiods               = undulator.number_of_periods(),
                K                                = undulator.K(),
                photon_energy                    = undulator._emin,
                EMAX                             = undulator._emax,
                NG_E                             = undulator._ng_e,
                MAXANGLE                         = undulator._maxangle,
                number_of_points                 = undulator._ng_t,
                NG_P                             = undulator._ng_p,
                number_of_trajectory_points      = undulator._ng_j,
                flag_size                        = undulator._flag_size,
                distance                         = undulator._distance,
                magnification                    = undulator._magnification,
                flag_backprop_recalculate_source = undulator._flag_backprop_recalculate_source,
                flag_backprop_weight             = undulator._flag_backprop_weight,
                weight_ratio                     = undulator._weight_ratio,
            )
        elif undulator.code_undul_phot == 'pysru' or  undulator.code_undul_phot == 'pySRU':
            undul_phot_dict = calculate_undulator_emission_pysru(
                electron_energy                  = syned_electron_beam.energy(),
                electron_current                 = syned_electron_beam.current(),
                undulator_period                 = undulator.period_length(),
                undulator_nperiods               = undulator.number_of_periods(),
                K                                = undulator.K(),
                photon_energy                    = undulator._emin,
                EMAX                             = undulator._emax,
                NG_E                             = undulator._ng_e,
                MAXANGLE                         = undulator._maxangle,
                number_of_points                 = undulator._ng_t,
                NG_P                             = undulator._ng_p,
                number_of_trajectory_points      = int(undulator._ng_j),
                flag_size                        = undulator._flag_size,
                distance                         = undulator._distance,
                magnification                    = undulator._magnification,
                flag_backprop_recalculate_source = undulator._flag_backprop_recalculate_source,
                flag_backprop_weight             = undulator._flag_backprop_weight,
                weight_ratio                     = undulator._weight_ratio,
                )
        elif undulator.code_undul_phot == 'srw' or  undulator.code_undul_phot == 'SRW':
            undul_phot_dict = calculate_undulator_emission_srw(
                electron_energy             = syned_electron_beam.energy(),
                electron_current            = syned_electron_beam.current(),
                undulator_period            = undulator.period_length(),
                undulator_nperiods          = undulator.number_of_periods(),
                K                           = undulator.K(),
                photon_energy               = undulator._emin,
                EMAX                        = undulator._emax,
                NG_E                        = undulator._ng_e,
                MAXANGLE                    = undulator._maxangle,
                number_of_points            = undulator._ng_t,
                NG_P                        = undulator._ng_p,
                number_of_trajectory_points = undulator._ng_j,
                flag_size                   = undulator._flag_size,
                distance                    = undulator._distance,
                srw_range                   = undulator._srw_range,
                srw_resolution              = undulator._srw_resolution,
                srw_semianalytical          = undulator._srw_semianalytical,
                )
        else:
            raise Exception("Not implemented undul_phot code: "+undulator.code_undul_phot)

        # add some info
        undul_phot_dict["code_undul_phot"] = undulator.code_undul_phot
        undul_phot_dict["info"] = self.info()

        self.__result_radiation = undul_phot_dict

    def __calculate_photon_size_distribution(self, sampled_photon_energy):
        # calculate (and stores) sizes of the photon undulator beam
        undulator = self.get_magnetic_structure()
        NRAYS = self.get_nrays()

        # see formulas 25 & 30 in Elleaume (Onaki & Elleaume)
        # sp_phot = 0.69*numpy.sqrt(lambda1/undulator_length) # todo: not needed now?
        lambda1 = codata.h * codata.c / codata.e / numpy.array(sampled_photon_energy).mean()
        s_phot = 2.740 / (4e0 * numpy.pi) * numpy.sqrt(undulator.length() * lambda1)

        dict_results = self.get_result_dictionary()

        # electron spread correction
        if undulator._flag_energy_spread and undulator._flag_size > 0:
            e = self.get_electron_beam()
            harmonic_number = int(sampled_photon_energy.mean() / undulator.resonance_energy(e.gamma(), harmonic=1))
            x = 2 * numpy.pi * harmonic_number * undulator.number_of_periods() * e._energy_spread
            q_s = S4UndulatorGaussianLightSource.q_s(x, factor=0.5)
            if is_debug():
                print(">>>>> energy spread correction x: %f; Qs(x)=%f" % (x, q_s))
                print("      n=%d, N=%f, sigma_delta=%g; Normalized energy spread = %g" %
                      (harmonic_number, undulator.number_of_periods(), e._energy_spread, x))
        else:
            q_s = 1.0
            if is_debug(): print(">>>>> NO energy spread correction")


        if undulator._flag_size == 0:
            x_photon = 0.0
            y_photon = 0.0
            z_photon = 0.0
            # for plot, a delta
            x = numpy.linspace(-1e-6, 1e-6, 101)
            y = numpy.zeros_like(x)
            y[y.size // 2] = 1.0

            # for plots, keep this result
            dict_results['BACKPROPAGATED_r']         = x
            dict_results['BACKPROPAGATED_radiation'] = numpy.outer(
                numpy.ones(undulator.get_number_of_energy_points()), y
                )
        elif undulator._flag_size == 1:
            # I added this correction to obtain the sigma in the RADIAL coordinate, not in x and z.
            # TODO: To be verified!
            self.__result_photon_size_sigma = s_phot
            s_phot_corrected = s_phot / numpy.sqrt(2)

            # now apply the electron spread correction
            s_phot_corrected *= q_s

            cov = [[s_phot_corrected**2, 0], [0, s_phot_corrected**2]]
            mean = [0.0,0.0]

            tmp = numpy.random.multivariate_normal(mean, cov, NRAYS)
            x_photon = tmp[:,0]
            y_photon = 0.0
            z_photon = tmp[:,1]

            # for plot, a Gaussian
            x = numpy.linspace(-5 * s_phot, 5 * s_phot, 101)
            y = numpy.exp(-x**2 / 2 / s_phot**2)


            # for plots, keep this result
            dict_results['BACKPROPAGATED_r']         = x
            dict_results['BACKPROPAGATED_radiation'] = numpy.outer(
                numpy.ones(undulator.get_number_of_energy_points()), y
                )

        elif undulator._flag_size == 2:
            # we need to retrieve the emission as a function of the angle
            # radiation, photon_energy, theta, phi
            radiation         = dict_results['radiation']
            photon_energy     = dict_results['photon_energy']
            theta             = dict_results['theta']
            phi               = dict_results['phi']
            e_amplitude_sigma = dict_results['e_amplitude_sigma']
            e_amplitude_pi    = dict_results['e_amplitude_pi']

            #
            # we propagated the emission at a long distance (fat field) to sample angles
            # Then it was propagated back to the source plane to sample the size.
            #
            if undulator.code_undul_phot == 'internal': # wofry1D
                dict1 = self.__result_radiation
                xx = dict1['BACKPROPAGATED_r']

                # now apply the electron spread correction
                xx *= q_s

                if self.get_magnetic_structure().is_monochromatic():
                    yy = dict1['BACKPROPAGATED_radiation'][0]  # todo something
                else:
                    yy = dict1['BACKPROPAGATED_radiation'].sum(axis=0)  # todo something

                sampler_radial = Sampler1D(yy * numpy.abs(xx), xx)
                r, hy, hx = sampler_radial.get_n_sampled_points_and_histogram(NRAYS, bins=101)
                angle = numpy.random.random(NRAYS) * 2 * numpy.pi

                x_photon = r / numpy.sqrt(2.0) * numpy.sin(angle)
                y_photon = 0.0
                z_photon = r / numpy.sqrt(2.0) * numpy.cos(angle)
            elif undulator.code_undul_phot == 'pysru':
                dict1 = self.__result_radiation
                if self.get_magnetic_structure().is_monochromatic():
                    i_prop = dict1['CART_BACKPROPAGATED_radiation'][0]  # todo something
                else:
                    i_prop = dict1['CART_BACKPROPAGATED_radiation'].sum(axis=0)  # todo something
                x = dict1['CART_BACKPROPAGATED_x']
                y = dict1['CART_BACKPROPAGATED_y']

                # now apply the electron spread correction
                x *= q_s / numpy.sqrt(2)
                y *= q_s / numpy.sqrt(2)

                s2d = Sampler2D(i_prop, x, y)
                sampled_x, sampled_z = s2d.get_n_sampled_points(NRAYS)

                x_photon = sampled_x
                y_photon = 0.0
                z_photon = sampled_z
            elif undulator.code_undul_phot == 'srw':
                dict1 = self.__result_radiation
                if self.get_magnetic_structure().is_monochromatic():
                    i_prop = dict1['CART_BACKPROPAGATED_radiation'][0] # todo something
                else:
                    i_prop = dict1['CART_BACKPROPAGATED_radiation'].sum(axis=0) # todo something
                x = dict1['CART_BACKPROPAGATED_x']
                y = dict1['CART_BACKPROPAGATED_y']

                # now apply the electron spread correction
                x *= q_s / numpy.sqrt(2)
                y *= q_s / numpy.sqrt(2)

                s2d = Sampler2D(i_prop, x, y)
                sampled_x, sampled_z = s2d.get_n_sampled_points(NRAYS)

                x_photon = sampled_x
                y_photon = 0.0
                z_photon = sampled_z

        return x_photon, y_photon, z_photon


    def __calculate_rays(self, user_unit_to_m=1.0, F_COHER=0):
        """
        compute the rays in SHADOW matrix (shape (npoints,18) )

        Parameters
        ----------
        user_unit_to_m
        F_COHER: float, optional
            set this flag for coherent beam.

        Returns
        -------
        numpy array
            rays matrix, a numpy.array((npoits,18))
        """

        if self.__result_radiation is None:
            self.__calculate_radiation()

        syned_electron_beam = self.get_electron_beam()
        undulator = self.get_magnetic_structure()
        sigmas = syned_electron_beam.get_sigmas_all()

        NRAYS = self.get_nrays()
        rays = numpy.zeros((NRAYS, 18))
        if self.get_seed() != 0: numpy.random.seed(self.get_seed())

        #
        # sample energies (col 11), theta and phi
        #
        sampled_photon_energy, sampled_theta, sampled_phi = self._sample_photon_energy_theta_and_phi(NRAYS)
        A2EV = 2.0 * numpy.pi / (codata.h * codata.c / codata.e*1e2)

        rays[:, 10] =  sampled_photon_energy * A2EV

        #
        # sample sizes (cols 1-3)
        #
        if undulator.get_flag_emittance():
            x_electron = numpy.random.normal(loc=0.0, scale=sigmas[0], size=NRAYS)
            y_electron = 0.0
            z_electron = numpy.random.normal(loc=0.0, scale=sigmas[2], size=NRAYS)
        else:
            x_electron = 0.0
            y_electron = 0.0
            z_electron = 0.0

        x_photon, y_photon, z_photon = self.__calculate_photon_size_distribution(sampled_photon_energy)

        rays[:, 0] = x_photon + x_electron
        rays[:, 1] = y_photon + y_electron
        rays[:, 2] = z_photon + z_electron

        if user_unit_to_m != 1.0:
            rays[:, 0] /= user_unit_to_m
            rays[:, 1] /= user_unit_to_m
            rays[:, 2] /= user_unit_to_m

        #
        # sample divergences (cols 4-6): the Shadow way
        #
        THETABM = sampled_theta
        PHI = sampled_phi
        A_Z = numpy.arcsin(numpy.sin(THETABM) * numpy.sin(PHI))
        A_X = numpy.arccos(numpy.cos(THETABM) / numpy.cos(A_Z))
        THETABM = A_Z
        PHI  = A_X
        # ! C Decide in which quadrant THETA and PHI are.
        myrand = numpy.random.random(NRAYS)
        THETABM[numpy.where(myrand < 0.5)] *= -1.0
        myrand = numpy.random.random(NRAYS)
        PHI[numpy.where(myrand < 0.5)] *= -1.0

        if undulator.get_flag_emittance():
            EBEAM1 = numpy.random.normal(loc=0.0, scale=sigmas[1], size=NRAYS)
            EBEAM3 = numpy.random.normal(loc=0.0, scale=sigmas[3], size=NRAYS)
            ANGLEX = EBEAM1 + PHI
            ANGLEV = EBEAM3 + THETABM
        else:
            ANGLEX = PHI # E_BEAM(1) + PHI
            ANGLEV = THETABM #  E_BEAM(3) + THETABM

        VX = numpy.tan(ANGLEX)
        VY = 1.0
        VZ = numpy.tan(ANGLEV) / numpy.cos(ANGLEX)
        VN = numpy.sqrt( VX * VX + VY * VY + VZ * VZ)
        VX /= VN
        VY /= VN
        VZ /= VN

        rays[:,3] = VX
        rays[:,4] = VY
        rays[:,5] = VZ

        #
        # electric field vectors (cols 7-9, 16-18) and phases (cols 14-15)
        #

        # ! C
        # ! C  ---------------------------------------------------------------------
        # ! C                 POLARIZATION
        # ! C
        # ! C   Generates the polarization of the ray. This is defined on the
        # ! C   source plane, so that A_VEC is along the X-axis and AP_VEC is along Z-axis.
        # ! C   Then care must be taken so that A will be perpendicular to the ray
        # ! C   direction.
        # ! C
        # ! C
        # A_VEC = [1,0,0]
        DIREC = rays[:, 3:6].copy()
        A_VEC = numpy.zeros_like(DIREC)
        A_VEC[:, 0] = 1.0

        # ! C
        # ! C   Rotate A_VEC so that it will be perpendicular to DIREC and with the
        # ! C   right components on the plane.
        # ! C
        A_TEMP = vector_cross(A_VEC, DIREC)
        A_VEC = vector_cross(DIREC, A_TEMP)
        A_VEC = vector_norm(A_VEC)
        AP_VEC = vector_cross(A_VEC, DIREC)
        AP_VEC = vector_norm(AP_VEC)

        #
        # obtain polarization for each ray (interpolation)
        #
        if undulator._ng_e == 1: # 2D interpolation
            sampled_photon_energy = numpy.array(sampled_photon_energy) # be sure is an array
            fn = interpolate.RegularGridInterpolator(
                (self.__result_radiation["theta"], self.__result_radiation["phi"]),
                self.__result_radiation["polarization"][0])

            pts = numpy.dstack( (sampled_theta, sampled_phi) )
            pts = pts[0]
            POL_DEG = fn(pts)
        else: # 3D interpolation
            fn = interpolate.RegularGridInterpolator(
                (self.__result_radiation["photon_energy"], self.__result_radiation["theta"], self.__result_radiation["phi"]),
                self.__result_radiation["polarization"])

            pts = numpy.dstack( (sampled_photon_energy,
                                 sampled_theta,
                                 sampled_phi) )
            pts = pts[0]
            POL_DEG = fn(pts)

        #     ! C
        #     ! C   WaNT A**2 = AX**2 + AZ**2 = 1 , instead of A_VEC**2 = 1 .
        #     ! C

        DENOM = numpy.sqrt(1.0 - 2.0 * POL_DEG + 2.0 * POL_DEG**2)
        AX = POL_DEG/DENOM
        for i in range(3):
            A_VEC[:, i] *= AX

        AZ = (1.0 - POL_DEG) / DENOM
        for i in range(3):
            AP_VEC[:, i] *= AZ

        rays[:, 6:9] =  A_VEC
        rays[:, 15:18] = AP_VEC

        # ! C
        # ! C Now the phases of A_VEC and AP_VEC.
        # ! C

        POL_ANGLE = 0.5 * numpy.pi

        if F_COHER == 1:
            PHASEX = 0.0
        else:
            PHASEX = numpy.random.random(NRAYS) * 2 * numpy.pi

        PHASEZ = PHASEX + POL_ANGLE * numpy.sign(ANGLEV)

        rays[:, 13] = PHASEX
        rays[:, 14] = PHASEZ

        # set flag (col 10)
        rays[:, 9] = 1.0



        # col 12 (ray index)
        rays[:, 11] =  1 + numpy.arange(NRAYS)

        # col 13 (optical path)
        rays[:, 12] = 0.0

        return rays

    def _sample_photon_energy_theta_and_phi(self, NRAYS):

        #
        # sample divergences
        #
        theta         = self.__result_radiation["theta"]
        phi           = self.__result_radiation["phi"]
        photon_energy = self.__result_radiation["photon_energy"]

        undulator = self.get_magnetic_structure()




        if undulator.is_monochromatic():
            # energy spread correction factor (for the moment for monochromatic only)
            if undulator._flag_energy_spread:
                e = self.get_electron_beam()
                harmonic_number = int(photon_energy.mean() / undulator.resonance_energy(e.gamma(), harmonic=1))
                x = 2 * numpy.pi * harmonic_number * undulator.number_of_periods() * e._energy_spread
                q_a = S4UndulatorGaussianLightSource.q_a(x)
                if is_debug():
                    print(">>>>> energy spread correction x: %f; Qa(x)=%f" % (x, q_a))
                    print("      n=%d, N=%f, sigma_delta=%g; Normalized energy spread = %g" %
                          (harmonic_number, undulator.number_of_periods(), e._energy_spread, x))
            else:
                q_a = 1.0
                if is_debug(): print(">>>>> NO energy spread correction")

            #2D case
            tmp = self.__result_radiation["radiation"][0, :, :].copy()
            tmp /= tmp.max()
            # correct radiation for DxDz / DthetaDphi
            tmp_theta = numpy.outer(theta, numpy.ones_like(phi))
            tmp_theta /= tmp_theta.max()
            tmp_theta += 1e-6 # to avoid zeros
            tmp *= tmp_theta
            # plot_image(tmp_theta,theta,phi,aspect='auto')

            theta *= q_a # apply energy spread correction

            s2d = Sampler2D(tmp, theta, phi)
            sampled_theta, sampled_phi = s2d.get_n_sampled_points(NRAYS)
            sampled_photon_energy = numpy.ones(NRAYS) * self.get_magnetic_structure()._emin
        else:
            # energy spread correction factor (for the moment for monochromatic only)
            if undulator._flag_energy_spread:
                raise Exception("** Error ** Energy spread cannot be calculated for polychromatic sources")

            #3D case
            tmp = self.__result_radiation["radiation"].copy()
            tmp /= tmp.max()
            # correct radiation for DxDz / DthetaDphi
            tmp_theta = numpy.outer(theta, numpy.ones_like(phi))
            tmp_theta /= tmp_theta.max()
            tmp_theta += 1e-6 # to avoid zeros
            for i in range(tmp.shape[0]):
                tmp[i,:,:] *= tmp_theta
            s3d = Sampler3D(tmp, photon_energy, theta, phi)
            sampled_photon_energy,sampled_theta,sampled_phi = s3d.get_n_sampled_points(NRAYS)

        return sampled_photon_energy, sampled_theta, sampled_phi



if __name__ == "__main__":

    import numpy

    do_plots = True

    su = S4Undulator(K_vertical=0.25,
                     period_length=0.032,
                     number_of_periods=50,
                     code_undul_phot='internal',
                     emin=10000.0,
                     emax=15000.0,
                     flag_energy_spread=0,
                     )

    ebeam = S4ElectronBeam(energy_in_GeV=6.04,
                 energy_spread = 0.0,
                 current = 0.2,
                 number_of_bunches = 0,
                 moment_xx=(400e-6)**2,
                 moment_xxp=0.0,
                 moment_xpxp=(10e-6)**2,
                 moment_yy=(10e-6)**2,
                 moment_yyp=0.0,
                 moment_ypyp=(4e-6)**2 )

    ls = S4UndulatorLightSource(name="", electron_beam=ebeam, magnetic_structure=su)

    # beam =  ls.get_beam()
    # #
    # # print(ls.info())
    # # # print(ls.to_python_code())
    # #
    # ls.set_energy_monochromatic_at_resonance(harmonic_number=1)
    # # print(ls.info())
    #
    # print("dict: ")
    # d = ls.get_result_dictionary()
    # for key in d.keys():
    #     print(key)
    #
    # # print("radiation: ", d["radiation"].shape)
    # # print("polarization: ", d["polarization"].shape)
    # # print("phi: ", d["phi"].shape)
    # # print("theta: ", d["theta"].shape)
    # traj = d['trajectory']
    # print(">>> traj ", traj.shape)
    # from srxraylib.plot.gol import plot
    # import scipy.constants as codata
    # plot(traj[3] * codata.c, traj[1] * codata.c, xtitle="Z", ytitle="X")
    # plot(traj[3] * codata.c, traj[4], xtitle="Z", ytitle="beta_x")

    print(ls.get_info())

