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
from shadow4.sources.undulator.source_undulator_factory import calculate_undulator_emission # SourceUndulatorFactory
from shadow4.sources.undulator.source_undulator_factory_srw import calculate_undulator_emission_SRW # SourceUndulatorFactorySrw
from shadow4.sources.undulator.source_undulator_factory_pysru import calculate_undulator_emission_pySRU # SourceUndulatorFactoryPysru
from shadow4.sources.undulator.s4_undulator import S4Undulator
from shadow4.sources.s4_light_source import S4LightSource
from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian
from shadow4.tools.arrayofvectors import vector_cross, vector_norm

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
        self.__result_photon_size_farfield = None



    def get_beam(self, F_COHER=0, verbose=1):
        """
        Creates the beam as emitted by the undulator.

        Parameters
        ----------
        F_COHER : int, optional
            A flag to indicate that the phase for the s-component is set to zero (coherent_beam=1) or is random for incoherent.
        verbose : int, optional
            Set to 1 for verbose output.

        Returns
        -------
        instance of S4Beam
        """
        if self.get_magnetic_structure().use_gaussian_approximation():
            return self.get_beam_in_gaussian_approximation()
        else:
            return S4Beam.initialize_from_array(self.__calculate_rays(
                user_unit_to_m=1.0, F_COHER=F_COHER, verbose=verbose))


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
        u.set_energy_monochromatic(
            u.resonance_energy(e.gamma(), harmonic=harmonic_number)
        )
        # take 3*sigma - _MAXANGLE is in RAD  TODO: why 0.69??
        u.set_maxangle(3 * 0.69 * u.gaussian_central_cone_aperture(e.gamma(), harmonic_number))

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
        u.set_energy_box(e0 - delta_e / 2, e0 + delta_e / 2)
        # take 3*sigma - _MAXANGLE is in RAD  TODO: why 0.69??
        u.set_maxangle(3 * 0.69 * u.gaussian_central_cone_aperture(e.gamma(), harmonic_number))



    #
    # get from results
    #
    def __calculate_photon_size_sigma(self, photon_energy=None):
        u = self.get_magnetic_structure()
        if photon_energy is None: photon_energy = u._EMIN
        lambda1 = codata.h * codata.c / codata.e / photon_energy
        s_phot = 2.740 / (4e0 * numpy.pi) * numpy.sqrt(u.length() * lambda1)
        self.__result_photon_size_sigma = s_phot

    def get_result_photon_size_sigma(self, photon_energy=None):
        """
        Returns the photon size in meters using Gaussian approximation at photon energy emin.

        Parameters
        ----------
        photon_energy : None or float, optional
            The phoron energy. If None it uses undulator emin.

        Returns
        -------
        float
        """
        if self.__result_photon_size_sigma is None: self.__calculate_photon_size_sigma(photon_energy=photon_energy)
        return self.__result_photon_size_sigma

    def get_result_dictionary(self):
        """
        Returns a dictionary with plenty of parameters.

        Returns
        -------
        dict
        """
        if self.__result_radiation is None: self.__calculate_radiation()
        return self.__result_radiation

    def get_result_radiation(self):
        """
        Returns the array with the radiation (intensity vs angles).

        Returns
        -------
        numpy array
            (1, ng_t, ng_p) shape 1 x  Number of points in angle theta x points in angle phi.
        """
        return self.get_result_dictionary()["radiation"]

    def get_result_polarization(self):
        """
        Returns the array with the polarization degree vs angles.

        Returns
        -------
        numpy array
            (1, ng_t, ng_p) shape 1 x  Number of points in angle theta x points in angle phi.
        """
        return self.get_result_dictionary()["polarization"]

    def get_result_polarisation(self):
        """
        Returns the array with the polarization degree vs angles.

        Returns
        -------
        numpy array
            (1, ng_t, ng_p) shape 1 x  Number of points in angle theta x points in angle phi.
        """
        return self.get_result_dictionary()["polarization"]

    def get_result_theta(self):
        """
        Returns the array with the theta angle.

        Returns
        -------
        numpy array
            (ng_t) size number of points in angle theta.
        """
        return self.get_result_dictionary()["theta"]

    def get_result_phi(self):
        """
        Returns the array with the phi angle.

        Returns
        -------
        numpy array
            (ng_p) size number of points in angle phi.
        """
        return self.get_result_dictionary()["phi"]

    def get_result_e_amplitudes(self):
        """
        Returns the arrays with the electric field amplitudes.

        Returns
        -------
        numpy array
            (ng_e) size number of points in photon energy.
        """
        return self.get_result_dictionary()["e_amplitude_sigma"], self.get_result_dictionary()["e_amplitude_pi"]

    def get_result_photon_energy(self):
        """
        Returns the array with the photon energy.

        Returns
        -------
        numpy array
            (ng_e) size number of points in photon energy.
        """
        return self.get_result_dictionary()["photon_energy"]

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

    def get_result_radiation_polar(self):
        """
        Returns all radiation arrays.

        Returns
        -------
        tuple
            (radiation, photon_energy, theta, phi)
        """
        return self.get_result_radiation(),\
               self.get_result_photon_energy(),\
               self.get_result_theta(),\
               self.get_result_phi()

    def get_radiation(self):
        """
        Returns all radiation arrays (the same as get_result_radiation_polar() ).

        Returns
        -------
        tuple
            (radiation, photon_energy, theta, phi)
        """
        return self.get_result_radiation_polar()

    def get_radiation_interpolated_cartesian(self, npointsx=100, npointsz=100, thetamax=None):
        """
        Interpolates the radiation array (in polar coordinates) to cartesian coordinates.

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
            (radiation, array_x, array_y) in units of W/rad^2.
        """

        radiation, photon_energy, thetabm,phi = self.get_result_radiation_polar()

        if thetamax is None:
            thetamax = thetabm.max()

        vx = numpy.linspace(-1.1 * thetamax, 1.1 * thetamax, npointsx)
        vz = numpy.linspace(-1.1 * thetamax, 1.1 * thetamax, npointsz)
        VX = numpy.outer(vx, numpy.ones_like(vz))
        VZ = numpy.outer(numpy.ones_like(vx), vz)
        VY = numpy.sqrt(1 - VX**2 - VZ**2)

        THETA = numpy.abs(numpy.arctan(numpy.sqrt(VX**2 + VZ**2) / VY))
        PHI = numpy.arctan2(numpy.abs(VZ), numpy.abs(VX))

        radiation_interpolated = numpy.zeros((radiation.shape[0], npointsx, npointsz))

        for i in range(radiation.shape[0]):
            interpolator_value = interpolate.RectBivariateSpline(thetabm, phi, radiation[i])
            radiation_interpolated[i] = interpolator_value.ev(THETA, PHI)

        return radiation_interpolated, photon_energy, vx, vz

    def get_power_density(self):
        """
        Returns the power density of the radiation stack (integrated oner photon energies).

        Returns
        -------
        tuple
            (power_density, theta, phi).
        """
        radiation = self.get_result_radiation().copy()
        theta = self.get_result_theta()
        phi = self.get_result_phi()
        photon_energy = self.get_result_photon_energy()

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
            (power density, array_x, array_y) in units W/m.

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
        radiation2 = self.get_result_radiation().copy()
        theta = self.get_result_theta()
        phi = self.get_result_phi()
        photon_energy = self.get_result_photon_energy()
        THETA = numpy.outer(theta, numpy.ones_like(phi))
        for i in range(radiation2.shape[0]):
            radiation2[i] *= THETA

        if INTEGRATION_METHOD == 0:
            flux = radiation2.sum(axis=2).sum(axis=1) * (1e-3 * photon_energy) # photons/eV -> photons/0.1%bw
            flux *= 4 * (theta[1] - theta[0]) * (phi[1] - phi[0]) # adding the four quadrants!
        else:
            flux = 4 * numpy.trapz(numpy.trapz(radiation2, phi, axis=2), theta, axis=1) * (1e-3 * photon_energy) # photons/eV -> photons/0.1%bw


        spectral_power = flux*codata.e*1e3

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
        flux, spectral_power, photon_energy = self.get_flux_and_spectral_power()
        return flux, photon_energy

    def get_spectral_power(self):
        """
        Return the integrated spectral power (W/eV) versus photon energy.

        Returns
        -------
        tuple
            (spectral_power, photon_energy).
        """
        flux, spectral_power, photon_energy = self.get_flux_and_spectral_power()
        return spectral_power, photon_energy

    def get_photon_size_distribution(self):
        """
        Returns the arrays of 1D photon size distribution.

        Returns
        -------
        tuple
            (array_x, array_z).
        """
        if self.__result_photon_size_distribution is None:
            raise Exception("Not yet calculated...")

        return self.__result_photon_size_distribution["x"], self.__result_photon_size_distribution["y"]

    def get_photon_size_farfield(self):
        """
        Returns the arrays of far field distribution.

        Returns
        -------
        tuple
            (array_x, array_z).
        """
        if self.__result_photon_size_distribution is None:
            raise Exception("Not yet calculated...")

        # theta, radial_flux, mean_photon_energy, distance, magnification
        return self.__result_photon_size_farfield["theta"],\
               self.__result_photon_size_farfield["radial_e_amplitude"], \
               self.__result_photon_size_farfield["mean_photon_energy"], \
               self.__result_photon_size_farfield["distance"], \
               self.__result_photon_size_farfield["magnification"]


    def get_beam_in_gaussian_approximation(self):
        """
        Returns the beam with the sampled rays using for an undulator modelled with the Gaussian approximation.

        Returns
        -------
        instance of S4Beam.
        """

        emin, emax, npoints = self.get_magnetic_structure().get_energy_box()

        NRAYS = self.get_nrays()

        if self.get_magnetic_structure().get_flag_emittance():
            sigma_x, sigdi_x, sigma_z, sigdi_z = self.get_electron_beam().get_sigmas_all()
            Sx, Sz, Spx, Spz = self.get_undulator_photon_beam_sizes_by_convolution(
                sigma_x=sigma_x,
                sigma_z=sigma_z,
                sigdi_x=sigdi_x,
                sigdi_z=sigdi_z,
                undulator_E0=0.5*(emin+emax),
                undulator_length=self.get_magnetic_structure().length(),
                )
            a = SourceGaussian.initialize_from_keywords(
                                                        sigmaX=Sx,
                                                        sigmaY=0,
                                                        sigmaZ=Sz,
                                                        sigmaXprime=Spx,
                                                        sigmaZprime=Spz,
                                                        real_space_center=[0.0, 0.0, 0.0],
                                                        direction_space_center=[0.0, 0.0],
                                                        nrays=NRAYS,
                                                        seed=self.get_seed()
                                                        )
        else:
            s, sp = self.get_undulator_photon_beam_sizes(
                                                        undulator_E0=0.5*(emin+emax),
                                                        undulator_length=self.get_magnetic_structure().length(),
                                                        )
            a = SourceGaussian.initialize_from_keywords(
                                                        sigmaX=s,
                                                        sigmaY=0,
                                                        sigmaZ=s,
                                                        sigmaXprime=sp,
                                                        sigmaZprime=sp,
                                                        real_space_center=[0.0, 0.0, 0.0],
                                                        direction_space_center=[0.0, 0.0],
                                                        nrays=NRAYS,
                                                        seed=self.get_seed()
                                                        )

        beam = a.get_beam()
        if emax == emin:
            e = numpy.zeros(NRAYS) + emin
        else:
            e = numpy.random.random(NRAYS) * (emax - emin) + emin

        beam.set_photon_energy_eV(e)

        return beam

    @classmethod
    def get_undulator_photon_beam_sizes(cls,
                                        undulator_E0=15000.0,
                                        undulator_length=4.0):
        """
        Returns the filament beam undulator size and divergence at resonance using Gaussian approximation.

        Parameters
        ----------
        undulator_E0 : float, optional
            The undulator resonance energy in eV.
        undulator_length : float, optional
            The undulator length in meters.

        Returns
        -------
        tuple
            (sigma, sigma').
        """
        m2ev = codata.c * codata.h / codata.e
        lambda1 = m2ev / undulator_E0
        # calculate sizes of the photon undulator beam
        # see formulas 25 & 30 in Elleaume (Onaki & Elleaume)
        s_phot = 2.740 / (4e0 * numpy.pi) * numpy.sqrt(undulator_length * lambda1)
        sp_phot = 0.69 * numpy.sqrt(lambda1 / undulator_length)
        return s_phot, sp_phot

    @classmethod
    def get_undulator_photon_beam_sizes_by_convolution(cls,
                                                       sigma_x=1e-4,
                                                       sigma_z=1e-4,
                                                       sigdi_x=1e-6,
                                                       sigdi_z=1e-6,
                                                       undulator_E0=15000.0,
                                                       undulator_length=4.0):
        """
        Return the photon beam sizes by convolution of the electron sizes and radiation sizes.

        Parameters
        ----------
        sigma_x : float, optional
            The sigma for X.
        sigma_z : float, optional
            The sigma for Z.
        sigdi_x : float, optional
            The sigma for X'.
        sigdi_z : float, optional
            The sigma for Z'.
        undulator_E0 : float, optional
            The undulator resonance energy in eV.
        undulator_length : float, optional
            The undulator length in meters.

        Returns
        -------
        tuple
            (SigmaX, SigmaZ, SigmaX', SigmaZ')
        """

        user_unit_to_m = 1.0

        s_phot, sp_phot = cls.get_undulator_photon_beam_sizes(undulator_E0=undulator_E0, undulator_length=undulator_length)

        photon_h = numpy.sqrt(numpy.power(sigma_x * user_unit_to_m, 2) + numpy.power(s_phot, 2) )
        photon_v = numpy.sqrt(numpy.power(sigma_z * user_unit_to_m, 2) + numpy.power(s_phot, 2) )
        photon_hp = numpy.sqrt(numpy.power(sigdi_x, 2) + numpy.power(sp_phot, 2) )
        photon_vp = numpy.sqrt(numpy.power(sigdi_z, 2) + numpy.power(sp_phot, 2) )

        m2ev = codata.c * codata.h / codata.e
        lambda1 = m2ev / undulator_E0
        txt = ""
        txt += "Photon single electron emission at wavelength %f A: \n" % (lambda1 * 1e10)
        txt += "    sigma_u:      %g um \n" % (1e6 * s_phot)
        txt += "    sigma_uprime: %g urad \n" % (1e6 * sp_phot)
        txt += "Electron sizes: \n"
        txt += "    sigma_x: %g um \n" % (1e6 * sigma_x * user_unit_to_m)
        txt += "    sigma_z: %g um \n" % (1e6 * sigma_z * user_unit_to_m)
        txt += "    sigma_x': %g urad \n" % (1e6 * sigdi_x)
        txt += "    sigma_z': %g urad \n" % (1e6 * sigdi_z)
        txt += "Photon source sizes (convolution): \n"
        txt += "    Sigma_x: %g um \n" % (1e6 * photon_h)
        txt += "    Sigma_z: %g um \n" % (1e6 * photon_v)
        txt += "    Sigma_x': %g urad \n" % (1e6 * photon_hp)
        txt += "    Sigma_z': %g urad \n" % (1e6 * photon_vp)

        print(txt)

        return (photon_h/user_unit_to_m, photon_v/user_unit_to_m, photon_hp, photon_vp)

    def get_info(self, debug=False):
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
        if undulator._FLAG_EMITTANCE:
            sigmas = syned_electron_beam.get_sigmas_all()
            txt += "        Electron sigmaX: %g [um]\n"%(1e6*sigmas[0])
            txt += "        Electron sigmaZ: %g [um]\n"%(1e6*sigmas[2])
            txt += "        Electron sigmaX': %f urad\n"%(1e6*sigmas[1])
            txt += "        Electron sigmaZ': %f urad\n"%(1e6*sigmas[3])
        txt += "Input Undulator parameters: \n"
        txt += "        period: %f m\n"%undulator.period_length()
        txt += "        number of periods: %d\n"%undulator.number_of_periods()
        txt += "        K-value: %f\n"%undulator.K_vertical()

        txt += "-----------------------------------------------------\n"

        txt += "Lorentz factor (gamma): %f\n"%syned_electron_beam.gamma()
        txt += "Electron velocity: %.12f c units\n"%(numpy.sqrt(1.0 - 1.0 / syned_electron_beam.gamma() ** 2))
        txt += "Undulator length: %f m\n"%(undulator.period_length()*undulator.number_of_periods())
        K_to_B = (2.0 * numpy.pi / undulator.period_length()) * codata.m_e * codata.c / codata.e

        txt += "Undulator peak magnetic field: %f T\n"%(K_to_B*undulator.K_vertical())
        txt += "Resonances: \n"
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

        txt += "-----------------------------------------------------\n"
        txt += "Grids: \n"
        if undulator._NG_E == 1:
            txt += "        photon energy %f eV\n"%(undulator._EMIN)
        else:
            txt += "        photon energy from %10.3f eV to %10.3f eV\n"%(undulator._EMIN,undulator._EMAX)
        txt += "        number of points for the trajectory: %d\n"%(undulator._NG_J)
        txt += "        number of energy points: %d\n"%(undulator._NG_E)
        txt += "        maximum elevation angle: %f urad\n"%(1e6*undulator._MAXANGLE)
        txt += "        number of angular elevation points: %d\n"%(undulator._NG_T)
        txt += "        number of angular azimuthal points: %d\n"%(undulator._NG_P)
        # txt += "        number of rays: %d\n"%(self.NRAYS)
        # txt += "        random seed: %d\n"%(self.SEED)
        txt += "-----------------------------------------------------\n"

        txt += "calculation code: %s\n"%undulator.code_undul_phot
        if self.__result_radiation is None:
            txt += "radiation: NOT YET CALCULATED\n"
        else:
            txt += "radiation: CALCULATED\n"
        txt += "Sampling: \n"
        if undulator._FLAG_SIZE == 0:
            flag = "point"
        elif undulator._FLAG_SIZE == 1:
            flag = "Gaussian"
        elif undulator._FLAG_SIZE == 2:
            flag = "Far field backpropagated"

        txt += "        Photon source size sampling flag: %d (%s)\n"%(undulator._FLAG_SIZE,flag)
        # if undulator._FLAG_SIZE == 1:
        #     if lightdource.get_result_photon_size_sigma() is not None:
        #         txt += "        Photon source size sigma (Gaussian): %6.3f um \n"%(1e6 * self.__result_photon_size_sigma)

        txt += "-----------------------------------------------------\n"
        return txt

    def to_python_code(self, **kwargs):
        """
        returns the python code for calculating the wiggler source.

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

        # script += "\n\n\n#beamline\nfrom shadow4.beamline.s4_beamline import S4Beamline"
        # script += "\nbeamline = S4Beamline(light_source=light_source)"
        script += "\nbeam = light_source.get_beam()"

        return script


    #
    # internal functions
    #
    def __calculate_radiation(self):
        """
        Calculates the radiation (emission) as a function of theta (elevation angle) and phi (azimuthal angle)
        This radiation will be sampled to create the source

        It calls undul_phot* in SourceUndulatorFactory

        :param code_undul_phot: 'internal' (calls undul_phot), 'pysru' (calls undul_phot_pysru) or
                'srw' (calls undul_phot_srw)
        :return: a dictionary (the output from undul_phot*)
        """

        # h = self.to_dictionary()
        # print(self.info())
        # os.system("rm -f xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")

        # if code_undul_phot != "internal" or code_undul_phot != "srw":
        #     dump_uphot_dot_dat = True

        syned_electron_beam = self.get_electron_beam()
        undulator = self.get_magnetic_structure()

        self.__result_radiation = None
        # undul_phot
        if undulator.code_undul_phot == 'internal':
            undul_phot_dict = calculate_undulator_emission(
                electron_energy             = syned_electron_beam.energy(),
                electron_current            = syned_electron_beam.current(),
                undulator_period            = undulator.period_length(),
                undulator_nperiods          = undulator.number_of_periods(),
                K                           = undulator.K(),
                photon_energy               = undulator._EMIN,
                EMAX                        = undulator._EMAX,
                NG_E                        = undulator._NG_E,
                MAXANGLE                    = undulator._MAXANGLE,
                number_of_points            = undulator._NG_T,
                NG_P                        = undulator._NG_P,
                number_of_trajectory_points = undulator._NG_J,
                flag_size                   = undulator._FLAG_SIZE,
                distance                   = undulator._distance,
                magnification              = undulator._magnification,
                )

        elif undulator.code_undul_phot == 'pysru' or  undulator.code_undul_phot == 'pySRU':
            undul_phot_dict = calculate_undulator_emission_pySRU(
                electron_energy            = syned_electron_beam.energy(),
                electron_current           = syned_electron_beam.current(),
                undulator_period           = undulator.period_length(),
                undulator_nperiods         = undulator.number_of_periods(),
                K                          = undulator.K(),
                photon_energy              = undulator._EMIN,
                EMAX                       = undulator._EMAX,
                NG_E                       = undulator._NG_E,
                MAXANGLE                   = undulator._MAXANGLE,
                number_of_points           = undulator._NG_T,
                NG_P                       = undulator._NG_P,
                number_of_trajectory_points= int(undulator._NG_J),
                flag_size                  = undulator._FLAG_SIZE,
                distance                   = undulator._distance,
                magnification              = undulator._magnification,
                pysru_source               = undulator._pysru_source,
                )
        elif undulator.code_undul_phot == 'srw' or  undulator.code_undul_phot == 'SRW':
            undul_phot_dict = calculate_undulator_emission_SRW(
                electron_energy             = syned_electron_beam.energy(),
                electron_current            = syned_electron_beam.current(),
                undulator_period            = undulator.period_length(),
                undulator_nperiods          = undulator.number_of_periods(),
                K                           = undulator.K(),
                photon_energy               = undulator._EMIN,
                EMAX                        = undulator._EMAX,
                NG_E                        = undulator._NG_E,
                MAXANGLE                    = undulator._MAXANGLE,
                number_of_points            = undulator._NG_T,
                NG_P                        = undulator._NG_P,
                number_of_trajectory_points = undulator._NG_J,
                flag_size                   = undulator._FLAG_SIZE,
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


        if undulator._FLAG_SIZE == 0:
            x_photon = 0.0
            y_photon = 0.0
            z_photon = 0.0
            # for plot, a delta
            x = numpy.linspace(-1e-6, 1e-6, 101)
            y = numpy.zeros_like(x)
            y[y.size // 2] = 1.0
            self.__result_photon_size_distribution = {"x":x, "y":y}

        elif undulator._FLAG_SIZE == 1:
            # TODO: I added this correction to obtain the sigma in the RADIAL coordinate, not in x and z.
            # TODO: TO be verified!

            self.__result_photon_size_sigma = s_phot
            s_phot_corrected = s_phot / numpy.sqrt(2)

            cov = [[s_phot_corrected**2, 0], [0, s_phot_corrected**2]]
            mean = [0.0,0.0]

            tmp = numpy.random.multivariate_normal(mean, cov, NRAYS)
            x_photon = tmp[:,0]
            y_photon = 0.0
            z_photon = tmp[:,1]

            # for plot, a Gaussian
            x = numpy.linspace(-5 * s_phot, 5 * s_phot, 101)
            y = numpy.exp(-x**2 / 2 / s_phot**2)
            self.__result_photon_size_distribution = {"x":x, "y":y}

        elif undulator._FLAG_SIZE == 2:
            # we need to retrieve the emission as a function of the angle
            radiation, photon_energy, theta, phi = self.get_result_radiation_polar()
            e_amplitude_sigma, e_amplitude_pi = self.get_result_e_amplitudes()


            #
            # we propagate the emission at a long distance (fat field) to sample angles
            # Then it was propagated back to the source plane to sample the size.
            #
            if undulator.code_undul_phot == 'internal': # wofry1D
                distance = 100.
                magnification = s_phot * 10 / (theta[-1] * distance)
                mean_photon_energy = numpy.array(sampled_photon_energy).mean()  # todo: use the weighted mean?

                if True:
                    THETA = numpy.concatenate((-theta[::-1], theta[1::]), axis=None)

                    radial_e_amplitude_sigma = e_amplitude_sigma[0, :, 0]
                    RADIAL_E_AMPLITUDE = numpy.concatenate((radial_e_amplitude_sigma[::-1], radial_e_amplitude_sigma[1::]), axis=None)
                    # doble the arrays for 1D propagation

                    # from srxraylib.plot.gol import plot
                    # plot(THETA, numpy.abs(RADIAL_E_AMPLITUDE)**2, title="WOFRY E=%f, Nener=%d" % (mean_photon_energy, photon_energy.size))

                    ############################### WOFRY ############################################
                    self.__result_photon_size_farfield = {
                                                        "theta": THETA,
                                                        "radial_e_amplitude": RADIAL_E_AMPLITUDE,
                                                        "mean_photon_energy": mean_photon_energy,
                                                        "distance": distance,
                                                        "magnification": magnification}

                    self._back_propagation_for_size_calculation_wofry()

                    # we sample rays following the resulting radial distribution
                    xx = self.__result_photon_size_distribution["x"]
                    yy = self.__result_photon_size_distribution["y"]

                    # plot(xx, yy)
                    # #################################################################################
                    sampler_radial = Sampler1D(yy * numpy.abs(xx), xx)
                    r, hy, hx = sampler_radial.get_n_sampled_points_and_histogram(NRAYS, bins=101)
                    angle = numpy.random.random(NRAYS) * 2 * numpy.pi

                    x_photon = r / numpy.sqrt(2.0) * numpy.sin(angle)
                    y_photon = 0.0
                    z_photon = r / numpy.sqrt(2.0) * numpy.cos(angle)
                else: # hankel
                    # print(">>>>>>>>>>>>>>>>>>>>>>>>>", e_amplitude_sigma.shape)
                    radial_e_amplitude_sigma =  e_amplitude_sigma[0, :, 0]

                    from srxraylib.plot.gol import plot
                    plot(theta, numpy.abs(radial_e_amplitude_sigma)**2, title="HANKEL E=%f, Nener=%d" % (mean_photon_energy, photon_energy.size))

                    self.__result_photon_size_farfield = {
                                                        "theta": theta,
                                                        "radial_e_amplitude": radial_e_amplitude_sigma,
                                                        "mean_photon_energy": mean_photon_energy,
                                                        "distance": distance,
                                                        "magnification": magnification}

                    self._back_propagation_for_size_calculation_hankel()

                    # we sample rays following the resulting radial distribution
                    xx = self.__result_photon_size_distribution["x"]
                    yy = self.__result_photon_size_distribution["y"]

                    plot(xx, yy)
                    # #################################################################################
                    sampler_radial = Sampler1D(yy * numpy.abs(xx), xx)
                    r, hy, hx = sampler_radial.get_n_sampled_points_and_histogram(NRAYS, bins=101)
                    angle = numpy.random.random(NRAYS) * 2 * numpy.pi

                    x_photon = r / numpy.sqrt(2.0) * numpy.sin(angle)
                    y_photon = 0.0
                    z_photon = r / numpy.sqrt(2.0) * numpy.cos(angle)
            elif undulator.code_undul_phot == 'pysru':
                self.__result_photon_size_farfield = {}
                dict1 = self.__result_radiation
                if self.get_magnetic_structure().is_monochromatic():
                    i_prop = dict1['CART_BACKPROPAGATED_radiation'][0]  # todo something
                else:
                    i_prop = dict1['CART_BACKPROPAGATED_radiation'].sum()  # todo something
                x = dict1['CART_BACKPROPAGATED_x']
                y = dict1['CART_BACKPROPAGATED_y']

                s2d = Sampler2D(i_prop, x, y)
                sampled_x, sampled_z = s2d.get_n_sampled_points(NRAYS)

                # from srxraylib.plot.gol import plot_scatter
                # plot_scatter(1e6 * sampled_x, 1e6 * sampled_z, title='>>>>>>>>>> SAMPLED (X,Z) in microns')

                x_photon = sampled_x
                y_photon = 0.0
                z_photon = sampled_z

            elif undulator.code_undul_phot == 'srw':
                self.__result_photon_size_farfield = {}
                dict1 = self.__result_radiation
                if self.get_magnetic_structure().is_monochromatic():
                    i_prop = dict1['CART_BACKPROPAGATED_radiation'][0] # todo something
                else:
                    i_prop = dict1['CART_BACKPROPAGATED_radiation'].sum() # todo something
                x = dict1['CART_BACKPROPAGATED_x']
                y = dict1['CART_BACKPROPAGATED_y']

                # from srxraylib.plot.gol import plot_image
                # plot_image(i_prop, x, y, title='backpropagated')
                # plot_image(tmp_theta,theta,phi,aspect='auto')

                s2d = Sampler2D(i_prop, x, y)
                sampled_x, sampled_z = s2d.get_n_sampled_points(NRAYS)

                # from srxraylib.plot.gol import plot_scatter
                # plot_scatter(1e6 * sampled_x, 1e6 * sampled_z, title='>>>>>>>>>> SAMPLED (X,Z) in microns')

                x_photon = sampled_x
                y_photon = 0.0
                z_photon = sampled_z

        return x_photon, y_photon, z_photon


    def __calculate_rays(self, user_unit_to_m=1.0, F_COHER=0, verbose=1):
        """
        compute the rays in SHADOW matrix (shape (npoints,18) )
        :param F_COHER: set this flag for coherent beam
        :param user_unit_to_m: default 1.0 (m)
        :return: rays, a numpy.array((npoits,18))
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
        if undulator._FLAG_EMITTANCE:
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

        if undulator._FLAG_EMITTANCE:
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

        # beam.rays[:,6] =  1.0

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
        # A_VEC(1) = 1.0D0
        # A_VEC(2) = 0.0D0
        # A_VEC(3) = 0.0D0

        DIREC = rays[:, 3:6].copy()
        A_VEC = numpy.zeros_like(DIREC)
        A_VEC[:, 0] = 1.0

        # ! C
        # ! C   Rotate A_VEC so that it will be perpendicular to DIREC and with the
        # ! C   right components on the plane.
        # ! C
        # CALL CROSS (A_VEC,DIREC,A_TEMP)
        A_TEMP = vector_cross(A_VEC, DIREC)
        # CALL CROSS (DIREC,A_TEMP,A_VEC)
        A_VEC = vector_cross(DIREC, A_TEMP)
        # CALL NORM (A_VEC,A_VEC)
        A_VEC = vector_norm(A_VEC)
        # CALL CROSS (A_VEC,DIREC,AP_VEC)
        AP_VEC = vector_cross(A_VEC, DIREC)
        # CALL NORM (AP_VEC,AP_VEC)
        AP_VEC = vector_norm(AP_VEC)

        #
        # obtain polarization for each ray (interpolation)
        #
        if undulator._NG_E == 1: # 2D interpolation
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
        #     DENOM = SQRT(1.0D0 - 2.0D0*POL_DEG + 2.0D0*POL_DEG**2)
        #     AX = POL_DEG/DENOM
        #     CALL SCALAR (A_VEC,AX,A_VEC)
        #     ! C
        #     ! C   Same procedure for AP_VEC
        #     ! C
        #     AZ = (1-POL_DEG)/DENOM
        #     CALL SCALAR  (AP_VEC,AZ,AP_VEC)

        DENOM = numpy.sqrt(1.0 - 2.0 * POL_DEG + 2.0 * POL_DEG**2)
        AX = POL_DEG/DENOM
        for i in range(3):
            A_VEC[:, i] *= AX

        AZ = (1.0 - POL_DEG) / DENOM
        for i in range(3):
            AP_VEC[:, i] *= AZ

        rays[:, 6:9] =  A_VEC
        rays[:, 15:18] = AP_VEC

        #
        # ! C
        # ! C Now the phases of A_VEC and AP_VEC.
        # ! C
        # IF (F_COHER.EQ.1) THEN
        #     PHASEX = 0.0D0
        # ELSE
        #     PHASEX = WRAN(ISTAR1) * TWOPI
        # END IF
        # PHASEZ = PHASEX + POL_ANGLE*I_CHANGE
        #
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


    def _back_propagation_for_size_calculation_wofry(self):
        """
        Calculate the radiation_flux vs theta at a "distance"
        Back propagate to -distance
        The result is the size distrubution

        :return: None; stores results in self._photon_size_distribution
        """

        from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
        from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
        from syned.beamline.beamline_element import BeamlineElement
        from syned.beamline.element_coordinates import ElementCoordinates
        from wofryimpl.propagator.propagators1D.fresnel_zoom import FresnelZoom1D
        from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

        theta          = self.__result_photon_size_farfield["theta"]
        photon_energy  = self.__result_photon_size_farfield["mean_photon_energy"]
        distance       = self.__result_photon_size_farfield["distance"]
        magnification  = self.__result_photon_size_farfield["magnification"]
        radial_e_amplitude = self.__result_photon_size_farfield["radial_e_amplitude"] # numpy.sqrt(self.__result_photon_size_farfield["radial_flux"]) + 0j

        input_wavefront = GenericWavefront1D().initialize_wavefront_from_arrays(theta * distance, radial_e_amplitude)
        input_wavefront.set_photon_energy(photon_energy)
        # input_wavefront.set_spherical_wave(radius=distance, complex_amplitude=radial_e_amplitude)
        # input_wavefront.save_h5_file("tmp2.h5","wfr")

        optical_element = WOScreen1D()
        #
        # propagating
        #
        #
        propagation_elements = PropagationElements()
        beamline_element = BeamlineElement(optical_element=optical_element,
                        coordinates=ElementCoordinates(p=0.0,q=-distance,
                        angle_radial=numpy.radians(0.000000),
                        angle_azimuthal=numpy.radians(0.000000)))
        propagation_elements.add_beamline_element(beamline_element)
        propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),propagation_elements = propagation_elements)
        propagation_parameters.set_additional_parameters('magnification_x', magnification)

        #
        propagator = PropagationManager.Instance()
        try:
            propagator.add_propagator(FresnelZoom1D())
        except:
            pass
        output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,handler_name='FRESNEL_ZOOM_1D')

        self.__result_photon_size_distribution = {"x":output_wavefront.get_abscissas(),
                                                  "y":output_wavefront.get_intensity()}

    def _back_propagation_for_size_calculation_hankel(self):
        """
        Calculate the radiation_flux vs theta at a "distance"
        Back propagate to -distance
        The result is the size distrubution

        :return: None; stores results in self._photon_size_distribution
        """

        theta          = self.__result_photon_size_farfield["theta"]
        photon_energy  = self.__result_photon_size_farfield["mean_photon_energy"]
        distance       = self.__result_photon_size_farfield["distance"]
        magnification  = self.__result_photon_size_farfield["magnification"]
        radial_e_amplitude = self.__result_photon_size_farfield["radial_e_amplitude"] # numpy.sqrt(self.__result_photon_size_farfield["radial_flux"]) + 0j

        n2 = theta.size
###############################################################################
        from pyhank import HankelTransform
        from srxraylib.plot.gol import plot, plot_image
        do_plot = 1

        r = theta * distance
        nr = r.size

        print(">>>>>>>> rminmax: ", r[0], r[-1], r.min(), r.max())
        lambda_ = codata.h * codata.c / codata.e / photon_energy #  1.5e-10  # 488e-9  # wavelength 488nm
        k0 = 2 * numpy.pi / lambda_  # Vacuum k vector

        # Set up a :class:`.HankelTransform` object, telling it the order (``0``) and
        # the radial grid.
        # def __init__(self, order: int, max_radius: float = None, n_points: int = None,
        #              radial_grid: np.ndarray = None, k_grid: np.ndarray = None):
        H = HankelTransform(order=0, radial_grid=r)

        # Set up the electric field profile at :math:`z = 0`, and resample onto the correct radial grid
        # (``transformer.r``) as required for the QDHT.
        Er = (radial_e_amplitude)     # Initial field
        ErH = H.to_transform_r(Er)  # Resampled field

        # # Now plot an image showing the intensity as a function of radius and propagation distance.

        if do_plot:
            plot(r   * 1e6, numpy.abs(Er) ** 2,
                 H.r * 1e6, numpy.abs(ErH) ** 2,
                 # xrange=[0, 1], yrange=[0, 1],
                 xtitle="r [um]", ytitle="Field intensity /arb", title="Initial electric field distribution",
                 legend=['$|E(r)|^2$', '$|E(H.r)|^2$', ],
                 marker=[None, None,], linestyle=[None, None,],
                 )

            plot(
                 r   * 1e6, numpy.unwrap(numpy.angle(Er)),
                 H.r * 1e6, numpy.unwrap(numpy.angle(ErH)),
                 # xrange=[0, 1], yrange=[0, 1],
                 xtitle="r [um]", ytitle="Phase /arb", title="Initial electric field distribution",
                 legend=['$\\phi(r)$', '$\\phi(H.r)$'],
                 marker=[None, '+'], linestyle=[None, ''],
                 )

        # Perform Hankel Transform
        # ------------------------
        # Convert from physical field to physical wavevector
        EkrH = H.qdht(ErH)

        if do_plot:
            plot(H.kr, numpy.abs(EkrH) ** 2,
                 # xrange=[0,1], yrange=[0,1],
                 xtitle=r'Radial wave-vector ($k_r$) /rad $m^{-1}$', ytitle='Field intensity /arb.',
                 title="Radial wave-vector distribution",
                 )


        # Propagate the beam - loop
        # -------------------------
        # Do the propagation in a loop over :math:`z`

        # Pre-allocate an array for field as a function of r and z
        # Erz = numpy.zeros((nr, Nz), dtype=complex)
        kz = numpy.sqrt(k0 ** 2 - H.kr ** 2)
        print(">>>>>000 kz", kz.shape)

        phi_z = kz * distance  # Propagation phase
        EkrHz = EkrH * numpy.exp(1j * phi_z)  # Apply propagation
        print("   >>>>", EkrHz.shape)
        ErHz = H.iqdht(EkrHz)  # iQDHT
        Erz = H.to_original_r(ErHz)  # Interpolate output
        Irz = numpy.abs(Erz) ** 2

###############################################################################


        self.__result_photon_size_distribution = {"x":r,
                                                  "y":numpy.abs(Erz) ** 2}



    def _sample_photon_energy_theta_and_phi(self, NRAYS):

        #
        # sample divergences
        #
        theta         = self.__result_radiation["theta"]
        phi           = self.__result_radiation["phi"]
        photon_energy = self.__result_radiation["photon_energy"]

        if self.get_magnetic_structure().is_monochromatic():
            #2D case
            tmp = self.__result_radiation["radiation"][0, :, :].copy()
            tmp /= tmp.max()
            # correct radiation for DxDz / DthetaDphi
            tmp_theta = numpy.outer(theta, numpy.ones_like(phi))
            tmp_theta /= tmp_theta.max()
            tmp_theta += 1e-6 # to avoid zeros
            tmp *= tmp_theta
            # plot_image(tmp_theta,theta,phi,aspect='auto')
            s2d = Sampler2D(tmp, theta, phi)
            sampled_theta,sampled_phi = s2d.get_n_sampled_points(NRAYS)
            sampled_photon_energy = self.get_magnetic_structure()._EMIN
        else:
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
                     use_gaussian_approximation=1,
                     emin=10000.0,
                     emax=15000.0,
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

    # # beam =  ls.get_beam()
    #
    # print(ls.info())
    # # print(ls.to_python_code())
    #
    ls.set_energy_monochromatic_at_resonance(harmonic_number=1)
    # print(ls.info())
    # print("source size: ", ls.get_result_photon_size_sigma())

    print("dict: ")
    d = ls.get_result_dictionary()
    for key in d.keys():
        print(key)

    # print("radiation: ", d["radiation"].shape)
    # print("polarization: ", d["polarization"].shape)
    # print("phi: ", d["phi"].shape)
    # print("theta: ", d["theta"].shape)
    traj = d['trajectory']
    from srxraylib.plot.gol import plot
    import scipy.constants as codata
    plot(traj[3] * codata.c, traj[1] * codata.c, xtitle="Z", ytitle="X")
    plot(traj[3] * codata.c, traj[4], xtitle="Z", ytitle="beta_x")

