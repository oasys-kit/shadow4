"""
Undulator light source.

Computes undulator radiation distributions and samples rays according to them.

The radiation divergences (far field) are computed in polar coordinates for a more efiicient sampling.
"""
import numpy


import scipy.constants as codata
from scipy.special import erf


from shadow4.sources.s4_electron_beam import S4ElectronBeam
from shadow4.beam.s4_beam import S4Beam
from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian
from shadow4.sources.s4_light_source import S4LightSource
from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian

class S4UndulatorGaussianLightSource(S4LightSource):
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
                         magnetic_structure=magnetic_structure if not magnetic_structure is None else S4UndulatorGaussian(),
                         nrays=nrays,
                         seed=seed,
                         )


    def get_beam(self):
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
        return self.get_beam_in_gaussian_approximation()

    def calculate_spectrum(self):
        """
        Calculates the spectrum (flux normalized to one).

        Returns
        -------
        tuple
            (e, f, w) numpy arrays with photon energy (e), photon flux (f) and spectral power (w).
        """
        emin, emax, npoints = self.get_magnetic_structure().get_energy_box()
        photon_energy = numpy.linspace(emin, emax, npoints)
        flux = numpy.ones_like(photon_energy) * self.get_magnetic_structure()._flux_peak
        spectral_power = flux * codata.e * 1e3
        return photon_energy, flux, spectral_power

    def get_flux(self):
        """
        Return the integrated flux (photons/s/0.1%bw) versus photon energy.

        Returns
        -------
        tuple
            (flux, photon_energy).
        """
        photon_energy, flux, spectral_power = self.calculate_spectrum()
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

    def get_scans_for_plots(self, photon_energy_array):

        u = self.get_magnetic_structure()

        s, sp = self.get_undulator_photon_beam_sizes(
            undulator_E0=photon_energy_array,
            undulator_length=u.length(),
        )

        if u.get_flag_emittance():
            sigma_x, sigdi_x, sigma_z, sigdi_z = self.get_electron_beam().get_sigmas_all()
        else:
            sigma_x, sigdi_x, sigma_z, sigdi_z = 0, 0, 0, 0

        Qa, Qs = self._get_q_a_and_q_s()

        Sx = numpy.sqrt( (s * Qs) ** 2 + sigma_x ** 2)
        Sz = numpy.sqrt( (s * Qs) ** 2 + sigma_z ** 2)
        Spx = numpy.sqrt((sp * Qa) ** 2 + sigdi_x ** 2)
        Spz = numpy.sqrt((sp * Qa) ** 2 + sigdi_z ** 2)

        return Sx, Spx, Sz, Spz

    def get_beam_in_gaussian_approximation(self):
        """
        Returns the beam with the sampled rays using for an undulator modelled with the Gaussian approximation.

        Returns
        -------
        instance of S4Beam.
        """
        E0, delta_e, npoints = self.get_magnetic_structure().get_energy_box()
        NRAYS = self.get_nrays()
        u = self.get_magnetic_structure()

        s, sp = self.get_undulator_photon_beam_sizes(
                                                    undulator_E0=E0,
                                                    undulator_length=u.length(),
                                                    )

        if u.get_flag_emittance():
            sigma_x, sigdi_x, sigma_z, sigdi_z = self.get_electron_beam().get_sigmas_all()
        else:
            sigma_x, sigdi_x, sigma_z, sigdi_z = 0, 0, 0, 0

        Qa, Qs = self._get_q_a_and_q_s()

        Sx = numpy.sqrt( (s  * Qs) ** 2 + sigma_x ** 2)
        Sz = numpy.sqrt( (s  * Qs) ** 2 + sigma_z ** 2)
        Spx = numpy.sqrt((sp * Qa) ** 2 + sigdi_x ** 2)
        Spz = numpy.sqrt((sp * Qa) ** 2 + sigdi_z ** 2)


        m2ev = codata.c * codata.h / codata.e
        lambda1 = m2ev / E0
        txt = ""
        txt += "Photon single electron emission at wavelength %f A: \n" % (lambda1 * 1e10)
        txt += "    sigma_u:      %g um \n" % (1e6 * s)
        txt += "    sigma_uprime: %g urad \n" % (1e6 * sp)
        if u.get_flag_energy_spread():
            txt += "Correction factor for energy spread (harmonic n=%d): \n" % (u._harmonic_number)
            txt += "    Qs:      %g \n" % (Qs)
            txt += "    Qa:      %g \n" % (Qa)
        else:
            txt += "No correction for electron energy spread.\n"

        txt += "Electron sizes: \n"
        txt += "    sigma_x: %g um \n" % (1e6 * sigma_x)
        txt += "    sigma_z: %g um \n" % (1e6 * sigma_z)
        txt += "    sigma_x': %g urad \n" % (1e6 * sigdi_x)
        txt += "    sigma_z': %g urad \n" % (1e6 * sigdi_z)
        txt += "Photon source sizes (convolution): \n"
        txt += "    Sigma_x: %g um \n" % (1e6 * Sx)
        txt += "    Sigma_z: %g um \n" % (1e6 * Sz)
        txt += "    Sigma_x': %g urad \n" % (1e6 * Spx)
        txt += "    Sigma_z': %g urad \n" % (1e6 * Spz)

        print(txt)

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

        beam = a.get_beam()
        if u.is_monochromatic():
            e = numpy.zeros(NRAYS) + E0
        else:
            e = numpy.random.random(NRAYS) * (delta_e) + E0

        beam.set_photon_energy_eV(e)

        return beam

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
        if undulator.get_flag_emittance():
            sigmas = syned_electron_beam.get_sigmas_all()
            txt += "        Electron sigmaX: %g [um]\n"%(1e6*sigmas[0])
            txt += "        Electron sigmaZ: %g [um]\n"%(1e6*sigmas[2])
            txt += "        Electron sigmaX': %f urad\n"%(1e6*sigmas[1])
            txt += "        Electron sigmaZ': %f urad\n"%(1e6*sigmas[3])
        txt += "Input Undulator parameters: \n"
        txt += "        period: %f m\n"%undulator.period_length()
        txt += "        number of periods: %d\n"%undulator.number_of_periods()


        if undulator.K_vertical() is not None:
            txt += "        K-value: %f\n"%undulator.K_vertical()

        txt += "-----------------------------------------------------\n"

        txt += "Lorentz factor (gamma): %f\n"%syned_electron_beam.gamma()
        txt += "Electron velocity: %.12f c units\n"%(numpy.sqrt(1.0 - 1.0 / syned_electron_beam.gamma() ** 2))
        txt += "Undulator length: %f m\n"%(undulator.period_length()*undulator.number_of_periods())

        txt += "-----------------------------------------------------\n"
        txt += "Grids: \n"
        if undulator.is_monochromatic():
            txt += "        photon energy %f eV\n"%(undulator._photon_energy)
        else:
            txt += "        photon energy from %10.3f eV to %10.3f eV\n" % (\
                undulator._photon_energy - 0.5 * undulator._delta_e, undulator._photon_energy - 0.5 * undulator._delta_e)
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


        script += "\n\n\n# light source\nfrom shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource"
        script += "\nlight_source = S4UndulatorGaussianLightSource(name='%s', electron_beam=electron_beam, magnetic_structure=source,nrays=%s,seed=%s)" % \
                                                          (self.get_name(),self.get_nrays(),self.get_seed())

        script += "\nbeam = light_source.get_beam()"

        return script

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

    def norm_energ_spr(self, verbose=1):
        """ Tanaka & Kitamura 2009 Normalized energy spread
        equation (13)"""
        u = self.get_magnetic_structure()
        if u.get_flag_energy_spread():
            e = self.get_electron_beam()
            out = 2 * numpy.pi * u._harmonic_number * u.number_of_periods() * e._energy_spread
            if verbose:
                print("n=%d, N=%f, sigma_delta=%g" % (u._harmonic_number, u.number_of_periods(), e._energy_spread) )
                print("Normalized energy spread = %g" % (out))
            return out
        else:
            return 0.0

    def _get_q_a_and_q_s(self):
        x = self.norm_energ_spr()
        return self.q_a(x), self.q_s(x, factor=0.5)

    @classmethod
    def q_a(cls, x):
        """ Tanaka & Kitamura 2009 equation (17), forzed to give
        1 in the limit close to zero"""
        if x > 1e-5:
            f_1 = -1 + numpy.exp(-2*numpy.square(x)) + numpy.sqrt(2*numpy.pi) * x * erf(numpy.sqrt(2)*x)
            value = numpy.sqrt(2*numpy.square(x)/f_1)
        elif x < 1e-5 and x >= 0:
            value = 1.0
        else:
            raise RuntimeError('ERROR: Please provide a positive energy spread')

        return value

    @classmethod
    def q_s(cls, x, factor=1.0):
        """ Tanaka & Kitamura 2009 equation (24), please noticed the correction factor
        which in our case of using Onuki&Elleaume should be factor = 0.5 """
        return 2 * numpy.power(cls.q_a(x/4), 2/3) * factor

if __name__ == "__main__":


    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam

    electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=1e-03, current=0.2)
    electron_beam.set_sigmas_all(sigma_x=0, sigma_y=0, sigma_xp=0, sigma_yp=0)

    # magnetic structure
    from shadow4.sources.undulator.s4_undulator import S4Undulator

    source = S4UndulatorGaussian(
        period_length=0.01999998,  # syned Undulator parameter
        number_of_periods=100,  # syned Undulator parameter
        photon_energy=10000.0,
        delta_e=0.0,
        ng_e=100,  # Photon energy scan number of points
        flag_emittance=0,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_energy_spread=1,
        harmonic_number=1,
        flux_peak=1e10,
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource

    light_source = S4UndulatorGaussianLightSource(name='GaussianUndulator', electron_beam=electron_beam,
                                          magnetic_structure=source, nrays=5000, seed=5676561)

    if True:
        print(light_source.info())
        beam = light_source.get_beam()

        print(beam)

        # test plot
        from srxraylib.plot.gol import plot_scatter

        rays = beam.get_rays()
        plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')
        plot_scatter(1e6 * rays[:, 3], 1e6 * rays[:, 5], title='(Xp,Zp) in microradians')


        print(light_source.get_info())

        ener = numpy.linspace(2000.0, 80000, 100)
        Sx, Spx, Sz, Spz  = light_source.get_scans_for_plots(ener)

        from srxraylib.plot.gol import plot
        plot(ener, Sx, ener, Sz, title="sizes")

        flux, ee = light_source.get_flux()
        print(flux.shape, flux)
        plot(ee, flux, title="flux")

    print(light_source.q_a(1.0))

    nes = light_source.norm_energ_spr(verbose=1)
    print("Normalized energy spread: ", nes )
    print("Qa: ", light_source.q_a(nes))
    print("Qs: ", light_source.q_s(nes, factor=0.5))


