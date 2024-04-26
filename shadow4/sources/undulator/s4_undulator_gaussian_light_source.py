"""
Undulator light source.

Computes undulator radiation distributions and samples rays according to them.

The radiation divergences (far field) are computed in polar coordinates for a more efiicient sampling.
"""
import numpy


import scipy.constants as codata
from scipy.special import erf


from shadow4.sources.s4_electron_beam import S4ElectronBeam
from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian
from shadow4.sources.s4_light_source import S4LightSource
from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian

from shadow4.tools.logger import is_verbose, is_debug

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

        if magnetic_structure._flag_autoset_flux_central_cone:
            magnetic_structure._flux_central_cone = self.get_flux_central_cone()


    def calculate_spectrum(self):
        """
        Calculates the spectrum.

        Returns
        -------
        tuple
            (e, f, w) numpy arrays with photon energy (e) in eV, photon flux (f) in (photons/s/0.1%bw) and
            spectral power (w) in W/eV.
        """

        E0, deltaE, npoints = self.get_magnetic_structure().get_energy_box()
        if deltaE == 0: deltaE = 1
        photon_energy = numpy.linspace(E0 - 0.5 * deltaE, E0 + 0.5 * deltaE, npoints)
        flux = numpy.ones_like(photon_energy) * self.get_magnetic_structure()._flux_central_cone
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

    def get_size_and_divergence_vs_photon_energy(self, Kmin=0, Kmax=5.0, max_harmonic_number=11):
        """
        Calculates the photon source size and divergences. It makes a scan in photon energy and considers
        electron emittance and electron energy spread. It uses Gaussian approximation.

        Parameters
        ----------
        Kmin : float, optional
            The minimum value of K to be used.
        Kmax : float, optional
            The maximum value of K to be used.
        max_harmonic_number : int, optional
            The maximum harmonic to be consideres in the calculation (1,3,5,... max_harmonic_number).

        Returns
        -------
        tuple
            (Energies, SizeH, SizeV, DivergenceH, DivergenceV, Labels) arrays (n,m) n: harmonic index.
            n=0 means no energy spread considered. m is the number of energy points.
        """

        u = self.get_magnetic_structure()

        if (max_harmonic_number % 2) == 0:  # even
            max_harmonic_number = (max_harmonic_number - 1)

        if u.get_flag_energy_spread:
            nmax = max_harmonic_number // 2 + 1
        else:
            nmax = 0

        nmax += 1  # add n=0 for the zero energy spread
        nmax = int(nmax) # just in case


        npoints = int(u._NG_E)
        syned_electron_beam = self.get_electron_beam()
        Gamma = syned_electron_beam.gamma()

        Energies = numpy.zeros((nmax, npoints))
        SizeH = numpy.zeros_like(Energies)
        SizeV = numpy.zeros_like(Energies)
        DivergenceH =  numpy.zeros_like(Energies)
        DivergenceV = numpy.zeros_like(Energies)
        Labels = ["Zero energy spread"]

        if u.get_flag_emittance():
            sigma_x, sigdi_x, sigma_z, sigdi_z = self.get_electron_beam().get_sigmas_all()
        else:
            sigma_x, sigdi_x, sigma_z, sigdi_z = 0, 0, 0, 0

        ii = 0
        if u.get_flag_energy_spread:
            for i in range(1, max_harmonic_number + 1, 2):
                ii += 1
                wavelength_min = (1 + 0.5 * Kmin**2) / (2 * Gamma**2) * u.period_length()
                wavelength_max = (1 + 0.5 * Kmax**2) / (2 * Gamma**2) * u.period_length()
                energy_max = codata.h * codata.c / codata.e / wavelength_min
                energy_min = codata.h * codata.c / codata.e / wavelength_max

                photon_energy = numpy.linspace(energy_min * i, energy_max * i, npoints)

                s, sp = self.get_undulator_photon_beam_sizes(
                    undulator_E0=photon_energy,
                    undulator_length=u.length(),
                )

                Qa, Qs = self._get_q_a_and_q_s(harmonic_number=i)

                Energies[ii, :] = photon_energy
                SizeH[ii, :] = numpy.sqrt( (s * Qs) ** 2 + sigma_x ** 2)
                SizeV[ii, :] = numpy.sqrt( (s * Qs) ** 2 + sigma_z ** 2)
                DivergenceH[ii, :] = numpy.sqrt((sp * Qa) ** 2 + sigdi_x ** 2)
                DivergenceV[ii, :] = numpy.sqrt((sp * Qa) ** 2 + sigdi_z ** 2)
                Labels.append("harmonic n=%d" % i)

        # add the zero spread in index 0
        if u.get_flag_energy_spread:
            emin = Energies[1:, :].min()
            emax = Energies[1:, :].max()
        else:
            wavelength_min = (1 + 0.5 * Kmin ** 2) / (2 * Gamma ** 2) * u.period_length()
            wavelength_max = (1 + 0.5 * Kmax ** 2) / (2 * Gamma ** 2) * u.period_length()
            emax = codata.h * codata.c / codata.e / wavelength_min
            emin = codata.h * codata.c / codata.e / wavelength_max

        photon_energy = numpy.linspace(emin, emax, npoints)

        s, sp = self.get_undulator_photon_beam_sizes(
            undulator_E0=photon_energy,
            undulator_length=u.length(),
        )

        Qa, Qs = self._get_q_a_and_q_s(harmonic_number=0)

        Energies[0, :] = photon_energy
        SizeH[0, :] = numpy.sqrt((s * Qs) ** 2 + sigma_x ** 2)
        SizeV[0, :] = numpy.sqrt((s * Qs) ** 2 + sigma_z ** 2)
        DivergenceH[0, :] = numpy.sqrt((sp * Qa) ** 2 + sigdi_x ** 2)
        DivergenceV[0, :] = numpy.sqrt((sp * Qa) ** 2 + sigdi_z ** 2)

        return Energies, SizeH, SizeV, DivergenceH, DivergenceV, Labels

    def get_beam(self):
        """
        Creates the beam as emitted by the undulator.

        Returns
        -------
        instance of S4Beam
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

        if is_verbose(): print(txt)

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
            e = (numpy.random.random(NRAYS) - 0.5) * delta_e + E0

        beam.set_photon_energy_eV(e)

        return beam

    def get_info(self, debug=False):
        """
        Returns the specific information for the Gaussian undulator light source.

        Returns
        -------
        str
        """
        syned_electron_beam = self.get_electron_beam()
        undulator = self.get_magnetic_structure()


        txt = ""
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

        References
        ----------
        Equations 25 & 30 Generalities on the synchrotron radiation, Pascal Elleaume
        in H. Onuki. P. Elleaume: “Undulators, Wigglers and their Applications Taylor & Francis, New York 2003.

        """
        m2ev = codata.c * codata.h / codata.e
        lambda1 = m2ev / undulator_E0
        # calculate sizes of the photon undulator beam
        # see formulas 25 & 30 in Elleaume (Onaki & Elleaume)
        s_phot = 2.740 / (4e0 * numpy.pi) * numpy.sqrt(undulator_length * lambda1)
        sp_phot = 0.69 * numpy.sqrt(lambda1 / undulator_length)
        return s_phot, sp_phot

    def norm_energ_spr(self, harmonic_number=None):
        """
        Calculate the "normalized" electron energy spread: 2 pi * harmonic_number * n_periods * energy spread DE/E.

        Parameters
        ----------
        harmonic_number : int or None, optional
            the harmonic number to be used. If None, use the one define in S4GaussianUndulator, otherwise use the
            one defined here. It must be odd.

        Returns
        -------
        float
            The normalized energy spread.

        References
        ----------
        Equation 13 in
        Tanaka, T. & Kitamura, H. (2009). Journal of Synchrotron Radiation, 16(3), 380–386

        """
        # Tanaka & Kitamura 2009 Normalized energy spread equation (13)
        u = self.get_magnetic_structure()
        if harmonic_number is None: harmonic_number = u._harmonic_number

        if u.get_flag_energy_spread():
            e = self.get_electron_beam()
            out = 2 * numpy.pi * harmonic_number * u.number_of_periods() * e._energy_spread
            if is_verbose():
                print("n=%d, N=%f, sigma_delta=%g; Normalized energy spread = %g" %
                      (harmonic_number, u.number_of_periods(), e._energy_spread, out) )
            return out
        else:
            return 0.0

    def _get_q_a_and_q_s(self, harmonic_number=None):
        x = self.norm_energ_spr(harmonic_number=harmonic_number)
        return self.q_a(x), self.q_s(x, factor=0.5)

    def get_flux_central_cone(self, K=None):
        """
        Calculate the flux in the central code.

        Parameters
        ----------
        K : float or None, optional
            If None, K is calculated for the photon energy and harmonic_number defined in S4UndulatorGaussian.
            Otherwise, use the value defined here.

        Returns
        -------
        float
            The flux in photons/s/0.1%bw.

        References
        ----------
        see pag 2.8 in X-RAY DATA BOOKLET https://xdb.lbl.gov/

        """
        # see X-ray data booklet pag 2.8
        u = self.get_magnetic_structure()
        syned_electron_beam = self.get_electron_beam()
        current = syned_electron_beam.current()
        if u._flag_energy_spread:
            wavelength = codata.c * codata.h / codata.e / (u._photon_energy / u._harmonic_number)
        else:
            wavelength = codata.c * codata.h / codata.e / u._photon_energy
        if K is None:
            K = numpy.sqrt( 2 * (2 * wavelength * syned_electron_beam.gamma()**2 / u.period_length() - 1))

            if numpy.isnan(K):
                print("\n\n*** Warning: impossible to get K ***\n")
                K = 0.0


        out = (numpy.pi * codata.alpha * 1e-3 / codata.e * u.number_of_periods() * current * self.Qn(K, u._harmonic_number) )

        if is_verbose():
            print("*** Flux Calculation ***")
            print("  Using target E=%f eV; n=%f periods, current=%f => K=%f" % (u._photon_energy, u._harmonic_number, current, K))
            print("  Calculated Flux: %g photons/s/0.1%%bw" % (out))
            # out2 = 1.431e14 * u.number_of_periods() * current * self.Qn(K, u._harmonic_number)
            # print("  Another way to get Flux: %g photons/s" % out2)

        return out

    @classmethod
    def q_a(cls, x):
        """
        Universal function Qa by Tanaka & Kitamura 2009 equation (17).

        Parameters
        ----------
        x : float
            The argument.

        Returns
        -------
        float
            Qa(x).

        References
        ----------
        Equation 17 in
        Tanaka, T. & Kitamura, H. (2009). Journal of Synchrotron Radiation, 16(3), 380–386
        """

        # Tanaka & Kitamura 2009 equation (17), forzed to give 1 in the limit close to zero.
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
        """
        Universal function Qs by Tanaka & Kitamura 2009 equation (24).

        Parameters
        ----------
        x : float
            The argument.
        factor : float, optional
            A multiplicative factor (set to 0.5 to force Qs(0)=1 ).

        Returns
        -------
        float
            Qs(x).

        References
        ----------
        Equation 24 in
        Tanaka, T. & Kitamura, H. (2009). Journal of Synchrotron Radiation, 16(3), 380–386
        """

        # Tanaka & Kitamura 2009 equation (24), please noticed the correction factor
        # which in our case of using Onuki&Elleaume should be factor = 0.5
        return 2 * numpy.power(cls.q_a(x/4), 2/3) * factor

    @classmethod
    def Fn(cls, K, n):
        """
        Calculate the universal function Fn(x).

        Parameters
        ----------
        K : float
            The K value
        n : int, optional
            The harmonic number (it must be odd).

        Returns
        -------
        float
            Fn(x)

        References
        ----------
        see fig 2.4 in X-RAY DATA BOOKLET https://xdb.lbl.gov/

        """
        from scipy.special import jn, yn, jv, yv
        if n%2 == 1 :
            cst1 = ((n * K)/(1. + (K**2) / 2.))**2
            cst2 = (n * K**2) / (4. + (2. * K**2))
            Fn = cst1*( jn(0.5 * (n - 1), cst2) - jn(0.5 * (n + 1), cst2))**2
        else :
            Fn = 0.0
        return Fn

    @classmethod
    def Qn(cls, K, n):
        """
        Calculate the universal function Qn(x).

        Parameters
        ----------
        K : float
            The K value
        n : int, optional
            The harmonic number (it must be odd).

        Returns
        -------
        float
            Qn(x)

        References
        ----------
        see fig 2.6 in X-RAY DATA BOOKLET https://xdb.lbl.gov/

        """

        if n == 0 :
            raise Exception(' the harmonic number can not be 0')
        res=(1. + 0.5 * K**2) * cls.Fn(K, n) / n
        return res

    @classmethod
    def theoretical_flux_integrated_central_cone(cls, N=100.0, n=1., current=0.2, K=1.68):
        """
        Calculate the theoretical flux in the central cone.

        Parameters
        ----------
        N : int, optional
            The number of periods of the undulator.
        n : int, optional
            The harmonic number (it must be odd).
        current: float, optional
            The storage ring electron current in A.
        K : float, optional
            The K value

        Returns
        -------
        float
            The flux in photons/s/0.1%bw.

        References
        ----------
        see pag. 2.8 in X-RAY DATA BOOKLET https://xdb.lbl.gov/

        """

        # see X-ray data booklet pag 2.8
        res=1.431e14 * N * current * cls.Qn(K, n)
        return res


if __name__ == "__main__":

    electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=1e-03, current=0.2)
    electron_beam.set_sigmas_all(sigma_x=3.01836e-05, sigma_y=3.63641e-06, sigma_xp=4.36821e-06, sigma_yp=1.37498e-06)

    # magnetic structure
    from shadow4.sources.undulator.s4_undulator import S4Undulator

    source = S4UndulatorGaussian(
        period_length=0.018,  # syned Undulator parameter
        number_of_periods=111,  # syned Undulator parameter
        photon_energy=10000.0,
        delta_e=0.0,
        ng_e=100,  # Photon energy scan number of points
        flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_energy_spread=1,
        harmonic_number=1,
        flag_autoset_flux_central_cone=1,
        flux_central_cone=1e10,
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource

    light_source = S4UndulatorGaussianLightSource(name='GaussianUndulator', electron_beam=electron_beam,
                                          magnetic_structure=source, nrays=5000, seed=5676561)

    # print(light_source.info())
    #
    # if True:
    #     beam = light_source.get_beam()
    #
    #     print(beam)
    #
    #     # test plot
    #     from srxraylib.plot.gol import plot_scatter
    #
    #     rays = beam.get_rays()
    #     plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')
    #     plot_scatter(1e6 * rays[:, 3], 1e6 * rays[:, 5], title='(Xp,Zp) in microradians')
    #
    # if True:
    #     from srxraylib.plot.gol import plot
    #
    #
    #     flux, ee = light_source.get_flux()
    #     plot(ee, flux, title="flux")
    #
    #     Energies, SizeH, SizeV, DivergenceH, DivergenceV, Labels =\
    #         light_source.get_size_and_divergence_vs_photon_energy(0.2, 1.479, max_harmonic_number=11)
    #
    #     plot(Energies[0, :], SizeH[0, :],
    #          Energies[1, :], SizeH[1, :],
    #          Energies[2, :], SizeH[2, :],
    #          Energies[3, :], SizeH[3, :],
    #          Energies[4, :], SizeH[4, :],
    #          Energies[5, :], SizeH[5, :],
    #          Energies[6, :], SizeH[6, :],
    #          legend=Labels, title="Size H", show=0)
    #
    #     plot(Energies[0, :], SizeV[0, :],
    #          Energies[1, :], SizeV[1, :],
    #          Energies[2, :], SizeV[2, :],
    #          Energies[3, :], SizeV[3, :],
    #          Energies[4, :], SizeV[4, :],
    #          Energies[5, :], SizeV[5, :],
    #          Energies[6, :], SizeV[6, :],
    #          legend=Labels, title="Size V", show=0)
    #
    #     plot(Energies[0, :], DivergenceH[0, :],
    #          Energies[1, :], DivergenceH[1, :],
    #          Energies[2, :], DivergenceH[2, :],
    #          Energies[3, :], DivergenceH[3, :],
    #          Energies[4, :], DivergenceH[4, :],
    #          Energies[5, :], DivergenceH[5, :],
    #          Energies[6, :], DivergenceH[6, :],
    #          legend=Labels, title="Divergence H", show=0)
    #
    #     plot(Energies[0, :], DivergenceV[0, :],
    #          Energies[1, :], DivergenceV[1, :],
    #          Energies[2, :], DivergenceV[2, :],
    #          Energies[3, :], DivergenceV[3, :],
    #          Energies[4, :], DivergenceV[4, :],
    #          Energies[5, :], DivergenceV[5, :],
    #          Energies[6, :], DivergenceV[6, :],
    #          legend=Labels, title="Divergence V")
    #
    #
    # if True:
    #     print(light_source.q_a(1.0))
    #     nes = light_source.norm_energ_spr()
    #     print("Normalized energy spread: ", nes )
    #     print("Qa: ", light_source.q_a(nes))
    #     print("Qs: ", light_source.q_s(nes, factor=0.5))
    #
    #     print("Flux: %g photons/s/0.1%%bw" % light_source.get_magnetic_structure()._flux_central_cone)

    print(light_source.get_info())