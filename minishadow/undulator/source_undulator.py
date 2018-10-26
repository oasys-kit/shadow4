__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "30-08-2018"

"""

Undulator code: computes undulator radiation distributions and samples rays according to them.

Fully replaces and upgrades the shadow3 undulator model.

The radiation is calculating using one of three cods: internal, pySRU and SRW. The internal
code is indeed based on pySRU.

The radiation divergences (far field) are computed in polar coordinates for a more efiicient sampling.

Usage:

su = SourceUndulator()        # use keywords to define parameters. It uses syned for electron beam and vertical undulator
rays = su.calculate_rays()    # sample rays. Result is a numpy.array of shape (NRAYS,18), exactly the same as in shadow3


"""


import numpy

from srxraylib.util.inverse_method_sampler import Sampler1D, Sampler2D, Sampler3D
import scipy.constants as codata
from scipy import interpolate

from syned.storage_ring.magnetic_structures.undulator import Undulator
from syned.storage_ring.electron_beam import ElectronBeam

from orangecontrib.shadow.util.undulator.source_undulator_factory import SourceUndulatorFactory
from orangecontrib.shadow.util.undulator.source_undulator_factory_srw import SourceUndulatorFactorySrw
from orangecontrib.shadow.util.undulator.source_undulator_factory_pysru import SourceUndulatorFactoryPysru


INTEGRATION_METHOD = 1 # 0=sum, 1=trapz

class SourceUndulator(object):
    def __init__(self,name="",
                 syned_electron_beam=None,
                 syned_undulator=None,
                 emin=10000.0,               # Photon energy scan from energy (in eV)
                 emax=11000.0,               # Photon energy scan to energy (in eV)
                 ng_e=11,                    # Photon energy scan number of points
                 maxangle=50e-6,             # Maximum radiation semiaperture in RADIANS
                 ng_t=31,                    # Number of points in angle theta
                 ng_p=21,                    # Number of points in angle phi
                 ng_j=20,                    # Number of points in electron trajectory (per period) for internal calculation only
                 code_undul_phot="internal", # internal, pysru, srw
                 flag_emittance=0,           # when sampling rays: Use emittance (0=No, 1=Yes)
                 flag_size=0,                # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
                 ):

        # # Machine
        if syned_electron_beam is None:
            self.syned_electron_beam = ElectronBeam()
        else:
            self.syned_electron_beam = syned_electron_beam

        # # Undulator
        if syned_undulator is None:
            self.syned_undulator = Undulator()
        else:
            self.syned_undulator = syned_undulator

        # Photon energy scan
        self._EMIN            = emin   # Photon energy scan from energy (in eV)
        self._EMAX            = emax   # Photon energy scan to energy (in eV)
        self._NG_E            = ng_e   # Photon energy scan number of points
        # Geometry
        self._MAXANGLE        = maxangle   # Maximum radiation semiaperture in RADIANS
        self._NG_T            = ng_t       # Number of points in angle theta
        self._NG_P            = ng_p       # Number of points in angle phi
        self._NG_J            = ng_j       # Number of points in electron trajectory (per period)
        # ray tracing
        # self.SEED            = SEED   # Random seed
        # self.NRAYS           = NRAYS  # Number of rays

        self.code_undul_phot = code_undul_phot

        self._FLAG_EMITTANCE  =  flag_emittance # Yes  # Use emittance (0=No, 1=Yes)
        self._FLAG_SIZE  =  flag_size # 0=point,1=Gaussian,2=backpropagate Divergences

        # results of calculations

        self._result_radiation = None
        self._result_photon_size_distribution = None
        self._result_photon_size_sigma = None


    def info(self,debug=False):
        """
        gets text info

        :param debug: if True, list the undulator variables (Default: debug=True)
        :return:
        """
        # list all non-empty keywords
        txt = ""


        txt += "-----------------------------------------------------\n"

        txt += "Input Electron parameters: \n"
        txt += "        Electron energy: %f geV\n"%self.syned_electron_beam._energy_in_GeV
        txt += "        Electron current: %f A\n"%self.syned_electron_beam._current
        if self._FLAG_EMITTANCE:
            sigmas = self.syned_electron_beam.get_sigmas_all()
            txt += "        Electron sigmaX: %g [um]\n"%(1e6*sigmas[0])
            txt += "        Electron sigmaZ: %g [um]\n"%(1e6*sigmas[2])
            txt += "        Electron sigmaX': %f urad\n"%(1e6*sigmas[1])
            txt += "        Electron sigmaZ': %f urad\n"%(1e6*sigmas[3])
        txt += "Input Undulator parameters: \n"
        txt += "        period: %f m\n"%self.syned_undulator.period_length()
        txt += "        number of periods: %d\n"%self.syned_undulator.number_of_periods()
        txt += "        K-value: %f\n"%self.syned_undulator.K_vertical()

        txt += "-----------------------------------------------------\n"

        txt += "Lorentz factor (gamma): %f\n"%self.syned_electron_beam.gamma()
        txt += "Electron velocity: %.12f c units\n"%(numpy.sqrt(1.0 - 1.0 / self.syned_electron_beam.gamma() ** 2))
        txt += "Undulator length: %f m\n"%(self.syned_undulator.period_length()*self.syned_undulator.number_of_periods())
        K_to_B = (2.0 * numpy.pi / self.syned_undulator.period_length()) * codata.m_e * codata.c / codata.e

        txt += "Undulator peak magnetic field: %f T\n"%(K_to_B*self.syned_undulator.K_vertical())
        txt += "Resonances: \n"
        txt += "        harmonic number [n]                   %10d %10d %10d \n"%(1,3,5)
        txt += "        wavelength [A]:                       %10.6f %10.6f %10.6f   \n"%(\
                                                                1e10*self.syned_undulator.resonance_wavelength(self.syned_electron_beam.gamma(),harmonic=1),
                                                                1e10*self.syned_undulator.resonance_wavelength(self.syned_electron_beam.gamma(),harmonic=3),
                                                                1e10*self.syned_undulator.resonance_wavelength(self.syned_electron_beam.gamma(),harmonic=5))
        txt += "        energy [eV]   :                       %10.3f %10.3f %10.3f   \n"%(\
                                                                self.syned_undulator.resonance_energy(self.syned_electron_beam.gamma(),harmonic=1),
                                                                self.syned_undulator.resonance_energy(self.syned_electron_beam.gamma(),harmonic=3),
                                                                self.syned_undulator.resonance_energy(self.syned_electron_beam.gamma(),harmonic=5))
        txt += "        frequency [Hz]:                       %10.3g %10.3g %10.3g   \n"%(\
                                                                1e10*self.syned_undulator.resonance_frequency(self.syned_electron_beam.gamma(),harmonic=1),
                                                                1e10*self.syned_undulator.resonance_frequency(self.syned_electron_beam.gamma(),harmonic=3),
                                                                1e10*self.syned_undulator.resonance_frequency(self.syned_electron_beam.gamma(),harmonic=5))
        txt += "        central cone 'half' width [urad]:     %10.6f %10.6f %10.6f   \n"%(\
                                                                1e6*self.syned_undulator.gaussian_central_cone_aperture(self.syned_electron_beam.gamma(),1),
                                                                1e6*self.syned_undulator.gaussian_central_cone_aperture(self.syned_electron_beam.gamma(),3),
                                                                1e6*self.syned_undulator.gaussian_central_cone_aperture(self.syned_electron_beam.gamma(),5))
        txt += "        first ring at [urad]:                 %10.6f %10.6f %10.6f   \n"%(\
                                                                1e6*self.get_resonance_ring(1,1),
                                                                1e6*self.get_resonance_ring(3,1),
                                                                1e6*self.get_resonance_ring(5,1))

        txt += "-----------------------------------------------------\n"
        txt += "Grids: \n"
        if self._NG_E == 1:
            txt += "        photon energy %f eV\n"%(self._EMIN)
        else:
            txt += "        photon energy from %10.3f eV to %10.3f eV\n"%(self._EMIN,self._EMAX)
        txt += "        number of points for the trajectory: %d\n"%(self._NG_J)
        txt += "        number of energy points: %d\n"%(self._NG_E)
        txt += "        maximum elevation angle: %f urad\n"%(1e6*self._MAXANGLE)
        txt += "        number of angular elevation points: %d\n"%(self._NG_T)
        txt += "        number of angular azimuthal points: %d\n"%(self._NG_P)
        # txt += "        number of rays: %d\n"%(self.NRAYS)
        # txt += "        random seed: %d\n"%(self.SEED)
        txt += "-----------------------------------------------------\n"

        txt += "calculation code: %s\n"%self.code_undul_phot
        if self._result_radiation is None:
            txt += "radiation: NOT YET CALCULATED\n"
        else:
            txt += "radiation: CALCULATED\n"
        txt += "Sampling: \n"
        if self._FLAG_SIZE == 0:
            flag = "point"
        elif self._FLAG_SIZE == 1:
            flag = "Gaussian"
        elif self._FLAG_SIZE == 2:
            flag = "Far field backpropagated"

        txt += "        Photon source size sampling flag: %d (%s)\n"%(self._FLAG_SIZE,flag)
        if self._FLAG_SIZE == 1:
            if self._result_photon_size_sigma is not None:
                txt += "        Photon source size sigma (Gaussian): %6.3f um \n"%(1e6*self._result_photon_size_sigma)

        txt += "-----------------------------------------------------\n"
        return txt

    def get_resonance_ring(self,harmonic_number=1, ring_order=1):
        return 1.0/self.syned_electron_beam.gamma()*numpy.sqrt( ring_order / harmonic_number * (1+0.5*self.syned_undulator.K_vertical()**2) )


    def set_energy_monochromatic_at_resonance(self,harmonic_number):

        self.set_energy_monochromatic(self.syned_undulator.resonance_energy(
            self.syned_electron_beam.gamma(),harmonic=harmonic_number))
        # take 3*sigma - _MAXANGLE is in RAD
        self._MAXANGLE = 3 * 0.69 * self.syned_undulator.gaussian_central_cone_aperture(self.syned_electron_beam.gamma(),harmonic_number)

    def set_energy_monochromatic(self,emin):
        """
        Sets a single energy line for the source (monochromatic)
        :param emin: the energy in eV
        :return:
        """
        self._EMIN = emin
        self._EMAX = emin
        self._NG_E = 1


    def set_energy_box(self,emin,emax,npoints=None):
        """
        Sets a box for photon energy distribution for the source
        :param emin:  Photon energy scan from energy (in eV)
        :param emax:  Photon energy scan to energy (in eV)
        :param npoints:  Photon energy scan number of points (optinal, if not set no changes)
        :return:
        """

        self._EMIN = emin
        self._EMAX = emax
        if npoints != None:
            self._NG_E = npoints

    def get_energy_box(self):
        """
        Gets the limits of photon energy distribution for the source
        :return: emin,emax,number_of_points
        """
        return self._EMIN,self._EMAX,self._NG_E


    def calculate_radiation(self):

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


        self._result_radiation = None

        # undul_phot
        if self.code_undul_phot == 'internal':
            undul_phot_dict = SourceUndulatorFactory.undul_phot(E_ENERGY  = self.syned_electron_beam.energy(),
                                         INTENSITY = self.syned_electron_beam.current(),
                                         LAMBDAU   = self.syned_undulator.period_length(),
                                         NPERIODS  = self.syned_undulator.number_of_periods(),
                                         K         = self.syned_undulator.K(),
                                         EMIN      = self._EMIN,
                                         EMAX      = self._EMAX,
                                         NG_E      = self._NG_E,
                                         MAXANGLE  = self._MAXANGLE,
                                         NG_T      = self._NG_T,
                                         NG_P      = self._NG_P,
                                         number_of_trajectory_points = self._NG_J)

        elif self.code_undul_phot == 'pysru' or  self.code_undul_phot == 'pySRU':
            undul_phot_dict = SourceUndulatorFactoryPysru.undul_phot(E_ENERGY  = self.syned_electron_beam.energy(),
                                         INTENSITY = self.syned_electron_beam.current(),
                                         LAMBDAU   = self.syned_undulator.period_length(),
                                         NPERIODS  = self.syned_undulator.number_of_periods(),
                                         K         = self.syned_undulator.K(),
                                         EMIN      = self._EMIN,
                                         EMAX      = self._EMAX,
                                         NG_E      = self._NG_E,
                                         MAXANGLE  = self._MAXANGLE,
                                         NG_T      = self._NG_T,
                                         NG_P      = self._NG_P,)
        elif self.code_undul_phot == 'srw' or  self.code_undul_phot == 'SRW':
            undul_phot_dict = SourceUndulatorFactorySrw.undul_phot(E_ENERGY  = self.syned_electron_beam.energy(),
                                         INTENSITY = self.syned_electron_beam.current(),
                                         LAMBDAU   = self.syned_undulator.period_length(),
                                         NPERIODS  = self.syned_undulator.number_of_periods(),
                                         K         = self.syned_undulator.K(),
                                         EMIN      = self._EMIN,
                                         EMAX      = self._EMAX,
                                         NG_E      = self._NG_E,
                                         MAXANGLE  = self._MAXANGLE,
                                         NG_T      = self._NG_T,
                                         NG_P      = self._NG_P,)
        else:
            raise Exception("Not implemented undul_phot code: "+self.code_undul_phot)

        # add some info
        undul_phot_dict["code_undul_phot"] = self.code_undul_phot
        undul_phot_dict["info"] = self.info()

        self._result_radiation = undul_phot_dict

    #
    # get from results
    #

    def get_result_dictionary(self):
        if self._result_radiation is None:
            self.calculate_radiation()
        return self._result_radiation

    def get_result_radiation(self):
        return self.get_result_dictionary()["radiation"]

    def get_result_polarization(self):
        return self.get_result_dictionary()["polarization"]

    def get_result_polarisation(self):
        return self.get_result_polarization()

    def get_result_theta(self):
        return self.get_result_dictionary()["theta"]

    def get_result_phi(self):
        return self.get_result_dictionary()["phi"]

    def get_result_photon_energy(self):
        return self.get_result_dictionary()["photon_energy"]


    def get_radiation_polar(self):
        return self.get_result_radiation(),self.get_result_photon_energy(),self.get_result_theta(),self.get_result_phi()

    def get_radiation(self):
        return self.get_radiation_polar()

    def get_radiation_interpolated_cartesian(self,npointsx=100,npointsz=100,thetamax=None):

        radiation,photon_energy, thetabm,phi = self.get_radiation_polar()

        if thetamax is None:
            thetamax = thetabm.max()

        vx = numpy.linspace(-1.1*thetamax,1.1*thetamax,npointsx)
        vz = numpy.linspace(-1.1*thetamax,1.1*thetamax,npointsz)
        VX = numpy.outer(vx,numpy.ones_like(vz))
        VZ = numpy.outer(numpy.ones_like(vx),vz)
        VY = numpy.sqrt(1 - VX**2 - VZ**2)

        THETA = numpy.abs(numpy.arctan( numpy.sqrt(VX**2+VZ**2)/VY))
        PHI = numpy.arctan2( numpy.abs(VZ),numpy.abs(VX))

        radiation_interpolated = numpy.zeros((radiation.shape[0],npointsx,npointsz))

        for i in range(radiation.shape[0]):
            interpolator_value = interpolate.RectBivariateSpline(thetabm, phi, radiation[i])
            radiation_interpolated[i] = interpolator_value.ev(THETA, PHI)

        return radiation_interpolated,photon_energy,vx,vz

    def get_power_density(self):

        radiation = self.get_result_radiation().copy()
        theta = self.get_result_theta()
        phi = self.get_result_phi()
        photon_energy = self.get_result_photon_energy()

        step_e = photon_energy[1]-photon_energy[0]

        for i in range(radiation.shape[0]):
            radiation[i] *= 1e-3 * photon_energy[i] # photons/eV/rad2 -> photons/0.1%bw/rad2

        if INTEGRATION_METHOD == 0:
            power_density = radiation.sum(axis=0) * step_e * codata.e * 1e3 # W/rad2
        else:
            power_density = numpy.trapz(radiation,photon_energy,axis=0) * codata.e * 1e3 # W/rad2

        return power_density,theta,phi


    def get_power_density_interpolated_cartesian(self,npointsx=100,npointsz=100,thetamax=None):

        power_density_polar,theta,phi = self.get_power_density()

        if thetamax is None:
            thetamax = theta.max()

        vx = numpy.linspace(-1.1*thetamax,1.1*thetamax,npointsx)
        vz = numpy.linspace(-1.1*thetamax,1.1*thetamax,npointsz)
        VX = numpy.outer(vx,numpy.ones_like(vz))
        VZ = numpy.outer(numpy.ones_like(vx),vz)
        VY = numpy.sqrt(1 - VX**2 - VZ**2)

        THETA = numpy.abs(numpy.arctan( numpy.sqrt(VX**2+VZ**2)/VY))
        PHI = numpy.arctan2( numpy.abs(VZ),numpy.abs(VX))

        interpolator_value = interpolate.RectBivariateSpline(theta, phi, power_density_polar)
        power_density_cartesian = interpolator_value.ev(THETA, PHI)

        return power_density_cartesian,vx,vz

    def get_flux_and_spectral_power(self):

        radiation2 = self.get_result_radiation().copy()
        theta = self.get_result_theta()
        phi = self.get_result_phi()
        photon_energy = self.get_result_photon_energy()
        THETA = numpy.outer(theta,numpy.ones_like(phi))
        for i in range(radiation2.shape[0]):
            radiation2[i] *= THETA

        if INTEGRATION_METHOD == 0:
            flux = radiation2.sum(axis=2).sum(axis=1) * (1e-3*photon_energy) # photons/eV -> photons/0.1%bw
            flux *= 4 * (theta[1]-theta[0]) * (phi[1]-phi[0]) # adding the four quadrants!
        else:
            flux = 4 * numpy.trapz(numpy.trapz(radiation2,phi,axis=2),theta,axis=1) * (1e-3*photon_energy) # photons/eV -> photons/0.1%bw


        spectral_power = flux*codata.e*1e3

        return flux,spectral_power,photon_energy

    def get_flux(self):
        flux,spectral_power,photon_energy = self.get_flux_and_spectral_power()
        return flux,photon_energy

    def get_spectral_power(self):
        flux,spectral_power,photon_energy = self.get_flux_and_spectral_power()
        return spectral_power,photon_energy

    def get_photon_size_distribution(self):
        return self._result_photon_size_distribution["x"],self._result_photon_size_distribution["y"]

    def calculate_rays(self,user_unit_to_m=1.0,F_COHER=0,NRAYS=5000,SEED=36255655452):
        """
        compute the rays in SHADOW matrix (shape (npoints,18) )
        :param F_COHER: set this flag for coherent beam
        :param user_unit_to_m: default 1.0 (m)
        :return: rays, a numpy.array((npoits,18))
        """

        if self._result_radiation is None:
            self.calculate_radiation()

        sampled_photon_energy,sampled_theta,sampled_phi = self._sample_photon_energy_theta_and_phi(NRAYS)

        if SEED != 0:
            numpy.random.seed(SEED)


        sigmas = self.syned_electron_beam.get_sigmas_all()

        rays = numpy.zeros((NRAYS,18))

        #
        # sample sizes (cols 1-3)
        #


        if self._FLAG_EMITTANCE:
            x_electron = numpy.random.normal(loc=0.0,scale=sigmas[0],size=NRAYS)
            y_electron = 0.0
            z_electron = numpy.random.normal(loc=0.0,scale=sigmas[2],size=NRAYS)
        else:
            x_electron = 0.0
            y_electron = 0.0
            z_electron = 0.0




        # calculate (and stores) sizes of the photon undulator beam
        # see formulas 25 & 30 in Elleaume (Onaki & Elleaume)
        # sp_phot = 0.69*numpy.sqrt(lambda1/undulator_length)
        undulator_length = self.syned_undulator.length()
        lambda1 = codata.h*codata.c/codata.e / numpy.array(sampled_photon_energy).mean()
        s_phot = 2.740/(4e0*numpy.pi)*numpy.sqrt(undulator_length*lambda1)
        self._result_photon_size_sigma = s_phot

        if self._FLAG_SIZE == 0:
            x_photon = 0.0
            y_photon = 0.0
            z_photon = 0.0
            # for plot, a delta
            x = numpy.linspace(-1e-6,1e-6,101)
            y = numpy.zeros_like(x)
            y[y.size//2] = 1.0
            self._result_photon_size_distribution = {"x":x,"y":y}
        elif self._FLAG_SIZE == 1:
            # TODO: I added this correction to obtain the sigma in the RADIAL coordinate, not in x and z.
            # RODO: TO be verified!
            s_phot_corrected = s_phot / numpy.sqrt(2)

            cov = [[s_phot_corrected**2, 0], [0, s_phot_corrected**2]]
            mean = [0.0,0.0]

            tmp = numpy.random.multivariate_normal(mean, cov, NRAYS)
            x_photon = tmp[:,0]
            y_photon = 0.0
            z_photon = tmp[:,1]

            # for plot, a Gaussian
            x = numpy.linspace(-5*s_phot,5*s_phot,101)
            y = numpy.exp(-x**2/2/s_phot**2)
            self._result_photon_size_distribution = {"x":x,"y":y}


        elif self._FLAG_SIZE == 2:
            # we need to retrieve the emission as a function of the angle
            radiation,photon_energy, theta,phi = self.get_radiation_polar()

            mean_photon_energy = numpy.array(sampled_photon_energy).mean() # todo: use the weighted mean?
            shape_radiation = radiation.shape
            radial_flux = radiation.sum(axis=2) / shape_radiation[2]
            radial_flux = radial_flux.sum(axis=0) / shape_radiation[0]
            # doble the arrays for 1D propagation
            THETA = numpy.concatenate((-theta[::-1],theta[1::]),axis=None)
            RADIAL_FLUX = numpy.concatenate( (radial_flux[::-1],radial_flux[1::]),axis=None)

            #
            # we propagate the emission at a long distance back to the source plane
            #
            distance = 100.

            magnification = s_phot*10 / (theta[-1]*distance)

            # do the propagation; result is stored in self._photon_size_distribution
            self._back_propagation_for_size_calculation(THETA,RADIAL_FLUX,
                                mean_photon_energy,distance=distance,magnification=magnification)


            # we sample rays following the resulting radial distribution
            xx = self._result_photon_size_distribution["x"]
            yy = self._result_photon_size_distribution["y"]


            # #########################################################
            # # for plot, a Gaussian
            # xx = numpy.linspace(-5*s_phot,5*s_phot,101)
            # yy = numpy.exp(-xx**2/2/s_phot**2)
            # self._result_photon_size_distribution = {"x":xx,"y":yy}
            # #########################################################

            sampler_radial = Sampler1D(yy*numpy.abs(xx),xx)
            r,hy,hx = sampler_radial.get_n_sampled_points_and_histogram(NRAYS,bins=101)
            angle = numpy.random.random(NRAYS) * 2 * numpy.pi

            x_photon = r / numpy.sqrt(2.0) * numpy.sin(angle)
            y_photon = 0.0
            z_photon = r / numpy.sqrt(2.0) * numpy.cos(angle)


        rays[:,0] = x_photon + x_electron
        rays[:,1] = y_photon + y_electron
        rays[:,2] = z_photon + z_electron


        if user_unit_to_m != 1.0:
            rays[:,0] /= user_unit_to_m
            rays[:,1] /= user_unit_to_m
            rays[:,2] /= user_unit_to_m

        #
        # sample divergences (cols 4-6): the Shadow way
        #
        THETABM = sampled_theta
        PHI = sampled_phi
        A_Z = numpy.arcsin(numpy.sin(THETABM)*numpy.sin(PHI))
        A_X = numpy.arccos(numpy.cos(THETABM)/numpy.cos(A_Z))
        THETABM = A_Z
        PHI  = A_X
        # ! C Decide in which quadrant THETA and PHI are.
        myrand = numpy.random.random(NRAYS)
        THETABM[numpy.where(myrand < 0.5)] *= -1.0
        myrand = numpy.random.random(NRAYS)
        PHI[numpy.where(myrand < 0.5)] *= -1.0

        if self._FLAG_EMITTANCE:
            EBEAM1 = numpy.random.normal(loc=0.0,scale=sigmas[1],size=NRAYS)
            EBEAM3 = numpy.random.normal(loc=0.0,scale=sigmas[3],size=NRAYS)
            ANGLEX = EBEAM1 + PHI
            ANGLEV = EBEAM3 + THETABM
        else:
            ANGLEX = PHI # E_BEAM(1) + PHI
            ANGLEV =THETABM #  E_BEAM(3) + THETABM

        VX = numpy.tan(ANGLEX)
        VY = 1.0
        VZ = numpy.tan(ANGLEV)/numpy.cos(ANGLEX)
        VN = numpy.sqrt( VX*VX + VY*VY + VZ*VZ)
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

        DIREC = rays[:,3:6].copy()
        A_VEC = numpy.zeros_like(DIREC)
        A_VEC[:,0] = 1.0

        # ! C
        # ! C   Rotate A_VEC so that it will be perpendicular to DIREC and with the
        # ! C   right components on the plane.
        # ! C
        # CALL CROSS (A_VEC,DIREC,A_TEMP)
        A_TEMP = self._cross(A_VEC,DIREC)
        # CALL CROSS (DIREC,A_TEMP,A_VEC)
        A_VEC = self._cross(DIREC,A_TEMP)
        # CALL NORM (A_VEC,A_VEC)
        A_VEC = self._norm(A_VEC)
        # CALL CROSS (A_VEC,DIREC,AP_VEC)
        AP_VEC = self._cross(A_VEC,DIREC)
        # CALL NORM (AP_VEC,AP_VEC)
        AP_VEC = self._norm(AP_VEC)

        #
        # obtain polarization for each ray (interpolation)
        #


        if self._NG_E == 1: # 2D interpolation
            sampled_photon_energy = numpy.array(sampled_photon_energy) # be sure is an array
            fn = interpolate.RegularGridInterpolator(
                (self._result_radiation["theta"],self._result_radiation["phi"]),
                self._result_radiation["polarization"][0])

            pts = numpy.dstack( (sampled_theta,
                                 sampled_phi) )
            pts = pts[0]
            POL_DEG = fn(pts)
        else: # 3D interpolation
            fn = interpolate.RegularGridInterpolator(
                (self._result_radiation["photon_energy"],self._result_radiation["theta"],self._result_radiation["phi"]),
                self._result_radiation["polarization"])

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
            A_VEC[:,i] *= AX

        AZ = (1.0-POL_DEG)/DENOM
        for i in range(3):
            AP_VEC[:,i] *= AZ

        rays[:,6:9] =  A_VEC
        rays[:,15:18] = AP_VEC

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

        rays[:,13] = PHASEX
        rays[:,14] = PHASEZ

        # set flag (col 10)
        rays[:,9] = 1.0

        #
        # photon energy (col 11)
        #

        A2EV = 2.0*numpy.pi/(codata.h*codata.c/codata.e*1e2)
        rays[:,10] =  sampled_photon_energy * A2EV

        # col 12 (ray index)
        rays[:,11] =  1 + numpy.arange(NRAYS)

        # col 13 (optical path)
        rays[:,11] = 0.0

        return rays

    def _back_propagation_for_size_calculation(self,theta,radiation_flux,photon_energy,
                                                distance=100.0,magnification=0.010000):
        """
        Calculate the radiation_flux vs theta at a "distance"
        Back propagate to -distance
        The result is the size distrubution

        :param theta:
        :param radiation_flux:
        :param photon_energy:
        :param distance:
        :param magnification:
        :return: None; stores results in self._photon_size_distribution
        """

        from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
        from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
        from syned.beamline.beamline_element import BeamlineElement
        from syned.beamline.element_coordinates import ElementCoordinates
        from wofry.propagator.propagators1D.fresnel_zoom import FresnelZoom1D
        from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen1D


        input_wavefront = GenericWavefront1D().initialize_wavefront_from_arrays(theta*distance,numpy.sqrt(radiation_flux)+0j)
        input_wavefront.set_photon_energy(photon_energy)
        input_wavefront.set_spherical_wave(radius=distance,complex_amplitude=numpy.sqrt(radiation_flux)+0j)
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

        self._result_photon_size_distribution = {"x":output_wavefront.get_abscissas(),"y":output_wavefront.get_intensity()}


    def _cross(self,u,v):
        # w = u X v
        # u = array (npoints,vector_index)

        w = numpy.zeros_like(u)
        w[:,0] = u[:,1] * v[:,2] - u[:,2] * v[:,1]
        w[:,1] = u[:,2] * v[:,0] - u[:,0] * v[:,2]
        w[:,2] = u[:,0] * v[:,1] - u[:,1] * v[:,0]

        return w

    def _norm(self,u):
        # w = u / |u|
        # u = array (npoints,vector_index)
        u_norm = numpy.zeros_like(u)
        uu = numpy.sqrt( u[:,0]**2 + u[:,1]**2 + u[:,2]**2)
        for i in range(3):
            u_norm[:,i] = uu
        return u / u_norm

    def _sample_photon_energy_theta_and_phi(self,NRAYS):

        #
        # sample divergences
        #

        theta = self._result_radiation["theta"]
        phi = self._result_radiation["phi"]
        photon_energy = self._result_radiation["photon_energy"]

        photon_energy_spectrum = 'polychromatic' # 'monochromatic' #
        if self._EMIN == self._EMAX:
            photon_energy_spectrum = 'monochromatic'
        if self._NG_E == 1:
            photon_energy_spectrum = 'monochromatic'


        if photon_energy_spectrum == 'monochromatic':

            #2D case
            tmp = self._result_radiation["radiation"][0,:,:].copy()
            tmp /= tmp.max()

            # correct radiation for DxDz / DthetaDphi
            tmp_theta = numpy.outer(theta,numpy.ones_like(phi))
            tmp_theta /= tmp_theta.max()
            tmp_theta += 1e-6 # to avoid zeros
            tmp *= tmp_theta
            # plot_image(tmp_theta,theta,phi,aspect='auto')

            s2d = Sampler2D(tmp,theta,phi)
            sampled_theta,sampled_phi = s2d.get_n_sampled_points(NRAYS)

            sampled_photon_energy = self._EMIN

        elif photon_energy_spectrum == "polychromatic":
            #3D case
            tmp = self._result_radiation["radiation"].copy()
            tmp /= tmp.max()
            # correct radiation for DxDz / DthetaDphi
            tmp_theta = numpy.outer(theta,numpy.ones_like(phi))
            tmp_theta /= tmp_theta.max()
            tmp_theta += 1e-6 # to avoid zeros
            for i in range(tmp.shape[0]):
                tmp[i,:,:] *= tmp_theta

            s3d = Sampler3D(tmp,photon_energy,theta,phi)

            sampled_photon_energy,sampled_theta,sampled_phi = s3d.get_n_sampled_points(NRAYS)


        return sampled_photon_energy,sampled_theta,sampled_phi


