__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "12/01/2017"


import numpy
import scipy.constants as codata

# import Shadow

from SourceUndulatorFactory import SourceUndulatorFactory
from SourceUndulatorFactorySrw import SourceUndulatorFactorySrw
from SourceUndulatorFactoryPysru import SourceUndulatorFactoryPysru

from syned.storage_ring.magnetic_structures.undulator import Undulator
from syned.storage_ring.electron_beam import ElectronBeam

import Shadow
from inverse_method_sampler import Sampler2D, Sampler3D

from scipy import interpolate


class SourceUndulator(object):
    def __init__(self,name="",
                 syned_electron_beam=ElectronBeam(),
                 syned_undulator=Undulator(),
                 FLAG_EMITTANCE=0,           # Use emittance (0=No, 1=Yes)
                 FLAG_SIZE=0,                # 0=point,1=Gaussian,2=FT(Divergences)
                 EMIN=10000.0,               # Photon energy scan from energy (in eV)
                 EMAX=11000.0,               # Photon energy scan to energy (in eV)
                 NG_E=11,                    # Photon energy scan number of points
                 MAXANGLE=0.5,               # Maximum radiation semiaperture in mrad # TODO: define it in rad, for consistency
                 NG_T=31,                    # Number of points in angle theta
                 NG_P=21,                    # Number of points in angle phi
                 NG_J=20,                    # Number of points in electron trajectory (per period)
                 SEED=36255655452,           # Random seed
                 NRAYS=5000,                 # Number of rays
                 code_undul_phot="internal", # internal, pysru, srw
                 ):

        # # Machine
        self.syned_electron_beam = syned_electron_beam

        # # Undulator
        self.syned_undulator = syned_undulator

        self.FLAG_EMITTANCE  =  FLAG_EMITTANCE # Yes  # Use emittance (0=No, 1=Yes)
        self.FLAG_SIZE  =  FLAG_SIZE # 0=point,1=Gaussian,2=FT(Divergences)

        # Photon energy scan
        self.EMIN            = EMIN   # Photon energy scan from energy (in eV)
        self.EMAX            = EMAX   # Photon energy scan to energy (in eV)
        self.NG_E            = NG_E   # Photon energy scan number of points
        # Geometry
        self.MAXANGLE        = MAXANGLE   # Maximum radiation semiaperture in mrad # TODO: define it in rad, for consistency
        self.NG_T            = NG_T       # Number of points in angle theta
        self.NG_P            = NG_P       # Number of points in angle phi
        self.NG_J            = NG_J       # Number of points in electron trajectory (per period)
        # ray tracing
        self.SEED            = SEED   # Random seed
        self.NRAYS           = NRAYS  # Number of rays

        self.code_undul_phot = code_undul_phot

        # results of calculations

        self.result_radiation = None


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
        if self.FLAG_EMITTANCE:
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
        txt += "        central cone 'half' width [mrad]:     %10.6f %10.6f %10.6f   \n"%(\
                                                                1e3*self.syned_undulator.gaussian_central_cone_aperture(self.syned_electron_beam.gamma(),1),
                                                                1e3*self.syned_undulator.gaussian_central_cone_aperture(self.syned_electron_beam.gamma(),3),
                                                                1e3*self.syned_undulator.gaussian_central_cone_aperture(self.syned_electron_beam.gamma(),5))
        txt += "        first ring at [mrad]:                 %10.6f %10.6f %10.6f   \n"%(\
                                                                1e3*self.get_resonance_ring(1,1),
                                                                1e3*self.get_resonance_ring(3,1),
                                                                1e3*self.get_resonance_ring(5,1))

        txt += "-----------------------------------------------------\n"
        txt += "Sampling: \n"
        if self.NG_E == 1:
            txt += "        photon energy %f eV\n"%(self.EMIN)
        else:
            txt += "        photon energy from %10.3f eV to %10.3f eV\n"%(self.EMIN,self.EMAX)
        txt += "        number of points for the trajectory %d\n"%(self.NG_J)
        txt += "        number of energy points %d\n"%(self.NG_E)
        txt += "        maximum elevation angle %f mrad\n"%(self.MAXANGLE)
        txt += "        number of angular elevation points %d\n"%(self.NG_T)
        txt += "        number of angular azimuthal points %d\n"%(self.NG_P)
        txt += "        number of rays %d\n"%(self.NRAYS)
        txt += "        random seed %d\n"%(self.SEED)
        txt += "-----------------------------------------------------\n"

        if self.result_radiation is None:
            txt += "        radiation: NOT YET CALCULATED\n"
        else:
            txt += "        radiation: CALCULATED\n"

        txt += "-----------------------------------------------------\n"
        return txt

    def get_resonance_ring(self,harmonic_number=1, ring_order=1):
        return 1.0/self.syned_electron_beam.gamma()*numpy.sqrt( ring_order / harmonic_number * (1+0.5*self.syned_undulator.K_vertical()**2) )

    # def set_harmonic(self,harmonic):

    def set_energy_monochromatic_at_resonance(self,harmonic_number):

        self.set_energy_monochromatic(self.syned_undulator.resonance_energy(
            self.syned_electron_beam.gamma(),harmonic=harmonic_number))
        # take 3*sigma - MAXANGLE is in mrad!!

        # self.MAXANGLE = 3 * 0.69 * 1e3 * self.get_resonance_central_cone(harmonic_number)
        self.MAXANGLE = 3 * 0.69 * 1e3 * self.syned_undulator.gaussian_central_cone_aperture(self.syned_electron_beam.gamma(),harmonic_number)

    def set_energy_monochromatic(self,emin):
        """
        Sets a single energy line for the source (monochromatic)
        :param emin: the energy in eV
        :return:
        """
        self.EMIN = emin
        self.EMAX = emin
        self.NG_E = 1


    def set_energy_box(self,emin,emax,npoints=None):
        """
        Sets a box energy distribution for the source (monochromatic)
        :param emin:  Photon energy scan from energy (in eV)
        :param emax:  Photon energy scan to energy (in eV)
        :param npoints:  Photon energy scan number of points (optinal, if not set no changes)
        :return:
        """

        self.EMIN = emin
        self.EMAX = emax
        if npoints != None:
            self.NG_E = npoints

    def get_radiation_polar(self):
        if self.result_radiation is None:
            self.calculate_radiation()
        return self.result_radiation["radiation"],self.result_radiation["theta"],self.result_radiation["phi"]

    def get_radiation_interpolated_cartesian(self,npointsx=100,npointsz=100,thetamax=None):

        radiation,thetabm,phi = self.get_radiation_polar()

        if thetamax is None:
            thetamax = thetabm.max()

        vx = numpy.linspace(-1.1*thetamax,1.1*thetamax,npointsx)
        vz = numpy.linspace(-1.1*thetamax,1.1*thetamax,npointsz)
        VX = numpy.outer(vx,numpy.ones_like(vz))
        VZ = numpy.outer(numpy.ones_like(vx),vz)
        VY = numpy.sqrt(1 - VX**2 - VZ**2)

        THETA = numpy.arctan( numpy.sqrt(VX**2+VZ**2)/VY)
        PHI = numpy.arctan(VZ/VX)

        radiation_interpolated = numpy.zeros((radiation.shape[0],npointsx,npointsz))

        for i in range(radiation.shape[0]):
            interpolator_value = interpolate.RectBivariateSpline(thetabm, phi, radiation[i])
            radiation_interpolated[i] = interpolator_value.ev(THETA, PHI)

        return radiation_interpolated,vx,vz


    def calculate_radiation(self):

        """
        Calculates the radiation (emission) as a function pf theta (elevation angle) and phi (azimuthal angle)
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


        self.result_radiation = None

        # undul_phot
        if self.code_undul_phot == 'internal':
            undul_phot_dict = SourceUndulatorFactory.undul_phot(E_ENERGY  = self.syned_electron_beam.energy(),
                                         INTENSITY = self.syned_electron_beam.current(),
                                         LAMBDAU   = self.syned_undulator.period_length(),
                                         NPERIODS  = self.syned_undulator.number_of_periods(),
                                         K         = self.syned_undulator.K(),
                                         EMIN      = self.EMIN,
                                         EMAX      = self.EMAX,
                                         NG_E      = self.NG_E,
                                         MAXANGLE  = self.MAXANGLE,
                                         NG_T      = self.NG_T,
                                         NG_P      = self.NG_P,
                                         number_of_trajectory_points = self.NG_J)

        elif self.code_undul_phot == 'pysru':
            undul_phot_dict = SourceUndulatorFactoryPysru.undul_phot(E_ENERGY  = self.syned_electron_beam.energy(),
                                         INTENSITY = self.syned_electron_beam.current(),
                                         LAMBDAU   = self.syned_undulator.period_length(),
                                         NPERIODS  = self.syned_undulator.number_of_periods(),
                                         K         = self.syned_undulator.K(),
                                         EMIN      = self.EMIN,
                                         EMAX      = self.EMAX,
                                         NG_E      = self.NG_E,
                                         MAXANGLE  = self.MAXANGLE,
                                         NG_T      = self.NG_T,
                                         NG_P      = self.NG_P,)
        elif self.code_undul_phot == 'srw':
            undul_phot_dict = SourceUndulatorFactorySrw.undul_phot(E_ENERGY  = self.syned_electron_beam.energy(),
                                         INTENSITY = self.syned_electron_beam.current(),
                                         LAMBDAU   = self.syned_undulator.period_length(),
                                         NPERIODS  = self.syned_undulator.number_of_periods(),
                                         K         = self.syned_undulator.K(),
                                         EMIN      = self.EMIN,
                                         EMAX      = self.EMAX,
                                         NG_E      = self.NG_E,
                                         MAXANGLE  = self.MAXANGLE,
                                         NG_T      = self.NG_T,
                                         NG_P      = self.NG_P,)
        else:
            raise Exception("Not implemented undul_phot code: "+code_undul_phot)

        # add some info
        undul_phot_dict["code_undul_phot"] = self.code_undul_phot
        undul_phot_dict["info"] = self.info()

        self.result_radiation = undul_phot_dict
        # return undul_phot_dict


    def calculate_shadow3_beam(self,user_unit_to_m=1.0,F_COHER=1):

        if self.result_radiation is None:
            self.calculate_radiation()

        sampled_photon_energy,sampled_theta,sampled_phi = self._sample_photon_beam()

        beam = self._sample_shadow3_beam(sampled_photon_energy,sampled_theta,sampled_phi,F_COHER=F_COHER)

        if user_unit_to_m != 1.0:
            beam.rays[:,0] /= user_unit_to_m
            beam.rays[:,1] /= user_unit_to_m
            beam.rays[:,2] /= user_unit_to_m

        return beam


    def _sample_shadow3_beam(self,sampled_photon_energy,sampled_theta,sampled_phi,F_COHER=1):

        beam = Shadow.Beam(N=self.NRAYS)

        sigmas = self.syned_electron_beam.get_sigmas_all()

        #
        # sample sizes (cols 1-3)
        #
        if self.FLAG_EMITTANCE:
            x_electron = numpy.random.normal(loc=0.0,scale=sigmas[0],size=self.NRAYS)
            y_electron = 0.0
            z_electron = numpy.random.normal(loc=0.0,scale=sigmas[2],size=self.NRAYS)
        else:
            x_electron = 0.0
            y_electron = 0.0
            z_electron = 0.0

        if self.FLAG_SIZE == 0:
            x_photon = 0.0
            y_photon = 0.0
            z_photon = 0.0
        elif self.FLAG_SIZE == 1:
            undulator_length = self.syned_undulator.length()
            lambda1 = codata.h*codata.c/codata.e / sampled_photon_energy.mean()

            # calculate sizes of the photon undulator beam
            # see formulas 25 & 30 in Elleaume (Onaki & Elleaume)
            # sp_phot = 0.69*numpy.sqrt(lambda1/undulator_length)
            s_phot = 2.740/(4e0*numpy.pi)*numpy.sqrt(undulator_length*lambda1)

            cov = [[s_phot**2, 0], [0, s_phot**2]]
            mean = [0.0,0.0]

            tmp = numpy.random.multivariate_normal(mean, cov, self.NRAYS)
            x_photon = tmp[:,0]
            y_photon = 0.0
            z_photon = tmp[:,1]
        elif self.FLAG_SIZE == 2:
            raise Exception("To be implemented")


        beam.rays[:,0] = x_photon + x_electron
        beam.rays[:,1] = y_photon + y_electron
        beam.rays[:,2] = z_photon + z_electron


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
        myrand = numpy.random.random(self.NRAYS)
        THETABM[numpy.where(myrand < 0.5)] *= -1.0
        myrand = numpy.random.random(self.NRAYS)
        PHI[numpy.where(myrand < 0.5)] *= -1.0

        if self.FLAG_EMITTANCE:
            EBEAM1 = numpy.random.normal(loc=0.0,scale=sigmas[1],size=self.NRAYS)
            EBEAM3 = numpy.random.normal(loc=0.0,scale=sigmas[3],size=self.NRAYS)
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

        beam.rays[:,3] = VX
        beam.rays[:,4] = VY
        beam.rays[:,5] = VZ


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

        DIREC = beam.rays[:,3:6].copy()
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


        if self.NG_E == 1: # 2D interpolation
            sampled_photon_energy = numpy.array(sampled_photon_energy) # be sure is an array
            fn = interpolate.RegularGridInterpolator(
                (self.result_radiation["theta"],self.result_radiation["phi"]),
                self.result_radiation["polarization"][0])

            pts = numpy.dstack( (sampled_theta,
                                 sampled_phi) )
            pts = pts[0]
            POL_DEG = fn(pts)
        else: # 3D interpolation
            fn = interpolate.RegularGridInterpolator(
                (self.result_radiation["photon_energy"],self.result_radiation["theta"],self.result_radiation["phi"]),
                self.result_radiation["polarization"])

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

        beam.rays[:,6:9] =  A_VEC
        beam.rays[:,15:18] = AP_VEC

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
            PHASEX = numpy.random.random(self.NRAYS) * 2 * numpy.pi

        PHASEZ = PHASEX + POL_ANGLE * numpy.sign(ANGLEV)

        beam.rays[:,13] = PHASEX
        beam.rays[:,14] = PHASEZ

        # set flag (col 10)
        beam.rays[:,9] = 1.0

        #
        # photon energy (col 11)
        #

        A2EV = 2.0*numpy.pi/(codata.h*codata.c/codata.e*1e2)
        beam.rays[:,10] =  sampled_photon_energy * A2EV

        # col 12 (ray index)
        beam.rays[:,11] =  1 + numpy.arange(self.NRAYS)

        # col 13 (optical path)
        beam.rays[:,11] = 0.0

        return beam

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

    def _sample_photon_beam(self):

        #
        # sample divergences
        #

        theta = self.result_radiation["theta"]
        phi = self.result_radiation["phi"]
        photon_energy = self.result_radiation["photon_energy"]

        photon_energy_spectrum = 'polychromatic' # 'monochromatic' #
        if self.EMIN == self.EMAX:
            photon_energy_spectrum = 'monochromatic'
        if self.NG_E == 1:
            photon_energy_spectrum = 'monochromatic'


        if photon_energy_spectrum == 'monochromatic':

            #2D case
            tmp = self.result_radiation["radiation"][0,:,:].copy()
            tmp /= tmp.max()

            # correct radiation for DxDz / DthetaDphi
            tmp_theta = numpy.outer(theta,numpy.ones_like(phi))
            tmp_theta /= tmp_theta.max()
            tmp_theta += 1e-6 # to avoid zeros
            tmp *= tmp_theta
            # plot_image(tmp_theta,theta,phi,aspect='auto')

            s2d = Sampler2D(tmp,theta,phi)
            sampled_theta,sampled_phi = s2d.get_n_sampled_points(self.NRAYS)

            sampled_photon_energy = self.EMIN

        elif photon_energy_spectrum == "polychromatic":
            #3D case
            tmp = self.result_radiation["radiation"].copy()
            tmp /= tmp.max()
            # correct radiation for DxDz / DthetaDphi
            tmp_theta = numpy.outer(theta,numpy.ones_like(phi))
            tmp_theta /= tmp_theta.max()
            tmp_theta += 1e-6 # to avoid zeros
            for i in range(tmp.shape[0]):
                tmp[i,:,:] *= tmp_theta

            s3d = Sampler3D(tmp,photon_energy,theta,phi)

            sampled_photon_energy,sampled_theta,sampled_phi = s3d.get_n_sampled_points(self.NRAYS)


        return sampled_photon_energy,sampled_theta,sampled_phi


