__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "30-08-2018"

"""

Wiggler code: computes wiggler radiation distributions and samples rays according to them.

Fully replaces and upgrades the shadow3 wiggler model.

The radiation is calculating using sr-xraylib


Usage:

sw = SourceWiggler()        # use keywords to define parameters. It uses syned for electron beam and wiggler
rays = sw.calculate_rays()    # sample rays. Result is a numpy.array of shape (NRAYS,18), exactly the same as in shadow3


"""


import numpy

from srxraylib.util.inverse_method_sampler import Sampler1D
from srxraylib.plot.gol import plot_scatter
import scipy.constants as codata

from syned.storage_ring.magnetic_structures.wiggler import Wiggler
from shadow4.syned.magnetic_structure_1D_field import MagneticStructure1DField

from syned.storage_ring.electron_beam import ElectronBeam


from srxraylib.sources.srfunc import wiggler_trajectory, wiggler_spectrum, wiggler_cdf, sync_f

import scipy
from scipy.interpolate import interp1d


# This is similar than sync_f in srxraylib but faster
def sync_f_sigma_and_pi(rAngle, rEnergy):
    r""" angular dependency of synchrotron radiation emission

      NAME:
            sync_f_sigma_and_pi

      PURPOSE:
            Calculates the function used for calculating the angular
         dependence of synchrotron radiation.

      CATEGORY:
            Mathematics.

      CALLING SEQUENCE:
            Result = sync_f_sigma_and_pi(rAngle,rEnergy)

      INPUTS:
            rAngle:  (array) the reduced angle, i.e., angle[rads]*Gamma. It can be a
             scalar or a vector.
            rEnergy:  (scalar) a value for the reduced photon energy, i.e.,
             energy/critical_energy.

      KEYWORD PARAMETERS:


      OUTPUTS:
            returns the value  of the sync_f for sigma and pi polarizations
             The result is an array of the same dimension as rAngle.

      PROCEDURE:
            The number of emitted photons versus vertical angle Psi is
         proportional to sync_f, which value is given by the formulas
         in the references.


         References:
             G K Green, "Spectra and optics of synchrotron radiation"
                 BNL 50522 report (1976)
             A A Sokolov and I M Ternov, Synchrotron Radiation,
                 Akademik-Verlag, Berlin, 1968

      OUTPUTS:
            returns the value  of the sync_f function

      PROCEDURE:
            Uses BeselK() function

      MODIFICATION HISTORY:
            Written by:     M. Sanchez del Rio, srio@esrf.fr, 2002-05-23
         2002-07-12 srio@esrf.fr adds circular polarization term for
             wavelength integrated spectrum (S&T formula 5.25)
         2012-02-08 srio@esrf.eu: python version
         2019-10-31 srio@lbl.gov  speed-up changes for shadow4

    """

    #
    # ; For 11 in Pag 6 in Green 1975
    #
    ji = numpy.sqrt((1.0 + rAngle**2)**3) * rEnergy / 2.0
    efe_sigma = scipy.special.kv(2.0 / 3.0, ji) * (1.0 + rAngle**2)
    efe_pi = rAngle * scipy.special.kv(1.0 / 3.0, ji) / numpy.sqrt(1.0 + rAngle ** 2) * (1.0 + rAngle ** 2)
    return efe_sigma**2,efe_pi**2



class SourceWiggler(object):
    def __init__(self,name="",
                 syned_electron_beam=None,
                 syned_wiggler=None,
                 emin=10000.0,               # Photon energy scan from energy (in eV)
                 emax=11000.0,               # Photon energy scan to energy (in eV)
                 ng_e=11,                    # Photon energy scan number of points
                 ng_j=20,                    # Number of points in electron trajectory (per period) for internal calculation only
                 flag_emittance=0,           # when sampling rays: Use emittance (0=No, 1=Yes)
                 ):

        # # Machine
        if syned_electron_beam is None:
            self.syned_electron_beam = ElectronBeam()
        else:
            self.syned_electron_beam = syned_electron_beam

        # # Undulator
        if syned_wiggler is None:
            self.syned_wiggler = Wiggler()
        else:
            self.syned_wiggler = syned_wiggler

        # Photon energy scan
        self._EMIN            = emin   # Photon energy scan from energy (in eV)
        self._EMAX            = emax   # Photon energy scan to energy (in eV)
        self._NG_E            = ng_e   # Photon energy scan number of points

        self._NG_J            = ng_j       # Number of points in electron trajectory (per period)

        self._FLAG_EMITTANCE  =  flag_emittance # Yes  # Use emittance (0=No, 1=Yes)

        # electron initial conditions for electron trahjectory calculations
        self.set_electron_initial_conditions()

        # results of calculations

        self._result_trajectory = None
        self._result_parameters = None
        self._result_cdf = None


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

        if isinstance(self.syned_wiggler,Wiggler):  # conventional wiggler
            txt += "Input Wiggler parameters: \n"
            txt += "        period: %f m\n"%self.syned_wiggler.period_length()
            txt += "        number of periods: %d\n"%self.syned_wiggler.number_of_periods()
            txt += "        K-value: %f\n"%self.syned_wiggler.K_vertical()

            txt += "-----------------------------------------------------\n"

            txt += "Lorentz factor (gamma): %f\n"%self.syned_electron_beam.gamma()
            txt += "Electron velocity: %.12f c units\n"%(numpy.sqrt(1.0 - 1.0 / self.syned_electron_beam.gamma() ** 2))
            txt += "Undulator length: %f m\n"%(self.syned_wiggler.period_length()*self.syned_wiggler.number_of_periods())
            K_to_B = (2.0 * numpy.pi / self.syned_wiggler.period_length()) * codata.m_e * codata.c / codata.e

            txt += "Wiggler peak magnetic field: %f T\n"%(K_to_B*self.syned_wiggler.K_vertical())

        if isinstance(self.syned_wiggler,MagneticStructure1DField):  # conventional wiggler
            txt += "Input Wiggler parameters: \n"
            txt += "        from external magnetic field \n"
            txt += self.syned_wiggler.info()


        # txt += "Resonances: \n"
        # txt += "        harmonic number [n]                   %10d %10d %10d \n"%(1,3,5)
        # txt += "        wavelength [A]:                       %10.6f %10.6f %10.6f   \n"%(\
        #                                                         1e10*self.syned_undulator.resonance_wavelength(self.syned_electron_beam.gamma(),harmonic=1),
        #                                                         1e10*self.syned_undulator.resonance_wavelength(self.syned_electron_beam.gamma(),harmonic=3),
        #                                                         1e10*self.syned_undulator.resonance_wavelength(self.syned_electron_beam.gamma(),harmonic=5))
        # txt += "        energy [eV]   :                       %10.3f %10.3f %10.3f   \n"%(\
        #                                                         self.syned_undulator.resonance_energy(self.syned_electron_beam.gamma(),harmonic=1),
        #                                                         self.syned_undulator.resonance_energy(self.syned_electron_beam.gamma(),harmonic=3),
        #                                                         self.syned_undulator.resonance_energy(self.syned_electron_beam.gamma(),harmonic=5))
        # txt += "        frequency [Hz]:                       %10.3g %10.3g %10.3g   \n"%(\
        #                                                         1e10*self.syned_undulator.resonance_frequency(self.syned_electron_beam.gamma(),harmonic=1),
        #                                                         1e10*self.syned_undulator.resonance_frequency(self.syned_electron_beam.gamma(),harmonic=3),
        #                                                         1e10*self.syned_undulator.resonance_frequency(self.syned_electron_beam.gamma(),harmonic=5))
        # txt += "        central cone 'half' width [urad]:     %10.6f %10.6f %10.6f   \n"%(\
        #                                                         1e6*self.syned_undulator.gaussian_central_cone_aperture(self.syned_electron_beam.gamma(),1),
        #                                                         1e6*self.syned_undulator.gaussian_central_cone_aperture(self.syned_electron_beam.gamma(),3),
        #                                                         1e6*self.syned_undulator.gaussian_central_cone_aperture(self.syned_electron_beam.gamma(),5))
        # txt += "        first ring at [urad]:                 %10.6f %10.6f %10.6f   \n"%(\
        #                                                         1e6*self.get_resonance_ring(1,1),
        #                                                         1e6*self.get_resonance_ring(3,1),
        #                                                         1e6*self.get_resonance_ring(5,1))
        #
        txt += "-----------------------------------------------------\n"
        txt += "Grids: \n"
        if self._NG_E == 1:
            txt += "        photon energy %f eV\n"%(self._EMIN)
        else:
            txt += "        photon energy from %10.3f eV to %10.3f eV\n"%(self._EMIN,self._EMAX)
        txt += "        number of energy points: %d\n"%(self._NG_E)
        txt += "        number of points for the trajectory: %d\n"%(self._NG_J)
        # txt += "        maximum elevation angle: %f urad\n"%(1e6*self._MAXANGLE)
        # txt += "        number of angular elevation points: %d\n"%(self._NG_T)
        # txt += "        number of angular azimuthal points: %d\n"%(self._NG_P)
        # # txt += "        number of rays: %d\n"%(self.NRAYS)
        # # txt += "        random seed: %d\n"%(self.SEED)
        # txt += "-----------------------------------------------------\n"
        #
        # txt += "calculation code: %s\n"%self.code_undul_phot
        # if self._result_radiation is None:
        #     txt += "radiation: NOT YET CALCULATED\n"
        # else:
        #     txt += "radiation: CALCULATED\n"
        # txt += "Sampling: \n"
        # if self._FLAG_SIZE == 0:
        #     flag = "point"
        # elif self._FLAG_SIZE == 1:
        #     flag = "Gaussian"
        # txt += "        sampling flag: %d (%s)\n"%(self._FLAG_SIZE,flag)

        txt += "-----------------------------------------------------\n"
        return txt


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


    def set_electron_initial_conditions(self,shift_x_flag=0,shift_x_value=0.0,shift_betax_flag=0,shift_betax_value=0.0):
        self.shift_x_flag      = shift_x_flag
        self.shift_x_value     = shift_x_value
        self.shift_betax_flag  = shift_betax_flag
        self.shift_betax_value = shift_betax_value

    def set_electron_initial_conditions_by_label(self,
                                        position_label="no_shift", # values are: no_shift, half_excursion, minimum, maximum, value_at_zero, user_value
                                        velocity_label="no_shift", # values are: no_shift, half_excursion, minimum, maximum, value_at_zero, user_value
                                        position_value=0.0,
                                        velocity_value=0.0,
                                        ):
        self.shift_x_value = 0.0
        self.shift_betax_value = 0.0

        if position_label == "no_shift":
            self.shift_x_flag = 0
        elif position_label == "half_excursion":
            self.shift_x_flag = 1
        elif position_label == "minimum":
            self.shift_x_flag = 2
        elif position_label == "maximum":
            self.shift_x_flag = 3
        elif position_label == "value_at_zero":
            self.shift_x_flag = 4
        elif position_label == "user_value":
            self.shift_x_flag = 5
            self.position_value = position_value
        else:
            raise Exception("Invalid value for keyword position_label")

        if velocity_label == "no_shift":
            self.shift_betax_flag = 0
        elif velocity_label == "half_excursion":
            self.shift_betax_flag = 1
        elif velocity_label == "minimum":
            self.shift_betax_flag = 2
        elif velocity_label == "maximum":
            self.shift_betax_flag = 3
        elif velocity_label == "value_at_zero":
            self.shift_betax_flag = 4
        elif velocity_label == "user_value":
            self.shift_betax_flag = 5
            self.shift_betax_value = velocity_value
        else:
            raise Exception("Invalid value for keyword velocity_label")

    def get_energy_box(self):
        """
        Gets the limits of photon energy distribution for the source
        :return: emin,emax,number_of_points
        """
        return self._EMIN,self._EMAX,self._NG_E


    def calculate_radiation(self):


        if isinstance(self.syned_wiggler,Wiggler):

            (traj, pars) = wiggler_trajectory(b_from=0,
                                                     inData="",
                                                     nPer=self.syned_wiggler.number_of_periods(),
                                                     nTrajPoints=self._NG_J,
                                                     ener_gev=self.syned_electron_beam._energy_in_GeV,
                                                     per=self.syned_wiggler.period_length(),
                                                     kValue=self.syned_wiggler.K_vertical(),
                                                     trajFile="",)

        elif isinstance(self.syned_wiggler,MagneticStructure1DField):

            print(">>>>>>>>>>>>>>>>>>>>>>",
                "shift_x_flag =      ",self.shift_x_flag,
                "shift_x_value =     ",self.shift_x_value,
                "shift_betax_flag =  ",self.shift_betax_flag,
                "shift_betax_value = ",self.shift_betax_value
                  )




            inData = numpy.vstack((self.syned_wiggler.get_abscissas(),self.syned_wiggler.get_magnetic_field())).T

            (traj, pars) = wiggler_trajectory(b_from=1,
                                                     inData=inData,
                                                     nPer=1,
                                                     nTrajPoints=self._NG_J,
                                                     ener_gev=self.syned_electron_beam._energy_in_GeV,
                                                     # per=self.syned_wiggler.period_length(),
                                                     # kValue=self.syned_wiggler.K_vertical(),
                                                     trajFile="",
                                                     shift_x_flag       = self.shift_x_flag     ,
                                                     shift_x_value      = self.shift_x_value    ,
                                                     shift_betax_flag   = self.shift_betax_flag ,
                                                     shift_betax_value  = self.shift_betax_value,)


        self._result_trajectory = traj
        self._result_parameters = pars

        # print(">>>>>>>>>> traj pars: ",traj.shape,pars)
        #
        # plot(traj[1, :], traj[0, :], xtitle="Y", ytitle="X")
        # plot(traj[1, :], traj[3, :], xtitle="Y", ytitle="BetaX")
        # plot(traj[1, :], traj[6, :], xtitle="Y", ytitle="Curvature")
        # plot(traj[1, :], traj[7, :], xtitle="Y", ytitle="B")


        # traj[0,ii] = yx[i]
        # traj[1,ii] = yy[i]+j * per - start_len
        # traj[2,ii] = 0.0
        # traj[3,ii] = betax[i]
        # traj[4,ii] = betay[i]
        # traj[5,ii] = 0.0
        # traj[6,ii] = curv[i]
        # traj[7,ii] = bz[i]


        #
        # calculate cdf and write file for Shadow/Source
        #

        print(">>>>>>>>>>>>>>>>>>>>  self._EMIN,self._EMAX,self._NG_E",self._EMIN,self._EMAX,self._NG_E)
        self._result_cdf = wiggler_cdf(self._result_trajectory,
                           enerMin=self._EMIN,
                           enerMax=self._EMAX,
                           enerPoints=self._NG_E,
                           outFile="tmp.cdf",
                           elliptical=False)


    def calculate_rays(self,user_unit_to_m=1.0,F_COHER=0,NRAYS=5000,SEED=123456,EPSI_DX=0.0,EPSI_DZ=0.0,
                       psi_interval_in_units_one_over_gamma=None,
                       psi_interval_number_of_points=1001,
                       verbose=True):
        """
        compute the rays in SHADOW matrix (shape (npoints,18) )
        :param F_COHER: set this flag for coherent beam
        :param user_unit_to_m: default 1.0 (m)
        :return: rays, a numpy.array((npoits,18))
        """

        if self._result_cdf is None:
            self.calculate_radiation()

        if verbose:
            print(">>>   Results of calculate_radiation")
            print(">>>       trajectory.shape: ",self._result_trajectory.shape)
            print(">>>       cdf: ", self._result_cdf.keys())


        sampled_photon_energy,sampled_theta,sampled_phi = self._sample_photon_energy_theta_and_phi(NRAYS)

        if verbose:
            print(">>> sampled sampled_photon_energy,sampled_theta,sampled_phi:  ",sampled_photon_energy,sampled_theta,sampled_phi)

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

        # traj[0,ii] = yx[i]
        # traj[1,ii] = yy[i]+j * per - start_len
        # traj[2,ii] = 0.0
        # traj[3,ii] = betax[i]
        # traj[4,ii] = betay[i]
        # traj[5,ii] = 0.0
        # traj[6,ii] = curv[i]
        # traj[7,ii] = bz[i]

        PATH_STEP = self._result_cdf["step"]
        X_TRAJ = self._result_cdf["x"]
        Y_TRAJ = self._result_cdf["y"]
        SEEDIN = self._result_cdf["cdf"]
        ANGLE  = self._result_cdf["angle"]
        CURV   = self._result_cdf["curv"]
        EPSI_PATH = numpy.arange(CURV.size) * PATH_STEP # self._result_trajectory[7,:]


        # ! C We define the 5 arrays:
        # ! C    Y_X(5,N)    ---> X(Y)
        # ! C    Y_XPRI(5,N) ---> X'(Y)
        # ! C    Y_CURV(5,N) ---> CURV(Y)
        # ! C    Y_PATH(5,N) ---> PATH(Y)
        # ! C    F(1,N) contains the array of Y values where the nodes are located.


        # CALL PIECESPL(SEED_Y, Y_TEMP,   NP_SY,   IER)
        # CALL CUBSPL (Y_X,    X_TEMP,   NP_TRAJ, IER)
        # CALL CUBSPL (Y_Z,    Z_TEMP,   NP_TRAJ, IER)
        # CALL CUBSPL (Y_XPRI, ANG_TEMP, NP_TRAJ, IER)
        # CALL CUBSPL (Y_ZPRI, ANG2_TEMP, NP_TRAJ, IER)
        # CALL CUBSPL (Y_CURV, C_TEMP,   NP_TRAJ, IER)
        # CALL CUBSPL (Y_PATH, P_TEMP,   NP_TRAJ, IER)

        SEED_Y = interp1d(SEEDIN,Y_TRAJ,kind='linear')
        Y_X    = interp1d(Y_TRAJ,X_TRAJ,kind='cubic')
        Y_XPRI = interp1d(Y_TRAJ,ANGLE,kind='cubic')
        Y_CURV = interp1d(Y_TRAJ,CURV,kind='cubic')
        Y_PATH = interp1d(Y_TRAJ,EPSI_PATH,kind='cubic')

        # ! C+++
        # ! C Compute the path length to the middle (origin) of the wiggler.
        # ! C We need to know the "center" of the wiggler coordinate.
        # ! C input:     Y_PATH  ---> spline array
        # ! C            NP_TRAJ ---> # of points
        # ! C            Y_TRAJ  ---> calculation point (ind. variable)
        # ! C output:    PATH0   ---> value of Y_PATH at X = Y_TRAJ. If
        # ! C                         Y_TRAJ = 0, then PATH0 = 1/2 length
        # ! C                         of trajectory.
        # ! C+++

        Y_TRAJ = 0.0
        # CALL SPL_INT (Y_PATH, NP_TRAJ, Y_TRAJ, PATH0, IER)
        PATH0 = Y_PATH(Y_TRAJ)


        # ! C
        # ! C These flags are set because of the original program structure.
        # ! C
        # F_PHOT  = 0
        # F_COLOR  = 3
        # FSOUR  = 3
        # FDISTR  = 4

        ws_ev,ws_f,tmp =  wiggler_spectrum(self._result_trajectory,
                                    enerMin=self._EMIN, enerMax=self._EMAX, nPoints=500,
                                    # per=self.syned_wiggler.period_length(),
                                    electronCurrent=self.syned_electron_beam._current,
                                    outFile="", elliptical=False)

        ws_flux_per_ev = ws_f / (ws_ev*1e-3)
        samplerE = Sampler1D(ws_flux_per_ev,ws_ev)

        sampled_energies,h,h_center = samplerE.get_n_sampled_points_and_histogram(NRAYS)


        ###############################################

        gamma = self.syned_electron_beam.gamma()
        m2ev = codata.c * codata.h / codata.e
        TOANGS = m2ev * 1e10


        #####################################################

        RAD_MIN = 1.0 / numpy.abs(self._result_cdf["curv"]).max()

        critical_energy = TOANGS * 3.0 * numpy.power(gamma, 3) / 4.0 / numpy.pi / 1.0e10 * (1.0 / RAD_MIN)

        if psi_interval_in_units_one_over_gamma is None:
            c = numpy.array([-0.3600382, 0.11188709])  # see file fit_psi_interval.py
            # x = numpy.log10(self._EMIN / critical_energy)
            x = numpy.log10(self._EMIN / (4 * critical_energy)) # the wiggler that does not have an unique
                                                                # Ec. To be safe, I use 4 times the
                                                                # Ec vale to make the interval wider than for the BM
            y_fit = c[1] + c[0] * x
            psi_interval_in_units_one_over_gamma = 10 ** y_fit  # this is the semi interval
            psi_interval_in_units_one_over_gamma *= 4  # doubled interval
            if psi_interval_in_units_one_over_gamma < 2:
                psi_interval_in_units_one_over_gamma = 2

        if verbose:
            print(">>> psi_interval_in_units_one_over_gamma: ",psi_interval_in_units_one_over_gamma)

        angle_array_mrad = numpy.linspace(-0.5*psi_interval_in_units_one_over_gamma * 1e3 / gamma,
                                          0.5*psi_interval_in_units_one_over_gamma * 1e3 / gamma,
                                          psi_interval_number_of_points)


        # a = numpy.linspace(-0.6,0.6,150)

        a = angle_array_mrad

        #####################################################################



        a8 = 1.0
        hdiv_mrad = 1.0
        # i_a = self.syned_electron_beam._current
        #
        # fm = sync_f(a*self.syned_electron_beam.gamma()/1e3,eene,polarization=0) * \
        #         numpy.power(eene,2)*a8*i_a*hdiv_mrad*numpy.power(self.syned_electron_beam._energy_in_GeV,2)
        #
        # plot(a,fm,title="sync_f")
        #
        # samplerAng = Sampler1D(fm,a)
        #
        # sampled_theta,hx,h = samplerAng.get_n_sampled_points_and_histogram(10*NRAYS)
        # plot(h,hx)



        for itik in range(NRAYS):

            #     ARG_Y = GRID(2,ITIK)
            #     CALL SPL_INT (SEED_Y, NP_SY,   ARG_Y,  Y_TRAJ,    IER)
            arg_y = numpy.random.random() # ARG_Y[itik]
            Y_TRAJ = SEED_Y(arg_y)


            #     ! srio@esrf.eu 2014-05-19
            #     ! in wiggler some problems arise because spl_int
            #     ! does not return a Y value in the correct range.
            #     ! In those cases, we make a linear interpolation instead.
            #     if ((y_traj.le.y_temp(1)).or.(y_traj.gt.y_temp(NP_SY))) then
            #         y_traj_old = y_traj
            #         CALL LIN_INT (SEED_Y, NP_SY,   ARG_Y,  Y_TRAJ,    IER)
            #         print*,'SOURCESYNC: bad y_traj from SPL_INT, corrected with LIN_SPL: ',y_traj_old,'=>',y_traj
            #     endif
            #
            #     CALL SPL_INT (Y_X,    NP_TRAJ, Y_TRAJ, X_TRAJ,    IER)
            #     CALL SPL_INT (Y_XPRI, NP_TRAJ, Y_TRAJ, ANGLE,     IER)
            #     CALL SPL_INT (Y_CURV, NP_TRAJ, Y_TRAJ, CURV,      IER)
            #     CALL SPL_INT (Y_PATH, NP_TRAJ, Y_TRAJ, EPSI_PATH, IER)
            # END IF

            X_TRAJ = Y_X(Y_TRAJ)
            ANGLE = Y_XPRI(Y_TRAJ)
            CURV = Y_CURV(Y_TRAJ)
            EPSI_PATH = Y_PATH(Y_TRAJ)

            # print("\n>>><<<",arg_y,Y_TRAJ,X_TRAJ,ANGLE,CURV,EPSI_PATH)


            # EPSI_PATH = EPSI_PATH - PATH0 ! now refer to wiggler's origin
            # IF (CURV.LT.0) THEN
            #     POL_ANGLE = 90.0D0  ! instant orbit is CW
            # ELSE
            #     POL_ANGLE = -90.0D0  !     CCW
            # END IF
            # IF (CURV.EQ.0) THEN
            #     R_MAGNET = 1.0D+20
            # ELSE
            #     R_MAGNET = ABS(1.0D0/CURV)
            # END IF
            # POL_ANGLE  = TORAD*POL_ANGLE

            EPSI_PATH = EPSI_PATH - PATH0 # now refer to wiggler's origin
            if CURV < 0:
                POL_ANGLE = 90.0 # instant orbit is CW
            else:
                POL_ANGLE = -90.0 # CCW

            if CURV == 0.0:
                R_MAGNET = 1.0e20
            else:
                R_MAGNET = numpy.abs(1.0/CURV)

            POL_ANGLE  = POL_ANGLE * numpy.pi / 180.0

            # ! C
            # ! C Compute the actual distance (EPSI_W*) from the orbital focus
            # ! C
            EPSI_WX = EPSI_DX + EPSI_PATH
            EPSI_WZ = EPSI_DZ + EPSI_PATH


            # ! BUG srio@esrf.eu found that these routine does not make the
            # ! calculation correctly. Changed to new one BINORMAL
            # !CALL GAUSS (SIGMAX, EPSI_X, EPSI_WX, XXX, E_BEAM(1), istar1)
            # !CALL GAUSS (SIGMAZ, EPSI_Z, EPSI_WZ, ZZZ, E_BEAM(3), istar1)
            # !
            # ! calculation of the electrom beam moments at the current position
            # ! (sX,sZ) = (epsi_wx,epsi_ez):
            # ! <x2> = sX^2 + sigmaX^2
            # ! <x x'> = sX sigmaXp^2
            # ! <x'2> = sigmaXp^2                 (same for Z)
            #
            # ! then calculate the new recalculated sigmas (rSigmas) and correlation rho of the
            # ! normal bivariate distribution at the point in the electron trajectory
            # ! rsigmaX  = sqrt(<x2>)
            # ! rsigmaXp = sqrt(<x'2>)
            # ! rhoX =  <x x'>/ (rsigmaX rsigmaXp)      (same for Z)
            #
            # if (abs(sigmaX) .lt. 1e-15) then  !no emittance
            #     sigmaXp = 0.0d0
            #     XXX = 0.0
            #     E_BEAM(1) = 0.0
            # else
            #     sigmaXp = epsi_Xold/sigmaX    ! true only at waist, use epsi_xOld as it has been redefined :(
            #     rSigmaX = sqrt( (epsi_wX**2) * (sigmaXp**2) + sigmaX**2 )
            #     rSigmaXp = sigmaXp
            #     if (abs(rSigmaX*rSigmaXp) .lt. 1e-15) then  !no emittance
            #         rhoX = 0.0
            #     else
            #         rhoX = epsi_wx * sigmaXp**2 / (rSigmaX * rSigmaXp)
            #     endif
            #
            #     CALL BINORMAL (rSigmaX, rSigmaXp, rhoX, XXX, E_BEAM(1), istar1)
            # endif
            #

            if self._FLAG_EMITTANCE:
                #     CALL BINORMAL (rSigmaX, rSigmaXp, rhoX, XXX, E_BEAM(1), istar1)
                #     [  c11  c12  ]     [  sigma1^2           rho*sigma1*sigma2   ]
                #     [  c21  c22  ]  =  [  rho*sigma1*sigma2  sigma2^2            ]
                sigmaX,sigmaXp,sigmaZ,sigmaZp = self.syned_electron_beam.get_sigmas_all()

                epsi_wX = sigmaX * sigmaXp
                rSigmaX = numpy.sqrt( (epsi_wX**2) * (sigmaXp**2) + sigmaX**2 )
                rSigmaXp = sigmaXp
                rhoX = epsi_wX * sigmaXp**2 / (rSigmaX * rSigmaXp)
                mean = [0, 0]
                cov = [[sigmaX**2, rhoX*sigmaX*sigmaXp], [rhoX*sigmaX*sigmaXp, sigmaXp**2]]  # diagonal covariance
                sampled_x, sampled_xp = numpy.random.multivariate_normal(mean, cov, 1).T
                # plot_scatter(sampled_x,sampled_xp)
                XXX = sampled_x
                E_BEAM1 = sampled_xp

                epsi_wZ = sigmaZ * sigmaZp
                rSigmaZ = numpy.sqrt( (epsi_wZ**2) * (sigmaZp**2) + sigmaZ**2 )
                rSigmaZp = sigmaZp
                rhoZ = epsi_wZ * sigmaZp**2 / (rSigmaZ * rSigmaZp)
                mean = [0, 0]
                cov = [[sigmaZ**2, rhoZ*sigmaZ*sigmaZp], [rhoZ*sigmaZ*sigmaZp, sigmaZp**2]]  # diagonal covariance
                sampled_z, sampled_zp = numpy.random.multivariate_normal(mean, cov, 1).T
                ZZZ = sampled_z
                E_BEAM3 = sampled_zp

            else:
                sigmaXp = 0.0
                XXX = 0.0
                E_BEAM1 = 0.0
                rhoX = 0.0
                sigmaZp = 0.0
                ZZZ = 0.0
                E_BEAM3 = 0.0

            #
            # ! C
            # ! C For normal wiggler, XXX is perpendicular to the electron trajectory at
            # ! C the point defined by (X_TRAJ,Y_TRAJ,0).
            # ! C
            # IF (F_WIGGLER.EQ.1) THEN   ! normal wiggler
            #     YYY = Y_TRAJ - XXX*SIN(ANGLE)
            #     XXX = X_TRAJ + XXX*COS(ANGLE)

            YYY = Y_TRAJ - XXX * numpy.sin(ANGLE)
            XXX = X_TRAJ + XXX * numpy.cos(ANGLE)

            rays[itik,0] = XXX
            rays[itik,1] = YYY
            rays[itik,2] = ZZZ

            #
            # directions
            #

            #     ! C
            #     ! C Synchrotron source
            #     ! C Note. The angle of emission IN PLANE is the same as the one used
            #     ! C before. This will give rise to a source curved along the orbit.
            #     ! C The elevation angle is instead characteristic of the SR distribution.
            #     ! C The electron beam emittance is included at this stage. Note that if
            #     ! C EPSI = 0, we'll have E_BEAM = 0.0, with no changes.
            #     ! C
            #     IF (F_WIGGLER.EQ.3) ANGLE=0        ! Elliptical Wiggler.
            #     ANGLEX =   ANGLE + E_BEAM(1)
            #     DIREC(1)  =   TAN(ANGLEX)
            #     IF (R_ALADDIN.LT.0.0D0) DIREC(1) = - DIREC(1)
            #     DIREC(2)  =   1.0D0
            #     ARG_ANG  =   GRID(6,ITIK)

            ANGLEX = ANGLE + E_BEAM1
            DIREC1 = numpy.tan(ANGLEX)
            DIREC2 = 1.0


            #     ! C
            #     ! C In the case of SR, we take into account the fact that the electron
            #     ! C trajectory is not orthogonal to the field. This will give a correction
            #     ! C to the photon energy.  We can write it as a correction to the
            #     ! C magnetic field strength; this will linearly shift the critical energy
            #     ! C and, with it, the energy of the emitted photon.
            #     ! C
            #     E_TEMP(3) =   TAN(E_BEAM(3))/COS(E_BEAM(1))
            #     E_TEMP(2) =   1.0D0
            #     E_TEMP(1) =   TAN(E_BEAM(1))
            #     CALL NORM (E_TEMP,E_TEMP)
            #     CORREC =   SQRT(1.0D0-E_TEMP(3)**2)
            #     4400 CONTINUE

            E_TEMP3 = numpy.tan(E_BEAM3)/numpy.cos(E_BEAM1)
            E_TEMP2 = 1.0
            E_TEMP1 = numpy.tan(E_BEAM1)

            e_temp_norm = numpy.sqrt( E_TEMP1**2 + E_TEMP2**2 + E_TEMP3**2)

            E_TEMP3 /= e_temp_norm
            E_TEMP2 /= e_temp_norm
            E_TEMP1 /= e_temp_norm


            CORREC = numpy.sqrt(1.0 - E_TEMP3**2)


            #     IF (FDISTR.EQ.6) THEN
            #         CALL ALADDIN1 (ARG_ANG,ANGLEV,F_POL,IER)
            #         Q_WAVE =   TWOPI*PHOTON(1)/TOCM*CORREC
            #         POL_DEG =   ARG_ANG
            #     ELSE IF (FDISTR.EQ.4) THEN
            #         ARG_ENER =   WRAN (ISTAR1)
            #         RAD_MIN =   ABS(R_MAGNET)
            #
            #         i1 = 1
            #         CALL WHITE  &
            #         (RAD_MIN,CORREC,ARG_ENER,ARG_ANG,Q_WAVE,ANGLEV,POL_DEG,i1)
            #     END IF

            RAD_MIN = numpy.abs(R_MAGNET)


            # CALL WHITE (RAD_MIN,CORREC,ARG_ENER,ARG_ANG,Q_WAVE,ANGLEV,POL_DEG,i1)
            ARG_ANG = numpy.random.random()
            ARG_ENER = numpy.random.random()


            # print("   >> R_MAGNET, DIREC",R_MAGNET,DIREC1,DIREC2)
            # print("   >> RAD_MIN,CORREC,ARG_ENER,ARG_ANG,",RAD_MIN,CORREC,ARG_ENER,ARG_ANG)

#######################################################################
            # gamma = self.syned_electron_beam.gamma()
            # m2ev = codata.c * codata.h / codata.e
            # TOANGS = m2ev * 1e10
            # critical_energy = TOANGS*3.0*numpy.power(gamma,3)/4.0/numpy.pi/1.0e10*(1.0/RAD_MIN)

            # sampled_photon_energy = sampled_energies[itik]
            # wavelength = codata.h * codata.c / codata.e /sampled_photon_energy
            # Q_WAVE = 2 * numpy.pi / (wavelength*1e2)
            # print("   >> PHOTON ENERGY, Ec, lambda, Q: ",sampled_photon_energy,critical_energy,wavelength*1e10,Q_WAVE)
###################################################################################
            sampled_photon_energy = sampled_energies[itik]
            # wavelength = codata.h * codata.c / codata.e /sampled_photon_energy
            critical_energy = TOANGS * 3.0 * numpy.power(gamma, 3) / 4.0 / numpy.pi / 1.0e10 * (1.0 / RAD_MIN)
            eene = sampled_photon_energy / critical_energy

            # TODO: remove old after testing...
            method = "new"

            if method == "old":
                # fm = sync_f(a*1e-3*self.syned_electron_beam.gamma(),eene,polarization=0) * \
                #     numpy.power(eene,2)*a8*self.syned_electron_beam._current*hdiv_mrad * \
                #     numpy.power(self.syned_electron_beam._energy_in_GeV,2)

                fm_s = sync_f(a*1e-3*self.syned_electron_beam.gamma(),eene,polarization=1) * \
                    numpy.power(eene,2)*a8*self.syned_electron_beam._current*hdiv_mrad * \
                    numpy.power(self.syned_electron_beam._energy_in_GeV,2)

                fm_p = sync_f(a*1e-3*self.syned_electron_beam.gamma(),eene,polarization=2) * \
                    numpy.power(eene,2)*a8*self.syned_electron_beam._current*hdiv_mrad * \
                    numpy.power(self.syned_electron_beam._energy_in_GeV,2)
            else:
                fm_s , fm_p = sync_f_sigma_and_pi(a*1e-3*self.syned_electron_beam.gamma(),eene)
                cte = eene ** 2 * a8 * self.syned_electron_beam._current * hdiv_mrad * self.syned_electron_beam._energy_in_GeV ** 2
                fm_s *= cte
                fm_p *= cte

            fm = fm_s + fm_p

            fm_pol = numpy.zeros_like(fm)
            for i in range(fm_pol.size):
                if fm[i] == 0.0:
                    fm_pol[i] = 0
                else:
                    fm_pol[i] = fm_s[i] / fm[i]

            fm.shape = -1
            fm_s.shape = -1
            fm_pol.shape = -1


            pol_deg_interpolator = interp1d(a*1e-3,fm_pol)

            samplerAng = Sampler1D(fm,a*1e-3)

            # samplerPol = Sampler1D(fm_s/fm,a*1e-3)

            # plot(a*1e-3,fm_s/fm)

            if fm.min() == fm.max():
                print("Warning: cannot compute divergence for ray index %d"%itik)
                sampled_theta = 0.0
            else:
                sampled_theta = samplerAng.get_sampled(ARG_ENER)

            sampled_pol_deg = pol_deg_interpolator(sampled_theta)


            # print("sampled_theta: ",sampled_theta, "sampled_energy: ",sampled_photon_energy, "sampled pol ",sampled_pol_deg)

            ANGLEV = sampled_theta
            ANGLEV += E_BEAM3
            #     IF (ANGLEV.LT.0.0) I_CHANGE = -1
            #     ANGLEV =   ANGLEV + E_BEAM(3)
            #     ! C
            #     ! C Test if the ray is within the specified limits
            #     ! C
            #     IF (FGRID.EQ.0.OR.FGRID.EQ.2) THEN
            #         IF (ANGLEV.GT.VDIV1.OR.ANGLEV.LT.-VDIV2) THEN
            #             ARG_ANG = WRAN(ISTAR1)
            #             ! C
            #             ! C If it is outside the range, then generate another ray.
            #             ! C
            #             GO TO 4400
            #         END IF
            #     END IF
            #     DIREC(3)  =   TAN(ANGLEV)/COS(ANGLEX)

            DIREC3 = numpy.tan(ANGLEV) / numpy.cos(ANGLEX)
            #     IF (F_WIGGLER.EQ.3) THEN
            #         CALL ROTATE (DIREC, ANGLE3,ANGLE2,ANGLE1,DIREC)
            #     END IF
            #     CALL NORM (DIREC,DIREC)

            direc_norm = numpy.sqrt(DIREC1**2 + DIREC2**2 + DIREC3**2)

            DIREC1 /= direc_norm
            DIREC2 /= direc_norm
            DIREC3 /= direc_norm

            rays[itik,3] = DIREC1 # VX
            rays[itik,4] = DIREC2 # VY
            rays[itik,5] = DIREC3 # VZ

        if user_unit_to_m != 1.0:
            rays[:,0] /= user_unit_to_m
            rays[:,1] /= user_unit_to_m
            rays[:,2] /= user_unit_to_m

        #
        # sample divergences (cols 4-6): the Shadow way
        #




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



        POL_DEG = sampled_pol_deg
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

        #
        POL_ANGLE = 0.5 * numpy.pi

        if F_COHER == 1:
            PHASEX = 0.0
        else:
            PHASEX = numpy.random.random(NRAYS) * 2 * numpy.pi

        # PHASEZ = PHASEX + POL_ANGLE * numpy.sign(ANGLEV)

        rays[:,13] = 0.0 # PHASEX
        rays[:,14] = 0.0 # PHASEZ

        # set flag (col 10)
        rays[:,9] = 1.0

        #
        # photon energy (col 11)
        #

        # A2EV = 2.0*numpy.pi/(codata.h*codata.c/codata.e*1e2)
        sampled_photon_energy = sampled_energies
        wavelength = codata.h * codata.c / codata.e /sampled_photon_energy
        Q_WAVE = 2 * numpy.pi / (wavelength*1e2)
        rays[:,10] =  Q_WAVE # sampled_photon_energy * A2EV

        # col 12 (ray index)
        rays[:,11] =  1 + numpy.arange(NRAYS)

        # col 13 (optical path)
        rays[:,11] = 0.0

        return rays

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
        return 0,0,0


if __name__ == "__main__":


    e_min = 5000.0 # 70490.0 #
    e_max = 100000.0 # 70510.0 #
    e_min = 70490.0 #
    e_max = 70510.0 #
    NRAYS = 5000
    use_emittances=True



    wigFile = "xshwig.sha"
    inData = ""

    nPer = 5 # 50
    nTrajPoints = 501
    ener_gev = 6.04
    per = 0.040
    kValue = 7.85
    trajFile = "tmp.traj"
    shift_x_flag = 0
    shift_x_value = 0.0
    shift_betax_flag = 0
    shift_betax_value = 0.0


    sw = SourceWiggler()

    #
    # syned
    #
    syned_wiggler = Wiggler(K_vertical=kValue,K_horizontal=0.0,period_length=per,number_of_periods=nPer)

    syned_electron_beam = ElectronBeam(energy_in_GeV=6.04,
                 energy_spread = 0.0,
                 current = 0.2,
                 number_of_bunches = 400,
                 moment_xx=(400e-6)**2,
                 moment_xxp=0.0,
                 moment_xpxp=(10e-6)**2,
                 moment_yy=(10e-6)**2,
                 moment_yyp=0.0,
                 moment_ypyp=(4e-6)**2 )

    sourcewiggler = SourceWiggler(name="test",syned_electron_beam=syned_electron_beam,
                    syned_wiggler=syned_wiggler,
                    flag_emittance=use_emittances,
                    emin=e_min,emax=e_max,ng_e=10, ng_j=nTrajPoints)



    print(sourcewiggler.info())


    rays = sourcewiggler.calculate_rays(NRAYS=NRAYS)

    plot_scatter(rays[:,1],rays[:,0],title="trajectory",show=False)
    plot_scatter(rays[:,0],rays[:,2],title="real space",show=False)
    plot_scatter(rays[:,3],rays[:,5],title="divergence space")

