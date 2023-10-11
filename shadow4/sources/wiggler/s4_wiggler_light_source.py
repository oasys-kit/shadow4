__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "30-08-2018"

"""

Wiggler code: computes wiggler radiation distributions and samples rays according to them.

Fully replaces and upgrades the shadow3 wiggler model.

The radiation is calculating using sr-xraylib

"""

import numpy

from srxraylib.util.inverse_method_sampler import Sampler1D
from srxraylib.sources.srfunc import wiggler_trajectory, wiggler_spectrum, wiggler_cdf, sync_f

import scipy
from scipy.interpolate import interp1d
import scipy.constants as codata

from shadow4.sources.s4_electron_beam import S4ElectronBeam
from shadow4.sources.s4_light_source import S4LightSource
from shadow4.sources.wiggler.s4_wiggler import S4Wiggler
from shadow4.beam.s4_beam import S4Beam

# This is similar to sync_f in srxraylib but faster
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


class S4WigglerLightSource(S4LightSource):

    def __init__(self,
                 name="Undefined",
                 electron_beam=None,
                 magnetic_structure=None,
                 nrays=5000,
                 seed=12345,
                 ):
        super().__init__(name,
                         electron_beam=electron_beam if not electron_beam is None else S4ElectronBeam(),
                         magnetic_structure=magnetic_structure if not magnetic_structure is None else S4Wiggler(),
                         nrays=nrays,
                         seed=seed)

        # results of calculations
        self.__result_trajectory = None
        self.__result_parameters = None
        self.__result_cdf = None


    def get_trajectory(self):
        return self.__result_trajectory, self.__result_parameters


    def __calculate_radiation(self):

        wiggler = self.get_magnetic_structure()
        electron_beam = self.get_electron_beam()

        if wiggler._magnetic_field_periodic == 1:

            (traj, pars) = wiggler_trajectory(b_from=0,
                                                     inData="",
                                                     nPer=wiggler.number_of_periods(),
                                                     nTrajPoints=wiggler._NG_J,
                                                     ener_gev=electron_beam._energy_in_GeV,
                                                     per=wiggler.period_length(),
                                                     kValue=wiggler.K_vertical(),
                                                     trajFile="",
                                                     shift_x_flag=wiggler._shift_x_flag,
                                                     shift_x_value=wiggler._shift_x_value,
                                                     shift_betax_flag=wiggler._shift_betax_flag,
                                                     shift_betax_value=wiggler._shift_betax_value, )

        elif wiggler._magnetic_field_periodic == 0:

            (traj, pars) = wiggler_trajectory(b_from=1,
                                                     inData=wiggler._file_with_magnetic_field,
                                                     nPer=1,
                                                     nTrajPoints=wiggler._NG_J,
                                                     ener_gev=electron_beam._energy_in_GeV,
                                                     # per=self.syned_wiggler.period_length(),
                                                     # kValue=self.syned_wiggler.K_vertical(),
                                                     trajFile="",
                                                     shift_x_flag       = wiggler._shift_x_flag     ,
                                                     shift_x_value      = wiggler._shift_x_value    ,
                                                     shift_betax_flag   = wiggler._shift_betax_flag ,
                                                     shift_betax_value  = wiggler._shift_betax_value,)


        self.__result_trajectory = traj
        self.__result_parameters = pars

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

        self.__result_cdf = wiggler_cdf(self.__result_trajectory,
                                        enerMin=wiggler._EMIN,
                                        enerMax=wiggler._EMAX,
                                        enerPoints=wiggler._NG_E,
                                        outFile="tmp.cdf",
                                        elliptical=False)


    def __calculate_rays(self,user_unit_to_m=1.0,F_COHER=0,EPSI_DX=0.0,EPSI_DZ=0.0,
                       psi_interval_in_units_one_over_gamma=None,
                       psi_interval_number_of_points=1001,
                       verbose=True):
        """
        compute the rays in SHADOW matrix (shape (npoints,18) )
        :param F_COHER: set this flag for coherent beam
        :param user_unit_to_m: default 1.0 (m)
        :return: rays, a numpy.array((npoits,18))
        """




        if self.__result_cdf is None:
            self.__calculate_radiation()

        if verbose:
            print(">>>   Results of calculate_radiation")
            print(">>>       trajectory.shape: ", self.__result_trajectory.shape)
            print(">>>       cdf: ", self.__result_cdf.keys())

        wiggler = self.get_magnetic_structure()
        syned_electron_beam = self.get_electron_beam()


        sampled_photon_energy,sampled_theta,sampled_phi = self._sample_photon_energy_theta_and_phi()

        if verbose:
            print(">>> sampled sampled_photon_energy,sampled_theta,sampled_phi:  ",sampled_photon_energy,sampled_theta,sampled_phi)

        if self.get_seed() != 0:
            numpy.random.seed(self.get_seed())


        sigmas = syned_electron_beam.get_sigmas_all()

        NRAYS = self.get_nrays()

        rays = numpy.zeros((NRAYS,18))

        #
        # sample sizes (cols 1-3)
        #

        #
        if wiggler._FLAG_EMITTANCE:
            if numpy.array(numpy.abs(sigmas)).sum() == 0:
                wiggler._FLAG_EMITTANCE = False

        if wiggler._FLAG_EMITTANCE:
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

        PATH_STEP = self.__result_cdf["step"]
        X_TRAJ = self.__result_cdf["x"]
        Y_TRAJ = self.__result_cdf["y"]
        SEEDIN = self.__result_cdf["cdf"]
        ANGLE  = self.__result_cdf["angle"]
        CURV   = self.__result_cdf["curv"]
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

        ws_ev,ws_f,tmp =  wiggler_spectrum(self.__result_trajectory,
                                           enerMin=wiggler._EMIN, enerMax=wiggler._EMAX, nPoints=500,
                                           # per=self.syned_wiggler.period_length(),
                                           electronCurrent=syned_electron_beam._current,
                                           outFile="", elliptical=False)

        ws_flux_per_ev = ws_f / (ws_ev*1e-3)
        samplerE = Sampler1D(ws_flux_per_ev,ws_ev)

        sampled_energies,h,h_center = samplerE.get_n_sampled_points_and_histogram(NRAYS)


        ###############################################

        gamma = syned_electron_beam.gamma()
        m2ev = codata.c * codata.h / codata.e
        TOANGS = m2ev * 1e10


        #####################################################

        RAD_MIN = 1.0 / numpy.abs(self.__result_cdf["curv"]).max()

        critical_energy = TOANGS * 3.0 * numpy.power(gamma, 3) / 4.0 / numpy.pi / 1.0e10 * (1.0 / RAD_MIN)

        if psi_interval_in_units_one_over_gamma is None:
            c = numpy.array([-0.3600382, 0.11188709])  # see file fit_psi_interval.py
            # x = numpy.log10(self._EMIN / critical_energy)
            x = numpy.log10(wiggler._EMIN / (4 * critical_energy)) # the wiggler that does not have an unique
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

            if wiggler._FLAG_EMITTANCE:
                #     CALL BINORMAL (rSigmaX, rSigmaXp, rhoX, XXX, E_BEAM(1), istar1)
                #     [  c11  c12  ]     [  sigma1^2           rho*sigma1*sigma2   ]
                #     [  c21  c22  ]  =  [  rho*sigma1*sigma2  sigma2^2            ]
                sigmaX,sigmaXp,sigmaZ,sigmaZp = syned_electron_beam.get_sigmas_all()

                rSigmaX = numpy.sqrt( (EPSI_WX**2) * (sigmaXp**2) + sigmaX**2 )
                rSigmaXp = sigmaXp
                rhoX = EPSI_WX * sigmaXp**2 / (rSigmaX * rSigmaXp)
                mean = [0, 0]
                cov = [[sigmaX**2, rhoX*sigmaX*sigmaXp], [rhoX*sigmaX*sigmaXp, sigmaXp**2]]  # diagonal covariance
                sampled_x, sampled_xp = numpy.random.multivariate_normal(mean, cov, 1).T
                XXX = sampled_x
                E_BEAM1 = sampled_xp

                rSigmaZ = numpy.sqrt( (EPSI_WZ**2) * (sigmaZp**2) + sigmaZ**2 )
                rSigmaZp = sigmaZp
                rhoZ = EPSI_WZ * sigmaZp**2 / (rSigmaZ * rSigmaZp)
                mean = [0, 0]
                cov = [[sigmaZ**2, rhoZ*sigmaZ*sigmaZp], [rhoZ*sigmaZ*sigmaZp, sigmaZp**2]]  # diagonal covariance
                sampled_z, sampled_zp = numpy.random.multivariate_normal(mean, cov, 1).T
                ZZZ = sampled_z
                E_BEAM3 = sampled_zp

            else:
                XXX = 0.0
                E_BEAM1 = 0.0
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
                fm_s , fm_p = sync_f_sigma_and_pi(a*1e-3*syned_electron_beam.gamma(),eene)
                cte = eene ** 2 * a8 * syned_electron_beam._current * hdiv_mrad * syned_electron_beam._energy_in_GeV ** 2
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
        rays[:,12] = 0.0

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

    def _sample_photon_energy_theta_and_phi(self):

        #
        # sample divergences
        #
        return 0,0,0

    ############################################################################
    #
    ############################################################################
    def get_beam(self):

        user_unit_to_m = 1.0
        F_COHER = 0
        EPSI_DX = self.get_magnetic_structure()._EPSI_DX
        EPSI_DZ = self.get_magnetic_structure()._EPSI_DZ
        psi_interval_in_units_one_over_gamma = None
        psi_interval_number_of_points = 1001
        verbose = True

        return S4Beam.initialize_from_array(self.__calculate_rays(
            user_unit_to_m=user_unit_to_m,
            F_COHER=F_COHER,
            EPSI_DX=EPSI_DX,
            EPSI_DZ=EPSI_DZ,
            psi_interval_in_units_one_over_gamma=psi_interval_in_units_one_over_gamma,
            psi_interval_number_of_points=psi_interval_number_of_points,
            verbose=verbose))

    def calculate_spectrum(self, output_file=""):
        traj, pars = self.get_trajectory()
        wig = self.get_magnetic_structure()
        e_min, e_max, ne = wig.get_energy_box()
        ring = self.get_electron_beam()
        if traj is not None:
            e, f, w = wiggler_spectrum(traj,
                          enerMin=e_min,
                          enerMax=e_max,
                          nPoints=ne,
                          electronCurrent=ring.current(),
                          outFile=output_file,
                          elliptical=False)
            return e,f,w
        else:
            raise Exception("Cannot compute spectrum")

    def info(self,debug=False):
        electron_beam = self.get_electron_beam()
        magnetic_structure = self.get_magnetic_structure()

        txt = ""
        txt += "-----------------------------------------------------\n"

        txt += "Input Electron parameters: \n"
        txt += "        Electron energy: %f geV\n"%electron_beam.energy()
        txt += "        Electron current: %f A\n"%electron_beam.current()
        if magnetic_structure._FLAG_EMITTANCE:
            sigmas = electron_beam.get_sigmas_all()
            txt += "        Electron sigmaX: %g [um]\n"%(1e6*sigmas[0])
            txt += "        Electron sigmaZ: %g [um]\n"%(1e6*sigmas[2])
            txt += "        Electron sigmaX': %f urad\n"%(1e6*sigmas[1])
            txt += "        Electron sigmaZ': %f urad\n"%(1e6*sigmas[3])

        txt += "Lorentz factor (gamma): %f\n"%electron_beam.gamma()

        txt2 = magnetic_structure.info()
        return (txt + "\n\n" + txt2)

    def to_python_code(self, **kwargs):
        script = ''
        try:
            script += self.get_electron_beam().to_python_code()
        except:
            script += "\n\n#Error retrieving electron_beam code"

        try:
            script += self.get_magnetic_structure().to_python_code()
        except:
            script += "\n\n#Error retrieving magnetic structure code"

        script += "\n\n\n#light source\nfrom shadow4.sources.wiggler.s4_wiggler_light_source import S4WigglerLightSource"
        script += "\nlight_source = S4WigglerLightSource(name='%s',electron_beam=electron_beam,magnetic_structure=source,nrays=%d,seed=%s)" % \
                                                          (self.get_name(),self.get_nrays(),self.get_seed())
        #
        # script += "\n\n\n#beamline\nfrom shadow4.beamline.s4_beamline import S4Beamline"
        # script += "\nbeamline = S4Beamline(light_source=light_source)"
        script += "\nbeam = light_source.get_beam()"
        return script
if __name__ == "__main__":
    from srxraylib.plot.gol import plot_scatter, set_qt
    set_qt()

    # e_min = 5000.0 # 70490.0 #
    # e_max = 100000.0 # 70510.0 #
    # e_min = 70490.0 #
    # e_max = 70510.0 #
    # NRAYS = 5000
    # use_emittances=True
    #
    #
    #
    # wigFile = "xshwig.sha"
    # inData = ""
    #
    # nPer = 5 # 50
    # nTrajPoints = 501
    # ener_gev = 6.04
    # per = 0.040
    # kValue = 7.85
    # trajFile = "tmp.traj"
    # shift_x_flag = 0
    # shift_x_value = 0.0
    # shift_betax_flag = 0
    # shift_betax_value = 0.0
    #
    #
    #
    #
    # #
    # # syned
    # #
    #
    # electron_beam = S4ElectronBeam(energy_in_GeV=6.04,
    #                                energy_spread = 0.0,
    #                                current = 0.2,
    #                                number_of_bunches = 400,
    #                                moment_xx=(400e-6)**2,
    #                                moment_xxp=0.0,
    #                                moment_xpxp=(10e-6)**2,
    #                                moment_yy=(10e-6)**2,
    #                                moment_yyp=0.0,
    #                                moment_ypyp=(4e-6)**2)
    #
    #
    # w = S4Wiggler(K_vertical=kValue,period_length=per,number_of_periods=nPer,
    #                             flag_emittance=use_emittances,
    #                             emin=e_min, emax=e_max,ng_e=10, ng_j=nTrajPoints)
    #
    #
    #
    # # print(w.info())
    #
    # ls = S4WigglerLightSource(name="Undefined", electron_beam=electron_beam, magnetic_structure=w,
    #                           nrays=NRAYS)
    #
    # print(ls.info())
    #
    #
    # beam = ls.get_beam()
    #
    # rays = beam.rays
    #
    # plot_scatter(rays[:,1],rays[:,0],title="trajectory",show=False)
    # plot_scatter(rays[:,0],rays[:,2],title="real space",show=False)
    # plot_scatter(rays[:,3],rays[:,5],title="divergence space")

    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam

    electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
    electron_beam.set_sigmas_all(sigma_x=2.377e-05, sigma_y=2.472e-05, sigma_xp=3.58e-06, sigma_yp=3.04e-06)

    # magnetic structure
    from shadow4.sources.wiggler.s4_wiggler import S4Wiggler

    source = S4Wiggler(
        magnetic_field_periodic=0,  # 0=external, 1=periodic
        file_with_magnetic_field="/nobackup/gurb1/srio/Oasys/SW_BM18_Joel.txt",  # used only if magnetic_field_periodic=0
        K_vertical=10.0,  # syned Wiggler pars: used only if magnetic_field_periodic=1
        period_length=0.1,  # syned Wiggler pars: used only if magnetic_field_periodic=1
        number_of_periods=10,  # syned Wiggler pars: used only if magnetic_field_periodic=1
        emin=100.0,  # Photon energy scan from energy (in eV)
        emax=200000.0,  # Photon energy scan to energy (in eV)
        ng_e=101,  # Photon energy scan number of points
        ng_j=501,  # Number of points in electron trajectory (per period) for internal calculation only
        flag_emittance=1,  # Use emittance (0=No, 1=Yes)
        shift_x_flag=4,  # 0="No shift", 1="Half excursion", 2="Minimum", 3="Maximum", 4="Value at zero", 5="User value"
        shift_x_value=0.001,  # used only if shift_x_flag=5
        shift_betax_flag=4,  # 0="No shift", 1="Half excursion", 2="Minimum", 3="Maximum", 4="Value at zero", 5="User value"
        shift_betax_value=0.0,  # used only if shift_betax_flag=5
    )

    # light source
    from shadow4.sources.wiggler.s4_wiggler_light_source import S4WigglerLightSource

    light_source = S4WigglerLightSource(name='wiggler', electron_beam=electron_beam, magnetic_structure=source, nrays=2000,
                                        seed=5676561)
    beam = light_source.get_beam()

    # test plot
    from srxraylib.plot.gol import plot_scatter

    rays = beam.get_rays()
    plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')