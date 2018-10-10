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

from srxraylib.util.inverse_method_sampler import Sampler1D, Sampler2D, Sampler3D
from srxraylib.plot.gol import plot,plot_scatter
import scipy.constants as codata
from scipy import interpolate

from syned.storage_ring.magnetic_structures.wiggler import Wiggler
from syned.storage_ring.electron_beam import ElectronBeam

from srfunc import wiggler_trajectory, wiggler_cdf

from scipy.interpolate import interp1d


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

        # ray tracing
        # self.SEED            = SEED   # Random seed
        # self.NRAYS           = NRAYS  # Number of rays


        self._FLAG_EMITTANCE  =  flag_emittance # Yes  # Use emittance (0=No, 1=Yes)

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

    def get_energy_box(self):
        """
        Gets the limits of photon energy distribution for the source
        :return: emin,emax,number_of_points
        """
        return self._EMIN,self._EMAX,self._NG_E


    def calculate_radiation(self):

        (traj, pars) = wiggler_trajectory(b_from=0,
                                                 inData="",
                                                 nPer=self.syned_wiggler.number_of_periods(),
                                                 nTrajPoints=self._NG_J,
                                                 ener_gev=self.syned_electron_beam._energy_in_GeV,
                                                 per=self.syned_wiggler.period_length(),
                                                 kValue=self.syned_wiggler.K_vertical(),
                                                 trajFile="",
                                                 shift_x_flag=0.0,
                                                 shift_x_value=0.0,
                                                 shift_betax_flag=0,
                                                 shift_betax_value=0.0)

        self._result_trajectory = traj
        self._result_parameters = pars


        # traj[0,ii] = yx[i]
        # traj[1,ii] = yy[i]+j * per - start_len
        # traj[2,ii] = 0.0
        # traj[3,ii] = betax[i]
        # traj[4,ii] = betay[i]
        # traj[5,ii] = 0.0
        # traj[6,ii] = curv[i]
        # traj[7,ii] = bz[i]

        # plot(traj[1,:],traj[0,:])
        # print(pars)

        #
        # calculate cdf and write file for Shadow/Source
        #

        self._result_cdf = wiggler_cdf(self._result_trajectory,
                           enerMin=self._EMIN,
                           enerMax=self._EMAX,
                           enerPoints=self._NG_E,
                           outFile="",
                           elliptical=False)


    def calculate_rays(self,user_unit_to_m=1.0,F_COHER=0,NRAYS=5000,SEED=123456):
        """
        compute the rays in SHADOW matrix (shape (npoints,18) )
        :param F_COHER: set this flag for coherent beam
        :param user_unit_to_m: default 1.0 (m)
        :return: rays, a numpy.array((npoits,18))
        """

        if self._result_cdf is None:
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

        x_photon = 0.0
        y_photon = 0.0 # numpy.random.random(NRAYS)
        z_photon = 0.0

        # traj[0,ii] = yx[i]
        # traj[1,ii] = yy[i]+j * per - start_len
        # traj[2,ii] = 0.0
        # traj[3,ii] = betax[i]
        # traj[4,ii] = betay[i]
        # traj[5,ii] = 0.0
        # traj[6,ii] = curv[i]
        # traj[7,ii] = bz[i]

        PATH_STEP = self._result_cdf["step"]
        print(">>>>>>",PATH_STEP)
        X_TRAJ = self._result_trajectory[0,:]
        Y_TRAJ = self._result_trajectory[1,:]
        ANGLE  = self._result_trajectory[3,:]
        CURV   = self._result_trajectory[6,:]
        EPSI_PATH = numpy.arange(CURV.size) * PATH_STEP # self._result_trajectory[7,:]

        # plot(Y_TRAJ,X_TRAJ)


        # ! C We define the 5 arrays:
        # ! C    Y_X(5,N)    ---> X(Y)
        # ! C    Y_XPRI(5,N) ---> X'(Y)
        # ! C    Y_CURV(5,N) ---> CURV(Y)
        # ! C    Y_PATH(5,N) ---> PATH(Y)
        # ! C    F(1,N) contains the array of Y values where the nodes are located.

        # print(">>>>>>>",self._result_trajectory.shape)
        # plot(self._result_trajectory[1,:],self._result_trajectory[0,:])

        # CALL PIECESPL(SEED_Y, Y_TEMP,   NP_SY,   IER)
        # CALL CUBSPL (Y_X,    X_TEMP,   NP_TRAJ, IER)
        # CALL CUBSPL (Y_Z,    Z_TEMP,   NP_TRAJ, IER)
        # CALL CUBSPL (Y_XPRI, ANG_TEMP, NP_TRAJ, IER)
        # CALL CUBSPL (Y_ZPRI, ANG2_TEMP, NP_TRAJ, IER)
        # CALL CUBSPL (Y_CURV, C_TEMP,   NP_TRAJ, IER)
        # CALL CUBSPL (Y_PATH, P_TEMP,   NP_TRAJ, IER)

        SEED_Y = interp1d(numpy.arange(Y_TRAJ.size)/(Y_TRAJ.size-1),Y_TRAJ,kind='linear')
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
        F_PHOT  = 0
        F_COLOR  = 3
        FSOUR  = 3
        FDISTR  = 4

        ARG_Y = numpy.random.random(NRAYS)

        for itik in range(NRAYS):


            # IF (F_WIGGLER.EQ.1) THEN
            #     ARG_Y = GRID(2,ITIK)
            #     CALL SPL_INT (SEED_Y, NP_SY,   ARG_Y,  Y_TRAJ,    IER)
            # arg_y = ARG_Y[itik]
            arg_y = itik/(NRAYS-1)
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

            print("\n>>><<<",arg_y,Y_TRAJ,X_TRAJ,ANGLE,CURV,EPSI_PATH)


            # ! C
            # ! C Compute the actual distance (EPSI_W*) from the orbital focus
            # ! C
            # EPSI_WX = EPSI_DX + EPSI_PATH
            # EPSI_WZ = EPSI_DZ + EPSI_PATH
            #

            # EPSI_DX = 0.0
            # EPSI_DZ = 0.0
            #
            # EPSI_WX = EPSI_DX + EPSI_PATH
            # EPSI_WZ = EPSI_DZ + EPSI_PATH


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
                #todo
                pass
            else:
                sigmaXp = 0.0
                XXX = 0.0
                E_BEAM1 = 0.0

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

            # print(">>>>> XXX,YYY,ZZZ,ANGLE",XXX,YYY,ZZZ,ANGLE)

            # plot_scatter(YYY,XXX,title="Number of rays: %d"%XXX.size)

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
            #     IF (F_WIGGLER.EQ.3) THEN
            #         CALL ROTATE (DIREC, ANGLE3,ANGLE2,ANGLE1,DIREC)
            #     END IF
            #     CALL NORM (DIREC,DIREC)
            # print*,">>>>DIREC,FGRID,R_ALADDIN: ",DIREC,FGRID,R_ALADDIN
            #     GO TO 1111


            rays[itik,3] = 0.0 # VX
            rays[itik,4] = 1.0 # VY
            rays[itik,5] = 0.0 # VZ

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




        # DENOM = numpy.sqrt(1.0 - 2.0 * POL_DEG + 2.0 * POL_DEG**2)
        # AX = POL_DEG/DENOM
        # for i in range(3):
        #     A_VEC[:,i] *= AX
        #
        # AZ = (1.0-POL_DEG)/DENOM
        # for i in range(3):
        #     AP_VEC[:,i] *= AZ

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

        A2EV = 2.0*numpy.pi/(codata.h*codata.c/codata.e*1e2)
        rays[:,10] =  sampled_photon_energy * A2EV

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


