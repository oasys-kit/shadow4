__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "30-08-2018"

"""


"""


import numpy

from srxraylib.util.inverse_method_sampler import Sampler1D, Sampler2D, Sampler3D
from srxraylib.plot.gol import plot,plot_scatter
import scipy.constants as codata
from scipy import interpolate

from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet
from syned.storage_ring.electron_beam import ElectronBeam

# from minishadow.bending_magnet.srfunc import wiggler_trajectory, wiggler_spectrum, wiggler_cdf, sync_g1, sync_f, sync_ene, sync_ang
#
# from scipy.interpolate import interp1d


class SourceBendingMagnet(object):
    def __init__(self,name="",
                 syned_electron_beam=None,
                 syned_bending_magnet=None,
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
        if syned_bending_magnet is None:
            self.syned_bending_magnet = BendingMagnet(radius,magnetic_field,length)
        else:
            self.syned_bending_magnet = syned_bending_magnet

        # Photon energy scan
        self._EMIN            = emin   # Photon energy scan from energy (in eV)
        self._EMAX            = emax   # Photon energy scan to energy (in eV)
        self._NG_E            = ng_e   # Photon energy scan number of points

        self._NG_J            = ng_j       # Number of points in electron trajectory (per period)


        self._FLAG_EMITTANCE  =  flag_emittance # Yes  # Use emittance (0=No, 1=Yes) #todo kw in calculate rays

        # # results of calculations
        #
        # self._result_trajectory = None
        # self._result_parameters = None
        # self._result_cdf = None


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
        txt += "        Electron energy: %f geV\n"%self.syned_electron_beam.energy()
        txt += "        Electron current: %f A\n"%self.syned_electron_beam.current()
        if self._FLAG_EMITTANCE:
            sigmas = self.syned_electron_beam.get_sigmas_all()
            txt += "        Electron sigmaX: %g [um]\n"%(1e6*sigmas[0])
            txt += "        Electron sigmaZ: %g [um]\n"%(1e6*sigmas[2])
            txt += "        Electron sigmaX': %f urad\n"%(1e6*sigmas[1])
            txt += "        Electron sigmaZ': %f urad\n"%(1e6*sigmas[3])
        txt += "Input Bending Magnet parameters: \n"
        txt += "        radius: %f m\n"%self.syned_bending_magnet._radius
        txt += "        magnetic field: %f T\n"%self.syned_bending_magnet._magnetic_field
        txt += "        length: %f m\n"%self.syned_bending_magnet._length


        txt += "-----------------------------------------------------\n"

        txt += "Lorentz factor (gamma): %f\n"%self.syned_electron_beam.gamma()
        txt += "Electron velocity: %.12f c units\n"%(numpy.sqrt(1.0 - 1.0 / self.syned_electron_beam.gamma() ** 2))

        #
        txt += "-----------------------------------------------------\n"
        txt += "Grids: \n"
        if self._NG_E == 1:
            txt += "        photon energy %f eV\n"%(self._EMIN)
        else:
            txt += "        photon energy from %10.3f eV to %10.3f eV\n"%(self._EMIN,self._EMAX)
        txt += "        number of energy points: %d\n"%(self._NG_E)
        txt += "        number of points for the trajectory: %d\n"%(self._NG_J)


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

    #
    # def calculate_radiation(self):
    #
    #     (traj, pars) = wiggler_trajectory(b_from=0,
    #                                              inData="",
    #                                              nPer=self.syned_wiggler.number_of_periods(),
    #                                              nTrajPoints=self._NG_J,
    #                                              ener_gev=self.syned_electron_beam._energy_in_GeV,
    #                                              per=self.syned_wiggler.period_length(),
    #                                              kValue=self.syned_wiggler.K_vertical(),
    #                                              trajFile="",
    #                                              shift_x_flag=0.0,
    #                                              shift_x_value=0.0,
    #                                              shift_betax_flag=0,
    #                                              shift_betax_value=0.0)
    #
    #     self._result_trajectory = traj
    #     self._result_parameters = pars
    #
    #
    #     # traj[0,ii] = yx[i]
    #     # traj[1,ii] = yy[i]+j * per - start_len
    #     # traj[2,ii] = 0.0
    #     # traj[3,ii] = betax[i]
    #     # traj[4,ii] = betay[i]
    #     # traj[5,ii] = 0.0
    #     # traj[6,ii] = curv[i]
    #     # traj[7,ii] = bz[i]
    #
    #
    #     #
    #     # calculate cdf and write file for Shadow/Source
    #     #
    #
    #     self._result_cdf = wiggler_cdf(self._result_trajectory,
    #                        enerMin=self._EMIN,
    #                        enerMax=self._EMAX,
    #                        enerPoints=self._NG_E,
    #                        outFile="",
    #                        elliptical=False)
    #
    #
    def calculate_rays(self,user_unit_to_m=1.0,F_COHER=0,NRAYS=5000,SEED=123456,
                       EPSI_DX=0.0,EPSI_DZ=0.0):
        """
        compute the rays in SHADOW matrix (shape (npoints,18) )
        :param F_COHER: set this flag for coherent beam
        :param user_unit_to_m: default 1.0 (m)
        :return: rays, a numpy.array((npoits,18))
        """
    #
    #     if self._result_cdf is None:
    #         self.calculate_radiation()
    #
    #     sampled_photon_energy,sampled_theta,sampled_phi = self._sample_photon_energy_theta_and_phi(NRAYS)
    #
        if SEED != 0:
            numpy.random.seed(SEED)


        sigmas = self.syned_electron_beam.get_sigmas_all()

        rays = numpy.zeros((NRAYS,18))

        RAD_MIN= numpy.abs(self.syned_bending_magnet._radius)
        RAD_MAX= numpy.abs(self.syned_bending_magnet._radius)

        # r_aladdin	=  bending magnet radius in units of length used for source size, CCW rings negative.

        r_aladdin = self.syned_bending_magnet._radius



        F_COHER = 0
        if r_aladdin < 0:
            POL_ANGLE = -90.0 * numpy.pi / 2
        else:
            POL_ANGLE = 90.0 * numpy.pi / 2

        # hdiv1	=  0.0000000000000000E+00 - horizontal divergence in +X (radians).
        # hdiv2	=  0.0000000000000000E+00 - horizontal divergence in -X (radians).

        HDIV1 = 0.5 * numpy.abs( self.syned_bending_magnet._length / self.syned_bending_magnet._radius )
        HDIV2 = HDIV1


        for itik in range(NRAYS):
            # ! Synchrontron depth
            ANGLE  =  numpy.random.random() * (HDIV1 + HDIV2) - HDIV2
            EPSI_PATH =  numpy.abs(r_aladdin) * ANGLE



            #
            # ! then calculate the new recalculated sigmas (rSigmas) and correlation rho of the
            # ! normal bivariate distribution at the point in the electron trajectory
            # ! rsigmaX  = sqrt(<x2>)
            # ! rsigmaXp = sqrt(<x'2>)
            # ! rhoX =  <x x'>/ (rsigmaX rsigmaXp)      (same for Z)
            #
            # if (abs(sigmaX) .lt. 1e-15) then  !no emittance
            # sigmaXp = 0.0d0
            # XXX = 0.0
            # E_BEAM(1) = 0.0
            # else
            # sigmaXp = epsi_Xold/sigmaX    ! true only at waist, use epsi_xOld as it has been redefined :(
            # rSigmaX = sqrt( (epsi_wX**2) * (sigmaXp**2) + sigmaX**2 )
            # rSigmaXp = sigmaXp
            # if (abs(rSigmaX*rSigmaXp) .lt. 1e-15) then  !no emittance
            # rhoX = 0.0
            # else
            # rhoX = epsi_wx * sigmaXp**2 / (rSigmaX * rSigmaXp)
            # endif


            if self._FLAG_EMITTANCE:
                sigma_x, sigma_xp, sigma_z, sigma_zp = self.syned_electron_beam.get_sigmas_all()

                # ! calculation of the electrom beam moments at the current position
                # ! (sX,sZ) = (epsi_wx,epsi_ez):
                # ! <x2> = sX^2 + sigmaX^2
                # ! <x x'> = sX sigmaXp^2
                # ! <x'2> = sigmaXp^2                 (same for Z)

                epsi_wX = sigma_x * sigma_xp
                rSigmaX = numpy.sqrt( (epsi_wX**2) * (sigma_xp**2) + sigma_x**2 )
                rSigmaXp = sigma_xp
                rhoX = epsi_wX * sigma_xp**2 / (rSigmaX * rSigmaXp)
                mean = [0, 0]
                cov = [[rSigmaX**2, rhoX*rSigmaX*rSigmaXp], [rhoX*rSigmaX*rSigmaXp, rSigmaXp**2]]  # diagonal covariance
                sampled_x, sampled_xp = numpy.random.multivariate_normal(mean, cov, 1).T
                # plot_scatter(sampled_x,sampled_xp,title="X")
                XXX = sampled_x
                E_BEAM1 = sampled_xp


                epsi_wZ = sigma_z * sigma_zp
                rSigmaZ = numpy.sqrt( (epsi_wZ**2) * (sigma_zp**2) + sigma_z**2 )
                rSigmaZp = sigma_zp
                rhoZ = epsi_wZ * sigma_zp**2 / (rSigmaZ * rSigmaZp)
                mean = [0, 0]
                cov = [[rSigmaZ**2, rhoZ*rSigmaZ*rSigmaZp], [rhoZ*rSigmaZ*rSigmaZp, rSigmaZp**2]]  # diagonal covariance
                sampled_z, sampled_zp = numpy.random.multivariate_normal(mean, cov, 1).T
                # plot_scatter(sampled_z,sampled_zp,title="Z")
                ZZZ = sampled_z
                E_BEAM3 = sampled_zp

                # print(">>>>>>>>>",sampled_x,sampled_z)
            else:
                sigma_x, sigma_xp, sigma_z, sigma_zp = (0.0, 0.0, 0.0, 0.0)

                rhoX = 0.0

                XXX = 0.0
                E_BEAM1 = 0.0

                ZZZ = 0.0
                E_BEAM3 = 0.0


            # ! C
            # ! C Synchrotron depth distribution
            # ! C
            # 440	CONTINUE
            # ! CC	R_ALADDIN NEGATIVE FOR COUNTER-CLOCKWISE SOURCE
            # IF (R_ALADDIN.LT.0) THEN
            # YYY = (ABS(R_ALADDIN) + XXX) * SIN(ANGLE)
            # ELSE
            # YYY = ( R_ALADDIN - XXX) * SIN(ANGLE)
            # END IF
            # XXX  =   COS(ANGLE) * XXX + R_ALADDIN * (1.0D0 - COS(ANGLE))


            # Synchrotron depth distribution
            # R_ALADDIN NEGATIVE FOR COUNTER-CLOCKWISE SOURCE
            if r_aladdin < 0:
                YYY = numpy.abs(r_aladdin + XXX) * numpy.sin(ANGLE)
            else:
                YYY = numpy.abs(r_aladdin - XXX) * numpy.sin(ANGLE)

            XXX = numpy.cos(ANGLE) * XXX + r_aladdin * (1.0 - numpy.cos(ANGLE))


            rays[itik,0] = XXX
            rays[itik,1] = YYY
            rays[itik,2] = ZZZ



    #         #
    #         # directions
    #         #
    #
    #         #     ! C
    #         #     ! C Synchrotron source
    #         #     ! C Note. The angle of emission IN PLANE is the same as the one used
    #         #     ! C before. This will give rise to a source curved along the orbit.
    #         #     ! C The elevation angle is instead characteristic of the SR distribution.
    #         #     ! C The electron beam emittance is included at this stage. Note that if
    #         #     ! C EPSI = 0, we'll have E_BEAM = 0.0, with no changes.
    #         #     ! C
    #         #     IF (F_WIGGLER.EQ.3) ANGLE=0        ! Elliptical Wiggler.
    #         #     ANGLEX =   ANGLE + E_BEAM(1)
    #         #     DIREC(1)  =   TAN(ANGLEX)
    #         #     IF (R_ALADDIN.LT.0.0D0) DIREC(1) = - DIREC(1)
    #         #     DIREC(2)  =   1.0D0
    #         #     ARG_ANG  =   GRID(6,ITIK)
    #
    #         ANGLEX = ANGLE + E_BEAM1
    #         DIREC1 = numpy.tan(ANGLEX)
    #         DIREC2 = 1.0
    #
    #
    #         #     ! C
    #         #     ! C In the case of SR, we take into account the fact that the electron
    #         #     ! C trajectory is not orthogonal to the field. This will give a correction
    #         #     ! C to the photon energy.  We can write it as a correction to the
    #         #     ! C magnetic field strength; this will linearly shift the critical energy
    #         #     ! C and, with it, the energy of the emitted photon.
    #         #     ! C
    #         #     E_TEMP(3) =   TAN(E_BEAM(3))/COS(E_BEAM(1))
    #         #     E_TEMP(2) =   1.0D0
    #         #     E_TEMP(1) =   TAN(E_BEAM(1))
    #         #     CALL NORM (E_TEMP,E_TEMP)
    #         #     CORREC =   SQRT(1.0D0-E_TEMP(3)**2)
    #         #     4400 CONTINUE
    #
    #         E_TEMP3 = numpy.tan(E_BEAM3)/numpy.cos(E_BEAM1)
    #         E_TEMP2 = 1.0
    #         E_TEMP1 = numpy.tan(E_BEAM1)
    #
    #         e_temp_norm = numpy.sqrt( E_TEMP1**2 + E_TEMP2**2 + E_TEMP3**2)
    #
    #         E_TEMP3 /= e_temp_norm
    #         E_TEMP2 /= e_temp_norm
    #         E_TEMP1 /= e_temp_norm
    #
    #
    #         CORREC = numpy.sqrt(1.0 - E_TEMP3**2)
    #
    #
    #         #     IF (FDISTR.EQ.6) THEN
    #         #         CALL ALADDIN1 (ARG_ANG,ANGLEV,F_POL,IER)
    #         #         Q_WAVE =   TWOPI*PHOTON(1)/TOCM*CORREC
    #         #         POL_DEG =   ARG_ANG
    #         #     ELSE IF (FDISTR.EQ.4) THEN
    #         #         ARG_ENER =   WRAN (ISTAR1)
    #         #         RAD_MIN =   ABS(R_MAGNET)
    #         #
    #         #         i1 = 1
    #         #         CALL WHITE  &
    #         #         (RAD_MIN,CORREC,ARG_ENER,ARG_ANG,Q_WAVE,ANGLEV,POL_DEG,i1)
    #         #     END IF
    #
    #         RAD_MIN = numpy.abs(R_MAGNET)
    #
    #
    #         # CALL WHITE (RAD_MIN,CORREC,ARG_ENER,ARG_ANG,Q_WAVE,ANGLEV,POL_DEG,i1)
    #         ARG_ANG = numpy.random.random()
    #         ARG_ENER = numpy.random.random()
    #
    #
    #         # print("   >> R_MAGNET, DIREC",R_MAGNET,DIREC1,DIREC2)
    #         # print("   >> RAD_MIN,CORREC,ARG_ENER,ARG_ANG,",RAD_MIN,CORREC,ARG_ENER,ARG_ANG)
    #
    #         #
    #         gamma = self.syned_electron_beam.gamma()
    #         m2ev = codata.c * codata.h / codata.e
    #         TOANGS = m2ev * 1e10
    #         critical_energy = TOANGS*3.0*numpy.power(gamma,3)/4.0/numpy.pi/1.0e10*(1.0/RAD_MIN)
    #
    #         sampled_photon_energy = sampled_energies[itik]
    #         wavelength = codata.h * codata.c / codata.e /sampled_photon_energy
    #         Q_WAVE = 2 * numpy.pi / (wavelength*1e2)
    #         # print("   >> PHOTON ENERGY, Ec, lambda, Q: ",sampled_photon_energy,critical_energy,wavelength*1e10,Q_WAVE)
    #
    #
    #         eene = sampled_photon_energy / critical_energy
    #
    #         fm = sync_f(a*1e-3*self.syned_electron_beam.gamma(),eene,polarization=0) * \
    #             numpy.power(eene,2)*a8*self.syned_electron_beam._current*hdiv_mrad * \
    #             numpy.power(self.syned_electron_beam._energy_in_GeV,2)
    #
    #         fm_s = sync_f(a*1e-3*self.syned_electron_beam.gamma(),eene,polarization=1) * \
    #             numpy.power(eene,2)*a8*self.syned_electron_beam._current*hdiv_mrad * \
    #             numpy.power(self.syned_electron_beam._energy_in_GeV,2)
    #
    #         fm_pol = numpy.zeros_like(fm)
    #         for i in range(fm_pol.size):
    #             if fm[i] == 0.0:
    #                 fm_pol[i] = 0
    #             else:
    #                 fm_pol[i] = fm_s[i] / fm[i]
    #
    #         fm.shape = -1
    #         fm_s.shape = -1
    #         fm_pol.shape = -1
    #         # print(">>>>",a.shape,fm.shape,fm_s.shape)
    #         pol_deg_interpolator = interp1d(a*1e-3,fm_pol)
    #
    #
    #         samplerAng = Sampler1D(fm,a*1e-3)
    #
    #         # samplerPol = Sampler1D(fm_s/fm,a*1e-3)
    #
    #         # plot(a*1e-3,fm_s/fm)
    #
    #         sampled_theta = samplerAng.get_sampled(ARG_ENER)
    #
    #         sampled_pol_deg = pol_deg_interpolator(sampled_theta)
    #
    #
    #         # print("sampled_theta: ",sampled_theta, "sampled_energy: ",sampled_photon_energy, "sampled pol ",sampled_pol_deg)
    #
    #         ANGLEV = sampled_theta
    #         ANGLEV += E_BEAM3
    #         #     IF (ANGLEV.LT.0.0) I_CHANGE = -1
    #         #     ANGLEV =   ANGLEV + E_BEAM(3)
    #         #     ! C
    #         #     ! C Test if the ray is within the specified limits
    #         #     ! C
    #         #     IF (FGRID.EQ.0.OR.FGRID.EQ.2) THEN
    #         #         IF (ANGLEV.GT.VDIV1.OR.ANGLEV.LT.-VDIV2) THEN
    #         #             ARG_ANG = WRAN(ISTAR1)
    #         #             ! C
    #         #             ! C If it is outside the range, then generate another ray.
    #         #             ! C
    #         #             GO TO 4400
    #         #         END IF
    #         #     END IF
    #         #     DIREC(3)  =   TAN(ANGLEV)/COS(ANGLEX)
    #
    #         DIREC3 = numpy.tan(ANGLEV) / numpy.cos(ANGLEX)
    #         #     IF (F_WIGGLER.EQ.3) THEN
    #         #         CALL ROTATE (DIREC, ANGLE3,ANGLE2,ANGLE1,DIREC)
    #         #     END IF
    #         #     CALL NORM (DIREC,DIREC)
    #
    #         direc_norm = numpy.sqrt(DIREC1**2 + DIREC2**2 + DIREC3**2)
    #
    #         DIREC1 /= direc_norm
    #         DIREC2 /= direc_norm
    #         DIREC3 /= direc_norm
    #
    #         rays[itik,3] = DIREC1 # VX
    #         rays[itik,4] = DIREC2 # VY
    #         rays[itik,5] = DIREC3 # VZ
    #
    #     if user_unit_to_m != 1.0:
    #         rays[:,0] /= user_unit_to_m
    #         rays[:,1] /= user_unit_to_m
    #         rays[:,2] /= user_unit_to_m
    #
    #     #
    #     # sample divergences (cols 4-6): the Shadow way
    #     #
    #
    #
    #
    #
    #     #
    #     # electric field vectors (cols 7-9, 16-18) and phases (cols 14-15)
    #     #
    #
    #     # ! C
    #     # ! C  ---------------------------------------------------------------------
    #     # ! C                 POLARIZATION
    #     # ! C
    #     # ! C   Generates the polarization of the ray. This is defined on the
    #     # ! C   source plane, so that A_VEC is along the X-axis and AP_VEC is along Z-axis.
    #     # ! C   Then care must be taken so that A will be perpendicular to the ray
    #     # ! C   direction.
    #     # ! C
    #     # ! C
    #     # A_VEC(1) = 1.0D0
    #     # A_VEC(2) = 0.0D0
    #     # A_VEC(3) = 0.0D0
    #
    #     DIREC = rays[:,3:6].copy()
    #     A_VEC = numpy.zeros_like(DIREC)
    #     A_VEC[:,0] = 1.0
    #
    #     # ! C
    #     # ! C   Rotate A_VEC so that it will be perpendicular to DIREC and with the
    #     # ! C   right components on the plane.
    #     # ! C
    #     # CALL CROSS (A_VEC,DIREC,A_TEMP)
    #     A_TEMP = self._cross(A_VEC,DIREC)
    #     # CALL CROSS (DIREC,A_TEMP,A_VEC)
    #     A_VEC = self._cross(DIREC,A_TEMP)
    #     # CALL NORM (A_VEC,A_VEC)
    #     A_VEC = self._norm(A_VEC)
    #     # CALL CROSS (A_VEC,DIREC,AP_VEC)
    #     AP_VEC = self._cross(A_VEC,DIREC)
    #     # CALL NORM (AP_VEC,AP_VEC)
    #     AP_VEC = self._norm(AP_VEC)
    #
    #     #
    #     # obtain polarization for each ray (interpolation)
    #     #
    #
    #
    #
    #     POL_DEG = sampled_pol_deg
    #     DENOM = numpy.sqrt(1.0 - 2.0 * POL_DEG + 2.0 * POL_DEG**2)
    #     AX = POL_DEG/DENOM
    #     for i in range(3):
    #         A_VEC[:,i] *= AX
    #
    #     AZ = (1.0-POL_DEG)/DENOM
    #     for i in range(3):
    #         AP_VEC[:,i] *= AZ
    #
    #
    #     rays[:,6:9] =  A_VEC
    #     rays[:,15:18] = AP_VEC
    #
    #     #
    #     # ! C
    #     # ! C Now the phases of A_VEC and AP_VEC.
    #     # ! C
    #
    #     #
    #     POL_ANGLE = 0.5 * numpy.pi
    #
    #     if F_COHER == 1:
    #         PHASEX = 0.0
    #     else:
    #         PHASEX = numpy.random.random(NRAYS) * 2 * numpy.pi
    #
    #     # PHASEZ = PHASEX + POL_ANGLE * numpy.sign(ANGLEV)
    #
    #     rays[:,13] = 0.0 # PHASEX
    #     rays[:,14] = 0.0 # PHASEZ
    #
    #     # set flag (col 10)
    #     rays[:,9] = 1.0
    #
    #     #
    #     # photon energy (col 11)
    #     #
    #
    #     # A2EV = 2.0*numpy.pi/(codata.h*codata.c/codata.e*1e2)
    #     sampled_photon_energy = sampled_energies
    #     wavelength = codata.h * codata.c / codata.e /sampled_photon_energy
    #     Q_WAVE = 2 * numpy.pi / (wavelength*1e2)
    #     rays[:,10] =  Q_WAVE # sampled_photon_energy * A2EV
    #
    #     # col 12 (ray index)
    #     rays[:,11] =  1 + numpy.arange(NRAYS)
    #
    #     # col 13 (optical path)
    #     rays[:,11] = 0.0
    #


        # plot_scatter(rays[:,0]*1e6,rays[:,2]*1e6,xtitle="X um",ytitle="Z um")
        # plot_scatter(rays[:,1],rays[:,0]*1e6,xtitle="Y m",ytitle="X um")
        # plot_scatter(rays[:,1],rays[:,2]*1e6,xtitle="Y m",ytitle="Z um")
        return rays
    #
    # def _cross(self,u,v):
    #     # w = u X v
    #     # u = array (npoints,vector_index)
    #
    #     w = numpy.zeros_like(u)
    #     w[:,0] = u[:,1] * v[:,2] - u[:,2] * v[:,1]
    #     w[:,1] = u[:,2] * v[:,0] - u[:,0] * v[:,2]
    #     w[:,2] = u[:,0] * v[:,1] - u[:,1] * v[:,0]
    #
    #     return w
    #
    # def _norm(self,u):
    #     # w = u / |u|
    #     # u = array (npoints,vector_index)
    #     u_norm = numpy.zeros_like(u)
    #     uu = numpy.sqrt( u[:,0]**2 + u[:,1]**2 + u[:,2]**2)
    #     for i in range(3):
    #         u_norm[:,i] = uu
    #     return u / u_norm
    #
    # def _sample_photon_energy_theta_and_phi(self,NRAYS):
    #
    #     #
    #     # sample divergences
    #     #
    #     return 0,0,0

if __name__ == "__main__":

    bm = SourceBendingMagnet()

