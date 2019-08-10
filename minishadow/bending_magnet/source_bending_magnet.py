__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "30-08-2018"

"""


"""


import numpy

from srxraylib.util.inverse_method_sampler import Sampler1D, Sampler2D

import scipy.constants as codata
from scipy import interpolate
from scipy.interpolate import interp1d,interp2d

from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet
from syned.storage_ring.electron_beam import ElectronBeam

from srxraylib.sources.srfunc import sync_ene, sync_ang

# from srxraylib.plot.gol import plot,plot_scatter,plot_image


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

    def is_monochromatic(self):
        if self._NG_E == 1:
            return True
        if self._EMAX == self._EMIN:
            return True
        return False

    def calculate_rays(self, F_COHER=0, NRAYS=5000, SEED=123456,
                       EPSI_DX=0.0, EPSI_DZ=0.0,
                       psi_interval_in_units_one_over_gamma=None,
                       psi_interval_number_of_points=1001,
                       verbose=False):
        """
        compute the rays in SHADOW matrix (shape (npoints,18) )
        :param F_COHER: set this flag for coherent beam
        :param user_unit_to_m: default 1.0 (m)
        :return: rays, a numpy.array((npoits,18))
        """

        if SEED != 0:
            numpy.random.seed(SEED)


        sigmas = self.syned_electron_beam.get_sigmas_all()

        rays = numpy.zeros((NRAYS,18))

        RAD_MIN= numpy.abs(self.syned_bending_magnet._radius)
        RAD_MAX= numpy.abs(self.syned_bending_magnet._radius)

        # r_aladdin	=  bending magnet radius in units of length used for source size, CCW rings negative.

        r_aladdin = self.syned_bending_magnet._radius

        if r_aladdin < 0:
            POL_ANGLE = -90.0 * numpy.pi / 2
        else:
            POL_ANGLE = 90.0 * numpy.pi / 2

        HDIV1 = 0.5 * self.syned_bending_magnet.horizontal_divergence()
        HDIV2 = HDIV1

        gamma = self.syned_electron_beam.gamma()
        critical_energy = self.syned_bending_magnet.get_critical_energy(self.syned_electron_beam.energy())

        if psi_interval_in_units_one_over_gamma is None:
            c = numpy.array([-0.3600382, 0.11188709])  # see file fit_psi_interval.py
            x = numpy.log10(self._EMIN / critical_energy)
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

        if self.is_monochromatic():
            if verbose:
                print(">>> calculate_rays: is monochromatic")
                print(">>> calculate_rays: sync_ang (s) E=%f GeV, I=%f A, D=%f mrad, R=%f m, PhE=%f eV, Ec=%f eV, PhE/Ec=%f "% ( \
                    self.syned_electron_beam.energy(),
                    self.syned_electron_beam.current(),
                    (HDIV1 + HDIV2) * 1e3,
                    self.syned_bending_magnet._radius,  # not needed anyway
                    self._EMIN,
                    critical_energy,
                    self._EMIN/critical_energy,
                      ))


            angular_distribution_s = sync_ang(1,#Flux at a given photon energy
                                            angle_array_mrad,
                                            polarization=1,#1 Parallel (l2=1, l3=0, in Sokolov&Ternov notation)
                                            e_gev=self.syned_electron_beam.energy(),
                                            i_a=self.syned_electron_beam.current(),
                                            hdiv_mrad=(HDIV1+HDIV2)*1e3,
                                            r_m=self.syned_bending_magnet._radius,#not needed anyway
                                            energy=self._EMIN,
                                            ec_ev=critical_energy)


            if verbose:
                print(">>> calculate_rays: sync_ang (p)")

            angular_distribution_p = sync_ang(1,#Flux at a given photon energy
                                            angle_array_mrad,
                                            polarization=2,#1 Parallel (l2=1, l3=0, in Sokolov&Ternov notation)
                                            e_gev=self.syned_electron_beam.energy(),
                                            i_a=self.syned_electron_beam.current(),
                                            hdiv_mrad=(HDIV1+HDIV2)*1e3,
                                            r_m=self.syned_bending_magnet._radius,#not needed anyway
                                            energy=self._EMIN,
                                            ec_ev=critical_energy)

            angular_distribution_s = angular_distribution_s.flatten()
            angular_distribution_p = angular_distribution_p.flatten()

            if verbose:
                from srxraylib.plot.gol import plot
                plot(angle_array_mrad,angular_distribution_s,
                     angle_array_mrad,angular_distribution_p,xtitle="angle / mrad",legend=["s","p"])

            sampler_angle = Sampler1D(angular_distribution_s+angular_distribution_p,angle_array_mrad*1e-3)
            if verbose:
                print(">>> calculate_rays: get_n_sampled_points (angle)")
            sampled_angle = sampler_angle.get_n_sampled_points(NRAYS)
            if verbose:
                print(">>> calculate_rays: DONE get_n_sampled_points (angle)  %d points"%(sampled_angle.size))

            pol_deg_interpolator = interp1d(angle_array_mrad*1e-3,
                    angular_distribution_s/(angular_distribution_s+angular_distribution_p))
            sampled_polarization = pol_deg_interpolator(sampled_angle)

            sampled_photon_energy = numpy.zeros_like(sampled_angle) + self._EMIN

        else: # polychromatic

            photon_energy_array = numpy.linspace(self._EMIN,self._EMAX,self._NG_E)

            if verbose:
                print(">>> sync_ene: calculating energy distribution")

            fm_s = sync_ene(4,photon_energy_array,
                          ec_ev=self.syned_bending_magnet.get_critical_energy(self.syned_electron_beam.energy()),
                          e_gev=self.syned_electron_beam.energy(),
                          i_a=self.syned_electron_beam.current(),
                          hdiv_mrad=1,
                          psi_min=angle_array_mrad.min(),
                          psi_max=angle_array_mrad.max(),
                          psi_npoints=angle_array_mrad.size,
                          polarization=1)

            fm_p = sync_ene(4,photon_energy_array,
                          ec_ev=self.syned_bending_magnet.get_critical_energy(self.syned_electron_beam.energy()),
                          e_gev=self.syned_electron_beam.energy(),
                          i_a=self.syned_electron_beam.current(),
                          hdiv_mrad=1,
                          psi_min=angle_array_mrad.min(),
                          psi_max=angle_array_mrad.max(),
                          psi_npoints=angle_array_mrad.size,
                          polarization=2)

            fm = fm_s + fm_p

            if verbose:
                print(">>> DONE sync_ene: calculating energy distribution",photon_energy_array.shape,fm.shape)
                from srxraylib.plot.gol import plot,plot_image
                plot(photon_energy_array,fm[fm.shape[0]//2,:],xtitle="Energy / eV",ytitle="Flux at zero elevation")
                plot(angle_array_mrad, fm[:,0], xtitle="Angle / mrad", ytitle="Flux at Emin="%(photon_energy_array[0]))
                print(">>>>>>>",fm.shape,angle_array_mrad.shape,photon_energy_array.shape)
                plot_image(fm,angle_array_mrad,photon_energy_array,aspect='auto',show=0,title="flux",xtitle="Psi / mrad",ytitle="Energy / eV")
                plot_image(fm_s/fm,angle_array_mrad,photon_energy_array,aspect='auto',title="polarization",xtitle="Psi / mrad",ytitle="Energy / eV")

            fm1 = numpy.zeros_like(fm)
            for i in range(fm.shape[0]):
                fm1[i,:] = fm[i,:] / (photon_energy_array*0.001)  # in photons/ev

            # plot_image(fm,angle_array_mrad,photon_energy_array,aspect='auto',show=0)
            # plot_image(fm_s/fm,angle_array_mrad,photon_energy_array,aspect='auto',title="polarization")


            sampler2 = Sampler2D(fm1,angle_array_mrad*1e-3,photon_energy_array)
            sampled_angle,sampled_photon_energy = sampler2.get_n_sampled_points(NRAYS)


            Angle_array_mrad = numpy.outer(angle_array_mrad,numpy.ones_like(photon_energy_array))
            Photon_energy_array = numpy.outer(numpy.ones_like(angle_array_mrad),photon_energy_array)
            Pi = numpy.array([Angle_array_mrad.flatten()*1e-3, Photon_energy_array.flatten()]).transpose()

            P = numpy.array([sampled_angle, sampled_photon_energy]).transpose()
            sampled_polarization = interpolate.griddata(Pi, (fm_s/fm).flatten(), P, method = "cubic")

        for itik in range(NRAYS):
            # ! Synchrontron depth
            ANGLE  =  numpy.random.random() * (HDIV1 + HDIV2) - HDIV2
            EPSI_PATH =  numpy.abs(r_aladdin) * ANGLE



            if self._FLAG_EMITTANCE:
                sigma_x, sigma_xp, sigma_z, sigma_zp = self.syned_electron_beam.get_sigmas_all()

                # ! calculation of the electrom beam moments at the current position
                # ! (sX,sZ) = (epsi_wx,epsi_ez):
                # ! <x2> = sX^2 + sigmaX^2
                # ! <x x'> = sX sigmaXp^2
                # ! <x'2> = sigmaXp^2                 (same for Z)

                epsi_wX = EPSI_DX + EPSI_PATH # sigma_x * sigma_xp


                # ! C
                # ! C Compute the actual distance (EPSI_W*) from the orbital focus
                # ! C
                # EPSI_WX = EPSI_DX + EPSI_PATH
                # EPSI_WZ = EPSI_DZ + EPSI_PATH

                rSigmaX = numpy.sqrt( (epsi_wX**2) * (sigma_xp**2) + sigma_x**2 )
                rSigmaXp = sigma_xp
                if rSigmaX * rSigmaXp != 0.0:
                    rhoX = epsi_wX * sigma_xp**2 / (rSigmaX * rSigmaXp)
                else:
                    rhoX = 0.0
                mean = [0, 0]
                cov = [[rSigmaX**2, rhoX*rSigmaX*rSigmaXp], [rhoX*rSigmaX*rSigmaXp, rSigmaXp**2]]  # diagonal covariance
                sampled_x, sampled_xp = numpy.random.multivariate_normal(mean, cov, 1).T
                # plot_scatter(sampled_x,sampled_xp,title="X")
                XXX = sampled_x
                E_BEAM1 = sampled_xp


                epsi_wZ = EPSI_DZ + EPSI_PATH # sigma_z * sigma_zp
                rSigmaZ = numpy.sqrt( (epsi_wZ**2) * (sigma_zp**2) + sigma_z**2 )
                rSigmaZp = sigma_zp
                if rSigmaZ * rSigmaZp != 0.0:
                    rhoZ = epsi_wZ * sigma_zp**2 / (rSigmaZ * rSigmaZp)
                else:
                    rhoZ = 0.0
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

            # ! C
            # ! C Synchrotron source
            # ! C Note. The angle of emission IN PLANE is the same as the one used
            # ! C before. This will give rise to a source curved along the orbit.
            # ! C The elevation angle is instead characteristic of the SR distribution.
            # ! C The electron beam emittance is included at this stage. Note that if
            # ! C EPSI = 0, we'll have E_BEAM = 0.0, with no changes.
            # ! C
            # ANGLEX =   ANGLE + E_BEAM(1)
            # DIREC(1)  =   TAN(ANGLEX)
            # IF (R_ALADDIN.LT.0.0D0) DIREC(1) = - DIREC(1)
            # DIREC(2)  =   1.0D0
            # ARG_ANG  =   GRID(6,ITIK)

            ANGLEX = ANGLE + E_BEAM1
            DIREC1 = numpy.tan(ANGLEX)
            if r_aladdin < 0:
                DIREC1 *= -1.0
            DIREC2 = 1.0
            ARG_ANG = numpy.random.random()

            # ! C
            # ! C In the case of SR, we take into account the fact that the electron
            # ! C trajectory is not orthogonal to the field. This will give a correction
            # ! C to the photon energy.  We can write it as a correction to the
            # ! C magnetic field strength; this will linearly shift the critical energy
            # ! C and, with it, the energy of the emitted photon.
            # ! C
            # E_TEMP(3) =   TAN(E_BEAM(3))/COS(E_BEAM(1))
            # E_TEMP(2) =   1.0D0
            # E_TEMP(1) =   TAN(E_BEAM(1))
            # CALL NORM (E_TEMP,E_TEMP)
            # CORREC =   SQRT(1.0D0-E_TEMP(3)**2)
            # 4400 CONTINUE
            E_TEMP3 = numpy.tan(E_BEAM3) / numpy.cos(E_BEAM1)
            E_TEMP2 = 1.0
            E_TEMP1 = numpy.tan(E_BEAM1)
            E_TEMP_MOD = numpy.sqrt(E_TEMP1**2 + E_TEMP2**2 + E_TEMP3**2)
            E_TEMP3 /= E_TEMP_MOD
            E_TEMP2 /= E_TEMP_MOD
            E_TEMP1 /= E_TEMP_MOD

            # IF (FDISTR.EQ.6) THEN ! exect synchtotron
            #     CALL ALADDIN1 (ARG_ANG,ANGLEV,F_POL,IER)
            #     Q_WAVE =   TWOPI*PHOTON(1)/TOCM*CORREC
            #     POL_DEG =   ARG_ANG
            # ELSE IF (FDISTR.EQ.4) THEN  ! synchrotron
            #     print*,"R_MAGNET, DIREC",R_MAGNET,DIREC
            #     ARG_ENER =   WRAN (ISTAR1)
            #     RAD_MIN =   ABS(R_MAGNET)
            #
            #     i1 = 1
            #     arg_ener = 0.5
            #     arg_ang = 0.5
            #     CALL WHITE (RAD_MIN,CORREC,ARG_ENER,ARG_ANG,Q_WAVE,ANGLEV,POL_DEG,i1)
            #
            #     print*,"RAD_MIN,CORREC,ARG_ENER,ARG_ANG,Q_WAVE,ANGLEV,POL_DEG",RAD_MIN,CORREC,ARG_ENER,ARG_ANG,Q_WAVE,ANGLEV,POL_DEG
            #     !Q_WAVE =   TWOPI*PHOTON(1)/TOCM*CORREC
            #     print*,"ENER,ANGLEV: ",Q_WAVE*TOCM/TWOPI,ANGLEV
            # END IF

            # interpolate for the photon energy,vertical angle,and the degree of polarization.

            wavelength = codata.h * codata.c / codata.e / sampled_photon_energy[itik]
            Q_WAVE = 2 * numpy.pi / (wavelength*1e2)
            ANGLEV = sampled_angle[itik]
            POL_DEG = sampled_polarization[itik]


            # IF (ANGLEV.LT.0.0) I_CHANGE = -1
            # ANGLEV =   ANGLEV + E_BEAM(3)
            if ANGLEV < 0:
                I_CHANGE = -1
            ANGLEV += E_BEAM3

            # ------ NOT LONGER DONE ------
            # ! C
            # ! C Test if the ray is within the specified limits
            # ! C
            # IF (FGRID.EQ.0.OR.FGRID.EQ.2) THEN
            #     IF (ANGLEV.GT.VDIV1.OR.ANGLEV.LT.-VDIV2) THEN
            #         ARG_ANG = WRAN(ISTAR1)
            #         ! C
            #         ! C If it is outside the range, then generate another ray.
            #         ! C
            #         GO TO 4400
            #     END IF
            # END IF


            # DIREC(3)  =   TAN(ANGLEV)/COS(ANGLEX)
            # CALL NORM (DIREC,DIREC)

            DIREC3 = numpy.tan(ANGLEV) / numpy.cos(ANGLEX)

            DIREC_MOD = numpy.sqrt(DIREC1**2 + DIREC2**2 + DIREC3**2)
            DIREC3 /= DIREC_MOD
            DIREC2 /= DIREC_MOD
            DIREC1 /= DIREC_MOD

            # print(">>>>DIREC,FGRID,R_ALADDIN: ",itik,DIREC1,DIREC2,DIREC3)

            rays[itik,3] = DIREC1
            rays[itik,4] = DIREC2
            rays[itik,5] = DIREC3


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


        POL_DEG = sampled_polarization
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
        # sampled_photon_energy = sampled_photon_energy
        wavelength = codata.h * codata.c / codata.e /sampled_photon_energy
        Q_WAVE = 2 * numpy.pi / (wavelength*1e2)
        rays[:,10] =  Q_WAVE # sampled_photon_energy * A2EV

        # col 12 (ray index)
        rays[:,11] =  1 + numpy.arange(NRAYS)

        # col 13 (optical path)
        rays[:,11] = 0.0

        POL_ANGLE = 0.5 * numpy.pi

        if F_COHER == 1:
            PHASEX = 0.0
        else:
            PHASEX = numpy.random.random(NRAYS) * 2 * numpy.pi

        PHASEZ = PHASEX + POL_ANGLE * numpy.sign(ANGLEV)

        rays[:,13] = PHASEX
        rays[:,14] = PHASEZ


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


if __name__ == "__main__":

    from srxraylib.plot.gol import plot_scatter

    syned_electron_beam = ElectronBeam(energy_in_GeV=6.04,current=0.2,
                                       moment_xx=(0.0078e-2)**2,
                                       moment_xpxp=(3.8e-07/0.0078)**2,
                                       moment_yy=(0.0036*1e-2)**2,
                                       moment_ypyp=(3.8e-09/0.0036)**2,
                                       )

    syned_bending_magnet = BendingMagnet(radius=25.1772,magnetic_field=0.8,length=25.1772*0.001)

    emin = 5000.0                # Photon energy scan from energy (in eV)
    emax = 100000.0              # Photon energy scan to energy (in eV)
    ng_e = 51                    # Photon energy scan number of points
    ng_j = 20                    # Number of points in electron trajectory (per period) for internal calculation only
    flag_emittance = 1           # when sampling rays: Use emittance (0=No, 1=Yes)



    bm = SourceBendingMagnet(syned_electron_beam=syned_electron_beam,
                 syned_bending_magnet=syned_bending_magnet,
                 emin=emin,               # Photon energy scan from energy (in eV)
                 emax=emax,               # Photon energy scan to energy (in eV)
                 ng_e=ng_e,                    # Photon energy scan number of points
                 ng_j=ng_j,                    # Number of points in electron trajectory (per period) for internal calculation only
                 flag_emittance=flag_emittance,           # when sampling rays: Use emittance (0=No, 1=Yes)
                )

    print(bm.info())

    rays = bm.calculate_rays(F_COHER=0,NRAYS=5000,SEED=123456,EPSI_DX=0.0,EPSI_DZ=0.0)


    plot_scatter(rays[:,0]*1e6,rays[:,2]*1e6,xtitle="X um",ytitle="Z um")
    plot_scatter(rays[:,1],rays[:,0]*1e6,xtitle="Y m",ytitle="X um")
    plot_scatter(rays[:,1],rays[:,2]*1e6,xtitle="Y m",ytitle="Z um")
    plot_scatter(rays[:,3]*1e6,rays[:,5]*1e6,xtitle="X' urad",ytitle="Z' urad")
