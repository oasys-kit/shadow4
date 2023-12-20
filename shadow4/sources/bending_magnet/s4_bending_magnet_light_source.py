"""
Bending magnet light source.
"""
import numpy

from srxraylib.util.inverse_method_sampler import Sampler1D, Sampler2D

import scipy.constants as codata
from scipy.interpolate import interp1d
from scipy.interpolate import griddata, CloughTocher2DInterpolator

from syned.storage_ring.electron_beam import ElectronBeam

from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet
from shadow4.beam.s4_beam import S4Beam
from shadow4.sources.s4_light_source import S4LightSource

from shadow4.tools.arrayofvectors import vector_cross, vector_norm

from shadow4.tools.sync_f_sigma_and_pi import sync_f_sigma_and_pi
import time

from srxraylib.sources.srfunc import sync_ene

class S4BendingMagnetLightSource(S4LightSource):
    """
    Defines a bending magnet light source and implements the mechanism of sampling rays.

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
                         electron_beam=electron_beam if not electron_beam is None else ElectronBeam(),
                         magnetic_structure=magnetic_structure if not magnetic_structure is None else S4BendingMagnet(),
                         nrays=nrays,
                         seed=seed,
                         )

    def get_beam(self, F_COHER=0,
                       EPSI_DX=0.0,
                       EPSI_DZ=0.0,
                       psi_interval_in_units_one_over_gamma=None,
                       psi_interval_number_of_points=1001,
                       verbose=False):
        """
        Creates the beam as emitted by the bending magnet.

        Parameters
        ----------
        F_COHER : int, optional
            A flag for coherent (1) or incoherent (0) rays.
        EPSI_DX : float, optional
            The distance from waist in X.
        EPSI_DZ : float, optional
            The distance from waist in Z.
        psi_interval_in_units_one_over_gamma : int, optional
            The sampling interval for the vertical divergence in units if 1/gamma.
            If None, it is calculated automatically.
        psi_interval_number_of_points : int, optional
            The number of points for the vertical divergency scan.
        verbose : int, optional
            set to 1 for verbose output.

        Returns
        -------
        instance od S4beam
        """
        return S4Beam.initialize_from_array(self.__calculate_rays(F_COHER=F_COHER,
                                                                  EPSI_DX=EPSI_DX,
                                                                  EPSI_DZ=EPSI_DZ,
                                                                  psi_interval_in_units_one_over_gamma=psi_interval_in_units_one_over_gamma,
                                                                  psi_interval_number_of_points=psi_interval_number_of_points,
                                                                  verbose=verbose))


    def get_info(self):
        """
        Returns the specific information for the bending magnet light source.

        Returns
        -------
        str
        """
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
        txt += "Input Bending Magnet parameters: \n"
        txt += "        radius: %f m\n"%magnetic_structure._radius
        txt += "        magnetic field: %f T\n"%magnetic_structure._magnetic_field
        txt += "        length: %f m\n"%magnetic_structure._length
        txt += "Lorentz factor (gamma): %f\n"%electron_beam.gamma()
        txt += "-----------------------------------------------------\n"
        return (txt)

    def __calculate_rays(self,
                         F_COHER=0,
                         EPSI_DX=0.0,
                         EPSI_DZ=0.0,
                         psi_interval_in_units_one_over_gamma=None,
                         psi_interval_number_of_points=1001,
                         verbose=0):
        #
        # compute the rays in SHADOW matrix (shape (npoints,18) )
        #
        t0 = time.time()

        # retrieve parameters
        NRAYS = self.get_nrays()

        if self.get_seed() != 0:
            numpy.random.seed(self.get_seed())

        r_aladdin = self.get_magnetic_structure().radius()

        if r_aladdin < 0:
            POL_ANGLE = -numpy.pi / 2
        else:
            POL_ANGLE = numpy.pi / 2

        HDIV1 = 0.5 * self.get_magnetic_structure().horizontal_divergence()
        HDIV2 = HDIV1

        gamma = self.get_electron_beam().gamma()
        critical_energy = self.get_magnetic_structure().get_critical_energy(self.get_electron_beam().energy())

        if psi_interval_in_units_one_over_gamma is None:
            c = numpy.array([-0.3600382, 0.11188709])  # see file fit_psi_interval.py
            x = numpy.log10(self.get_magnetic_structure()._EMIN / critical_energy)
            y_fit = c[1] + c[0] * x
            psi_interval_in_units_one_over_gamma = 10**y_fit  # this is the semi interval
            psi_interval_in_units_one_over_gamma *= 4  # doubled interval
            if psi_interval_in_units_one_over_gamma < 2:
                psi_interval_in_units_one_over_gamma = 2

        if verbose:
            print(">>> psi_interval_in_units_one_over_gamma: ",psi_interval_in_units_one_over_gamma)


        angle_array_mrad = numpy.linspace(-0.5 * psi_interval_in_units_one_over_gamma * 1e3 / gamma,
                                          0.5 * psi_interval_in_units_one_over_gamma * 1e3 / gamma,
                                          psi_interval_number_of_points)


        # initialize arrays
        rays = numpy.zeros((NRAYS,18))
        anglev_sign = numpy.zeros(NRAYS)

        # calculate the sampled_angle and sampled_photon_energy.
        # separate the calculation for monochromatic case (using Sampler1D for angles) and
        # polychromatic case (use Sampler2D for angles and energies).
        t1 = time.time()

        if self.get_magnetic_structure().is_monochromatic():
            if verbose:
                print(">>> calculate_rays: is monochromatic")
                print(">>> calculate_rays: (s) E=%f GeV, I=%f A, D=%f mrad, R=%f m, PhE=%f eV, Ec=%f eV, PhE/Ec=%f "% ( \
                    self.get_electron_beam().energy(),
                    self.get_electron_beam().current(),
                    (HDIV1 + HDIV2) * 1e3,
                    self.get_magnetic_structure().radius(),  # not needed anyway
                    self.get_magnetic_structure()._EMIN,
                    critical_energy,
                    self.get_magnetic_structure()._EMIN/critical_energy,
                      ))
            t2 = time.time()

            e_gev = self.get_electron_beam().energy()
            i_a = self.get_electron_beam().current()
            hdiv_mrad = (HDIV1 + HDIV2) * 1e3
            energy = self.get_magnetic_structure()._EMIN
            ec_ev = critical_energy

            angle_mrad = angle_array_mrad
            codata_mee = 1e-6 * codata.m_e * codata.c ** 2 / codata.e

            eene = energy / ec_ev
            gamma = e_gev * 1e3 / codata_mee
            angular_distribution_s, angular_distribution_p = sync_f_sigma_and_pi(angle_mrad * gamma / 1e3, eene)

            a8 = codata.e / numpy.power(codata_mee, 2) / codata.h * (9e-2 / 2 / numpy.pi) # a8 = 1.3264d13
            angular_distribution_s *= eene**2 * a8 * i_a * hdiv_mrad * e_gev**2
            angular_distribution_p *= eene**2 * a8 * i_a * hdiv_mrad * e_gev**2

            if verbose:
                from srxraylib.plot.gol import plot
                plot(angle_array_mrad,angular_distribution_s,
                     angle_array_mrad,angular_distribution_p,
                     xtitle="angle / mrad", legend=["s","p"])

            t3 = time.time()
            sampler_angle = Sampler1D(angular_distribution_s + angular_distribution_p, angle_array_mrad * 1e-3)

            if verbose: print(">>> calculate_rays: get_n_sampled_points (angle)")
            sampled_angle = sampler_angle.get_n_sampled_points(NRAYS)
            if verbose: print(">>> calculate_rays: DONE get_n_sampled_points (angle)  %d points"%(sampled_angle.size))

            sampled_photon_energy = numpy.zeros_like(sampled_angle) + self.get_magnetic_structure()._EMIN

        else: # polychromatic

            photon_energy_array = numpy.linspace(self.get_magnetic_structure()._EMIN,
                                                 self.get_magnetic_structure()._EMAX,
                                                 self.get_magnetic_structure()._NG_E)

            t2 = time.time()
            # energy array
            energy_ev = numpy.array(photon_energy_array)
            ec_ev = self.get_magnetic_structure().get_critical_energy(
                self.get_electron_beam().energy())
            eene = energy_ev / ec_ev

            # angle array
            e_gev = self.get_electron_beam().energy()
            codata_mee = 1e-6 * codata.m_e * codata.c ** 2 / codata.e
            gamma = e_gev * 1e3 / codata_mee
            angle_mrad = numpy.array(angle_array_mrad) * gamma / 1e3


            eene2 = numpy.outer(numpy.ones_like(angle_mrad), eene)
            angle_mrad2 = numpy.outer(angle_mrad, numpy.ones_like(eene))

            a5_s, a5_p = sync_f_sigma_and_pi(angle_mrad2, eene2)

            i_a = self.get_electron_beam().current()
            hdiv_mrad = 1
            a8 = codata.e / numpy.power(codata_mee, 2) / codata.h * (9e-2 / 2 / numpy.pi)
            fm_s = a5_s * eene2**2 * a8 * i_a * hdiv_mrad * e_gev**2
            fm_p = a5_p * eene2**2 * a8 * i_a * hdiv_mrad * e_gev**2

            fm = fm_s + fm_p

            # todo: remove? (we store it for external plots)
            self.angle_array_mrad = angle_array_mrad
            self.fm = fm
            self.fm_s = fm_s
            self.fm_p = fm_p

            if verbose:
                print(angle_array_mrad.shape, photon_energy_array.shape, fm_s.shape, fm_p.shape)
                print(">>> DONE : calculating energy distribution",photon_energy_array.shape,fm.shape)
                from srxraylib.plot.gol import plot,plot_image
                plot(photon_energy_array, fm[fm.shape[0] // 2, :],
                     xtitle="Energy / eV", ytitle="Flux at zero elevation")
                plot(
                     angle_array_mrad, fm[:, 0],
                     angle_array_mrad, fm_s[:, 0],
                     angle_array_mrad, fm_p[:, 0],
                     xtitle="Angle / mrad", ytitle="Flux at Emin=%f eV"%(photon_energy_array[0]))

                plot_image(fm,angle_array_mrad,photon_energy_array, aspect='auto', show=0, title="flux (total)",
                           xtitle="Psi / mrad", ytitle="Energy / eV")
                plot_image(fm_s, angle_array_mrad, photon_energy_array, aspect='auto', show=0, title="flux (sigma)",
                           xtitle="Psi / mrad", ytitle="Energy / eV")
                plot_image(fm_p, angle_array_mrad, photon_energy_array, aspect='auto', show=0, title="flux (pi)",
                           xtitle="Psi / mrad", ytitle="Energy / eV")
                plot_image(fm*0+1-fm_s/fm,angle_array_mrad,photon_energy_array,aspect='auto',title="polarization-p",xtitle="Psi / mrad",ytitle="Energy / eV")

            fm1 = numpy.zeros_like(fm)
            for i in range(fm.shape[0]):
                fm1[i,:] = fm[i, :] / (photon_energy_array * 0.001)  # in photons/ev


            sampler2 = Sampler2D(fm1, angle_array_mrad * 1e-3, photon_energy_array)
            sampled_angle, sampled_photon_energy = sampler2.get_n_sampled_points(NRAYS)

            # Angle_array_mrad = numpy.outer(angle_array_mrad,numpy.ones_like(photon_energy_array))
            # Photon_energy_array = numpy.outer(numpy.ones_like(angle_array_mrad),photon_energy_array)
            # Pi = numpy.array([Angle_array_mrad.flatten() * 1e-3, Photon_energy_array.flatten()]).transpose()
            t3 = time.time()

        t4 = time.time()

        # sample points in the electron phase space
        E_BEAM1_array = numpy.zeros(NRAYS)
        E_BEAM3_array = numpy.zeros(NRAYS)
        E_BEAMXXX_array = numpy.zeros(NRAYS)
        E_BEAMZZZ_array = numpy.zeros(NRAYS)

        ANGLE_array = numpy.random.random(NRAYS) * (HDIV1 + HDIV2) - HDIV2

        if self.get_magnetic_structure()._FLAG_EMITTANCE:
            sigma_x, sigma_xp, sigma_z, sigma_zp = self.get_electron_beam().get_sigmas_all()

            # emittance loop

            # unchanged values of the covariance matrix
            meanX = [0, 0]
            meanZ = [0, 0]
            rSigmaXp = sigma_xp
            rSigmaZp = sigma_zp
            for itik in range(NRAYS):
                EPSI_PATH = numpy.abs(r_aladdin) * ANGLE_array[itik]

                # ! C Compute the actual distance (EPSI_W*) from the orbital focus
                epsi_wX = EPSI_DX + EPSI_PATH
                epsi_wZ = EPSI_DZ + EPSI_PATH

                # ! calculation of the electrom beam moments at the current position
                # ! (sX,sZ) = (epsi_wx,epsi_ez):
                # ! <x2> = sX^2 + sigmaX^2
                # ! <x x'> = sX sigmaXp^2
                # ! <x'2> = sigmaXp^2                 (same for Z)
                rSigmaX = numpy.sqrt( (epsi_wX**2) * (sigma_xp**2) + sigma_x**2 )

                rSigmaZ = numpy.sqrt( (epsi_wZ**2) * (sigma_zp**2) + sigma_z**2 )


                if rSigmaX * rSigmaXp != 0.0:
                    rhoX = epsi_wX * sigma_xp**2 / (rSigmaX * rSigmaXp)
                else:
                    rhoX = 0.0

                if rSigmaZ * rSigmaZp != 0.0:
                    rhoZ = epsi_wZ * sigma_zp**2 / (rSigmaZ * rSigmaZp)
                else:
                    rhoZ = 0.0

                covX = [[rSigmaX**2, rhoX * rSigmaX * rSigmaXp], [rhoX * rSigmaX * rSigmaXp, rSigmaXp**2]]  # diagonal covariance
                covZ = [[rSigmaZ**2, rhoZ * rSigmaZ * rSigmaZp], [rhoZ * rSigmaZ * rSigmaZp, rSigmaZp**2]]  # diagonal covariance

                # sampling using a multivariare (2) normal distribution
                # multivariate_normal is very slow.
                # See https://github.com/numpy/numpy/issues/14193 and shadow4-tests/shadow4tests/devel/check_multivariate_normal.py
                # for acceleration recipee. For 5k rays this emittance loop passed from 31% to 8% of the total time.
                # sampled_x, sampled_xp = numpy.random.multivariate_normal(meanX, covX, 1).T
                # sampled_z, sampled_zp = numpy.random.multivariate_normal(meanZ, covZ, 1).T
                sampled_x, sampled_xp = (meanX + numpy.linalg.cholesky(covX) @ numpy.random.standard_normal(len(meanX))).T
                sampled_z, sampled_zp = (meanZ + numpy.linalg.cholesky(covZ) @ numpy.random.standard_normal(len(meanZ))).T

                XXX = sampled_x
                E_BEAM1 = sampled_xp

                ZZZ = sampled_z
                E_BEAM3 = sampled_zp

                E_BEAM1_array[itik] = E_BEAM1
                E_BEAM3_array[itik] = E_BEAM3
                E_BEAMXXX_array[itik] = XXX
                E_BEAMZZZ_array[itik] = ZZZ


            if verbose:
                from srxraylib.plot.gol import plot_scatter
                plot_scatter(E_BEAMXXX_array,E_BEAM1_array,title="sampled electron X Xp")
                plot_scatter(E_BEAMZZZ_array,E_BEAM3_array,title="sampled electron Z Zp")

        t5 = time.time()

        # spatial coordinates
        do_loop = 0 # todo: clean the part with the loop
        if do_loop:
            for itik in range(NRAYS):
                # retrieve the electron sampled values
                ANGLE = ANGLE_array[itik]
                XXX = E_BEAMXXX_array[itik]
                E_BEAM1 = E_BEAM1_array[itik]
                ZZZ = E_BEAMZZZ_array[itik]
                E_BEAM3 = E_BEAM3_array[itik]

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
        else:
            XXX = E_BEAMXXX_array
            ZZZ = E_BEAMZZZ_array
            # Synchrotron depth distribution
            # R_ALADDIN NEGATIVE FOR COUNTER-CLOCKWISE SOURCE
            if r_aladdin < 0:
                YYY = numpy.abs(r_aladdin + XXX) * numpy.sin(ANGLE_array)
            else:
                YYY = numpy.abs(r_aladdin - XXX) * numpy.sin(ANGLE_array)

            XXX = numpy.cos(ANGLE_array) * XXX + r_aladdin * (1.0 - numpy.cos(ANGLE_array))

            rays[:, 0] = XXX
            rays[:, 1] = YYY
            rays[:, 2] = ZZZ

        # angular coordinates
        if do_loop:
            for itik in range(NRAYS):
                ANGLE = ANGLE_array[itik]
                E_BEAM1 = E_BEAM1_array[itik]
                E_BEAM3 = E_BEAM3_array[itik]

                # ! C Synchrotron source
                # ! C Note. The angle of emission IN PLANE is the same as the one used
                # ! C before. This will give rise to a source curved along the orbit.
                # ! C The elevation angle is instead characteristic of the SR distribution.
                # ! C The electron beam emittance is included at this stage. Note that if
                # ! C EPSI = 0, we'll have E_BEAM = 0.0, with no changes.
                ANGLEX = ANGLE + E_BEAM1
                DIREC1 = numpy.tan(ANGLEX)
                if r_aladdin < 0:
                    DIREC1 *= -1.0
                DIREC2 = 1.0
                ARG_ANG = numpy.random.random()

                # ! C In the case of SR, we take into account the fact that the electron
                # ! C trajectory is not orthogonal to the field. This will give a correction
                # ! C to the photon energy.  We can write it as a correction to the
                # ! C magnetic field strength; this will linearly shift the critical energy
                # ! C and, with it, the energy of the emitted photon.
                E_TEMP3 = numpy.tan(E_BEAM3) / numpy.cos(E_BEAM1)
                E_TEMP2 = 1.0
                E_TEMP1 = numpy.tan(E_BEAM1)
                E_TEMP_MOD = numpy.sqrt(E_TEMP1**2 + E_TEMP2**2 + E_TEMP3**2)
                E_TEMP3 /= E_TEMP_MOD
                E_TEMP2 /= E_TEMP_MOD
                E_TEMP1 /= E_TEMP_MOD

                # interpolate for the photon energy and vertical angle.
                wavelength = codata.h * codata.c / codata.e / sampled_photon_energy[itik]
                Q_WAVE = 2 * numpy.pi / (wavelength*1e2)
                ANGLEV = sampled_angle[itik]


                # POL_DEG = sampled_polarization[itik]
                # fm_s, fm_p = sync_f_sigma_and_pi(ANGLEV * gamma, sampled_photon_energy[itik] / critical_energy)
                # POL_DEG = numpy.sqrt(fm_s) / (numpy.sqrt(fm_s) + numpy.sqrt(fm_p))

                if ANGLEV < 0: I_CHANGE = -1
                ANGLEV += E_BEAM3

                anglev_sign[itik] = numpy.sign(ANGLEV)

                DIREC3 = numpy.tan(ANGLEV) / numpy.cos(ANGLEX)

                DIREC_MOD = numpy.sqrt(DIREC1**2 + DIREC2**2 + DIREC3**2)
                DIREC3 /= DIREC_MOD
                DIREC2 /= DIREC_MOD
                DIREC1 /= DIREC_MOD

                rays[itik,3] = DIREC1
                rays[itik,4] = DIREC2
                rays[itik,5] = DIREC3
        else:
            ANGLE = ANGLE_array
            E_BEAM1 = E_BEAM1_array
            E_BEAM3 = E_BEAM3_array

            # ! C Synchrotron source
            # ! C Note. The angle of emission IN PLANE is the same as the one used
            # ! C before. This will give rise to a source curved along the orbit.
            # ! C The elevation angle is instead characteristic of the SR distribution.
            # ! C The electron beam emittance is included at this stage. Note that if
            # ! C EPSI = 0, we'll have E_BEAM = 0.0, with no changes.
            ANGLEX = ANGLE + E_BEAM1
            DIREC1 = numpy.tan(ANGLEX)
            if r_aladdin < 0:
                DIREC1 *= -1.0
            DIREC2 = numpy.ones_like(DIREC1) # 1.0
            ARG_ANG = numpy.random.random(NRAYS)

            # ! C In the case of SR, we take into account the fact that the electron
            # ! C trajectory is not orthogonal to the field. This will give a correction
            # ! C to the photon energy.  We can write it as a correction to the
            # ! C magnetic field strength; this will linearly shift the critical energy
            # ! C and, with it, the energy of the emitted photon.
            E_TEMP3 = numpy.tan(E_BEAM3) / numpy.cos(E_BEAM1)
            E_TEMP2 = numpy.ones_like(E_TEMP3) # 1.0
            E_TEMP1 = numpy.tan(E_BEAM1)
            E_TEMP_MOD = numpy.sqrt(E_TEMP1 ** 2 + E_TEMP2 ** 2 + E_TEMP3 ** 2)
            E_TEMP3 /= E_TEMP_MOD
            E_TEMP2 /= E_TEMP_MOD
            E_TEMP1 /= E_TEMP_MOD

            # interpolate for the photon energy, vertical angle, and the degree of polarization.
            wavelength = codata.h * codata.c / codata.e / sampled_photon_energy
            Q_WAVE = 2 * numpy.pi / (wavelength * 1e2)
            ANGLEV = sampled_angle

            # POL_DEG = sampled_polarization

            ANGLEV = ANGLEV + E_BEAM3

            anglev_sign = numpy.sign(ANGLEV)

            DIREC3 = numpy.tan(ANGLEV) / numpy.cos(ANGLEX)

            DIREC_MOD = numpy.sqrt(DIREC1 ** 2 + DIREC2 ** 2 + DIREC3 ** 2)
            DIREC3 /= DIREC_MOD
            DIREC2 /= DIREC_MOD
            DIREC1 /= DIREC_MOD

            rays[:, 3] = DIREC1
            rays[:, 4] = DIREC2
            rays[:, 5] = DIREC3

        t6 = time.time()
        #
        # electric field vectors (cols 7-9, 16-18) and phases (cols 14-15)
        #
        # ! C                 POLARIZATION
        # ! C   Generates the polarization of the ray. This is defined on the
        # ! C   source plane, so that A_VEC is along the X-axis and AP_VEC is along Z-axis.
        # ! C   Then care must be taken so that A will be perpendicular to the ray
        # ! C   direction.

        DIREC = rays[:,3:6].copy()
        A_VEC = numpy.zeros_like(DIREC)
        A_VEC[:,0] = 1.0

        # ! C   Rotate A_VEC so that it will be perpendicular to DIREC and with the
        # ! C   right components on the plane.
        A_TEMP = vector_cross(A_VEC, DIREC)
        A_VEC = vector_cross(DIREC, A_TEMP)
        A_VEC = vector_norm(A_VEC)
        AP_VEC = vector_cross(A_VEC, DIREC)
        AP_VEC = vector_norm(AP_VEC)

        #
        # obtain polarization for each ray (interpolation)
        #
        fm_s, fm_p = sync_f_sigma_and_pi(rays[:,5] * gamma, sampled_photon_energy / critical_energy)
        POL_DEG = numpy.sqrt(fm_s) / (numpy.sqrt(fm_s) + numpy.sqrt(fm_p))
        DENOM = numpy.sqrt(1.0 - 2.0 * POL_DEG + 2.0 * POL_DEG**2)
        AX = POL_DEG/DENOM
        for i in range(3):
            A_VEC[:, i] *= AX

        AZ = (1.0 - POL_DEG) / DENOM
        for i in range(3):
            AP_VEC[:, i] *= AZ

        rays[:, 6:9] =  A_VEC
        rays[:, 15:18] = AP_VEC


        # set flag (col 10)
        rays[:, 9] = 1.0

        #
        # photon energy (col 11)
        #
        wavelength = codata.h * codata.c / codata.e /sampled_photon_energy
        Q_WAVE = 2 * numpy.pi / (wavelength*1e2)
        rays[:, 10] =  Q_WAVE # sampled_photon_energy * A2EV

        # col 12 (ray index)
        rays[:, 11] =  1 + numpy.arange(NRAYS)

        # col 13 (optical path)
        rays[:, 12] = 0.0

        # ! C Now the phases of A_VEC and AP_VEC.

        if F_COHER == 1:
            PHASEX = 0.0
        else:
            PHASEX = numpy.random.random(NRAYS) * 2 * numpy.pi

        PHASEZ = PHASEX + POL_ANGLE * anglev_sign

        rays[:, 13] = PHASEX
        rays[:, 14] = PHASEZ
        t7 = time.time()

        if verbose:
            print("------------ timing---------")

            t = t7-t0
            print("            Total: ", t)
            print("            Pre (t1-t0)  ",    (t1-t0), 100 * (t1-t0) / t)
            print("            Pre (t2-t1): ",    (t2-t1), 100 * (t2-t1) / t)
            print("            Pre (t3-t2): ",    (t3-t2), 100 * (t3-t2) / t)
            print("            Pre (t4-t3): ",    (t4-t3), 100 * (t4-t3) / t)
            print("            emittance:  ",     (t5-t4), 100 * (t5-t4) / t)
            print("            coor+diverg: ",    (t6-t5), 100 * (t6-t5) / t)
            print("            Post: ",           (t7-t6), 100 * (t7-t6) / t)

        return rays



    def calculate_spectrum(self, output_file=""):
        """
                Calculates the spectrum.

        Parameters
        ----------
        output_file : str, optional
            Name of the file to write the spectrom (use "" for not writing file).

        Returns
        -------
        tuple
            (e, f, w) numpy arrays with photon energy (e), photon flux (f) and spectral power (w).
        """
        try:
            bm = self.get_magnetic_structure()
            ring = self.get_electron_beam()

            if bm.is_monochromatic():
                photon_energy_array = numpy.array([bm._EMIN])
            else:
                photon_energy_array = numpy.linspace(bm._EMIN,
                                                     bm._EMAX,
                                                     bm._NG_E)

            hdiv = numpy.abs(bm.length() / bm.radius())

            ec_ev = bm.get_critical_energy(ring.energy())
            i_a = ring.current()
            e_gev = ring.energy()

            f = sync_ene(0,
                         photon_energy_array,
                         ec_ev=ec_ev,
                         polarization=0,
                         e_gev=e_gev,
                         i_a=i_a,
                         hdiv_mrad=1e3*hdiv,
                         psi_min=0.0, # not needed
                         # psi_max=0.0, # not needed
                         # psi_npoints=1, # not needed
                        )
            return photon_energy_array, f, f * codata.e * 1e3
        except:
            raise Exception("Cannot compute spectrum")

    def to_python_code(self, **kwargs):
        """
        returns the python code for calculating the bending magnet source.

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

        script += "\n\n\n#light source\nfrom shadow4.sources.bending_magnet.s4_bending_magnet_light_source import S4BendingMagnetLightSource"
        script += "\nlight_source = S4BendingMagnetLightSource(name='%s', electron_beam=electron_beam, magnetic_structure=source, nrays=%d, seed=%s)" % \
                                                          (self.get_name(),self.get_nrays(),self.get_seed())

        script += "\nbeam = light_source.get_beam()"
        return script


if __name__ == "__main__":
    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam

    electron_beam = S4ElectronBeam(energy_in_GeV=6.04, energy_spread=0.001, current=0.2)
    electron_beam.set_sigmas_all(sigma_x=7.8e-05, sigma_y=4.87179e-05, sigma_xp=3.6e-05, sigma_yp=1.05556e-06)

    # magnetic structure
    # from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet

    source = S4BendingMagnet(
        radius=25.18408918746048,
        # from syned BM, can be obtained as S4BendingMagnet.calculate_magnetic_radius(0.8, electron_beam.energy())
        magnetic_field=0.8,  # from syned BM
        length=0.02518408918746048,  # from syned BM = abs(BM divergence * magnetic_field)
        emin=100.0,  # Photon energy scan from energy (in eV)
        emax=100000,  # Photon energy scan to energy (in eV)
        ng_e=500,  # Photon energy scan number of points
        flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
        )

    # light source
    from shadow4.sources.bending_magnet.s4_bending_magnet_light_source import S4BendingMagnetLightSource

    light_source = S4BendingMagnetLightSource(name='BendingMagnet',
                                              electron_beam=electron_beam,
                                              magnetic_structure=source,
                                              nrays=25000,
                                              seed=5676561)

    e, f, w = light_source.calculate_spectrum()

    from srxraylib.plot.gol import plot
    plot(e, f, xlog=1, ylog=1, marker='+')
    plot(e, w, xlog=1, ylog=1, marker='+')
