__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "12/01/2017"

import json
import os
import numpy

import scipy.constants as codata

# import Shadow
from SourceUndulatorFactory import undul_cdf, undul_phot, undul_phot_srw,  undul_phot_pysru

from SourceUndulatorInputOutput import write_file_undul_cdf, load_file_undul_cdf
from SourceUndulatorInputOutput import write_file_undul_phot_h5, write_file_undul_cdf_h5
# from SourceUndulatorInputOutput import load_file_undul_phot,write_file_undul_phot
# from SourceUndulatorInputOutput import load_file_undul_cdf,write_file_undul_sha
from syned.storage_ring.magnetic_structures.undulator import Undulator
from syned.storage_ring.electron_beam import ElectronBeam

import platform


class SampleUndulator(object):
    def __init__(self,name="",
                 syned_electron_beam=ElectronBeam(),
                 syned_undulator=Undulator(),
                 FLAG_EMITTANCE=0,EMIN=10000.0,EMAX=11000.0,NG_E=11,MAXANGLE=0.5,NG__T=31,NG_P=21,NG_J=20,SEED=36255655452,NRAYS=5000):

        # # Machine
        self.syned_electron_beam = syned_electron_beam

        # # Undulator
        self.syned_undulator = syned_undulator

        self.FLAG_EMITTANCE  =  FLAG_EMITTANCE# Yes  # Use emittance (0=No, 1=Yes)
        # Photon energy scan
        self.EMIN            = EMIN   # Photon energy scan from energy (in eV)
        self.EMAX            = EMAX   # Photon energy scan to energy (in eV)
        self.NG_E            = NG_E        # Photon energy scan number of points
        # Geometry
        self.MAXANGLE        = MAXANGLE      # Maximum radiation semiaperture in mrad # TODO: define it in rad, for consistency
        self.NG_T            = NG__T       # Number of points in angle theta
        self.NG_P            = NG_P       # Number of points in angle phi
        self.NG_J            = NG_J       # Number of points in electron trajectory (per period)
        # ray tracing
        self.SEED            = SEED  # Random seed
        self.NRAYS           = NRAYS         # Number of rays


        self.result_radiation = None
        self.result_cdf = None

        # if platform.system() == "Linux":
        #     self.SHADOW3_BINARY = "/users/srio/OASYS1.1/shadow3/shadow3"
        # else:
        #     self.SHADOW3_BINARY = "/Users/srio/Oasys/OASYS1.1/shadow3/shadow3"


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

    def calculate_radiation(self,code_undul_phot='internal'):
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


        # undul_phot
        if code_undul_phot == 'internal':
            undul_phot_dict = undul_phot(E_ENERGY  = self.syned_electron_beam.energy(),
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

        elif code_undul_phot == 'pysru':
            undul_phot_dict = undul_phot_pysru(E_ENERGY  = self.syned_electron_beam.energy(),
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
        elif code_undul_phot == 'srw':
            undul_phot_dict = undul_phot_srw(E_ENERGY  = self.syned_electron_beam.energy(),
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
        undul_phot_dict["code_undul_phot"] = code_undul_phot
        undul_phot_dict["info"] = self.info()

        self.result_radiation = undul_phot_dict
        return undul_phot_dict


    def calculate_cdf(self):
                    #,code_undul_phot='internal'):
                    #,use_existing_undul_phot_output=None,
                    #dump_undul_phot_file=False):
        """
        calculates cdf (cumulative distribution functions)

        it first computes the radiation (calculate_radiation) and the integrate using
        SourceUndulatorFactory.undul_cdf

        :param code_undul_phot: code_undul_phot: 'internal' (calls undul_phot), 'pysru'
                (calls undul_phot_pysru) or 'srw' (calls undul_phot_srw)
        :param use_existing_undul_phot_output: set to a file name or dictionary to use this
            particular output from undul_phot
        :param dump_undul_phot_file: if True writes the undul_phot output in uphot.dat file
        :return: a dictionary with the undul_cdf output. It also dumps this output to a file
            called xshundul.sha file.
        """

        #
        # undul_phot
        #

        # undul_phot_dict = self.calculate_radiation(code_undul_phot=code_undul_phot)

        # if use_existing_undul_phot_output is None:
        #     undul_phot_dict = self.calculate_radiation(code_undul_phot=code_undul_phot)
        # else:
        #     if isinstance(use_existing_undul_phot_output,str):
        #         undul_phot_dict = load_file_undul_phot(use_existing_undul_phot_output)
        #     elif isinstance(use_existing_undul_phot_output,dict):
        #         undul_phot_dict = use_existing_undul_phot_output
        #         #TODO: import parameters from external file E_MIN, E_MAX, MAXANGLE, N_*
        #     else:
        #         raise Exception("Bad undul_phot data.")
        #
        #
        #
        # if dump_undul_phot_file:
        #     write_file_undul_phot(undul_phot_dict,file_out="uphot.dat")



        # undul_cdf

        undul_cdf_dict = undul_cdf(self.result_radiation,method='trapz')
        # write_file_undul_sha(undul_cdf_dict,file_out="xshundul.sha",)

        self.result_cdf = undul_cdf_dict
        return undul_cdf_dict

    def load_from_file_cdf(self,file_in="xshundul.sha"):
        self.result_radiation = None
        self.result_radiation = None
        self.result_radiation = load_file_undul_cdf(file_in)

    def write_file_h5(self,file_out):
        write_file_undul_phot_h5(self.result_radiation,file_out=file_out,mode='w',entry_name="radiation")
        write_file_undul_cdf_h5(self.result_cdf,file_out=file_out,mode='a',entry_name="cdf")




    # TODO: remove in far future
    def get_shadow3_source_object(self,m_to_user_unit=1e2):
        """

        creates a Shadow.Source object with the undulator parameters inside (proprocessor file: xshundul.sha)

        :return:
        """
        # initialize shadow3 source (oe0) and beam

        import Shadow

        oe0 = Shadow.Source()

        if self.FLAG_EMITTANCE:
            sigmas = self.syned_electron_beam.get_sigmas_all()
            oe0.EPSI_X = m_to_user_unit*sigmas[0]*sigmas[2]
            oe0.EPSI_Z = m_to_user_unit*sigmas[1]*sigmas[3]
            oe0.SIGDIX = 0.0
            oe0.SIGDIZ = 0.0
            oe0.SIGMAX = m_to_user_unit*sigmas[0]
            oe0.SIGMAY = 0.0
            oe0.SIGMAZ = m_to_user_unit*sigmas[1]
        else:
            oe0.EPSI_X = 0.0
            oe0.EPSI_Z = 0.0
            oe0.SIGDIX = 0.0
            oe0.SIGDIZ = 0.0
            oe0.SIGMAX = 0.0
            oe0.SIGMAY = 0.0
            oe0.SIGMAZ = 0.0

        oe0.FILE_TRAJ = b'xshundul.sha'
        oe0.ISTAR1 = self.SEED
        oe0.NPOINT = self.NRAYS
        oe0.F_WIGGLER = 2

        return oe0


    def calculate_shadow3_beam(self,dump_start_files=False):
        import Shadow

        # # using existing run
        # oe0 = Shadow.Source()
        # oe0.load("start.00")
        # beam = Shadow.Beam()
        # beam.genSource(oe0)


        # # creating start.00
        # oe0 = self.get_shadow3_source_object(m_to_user_unit=1e2)
        # oe0.write("start.00")
        # beam = Shadow.Beam()
        # beam.genSource(oe0)
        # beam.write("begin.dat")


        # # creating preprocessor file xshundul.sha
        # self.calculate_radiation(code_undul_phot='srw')
        # self.calculate_cdf()
        # write_file_undul_cdf(self.result_cdf,file_out="xshundul.sha")
        # # initialize shadow3 source (oe0) and beam
        # oe0 = self.get_shadow3_source_object(m_to_user_unit=1e2)
        # oe0.write("start.00")
        # beam = Shadow.Beam()
        # beam.genSource(oe0)
        # beam.write("begin.dat")


        # loading preprocessor file xshundul.sha
        cdf = load_file_undul_cdf("xshundul.sha")
        self.result_radiation = None
        self.result_cdf = cdf
        # initialize shadow3 source (oe0) and beam
        beam = self.sample() # Shadow.Beam(N=self.NRAYS)



        return beam

    def sample(self):

        import Shadow
        from inverse_method_sampler import Sampler2D, Sampler3D
        from srxraylib.plot.gol import plot_scatter, plot_image, plot

        beam = Shadow.Beam(N=self.NRAYS)

        #
        # sample sizes
        #
        sigmas = self.syned_electron_beam.get_sigmas_all()

        beam.rays[:,9] = 1.0



        #
        # sample divergences
        #
        self.calculate_radiation()
        for k in self.result_radiation.keys():
            print("radiation>>>>",k)
        # self.calculate_cdf()
        # for k in self.result_cdf.keys():
        #     print("cdf>>>>",k)
        # self.write_file_h5("tmp.h5")


        theta = self.result_radiation["theta"]
        phi = self.result_radiation["phi"]
        photon_energy = self.result_radiation["photon_energy"]


        photon_energy_spectrum = 'polychromatic' # 'monochromatic' #
        #
        # monochromatic case
        #
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
            # plot_image(tmp_theta,theta,phi,aspect='auto')

            s3d = Sampler3D(tmp,photon_energy,theta,phi)

            sampled_photon_energy,sampled_theta,sampled_phi = s3d.get_n_sampled_points(self.NRAYS)
            # print(sampled_photon_energy)


        #
        # the Shadow way
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
        # photon energy
        #

        # xx = self.result_radiation["photon_energy"]
        # yy = self.result_radiation["radiation"].sum(axis=2).sum(axis=1)
        #
        #
        # print(xx.shape,yy.shape)
        # s3d._cdf_calculate()
        # plot(xx,yy,title="energy",show=0)
        # plot(xx,s3d._cdf1,title="energy cdf")
        # plot_image(s3d._cdf2)

        A2EV = 2.0*numpy.pi/(codata.h*codata.c/codata.e*1e2)
        beam.rays[:,10] =  sampled_photon_energy * A2EV


        #
        # electric vectors
        #

        beam.rays[:,6] =  1.0



        #
        # # plot_image(s2d.pdf(),cmap='binary',title="pdf")
        #
        # cdf2,cdf1 = s2d.cdf()
        # plot_image(cdf2,cmap='binary',title="cdf")
        # # plot(s2d.abscissas()[0],s2d.cdf()[0][:,-1])
        # plot(s2d.abscissas()[0],cdf1)
        #
        # x0s,x1s = s2d.get_n_sampled_points(100000)
        # plot_scatter(x0s,x1s)


        #
        # write output
        #

        beam.write("begin.dat")



        return beam


    # def to_dictionary(self):
    #     """
    #     returns a python dictionary of the Shadow.Source instance
    #     :return: a dictionary
    #     """
    #     h = {}
    #     h["E_ENERGY"]        = self.E_ENERGY
    #     h["INTENSITY"]       = self.INTENSITY
    #     h["SX"]              = self.SX
    #     h["SZ"]              = self.SZ
    #     h["SXP"]             = self.SXP
    #     h["SZP"]             = self.SZP
    #     h["FLAG_EMITTANCE"]  = self.FLAG_EMITTANCE
    #     h["LAMBDAU"]         = self.LAMBDAU
    #     h["NPERIODS"]        = self.NPERIODS
    #     h["K"]               = self.K
    #     h["EMIN"]            = self.EMIN
    #     h["EMAX"]            = self.EMAX
    #     h["NG_E"]            = self.NG_E
    #     h["MAXANGLE"]        = self.MAXANGLE
    #     h["NG_T"]            = self.NG_T
    #     h["NG_P"]            = self.NG_P
    #     h["N_J"]             = self.N_J
    #     h["SEED"]            = self.SEED
    #     h["NRAYS"]           = self.NRAYS
    #
    #     return h
    #
    # def duplicate(self):
    #     """
    #     makes a copy of the source
    #     :return: new instance of Shadow.Source()
    #     """
    #     out = SourceUndulator()
    #     out.set_from_dictionary(self.to_dictionary())
    #     return out
    #
    # def set_from_keywords(self,E_ENERGY=6.04,INTENSITY=0.2,
    #             SX=0.04,SZ=0.001,SXP=10e-6,SZP=4e-6,FLAG_EMITTANCE=1,
    #             LAMBDAU=0.032,NPERIODS=50,K=0.25,
    #             EMIN= 10498.0000,EMAX= 10499.0000,NG_E=101,MAXANGLE=0.1,NG_T=51,NG_P=11,N_J=10,
    #             SEED=36255,NRAYS=15000,):
    #     """
    #     Sets the undulator parameters using keywords:
    #
    #     :param E_ENERGY: Electron energy in GeV
    #     :param INTENSITY: Electron current in A
    #     :param SX: Electron sigma X
    #     :param SZ: Electron sigma Z
    #     :param SXP: Electron sigma X'
    #     :param SZP: Electron sigma Z'
    #     :param FLAG_EMITTANCE: 0=No emittance (S* ignored), 1=Use emittance
    #     :param LAMBDAU: Undulator period in m
    #     :param NPERIODS: Undulator number of periodas
    #     :param K: Undulator K-value
    #     :param EMIN: Photon energy scan from energy (in eV)
    #     :param EMAX: Photon energy scan to energy (in eV)
    #     :param NG_E: Photon energy scan number of points
    #     :param MAXANGLE: Maximum radiation semiaperture in mrad # TODO: define it in rad, for consistency
    #     :param NG_T: Number of points in angle theta
    #     :param NG_P: Number of points in angle phi
    #     :param N_J: Number of points in electron trajectory (per period)
    #     :param SEED: Random seed
    #     :param NRAYS: Number of rays
    #     :return:
    #     """
    #     self.E_ENERGY          = E_ENERGY
    #     self.INTENSITY         = INTENSITY
    #     self.SX                = SX
    #     self.SZ                = SZ
    #     self.SXP               = SXP
    #     self.SZP               = SZP
    #     self.FLAG_EMITTANCE    = FLAG_EMITTANCE
    #     self.LAMBDAU           = LAMBDAU
    #     self.NPERIODS          = NPERIODS
    #     self.K                 = K
    #     self.EMIN              = EMIN
    #     self.EMAX              = EMAX
    #     self.NG_E              = NG_E
    #     self.MAXANGLE          = MAXANGLE
    #     self.NG_T              = NG_T
    #     self.NG_P              = NG_P
    #     self.N_J               = N_J
    #     self.SEED              = SEED
    #     self.NRAYS             = NRAYS
    #



    #
    #
    # def set_from_dictionary(self,h):
    #     """
    #     set undulator variables from a dictionary. The following keys must be present:
    #       E_ENERGY
    #       INTENSITY
    #       SX
    #       SZ
    #       SXP
    #       SZP
    #       FLAG_EMITTANCE
    #       LAMBDAU
    #       NPERIODS
    #       K
    #       EMIN
    #       EMAX
    #       NG_E
    #       MAXANGLE
    #       NG_T
    #       NG_P
    #       N_J
    #       SEED
    #       NRAYS
    #
    #     :param h:
    #     :return:
    #     """
    #     self.E_ENERGY        = h["E_ENERGY"]
    #     self.INTENSITY       = h["INTENSITY"]
    #     self.SX              = h["SX"]
    #     self.SZ              = h["SZ"]
    #     self.SXP             = h["SXP"]
    #     self.SZP             = h["SZP"]
    #     self.FLAG_EMITTANCE  = h["FLAG_EMITTANCE"]
    #     self.LAMBDAU         = h["LAMBDAU"]
    #     self.NPERIODS        = h["NPERIODS"]
    #     self.K               = h["K"]
    #     self.EMIN            = h["EMIN"]
    #     self.EMAX            = h["EMAX"]
    #     self.NG_E            = h["NG_E"]
    #     self.MAXANGLE        = h["MAXANGLE"]
    #     self.NG_T            = h["NG_T"]
    #     self.NG_P            = h["NG_P"]
    #     self.N_J             = h["N_J"]
    #     self.SEED            = h["SEED"]
    #     self.NRAYS           = h["NRAYS"]
    #
    #
    # #
    # # useful getters
    # #
    # def get_lorentz_factor(self):
    #     """
    #     Calculates the Lorentz factor (gamma) from electron energy
    #     :return:
    #     """
    #     return ( self.E_ENERGY * 1e9) / (codata.m_e *  codata.c**2 / codata.e)
    #
    # def get_resonance_wavelength(self, harmonic_number=1, theta_x=0.0, theta_z=0.0):
    #     """
    #     sets photon energy to monochromatic line at resonance
    #     :param harmonic_number:
    #     :param theta_x:
    #     :param theta_z:
    #     :return:
    #     """
    #     gamma = self.get_lorentz_factor()
    #     wavelength = self.LAMBDAU / (2.0*gamma **2) * (1 + self.K**2 / 2.0 + gamma**2 * (theta_x**2 + theta_z ** 2))
    #     return wavelength / harmonic_number
    #
    # def get_resonance_frequency(self, harmonic_number=1, theta_x=0, theta_z=0):
    #     """
    #     gets the resonance wavelength (in m)
    #     :param harmonic_number:
    #     :param theta_x:
    #     :param theta_z:
    #     :return:
    #     """
    #     return codata.c / self.get_resonance_wavelength(harmonic_number,theta_x,theta_z)
    #
    # def get_resonance_energy(self, harmonic_number=1, theta_x=0, theta_z=0):
    #     """
    #     gets the resonance photon energy (in eV)
    #     :param harmonic_number:
    #     :param theta_x:
    #     :param theta_z:
    #     :return:
    #     """
    #     return codata.h*codata.c/codata.e*1e10 / (1e10*self.get_resonance_wavelength(harmonic_number, theta_x, theta_z))
    #
    # def get_resonance_central_cone(self,harmonic_number=1):
    #     """
    #     gets the angular aperture of the central cone (in rad).
    #     It is calculated as:
    #         CC = (1/gamma) * sqrt(  (1+0.5*K^2) / ( 2 N n) )
    #         with N = number of periods, n = harmonic number
    #
    #         Note that the standard deviation is ~ 0.68*CC and the FWHM is 2.35*0.68*CC
    #
    #     :param harmonic_number:
    #     :return:
    #     """
    #     return 1.0/self.get_lorentz_factor()*numpy.sqrt( (1+0.5*self.K**2)/(2*self.NPERIODS*harmonic_number) )
    #

    #

    #
    #
    # def set_from_dictionary(self,h):
    #     """
    #     set undulator variables from a dictionary. The following keys must be present:
    #       E_ENERGY
    #       INTENSITY
    #       SX
    #       SZ
    #       SXP
    #       SZP
    #       FLAG_EMITTANCE
    #       LAMBDAU
    #       NPERIODS
    #       K
    #       EMIN
    #       EMAX
    #       NG_E
    #       MAXANGLE
    #       NG_T
    #       NG_P
    #       N_J
    #       SEED
    #       NRAYS
    #
    #     :param h:
    #     :return:
    #     """
    #     self.E_ENERGY        = h["E_ENERGY"]
    #     self.INTENSITY       = h["INTENSITY"]
    #     self.SX              = h["SX"]
    #     self.SZ              = h["SZ"]
    #     self.SXP             = h["SXP"]
    #     self.SZP             = h["SZP"]
    #     self.FLAG_EMITTANCE  = h["FLAG_EMITTANCE"]
    #     self.LAMBDAU         = h["LAMBDAU"]
    #     self.NPERIODS        = h["NPERIODS"]
    #     self.K               = h["K"]
    #     self.EMIN            = h["EMIN"]
    #     self.EMAX            = h["EMAX"]
    #     self.NG_E            = h["NG_E"]
    #     self.MAXANGLE        = h["MAXANGLE"]
    #     self.NG_T            = h["NG_T"]
    #     self.NG_P            = h["NG_P"]
    #     self.N_J             = h["N_J"]
    #     self.SEED            = h["SEED"]
    #     self.NRAYS           = h["NRAYS"]
    #
    # #
    # # load/writers
    # #
    # def load_json_shadowvui_file(self,inFileTxt):
    #     """
    #     sets undulator parameters from a json file from shadowVUI and
    #
    #     :param inFileTxt:
    #     :return:
    #     """
    #     #TODO: remove this routine in the future
    #     #
    #     # read inputs from a file created by ShadowVUI ----------------------------
    #     #
    #
    #     with open(inFileTxt, mode='r') as f1:
    #         h = json.load(f1)
    #     self.load_json_shadowvui_dictionary(h)
    #
    #
    # def load_json_shadowvui_dictionary(self,h):
    #     """
    #     sets undulator parameters from a dictionary loaded from a shadowVUI json file
    #     :param h:
    #     :return:
    #     """
    #     # self.dict = {
    #     #     "LAMBDAU":     0.0320000015,
    #     #     "K":      0.250000000,
    #     #     "E_ENERGY":       6.03999996,
    #     #     "E_ENERGY_SPREAD":    0.00100000005,
    #     #     "NPERIODS": 50,
    #     #     "EMIN":       10498.0000,
    #     #     "EMAX":       10499.0000,
    #     #     "INTENSITY":      0.200000003,
    #     #     "MAXANGLE":      0.100000001,
    #     #     "NG_E": 101,
    #     #     "NG_T": 51,
    #     #     "NG_P": 11,
    #     #     "NG_PLOT(1)":"1",
    #     #     "NG_PLOT(2)":"No",
    #     #     "NG_PLOT(3)":"Yes",
    #     #     "UNDUL_PHOT_FLAG(1)":"4",
    #     #     "UNDUL_PHOT_FLAG(2)":"Shadow code",
    #     #     "UNDUL_PHOT_FLAG(3)":"Urgent code",
    #     #     "UNDUL_PHOT_FLAG(4)":"SRW code",
    #     #     "UNDUL_PHOT_FLAG(5)":"Gaussian Approximation",
    #     #     "UNDUL_PHOT_FLAG(6)":"ESRF python code",
    #     #     "SEED": 36255,
    #     #     "SX":     0.0399999991,
    #     #     "SZ":    0.00100000005,
    #     #     "EX":   4.00000005E-07,
    #     #     "EZ":   3.99999989E-09,
    #     #     "FLAG_EMITTANCE(1)":"0",
    #     #     "FLAG_EMITTANCE(2)":"No",
    #     #     "FLAG_EMITTANCE(3)":"Yes",
    #     #     "NRAYS": 15000,
    #     #     "F_BOUND_SOUR": 0,
    #     #     "FILE_BOUND":"NONESPECIFIED",
    #     #     "SLIT_DISTANCE":       1000.00000,
    #     #     "SLIT_XMIN":      -1.00000000,
    #     #     "SLIT_XMAX":       1.00000000,
    #     #     "SLIT_ZMIN":      -1.00000000,
    #     #     "SLIT_ZMAX":       1.00000000,
    #     #     "NTOTALPOINT": 10000000,
    #     #     "JUNK4JSON":0
    #     #     }
    #
    #     self.E_ENERGY        = h["E_ENERGY"]
    #     self.N_J             = 20
    #     self.INTENSITY       = h["INTENSITY"]
    #     self.SX              = h["SX"]
    #     self.SZ              = h["SZ"]
    #     self.SXP             = h["EX"] / h["SX"]
    #     self.SZP             = h["EZ"] / h["SZ"]
    #     self.FLAG_EMITTANCE  = int(h["FLAG_EMITTANCE(1)"])
    #     self.LAMBDAU         = h["LAMBDAU"]
    #     self.NPERIODS        = h["NPERIODS"]
    #     self.K               = h["K"]
    #     self.EMIN            = h["EMIN"]
    #     self.EMAX            = h["EMAX"]
    #     self.NG_E            = h["NG_E"]
    #     self.MAXANGLE        = h["MAXANGLE"]
    #     self.NG_T            = h["NG_T"]
    #     self.NG_P            = h["NG_P"]
    #     self.SEED            = h["SEED"]
    #     self.NRAYS           = h["NRAYS"]
    #
    #
    # def load(self,inFileTxt):
    #     """
    #     load undulator parameters from a json file (type startj.00)
    #     :param inFileTxt:
    #     :return:
    #     """
    #     with open(inFileTxt, mode='r') as f1:
    #         h = json.load(f1)
    #     self.set_from_dictionary(h)
    #
    # def write(self,file_out='startj.00'):
    #     """
    #     write undulator parameters to a json file (type startj.00)
    #     :param file_out:
    #     :return:
    #     """
    #     data = self.to_dictionary()
    #     with open(file_out, 'w') as outfile:
    #         json.dump(data, outfile, indent=4, sort_keys=True, separators=(',', ':'))
    #
    # #
    # # util
    # #
    #
    # def sourcinfo(self,title=None):
    #     '''
    #     mimics SHADOW sourcinfo postprocessor. Returns a text array.
    #     :return: a text string
    #     '''
    #
    #     txt = ''
    #     TOPLIN = '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
    #
    #
    #     txt += TOPLIN
    #     txt += '**************  S O U R C E       D E S C R I P T I O N  **************\n'
    #     if title == None:
    #         txt += '\n\n'
    #     else:
    #         txt += title+'\n'
    #     txt += TOPLIN
    #
    #     # TODO improve look
    #     txt += self.info()
    #
    #     txt += TOPLIN
    #     txt += '***************                 E N D                  ***************\n'
    #     txt += TOPLIN
    #     return (txt)
    #
    #

    #
    # def _shadow3_commands(self,commands="exit\n",input_file="shadow3_tmp.inp"):
    #     # for internal use
    #     f = open(input_file,'w')
    #     f.write(commands)
    #     f.close()
    #     os.system(self.SHADOW3_BINARY+" < "+input_file)
    #
    # #
    # # call calculators
    # #
    # def calculate_shadow3_beam_using_preprocessors(self):
    #     """
    #     Calculates rays (source) using preprocessors in binary shadow3 (for
    #         comparison purposes, from python use calculate_beam() instead)
    #
    #     :return:  a Shadow.Beam object
    #     """
    #     jsn = self.to_dictionary()
    #     # jsn = self.dict
    #     # print(self.info())
    #     os.system("rm -f start.00 systemfile.dat begin.dat xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")
    #
    #     # epath
    #     commands = "epath\n2\n%f \n%f \n%f \n%d \n1.\nxshundul.par\nxshundul.traj\n1\nxshundul.plt\nexit\n"% \
    #     (jsn["LAMBDAU"],jsn["K"],jsn["E_ENERGY"],101)
    #     self._shadow3_commands(commands=commands,input_file="shadow3_epath.inp")
    #
    #     # undul_set
    #     NG_E = jsn["NG_E"]
    #     # TODO: is seems a bug in shadow3: undul_set must be NG_E>1 (otherwise nosense)
    #     # but if emin=emaxthe  resulting uphot.nml has NG_E=1
    #     if NG_E == 1: NG_E = 2
    #
    #     commands = "undul_set\n0\n0\n%d \n%d \n%d \nxshundul.traj\n%d\n%f\n%f\n%f\n%f\n0\n1000\nexit\n"% \
    #         (NG_E,jsn["NG_T"],jsn["NG_P"],
    #          jsn["NPERIODS"],jsn["EMIN"],jsn["EMAX"],jsn["INTENSITY"],jsn["MAXANGLE"])
    #     self._shadow3_commands(commands=commands,input_file="shadow3_undul_set.inp")
    #
    #     # undul_phot
    #     self._shadow3_commands(commands="undul_phot\nexit\n",input_file="shadow3_undul_phot.inp")
    #     self._shadow3_commands(commands="undul_phot_dump\nexit\n",input_file="shadow3_undul_phot_dump.inp")
    #
    #     # undul_cdf
    #     self._shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
    #                     input_file="shadow3_undul_cdf.inp")
    #
    #
    #     # input source
    #     if self.FLAG_EMITTANCE:
    #         commands = "input_source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
    #         (jsn["NRAYS"],jsn["SEED"],jsn["SX"],jsn["SZ"],jsn["SX"]*jsn["SXP"],0,jsn["SZ"]*jsn["SZP"],0,3,1,1)
    #     else:
    #         commands = "input_source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
    #         (jsn["NRAYS"],jsn["SEED"],0,0,0,0,0,0,3,1,1)
    #
    #     self._shadow3_commands(commands=commands,input_file="shadow3_input_source.inp")
    #
    #     # run source
    #     commands = "source\nsystemfile\nexit\n"
    #     self._shadow3_commands(commands=commands,input_file="shadow3_source.inp")
    #
    #
    #     # return shadow3 beam
    #     beam = Shadow.Beam()
    #     beam.load("begin.dat")
    #     return beam
    #

    #
    #
    #


