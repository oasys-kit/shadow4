__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "12/01/2017"

import json
import os
import numpy

import scipy.constants as codata

import Shadow
from SourceUndulatorFactory import undul_cdf, undul_phot, undul_phot_srw,  undul_phot_pysru
from SourceUndulatorInputOutput import load_file_undul_phot,write_file_undul_phot
from SourceUndulatorInputOutput import load_file_undul_cdf,write_file_undul_sha

import platform


class SourceUndulator(object):
    def __init__(self):

        # Machine
        self.E_ENERGY        = 1.0   # Electron energy in GeV
        self.INTENSITY       = 1.0   # Electron current in A
        self.SX              = 1e-6  # Electron sigma X
        self.SZ              = 1e-6  # Electron sigma Z
        self.SXP             = 1e-6  # Electron sigma X'
        self.SZP             = 1e-6  # Electron sigma Z'
        # self.EX              = 4.00000005E-07
        # self.EZ              = 3.99999989E-09
        self.FLAG_EMITTANCE  = 1 # Yes  # Use emittance (0=No, 1=Yes)
        # Undulator
        self.LAMBDAU         = 0.01  # Undulator period in m
        self.NPERIODS        = 100   # Undulator number of periodas
        self.K               = 0.8   # Undulator K-value
        # Photon energy scan
        self.EMIN            = 10000.0   # Photon energy scan from energy (in eV)
        self.EMAX            = 11000.0   # Photon energy scan to energy (in eV)
        self.NG_E            = 11        # Photon energy scan number of points
        # Geometry
        self.MAXANGLE        = 0.5      # Maximum radiation semiaperture in mrad # TODO: define it in rad, for consistency
        self.NG_T            = 31       # Number of points in angle theta
        self.NG_P            = 21       # Number of points in angle phi
        self.NG_J            = 20       # Number of points in electron trajectory (per period)
        # ray tracing
        self.SEED            = 36255655452  # Random seed
        self.NRAYS           = 5000         # Number of rays


        # # Machine
        # self.E_ENERGY        = 6.04
        # self.E_ENERGY_SPREAD = 0.001
        # self.INTENSITY       = 0.2
        # self.SX              = 0.0399999991
        # self.SZ              = 0.00100000005
        # self.SXP             = 4.00000005E-07 / 0.0399999991
        # self.SZP             = 3.99999989E-09 / 0.00100000005
        # # self.EX              = 4.00000005E-07
        # # self.EZ              = 3.99999989E-09
        # self.FLAG_EMITTANCE  = 1 # Yes
        # # Undulator
        # self.LAMBDAU         = 0.0320000015
        # self.NPERIODS        = 50
        # self.K               = 0.250000000
        # # Photon energy scan
        # self.EMIN            = 10498.0000
        # self.EMAX            = 10499.0000
        # self.NG_E            = 101
        # # Geometry
        # self.MAXANGLE        = 0.1
        # self.NG_T            = 51
        # self.NG_P            = 11
        # # ray tracing
        # self.SEED            = 36255
        # self.NRAYS           = 15000


        if platform.system() == "Linux":
            self.SHADOW3_BINARY = "/users/srio/OASYS1.1/shadow3/shadow3"
        else:
            self.SHADOW3_BINARY = "/Users/srio/Oasys/OASYS1.1/shadow3/shadow3"

    def to_dictionary(self):
        """
        returns a python dictionary of the Shadow.Source instance
        :return: a dictionary
        """
        h = {}
        h["E_ENERGY"]        = self.E_ENERGY
        h["INTENSITY"]       = self.INTENSITY
        h["SX"]              = self.SX
        h["SZ"]              = self.SZ
        h["SXP"]             = self.SXP
        h["SZP"]             = self.SZP
        h["FLAG_EMITTANCE"]  = self.FLAG_EMITTANCE
        h["LAMBDAU"]         = self.LAMBDAU
        h["NPERIODS"]        = self.NPERIODS
        h["K"]               = self.K
        h["EMIN"]            = self.EMIN
        h["EMAX"]            = self.EMAX
        h["NG_E"]            = self.NG_E
        h["MAXANGLE"]        = self.MAXANGLE
        h["NG_T"]            = self.NG_T
        h["NG_P"]            = self.NG_P
        h["N_J"]             = self.N_J
        h["SEED"]            = self.SEED
        h["NRAYS"]           = self.NRAYS

        return h

    def duplicate(self):
        """
        makes a copy of the source
        :return: new instance of Shadow.Source()
        """
        out = SourceUndulator()
        out.set_from_dictionary(self.to_dictionary())
        return out

    def set_from_keywords(self,E_ENERGY=6.04,INTENSITY=0.2,
                SX=0.04,SZ=0.001,SXP=10e-6,SZP=4e-6,FLAG_EMITTANCE=1,
                LAMBDAU=0.032,NPERIODS=50,K=0.25,
                EMIN= 10498.0000,EMAX= 10499.0000,NG_E=101,MAXANGLE=0.1,NG_T=51,NG_P=11,N_J=10,
                SEED=36255,NRAYS=15000,):
        """
        Sets the undulator parameters using keywords:

        :param E_ENERGY: Electron energy in GeV
        :param INTENSITY: Electron current in A
        :param SX: Electron sigma X
        :param SZ: Electron sigma Z
        :param SXP: Electron sigma X'
        :param SZP: Electron sigma Z'
        :param FLAG_EMITTANCE: 0=No emittance (S* ignored), 1=Use emittance
        :param LAMBDAU: Undulator period in m
        :param NPERIODS: Undulator number of periodas
        :param K: Undulator K-value
        :param EMIN: Photon energy scan from energy (in eV)
        :param EMAX: Photon energy scan to energy (in eV)
        :param NG_E: Photon energy scan number of points
        :param MAXANGLE: Maximum radiation semiaperture in mrad # TODO: define it in rad, for consistency
        :param NG_T: Number of points in angle theta
        :param NG_P: Number of points in angle phi
        :param N_J: Number of points in electron trajectory (per period)
        :param SEED: Random seed
        :param NRAYS: Number of rays
        :return:
        """
        self.E_ENERGY          = E_ENERGY
        self.INTENSITY         = INTENSITY
        self.SX                = SX
        self.SZ                = SZ
        self.SXP               = SXP
        self.SZP               = SZP
        self.FLAG_EMITTANCE    = FLAG_EMITTANCE
        self.LAMBDAU           = LAMBDAU
        self.NPERIODS          = NPERIODS
        self.K                 = K
        self.EMIN              = EMIN
        self.EMAX              = EMAX
        self.NG_E              = NG_E
        self.MAXANGLE          = MAXANGLE
        self.NG_T              = NG_T
        self.NG_P              = NG_P
        self.N_J               = N_J
        self.SEED              = SEED
        self.NRAYS             = NRAYS

    def set_energy_monochromatic(self,emin):
        """
        Sets a single energy line for the source (monochromatic)
        :param emin: the energy in eV
        :return:
        """
        self.EMIN = emin
        self.EMAX = emin
        self.NG_E = 1

    def set_energy_monochromatic_at_resonance(self,harmonic_number):

        self.set_energy_monochromatic(self.get_resonance_energy(harmonic_number=harmonic_number))
        # take 3*sigma - MAXANGLE is in mrad!!
        self.MAXANGLE = 3 * 0.69 * 1e3 * self.get_resonance_central_cone(harmonic_number)

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


    def set_from_dictionary(self,h):
        """
        set undulator variables from a dictionary. The following keys must be present:
          E_ENERGY
          INTENSITY
          SX
          SZ
          SXP
          SZP
          FLAG_EMITTANCE
          LAMBDAU
          NPERIODS
          K
          EMIN
          EMAX
          NG_E
          MAXANGLE
          NG_T
          NG_P
          N_J
          SEED
          NRAYS

        :param h:
        :return:
        """
        self.E_ENERGY        = h["E_ENERGY"]
        self.INTENSITY       = h["INTENSITY"]
        self.SX              = h["SX"]
        self.SZ              = h["SZ"]
        self.SXP             = h["SXP"]
        self.SZP             = h["SZP"]
        self.FLAG_EMITTANCE  = h["FLAG_EMITTANCE"]
        self.LAMBDAU         = h["LAMBDAU"]
        self.NPERIODS        = h["NPERIODS"]
        self.K               = h["K"]
        self.EMIN            = h["EMIN"]
        self.EMAX            = h["EMAX"]
        self.NG_E            = h["NG_E"]
        self.MAXANGLE        = h["MAXANGLE"]
        self.NG_T            = h["NG_T"]
        self.NG_P            = h["NG_P"]
        self.N_J             = h["N_J"]
        self.SEED            = h["SEED"]
        self.NRAYS           = h["NRAYS"]


    #
    # useful getters
    #
    def get_lorentz_factor(self):
        """
        Calculates the Lorentz factor (gamma) from electron energy
        :return:
        """
        return ( self.E_ENERGY * 1e9) / (codata.m_e *  codata.c**2 / codata.e)

    def get_resonance_wavelength(self, harmonic_number=1, theta_x=0.0, theta_z=0.0):
        """
        sets photon energy to monochromatic line at resonance
        :param harmonic_number:
        :param theta_x:
        :param theta_z:
        :return:
        """
        gamma = self.get_lorentz_factor()
        wavelength = self.LAMBDAU / (2.0*gamma **2) * (1 + self.K**2 / 2.0 + gamma**2 * (theta_x**2 + theta_z ** 2))
        return wavelength / harmonic_number

    def get_resonance_frequency(self, harmonic_number=1, theta_x=0, theta_z=0):
        """
        gets the resonance wavelength (in m)
        :param harmonic_number:
        :param theta_x:
        :param theta_z:
        :return:
        """
        return codata.c / self.get_resonance_wavelength(harmonic_number,theta_x,theta_z)

    def get_resonance_energy(self, harmonic_number=1, theta_x=0, theta_z=0):
        """
        gets the resonance photon energy (in eV)
        :param harmonic_number:
        :param theta_x:
        :param theta_z:
        :return:
        """
        return codata.h*codata.c/codata.e*1e10 / (1e10*self.get_resonance_wavelength(harmonic_number, theta_x, theta_z))

    def get_resonance_central_cone(self,harmonic_number=1):
        """
        gets the angular aperture of the central cone (in rad).
        It is calculated as:
            CC = (1/gamma) * sqrt(  (1+0.5*K^2) / ( 2 N n) )
            with N = number of periods, n = harmonic number

            Note that the standard deviation is ~ 0.68*CC and the FWHM is 2.35*0.68*CC

        :param harmonic_number:
        :return:
        """
        return 1.0/self.get_lorentz_factor()*numpy.sqrt( (1+0.5*self.K**2)/(2*self.NPERIODS*harmonic_number) )

    def get_resonance_ring(self,harmonic_number=1, ring_order=1):
        return 1.0/self.get_lorentz_factor()*numpy.sqrt( ring_order / harmonic_number * (1+0.5*self.K**2) )

    # TODO: remove in far future
    def get_shadow3_source_object(self):
        """

        creates a Shadow.Source object with the undulator parameters inside (proprocessor file: xshundul.sha)

        :return:
        """
        # initialize shadow3 source (oe0) and beam
        h = self.to_dictionary()

        oe0 = Shadow.Source()

        if self.FLAG_EMITTANCE:
            oe0.EPSI_X = h["SX"] * h["SXP"]
            oe0.EPSI_Z = h["SZ"] * h["SZP"]
            oe0.SIGDIX = 0.0
            oe0.SIGDIZ = 0.0
            oe0.SIGMAX = h["SX"]
            oe0.SIGMAY = 0.0
            oe0.SIGMAZ = h["SZ"]
        else:
            oe0.EPSI_X = 0.0
            oe0.EPSI_Z = 0.0
            oe0.SIGDIX = 0.0
            oe0.SIGDIZ = 0.0
            oe0.SIGMAX = 0.0
            oe0.SIGMAY = 0.0
            oe0.SIGMAZ = 0.0

        oe0.FILE_TRAJ = b'xshundul.sha'
        oe0.ISTAR1 = h["SEED"]
        oe0.NPOINT = h["NRAYS"]
        oe0.F_WIGGLER = 2

        return oe0


    def set_from_dictionary(self,h):
        """
        set undulator variables from a dictionary. The following keys must be present:
          E_ENERGY
          INTENSITY
          SX
          SZ
          SXP
          SZP
          FLAG_EMITTANCE
          LAMBDAU
          NPERIODS
          K
          EMIN
          EMAX
          NG_E
          MAXANGLE
          NG_T
          NG_P
          N_J
          SEED
          NRAYS

        :param h:
        :return:
        """
        self.E_ENERGY        = h["E_ENERGY"]
        self.INTENSITY       = h["INTENSITY"]
        self.SX              = h["SX"]
        self.SZ              = h["SZ"]
        self.SXP             = h["SXP"]
        self.SZP             = h["SZP"]
        self.FLAG_EMITTANCE  = h["FLAG_EMITTANCE"]
        self.LAMBDAU         = h["LAMBDAU"]
        self.NPERIODS        = h["NPERIODS"]
        self.K               = h["K"]
        self.EMIN            = h["EMIN"]
        self.EMAX            = h["EMAX"]
        self.NG_E            = h["NG_E"]
        self.MAXANGLE        = h["MAXANGLE"]
        self.NG_T            = h["NG_T"]
        self.NG_P            = h["NG_P"]
        self.N_J             = h["N_J"]
        self.SEED            = h["SEED"]
        self.NRAYS           = h["NRAYS"]

    #
    # load/writers
    #
    def load_json_shadowvui_file(self,inFileTxt):
        """
        sets undulator parameters from a json file from shadowVUI and

        :param inFileTxt:
        :return:
        """
        #TODO: remove this routine in the future
        #
        # read inputs from a file created by ShadowVUI ----------------------------
        #

        with open(inFileTxt, mode='r') as f1:
            h = json.load(f1)
        self.load_json_shadowvui_dictionary(h)


    def load_json_shadowvui_dictionary(self,h):
        """
        sets undulator parameters from a dictionary loaded from a shadowVUI json file
        :param h:
        :return:
        """
        # self.dict = {
        #     "LAMBDAU":     0.0320000015,
        #     "K":      0.250000000,
        #     "E_ENERGY":       6.03999996,
        #     "E_ENERGY_SPREAD":    0.00100000005,
        #     "NPERIODS": 50,
        #     "EMIN":       10498.0000,
        #     "EMAX":       10499.0000,
        #     "INTENSITY":      0.200000003,
        #     "MAXANGLE":      0.100000001,
        #     "NG_E": 101,
        #     "NG_T": 51,
        #     "NG_P": 11,
        #     "NG_PLOT(1)":"1",
        #     "NG_PLOT(2)":"No",
        #     "NG_PLOT(3)":"Yes",
        #     "UNDUL_PHOT_FLAG(1)":"4",
        #     "UNDUL_PHOT_FLAG(2)":"Shadow code",
        #     "UNDUL_PHOT_FLAG(3)":"Urgent code",
        #     "UNDUL_PHOT_FLAG(4)":"SRW code",
        #     "UNDUL_PHOT_FLAG(5)":"Gaussian Approximation",
        #     "UNDUL_PHOT_FLAG(6)":"ESRF python code",
        #     "SEED": 36255,
        #     "SX":     0.0399999991,
        #     "SZ":    0.00100000005,
        #     "EX":   4.00000005E-07,
        #     "EZ":   3.99999989E-09,
        #     "FLAG_EMITTANCE(1)":"0",
        #     "FLAG_EMITTANCE(2)":"No",
        #     "FLAG_EMITTANCE(3)":"Yes",
        #     "NRAYS": 15000,
        #     "F_BOUND_SOUR": 0,
        #     "FILE_BOUND":"NONESPECIFIED",
        #     "SLIT_DISTANCE":       1000.00000,
        #     "SLIT_XMIN":      -1.00000000,
        #     "SLIT_XMAX":       1.00000000,
        #     "SLIT_ZMIN":      -1.00000000,
        #     "SLIT_ZMAX":       1.00000000,
        #     "NTOTALPOINT": 10000000,
        #     "JUNK4JSON":0
        #     }

        self.E_ENERGY        = h["E_ENERGY"]
        self.N_J             = 20
        self.INTENSITY       = h["INTENSITY"]
        self.SX              = h["SX"]
        self.SZ              = h["SZ"]
        self.SXP             = h["EX"] / h["SX"]
        self.SZP             = h["EZ"] / h["SZ"]
        self.FLAG_EMITTANCE  = int(h["FLAG_EMITTANCE(1)"])
        self.LAMBDAU         = h["LAMBDAU"]
        self.NPERIODS        = h["NPERIODS"]
        self.K               = h["K"]
        self.EMIN            = h["EMIN"]
        self.EMAX            = h["EMAX"]
        self.NG_E            = h["NG_E"]
        self.MAXANGLE        = h["MAXANGLE"]
        self.NG_T            = h["NG_T"]
        self.NG_P            = h["NG_P"]
        self.SEED            = h["SEED"]
        self.NRAYS           = h["NRAYS"]


    def load(self,inFileTxt):
        """
        load undulator parameters from a json file (type startj.00)
        :param inFileTxt:
        :return:
        """
        with open(inFileTxt, mode='r') as f1:
            h = json.load(f1)
        self.set_from_dictionary(h)

    def write(self,file_out='startj.00'):
        """
        write undulator parameters to a json file (type startj.00)
        :param file_out:
        :return:
        """
        data = self.to_dictionary()
        with open(file_out, 'w') as outfile:
            json.dump(data, outfile, indent=4, sort_keys=True, separators=(',', ':'))

    #
    # util
    #

    def sourcinfo(self,title=None):
        '''
        mimics SHADOW sourcinfo postprocessor. Returns a text array.
        :return: a text string
        '''

        txt = ''
        TOPLIN = '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'


        txt += TOPLIN
        txt += '**************  S O U R C E       D E S C R I P T I O N  **************\n'
        if title == None:
            txt += '\n\n'
        else:
            txt += title+'\n'
        txt += TOPLIN

        # TODO improve look
        txt += self.info()

        txt += TOPLIN
        txt += '***************                 E N D                  ***************\n'
        txt += TOPLIN
        return (txt)


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
        txt += "        Electron energy: %f geV\n"%self.E_ENERGY
        txt += "        Electron current: %f A\n"%self.INTENSITY
        if self.FLAG_EMITTANCE:
            txt += "        Electron sigmaX: %g [user units]\n"%self.SX
            txt += "        Electron sigmaZ: %g [user units]\n"%self.SZ
            txt += "        Electron sigmaX': %f urad\n"%(1e6*self.SXP)
            txt += "        Electron sigmaZ': %f urad\n"%(1e6*self.SZP)
        txt += "Input Undulator parameters: \n"
        txt += "        period: %f m\n"%self.LAMBDAU
        txt += "        number of periods: %d\n"%self.NPERIODS
        txt += "        K-value: %f\n"%self.K

        txt += "-----------------------------------------------------\n"

        txt += "Lorentz factor (gamma): %f\n"%self.get_lorentz_factor()
        txt += "Electron velocity: %.12f c units\n"%(numpy.sqrt(1.0 - 1.0 / self.get_lorentz_factor() ** 2))
        txt += "Undulator length: %f m\n"%(self.LAMBDAU*self.NPERIODS)
        K_to_B = (2.0 * numpy.pi / self.LAMBDAU) * codata.m_e * codata.c / codata.e

        txt += "Undulator peak magnetic field: %f T\n"%(K_to_B*self.K)
        txt += "Resonances: \n"
        txt += "        harmonic number [n]                   %10d %10d %10d \n"%(1,3,5)
        txt += "        wavelength [A]:                       %10.6f %10.6f %10.6f   \n"%(\
                                                                1e10*self.get_resonance_wavelength(1),
                                                                1e10*self.get_resonance_wavelength(3),
                                                                1e10*self.get_resonance_wavelength(5))
        txt += "        energy [eV]   :                       %10.3f %10.3f %10.3f   \n"%(\
                                                                self.get_resonance_energy(1),
                                                                self.get_resonance_energy(3),
                                                                self.get_resonance_energy(5))
        txt += "        frequency [Hz]:                       %10.3g %10.3g %10.3g   \n"%(\
                                                                self.get_resonance_frequency(1),
                                                                self.get_resonance_frequency(3),
                                                                self.get_resonance_frequency(5))
        txt += "        central cone 'half' width [mrad]:     %10.6f %10.6f %10.6f   \n"%(\
                                                                1e3*self.get_resonance_central_cone(1),
                                                                1e3*self.get_resonance_central_cone(3),
                                                                1e3*self.get_resonance_central_cone(5))
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
        if debug:
            txt += "E_ENERGY        = "+repr(self.E_ENERGY        ) + "\n"
            txt += "INTENSITY       = "+repr(self.INTENSITY       ) + "\n"
            txt += "SX              = "+repr(self.SX              ) + "\n"
            txt += "SZ              = "+repr(self.SZ              ) + "\n"
            txt += "SXP             = "+repr(self.SXP             ) + "\n"
            txt += "SZP             = "+repr(self.SZP             ) + "\n"
            txt += "  SX*SXP        = "+repr(self.SX * self.SXP   ) + "\n"
            txt += "  SZ*SZP        = "+repr(self.SZ * self.SZP   ) + "\n"
            txt += "FLAG_EMITTANCE  = "+repr(self.FLAG_EMITTANCE  ) + "\n"
            txt += "LAMBDAU         = "+repr(self.LAMBDAU         ) + "\n"
            txt += "NPERIODS        = "+repr(self.NPERIODS        ) + "\n"
            txt += "K               = "+repr(self.K               ) + "\n"
            txt += "EMIN            = "+repr(self.EMIN            ) + "\n"
            txt += "EMAX            = "+repr(self.EMAX            ) + "\n"
            txt += "NG_E            = "+repr(self.NG_E            ) + "\n"
            txt += "MAXANGLE        = "+repr(self.MAXANGLE        ) + "\n"
            txt += "NG_T            = "+repr(self.NG_T            ) + "\n"
            txt += "NG_P            = "+repr(self.NG_P            ) + "\n"
            txt += "N_J             = "+repr(self.N_J ) + "\n"
            txt += "SEED            = "+repr(self.SEED            ) + "\n"
            txt += "NRAYS           = "+repr(self.NRAYS           ) + "\n"
            txt += "-----------------------------------------------------\n"

        return txt

    def _shadow3_commands(self,commands="exit\n",input_file="shadow3_tmp.inp"):
        # for internal use
        f = open(input_file,'w')
        f.write(commands)
        f.close()
        os.system(self.SHADOW3_BINARY+" < "+input_file)

    #
    # call calculators
    #
    def calculate_shadow3_beam_using_preprocessors(self):
        """
        Calculates rays (source) using preprocessors in binary shadow3 (for
            comparison purposes, from python use calculate_beam() instead)

        :return:  a Shadow.Beam object
        """
        jsn = self.to_dictionary()
        # jsn = self.dict
        # print(self.info())
        os.system("rm -f start.00 systemfile.dat begin.dat xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")

        # epath
        commands = "epath\n2\n%f \n%f \n%f \n%d \n1.\nxshundul.par\nxshundul.traj\n1\nxshundul.plt\nexit\n"% \
        (jsn["LAMBDAU"],jsn["K"],jsn["E_ENERGY"],101)
        self._shadow3_commands(commands=commands,input_file="shadow3_epath.inp")

        # undul_set
        NG_E = jsn["NG_E"]
        # TODO: is seems a bug in shadow3: undul_set must be NG_E>1 (otherwise nosense)
        # but if emin=emaxthe  resulting uphot.nml has NG_E=1
        if NG_E == 1: NG_E = 2

        commands = "undul_set\n0\n0\n%d \n%d \n%d \nxshundul.traj\n%d\n%f\n%f\n%f\n%f\n0\n1000\nexit\n"% \
            (NG_E,jsn["NG_T"],jsn["NG_P"],
             jsn["NPERIODS"],jsn["EMIN"],jsn["EMAX"],jsn["INTENSITY"],jsn["MAXANGLE"])
        self._shadow3_commands(commands=commands,input_file="shadow3_undul_set.inp")

        # undul_phot
        self._shadow3_commands(commands="undul_phot\nexit\n",input_file="shadow3_undul_phot.inp")
        self._shadow3_commands(commands="undul_phot_dump\nexit\n",input_file="shadow3_undul_phot_dump.inp")

        # undul_cdf
        self._shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
                        input_file="shadow3_undul_cdf.inp")


        # input source
        if self.FLAG_EMITTANCE:
            commands = "input_source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
            (jsn["NRAYS"],jsn["SEED"],jsn["SX"],jsn["SZ"],jsn["SX"]*jsn["SXP"],0,jsn["SZ"]*jsn["SZP"],0,3,1,1)
        else:
            commands = "input_source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
            (jsn["NRAYS"],jsn["SEED"],0,0,0,0,0,0,3,1,1)

        self._shadow3_commands(commands=commands,input_file="shadow3_input_source.inp")

        # run source
        commands = "source\nsystemfile\nexit\n"
        self._shadow3_commands(commands=commands,input_file="shadow3_source.inp")


        # return shadow3 beam
        beam = Shadow.Beam()
        beam.load("begin.dat")
        return beam

    def calculate_radiation(self,code_undul_phot='internal'):
        """
        Calculates the radiation (emission) as a function pf theta (elevation angle) and phi (azimuthal angle)
        This radiation will be sampled to create the source

        It calls undul_phot* in SourceUndulatorFactory

        :param code_undul_phot: 'internal' (calls undul_phot), 'pysru' (calls undul_phot_pysru) or
                'srw' (calls undul_phot_srw)
        :return: a dictionary (the output from undul_phot*)
        """

        h = self.to_dictionary()
        # print(self.info())
        # os.system("rm -f xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")

        # if code_undul_phot != "internal" or code_undul_phot != "srw":
        #     dump_uphot_dot_dat = True



        # undul_phot
        if code_undul_phot == 'internal':
            undul_phot_dict = undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                    LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                    EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                    MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                    NG_P = h["NG_P"],
                                    number_of_trajectory_points=h["N_J"])
        elif code_undul_phot == 'pysru':
            undul_phot_dict = undul_phot_pysru(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                    LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                    EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                    MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                    NG_P = h["NG_P"])
        elif code_undul_phot == 'srw':
            undul_phot_dict = undul_phot_srw(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                    LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                    EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                    MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                    NG_P = h["NG_P"])
        else:
            raise Exception("Not implemented undul_phot code: "+code_undul_phot)


        return undul_phot_dict


    def calculate_cdf(self,code_undul_phot='internal',use_existing_undul_phot_output=None,
                       dump_undul_phot_file=False):
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

        if use_existing_undul_phot_output is None:
            undul_phot_dict = self.calculate_radiation(code_undul_phot=code_undul_phot)
        else:
            if isinstance(use_existing_undul_phot_output,str):
                undul_phot_dict = load_file_undul_phot(use_existing_undul_phot_output)
            elif isinstance(use_existing_undul_phot_output,dict):
                undul_phot_dict = use_existing_undul_phot_output
                #TODO: import parameters from external file E_MIN, E_MAX, MAXANGLE, N_*
            else:
                raise Exception("Bad undul_phot data.")



        if dump_undul_phot_file:
            write_file_undul_phot(undul_phot_dict,file_out="uphot.dat")



        # undul_cdf

        undul_cdf_dict = undul_cdf(undul_phot_dict,method='trapz')
        write_file_undul_sha(undul_cdf_dict,file_out="xshundul.sha",)

        return undul_cdf_dict



    def calculate_shadow3_beam(self,code_undul_phot='internal',use_existing_undul_phot_output=None,
                       dump_undul_phot_file=False,dump_start_files=False):
        """
        Calculates rays (source)

        :param code_undul_phot: 'internal' (calls undul_phot), 'pysru' (calls undul_phot_pysru) or
                'srw' (calls undul_phot_srw)
        :param use_existing_undul_phot_output: set to a file name or dictionary to use this
            particular output from undul_phot
        :param dump_undul_phot_file: if True writes uphot.dat with output from undul_phot
        :param dump_start_files: if Truie writes start.00 and end.00
        :return: a Shadow.Beam object
        """

        # create preprocessor file xshundul.sha
        tmp = self.calculate_cdf(code_undul_phot=code_undul_phot,
                                 use_existing_undul_phot_output=use_existing_undul_phot_output,
                                 dump_undul_phot_file=dump_undul_phot_file)

        # initialize shadow3 source (oe0) and beam
        oe0 = self.get_shadow3_source_object()

        if dump_start_files: oe0.write("start.00")
        beam = Shadow.Beam()
        beam.genSource(oe0)
        if dump_start_files: oe0.write("end.00")

        return beam


