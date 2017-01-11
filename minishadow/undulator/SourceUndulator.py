import json
import os
import numpy

import scipy.constants as codata

import Shadow
from SourceUndulatorFactory import undul_cdf, undul_phot, undul_phot_srw,  undul_phot_pysru
from SourceUndulatorInputOutput import load_file_undul_phot,write_file_undul_phot
from SourceUndulatorInputOutput import load_fule_undul_cdf,write_file_undul_sha

import platform


class SourceUndulator(object):
    def __init__(self):

        # Machine
        self.E_ENERGY        = 1.0
        self.E_ENERGY_SPREAD = 0.0
        self.INTENSITY       = 1.0
        self.SX              = 1e-6
        self.SZ              = 1e-6
        self.SXP             = 1e-6
        self.SZP             = 1e-6
        # self.EX              = 4.00000005E-07
        # self.EZ              = 3.99999989E-09
        self.FLAG_EMITTANCE  = 1 # Yes
        # Undulator
        self.LAMBDAU         = 0.01
        self.NPERIODS        = 100
        self.K               = 0.8
        # Photon energy scan
        self.EMIN            = 10000.0
        self.EMAX            = 11000.0
        self.NG_E            = 11
        # Geometry
        self.MAXANGLE        = 0.5
        self.NG_T            = 31
        self.NG_P            = 21
        # ray tracing
        self.SEED            = 36255655452
        self.NRAYS           = 5000


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
            self.SHADOW3_BINARY = "/users/srio/OASYS_VE/shadow3/shadow3"
        else:
            self.SHADOW3_BINARY = "/Users/srio/Oasys/OASYS_VE/shadow3/shadow3"

    def to_dictionary(self):
        """
        returns a python dictionary of the Shadow.Source instance
        :return: a dictionary
        """
        h = {}
        h["E_ENERGY"]        = self.E_ENERGY
        h["E_ENERGY_SPREAD"] = self.E_ENERGY_SPREAD
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

    def set_from_keywords(self,E_ENERGY=6.04,E_ENERGY_SPREAD=0.001,INTENSITY=0.2,
                SX=0.04,SZ=0.001,SXP=10e-6,SZP=4e-6,FLAG_EMITTANCE=1,
                LAMBDAU=0.032,NPERIODS=50,K=0.25,
                EMIN= 10498.0000,EMAX= 10499.0000,NG_E=101,MAXANGLE=0.1,NG_T=51,NG_P=11,
                SEED=36255,NRAYS=15000,):
        self.E_ENERGY          = E_ENERGY
        self.E_ENERGY_SPREAD   = E_ENERGY_SPREAD
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
        self.SEED              = SEED
        self.NRAYS             = NRAYS

    def set_energy_monochromatic(self,emin):
        """
        Sets a single energy line for the source (monochromatic)
        :param emin: the energy in eV
        :return: self
        """
        self.EMIN = emin
        self.EMAX = emin
        self.NG_E = 1

    def set_energy_monochromatic_at_resonance(self,harmonic_number):

        self.EMIN = self.get_resonance_energy(harmonic_number=harmonic_number)
        self.EMAX = self.EMIN
        self.NG_E = 1

    def set_energy_box(self,emin,emax):
        """
        Sets a box energy distribution for the source (monochromatic)
        :param emin: minimum energy in eV
        :param emax: maximum energy in eV
        :return: self
        """
        self.EMIN = emin
        self.EMAX = emax

    def set_energy_number_of_points(self,npoints):
        self.NG_E = npoints

    #
    # useful getters
    #
    def get_lorentz_factor(self):
        return ( self.E_ENERGY * 1e9) / (codata.m_e *  codata.c**2 / codata.e)

    def get_resonance_wavelength(self, harmonic_number=1, theta_x=0.0, theta_z=0.0):
        gamma = self.get_lorentz_factor()
        wavelength = self.LAMBDAU / (2.0*gamma **2) * (1 + self.K**2 / 2.0 + gamma**2 * (theta_x**2 + theta_z ** 2))
        return wavelength / harmonic_number

    def get_resonance_frequency(self, harmonic_number=1, theta_x=0, theta_z=0):
        return codata.c / self.get_resonance_wavelength(harmonic_number,theta_x,theta_z)

    def get_resonance_energy(self, harmonic_number=1, theta_x=0, theta_z=0):
        return codata.h*codata.c/codata.e*1e10 / (1e10*self.get_resonance_wavelength(harmonic_number, theta_x, theta_z))

    def get_resonance_central_cone(self,harmonic_number=1):
        return 1.0/self.get_lorentz_factor()*numpy.sqrt( (1+0.5*self.K**2)/(2*self.NPERIODS*harmonic_number) )

    def get_resonance_ring(self,harmonic_number=1, ring_order=1):
        return 1.0/self.get_lorentz_factor()*numpy.sqrt( ring_order / harmonic_number * (1+0.5*self.K**2) )

    #     codata = scipy.constants.codata.physical_constants
    #     codata_c = codata["speed of light in vacuum"][0]
    #
    #     frequency = codata_c / self.resonanceWavelength(gamma, theta_x, theta_z)
    #     return frequency
    #
    #
    #
    # def resonanceEnergy(self, gamma, theta_x, theta_y, harmonic=1):
    #     codata = scipy.constants.codata.physical_constants
    #     energy_in_ev = codata["Planck constant"][0] * self.resonanceFrequency(gamma, theta_x, theta_y) / codata["elementary charge"][0]
    #     return energy_in_ev*harmonic

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


    def set_from_dictionary(self,h):
        self.E_ENERGY        = h["E_ENERGY"]
        self.E_ENERGY_SPREAD = h["E_ENERGY_SPREAD"]
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
        self.SEED            = h["SEED"]
        self.NRAYS           = h["NRAYS"]


    def load_json_shadowvui_file(self,inFileTxt):
        #
        # read inputs from a file created by ShadowVUI ----------------------------
        #

        with open(inFileTxt, mode='r') as f1:
            h = json.load(f1)
        self.load_json_shadowvui_dictionary(h)


    def load_json_shadowvui_dictionary(self,h):
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
        self.E_ENERGY_SPREAD = h["E_ENERGY_SPREAD"]
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
        with open(inFileTxt, mode='r') as f1:
            h = json.load(f1)
        self.set_from_dictionary(h)

    def write(self,file_out='startj.00'):
        data = self.to_dictionary()
        with open(file_out, 'w') as outfile:
            json.dump(data, outfile, indent=4, sort_keys=True, separators=(',', ':'))


    def info(self,debug=False):
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
        txt += "        central cone half width [mrad]:       %10.6f %10.6f %10.6f   \n"%(\
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
        txt += "        number of energy points %d\n"%(self.NG_E)
        txt += "        number of angular elevation points %d\n"%(self.NG_T)
        txt += "        number of angular azimuthal points %d\n"%(self.NG_P)
        txt += "        number of rays %d\n"%(self.NRAYS)
        txt += "        random seed %d\n"%(self.SEED)

        txt += "-----------------------------------------------------\n"
        if debug:
            txt += "E_ENERGY        = "+repr(self.E_ENERGY        ) + "\n"
            txt += "E_ENERGY_SPREAD = "+repr(self.E_ENERGY_SPREAD ) + "\n"
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
            txt += "SEED            = "+repr(self.SEED            ) + "\n"
            txt += "NRAYS           = "+repr(self.NRAYS           ) + "\n"
            txt += "-----------------------------------------------------\n"

        return txt

    def shadow3_commands(self,commands="exit\n",input_file="shadow3_tmp.inp"):
        f = open(input_file,'w')
        f.write(commands)
        f.close()
        os.system(self.SHADOW3_BINARY+" < "+input_file)

    def run_using_preprocessors(self):
        jsn = self.to_dictionary()
        # jsn = self.dict
        # print(self.info())
        os.system("rm -f start.00 systemfile.dat begin.dat xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")

        # epath
        commands = "epath\n2\n%f \n%f \n%f \n%d \n1.\nxshundul.par\nxshundul.traj\n1\nxshundul.plt\nexit\n"% \
        (jsn["LAMBDAU"],jsn["K"],jsn["E_ENERGY"],101)
        self.shadow3_commands(commands=commands,input_file="shadow3_epath.inp")

        # undul_set
        NG_E = jsn["NG_E"]
        # TODO: is seems a bug in shadow3: undul_set must be NG_E>1 (otherwise nosense)
        # but if emin=emaxthe  resulting uphot.nml has NG_E=1
        if NG_E == 1: NG_E = 2

        commands = "undul_set\n0\n0\n%d \n%d \n%d \nxshundul.traj\n%d\n%f\n%f\n%f\n%f\n0\n1000\nexit\n"% \
            (NG_E,jsn["NG_T"],jsn["NG_P"],
             jsn["NPERIODS"],jsn["EMIN"],jsn["EMAX"],jsn["INTENSITY"],jsn["MAXANGLE"])
        self.shadow3_commands(commands=commands,input_file="shadow3_undul_set.inp")

        # undul_phot
        self.shadow3_commands(commands="undul_phot\nexit\n",input_file="shadow3_undul_phot.inp")
        self.shadow3_commands(commands="undul_phot_dump\nexit\n",input_file="shadow3_undul_phot_dump.inp")

        # undul_cdf
        self.shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
                        input_file="shadow3_undul_cdf.inp")


        # input source
        if self.FLAG_EMITTANCE:
            commands = "input_source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
            (jsn["NRAYS"],jsn["SEED"],jsn["SX"],jsn["SZ"],jsn["SX"]*jsn["SXP"],0,jsn["SZ"]*jsn["SZP"],0,3,1,1)
        else:
            commands = "input_source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
            (jsn["NRAYS"],jsn["SEED"],0,0,0,0,0,0,3,1,1)

        self.shadow3_commands(commands=commands,input_file="shadow3_input_source.inp")

        # run source
        commands = "source\nsystemfile\nexit\n"
        self.shadow3_commands(commands=commands,input_file="shadow3_source.inp")


        # return shadow3 beam
        beam = Shadow.Beam()
        beam.load("begin.dat")
        return beam

    def run(self,code_undul_phot='internal',dump_uphot_dot_dat=False,dump_start_files=False):

        h = self.to_dictionary()
        # print(self.info())
        # os.system("rm -f xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")

        if code_undul_phot != "internal" or code_undul_phot != "srw":
            dump_uphot_dot_dat = True



        # undul_phot
        if code_undul_phot == 'internal':
            undul_phot_dict = undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                    LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                    EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                    MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                    NG_P = h["NG_P"])
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

        if dump_uphot_dot_dat:
            write_file_undul_phot(undul_phot_dict,file_out="uphot.dat")


        # undul_cdf

        undul_cdf_dict = undul_cdf(undul_phot_dict,method='trapz',do_plot=False)
        write_file_undul_sha(undul_cdf_dict,file_out="xshundul.sha",)

        # # TODO: remove this cannot be done here as xshundul.traj?? is needed but not existing!!
        # elif code_undul_cdf == 'preprocessor':
        #     self.shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
        #                     input_file="shadow3_undul_cdf.inp")
        # else:
        #     raise Exception("Not implemented undul_cdf code: "+code_undul_cdf)

        # run source

        # initialize shadow3 source (oe0) and beam
        oe0 = Shadow.Source()
        beam = Shadow.Beam()
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

        if dump_start_files: oe0.write("start.00")
        beam.genSource(oe0)
        if dump_start_files: oe0.write("end.00")
        # if dump_begin_dot_dat: beam.write("begin.dat")

        return beam


