
import numpy
import json
import os

from numpy.testing import assert_equal, assert_almost_equal

from srxraylib.plot.gol import plot,plot_image,plot_show


import Shadow
from SourceUndulatorFactory import undul_cdf, undul_phot, undul_phot_srw
from SourceUndulatorInputOutput import load_uphot_dot_dat,write_uphot_dot_dat
from SourceUndulatorInputOutput import load_xshundul_dot_sha,write_xshundul_dot_sha

import platform


class SourceUndulator(object):
    def __init__(self):

        # Machine
        self.E_ENERGY        = 6.04
        self.E_ENERGY_SPREAD = 0.001
        self.INTENSITY       = 0.2
        self.SX              = 0.0399999991
        self.SZ              = 0.00100000005
        self.EX              = 4.00000005E-07
        self.EZ              = 3.99999989E-09
        self.FLAG_EMITTANCE  = 1 # Yes
        # Undulator
        self.LAMBDAU         = 0.0320000015
        self.NPERIODS        = 50
        self.K               = 0.250000000
        # Photon energy scan
        self.EMIN            = 10498.0000
        self.EMAX            = 10499.0000
        self.NG_E            = 101
        # Geometry
        self.MAXANGLE        = 0.1
        self.NG_T            = 51
        self.NG_P            = 11
        # ray tracing
        self.SEED            = 36255
        self.NRAYS           = 15000


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
        h["EX"]              = self.EX
        h["EZ"]              = self.EZ
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
                SX=   0.0399999991,SZ=  0.00100000005,EX= 4.00000005E-07,EZ= 3.99999989E-09,FLAG_EMITTANCE=1,
                LAMBDAU=0.032,NPERIODS=50,K=0.25,
                EMIN= 10498.0000,EMAX= 10499.0000,NG_E=101,MAXANGLE=0.1,NG_T=51,NG_P=11,
                SEED=36255,NRAYS=15000,):
        self.E_ENERGY          = E_ENERGY
        self.E_ENERGY_SPREAD   = E_ENERGY_SPREAD
        self.INTENSITY         = INTENSITY
        self.SX                = SX
        self.SZ                = SZ
        self.EX                = EX
        self.EZ                = EZ
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
        self.EX              = h["EX"]
        self.EZ              = h["EZ"]
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
        self.EX              = h["EX"]
        self.EZ              = h["EZ"]
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

    def write(self,file_out):
        pass

    def info(self):
        # list all non-empty keywords
        txt = "-----------------------------------------------------\n"
        txt += "E_ENERGY        = "+repr(self.E_ENERGY        ) + "\n"
        txt += "E_ENERGY_SPREAD = "+repr(self.E_ENERGY_SPREAD ) + "\n"
        txt += "INTENSITY       = "+repr(self.INTENSITY       ) + "\n"
        txt += "SX              = "+repr(self.SX              ) + "\n"
        txt += "SZ              = "+repr(self.SZ              ) + "\n"
        txt += "EX              = "+repr(self.EX              ) + "\n"
        txt += "EZ              = "+repr(self.EZ              ) + "\n"
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
        os.system("rm -f xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")

        # epath
        commands = "epath\n2\n%f \n%f \n%f \n%d \n1.\nxshundul.par\nxshundul.traj\n1\nxshundul.plt\nexit\n"% \
        (jsn["LAMBDAU"],jsn["K"],jsn["E_ENERGY"],jsn["NG_E"])
        self.shadow3_commands(commands=commands,input_file="shadow3_epath.inp")

        # undul_set
        commands = "undul_set\n0\n0\n%d \n%d \n%d \nxshundul.traj\n%d\n%f\n%f\n%f\n%f\n0\n1000\nexit\n"% \
            (jsn["NG_E"],jsn["NG_T"],jsn["NG_P"],
             jsn["NPERIODS"],jsn["EMIN"],jsn["EMAX"],jsn["INTENSITY"],jsn["MAXANGLE"])
        self.shadow3_commands(commands=commands,input_file="shadow3_undul_set.inp")

        # undul_phot
        self.shadow3_commands(commands="undul_phot\nexit\n",input_file="shadow3_undul_phot.inp")
        self.shadow3_commands(commands="undul_phot_dump\nexit\n",input_file="shadow3_undul_phot_dump.inp")

        # undul_cdf
        self.shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
                        input_file="shadow3_undul_cdf.inp")


        # input source
        commands = "input source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
        (jsn["NRAYS"],jsn["SEED"],jsn["SX"],jsn["SZ"],jsn["EX"],0,jsn["EZ"],0,3,1,1)


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
        print(self.info())
        os.system("rm -f xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")

        if code_undul_phot != "internal" or code_undul_phot != "srw":
            dump_uphot_dot_dat = True



        # undul_phot
        if code_undul_phot == 'internal':
            undul_phot_dict = undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
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
            write_uphot_dot_dat(undul_phot_dict,file_out="uphot.dat")


        # undul_cdf

        undul_cdf_dict = undul_cdf(undul_phot_dict,method='trapz',do_plot=False)
        write_xshundul_dot_sha(undul_cdf_dict,file_out="xshundul.sha",)

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
        oe0.EPSI_X = h["EX"]
        oe0.EPSI_Z = h["EZ"]
        oe0.SIGDIX = 0.0
        oe0.SIGDIZ = 0.0
        oe0.SIGMAX = h["SX"]
        oe0.SIGMAY = 0.0
        oe0.SIGMAZ = h["SZ"]
        oe0.FILE_TRAJ = b'xshundul.sha'
        oe0.ISTAR1 = h["SEED"]
        oe0.NPOINT = h["NRAYS"]
        oe0.F_WIGGLER = 2
        if dump_start_files: oe0.write("start.00")
        beam.genSource(oe0)
        if dump_start_files: oe0.write("end.00")
        # if dump_begin_dot_dat: beam.write("begin.dat")

        return beam


def compare_shadow3_files(file1,file2,do_assert=True):

    Shadow.ShadowTools.plotxy(file1,4,6,nbins=101,nolost=1,title=file1)
    Shadow.ShadowTools.plotxy(file2,4,6,nbins=101,nolost=1,title=file2)

    if do_assert:
        begin1 = Shadow.Beam()
        begin1.load(file1)
        begin2     = Shadow.Beam()
        begin2.load(file2)
        assert_almost_equal(begin1.rays[:,0:6],begin2.rays[:,0:6],3)



if __name__ == "__main__":

    u = SourceUndulator()

    u.load_json_shadowvui_file("xshundul.json")
    print(u.to_dictionary())
    print(u.info())

    #
    # run using binary shadow3 (with preprocessors)
    #
    u.run_using_preprocessors()

    os.system("mv begin.dat begin_shadow3.dat")
    os.system("mv uphot.dat uphot_shadow3.dat")
    os.system("mv xshundul.sha xshundul_shadow3.sha")
    # make script
    # oe0 = Shadow.Source()
    # oe0.load("start.00")
    # Shadow.ShadowTools.make_python_script_from_list([oe0],script_file="script_undulator.py")

    #
    # run using python
    #
    beam = u.run(code_undul_phot='internal',
            dump_uphot_dot_dat=True,dump_start_files=True)

    beam.write("begin.dat")
    os.system("mv begin.dat begin_minishadow.dat")
    os.system("mv uphot.dat uphot_minishadow.dat")
    os.system("mv xshundul.sha xshundul_minishadow.sha")


    #
    # compare radiation
    #

    #
    # load_uphot_dot_dat("uphot_shadow3.dat",do_plot=True)
    # load_uphot_dot_dat("uphot_minishadow.dat",do_plot=True)
    #

    #
    # compare CDF
    #
    # load_xshundul_dot_sha("xshundul_shadow3.sha",do_plot=True)
    # load_xshundul_dot_sha("xshundul_minishadow.sha",do_plot=True)

    #
    # compare binary files
    #
    compare_shadow3_files("begin_shadow3.dat","begin_minishadow.dat",do_assert=True)


    plot_show()