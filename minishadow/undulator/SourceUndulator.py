
import numpy
import json
import os

from numpy.testing import assert_equal, assert_almost_equal

from srxraylib.plot.gol import plot,plot_image,plot_show


import Shadow
from SourceUndulatorFactory import undul_cdf, undul_phot, undul_phot_srw

import platform


class SourceUndulator(object):
    def __init__(self):
        self.dict = {
            "LAMBDAU":     0.0320000015,
            "K":      0.250000000,
            "E_ENERGY":       6.03999996,
            "E_ENERGY_SPREAD":    0.00100000005,
            "NPERIODS": 50,
            "EMIN":       10498.0000,
            "EMAX":       10499.0000,
            "INTENSITY":      0.200000003,
            "MAXANGLE":      0.100000001,
            "NG_E": 101,
            "NG_T": 51,
            "NG_P": 11,
            "NG_PLOT(1)":"1",
            "NG_PLOT(2)":"No",
            "NG_PLOT(3)":"Yes",
            "UNDUL_PHOT_FLAG(1)":"4",
            "UNDUL_PHOT_FLAG(2)":"Shadow code",
            "UNDUL_PHOT_FLAG(3)":"Urgent code",
            "UNDUL_PHOT_FLAG(4)":"SRW code",
            "UNDUL_PHOT_FLAG(5)":"Gaussian Approximation",
            "UNDUL_PHOT_FLAG(6)":"ESRF python code",
            "SEED": 36255,
            "SX":     0.0399999991,
            "SZ":    0.00100000005,
            "EX":   4.00000005E-07,
            "EZ":   3.99999989E-09,
            "FLAG_EMITTANCE(1)":"0",
            "FLAG_EMITTANCE(2)":"No",
            "FLAG_EMITTANCE(3)":"Yes",
            "NRAYS": 15000,
            "F_BOUND_SOUR": 0,
            "FILE_BOUND":"NONESPECIFIED",
            "SLIT_DISTANCE":       1000.00000,
            "SLIT_XMIN":      -1.00000000,
            "SLIT_XMAX":       1.00000000,
            "SLIT_ZMIN":      -1.00000000,
            "SLIT_ZMAX":       1.00000000,
            "NTOTALPOINT": 10000000,
            "JUNK4JSON":0
            }

        if platform.system() == "Linux":
            self.SHADOW3_BINARY = "/users/srio/OASYS_VE/shadow3/shadow3"
        else:
            self.SHADOW3_BINARY = "/Users/srio/Oasys/OASYS_VE/shadow3/shadow3"

    def to_dictionary(self):
        """
        returns a python dictionary of the Shadow.Source instance
        :return: a dictionary
        """
        return(self.dict)

    def duplicate(self):
        """
        makes a copy of the source
        :return: new instance of Shadow.Source()
        """
        pass


    def set_energy_monochromatic(self,emin):
        """
        Sets a single energy line for the source (monochromatic)
        :param emin: the energy in eV
        :return: self
        """
        pass

    def set_energy_box(self,emin,emax):
        """
        Sets a box energy distribution for the source (monochromatic)
        :param emin: minimum energy in eV
        :param emax: maximum energy in eV
        :return: self
        """
        pass


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

        txt += TOPLIN
        txt += '***************                 E N D                  ***************\n'
        txt += TOPLIN
        return (txt)


    def set_from_dictionary(self,myinput):
        #TODO copy one by one
        self.dict = myinput


    def load_json_shadowvui(self,inFileTxt):
        #
        # read inputs from a file created by ShadowVUI ----------------------------
        #
        with open(inFileTxt, mode='r') as f1:
            h = json.load(f1)
        self.dict = h

    def load(self,inFileTxt):
        pass

    def write(self,file_out):
        pass

    def info(self):

        # list all non-empty keywords
        txt = "-----------------------------------------------------\n"
        for i,j in self.dict.items():
            if (j != None):
                txt += "%s = %s\n" %(i,j)
        txt += "-----------------------------------------------------\n"
        txt += "k: %f \n"%(self.dict['K'])
        return txt

    def shadow3_commands(self,commands="exit\n",input_file="shadow3_tmp.inp"):
        f = open(input_file,'w')
        f.write(commands)
        f.close()
        os.system(self.SHADOW3_BINARY+" < "+input_file)

    def run_using_preprocessors(self):
        jsn = self.dict
        print(self.info())
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

    def run(self,code_undul_phot='internal',code_undul_cdf='internal',
            dump_uphot_dot_dat=False,dump_xshundul_dot_sha=False,
            dump_start_files=False,dump_begin_dot_dat=False):

        jsn = self.dict
        print(self.info())
        os.system("rm -f xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")

        if code_undul_phot != "internal" or code_undul_phot != "srw":
            dump_uphot_dot_dat = True

        if code_undul_cdf == "preprocessor":
            dump_uphot_dot_dat = True


        # undul_phot
        if code_undul_phot == 'internal':
            undul_phot_dict = undul_phot(jsn)
        elif code_undul_phot == 'srw':
            undul_phot_srw(jsn)
        else:
            raise Exception("Not implemented undul_phot code: "+code_undul_phot)

        if dump_uphot_dot_dat:
            write_uphot_dot_dat(undul_phot_dict,file_out="uphot.dat")


        # undul_cdf
        if code_undul_cdf == 'internal':
            # undul_phot_dict = load_uphot_dot_dat(file_in="uphot.dat")
            undul_cdf_dict = undul_cdf(undul_phot_dict,method='trapz',do_plot=False)
            write_xshundun_dot_sha(undul_cdf_dict,file_out="xshundul.sha",)

            if dump_xshundul_dot_sha:
                write_xshundun_dot_sha(undul_cdf_dict)

        elif code_undul_cdf == 'preprocessor':
            self.shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
                            input_file="shadow3_undul_cdf.inp")
        else:
            raise Exception("Not implemented undul_cdf code: "+code_undul_cdf)

        # run source

        # initialize shadow3 source (oe0) and beam
        oe0 = Shadow.Source()
        beam = Shadow.Beam()
        oe0.EPSI_X = jsn["EX"]
        oe0.EPSI_Z = jsn["EZ"]
        oe0.SIGDIX = 0.0
        oe0.SIGDIZ = 0.0
        oe0.SIGMAX = jsn["SX"]
        oe0.SIGMAY = 0.0
        oe0.SIGMAZ = jsn["SZ"]
        oe0.FILE_TRAJ = b'xshundul.sha'
        oe0.ISTAR1 = jsn["SEED"]
        oe0.NPOINT = jsn["NRAYS"]
        oe0.F_WIGGLER = 2
        if dump_start_files: oe0.write("start.00")
        beam.genSource(oe0)
        if dump_start_files: oe0.write("end.00")
        if dump_begin_dot_dat: beam.write("begin.dat")

        return beam

#
# auxiliar functoins
#


def load_uphot_dot_dat(file_in="uphot.dat",do_plot=False,verbose=True):
    """
    read uphot.dat file (like in SHADOW undul_phot_dump)

    :param file_in: name of the file to be read
    :param do_plot: flag to plot
    :param verbose: flag to print some info
    :return: a dictionary {'radiation':RN0, 'polarization':POL_DEG, 'photon_energy':E, 'theta':TT, 'phi':PP}
    """

    f = open(file_in,'r')
    firstline = f.readline()
    f.close()

    NG_E,NG_T,NG_P = numpy.fromstring(firstline,dtype=int,sep=" ")
    if verbose:
        print("\nload_uphot_dot_dat: file is %s"%(file_in))
        print("load_uphot_dot_dat: NG_E,NG_T,NG_P  %d  %d  %d"%(NG_E,NG_T,NG_P ))

    tmp = numpy.loadtxt(file_in,skiprows=1)

    if tmp.size != 3*(NG_E*NG_T*NG_P)+NG_E+NG_E*NG_T:
        raise Exception("load_uphot_dot_dat: File not understood")

    E = numpy.zeros(NG_E)
    T = numpy.zeros((NG_E,NG_T))
    P = numpy.zeros((NG_E,NG_T,NG_P))


    itmp = 0
    for ie in range(NG_E):
        E[ie] = tmp[itmp]
        itmp += 1

    for ie in range(NG_E):
        for it in range(NG_T):
            T[ie,it] = tmp[itmp]
            itmp += 1

    for ie in range(NG_E):
        for it in range(NG_T):
            for ip in range(NG_P):
                P[ie,it,ip] = tmp[itmp]
                itmp += 1

    RN0 = numpy.zeros((NG_E,NG_T,NG_P))
    POL_DEG = numpy.zeros((NG_E,NG_T,NG_P))

    for e in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                RN0[e,t,p] = tmp[itmp]
                itmp += 1

    for e in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                POL_DEG[e,t,p] = tmp[itmp]
                itmp += 1

    TT = T.flatten()[0:NG_T].copy()
    PP = P.flatten()[0:NG_P].copy()
    if verbose:
        print("load_uphot_dot_dat: Step in E: %f, Interval in E: %f"%( (E[1]-E[0]), (E[-1]-E[0]) ))
        print("load_uphot_dot_dat: Step in P: %f, Interval in P: %f"%( (PP[1]-PP[0]), (PP[-1]-PP[0]) ))
        print("load_uphot_dot_dat: Step in T: %f, Interval in T: %f"%( (TT[1]-TT[0]), (TT[-1]-TT[0]) ))

        print("load_uphot_dot_dat: RN0 max: %f min: %f"%(RN0.max(),RN0.min()) )
        print("load_uphot_dot_dat: POL_DEG max: %f min: %f"%(POL_DEG.max(),POL_DEG.min()) )

    if do_plot:
        plot_image(RN0[0,:,:],TT*1e6,PP*180/numpy.pi,title=file_in+" RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
        plot_image(POL_DEG[0,:,:],TT*1e6,PP*180/numpy.pi,title=file_in+" POL_DEG[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)

    return {'radiation':RN0, 'polarization':POL_DEG, 'photon_energy':E, 'theta':TT, 'phi':PP}



def write_uphot_dot_dat(undul_phot_dict,file_out="uphot.dat"):

    Z2      = undul_phot_dict['radiation']
    POL_DEG = undul_phot_dict['polarization']
    e       = undul_phot_dict['photon_energy']
    theta   = undul_phot_dict['theta']
    phi     = undul_phot_dict['phi']

    NG_E = e.size
    NG_T = theta.size
    NG_P = phi.size

    f = open(file_out,'w')
    f.write("%d  %d  %d \n"%(NG_E,NG_T,NG_P))
    for ie in range(NG_E):
        f.write("%20.10f \n"%(e[ie]))

    for ie in range(NG_E):
        for t in range(NG_T):
            f.write("%20.10f \n"%(theta[t]))

    for ie in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                f.write("%20.10f \n"%(phi[p]))


    for ie in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                f.write("%20.10f \n"%Z2[ie,t,p])

    for ie in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                f.write("%20.10f \n"%(POL_DEG[ie,t,p]))

    f.close()
    print("File written to disk: %s"%file_out)




def load_xshundun_dot_sha(file_in="xshundul.sha",do_plot=False,verbose=True):
    #
    # read uphot.dat file (like in SHADOW undul_phot_dump)
    #
    f = open(file_in,'r')
    firstline = f.readline()
    f.close()

    NG_E,NG_T,NG_P, IANGLE = numpy.fromstring(firstline,dtype=int,sep=" ")
    if verbose: print("NG_E,NG_T,NG_P, IANGLE  %d  %d  %d %d \n"%(NG_E,NG_T,NG_P,IANGLE ))

    tmp = numpy.loadtxt(file_in,skiprows=1)

    if tmp.size != 2*(NG_E + NG_E*NG_T + NG_E*NG_T*NG_P) + NG_E*NG_T*NG_P:
        raise Exception("File not understood")

    E = numpy.zeros(NG_E)
    T = numpy.zeros((NG_E,NG_T))
    P = numpy.zeros((NG_E,NG_T,NG_P))


    itmp = 0
    for ie,e in enumerate(E):
        E[ie] = tmp[itmp]
        itmp += 1

    for ie in range(NG_E):
        for it in range(NG_T):
            T[ie,it] = tmp[itmp]
            itmp += 1

    for ie in range(NG_E):
        for it in range(NG_T):
            for ip in range(NG_P):
                P[ie,it,ip] = tmp[itmp]
                itmp += 1

    TWO = numpy.zeros((NG_E))
    ONE = numpy.zeros((NG_E,NG_T))
    ZERO = numpy.zeros((NG_E,NG_T,NG_P))
    POL_DEGREE = numpy.zeros_like(ZERO)


    for e in range(NG_E):
        TWO[e] = tmp[itmp]
        itmp += 1

    for e in range(NG_E):
        for t in range(NG_T):
            ONE[e,t] = tmp[itmp]
            itmp += 1

    for e in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                ZERO[e,t,p] = tmp[itmp]
                itmp += 1

    for e in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                POL_DEGREE[e,t,p] = tmp[itmp]
                itmp += 1

    if do_plot:
        plot(E,TWO,title="TWO %s"%file_in,xtitle="E",ytitle="TWO",show=False)
        plot_image(ONE,numpy.arange(NG_E),numpy.arange(NG_T),title="ONE %s "%file_in,xtitle="index Energy",ytitle="index Theta",show=False)
        plot_image(ZERO[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="ZERO[0] %s"%file_in,xtitle="index Theta",ytitle="index Phi",show=False)
        # plot_image(POL_DEGREE[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="POL_DEGREE[0]",xtitle="index Theta",ytitle="index Phi",show=0)


    return {'cdf_EnergyThetaPhi':TWO,'cdf_EnergyTheta':ONE,'cdf_Energy':ZERO,'energy':E,'theta':T,'phi':P,'polarization':POL_DEGREE}

def write_xshundun_dot_sha(dict,file_out="xshundul.sha"):
    #
    # create xshundul.sha file (like in SHADOW undul_cdf)
    #
    TWO =     dict['cdf_EnergyThetaPhi']
    ONE =     dict['cdf_EnergyTheta']
    ZERO =    dict['cdf_Energy']
    E =       dict['energy']
    T =       dict['theta']
    P =       dict['phi']
    POL_DEG = dict['polarization']

    NG_E = E.size
    NG_T = T.size
    NG_P = P.size

    if file_out is not None:
        f = open(file_out,'w')
        f.write("%d  %d  %d 1 \n"%(NG_E,NG_T,NG_P))

        for e in E:
            f.write("%g \n"%(e))

        for e in E:
            for t in T:
                f.write("%g \n"%t)

        for e in E:
            for t in T:
                for p in P:
                    f.write("%g \n"%p)

        for e in numpy.arange(NG_E):
            f.write("%g \n"%(TWO[e]))

        for e in numpy.arange(NG_E):
            for t in numpy.arange(NG_T):
                f.write("%g \n"%(ONE[e,t]))

        for e in numpy.arange(NG_E):
            for t in numpy.arange(NG_T):
                for p in numpy.arange(NG_P):
                    f.write("%g \n"%(ZERO[e,t,p]))

        for e in numpy.arange(NG_E):
            for t in numpy.arange(NG_T):
                for p in numpy.arange(NG_P):
                    f.write("%g \n"%(POL_DEG[e,t,p]))


        f.close()
        print("File written to disk: %s"%file_out)


#
# comparisons/tests
#
# def epath_compare():
#     a = numpy.loadtxt("xshundul.plt").T
#     # TTT = ETOFZ(I)/C  # X/c  Y=0
#     # BBB = 1.0D0 - EBETAZ(I)
#     # WRITE(33,1010) EXOFZ(I), EBETAX(I), EZ(I), BBB, TTT
#     # X, BetaX, Z, 1-betaZ, X/c
#     plot(a[2],a[0],xtitle="Z(along)",ytitle="X(Horizontal)",title="trajectory at the edges??")


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

    u.load_json_shadowvui("xshundul.json")
    print(u.to_dictionary())
    print(u.info())


    #
    # run using binary shadow3 (with preprcessors)
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
    # run using python preprocessors
    #
    u.run(code_undul_phot='internal',code_undul_cdf='internal',
            dump_uphot_dot_dat=True,dump_xshundul_dot_sha=True,
            dump_start_files=True,dump_begin_dot_dat=True)

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
    # xshundul.sha
    # TWO,ONE,ZERO,E,T,P,POL_DEGREE = load_xshundun_dot_sha(do_plot=True)
    # print("Shadow TWO",TWO.shape)
    # print("Shadow ONE",ONE.shape)
    # print("Shadow ZERO",ZERO.shape)
    #
    # print("Total TWO",  TWO.sum())
    # print("Total ONE",  ONE.sum())
    # print("Total ZERO",ZERO.sum())
    # NG_E,NG_T,NG_P = ZERO.shape

    # new undul_cdf
    # TWO,ONE,ZERO,E,T,P,POL_DEGREE = undul_cdf(do_plot=True)


    #
    # compare binary files
    #
    compare_shadow3_files("begin_shadow3.dat","begin_minishadow.dat",do_assert=True)


    #===============================================================

    # uphot.dat

    # RN0,POL_DEG,E,T,P = load_uphot_dot_dat(do_plot=True)
    # print("RN0,POL_DEG",RN0.shape,POL_DEG.shape)
    # print("E",E)
    # print("T",T)
    # print("P",P)
    #
    plot_show()