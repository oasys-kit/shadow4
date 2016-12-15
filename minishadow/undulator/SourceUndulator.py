
import numpy
import scipy.constants as codata
import json
import os

from numpy.testing import assert_equal, assert_almost_equal
from srxraylib.plot.gol import plot,plot_image,plot_show

import Shadow

SHADOW3_BINARY = "/Users/srio/Oasys/OASYS_VE/shadow3/shadow3"

angstroms_to_eV = codata.h*codata.c/codata.e*1e10

def read_json(myinput="xshundul.json"):
    #
    # read inputs from a file created by ShadowVUI ----------------------------
    #
    if isinstance(myinput,str):
        inFileTxt = myinput # "xshundul.json"
        with open(inFileTxt, mode='r') as f1:
            h = json.load(f1)

    elif isinstance(myinput,dict):
        h = myinput
    else:
        raise Exception("Unknown input")

    return h


def info_dictionary(h):

    # list all non-empty keywords
    txt = "-----------------------------------------------------\n"
    for i,j in h.items():
        if (j != None):
            txt += "%s = %s\n" %(i,j)
    txt += "-----------------------------------------------------\n"
    txt += "k: %f \n"%(h['K'])
    return txt



def epath_compare():
    a = numpy.loadtxt("xshundul.plt").T
    # TTT = ETOFZ(I)/C  # X/c  Y=0
    # BBB = 1.0D0 - EBETAZ(I)
    # WRITE(33,1010) EXOFZ(I), EBETAX(I), EZ(I), BBB, TTT
    # X, BetaX, Z, 1-betaZ, X/c
    from srxraylib.plot.gol import plot
    plot(a[2],a[0],xtitle="Z(along)",ytitle="X(Horizontal)",title="trajectory at the edges??")


def shadow3_commands(commands="exit\n",input_file="shadow3_tmp.inp"):

    f = open(input_file,'w')
    f.write(commands)
    f.close()

    os.system(SHADOW3_BINARY+" < "+input_file)



def compare_results(do_assert=True):
    Shadow.ShadowTools.plotxy("begin_shadow3.dat",3,6,nbins=101,nolost=1,title="Undulator (SHADOW)")

    Shadow.ShadowTools.plotxy("begin_minishadow.dat",3,6,nbins=101,nolost=1,title="Undulator (minishadow)")

    if do_assert:

        begin_shadow3 = Shadow.Beam()
        begin_shadow3.load("begin_shadow3.dat")
        begin_minishadow     = Shadow.Beam()
        begin_minishadow.load("begin_minishadow.dat")
        assert_almost_equal(begin_shadow3.rays[:,0:6],begin_minishadow.rays[:,0:6],3)


def load_uphot_dot_dat(file_in="uphot.dat",do_plot=False):
    #
    # read uphot.dat file (like in SHADOW undul_phot_dump)
    #
    f = open(file_in,'r')
    firstline = f.readline()
    f.close()

    NG_E,NG_T,NG_P = numpy.fromstring(firstline,dtype=int,sep=" ")
    print("NG_E,NG_T,NG_P  %d  %d  %d \n"%(NG_E,NG_T,NG_P ))

    tmp = numpy.loadtxt(file_in,skiprows=1)

    if tmp.size != 3*(NG_E*NG_T*NG_P)+NG_E+NG_E*NG_T:
        raise Exception("File not understood")

    E = numpy.zeros(NG_E)
    T = numpy.zeros((NG_E,NG_T))
    P = numpy.zeros((NG_E,NG_T,NG_P))


    # for e in E:
    #     f.write("%g \n"%(e))
    #
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

    if do_plot:
        # with anscissas
        # plot_image(RN0[0,:,:],T.flatten()[0:NG_T].copy()*1e6,P.flatten()[0:NG_P].copy(),title="RN0[0]",xtitle="Theta [urad]",ytitle="Phi [rad]",aspect='auto',show=0)
        # plot_image(POL_DEG[0,:,:],T.flatten()[0:NG_T].copy()*1e6,P.flatten()[0:NG_P].copy(),title="POL_DEG[0]",xtitle="Theta [urad]",ytitle="Phi [rad]",aspect='auto',show=1)

        # with indices
        plot_image(RN0[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="RN0[0]",xtitle="Theta [index]",ytitle="Phi [index]",aspect=None,show=0)
        plot_image(POL_DEG[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="POL_DEG[0]",xtitle="Theta [index]",ytitle="Phi [index]",aspect=None,show=0)


    return RN0, POL_DEG, E, T.flatten()[0:NG_T].copy(), P.flatten()[0:NG_P].copy()


def load_xshundun_dot_sha(file_in="xshundul.sha",do_plot=False):
    #
    # read uphot.dat file (like in SHADOW undul_phot_dump)
    #
    f = open(file_in,'r')
    firstline = f.readline()
    f.close()

    NG_E,NG_T,NG_P, IANGLE = numpy.fromstring(firstline,dtype=int,sep=" ")
    print("NG_E,NG_T,NG_P, IANGLE  %d  %d  %d %d \n"%(NG_E,NG_T,NG_P,IANGLE ))

    tmp = numpy.loadtxt(file_in,skiprows=1)

    if tmp.size != 2*(NG_E + NG_E*NG_T + NG_E*NG_T*NG_P) + NG_E*NG_T*NG_P:
        raise Exception("File not understood")

    E = numpy.zeros(NG_E)
    T = numpy.zeros((NG_E,NG_T))
    P = numpy.zeros((NG_E,NG_T,NG_P))


    # for e in E:
    #     f.write("%g \n"%(e))
    #
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


    plot(E,TWO,title="TWO %s"%file_in,xtitle="E",ytitle="TWO",show=0)
    plot_image(ONE,numpy.arange(NG_E),numpy.arange(NG_T),title="ONE %s "%file_in,xtitle="index Energy",ytitle="index Theta",show=0)
    plot_image(ZERO[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="ZERO[0] %s"%file_in,xtitle="index Theta",ytitle="index Phi",show=0)
    # plot_image(POL_DEGREE[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="POL_DEGREE[0]",xtitle="index Theta",ytitle="index Phi",show=0)


    return TWO,ONE,ZERO,E,T,P,POL_DEGREE

def undul_cdf(file_out=None,do_plot=False):
    #
    # read uphot.dat file (like in SHADOW undul_phot_dump)
    #
    RN0,POL_DEG,E,T,P = load_uphot_dot_dat(do_plot=False)


    NG_E,NG_T,NG_P = RN0.shape
    print("NG_E,NG_T,NG_P, %d  %d %d \n"%(NG_E,NG_T,NG_P))

    # TWO = numpy.zeros(NG_E)
    # ONE = numpy.zeros((NG_E,NG_T))
    # ZERO = numpy.zeros((NG_E,NG_T,NG_P))

    RN1 = RN0.sum(axis=2) * 2 * numpy.pi * (P[1]-P[0])  # RN1(e,t)
    RN2 = RN1.sum(axis=1) * (T[1]-T[0])  # RN2(e)



    ZERO  = numpy.cumsum(RN0,axis=2) # CDF(e,t,p)
    ONE   = numpy.cumsum(RN1,axis=1) # CDF(e,t)
    TWO   = numpy.cumsum(RN2)        # CDF(e)

    print("Shadow ZERO,ONE,TWO: ",ZERO.shape,ONE.shape,TWO.shape)

    if do_plot:
        plot(E,TWO,title="NEW TWO",xtitle="E",ytitle="TWO",show=0)
        plot_image(ONE,numpy.arange(NG_E),numpy.arange(NG_T),title="NEW ONE",xtitle="index Energy",ytitle="index Theta",show=0)
        plot_image(ZERO[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="NEW ZERO[0]",xtitle="index Theta",ytitle="index Phi",show=0)
        #plot_image(POL_DEGREE[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="POL_DEGREE[0]",xtitle="index Theta",ytitle="index Phi",show=0)



    #
    # create xshundul.sha file (like in SHADOW undul_cdf)
    #
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




    return TWO,ONE,ZERO,E,T,P


def undulator(code='shadow3'):
    jsn = read_json("xshundul.json")
    os.system("rm -f xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")
    print(info_dictionary(jsn))

    if code == 'shadow3':

        # epath

        commands = "epath\n2\n%f \n%f \n%f \n%d \n1.\nxshundul.par\nxshundul.traj\n1\nxshundul.plt\nexit\n"% \
        (jsn["LAMBDAU"],jsn["K"],jsn["E_ENERGY"],jsn["NG_E"])
        shadow3_commands(commands=commands,input_file="shadow3_epath.inp")

        # undul_set
        commands = "undul_set\n0\n0\n%d \n%d \n%d \nxshundul.traj\n%d\n%f\n%f\n%f\n%f\n0\n1000\nexit\n"% \
            (jsn["NG_E"],jsn["NG_T"],jsn["NG_P"],
             jsn["NPERIODS"],jsn["EMIN"],jsn["EMAX"],jsn["INTENSITY"],jsn["MAXANGLE"])
        shadow3_commands(commands=commands,input_file="shadow3_undul_set.inp")

        # undul_phot
        shadow3_commands(commands="undul_phot\nexit\n",input_file="shadow3_undul_phot.inp")
        shadow3_commands(commands="undul_phot_dump\nexit\n",input_file="shadow3_undul_phot_dump.inp")

        # undul_cdf


        shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
                        input_file="shadow3_undul_cdf.inp")


        # input source
        commands = "input source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
        (jsn["NRAYS"],jsn["SEED"],jsn["SX"],jsn["SZ"],jsn["EX"],0,jsn["EZ"],0,3,1,1)


        shadow3_commands(commands=commands,input_file="shadow3_input_source.inp")

        # run source
        commands = "source\nsystemfile\nexit\n"
        shadow3_commands(commands=commands,input_file="shadow3_source.inp")

        os.system("cp begin.dat begin_shadow3.dat")
    else:

        from undul_phot import undul_phot

        # epath

        commands = "epath\n2\n%f \n%f \n%f \n%d \n1.\nxshundul.par\nxshundul.traj\n1\nxshundul.plt\nexit\n"% \
        (jsn["LAMBDAU"],jsn["K"],jsn["E_ENERGY"],jsn["NG_E"])
        shadow3_commands(commands=commands,input_file="shadow3_epath.inp")


        # undul_phot
        undul_phot(jsn)
        shadow3_commands(commands="undul_phot_dump\nexit\n",input_file="shadow3_undul_phot_dump.inp")


        # undul_cdf


        shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
                        input_file="shadow3_undul_cdf.inp")


        # input source
        commands = "input source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
        (jsn["NRAYS"],jsn["SEED"],jsn["SX"],jsn["SZ"],jsn["EX"],0,jsn["EZ"],0,3,1,1)


        shadow3_commands(commands=commands,input_file="shadow3_input_source.inp")

        # run source
        commands = "source\nsystemfile\nexit\n"
        shadow3_commands(commands=commands,input_file="shadow3_source.inp")

        os.system("cp begin.dat begin_minishadow.dat")





if __name__ == "__main__":

    # undulator(code='shadow3')
    # undulator(code='minishadow')
    #
    # compare_results()

    # uphot.dat

    # RN0,POL_DEG,E,T,P = load_uphot_dot_dat(do_plot=True)
    # print("RN0,POL_DEG",RN0.shape,POL_DEG.shape)
    # print("E",E)
    # print("T",T)
    # print("P",P)
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
    TWO,ONE,ZERO,E,T,P = undul_cdf(do_plot=False,file_out="xshundul2.sha")
    load_xshundun_dot_sha(file_in="xshundul.sha",do_plot=True)
    load_xshundun_dot_sha(file_in="xshundul2.sha",do_plot=True)

    plot_show()