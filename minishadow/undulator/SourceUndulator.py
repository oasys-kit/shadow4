
import numpy
import scipy.constants as codata
import json
import os

from numpy.testing import assert_equal, assert_almost_equal

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


def load_uphot_dot_dat(file_in="uphot.dat"):
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
    T = numpy.zeros(NG_T)
    P = numpy.zeros(NG_P)


    # for e in E:
    #     f.write("%g \n"%(e))
    #
    itmp = 0
    for ie,e in enumerate(E):
        E[ie] = tmp[itmp]
        itmp += 1

    for ie in range(NG_E):
        for it in range(NG_T):
            T[it] = tmp[itmp]
            itmp += 1

    for ie in range(NG_E):
        for it in range(NG_T):
            for ip in range(NG_P):
                P[ip] = tmp[itmp]
                itmp += 1

    Z2 = numpy.zeros((NG_E,NG_T,NG_P))
    Z3 = numpy.zeros((NG_E,NG_T,NG_P))

    for e in range(E.size):
        for t in range(T.size):
            for p in range(P.size):
                Z2[e,t,p] = tmp[itmp]
                itmp += 1

    for e in range(E.size):
        for t in range(T.size):
            for p in range(P.size):
                Z3[e,t,p] = tmp[itmp]
                itmp += 1


    return Z2,Z3,E,T,P


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
    Z2,Z3,E,T,P = load_uphot_dot_dat()
    print("Z2,Z3",Z2.shape,Z3.shape)
    print("E",E)
    print("T",T)
    print("P",P)

    from srxraylib.plot.gol import plot_image
    plot_image(Z2[0,:,:],T*1e6,P,xtitle="Theta [urad]",ytitle="Phi [rad]",aspect='auto')