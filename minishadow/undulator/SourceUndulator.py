
import numpy
import scipy.constants as codata
import json
import os

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



def compare_results(do_assert=False):
    Shadow.ShadowTools.plotxy("begin.dat",3,6,nbins=101,nolost=1,title="Undulator (SHADOW)")



    if do_assert:

        minimirr = Shadow.Beam()
        minimirr.load("minimirr.01")
        mirr     = Shadow.Beam()
        mirr.load("mirr.01")
        assert_almost_equal(minimirr.rays[:,0:6],mirr.rays[:,0:6])


        ministar = Shadow.Beam()
        ministar.load("ministar.01")
        star     = Shadow.Beam()
        star.load("star.01")
        assert_almost_equal(ministar.rays[:,0:6],star.rays[:,0:6])

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

        compare_results()


if __name__ == "__main__":

    undulator(code='shadow3')
